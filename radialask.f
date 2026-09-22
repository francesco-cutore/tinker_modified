c     ##############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong and Jay William Ponder  ##
c     ##              Modified by Francesco Cutore                ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine radialask - manages on-the-fly structure factor calculation  ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c
      subroutine radialask
      use argue
      use atomid
      use atoms
      use bound
      use boxes
      use files
      use inform
      use iounit
      use keys
      use limits
      use math
      use molcul
      use potent
      use rdfparams
      use factors

      implicit none
      integer slw                                     ! Sliding window size
      integer i, j, next
      real*8 rmax
      logical exist,query
      logical read_file
      integer freeunit
      integer unit
      character*20 keyword
      character*240 record
      character*240 string
      real*8 q_min
      real*8 kappa
      character*20 kstr

      real*8   q_max              ! Maximum q value
      real*8   q_spat             ! Spatial bin size change also q step
      real*8   q_threshold        ! Threshold for decimation
      real*8   decimation_power   ! Power for decimation

      logical, parameter :: debug_ra = .true.         ! Debug flag writing additional files
      logical, parameter :: verbose_ra = .true.       ! Verbose flag for printing info to console

      integer  normalization_mode ! Normalization mode for S(q) calculation refers
      real*8   manual_factor      ! to table on the bottom.
      logical  normmode_set       ! .true. if SCATTER-NORMMODE was set in keyfile
      logical  expminusone_set    ! .true. if SCATTER-EXPMINUSONE was set in keyfile

      integer :: num_unique
      logical :: is_new
      integer, allocatable :: unique_z(:)
c
c     Allocate global variables and set parameters
c
      verbose_global = verbose_ra
      debug_global = debug_ra
c
c     Default values for the scattering-restraint parameters, then
c     allow override via keyfile (SCATTER-QMAX, SCATTER-QSPAT,
c     SCATTER-QTHRESH, SCATTER-DECPOWER, SCATTER-NORMMODE,
c     SCATTER-MANUALFAC, SCATTER-EXPMINUSONE, SCATTER-AUTOTUNE) so
c     changing them doesn't require a rebuild.
c
      q_max = 4.0d0
      q_spat = 0.1d0
      q_threshold = 10.0d0
      decimation_power = 2.0d0
      normalization_mode = 1
      manual_factor = 0.9d0
      scatter_autotune = .false.
      normmode_set = .false.
      expminusone_set = .false.
      scatter_normalize = .false.
      scatter_peakemph = .false.
      peak_qlo = 0.0d0
      peak_qhi = 0.0d0
      peak_factor = 1.0d0

      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:13) .eq. 'SCATTER-QMAX ') then
            read (string,*,err=1,end=1)  q_max
    1       continue
         else if (keyword(1:14) .eq. 'SCATTER-QSPAT ') then
            read (string,*,err=2,end=2)  q_spat
    2       continue
         else if (keyword(1:16) .eq. 'SCATTER-QTHRESH ') then
            read (string,*,err=3,end=3)  q_threshold
    3       continue
         else if (keyword(1:17) .eq. 'SCATTER-DECPOWER ') then
            read (string,*,err=4,end=4)  decimation_power
    4       continue
         else if (keyword(1:17) .eq. 'SCATTER-NORMMODE ') then
            read (string,*,err=5,end=5)  normalization_mode
            normmode_set = .true.
    5       continue
         else if (keyword(1:18) .eq. 'SCATTER-MANUALFAC ') then
            read (string,*,err=6,end=6)  manual_factor
    6       continue
         else if (keyword(1:20) .eq. 'SCATTER-EXPMINUSONE ') then
            read (string,*,err=7,end=7)  scatter_expminusone
            expminusone_set = .true.
    7       continue
         else if (keyword(1:17) .eq. 'SCATTER-AUTOTUNE ') then
            scatter_autotune = .true.
         else if (keyword(1:18) .eq. 'SCATTER-NORMALIZE ') then
            scatter_normalize = .true.
         else if (keyword(1:17) .eq. 'SCATTER-PEAKEMPH ') then
            read (string,*,err=8,end=8)  peak_qlo, peak_qhi,
     &                                   peak_factor
            scatter_peakemph = .true.
    8       continue
         end if
      end do

      q_min = 2 * pi / (xbox)
      q_step = q_min
      
      rdf_width = q_spat
      NN = int(q_max / rdf_width) + 1

      if (verbose_global) then

      print *, '---------------------------------------------'
      print *, ' RESOLUTION SETUP:'
      print *, ' q_spat (Bin Width): ', rdf_width
      print *, ' q_step (Lattice):   ', q_step
      print *, ' Total Bins (NN):    ', NN
      print *, '---------------------------------------------'

      end if

      allocate(unique_z(n))
      
      ewma_initialized = .false.
      current_md_step = 0
      read_file = .true.
      query = .true.
c
c     Get sliding window size from input or user
c
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  slw
  10    continue
         query = .false.
      end if
      if (query) then
         write (iout,20)
   20    format (/,' Enter dimension of the',
     &              ' sliding window:  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  slw
   40    continue
      end if 

      rdf_mean = slw
      if (rdf_mean .lt. 1) then
         write (iout,45)  rdf_mean
   45    format (/,' RADIALASK  --  Sliding Window Size of',i8,
     &              ' is Invalid, Value of 1 will be Used')
         rdf_mean = 1
      end if
c
c     Convert the requested window size into the smoothing factor
c     "alpha" of an exponentially-decaying memory function. The
c     restraint acts on a running exponential average of S(q),
c     <S>(t) = alpha*S_inst(t) + (1-alpha)*<S>(t-1) (see
c     structfactor_forces), rather than either the raw instantaneous
c     S(q) or a plain rectangular (box-car) moving average -- the
c     latter has a small force discontinuity every time the oldest
c     sample falls out of the window, since an EWMA has no fixed
c     window edge for a sample to fall out of.
c
c     alpha = 2/(span+1) is the standard "span" parameterization of an
c     exponential moving average (e.g. pandas' ewm(span=...)): it is
c     chosen so the EWMA has about the same memory depth as a simple
c     moving average over "span" points. A span of 1 gives alpha = 1,
c     i.e. <S>(t) = S_inst(t) exactly -- the undamped, purely
c     instantaneous restraint -- so SCATTER-RESTRAIN runs that used to
c     rely on window size 1 behave identically to before.
c
      scatter_alpha = 2.0d0 / (dble(rdf_mean) + 1.0d0)
c
c     Get force constant k value from input or user
c      
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=50,end=50)  kappa
  50    continue
         query = .false.
      end if
      if (query) then
         write (iout,60)
   60    format (/,' Enter force constant (K):  ',$)
         read (input,70)  record
   70    format (a240)
         read (record,*,err=80,end=80)  kappa
   80    continue
      end if
      rdf_kappa = kappa
c
c     Build a filename tag ("<base>_k<value>") so that runs with a
c     different restraint force constant k on the same input file
c     don't overwrite each other's diagnostic output, while keeping
c     Tinker's usual filename(1:leng) naming convention
c
      write (kstr,'(f0.4)')  rdf_kappa
      if (kstr(1:1) .eq. '.')  kstr = '0'//kstr
      call trim_trailing_zeros (kstr)
      scatter_tag = filename(1:leng)//'_k'//trim(kstr)
      scatter_basetag = filename(1:leng)
c
c     Count number of different occurrences of atom names
c
      num_unique = 0
      do i = 1, n
          is_new = .true.
          do j = 1, num_unique
              if (atomic(i) == unique_z(j)) then
                  is_new = .false.
                  exit
              end if
          end do
          if (is_new) then
             num_unique = num_unique + 1
             unique_z(num_unique) = atomic(i)
          end if
      end do
      
      n_species = num_unique
c
c     Allocate global atom_name array and fill it from atomic numbers
c
      allocate(atom_name(n_species))
      do i = 1, n_species
          atom_name(i) = get_name_from_atomic_number(unique_z(i))
      end do
c
c     Map atoms to type index AND count occurrences 
c
      allocate(atype(n))
      allocate(mole_fractions(n_species))
      mole_fractions = 0.0d0

      do i = 1, n
          do j = 1, n_species
              if (atomic(i) == unique_z(j)) then                  ! If the atom's Z-number matches a number in the unique list
                  atype(i) = j                                    ! store its type index
                  mole_fractions(j) = mole_fractions(j) + 1.0d0   ! and increment the count for that type
                  exit
              end if
          end do
      end do
c
c     Normalize counts to get mole fractions
c
      mole_fractions = mole_fractions / dble(n)
c
c     Get maximum distance from input or minimum image convention
c
      if (.not. use_bounds) then
         rmax = -1.0d0
        query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=90,end=90)  rmax
            query = .false.
         end if
   90    continue
         if (query) then
            write (iout,100)
  100       format (/,' Enter Maximum Distance to Accumulate',
     &                 ' [10.0 Ang] :  ',$)
            read (input,110)  rmax
  110       format (f20.0)
         end if
         if (rmax .le. 0.0d0)  rmax = 10.0d0
      else if (octahedron) then
         rmax = (sqrt(3.0d0)/4.0d0) * xbox
         rmax = 0.95d0 * rmax
      else
         rmax = min(xbox2*beta_sin*gamma_sin,ybox2*gamma_sin,
     &                         zbox2*beta_sin)
         rmax = 0.95d0 * rmax
      end if
c
c     Overwrite output files (named from the run's base filename plus
c     the restraint force constant k, so concurrent runs on different
c     systems, or different k values for the same system, in the same
c     directory don't clobber each other's diagnostics)
c
      unit = freeunit ()
      open(unit,file=trim(scatter_tag)//'_st.csv',status='replace')
      write(unit, 180)
  180 format('q_value(A^-1);S_q_avg')
      close(unit)

      unit = freeunit ()
      open(unit,file=trim(scatter_basetag)//'_st_init.csv',
     &     status='replace')
      write(unit, 190)
  190 format('q_value(A^-1);S_q_avg')
      close(unit)
c
c     S(q) time series (item 27).  _st.csv is rewritten in place at every
c     dump and so only ever holds the LAST snapshot; this companion file
c     is appended to instead, keeping every dump, so S(q,t) can be
c     reconstructed after the run.  Both the instantaneous per-step S(q)
c     and the EWMA (the quantity the restraint actually biases) are kept,
c     so the lag between them is visible.  Created (truncated) here once
c     at setup; structfactor_direct.f only ever appends.
c
      unit = freeunit ()
      open(unit,file=trim(scatter_tag)//'_st_traj.csv',
     &     status='replace')
      write(unit, 195)
  195 format('Step;q_value(A^-1);S_inst;S_ewma')
      close(unit)

      unit = freeunit()
      open(unit, file=trim(scatter_tag)//'_energy.csv',
     &     status='replace')
c     energy.csv header, temporarily reduced to vdW + restraint RMS(dS)
c     only (keep in sync with the write block in gradient.f):
c     write(unit, '(A)') 'Step,Ebond,Eangle,EUrey-Bradley,'//
c    &         'EvdW,Echarge-charge,Ex,Total'
      write(unit, '(A)') 'Step,EvdW,RMS_dS'
      close(unit)
c
c     Perform array allocations and initializations
c
      allocate (S_bin_ewma(NN))
      allocate (q_magnitude(NN))

      S_bin_ewma = 0.0d0

      do i = 1, NN
            q_magnitude(i) =  dble(i) * rdf_width
      end do
c
c     Per-bin restraint weights: 1.0 everywhere by default, so with
c     SCATTER-PEAKEMPH absent the weighted energy reduces exactly to the
c     unweighted one. When present, bins whose q falls in [peak_qlo,
c     peak_qhi] are scaled by peak_factor, concentrating the restraint
c     force on that q-window (e.g. the first-peak region).
c
      allocate (S_exp_weight(NN))
      S_exp_weight = 1.0d0
      if (scatter_peakemph) then
         do i = 1, NN
            if (q_magnitude(i) .ge. peak_qlo .and.
     &          q_magnitude(i) .le. peak_qhi) then
               S_exp_weight(i) = peak_factor
            end if
         end do
         if (verbose_global) then
            print *, '---------------------------------------------'
            print *, ' SCATTER-PEAKEMPH active:'
            print *, '   q window (A^-1): ', peak_qlo, peak_qhi
            print *, '   weight factor:   ', peak_factor
            print *, '---------------------------------------------'
         end if
      end if
c
c   Generate q-vector lattice
c
      call generate_lattice (q_max, q_min, q_step,
     & q_spat)
c
c   Apply decimation to q-vectors
c
      call decimation (q_threshold, decimation_power,q_max)
c
c   Load experimental data
c
      call load_experimental_data()
c
c   Optionally auto-select the best normalization mode (1 or 2) and
c   whether to subtract 1 from the experimental S(q), by comparing the
c   calculated S(q) at the starting configuration against experiment
c
      if (scatter_autotune) then
         if (normmode_set) then
            print *, 'WARNING -- SCATTER-NORMMODE was set in the',
     &         ' keyfile, but SCATTER-AUTOTUNE is also active;',
     &         ' the keyfile value will be overridden by the',
     &         ' auto-selected normalization mode.'
         end if
         if (expminusone_set) then
            print *, 'WARNING -- SCATTER-EXPMINUSONE was set in the',
     &         ' keyfile, but SCATTER-AUTOTUNE is also active;',
     &         ' the keyfile value will be overridden by the',
     &         ' auto-selected shift setting.'
         end if
         call autotune_scatter_settings (normalization_mode,
     &                                    manual_factor)
      end if
c
c   Precompute form factors
c
      call precompute_form_factors(normalization_mode,
     & manual_factor)

      end subroutine radialask

c     ################################################################
c     ##  Strip trailing zeros (and a bare trailing decimal point)  ##
c     ##  left over from an F-edit-descriptor write, e.g. so that   ##
c     ##  50.0000 -> 50 and 12.5000 -> 12.5                         ##
c     ################################################################
      subroutine trim_trailing_zeros (str)
      implicit none
      character(len=*), intent(inout) :: str
      integer :: i

      i = len_trim(str)
      do while (i .gt. 1 .and. str(i:i) .eq. '0')
         i = i - 1
      end do
      if (str(i:i) .eq. '.')  i = i - 1
      str = str(1:i)

      end subroutine trim_trailing_zeros

c     ################################################################
c     ##  Load and Re-bin Experimental Structure Factor             ##
c     ################################################################
      subroutine load_experimental_data()
      use rdfparams
      use iounit
      implicit none
      
      integer :: i, ios, bin_idx
      integer :: file_unit, debug_unit
      integer :: bin_counts(NN)
      real*8 :: q_val, s_val
      logical :: file_exists
      exp_present = .true.
c
c   Check if file exists
c
      inquire(file='water_sfac.dat', exist=file_exists)
      if (.not. file_exists) then
          write(*,*) 'EXPERIMENTAL DATA FILE NOT FOUND.'
          write(*,*) 'Performing base calculation', 
     &       'without experimental comparison.'
          exp_present = .false.
          return
      end if
c
c   Allocate arrays, read data, and bin it according to q values
c
      if (allocated(S_exp_binned)) deallocate(S_exp_binned)
      allocate(S_exp_binned(NN))

      S_exp_binned = 0.0d0
      bin_counts = 0

      open(newunit=file_unit, file='water_sfac.dat', status='old',
     &     action='read')
      do
          read(file_unit, *, iostat=ios) q_val, s_val
          if (ios /= 0) exit 

          bin_idx = nint(q_val / rdf_width)

          if (bin_idx .ge. 1 .and. bin_idx .le. NN) then
              S_exp_binned(bin_idx) = S_exp_binned(bin_idx) + s_val
              bin_counts(bin_idx) = bin_counts(bin_idx) + 1
          end if
      end do
      close(file_unit)

      do i = 1, NN
          if (bin_counts(i) .gt. 0) then
              S_exp_binned(i) = S_exp_binned(i) / dble(bin_counts(i))
          else
              S_exp_binned(i) = 0.0d0
              print *, 'Missing experimental S(q) bin ', i
          end if
      end do
c
c   If the flag is set, subtract 1 from all S_exp values
c   (SCATTER-EXPMINUSONE keyfile-configurable; see radialask())
c
      if (scatter_expminusone) then
          do i = 1, NN
              S_exp_binned(i) = S_exp_binned(i) - 1.0d0
          end do
      end if

      if (debug_global) then

        open(newunit=debug_unit,
     &     file=trim(scatter_basetag)//'_debug_experimental_sfac.csv',
     &     status='replace')
        
        write(debug_unit, '(a)') 'q_value,S_exp_rebinned,count'

        do i = 1, NN
            if (bin_counts(i) .gt. 0) then
                write(debug_unit, '(f12.6,a,f12.6,a,i8)') 
     &              q_magnitude(i), ',', 
     &              S_exp_binned(i), ',', 
     &              bin_counts(i)
            end if
        end do

        close(debug_unit)

      end if

      if (verbose_global) then
        print *, '---------------------------------------------'

        print *, ' LOADED EXPERIMENTAL S(Q):'
        print *, ' Total Bins Loaded: ', NN
        print *, ' Minus One Subtracted: ', scatter_expminusone
        print *, '---------------------------------------------'
      end if

      return
      end subroutine load_experimental_data

c     ################################################################
c     ##  Auto-select the best normalization mode (1 or 2, manual   ##
c     ##  mode 3 is never tried) and whether to subtract 1 from the ##
c     ##  experimental S(q), by comparing the S(q) calculated from  ##
c     ##  the starting configuration against experiment for all     ##
c     ##  four combinations and keeping the one with the lowest     ##
c     ##  mean squared deviation. Opt-in via SCATTER-AUTOTUNE.      ##
c     ################################################################
      subroutine autotune_scatter_settings (normalization_mode,
     &                                       manual_factor)
      use atomid
      use atoms
      use factors
      use rdfparams
      use iounit
      implicit none

      integer, intent(inout) :: normalization_mode
      real*8, intent(in) :: manual_factor

      integer :: mode, ishift, i, j, a, bin_idx
      integer :: best_mode, ncomp
      logical :: best_shift
      real*8 :: qx, qy, qz, q_dot_r
      real*8 :: cos_sum_w, sin_sum_w, S_raw, S_normalized
      real*8 :: sq_diff, best_diff
      real*8, dimension(n_species) :: C_species, S_species
      real*8, dimension(NN) :: S_bin_avg
      integer, dimension(NN) :: bin_counts

      if (.not. exp_present) then
         if (verbose_global) then
            print *, 'SCATTER-AUTOTUNE: no experimental data loaded,',
     &               ' skipping automatic mode/shift selection'
         end if
         return
      end if

      best_diff = huge(1.0d0)
      best_mode = normalization_mode
      best_shift = scatter_expminusone

c     Try both candidate normalization modes (never manual mode 3)
      do mode = 1, 2
         call precompute_form_factors(mode, manual_factor)

c        Calculate S(q) per bin from the starting configuration
c        (single unweighted pass, same math as structfactor_forces'
c        Pass 1, but self-contained here since it only runs a
c        handful of times at setup, never during the MD loop)
         S_bin_avg = 0.0d0
         bin_counts = 0
         do j = 1, max_vect
            qx = q_vec(j,1)
            qy = q_vec(j,2)
            qz = q_vec(j,3)
            bin_idx = int(q_vec(j,4))
            if (bin_idx < 1 .or. bin_idx > NN) cycle

            C_species = 0.0d0
            S_species = 0.0d0
            do i = 1, n
               q_dot_r = qx*x(i) + qy*y(i) + qz*z(i)
               C_species(atype(i)) = C_species(atype(i))
     &                                + cos(q_dot_r)
               S_species(atype(i)) = S_species(atype(i))
     &                                + sin(q_dot_r)
            end do

            cos_sum_w = 0.0d0
            sin_sum_w = 0.0d0
            do a = 1, n_species
               cos_sum_w = cos_sum_w + pc_fa(j,a) * C_species(a)
               sin_sum_w = sin_sum_w + pc_fa(j,a) * S_species(a)
            end do

            S_raw = (cos_sum_w**2 + sin_sum_w**2) / dble(n)
            if (pc_favg_sq(j) .ne. 0.0d0) then
               S_normalized = S_raw / pc_favg_sq(j) - pc_self_term(j)
            else
               S_normalized = 0.0d0
            end if

            S_bin_avg(bin_idx) = S_bin_avg(bin_idx) + S_normalized
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1
         end do

         do i = 1, NN
            if (bin_counts(i) .gt. 0) then
               S_bin_avg(i) = S_bin_avg(i) / dble(bin_counts(i))
            end if
         end do

c        Try both candidate values of the -1 shift against this mode
         do ishift = 0, 1
            scatter_expminusone = (ishift .eq. 1)
            call load_experimental_data ()

            sq_diff = 0.0d0
            ncomp = 0
            do i = 1, NN
               if (bin_counts(i) .gt. 0 .and.
     &             S_exp_binned(i) .ne. 0.0d0) then
                  sq_diff = sq_diff
     &                      + (S_bin_avg(i) - S_exp_binned(i))**2
                  ncomp = ncomp + 1
               end if
            end do

            if (ncomp .gt. 0) then
               sq_diff = sq_diff / dble(ncomp)
               if (verbose_global) then
                  print *, ' SCATTER-AUTOTUNE: mode', mode,
     &                     ' shift', scatter_expminusone,
     &                     ' MSE', sq_diff, ' over', ncomp, ' bins'
               end if
               if (sq_diff .lt. best_diff) then
                  best_diff = sq_diff
                  best_mode = mode
                  best_shift = scatter_expminusone
               end if
            end if
         end do
      end do

c     Lock in the winning combination (form factors get their final
c     recompute back in radialask() right after this routine returns)
      normalization_mode = best_mode
      scatter_expminusone = best_shift
      call load_experimental_data ()

      if (verbose_global) then
         print *, '---------------------------------------------'
         print *, ' SCATTER-AUTOTUNE RESULT:'
         print *, ' Normalization Mode: ', normalization_mode
         print *, ' Minus One Subtracted: ', scatter_expminusone
         print *, ' Mean Sq. Deviation:  ', best_diff
         print *, '---------------------------------------------'
      end if

      return
      end subroutine autotune_scatter_settings

c     ==================================================================
c     NORMALIZATION SETTINGS
c     1 = Average Atom  (C++ Case 0): Divides by <f>^2.
c     2 = Faber-Ziman   (C++ Case 1): Divides by <f^2>.
c     3 = Manual        (C++ Case 2): Divides by manual_factor.
c     4 = None          (C++ Case 3): Divides by 1.
c     ==================================================================