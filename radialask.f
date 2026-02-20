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
      use limits
      use math
      use molcul
      use potent
      use rdfparams
      use factors
      
      implicit none
      integer slw                                     ! Sliding window size
      integer i, j
      real*8 rmax
      logical exist,query
      logical read_file
      integer freeunit
      integer unit
      character*240 record
      character*240 string
      real*8 q_min
      real*8 kappa

      real*8, parameter :: q_max = 3.0d0              ! Maximum q value
      real*8, parameter :: q_spat = 0.1               ! Spatial bin size change also q step
      real*8, parameter :: q_threshold = 10.d0        ! Threshold for decimation
      real*8, parameter :: decimation_power = 2.0d0   ! Power for decimation 

      logical, parameter :: debug_ra = .true.         ! Debug flag writing additional files
      logical, parameter :: verbose_ra = .true.       ! Verbose flag for printing info to console

      integer, parameter :: normalization_mode = 1    ! Normalization mode for S(q) calculation refers 
      real*8, parameter  :: manual_factor = 0.9d0     ! to table on the bottom. 

      integer :: num_unique
      logical :: is_new
      integer, allocatable :: unique_z(:)
c
c     Allocate global variables and set parameters 
c
      verbose_global = verbose_ra
      debug_global = debug_ra

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
      
      ctn = 1
      savelock = 0
      current_md_step = 0
      rdf_num = 1
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
c     Overwrite output files
c
      unit = freeunit ()
      open(unit,file='st.csv',status='replace')
      write(unit, 180) 
  180 format('q_value(A^-1);S_q_avg')
      close(unit)

      unit = freeunit ()
      open(unit,file='st_init.csv',status='replace')
      write(unit, 190) 
  190 format('q_value(A^-1);S_q_avg')
      close(unit)

      unit = freeunit()
      open(unit, file='energy.csv', status='replace')
      write(unit, '(A)') 'Step,Ebond,Eangle,EUrey-Bradley,'//
     &         'EvdW,Echarge-charge,Ex,Total'
      close(unit)
c
c     Perform array allocations and initializations
c
      allocate (S_hist(NN,rdf_mean,rdf_num))
      allocate (q_magnitude(NN))
      
      S_hist = 0.0d0

      do i = 1, NN
            q_magnitude(i) =  dble(i) * rdf_width
      end do
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
c   Precompute form factors
c
      call precompute_form_factors(normalization_mode, 
     & manual_factor)

      end subroutine radialask

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
      logical, parameter :: expminusone = .true.            ! Flag to subtract 1 from experimental S(q) values
      exp_present = .true.
c
c   Check if file exists
c
      inquire(file='water_sfact.dat', exist=file_exists)
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

      open(newunit=file_unit, file='water_sface.dat', status='old', 
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
c
      if (expminusone) then
          do i = 1, NN
              S_exp_binned(i) = S_exp_binned(i) - 1.0d0
          end do
      end if

      if (debug_global) then

        open(newunit=debug_unit, file='debug_experimental_sfac.csv', 
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
        print *, ' Minus One Subtracted: ', expminusone
        print *, '---------------------------------------------'
      end if

      return
      end subroutine load_experimental_data

c     ==================================================================
c     NORMALIZATION SETTINGS 
c     1 = Average Atom  (C++ Case 0): Divides by <f>^2. 
c     2 = Faber-Ziman   (C++ Case 1): Divides by <f^2>. 
c     3 = Manual        (C++ Case 2): Divides by manual_factor.
c     4 = None          (C++ Case 3): Divides by 1.
c     ==================================================================