c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong and Jay William Ponder  ##
c     ##                      Modified by FC                      ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine radial  --  compute radial distribution function on the fly  ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "radial" finds the radial distribution function for a specified
c     pair of atom types via analysis of a set of coordinate frames
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
      use atom_sort
      
      implicit none
      integer slw                                   ! sliding window size
      integer i, j
      real*8 rmax
      logical exist,query
      logical read_file
      integer freeunit
      integer unit
      character*240 record
      character*240 string

      real*8 q_min
      real*8, parameter :: q_max = 10.0d0          ! Maximum q value
      real*8, parameter :: q_spat = 0.1d0             ! Spatial bin size
      real*8, parameter :: q_aggr = 1.5d0             ! Overgeneration maximum q value
      real*8, parameter :: force = 1.0d0              ! Force of overgeneration
      real*8, parameter :: q_threshold = 2.0d0        ! Threshold for decimation
      real*8, parameter :: decimation_power = 2.0d0   ! Power for decimation
      logical, parameter :: sort = .false.

c     for the counter
      integer num_unique
      logical is_new 
      character :: unique(n)

c     global variables

      rdf_sort = sort      

      q_step = 0.1d0
      q_min = 4 * pi / (xbox)
      NN = int(q_max * (10.0d0 / q_step))
      ctn = 1
      savelock = 0
      rdf_num = 1
      read_file = .true.
      query = .true.
      if (query) then
         write (iout,20)
   20    format (/,' Enter dimension of the',
     &              ' sliding window:  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  slw
   40    continue
      end if

      print *, 'names' 
      do i = 1, n
          print *, trim(name(i))
      end do

      if (sort) then
            call prepare_atom_sort() 
            call apply_sorted_order()
      endif      
c             
c     count number of different occurrences of atom names
c
      num_unique = 0
      do i = 1, n
      is_new = .true.
            do j = 1, num_unique
                  if (trim(name(i)) == trim(unique(j))) then
                        is_new = .false.
                        exit
                  end if
            end do
            if (is_new) then
               num_unique = num_unique + 1
               unique(num_unique) = trim(name(i))
            end if
      end do
      
      n_species = num_unique
c
c     Map atoms to type index AND count occurrences 
c
      allocate(atype(n))
      allocate(mole_fractions(n_species))
      mole_fractions = 0.0d0

      do i = 1, n
          do j = 1, n_species
              ! If the name from the original list matches a name in our unique key...
              if (trim(name(i)) == trim(unique(j))) then
                  atype(i) = j  ! ...store its type index
                  
                  ! ...and increment the count for that type
                  mole_fractions(j) = mole_fractions(j) + 1.0d0 
                  
                  exit ! Move to the next atom (i)
              end if
          end do
      end do
      if (sort) then
            call apply_original_order()
      endif
c
c     Normalize counts to get mole fractions ---
c
      mole_fractions = mole_fractions / dble(n)

c
c     get maximum distance from input or minimum image convention
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

      rdf_mean = slw
c
c     overwrite rdf.txt file
c
      unit = freeunit ()
      open(unit,file='st.txt',status='unknown')
      write(unit, 180) 
  180 format(/,' Direct Structure Factor - Sliding Window Average',
     &          /,' q-value (Å⁻¹)    S(q) (sliding avg)',/)
      close(unit)      

      allocate (S_hist(NN,rdf_mean,rdf_num))
      allocate (q_magnitude(NN))
      
      S_hist = 0.0d0

      ! Only compute q_magnitude once, outside this routine if possible
      do i = 1, NN
            q_magnitude(i) =  (i) * q_step / 10.0d0
      end do

      call generate_qshell_imp (q_max, q_min, q_step,
     & q_spat, q_aggr, force)
      call decimation (q_threshold, decimation_power,q_max)

      call precompute_form_factors()

      end
