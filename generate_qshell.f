      subroutine generate_lattice(q_max, q_mi, dq, q_spat)

      use rdfparams
      implicit none 

c     Input parameters
      real*8, intent(in) :: q_mi      ! Min q magnitude
      real*8, intent(in) :: q_max     ! Max q magnitude
      real*8, intent(in) :: dq        ! Fundamental frequency (2pi/L)
      real*8, intent(in) :: q_spat    ! Bin width

c     Local variables
      real*8 :: qx, qy, qz, q_mag
      real*8, allocatable :: bin_edges(:)
      integer :: bin_idx
      integer :: qtint
      integer :: n_sbins
      integer :: n_vectors      ! Total count for allocation
      integer :: idx_vec        ! Index for filling
      integer :: nx, ny, nz, n_limit
      integer :: i
      real*8 :: dq_l, q_spat_l
        
c     File handling
      integer :: file_unit

      dq_l = dq 
      q_spat_l = q_spat 
      
c     Initialize
      n_vectors = 0
      qtint = int(10.0d0 / dq)

c     Calculate limits for integer loops
      n_limit = int(q_max / dq_l) + 1

c     Setup Binning
      n_sbins = int((q_max - q_mi) / q_spat_l) + 1
      sfac_nbins = n_sbins - 1
      allocate(bin_edges(n_sbins))

      do i = 1, n_sbins
        bin_edges(i) = q_mi + real(i-1) * q_spat_l
      end do 

c     ###########################################
c     ##  PASS 1: COUNT VALID LATTICE VECTORS
c     ###########################################
      write(*,*) 'Pass 1: Counting lattice vectors...'
      
      do nx = -n_limit, n_limit
        do ny = -n_limit, n_limit
            do nz = -n_limit, n_limit
                
                qx = dble(nx) * dq_l
                qy = dble(ny) * dq_l
                qz = dble(nz) * dq_l
                q_mag = sqrt(qx*qx + qy*qy + qz*qz)

                if (q_mag .ge. q_mi .and. q_mag .le. q_max) then
                    n_vectors = n_vectors + 1
                end if

            end do
        end do
      end do

      write(*,*) 'Pass 1 complete. Vectors found:', n_vectors

c     ###########################################
c     ##  ALLOCATE Q_VEC
c     ###########################################
      if (n_vectors == 0) then
          write(*,*) 'Error: No q-vectors found in shell.'
          stop
      end if

      if (allocated(q_vec)) deallocate(q_vec)
      allocate(q_vec(n_vectors, 4))
      write(*,*) 'Allocated q_vec memory.'

c     ###########################################
c     ##  PASS 2: FILL Q_VEC ARRAY
c     ###########################################
      write(*,*) 'Pass 2: Filling array...'
      
      idx_vec = 0

      do nx = -n_limit, n_limit
        do ny = -n_limit, n_limit
            do nz = -n_limit, n_limit
                
                qx = dble(nx) * dq_l
                qy = dble(ny) * dq_l
                qz = dble(nz) * dq_l
                q_mag = sqrt(qx*qx + qy*qy + qz*qz)

                if (q_mag .ge. q_mi .and. q_mag .le. q_max) then
                    
                    idx_vec = idx_vec + 1
                    
                    if (idx_vec .gt. n_vectors) then 
                        write(*,*) 'Error: Pass 2 overflow'
                        stop
                    endif

                    ! Store coordinates
                    q_vec(idx_vec, 1) = qx
                    q_vec(idx_vec, 2) = qy
                    q_vec(idx_vec, 3) = qz

                    bin_idx = nint(q_mag / q_spat)
                    
                    ! Clamp to ensure we don't crash
                    if (bin_idx < 1) bin_idx = 1
                    if (bin_idx > int(q_max/q_spat) + 1) 
     &     bin_idx = int(q_max/q_spat) + 1
                    
                    q_vec(idx_vec, 4) = dble(bin_idx)

                end if
            end do
        end do
      end do

c     ###########################################
c     ##  FINALIZE
c     ###########################################
      
      write(*,*) 'Pass 2 complete.'
      max_vect = idx_vec

      ! Optional CSV dump
c      if (savelock .eq. 1000) then                  
        file_unit = 10
        open(unit=file_unit, file='qmatrix_imp.csv', status='replace')
        write(file_unit, '(a)') 'x,y,z,bin_number'
        do i = 1, n_vectors
            ! Changed i4 to i8 to prevent **** for large bin numbers
            write(file_unit, '(f12.6,a,f12.6,a,f12.6,a,i8)') 
     &            q_vec(i,1), ',', 
     &            q_vec(i,2), ',', 
     &            q_vec(i,3), ',', 
     &            int(q_vec(i,4))
        end do
        close(file_unit)
c      endif

      deallocate(bin_edges)    
      
      return
      end subroutine generate_lattice

c ################################################################
c ###      Improved qshellgeneration                           ###
c ################################################################
                                                              
      subroutine generate_qshell_imp(q_max, q_mi, dq, q_spat, 
     &      q_aggr, force)
c     Generate 3D q-vectors in spherical shell and assign bin numbers
c     Write results to CSV file with format: x, y, z, bin_number
c
c     MODIFIED: Uses a two-pass method to allocate q_vec to its
c               exact size, saving significant memory.
c
      use rdfparams
      implicit none

c     Input parameters
      real*8, intent(in) :: q_mi, q_max, dq, q_spat, q_aggr, force

c     Local variables
      real*8 :: qx, qy, qz, q_mag
      real*8, allocatable :: bin_edges(:)
      real*8 :: mean
      integer :: bin_idx
      integer :: qtint
      integer :: n_sbins
      integer :: n_vectors  ! Total number of vectors (from Pass 1)
      integer :: idx_vec    ! Index counter (for Pass 2)
      integer :: i, j, k, idx
      integer :: n_grid
      integer :: q_count

c     File handling
      integer :: file_unit, freeunit

c     Initialize counters
      n_vectors = 0
      q_count = 0
      qtint = int(10 / dq)

c     Calculate number of bins (for q_spat)
      n_sbins = int((q_max - q_mi) / q_spat) + 1
      sfac_nbins = n_sbins - 1
      allocate(bin_edges(n_sbins))

c     Create bin edges
      do i = 1, n_sbins
        bin_edges(i) = q_mi + real(i-1) * q_spat
      end do

c     ###########################################
c     ##  PASS 1: COUNT Q-VECTORS
c     ###########################################
      write(*,*) 'Pass 1: Counting q-vectors...'

c     --- Loop 1: Normal generation (q_aggr to q_max) ---
      n_grid = int(2.0d0 * q_max / dq) + 1
      do i = 1, n_grid
        qx = -q_max + real(i-1) * dq
        do j = 1, n_grid
            qy = -q_max + real(j-1) * dq
            do k = 1, n_grid
                 qz = -q_max + real(k-1) * dq
                 q_mag = sqrt(qx*qx + qy*qy + qz*qz)

c                Check if vector is in spherical shell
                 if (q_mag .ge. q_aggr .and. q_mag .le. q_max) then
                    n_vectors = n_vectors + 1
                 end if
            end do
        end do
      end do

c     --- Loop 2: Aggressive generation (q_mi to q_aggr) ---
      n_grid = int(2.0d0 * q_aggr / (dq/force)) + 1
      do i = 1, n_grid
        qx = -q_aggr + real(i-1) * dq/force
        do j = 1, n_grid
            qy = -q_aggr + real(j-1) * dq/force
            do k = 1, n_grid
                qz = -q_aggr + real(k-1) * dq/force
                 q_mag = sqrt(qx*qx + qy*qy + qz*qz)

c                Check if vector is in spherical shell
                 if (q_mag .ge. q_mi .and. q_mag .lt. q_aggr) then
                    n_vectors = n_vectors + 1
                    q_count = q_count + 1
                 end if
            end do
        end do
      end do

      write(*,*) 'Pass 1 complete.'
      write(*,*) 'Total q-vectors to generate:', n_vectors
      write(*,*) 'Aggressive q-vectors:', q_count

c     ###########################################
c     ##  ALLOCATE Q_VEC TO EXACT SIZE
c     ###########################################

      if (n_vectors == 0) then
          write(*,*) 'Error: No q-vectors were found. Stopping.'
          stop
      end if


      allocate(q_vec(n_vectors, 4))
      write(*,*) 'Allocated q_vec with size:', n_vectors

c     ###########################################
c     ##  PASS 2: FILL Q_VEC ARRAY
c     ###########################################
      write(*,*) 'Pass 2: Filling q-vector array...'

      idx_vec = 0 ! Use a new counter for the array index
      q_count = 0 ! Reset for the final print

c     --- Loop 1: Normal generation (q_aggr to q_max) ---
      n_grid = int(2.0d0 * q_max / dq) + 1
      do i = 1, n_grid
        qx = -q_max + real(i-1) * dq
        do j = 1, n_grid
            qy = -q_max + real(j-1) * dq
            do k = 1, n_grid
                 qz = -q_max + real(k-1) * dq
                 q_mag = sqrt(qx*qx + qy*qy + qz*qz)

                 if (q_mag .ge. q_aggr .and. q_mag .le. q_max) then
                    idx_vec = idx_vec + 1

c                   Check bounds (safety check, should not fail)
                    if (idx_vec > n_vectors) then
                        write(*,*) 'Error: Mismatch in vector count!'
                        stop
                    end if

c                   Store coordinates
                    q_vec(idx_vec, 1) = qx
                    q_vec(idx_vec, 2) = qy
                    q_vec(idx_vec, 3) = qz

c                   Assign bin number
                    bin_idx = 0
                    do idx = 1, n_sbins - 1
                       if (q_mag .ge. bin_edges(idx) .and.
     &                     q_mag .lt. bin_edges(idx+1)) then
                           mean = (bin_edges(idx) +
     &                             bin_edges(idx+1)) / 2.0d0
                           bin_idx = int(mean * qtint)
                           exit
                       end if
                    end do
                    q_vec(idx_vec, 4) = bin_idx
                 end if
            end do
        end do
      end do

c     --- Loop 2: Aggressive generation (q_mi to q_aggr) ---
      n_grid = int(2.0d0 * q_aggr / (dq/force)) + 1
      do i = 1, n_grid
        qx = -q_aggr + real(i-1) * dq/force
        do j = 1, n_grid
            qy = -q_aggr + real(j-1) * dq/force
            do k = 1, n_grid
                qz = -q_aggr + real(k-1) * dq/force
                 q_mag = sqrt(qx*qx + qy*qy + qz*qz)

                 if (q_mag .ge. q_mi .and. q_mag .lt. q_aggr) then
                    idx_vec = idx_vec + 1
                    q_count = q_count + 1

                    if (idx_vec > n_vectors) then
                        write(*,*) 'Error: Mismatch in vector count!'
                        stop
                    end if

c                   Store coordinates
                    q_vec(idx_vec, 1) = qx
                    q_vec(idx_vec, 2) = qy
                    q_vec(idx_vec, 3) = qz

c                   Assign bin number
                    bin_idx = 0
                    do idx = 1, n_sbins - 1
                        if (q_mag .ge. bin_edges(idx) .and.
     &                      q_mag .lt. bin_edges(idx+1)) then
                            mean = (bin_edges(idx) +
     &                              bin_edges(idx+1)) / 2.0d0
                            bin_idx = int(mean * qtint)
                            exit
                        end if
                    end do

c                   Handle edge case for exact q_aggr
                    if (q_mag .eq. q_aggr) then
                        bin_idx = int(q_aggr * qtint)
                    end if

                    q_vec(idx_vec, 4) = bin_idx
                 end if
            end do
        end do
      end do

c     ###########################################
c     ##  FINALIZE AND WRITE FILE
c     ###########################################

      write(*,*) 'Pass 2 complete.'
      write(*,*) 'Total q-vectors generated:', idx_vec
      print *, 'Aggressive q-vectors generated:', q_count

      max_vect = idx_vec ! Set the global variable

      if (savelock .eq. 1000) then
         file_unit = freeunit()
         open(unit=file_unit, file='qmatrix.csv', status='replace')

c        Write header
         write(file_unit, '(a)') 'x,y,z,bin_number'

c        Write data
         do i = 1, n_vectors ! n_vectors is the total size
           write(file_unit, '(f12.6,a,f12.6,a,f12.6,a,i4)')
     &           q_vec(i,1), ',',
     &           q_vec(i,2), ',',
     &           q_vec(i,3), ',',
     &           int(q_vec(i,4))
         end do

         close(file_unit)
      endif

      deallocate(bin_edges)

      return
      end       