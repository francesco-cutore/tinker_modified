      subroutine generate_lattice(q_max, q_mi, dq, q_spat)

      use rdfparams
      implicit none

c     Input parameters
      real*8, intent(in) :: q_mi      ! Min q magnitude
      real*8, intent(in) :: q_max     ! Max q magnitude
      real*8, intent(in) :: dq        ! Fundamental frequency (2pi/L)
      real*8, intent(in) :: q_spat    ! Spatial bin width

c     Local variables
      real*8 :: qx, qy, qz, q_mag
      integer :: bin_idx
      integer :: n_vectors                ! Total count for allocation
      integer :: idx_vec                  ! Index for filling
      integer :: nx, ny, nz, n_limit
      integer :: i

c     File handling
      integer :: file_unit

c     Initialize
      n_vectors = 0                       ! Counter for number of q-vectors

c     Calculate limits for integer loops
      n_limit = int(q_max / dq) + 1

c     Number of bins matches the bin_idx clamp used below (absolute,
c     from zero -- consistent with q_magnitude(i) = i*rdf_width in
c     radialask.f), not offset by q_mi
      sfac_nbins = int(q_max / q_spat) + 1

c     Start Pass 1:
c     Count number of q-vectors in the shell

      do nx = -n_limit, n_limit
        do ny = -n_limit, n_limit
            do nz = -n_limit, n_limit
                
                qx = dble(nx) * dq
                qy = dble(ny) * dq
                qz = dble(nz) * dq
                q_mag = sqrt(qx*qx + qy*qy + qz*qz)

                if (q_mag .ge. q_mi .and. q_mag .le. q_max) then
                    n_vectors = n_vectors + 1
                end if

            end do
        end do
      end do

c     Start Pass 2:
c     Allocate the q_vec array

      if (allocated(q_vec)) deallocate(q_vec)
      allocate(q_vec(n_vectors, 4))

c     Start Pass 3:
c     Store the q-vectors and their bin numbers
      
      idx_vec = 0

      do nx = -n_limit, n_limit
        do ny = -n_limit, n_limit
            do nz = -n_limit, n_limit
                
                qx = dble(nx) * dq
                qy = dble(ny) * dq
                qz = dble(nz) * dq
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
      
      max_vect = idx_vec

      if (verbose_global) then
      print *, '---------------------------------------------'
      print *, ' GENERATED Q-VECTORS:'
      print *, ' Total Vectors: ', n_vectors
      print *, '---------------------------------------------'
      end if

      if (debug_global) then   

        file_unit = 10
        open(unit=file_unit,
     &        file=trim(scatter_basetag)//'_qmatrix_imp.csv',
     &        status='replace')
        write(file_unit, '(a)') 'x,y,z,bin_number'
        do i = 1, n_vectors
            write(file_unit, '(f12.6,a,f12.6,a,f12.6,a,i8)') 
     &            q_vec(i,1), ',', 
     &            q_vec(i,2), ',', 
     &            q_vec(i,3), ',', 
     &            int(q_vec(i,4))
        end do
        close(file_unit)
        
      endif

      return
      end subroutine generate_lattice
