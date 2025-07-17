      subroutine generate_qshell(q_max, q_mi, dq, q_max_u, dq_u)
c     Generate 3D q-vectors in spherical shell and assign bin numbers
c     Write results to CSV file with format: x, y, z, bin_number

      use rdfparams
      implicit none            
c     Input parameters
      real*8, intent(in) :: q_mi, q_max, dq
      real*8, intent(in) :: q_max_u, dq_u
                  
c     Local variables
      integer, parameter :: max_vectors = 5000000
      real*8 :: qx, qy, qz, q_mag
      real*8 :: bin_edges(200)
      real*8 :: bin_edges_u(1000)
      integer :: n_bins, bin_idx
      integer :: n_bins_u
      integer :: n_vectors
      integer :: i, j, k, idx, step
      integer :: n_grid, n_steps
      real*8 :: qval
                  
c     File handling
      integer :: file_unit
                  
c     Initialize counters
      n_vectors = 0
c     Allocate array for q-vectors
      allocate(q_vec(max_vectors, 4))  ! 4 columns: x, y, z, bin_number      
                  
c     Calculate number of bins
      n_bins = int((q_max - q_mi) / dq) + 1
c      print *, 'Number of bins:', n_bins

c     Calculate the other number of bins 
      n_bins_u = int((q_max_u - q_mi) / dq_u) + 1
c      print *, 'Number of bins up:', n_bins_u
                       
c     Create bin edges
      do i = 1, n_bins + 1
        bin_edges(i) = q_mi + real(i-1) * dq
      end do
      
c     other bin edges
      do i = 1, n_bins_u + 1
        bin_edges_u(i) = q_max + real(i-1) * dq_u
      end do      
                  
c     Calculate grid size
      n_grid = int(2.0d0 * q_max / dq) + 1
                  
c      write(*,*) 'Generating q-vectors...'
c      write(*,*) 'Grid size:', n_grid
c      write(*,*) 'q_min =', q_mi
c      write(*,*) 'q_max =', q_max
c      write(*,*) 'dq =', dq
c      write(*,*) 'Number of bins:', n_bins
                  
c     Generate q-vectors on Cartesian grid
      do i = 1, n_grid
        qx = -q_max + real(i-1) * dq
        do j = 1, n_grid
            qy = -q_max + real(j-1) * dq
            do k = 1, n_grid
                qz = -q_max + real(k-1) * dq
                              
c                Calculate magnitude
                 q_mag = sqrt(qx*qx + qy*qy + qz*qz)
                              
c                Check if vector is in spherical shell
                 if (q_mag .ge. q_mi .and. q_mag .le. q_max) then
                    n_vectors = n_vectors + 1
                                  
c                   Check array bounds
                    if (n_vectors .gt. max_vectors) then
                        write(*,*) 'Error: Too many vectors!'
                        write(*,*) 'Increase max_vectors parameter'
                             stop
                    end if
                                  
c                   Store coordinates
                    q_vec(n_vectors, 1) = qx
                    q_vec(n_vectors, 2) = qy
                    q_vec(n_vectors, 3) = qz
                                  
c                   Assign bin number
                    bin_idx = 0
                    do idx = 1, n_bins
                        if (q_mag .ge. bin_edges(idx) .and. 
     &                        q_mag .lt. bin_edges(idx+1)) then
                                bin_idx = idx*int(dq*100.0d0) 
c                                print *, 'bin_idx:', bin_idx    
                        exit
                        end if
                    end do
                                  
c                   Handle edge case for maximum q
                    if (q_mag .eq. q_max) then
                        bin_idx = n_bins*int(dq*100.0d0)
                    end if
                                  
                q_vec(n_vectors, 4) = bin_idx
c                print *, 'bin:', bin_idx
                end if
            end do
        end do
      end do
      
c     Add extended q-vectors along axes
      n_steps = int((q_max_u - q_max) / dq_u)
      do step = 1, n_steps
        qval = q_max + step * dq_u

       do i = 1, 3
          n_vectors = n_vectors + 1
          if (n_vectors .gt. max_vectors) then
              write(*,*) 'Error: Too many vectors in extended region!'
              stop
          end if

c         Set qx, qy, qz based on axis
          qx = 0.0d0
          qy = 0.0d0
          qz = 0.0d0
          if (i .eq. 1) qx = qval
          if (i .eq. 2) qy = qval
          if (i .eq. 3) qz = qval
          q_vec(n_vectors,1) = qx
          q_vec(n_vectors,2) = qy
          q_vec(n_vectors,3) = qz

          q_mag = sqrt(qx*qx + qy*qy + qz*qz)
          bin_idx = 0
          do idx = 1, n_bins_u
            if (q_mag .ge. bin_edges_u(idx) .and. 
     &          q_mag .lt. bin_edges_u(idx+1)) then
                bin_idx = (idx*int(dq_u*100.0d0) +
     &            n_bins*int(dq*100.0d0) -1) 
                exit
            end if
          end do
          if (q_mag .eq. q_max) bin_idx = n_bins_u
          q_vec(n_vectors,4) = real(bin_idx)
        end do
      end do


      write(*,*) 'Total q-vectors generated:', n_vectors

      max_vect = n_vectors
      if (savelock .eq. 1000) then                  
c     Write to CSV file
      file_unit = 10
      open(unit=file_unit, file='qmatrix.csv', status='replace')
                  
c     Write header
      write(file_unit, '(a)') 'x,y,z,bin_number'
                  
c     Write data
      do i = 1, n_vectors
        write(file_unit, '(f12.6,a,f12.6,a,f12.6,a,i4)') 
     &        q_vec(i,1), ',', 
     &        q_vec(i,2), ',', 
     &        q_vec(i,3), ',', 
     &        int(q_vec(i,4))
      end do
           
      close(file_unit)
      endif

      sfac_nbins = n_bins
      sfac_nbins_u = n_bins_u
      tot_bins = sfac_nbins + sfac_nbins_u
           
      return
      end 
            
        