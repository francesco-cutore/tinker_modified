      subroutine generate_qshell(q_max, q_mi, dq, q_spat)
c     Generate 3D q-vectors in spherical shell and assign bin numbers
c     Write results to CSV file with format: x, y, z, bin_number

      use rdfparams
      implicit none 

c     Input parameters
      real*8, intent(in) :: q_mi, q_max, dq
                  
c     Local variables
      integer, parameter :: max_vectors = 5000000
      real*8 :: qx, qy, qz, q_mag
      real*8, allocatable :: bin_edges(:)
      real*8 :: q_spat
      real*8 :: mean
      integer :: n_bins, bin_idx
      integer :: n_sbins
      integer :: n_vectors
      integer :: i, j, k, idx
      integer :: n_grid

                  
c     File handling
      integer :: file_unit
                  
c     Initialize counters
      n_vectors = 0

c     Allocate array for q-vectors
      allocate(q_vec(max_vectors, 4))  ! 4 columns: x, y, z, bin_number      
                  
c     Calculate number of bins
      n_bins = int((q_max - q_mi) / dq) + 1

c      print *, 'Number of bins:', n_bins

c     calculate number of edges

      n_sbins = int((q_max - q_mi) / q_spat) + 1

c      print *, 'Number of sbins:', n_sbins

      sfac_nbins = n_sbins - 1

      allocate(bin_edges(n_sbins))
            
c     Create bin edges
      do i = 1, n_sbins
        bin_edges(i) = q_mi + real(i-1) * q_spat
      end do 

      do i = 1, n_sbins
        print *, 'bin_edges:', bin_edges(i)
      end do

      
c      print *, size(bin_edges)
                  
c     Calculate grid size
      n_grid = int(2.0d0 * q_max / dq) + 1
                  
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
                    do idx = 1, n_sbins -1
                        if (q_mag .ge. bin_edges(idx) .and. 
     &                        q_mag .lt. bin_edges(idx+1)) then
                                mean = (bin_edges(idx) + 
     &                           bin_edges(idx+1)) / 2.0d0
                                bin_idx = int(mean * 100)
                                exit   
                        end if
                    end do
c                    print *, 'bin_edges:', bin_edges(n_sbins)             
c                   Handle edge case for maximum q
                    if (q_mag .eq. q_max .or. 
     &               q_mag .gt. bin_edges(n_sbins)) then
                        bin_idx = int(q_max * 100)
                    end if
                                  
                q_vec(n_vectors, 4) = bin_idx
c                print *, 'bin:', bin_idx
                end if
            end do
        end do
      end do

      write(*,*) 'Total q-vectors generated:', n_vectors

      max_vect = n_vectors
c      if (savelock .eq. 1000) then                  
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
c      endif

c      sfac_nbins = n_bins
c      sfac_nbins_u = n_bins_u
c      tot_bins = n_bins

      deallocate(bin_edges)    
      
      return
      end 
            
        