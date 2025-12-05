      subroutine decimation(q_threshold, decimation_power, q_max)
      use rdfparams
      implicit none
      real*8, intent(in) :: q_threshold
      real*8, intent(in) :: decimation_power
      real*8, intent(inout) :: q_max
      
      integer :: i, j, n_kept
      real*8 :: q_mag, q_ratio
      integer :: decimation_factor
      logical :: keep_vector
      integer :: file_unit
      
c      write(*,*) 'Applying q-vector reduction:'
c      write(*,*) '  Threshold q =', q_threshold
c      write(*,*) '  Decimation power =', decimation_power
c      write(*,*) '  Initial vectors:', max_vect
      
      n_kept = 0
      j = 0
      
      do i = 1, max_vect
        q_mag = sqrt(q_vec(i,1)**2 + q_vec(i,2)**2 + q_vec(i,3)**2)
        if (q_mag .le. q_threshold) then
c           Keep all vectors below threshold
            keep_vector = .true.
        else if (q_mag .gt. q_max) then
c           Keep all vectors above q_max
            keep_vector = .true.
        else
c            print * , 'reducing'
c           Apply tunable decimation      
            q_ratio = q_mag / q_threshold
            decimation_factor = max(1, int(q_ratio**decimation_power))           
            j = j + 1
            keep_vector = (mod(j, decimation_factor) .eq. 0)
        end if
        
        if (keep_vector) then
            n_kept = n_kept + 1
            if (n_kept .ne. i) then
                q_vec(n_kept, 1:4) = q_vec(i, 1:4)
            end if
        end if
      end do
      
      write(*,*) '  Final vectors:', n_kept
      write(*,*) '  Reduction factor:', real(max_vect)/real(n_kept)
      
      max_vect = n_kept
      if (savelock .eq. 1000) then
    
c       Write to CSV file
        file_unit = 11
        open(unit=file_unit, file='qmatrix_red.csv', status='replace')
                  
c       Write header
        write(file_unit, '(a)') 'x,y,z,bin_number'
                  
c       Write data
        do i = 1, max_vect
            write(file_unit, '(f12.6,a,f12.6,a,f12.6,a,i4)') 
     &          q_vec(i,1), ',', 
     &          q_vec(i,2), ',', 
     &          q_vec(i,3), ',', 
     &          int(q_vec(i,4))
        end do

                
        close(file_unit)
      endif
      
      return
      end