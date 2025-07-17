c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine structfactor_direct  --  compute structure factor directly  ##
c     ##                        using atomic positions and scattering formula    ##
c     ##                                                                          ##
c     ##############################################################################
c

      subroutine structfactor_direct(istep)
      use atomid
      use atoms
      use rdfparams
      implicit none

      integer i, j, f, bin_idx
      integer istep
      integer freeunit
      integer unit
      integer block_size, full_blocks, remainder
      integer block, start_idx, end_idx
      integer bin_offset
      
      ! Arrays to store structure factor for sliding window average
      real*8 S_q_current(NN)     ! Current step structure factor
      real*8 S_q_mean(NN)        ! Mean structure factor over sliding window
      
      ! Working variables
      real*8 :: qx, qy, qz
      real*8 :: cos_sum, sin_sum, S_direct, q_dot_r
      real*8 :: q_magnitude

      integer :: bin_counts(NN)

      if (ctn > rdf_mean) then
            ctn = 1
      end if

      ! Initialize arrays
c      print *, 'Starting structfactor_direct_spherical...'
c      print *, 'Step:', istep   
      S_q_current = 0.0d0
      S_q_mean = 0.0d0
      f = 1

c      print *, max_vect, 'q-vectors to process'
      
      ! Initialize bin counts
!$OMP SIMD
      do i = 1, NN
          bin_counts(i) = 0
      end do
  
      f = 1
      
      ! Set block size for cache optimization
      block_size = 256
      full_blocks = n / block_size
      remainder = mod(n, block_size)
         
      do j = 1, max_vect
         qx = q_vec(j, 1)
         qy = q_vec(j, 2)
         qz = q_vec(j, 3)
         bin_idx = int(q_vec(j, 4))

         if (bin_idx < 1 .or. bin_idx > tot_bins) then
            print *, 'Warning: out of range:', bin_idx
            cycle
         end if
   
         ! Initialize sums for this q-vector
         cos_sum = 0.0d0
         sin_sum = 0.0d0
         
         ! Process atoms in blocks for better cache usage
         do block = 0, full_blocks - 1
            start_idx = block * block_size + 1
            end_idx = start_idx + block_size - 1
            
            ! Vectorized computation of dot products and trigonometric functions
!$OMP SIMD REDUCTION(+:cos_sum,sin_sum)
            do i = start_idx, end_idx
               q_dot_r = qx * x(i) + qy * y(i) + qz * z(i)
               cos_sum = cos_sum + cos(q_dot_r)
               sin_sum = sin_sum + sin(q_dot_r)
            end do
!$OMP END SIMD
         end do
         
         ! Process remaining atoms
         if (remainder > 0) then
            start_idx = full_blocks * block_size + 1
!$OMP SIMD REDUCTION(+:cos_sum,sin_sum)
            do i = start_idx, n
               q_dot_r = qx * x(i) + qy * y(i) + qz * z(i)
               cos_sum = cos_sum + cos(q_dot_r)
               sin_sum = sin_sum + sin(q_dot_r)
            end do
!$OMP END SIMD
         end if
               
         ! Calculate structure factor for this q-vector
         S_direct = (cos_sum**2 + sin_sum**2) / n
         bin_offset = bin_idx + 40
         S_q_current(bin_offset) = S_q_current(bin_offset) + S_direct
         bin_counts(bin_offset) = bin_counts(bin_offset) + 1
      end do    
      
      ! Divide by number of q-vectors in each bin - vectorized
!$OMP SIMD
      do i = 1, NN
         if (bin_counts(i) .gt. 0) then
            S_q_current(i) = S_q_current(i) / dble(bin_counts(i))
         else
            S_q_current(i) = 0.0d0
         end if
      end do      
       
      ! Store current step structure factor in sliding window array - vectorized
!$OMP SIMD
      do i = 1, NN
         S_hist(i, ctn, f) = S_q_current(i)
      end do                      
      
      ! Check if we need to calculate and save sliding window average
      if (ctn .eq. rdf_mean) then
         savelock = 1
      endif
      
      ! Calculate and write sliding window average for each q-value
      if (savelock .eq. 1) then
         print *, 'Calculating sliding window average...'
         print *, 'Step:', istep
         unit = freeunit ()
         open(unit, file='st.txt', status='unknown',  
     &        position='append')
         write(unit, 301) istep
 301     format(/,' Step: ',i8,/)
         
         do i = 1, NN
            ! Calculate sliding window average across all atom pairs
            S_q_mean(i) = 0.0d0
            do f = 1, rdf_num
               ! Calculate mean over the sliding window (last rdf_mean steps)
               do j = 1, rdf_mean
                  S_q_mean(i) = S_q_mean(i) + S_hist(i,j,f)
               end do
               S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)
            end do
            S_q_mean(i) = S_q_mean(i) / dble(rdf_num)  ! Average over atom pairs
            q_magnitude =  (i) * 0.01d0

            ! Write q-value and sliding window averaged structure factor
            if (S_q_mean(i) .eq. 0.0d0) then
               cycle
            end if
            write(unit, 310) q_magnitude, S_q_mean(i)
 310        format(3x, f12.4, 3x, f12.4)
         end do       
         close(unit)
      endif
      
      ! Increment counter (matching original radialsub logic)
      ctn = ctn + 1
      
      return      

      end

c     ##############################################################################
c     ##                                                                          ##
c     ##  Alternative version with spherical averaging over q-directions          ##
c     ##                                                                          ##
c     ##############################################################################
c
c     This version averages over multiple q-directions for more accurate results
c
      subroutine structfactor_direct_spherical(istep)
      use atomid
      use atoms
c      use bound
c      use boxes
c      use limits
c      use math
c      use molcul
c      use potent
      use rdfparams
c      use factors
      implicit none

      integer i, j, f, bin_idx
      integer istep
c      integer numj, numk
      integer freeunit
      integer unit
c      character*3 namej, namek
c      integer typej, typek
      
c     Arrays to store structure factor for sliding window average
      real*8 S_q_current(NN)     ! Current step structure factor
      real*8 S_q_mean(NN)        ! Mean structure factor over sliding window
      
      
c     Working variables
      real*8 :: qx, qy, qz
      real*8 :: cos_sum, sin_sum, S_direct, q_dot_r
      real*8 :: q_magnitude

      integer :: bin_counts(NN)

c     Initialize arrays

      print *, 'Starting structfactor_direct_spherical...'
      print *, 'Step:', istep   
      S_q_current = 0.0d0
      S_q_mean = 0.0d0

      f=1

      
c      namej = rdf_namej(f)
c      namek = rdf_namek(f)
c      typej = rdf_typej(f)
c      typek = rdf_typek(f)

cc        Count atoms of each type for normalization
c         numj = 0
c         numk = 0
c         do i = 1, n
            ! you can do this with count
c            if (name(i).eq.namej .or. type(i).eq.typej) then
c               numj = numj + 1
c            end if
c            if (name(i).eq.namek .or. type(i).eq.typek) then
c               numk = numk + 1
c            end if
c         end do      

      print *, max_vect, 'q-vectors to process'
c     Initialize arrays
      do i = 1, NN
          bin_counts(i) = 0
      end do
  
      f = 1
         
      do j = 1, max_vect
         qx = q_vec(j, 1)
         qy = q_vec(j, 2)
         qz = q_vec(j, 3)
         bin_idx = int(q_vec(j, 4))
c         print *, 'Processing q-vector:', qx, qy, qz, 'Bin:', bin_idx

         if (bin_idx < 1 .or. bin_idx > tot_bins) then
            print *, 'Warning: out of range:', bin_idx
            cycle
         end if
   
c         Initialize sums for this q-vector
          cos_sum = 0.0d0
          sin_sum = 0.0d0
               
c         Sum over all atoms
            do i = 1, n
c              Calculate q·r for atom i
               q_dot_r = qx * x(i) + qy * y(i) + qz * z(i)
c              print *, 'q_dot_r:', q_dot_r, 'for atom:', i
                  
c              Add contributions to cosine and sine sums
               cos_sum = cos_sum + cos(q_dot_r)
c               print *, 'cos(q_dot_r):', cos(q_dot_r), 'for atom:', i
               sin_sum = sin_sum + sin(q_dot_r)
c               print *, 'sin(q_dot_r):', sin(q_dot_r), 'for atom:', i
            end do
               
c           Calculate structure factor for this q-vector
            S_direct = (cos_sum**2 + sin_sum**2) 
     &          / n
c            print *, 'S_direct:', S_direct
c            print *, 'S_direct for q-vector:', qx, qy, qz, 
c     &      'is:', S_direct
            S_q_current(bin_idx + 40) = S_q_current(bin_idx + 40) 
     &      + S_direct
c           Increment count for this bin
            bin_counts(bin_idx +40) = bin_counts(bin_idx + 40) + 1
      end do 
c      print *, 'Finished processing q-vectors'
c      print*, 'bin counts:'
      do i = 1, NN
c         print *, 'Bin', i, ':', bin_counts(i)
      end do     
c     divide by number of q-vectors in each bin
      do i = 1, NN
         if (bin_counts(i) .gt. 0) then
            S_q_current(i) = S_q_current(i) / dble(bin_counts(i))
         else
            S_q_current(i) = 0.0d0
         end if
      end do      
c     store current step structure factor in sliding window array
      do i = 1, NN
         S_hist(i, ctn, f) = S_q_current(i)
c         print *, 'S_hist(', i, ',', ctn, ',', 
c     &      f, ') = ', S_hist(i, ctn, f)
      end do                      
      
c     Check if we need to calculate and save sliding window average
      if (ctn .eq. rdf_mean) then
         savelock = 1
      endif
      
c     Calculate and write sliding window average for each q-value
      if (savelock .eq. 1) then
c         print *, 'Calculating sliding window average...'
c         print *, 'Step:', istep
         unit = freeunit ()
         open(unit, file='st.txt', status='unknown',  
     &        position='append')
         write(unit, 301) istep
 301     format(/,' Step: ',i8,/)
         
         do i = 1, NN
c           Calculate sliding window average across all atom pairs
            S_q_mean(i) = 0.0d0
            do f = 1, rdf_num
c              Calculate mean over the sliding window (last rdf_mean steps)
               do j = 1, rdf_mean
                  S_q_mean(i) = S_q_mean(i) + S_hist(i,j,f)
               end do
               S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)
            end do
            S_q_mean(i) = S_q_mean(i) / dble(rdf_num)  ! Average over atom pairs
            q_magnitude =  (i) * 0.01d0

            
c           Write q-value and sliding window averaged structure factor
            if (S_q_mean(i) .eq. 0.0d0) then
               cycle
            end if
            write(unit, 310) q_magnitude, S_q_mean(i)
 310        format(3x, f12.4, 3x, f12.4)
         end do       
         close(unit)
      endif
      
c     Increment counter (matching original radialsub logic)
      ctn = ctn + 1
      
      return      

      end