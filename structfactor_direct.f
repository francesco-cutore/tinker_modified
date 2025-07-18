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

         if (bin_idx < 1 .or. bin_idx > NN) then
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
      use factors
      use rdfparams
      implicit none

      integer i, j, f, bin_idx
      integer istep
      integer freeunit
      integer unit
      integer bin_offset
      
      ! Arrays to store structure factor for sliding window average
      real*8 S_q_current(NN)     ! Current step structure factor
      real*8 S_q_mean(NN)        ! Mean structure factor over sliding window
      
      ! Working variables
      real*8 :: qx, qy, qz
      real*8 :: cos_sum, sin_sum, S_direct, q_dot_r
      real*8 :: q_magnitude(NN)
      integer :: bin_counts(NN)

      ! Variables for weighting the structure factor
      real*8 :: factor(2)
      real*8 :: sumquad

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
      do i = 1, NN
          bin_counts(i) = 0
      end do

      do i = 1, NN
         q_magnitude(i) = 0
      end do
  
      f = 1     
      do j = 1, max_vect
         qx = q_vec(j, 1)
         qy = q_vec(j, 2)
         qz = q_vec(j, 3)
         bin_idx = int(q_vec(j, 4))

c         if (bin_idx < 1 .or. bin_idx > NN) then
c            print *, 'Warning: out of range:', bin_idx
c            cycle
c         end if
   
         ! Initialize sums for this q-vector
         cos_sum = 0.0d0
         sin_sum = 0.0d0
            do i = 1, n
               q_dot_r = qx * x(i) + qy * y(i) + qz * z(i)
               cos_sum = cos_sum + cos(q_dot_r)
               sin_sum = sin_sum + sin(q_dot_r)
            end do
         
               
         ! Calculate structure factor for this q-vector
         S_direct = (cos_sum**2 + sin_sum**2) / n
         bin_offset = bin_idx
         S_q_current(bin_offset) = S_q_current(bin_offset) + S_direct
         bin_counts(bin_offset) = bin_counts(bin_offset) + 1
      end do  
      
c      print *, 'bin_counts:', bin_counts(:)
      
      ! Divide by number of q-vectors in each bin - vectorized
      do i = 1, NN
         if (bin_counts(i) .gt. 0) then
            S_q_current(i) = S_q_current(i) / dble(bin_counts(i))
         else
c            print *, 'Warning: no vectors in bin', i
            S_q_current(i) = 0.0d0
         end if
      end do      
c      do i = 1, NN
c         print *, 'S_q_current:', S_q_current(i)
c      end do

      ! Store current step structure factor in sliding window array - vectorized

      do i = 1, NN
         if (bin_counts(i) .gt. 0) then
            S_hist(i, ctn, f) = S_q_current(i)
         else
            S_hist(i, ctn, f) = 0.0d0
         end if
         S_hist(i, ctn, f) = S_q_current(i)
      end do
      
      ! Check if we need to calculate and save sliding window average
      if (ctn .eq. rdf_mean) then
         savelock = 1
      endif
      
      ! Calculate and write sliding window average for each q-value
      if (savelock .eq. 1) then
c         print *, 'Calculating sliding window average...'
c         print *, 'Step:', istep
         unit = freeunit ()
         open(unit, file='st.txt', status='unknown',  
     &        position='append')
         write(unit, 301) istep
 301     format(/,' Step: ',i8,/)
         
         do i = 1, NN
            ! Calculate sliding window average across all atom pairs
            S_q_mean(i) = 0.0d0

            ! Calculate mean over the sliding window (last rdf_mean steps)
            do j = 1, rdf_mean
                  S_q_mean(i) = S_q_mean(i) + S_hist(i,j,f)
            end do
            S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)

            S_q_mean(i) = S_q_mean(i) / dble(rdf_num)  ! Average over atom pairs
            
            q_magnitude(i) =  (i) * 0.01d0

c            print *, 'sq_mean:', S_q_mean(i)
c            print *, 'q_magnitude:', q_magnitude(i)


            ! Add atomic form factors to structure factor

c            factor(1) = atomF('H ', q_magnitude(i))
c            print *, 'H factor:', factor(1)
c            factor(2) = atomF('O ', q_magnitude(i))
c            print *, 'O factor:', factor(2)

            ! Calculate sum of squares for normalization

            sumquad = (2.0d0 / 3.0d0) * factor(1)**2 +
     &                (1.0d0 / 3.0d0) * factor(2)**2
c           S_q_mean(i) = S_q_mean(i) * d / sumquad
c            S_q_mean(i) = S_q_mean(i) * (factor(1)*factor(2)) / sumquad


            ! Write q-value and sliding window averaged structure factor
            if (S_q_mean(i) .gt. 3.0d0) then
               S_q_mean(i) = 3.0d0
            else if (S_q_mean(i) .eq. 0.0d0) then
               cycle
            end if
            write(unit, 310) q_magnitude(i), S_q_mean(i)
 310        format( 3x, f12.4, 3x, f12.4)
         end do       
         close(unit)
      endif
      
      ! Increment counter (matching original radialsub logic)
      ctn = ctn + 1
      
      return      
      
      end

