
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

         integer :: i, j, f, a, bin_idx
         integer :: istep, freeunit, unit
         real*8 :: qx, qy, qz, q_dot_r
         real*8, dimension(NN) :: S_q_current, S_q_mean
         integer, dimension(NN) :: bin_counts
         real*8, dimension(n_species) :: C_species, S_species
         
         real*8 :: cos_sum_w, sin_sum_w, S_raw, S_normalized 
         
         print *, 'step',istep

         if (ctn > rdf_mean) ctn = 1

         S_q_current = 0.0d0
         S_q_mean = 0.0d0
         bin_counts = 0
         f = 1
c
c        Main loop over q-vectors (max_vect) 
c
!$OMP PARALLEL DO DEFAULT(NONE) NUM_THREADS(12)
!$OMP& SHARED(max_vect, q_vec, n,pc_self_term)
!$OMP& SHARED(x, y, z, atype, pc_fa, pc_favg_sq, n_species) 
!$OMP& PRIVATE(j, qx, qy, qz, bin_idx, i, q_dot_r, a)  
!$OMP& PRIVATE(C_species, S_species, cos_sum_w, sin_sum_w)
!$OMP& PRIVATE(S_raw, S_normalized)
!$OMP& REDUCTION(+:S_q_current) REDUCTION(+:bin_counts)
         do j = 1, max_vect
            qx = q_vec(j, 1)
            qy = q_vec(j, 2)
            qz = q_vec(j, 3)
            bin_idx = int(q_vec(j, 4))
            if (bin_idx < 1) then
c               print *, 'Warning: bin_idx < 1, skipping q-vector', j
               cycle
            end if

c           Step 1: Accumulate cos/sin sums per species            
            C_species = 0.0d0
            S_species = 0.0d0

            do i = 1, n
               q_dot_r = qx*x(i) + qy*y(i) + qz*z(i)
               C_species(atype(i)) = C_species(atype(i)) + cos(q_dot_r)
               S_species(atype(i)) = S_species(atype(i)) + sin(q_dot_r)
            end do

c           Step 2: Apply form factors to get weighted sums
            cos_sum_w = 0.0d0
            sin_sum_w = 0.0d0

            do a = 1, n_species
                cos_sum_w = cos_sum_w + pc_fa(j, a) * C_species(a)
                sin_sum_w = sin_sum_w + pc_fa(j, a) * S_species(a)
            end do

c           Step 3: Compute RAW (unnormalized) structure factor
            S_raw = (cos_sum_w**2 + sin_sum_w**2) / dble(n)

c           Step 4: Apply normalization (from precomputed denominator)
            if (pc_favg_sq(j) == 0.0d0) then
                S_normalized = 0.0d0
                print *, 'warning: precomputed_favg_sq = 0', 
     &                   ' for q-vector', j
            else
                S_normalized = S_raw / pc_favg_sq(j)
                S_normalized = S_normalized - pc_self_term(j)
            end if

            S_q_current(bin_idx) = S_q_current(bin_idx) + S_normalized
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1

         end do
!$OMP END PARALLEL DO

         ! Normalize by bin count
         do i = 1, NN
            if (bin_counts(i) > 0) then
                  S_q_current(i) = S_q_current(i) / dble(bin_counts(i))
c                  print *, 'bin ', i, ' count ', bin_counts(i)
            else
                  S_q_current(i) = 0.0d0
            end if
         end do

         ! Store in history
         do i = 1, NN
            S_hist(i, ctn, f) = S_q_current(i)
         end do

         ! Sliding window logic
         if (ctn == rdf_mean) savelock = 1

         if (savelock == 1) then
            unit = freeunit()
            ! Ensure this matches the filename above
            open(unit, file='st.csv', status='unknown'
     &       , position='append')

            do i = 1, (NN)
                  S_q_mean(i) = 0.0d0
                  do j = 1, rdf_mean
                     S_q_mean(i) = S_q_mean(i) + S_hist(i, j, f)
                  end do
                  S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)

                  if (S_q_mean(i) == 0.0d0) cycle
                  
                  write(unit, 310) q_magnitude(i), S_q_mean(i)
                  
310               format(f10.4, ';', f10.4) 
            end do
            close(unit)
         end if

         ctn = ctn + 1
         return            
      return      
      
      end