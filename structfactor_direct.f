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
         use factors
         use rdfparams
         implicit none

         integer :: i, j, f, bin_idx
         integer :: istep, freeunit, unit
         real*8 :: qx, qy, qz, cos_sum, sin_sum, S_direct, q_dot_r
         real*8, dimension(NN) :: S_q_current, S_q_mean
         real*8 :: numerator
         integer, dimension(NN) :: bin_counts
         real*8 :: factor(2), sumquad


         if (ctn > rdf_mean) ctn = 1

         S_q_current = 0.0d0
         S_q_mean = 0.0d0
         bin_counts = 0
         f = 1

         ! Main loop over q-vectors (max_vect)
         do j = 1, max_vect
            qx = q_vec(j, 1)
            qy = q_vec(j, 2)
            qz = q_vec(j, 3)
            bin_idx = int(q_vec(j, 4))
            if (bin_idx < 1) then
               print *, 'Warning: bin_idx < 1, skipping q-vector', j
               cycle
            end if
            cos_sum = 0.0d0
            sin_sum = 0.0d0

            ! Parallelize over atoms

            do i = 1, n
                  q_dot_r = qx * x(i) + qy * y(i) + qz * z(i)
                  cos_sum = cos_sum + cos(q_dot_r)
                  sin_sum = sin_sum + sin(q_dot_r)
            end do


            ! Structure factor for this q-vector
            S_direct = (cos_sum**2 + sin_sum**2) / n
            S_q_current(bin_idx) = S_q_current(bin_idx) + S_direct
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1
         end do
         
         ! Normalize by bin count
         do i = 1, NN
            if (bin_counts(i) > 0) then
                  S_q_current(i) = S_q_current(i) / dble(bin_counts(i))
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
            open(unit, file='st.txt', status='unknown'
     &       , position='append')
            write(unit, 301) istep
301     format(/,' Step: ',i8,/)

            do i = 1, (NN)
                  S_q_mean(i) = 0.0d0
                  do j = 1, rdf_mean
                     S_q_mean(i) = S_q_mean(i) + S_hist(i, j, f)
                  end do
                  S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)

                  ! Add atomic form factors if needed (uncomment if you have atomF)
                  factor(1) = atomF('H ', q_magnitude(i))
                  factor(2) = atomF('O ', q_magnitude(i))
                  sumquad = (2.0d0 / 3.0d0) * factor(1)**2 
     &            + (1.0d0 / 3.0d0) * factor(2)**2
                  numerator = factor(1) * factor(2)

c                 S_q_mean(i) = S_q_mean(i) * (numerator / sumquad)

                  ! Clipping
c                  if (S_q_mean(i) > 3.0d0) S_q_mean(i) = 3.0d0
                  if (S_q_mean(i) == 0.0d0) cycle

                  write(unit, 310) q_magnitude(i), S_q_mean(i)
310         format(3x, f12.4, 3x, f12.4)
            end do
            close(unit)
         end if

         ctn = ctn + 1
         return
      end subroutine structfactor_direct

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
         use atom_sort
         implicit none

         integer :: i, j, f, a, bin_idx
         integer :: istep, freeunit, unit
         real*8 :: qx, qy, qz, S_direct, q_dot_r
         real*8, dimension(NN) :: S_q_current, S_q_mean
         integer, dimension(NN) :: bin_counts
         real*8, dimension(n_species) :: C_species, S_species
         
         real*8 :: cos_sum_w, sin_sum_w 
!        real*8 :: favg, qmg, fa
         
         print *, 'step',istep
         if (rdf_sort) then
            call apply_sorted_order() 
         end if


         if (ctn > rdf_mean) ctn = 1

         S_q_current = 0.0d0
         S_q_mean = 0.0d0
         bin_counts = 0
         f = 1
c
c        Main loop over q-vectors (max_vect) 
c
!$OMP PARALLEL DO DEFAULT(NONE) NUM_THREADS(12)
!$OMP& SHARED(max_vect, q_vec, n)
!$OMP& SHARED(x, y, z, atype, pc_fa, pc_favg_sq, n_species) 
!$OMP& PRIVATE(j, qx, qy, qz, bin_idx, i, q_dot_r, a)  
!$OMP& PRIVATE(C_species, S_species, cos_sum_w, sin_sum_w, S_direct)
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

            C_species = 0.0d0
            S_species = 0.0d0

            do i = 1, n
               q_dot_r = qx*x(i) + qy*y(i) + qz*z(i)
               C_species(atype(i)) = C_species(atype(i)) + cos(q_dot_r)
               S_species(atype(i)) = S_species(atype(i)) + sin(q_dot_r)
            end do

            cos_sum_w = 0.0d0
            sin_sum_w = 0.0d0


            do a = 1, n_species
                ! Look up pre-computed form factor
                cos_sum_w = cos_sum_w + pc_fa(j, a) * C_species(a)
                sin_sum_w = sin_sum_w + pc_fa(j, a) * S_species(a)
            end do

            if (pc_favg_sq(j) == 0.0d0) then
                S_direct = 0.0d0
c                print *, 'warning: precomputed_favg_sq = 0', 
c     &                   ' for q-vector', j
            else
c               This is (Total Scattering) / (N * <f^2>)
                S_direct = (cos_sum_w**2 + sin_sum_w**2)  
     &                 / ( dble(n) * pc_favg_sq(j) )
            end if

            S_q_current(bin_idx) = S_q_current(bin_idx) + S_direct
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
            open(unit, file='st.txt', status='unknown'
     &       , position='append')
            write(unit, 301) istep
301     format(/,' Step: ',i8,/)

            do i = 1, (NN)
                  S_q_mean(i) = 0.0d0
                  do j = 1, rdf_mean
                     S_q_mean(i) = S_q_mean(i) + S_hist(i, j, f)
                  end do
                  S_q_mean(i) = S_q_mean(i) / dble(rdf_mean)

                  ! Clipping
c                  if (S_q_mean(i) > 3.0d0) S_q_mean(i) = 3.0d0
                  if (S_q_mean(i) == 0.0d0) cycle
                  write(unit, 310) q_magnitude(i), S_q_mean(i)
310         format(3x, f12.4, 3x, f12.4)
            end do
            close(unit)
         end if

         ctn = ctn + 1
         if (rdf_sort) then
         call apply_original_order()
         end if
         return            
      return      
      
      end



!            favg = 0.0d0
            
!            qmg = sqrt(qx*qx + qy*qy + qz*qz)         ! Just needed for atomic form factors

!            do a = 1, n_species
!               fa = atomF_by_index(a, qmg)
!               cos_sum_w = cos_sum_w + fa * C_species(a)
!               sin_sum_w = sin_sum_w + fa * S_species(a)
!               favg = favg + mole_fractions(a) * fa**2    ! x_a * f_a(q)
!            end do

!            if (favg == 0.0d0) then
!                S_direct = 0.0d0
!                print *, 'warning: favg = 0 at q=', qmg
!            else
!               S_direct = (cos_sum_w**2 + sin_sum_w**2) 
!     &      / ( dble(n) * favg )
!            end if      