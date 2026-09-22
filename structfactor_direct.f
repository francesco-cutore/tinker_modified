      subroutine structfactor_forces(istep, fx_s, fy_s, fz_s,
     &   total_energy)
         use atomid
         use atoms
         use factors
         use rdfparams
         implicit none

         ! Arguments
         integer, intent(in) :: istep
         real*8, dimension(n), intent(inout) :: fx_s, fy_s, fz_s
         real*8, intent(inout) :: total_energy

         ! Local variables
         integer :: i, j, a, bin_idx
         integer :: freeunit, unit
         real*8 :: qx, qy, qz, q_dot_r

         ! Arrays for Two-Pass Algorithm
         real*8, allocatable :: saved_Eq(:), saved_Gq(:)
         real*8, dimension(NN) :: S_bin_avg, S_diff
         integer, dimension(NN) :: bin_counts
         real*8, dimension(n_species) :: C_species, S_species

         ! Scalar temps
         real*8 :: cos_sum_w, sin_sum_w, S_raw, S_normalized
         real*8 :: S_exp, prefactor, term_energy
         real*8 :: cos_i, sin_i, f_R_i
         real*8 :: dE_dri_x, dE_dri_y, dE_dri_z
         real*8 :: dG_dri_x, dG_dri_y, dG_dri_z
         real*8 :: dS_dri_x, dS_dri_y, dS_dri_z
         real*8 :: deriv_weight

         ! Weighted-restraint normalization (see SCATTER-NORMALIZE)
         real*8 :: Wnorm, esum, norm_denom
         integer :: n_act

         ! Parameters
         real*8 :: k_force
         integer, parameter :: savecycles = 100

         ! Output History Array
         real*8, dimension(NN) :: S_q_mean_out

         ! Force Constant
         k_force = rdf_kappa

         ! Initialize Outputs
         fx_s = 0.0d0
         fy_s = 0.0d0
         fz_s = 0.0d0
         total_energy = 0.0d0
         
         ! Initialize Bins
         S_bin_avg = 0.0d0
         bin_counts = 0
         
         ! Allocate Temporary Storage for Pass 1 -> Pass 2

         allocate(saved_Eq(max_vect))
         allocate(saved_Gq(max_vect))

         saved_Eq = 0.0d0
         saved_Gq = 0.0d0


!$OMP PARALLEL DO DEFAULT(NONE) NUM_THREADS(16)
!$OMP& SHARED(max_vect, q_vec, n, x, y, z, atype, n_species, NN)
!$OMP& SHARED(pc_fa, pc_favg_sq, pc_self_term)
!$OMP& SHARED(saved_Eq, saved_Gq)
!$OMP& PRIVATE(j, qx, qy, qz, bin_idx, i, q_dot_r, a)
!$OMP& PRIVATE(C_species, S_species, cos_sum_w, sin_sum_w)
!$OMP& PRIVATE(S_raw, S_normalized)
!$OMP& REDUCTION(+:S_bin_avg) REDUCTION(+:bin_counts)

         ! ---------------------------------------------------------
         ! PASS 1: Calculate S(q) and Accumulate Bin Averages
         ! ---------------------------------------------------------         
         
         do j = 1, max_vect                        ! For each q-vector                  
            qx = q_vec(j, 1)
            qy = q_vec(j, 2)
            qz = q_vec(j, 3)
            bin_idx = int(q_vec(j, 4))
            
            if (bin_idx < 1 .or. bin_idx > NN) cycle

            ! 1. Calculate Cos/Sin Sums 
            C_species = 0.0d0
            S_species = 0.0d0
            do i = 1, n
               q_dot_r = qx*x(i) + qy*y(i) + qz*z(i)
               C_species(atype(i)) = C_species(atype(i)) + cos(q_dot_r) 
               S_species(atype(i)) = S_species(atype(i)) + sin(q_dot_r)
            end do

            ! 2. Apply Form Factors 
            cos_sum_w = 0.0d0
            sin_sum_w = 0.0d0
            do a = 1, n_species
                cos_sum_w = cos_sum_w + pc_fa(j, a) * C_species(a)
                sin_sum_w = sin_sum_w + pc_fa(j, a) * S_species(a)
            end do

            ! 3. Store for Pass 2
            saved_Eq(j) = cos_sum_w
            saved_Gq(j) = sin_sum_w
            
            ! 4. Compute S_normalized
            S_raw = (cos_sum_w**2 + sin_sum_w**2) / dble(n)
            
            if (pc_favg_sq(j) .ne. 0.0d0) then
                S_normalized = S_raw / pc_favg_sq(j)
                S_normalized = S_normalized - pc_self_term(j)
            else
                S_normalized = 0.0d0
            end if

            ! 5. Accumulate into Bins
            S_bin_avg(bin_idx) = S_bin_avg(bin_idx) + S_normalized
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1
         end do
!$OMP END PARALLEL DO

      
         

         ! ---------------------------------------------------------
         ! PASS 2: Normalize This Step's Instantaneous S(q) per Bin,
         ! then Update the Exponentially-Decaying Running Average the
         ! Restraint Actually Acts On
         ! ---------------------------------------------------------
c
c     Rather than a plain rectangular (box-car) moving average -- which
c     has a small force discontinuity every time the oldest sample
c     falls out of a fixed-width window -- the restraint uses an
c     exponential moving average (EWMA): a memory function that decays
c     smoothly and indefinitely into the past instead of cutting off
c     sharply at a window edge (Torda, Scheek & van Gunsteren, 1989
c     -style time-averaged restraint, here generalized from distance
c     restraints to S(q)):
c
c         <S>(t) = alpha*S_inst(t) + (1-alpha)*<S>(t-1)
c
c     with alpha = scatter_alpha derived from the requested window
c     size in radialask.f. There is no <S>(t-1) on the very first
c     restrained step, so <S> is seeded directly from S_inst there
c     (equivalent to using alpha=1 for that one step only).
c
c     Differentiating the recursion with respect to the CURRENT
c     positions: <S>(t-1) was computed at an earlier step from
c     positions that no longer depend on r_i(t), so only the S_inst(t)
c     term survives the chain rule, each with weight "deriv_weight"
c     below (alpha normally, or 1 on the seeding step):
c
c         d<S>(t)/dr_i(t) = deriv_weight * dS_inst(t)/dr_i(t)
c
         do i = 1, NN
             if (bin_counts(i) > 0) then
                 S_bin_avg(i) = S_bin_avg(i) / dble(bin_counts(i))
             endif
         end do

         if (.not. ewma_initialized) then
            S_bin_ewma = S_bin_avg
            deriv_weight = 1.0d0
            ewma_initialized = .true.
         else
            S_bin_ewma = scatter_alpha * S_bin_avg
     &                   + (1.0d0 - scatter_alpha) * S_bin_ewma
            deriv_weight = scatter_alpha
         end if

         ! ---------------------------------------------------------
         ! Compute Deviations of the Running Average from Experiment
         ! ---------------------------------------------------------
         S_diff = 0.0d0
         esum = 0.0d0
         Wnorm = 0.0d0
         n_act = 0
         do i = 1, NN
             if (bin_counts(i) > 0) then
                 if (exp_present) then
                 ! Calculate Difference (Running Average - Experiment)
                 S_exp = S_exp_binned(i)
                 if (S_exp .ne. 0.0d0) then
                     S_diff(i) = S_bin_ewma(i) - S_exp

                     ! Accumulate the weighted sum of squares and the
                     ! weight total W = sum(w_i) over contributing bins.
                     ! n_act counts them (= W when all weights are 1).
                     esum = esum + S_exp_weight(i)*(S_diff(i)**2)
                     Wnorm = Wnorm + S_exp_weight(i)
                     n_act = n_act + 1
                 end if
                 end if
             endif
         end do

         ! Restraint normalization: with SCATTER-NORMALIZE on, divide by
         ! W = sum(w_i) so k is a penalty per unit mean-square deviation,
         ! independent of the number of bins / q-range / weighting. Off
         ! by default (norm_denom = 1) -> exact original 0.5*k*sum(...).
         if (scatter_normalize) then
             norm_denom = Wnorm
             if (norm_denom .le. 0.0d0) norm_denom = 1.0d0
         else
             norm_denom = 1.0d0
         end if

         ! Energy V = 0.5 * (k / W) * sum_i w_i * (S_diff_i)^2
         total_energy = 0.5d0 * (k_force / norm_denom) * esum

         ! Weighted RMS deviation from experiment (independent of k and of
         ! the normalization switch) -- the restraint convergence metric
         ! gradient.f logs to *_energy.csv.
         if (Wnorm .gt. 0.0d0) then
            scatter_rms = sqrt(esum / Wnorm)
         else
            scatter_rms = 0.0d0
         end if

         if (verbose_global .and. istep == 1) then
            print *, ' SCATTER restraint normalization:'
            print *, '   SCATTER-NORMALIZE : ', scatter_normalize
            print *, '   active bins N_act : ', n_act
            print *, '   weight sum W      : ', Wnorm
            print *, '   effective k (k/W) : ', k_force / norm_denom
         end if

      if (exp_present) then
!$OMP PARALLEL DO DEFAULT(NONE) NUM_THREADS(16)
!$OMP& SHARED(max_vect, q_vec, n, x, y, z, atype, NN)
!$OMP& SHARED(pc_fa, pc_favg_sq, k_force, deriv_weight)
!$OMP& SHARED(saved_Eq, saved_Gq, S_diff, bin_counts)
!$OMP& SHARED(S_exp_weight, norm_denom)
!$OMP& PRIVATE(j, qx, qy, qz, bin_idx, i, a, q_dot_r)
!$OMP& PRIVATE(cos_sum_w, sin_sum_w, prefactor)
!$OMP& PRIVATE(cos_i, sin_i, f_R_i)
!$OMP& PRIVATE(dE_dri_x, dE_dri_y, dE_dri_z)
!$OMP& PRIVATE(dG_dri_x, dG_dri_y, dG_dri_z)
!$OMP& PRIVATE(dS_dri_x, dS_dri_y, dS_dri_z)
!$OMP& REDUCTION(+:fx_s, fy_s, fz_s)
         

         ! ---------------------------------------------------------
         ! PASS 2: Calculate Forces Distributed to Atoms
         ! ---------------------------------------------------------
         do j = 1, max_vect
            bin_idx = int(q_vec(j, 4))
            
            ! Skip if invalid bin or if S_exp was zero 
            if (bin_idx < 1 .or. bin_idx > NN) cycle
            if (S_diff(bin_idx) == 0.0d0) cycle
            
            qx = q_vec(j, 1)
            qy = q_vec(j, 2)
            qz = q_vec(j, 3)

            ! Retrieve Saved E and G
            cos_sum_w = saved_Eq(j)
            sin_sum_w = saved_Gq(j)

            ! dE/dr = k * (<S> - S_exp) * (1/N_count)
            !         * deriv_weight * dS_vector/dr
            ! (Tinker convention: dex holds dE/dr, not the force -dE/dr;
            ! the sign flip to get the force is applied later by the
            ! integrator, e.g. beeman.f's a = -ekcal*derivs/mass)
            !
            ! The deriv_weight factor comes from the chain rule through
            ! the exponential running average -- see the comment above
            ! Pass 2's EWMA update for the derivation.
            prefactor = k_force * S_diff(bin_idx) * deriv_weight

            ! Same per-bin weight and 1/W normalization applied to the
            ! energy above -- keeps the force an exact gradient of E.
            prefactor = prefactor*S_exp_weight(bin_idx) / norm_denom

            if (bin_counts(bin_idx) > 0) then
               prefactor = prefactor / dble(bin_counts(bin_idx))
            endif

            ! Normalization from S definition
            if (pc_favg_sq(j) /= 0.0d0) then
                 prefactor = prefactor / (dble(n) * pc_favg_sq(j))
            endif

            ! Atom Loop for Derivatives
            do i = 1, n
               q_dot_r = qx*x(i) + qy*y(i) + qz*z(i)
               cos_i = cos(q_dot_r)
               sin_i = sin(q_dot_r)
               a = atype(i)
               f_R_i = pc_fa(j, a)

               ! Derivatives
               dE_dri_x = -qx * f_R_i * sin_i
               dE_dri_y = -qy * f_R_i * sin_i
               dE_dri_z = -qz * f_R_i * sin_i

               dG_dri_x = qx * f_R_i * cos_i
               dG_dri_y = qy * f_R_i * cos_i
               dG_dri_z = qz * f_R_i * cos_i

               ! Chain Rule
               dS_dri_x = 2.0d0*(cos_sum_w * dE_dri_x + 
     &                           sin_sum_w * dG_dri_x)
               dS_dri_y = 2.0d0*(cos_sum_w * dE_dri_y + 
     &                           sin_sum_w * dG_dri_y)
               dS_dri_z = 2.0d0*(cos_sum_w * dE_dri_z + 
     &                           sin_sum_w * dG_dri_z)

               ! Accumulate derivative contribution from this q-vector
               fx_s(i) = fx_s(i) + prefactor * dS_dri_x
               fy_s(i) = fy_s(i) + prefactor * dS_dri_y
               fz_s(i) = fz_s(i) + prefactor * dS_dri_z
            end do
         end do
!$OMP END PARALLEL DO

         end if

         ! Cleanup Temps
         deallocate(saved_Eq, saved_Gq)
         if (debug_global) then

         if (istep == 1) then
            print *, 'Writing initial S(q)'

            ! The EWMA is seeded directly from the instantaneous S(q)
            ! on step 1 (see the Pass 2 comment above), so this is
            ! exactly the starting configuration's S(q)
            S_q_mean_out = S_bin_ewma

            ! 3. Write
            unit = freeunit()
            open(unit, file=trim(scatter_basetag)//'_st_init.csv',
     &           status='unknown',
     &           position='append', action='write')

            do i = 1, NN
               if (S_q_mean_out(i) .ne. 0.0d0) then
                  write(unit, 310) q_magnitude(i), S_q_mean_out(i)
               end if
            end do

 310        format(f10.4, ';', f10.4)
            flush(unit)
            close(unit)
         end if

          if (mod(istep, savecycles) == 0) then
            print *, 'Writing S(q) at step', istep

            ! Reuse the running average already computed above
            S_q_mean_out = S_bin_ewma

            ! 3. Write
            unit = freeunit()
c            open(unit, file=trim(scatter_tag)//'_st.csv',
c     &           status='unknown',
c     &           position='append', action='write')
            open(unit, file=trim(scatter_tag)//'_st.csv',
     &           status='unknown',
     &           position='asis', action='write')
             write(unit, 181) 
181          format('q_value(A^-1);S_q_avg')
            
            do i = 1, NN
               if (S_q_mean_out(i) .ne. 0.0d0) then
                  write(unit, 320) q_magnitude(i), S_q_mean_out(i)
               end if
            end do

 320        format(f10.4, ';', f10.4)
            flush(unit)
            close(unit)

            ! S(q,t) time series (item 27): the _st.csv write above is a
            ! rewrite-in-place, so it only ever preserves the LAST dump.
            ! Append the same snapshot here instead, tagged with the step,
            ! keeping both the instantaneous S(q) (S_bin_avg, this step's
            ! value) and the EWMA (S_bin_ewma, what the restraint biases).
            unit = freeunit()
            open(unit, file=trim(scatter_tag)//'_st_traj.csv',
     &           status='unknown',
     &           position='append', action='write')

            do i = 1, NN
               if (S_bin_ewma(i) .ne. 0.0d0) then
                  write(unit, 330) istep, q_magnitude(i),
     &                             S_bin_avg(i), S_bin_ewma(i)
               end if
            end do

 330        format(i10, ';', f10.4, ';', f12.4, ';', f12.4)
            flush(unit)
            close(unit)
         end if
      end if

         return
      end subroutine structfactor_forces
