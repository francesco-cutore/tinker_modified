c     ##############################################################################
c     ##                                                                          ##
c     ##  Optimized & Corrected structfactor_forces Routine                       ##
c     ##                                                                          ##
c     ##############################################################################

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
         integer :: i, j, f, a, bin_idx
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
         
         ! Parameters
         real*8 :: k_force
         integer, parameter :: savecycles = 100
         
         ! Output History Array
         real*8, dimension(NN) :: S_q_mean_out

         if (ctn > rdf_mean) ctn = 1
         f = 1
         
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
         ! PASS 2: Average Bins and Compute Deviations
         ! ---------------------------------------------------------
         S_diff = 0.0d0
         do i = 1, NN
             if (bin_counts(i) > 0) then
                 ! Normalize sum to get Average S for this bin
                 S_bin_avg(i) = S_bin_avg(i) / dble(bin_counts(i))
                 if (exp_present) then 
                 ! Calculate Difference (Average - Experiment)
                 S_exp = S_exp_binned(i)
                 if (S_exp .ne. 0.0d0) then
                     S_diff(i) = S_bin_avg(i) - S_exp 
                     
                     ! Calculate Energy (V = 0.5 * k * Diff^2)
                     term_energy = 0.5d0 * k_force * (S_diff(i)**2)
                     total_energy = total_energy + term_energy
                 end if
                 end if
             endif
         end do

      
        
         ! Store in History
         do i = 1, NN
             S_hist(i, ctn, f) = S_bin_avg(i)
         end do

      if (exp_present) then 
!$OMP PARALLEL DO DEFAULT(NONE) NUM_THREADS(16)
!$OMP& SHARED(max_vect, q_vec, n, x, y, z, atype, NN)
!$OMP& SHARED(pc_fa, pc_favg_sq, k_force)
!$OMP& SHARED(saved_Eq, saved_Gq, S_diff, bin_counts)
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

            ! Force F = -k * (S_avg - S_exp) * (1/N_count) * dS_vector/dr
            prefactor = -k_force * S_diff(bin_idx)
            
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

               ! Accumulate Force
               fx_s(i) = -fx_s(i) + prefactor * dS_dri_x
               fy_s(i) = -fy_s(i) + prefactor * dS_dri_y
               fz_s(i) = -fz_s(i) + prefactor * dS_dri_z
            end do
         end do
!$OMP END PARALLEL DO

         end if

         ! Cleanup Temps
         deallocate(saved_Eq, saved_Gq)
         if (debug_global) then     
         if (ctn == rdf_mean) savelock = 1

         if (savelock == 1 .and. istep == 1) then
            print *, 'Writing initial S(q)'
            
            ! 1. Initialize
            S_q_mean_out = 0.0d0

            ! 2. Compute Average 
            do j = 1, rdf_mean
               do i = 1, NN
                  S_q_mean_out(i) = S_q_mean_out(i) + S_hist(i, j, f)
               end do
            end do
            
            ! Normalize
            if (rdf_mean .gt. 0) then
               S_q_mean_out = S_q_mean_out / dble(rdf_mean)
            endif

            ! 3. Write
            unit = freeunit()
            open(unit, file='st_init.csv', status='unknown', 
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

          if (savelock == 1 .and. mod(istep, savecycles) == 0) then
            print *, 'Writing S(q) at step', istep
            
            ! 1. Initialize
            S_q_mean_out = 0.0d0

            ! 2. Compute Average 
            do j = 1, rdf_mean
               do i = 1, NN
                  S_q_mean_out(i) = S_q_mean_out(i) + S_hist(i, j, f)
               end do
            end do
            
            ! Normalize
            if (rdf_mean .gt. 0) then
               S_q_mean_out = S_q_mean_out / dble(rdf_mean)
            endif

            ! 3. Write
            unit = freeunit()
c            open(unit, file='st.csv', status='unknown', 
c     &           position='append', action='write')
            open(unit, file='st.csv', status='unknown', 
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
         end if
      end if
         ctn = ctn + 1
          
         return      
      end subroutine structfactor_forces