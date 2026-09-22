c
c ######################################################################
c ##                                                                  ##
c ##  Pre-computes form factors, denominators, and SELF-TERMS.        ##
c ##                                                                  ##
c ##  Updated to match TRAVIS "Distinct" Structure Factor output      ##
c ##  by calculating the self-scattering term for subtraction.        ##
c ##                                                                  ##
c ######################################################################
c
      subroutine precompute_form_factors(norm_mode, manual_val)
        use factors      
        use rdfparams    
        implicit none
        
        integer, intent(in) :: norm_mode
        real*8, intent(in)  :: manual_val
        
        integer :: j, a
        real*8 :: qx, qy, qz, qmg, fa
        real*8 :: sum_xf, sum_xf2
        
        ! Allocate global arrays (guard against re-allocation --
        ! this routine may be called more than once per run, e.g.
        ! by SCATTER-AUTOTUNE trying multiple normalization modes)

        if (allocated(pc_fa)) deallocate(pc_fa)
        if (allocated(pc_favg_sq)) deallocate(pc_favg_sq)
        if (allocated(pc_self_term)) deallocate(pc_self_term)

        allocate(pc_fa(max_vect, n_species))
        allocate(pc_favg_sq(max_vect))
        allocate(pc_self_term(max_vect))

        pc_fa = 0.0d0
        pc_favg_sq = 0.0d0
        pc_self_term = 0.0d0

        ! Loop over every q-vector
        do j = 1, max_vect
            qx = q_vec(j, 1)
            qy = q_vec(j, 2)
            qz = q_vec(j, 3)
            qmg = sqrt(qx*qx + qy*qy + qz*qz)

            sum_xf  = 0.0d0  ! Accumulator for <f>
            sum_xf2 = 0.0d0  ! Accumulator for <f^2>

            do a = 1, n_species
               fa = atomF_by_index(a, qmg)

               ! Store individual form factor
               pc_fa(j, a) = fa
               
               ! Accumulate weighted sums
               sum_xf  = sum_xf  + mole_fractions(a) * fa
               sum_xf2 = sum_xf2 + mole_fractions(a) * fa**2
               
            end do
            
            select case (norm_mode)
            
            case (1) 
               ! C++ Case 0: Average Atom (Ashcroft-Langreth)
               ! Denominator = <f>^2
               if (sum_xf .ne. 0.0d0) then
                   pc_favg_sq(j) = sum_xf**2
                   ! Self-term = <f^2> / <f>^2
                   pc_self_term(j) = sum_xf2 / (sum_xf**2)
               else
                   pc_favg_sq(j) = 1.0d0
                   pc_self_term(j) = sum_xf2
               endif

            case (2)
               ! C++ Case 1: Faber-Ziman (Standard)
               ! Denominator = <f^2>
               if (sum_xf2 .ne. 0.0d0) then
                   pc_favg_sq(j) = sum_xf2
                   ! Self-term = <f^2> / <f^2> = 1.0
                   pc_self_term(j) = 1.0d0
               else
                   pc_favg_sq(j) = 1.0d0
                   pc_self_term(j) = 1.0d0
               endif

            case (3)
               ! Manual Normalization
               if (manual_val .ne. 0.0d0) then
                   pc_favg_sq(j) = manual_val
                   pc_self_term(j) = sum_xf2 / manual_val
               else
                   pc_favg_sq(j) = 1.0d0
                   pc_self_term(j) = 1.0d0
               endif

            case (4)
               ! No Normalization
               pc_favg_sq(j) = 1.0d0
               pc_self_term(j) = sum_xf2
               
            case default
               print *, 'Error: Unknown normalization mode'
               stop
            end select

        end do

        if (verbose_global) then
            print *, '---------------------------------------------'
            print *, ' PRECOMPUTE FORM FACTORS SUMMARY '
            print *, ' Normalization Mode: ', norm_mode
            print *, '---------------------------------------------'
        end if
      end subroutine precompute_form_factors