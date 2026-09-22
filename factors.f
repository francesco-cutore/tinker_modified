! Subroutine factors: This routine computes the structure factor as
!                    a function of q for each element in the box.
!

        module factors
            use math
            implicit none
            ! Array to store form factor for each q-vector (j) and species (a)
            real*8, allocatable, dimension(:,:) :: pc_fa
    
            ! Array to store the denominator <f^2> for each q-vector (j)
            real*8, allocatable, dimension(:) :: pc_favg_sq

            real*8, allocatable :: pc_self_term(:)
        
        contains
        real(8) function atomF(nome,k)
            character nome*2
            real(8) k
c           q02's low-q Cromer-Mann fit is used for the full range;
c           see q_max in radialask.f -- kt=k/(4*pi) stays well within
c           q02's valid domain there, so no high-q blend is needed.
            atomF = q02(nome,k)
        return
        end function atomF

        integer function get_pair_index(alpha, beta, n_species)
        implicit none
        integer, intent(in) :: alpha, beta, n_species
        get_pair_index = (alpha - 1) * n_species + beta
        end function
        
        !==============================================================
        ! Wrapper: get atomic form factor by species index
        !==============================================================
        real*8 function atomF_by_index(a, q)
        use atomid
        implicit none
        integer, intent(in) :: a
        real*8, intent(in) :: q
        character :: a_name*2

        a_name = atom_name(a)
        atomF_by_index = atomF(a_name, q)
        end function atomF_by_index

!
!==============================================================
! NEW FUNCTION: Atom Number to Name Conversion Table
!==============================================================
      character*2 function get_name_from_atomic_number(atomic_num)
          implicit none
          integer, intent(in) :: atomic_num
          
c         Covers every element with q02() Cromer-Mann coefficients
c         (atomic numbers 1-30, H through Zn); q02 is the only form-
c         factor branch actually reachable from atomF (see factors.f).
          SELECT CASE (atomic_num)
              CASE (1)
                  get_name_from_atomic_number = 'H '
              CASE (2)
                  get_name_from_atomic_number = 'He'
              CASE (3)
                  get_name_from_atomic_number = 'Li'
              CASE (4)
                  get_name_from_atomic_number = 'Be'
              CASE (5)
                  get_name_from_atomic_number = 'B '
              CASE (6)
                  get_name_from_atomic_number = 'C '
              CASE (7)
                  get_name_from_atomic_number = 'N '
              CASE (8)
                  get_name_from_atomic_number = 'O '
              CASE (9)
                  get_name_from_atomic_number = 'F '
              CASE (10)
                  get_name_from_atomic_number = 'Ne'
              CASE (11)
                  get_name_from_atomic_number = 'Na'
              CASE (12)
                  get_name_from_atomic_number = 'Mg'
              CASE (13)
                  get_name_from_atomic_number = 'Al'
              CASE (14)
                  get_name_from_atomic_number = 'Si'
              CASE (15)
                  get_name_from_atomic_number = 'P '
              CASE (16)
                  get_name_from_atomic_number = 'S '
              CASE (17)
                  get_name_from_atomic_number = 'Cl'
              CASE (18)
                  get_name_from_atomic_number = 'Ar'
              CASE (19)
                  get_name_from_atomic_number = 'K '
              CASE (20)
                  get_name_from_atomic_number = 'Ca'
              CASE (21)
                  get_name_from_atomic_number = 'Sc'
              CASE (22)
                  get_name_from_atomic_number = 'Ti'
              CASE (23)
                  get_name_from_atomic_number = 'V '
              CASE (24)
                  get_name_from_atomic_number = 'Cr'
              CASE (25)
                  get_name_from_atomic_number = 'Mn'
              CASE (26)
                  get_name_from_atomic_number = 'Fe'
              CASE (27)
                  get_name_from_atomic_number = 'Co'
              CASE (28)
                  get_name_from_atomic_number = 'Ni'
              CASE (29)
                  get_name_from_atomic_number = 'Cu'
              CASE (30)
                  get_name_from_atomic_number = 'Zn'
              CASE DEFAULT
                  get_name_from_atomic_number = '??'
                  print *, 'Error: Unknown atomic number:', atomic_num
                  stop
          END SELECT
          
      end function get_name_from_atomic_number



        real(8) function q02(nome,k)
            real(8) a1, b1, a2, b2, a3, b3, a4, b4, c, k
            real(8) kt
            character nome*2

            SELECT CASE (nome)
            CASE ('H')
                a1 = 0.489918
                b1 = 20.65930
                a2 = 0.262003
                b2 = 7.740390
                a3 = 0.196767
                b3 = 49.55190
                a4 = 0.049879
                b4 = 2.201590
                c  = 0.001305
            CASE ('He')
                a1 = 0.87340
                b1 = 9.10370
                a2 = 0.63090
                b2 = 3.35680
                a3 = 0.31120
                b3 = 22.9276
                a4 = 0.17800
                b4 = 0.98210
                c  = 0.00640
            CASE ('Li')
                a1 = 1.12820
                b1 = 3.95460
                a2 = 0.75080
                b2 = 1.05240
                a3 = 0.61750
                b3 = 85.39050
                a4 = 0.46530
                b4 = 168.26100
                c  = 0.03770
            CASE ('Be')
                a1 = 1.59190
                b1 = 43.64270
                a2 = 1.12780
                b2 = 1.86230
                a3 = 0.53910
                b3 = 103.48300
                a4 = 0.70290
                b4 = 0.54200
                c  = 0.03850
            CASE ('B')
                a1 = 2.05450
                b1 = 23.21850
                a2 = 1.33260
                b2 = 1.02100
                a3 = 1.09790
                b3 = 60.34980
                a4 = 0.70680
                b4 = 0.14030
                c  = -0.19320
            CASE ('C')
                a1 = 2.31000
                b1 = 20.84390
                a2 = 1.02000
                b2 = 10.20750
                a3 = 1.58860
                b3 = 0.56870
                a4 = 0.86500
                b4 = 51.65120
                c  = 0.21560
            CASE ('N')
                a1 = 12.21260
                b1 = 0.00570
                a2 = 3.13220
                b2 = 9.89330
                a3 = 2.01250
                b3 = 28.99750
                a4 = 1.16630
                b4 = 0.58260
                c  = -11.52900 
            CASE ('O')
                a1 = 3.04850
                b1 = 13.27710
                a2 = 2.28680
                b2 = 5.70110
                a3 = 1.54630
                b3 = 0.32390
                a4 = 0.86700
                b4 = 32.90890
                c  = 0.25080
            CASE ('F')
                a1 = 3.53920
                b1 = 10.28250
                a2 = 2.64120
                b2 = 4.29440
                a3 = 1.51700
                b3 = 0.26150
                a4 = 1.02430
                b4 = 26.14760
                c  = 0.27760 
            CASE ('Ne')
                a1 = 3.95530
                b1 = 8.40420
                a2 = 3.11250
                b2 = 3.42620
                a3 = 1.45460
                b3 = 0.23060
                a4 = 1.12510
                b4 = 21.71840
                c  = 0.35150
            CASE ('Na')
                a1 = 3.25650
                b1 = 2.66710
                a2 = 3.93620
                b2 = 6.11530
                a3 = 1.39980
                b3 = 0.20010
                a4 = 1.00320
                b4 = 14.03900
                c  = 0.40400
            CASE ('Mg')
                a1 = 3.49880
                b1 = 2.16760
                a2 = 3.83780
                b2 = 4.75420
                a3 = 1.32840
                b3 = 0.18500
                a4 = 0.84970
                b4 = 10.14110
                c  = 0.48530
            CASE ('Al')
                a1 = 6.42020
                b1 = 3.03870
                a2 = 1.90020
                b2 = 0.74260
                a3 = 1.59360
                b3 = 31.54720
                a4 = 1.96460
                b4 = 85.08860
                c  = 1.11510
            CASE ('Si')
                a1 = 6.29150
                b1 = 2.43860
                a2 = 3.03530
                b2 = 32.33370
                a3 = 1.98910
                b3 = 0.67850
                a4 = 1.54100
                b4 = 81.69370
                c  = 1.14070
            CASE ('P')
                a1 = 6.43450
                b1 = 1.90670
                a2 = 4.17910
                b2 = 27.15700
                a3 = 1.78000
                b3 = 0.52600
                a4 = 1.49080
                b4 = 68.16450
                c  = 1.11490
            CASE ('S')
                a1 = 6.90530
                b1 = 1.46790
                a2 = 5.20340
                b2 = 22.21510
                a3 = 1.43790
                b3 = 0.25360
                a4 = 1.58630
                b4 = 56.17200
                c  = 0.86690
            CASE ('Cl')
                a1 = 11.46040
                b1 = 0.01040
                a2 = 7.19640
                b2 = 1.66200
                a3 = 6.25560
                b3 = 18.51940
                a4 = 1.64550
                b4 = 47.77840
                c  = -9.55740
            CASE ('Ar')
                a1 = 7.48450
                b1 = 0.90720
                a2 = 6.77230
                b2 = 14.84070
                a3 = 0.65390
                b3 = 43.89830
                a4 = 1.64420
                b4 = 33.39290
                c  = 1.44450
            CASE ('K')
                a1 = 8.21860
                b1 = 12.79490
                a2 = 7.43980
                b2 = 0.77480
                a3 = 1.05190
                b3 = 213.18700
                a4 = 0.86590
                b4 = 41.68410
                c  = 1.42280
            CASE ('Ca')
                a1 = 8.62660
                b1 = 10.44210
                a2 = 7.38730
                b2 = 0.65990
                a3 = 1.58990
                b3 = 85.74840
                a4 = 1.02110
                b4 = 178.43700
                c  = 1.37510
            CASE ('Sc')
                a1 = 9.18900
                b1 = 9.02130
                a2 = 7.36790
                b2 = 0.57290
                a3 = 1.64090
                b3 = 136.10800
                a4 = 1.46800
                b4 = 51.35310
                c  = 1.33290
            CASE ('Ti')
                a1 = 9.75950
                b1 = 7.85080
                a2 = 7.35580
                b2 = 0.50000
                a3 = 1.69910
                b3 = 35.63380
                a4 = 1.90210
                b4 = 116.10500
                c  = 1.28070
            CASE ('V')
                a1 = 10.29710
                b1 = 6.86570
                a2 = 7.35110
                b2 = 0.43850
                a3 = 2.07030
                b3 = 26.89380
                a4 = 2.05710
                b4 = 102.47800
                c  = 1.21990
            CASE ('Cr')
                a1 = 10.64060
                b1 = 6.10380
                a2 = 7.35370
                b2 = 0.39200
                a3 = 3.32400
                b3 = 20.26260
                a4 = 1.49220
                b4 = 98.73990
                c  = 1.18320
            CASE ('Mn')
                a1 = 11.28190
                b1 = 5.34090
                a2 = 7.35730
                b2 = 0.34320
                a3 = 3.01930
                b3 = 17.86740
                a4 = 2.24410
                b4 = 83.75430
                c  = 1.08960
            CASE ('Fe')
                a1 = 11.76950
                b1 = 4.76110
                a2 = 7.35730
                b2 = 0.30720
                a3 = 3.52220
                b3 = 15.35350
                a4 = 2.30450
                b4 = 76.88050
                c  = 1.03690
            CASE ('Co')
                a1 = 12.28410
                b1 = 4.27910
                a2 = 7.34090
                b2 = 0.27840
                a3 = 4.00340
                b3 = 13.53590
                a4 = 2.34880
                b4 = 71.16920
                c  = 1.01180
            CASE ('Ni')
                a1 = 12.83760
                b1 = 3.87850
                a2 = 7.29200
                b2 = 0.25650
                a3 = 4.44380
                b3 = 12.17630
                a4 = 2.38000
                b4 = 66.34210
                c  = 1.03410
            CASE ('Cu')
                a1 = 13.33800
                b1 = 3.58280
                a2 = 7.16760
                b2 = 0.24700
                a3 = 5.61580
                b3 = 11.39660
                a4 = 1.67350
                b4 = 64.81260
                c  = 1.19100
            CASE ('Zn')
                a1 = 14.07430
                b1 = 3.26550
                a2 = 7.03180
                b2 = 0.23330
                a3 = 5.16520
                b3 = 10.31630
                a4 = 2.41000
                b4 = 58.70970
                c  = 1.30410
            CASE DEFAULT
                q02 = 0
                print *, 'Error: Element not recognized', 
     &            ' in q02 function:', nome
                RETURN
            END SELECT
            !CALCULA O TERMO
            q02 = 0
            kt = k / (4*pi)
            q02 = a1 * exp(-b1 * kt ** 2)
            q02 = q02 + a2 * exp(-b2 * kt ** 2)
            q02 = q02 + a3 * exp(-b3 * kt ** 2)
            q02 = q02 + a4 * exp(-b4 * kt ** 2)
            q02 = q02 + c
            RETURN
            end function
      end module factors
    