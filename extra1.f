c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################

      subroutine extra1
      use atoms
      use deriv
      use energi
      use rdfparams
      implicit none
      integer i

      real*8, allocatable :: fx_s(:), fy_s(:), fz_s(:)
      real*8 :: energy_s

c     zero out the extra energy term and first derivatives
c
      ex = 0.0d0
      do i = 1, n
         dex(1,i) = 0.0d0
         dex(2,i) = 0.0d0
         dex(3,i) = 0.0d0
      end do

      if (.not. use_scatter) return
      if (current_md_step == 0) return

c     1. Allocate local arrays to size 'n' (number of atoms)
      allocate (fx_s(n))
      allocate (fy_s(n))
      allocate (fz_s(n))

c     2. Initialize them to zero
      fx_s = 0.0d0
      fy_s = 0.0d0
      fz_s = 0.0d0
      energy_s = 0.0d0

c     3. Call your routine
      call structfactor_forces(current_md_step,fx_s, fy_s,
     &      fz_s, energy_s)

c     4. Accumulate Energy into global 'ex'
      ex = ex + energy_s

c     5. Accumulate Forces into global 'dex'
      do i = 1, n
         dex(1,i) = dex(1,i) + fx_s(i)
         dex(2,i) = dex(2,i) + fy_s(i)
         dex(3,i) = dex(3,i) + fz_s(i)
      end do

c     6. Clean up
      deallocate (fx_s, fy_s, fz_s)

      return
      end