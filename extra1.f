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
      use files
      implicit none
      integer i

      real*8, allocatable :: fx_s(:), fy_s(:), fz_s(:)
      real*8 :: energy_s
      real*8 :: f_sq_sum, f_rms
      integer log_unit
      integer :: freeunit
      logical :: file_exists
      
c      return
c
c     zero out the extra energy term and first derivatives
c
      ex = 0.0d0
      do i = 1, n
         dex(1,i) = 0.0d0
         dex(2,i) = 0.0d0
         dex(3,i) = 0.0d0
      end do

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
      
c     Calculate RMS force for logging
c     sqrt(sum_i (fx_i^2 + fy_i^2 + fz_i^2) / N)
      if (.false.) then
      f_sq_sum = 0.0d0
      do i = 1, n
         f_sq_sum = f_sq_sum + fx_s(i)**2 + fy_s(i)**2 + fz_s(i)**2
      end do
      
      if (n .gt. 0) then
         f_rms = sqrt(f_sq_sum / (dble(n)))
      else
         f_rms = 0.0d0
      end if
      if (mod(current_md_step, 10) == 0) then
      ! Print to standard output (screen)
      write(*, '("Scattering RMS Force: ", F15.6, " Energy: ", F15.6)') 
     &      f_rms, energy_s
      log_unit = freeunit()
      
      inquire(file='scattering_log.csv', exist=file_exists)
      
      open(unit=log_unit, file='scattering_log.csv', status='unknown', 
     &     position='append')
      
      if (.not. file_exists) then
         write(log_unit, '(A)') 'RMS_Force,Scattering_Energy'
      end if
      
      ! Write data 
      write(log_unit, '(",", F15.6, ",", F15.6)') 
     &       f_rms, energy_s
      close(log_unit)
      end if
      end if

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