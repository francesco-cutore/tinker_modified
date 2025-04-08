c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong and Jay William Ponder  ##
c     ##                      Modified by FC                      ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################################
c     ##                                                                          ##
c     ##  subroutine radial  --  compute radial distribution function on the fly  ##
c     ##                                                                          ##
c     ##############################################################################
c
c
c     "radialsub" finds the radial distribution function for a specified
c     pair of atom types via analysis of a set of coordinate frames
c
c
      subroutine radialsub (istep)
      use argue
      use atomid
      use atoms
      use bound
      use boxes
      use files
      use inform
      use iounit
      use limits
      use math
      use molcul
      use potent
      use rdfparams
      implicit none
      integer i,j,k,f
      integer nframe
      integer molj,molk
      integer numj,numk
      integer typej,typek
      integer nbin,bin
      integer freeunit
      integer unit
      integer istep
      real*8 xj,yj,zj
      real*8 dx,dy,dz
      real*8 rjk,rmax,width
      real*8 rlower,rupper
      real*8 factor,pairs
      real*8 volume,expect
      logical intramol
      character*3 namej,namek
      character*6 labelj,labelk

      logical use_vdw_temp
      logical use_charge_temp
      logical use_dipole_temp
      logical use_mpole_temp
      logical use_ewald_temp
      real*8 vdwcut_temp
c
c     store current cutoffs
c
         use_vdw_temp = use_vdw
         use_charge_temp = use_charge
         use_dipole_temp = use_dipole
         use_mpole_temp = use_mpole
         use_ewald_temp = use_ewald
         vdwcut_temp = vdwcut
c
c     set cutoffs small to enforce use of minimum images
c
      use_vdw = .true.
      use_charge = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_ewald = .false.
      vdwcut = 0.01d0
c!!      call lattice
c
c     get maximum distance from input or minimum image convention
c
      rmax = rdf_rmax
c
c     get the desired width of the radial distance bins
c
      width = rdf_width
c
c     decide whether to restrict to intermolecular atom pairs
c
      intramol = rdf_intramol
c
c     set the number of distance bins
c
      nbin = rdf_nbin
c
c     if ctn is greater than rdf_mean then ctn = 1
c
      if (ctn > rdf_mean) then
         ctn = 1
      end if
c          
c     mypart :):)
c
      istep = istep
      do f = 1, rdf_num
         namej = rdf_namej(f)
         namek = rdf_namek(f)
         typej = rdf_typej(f)
         typek = rdf_typek(f)
         labelj = rdf_labelj(f)
         labelk = rdf_labelk(f)

         do j = 1, n
            if (name(j).eq.namej .or. 
     &          type(j).eq.typej) then
               xj = x(j)
               yj = y(j)
               zj = z(j)
               molj = molcule(j)
               do k = 1, n
                  if (name(k).eq.namek .or. 
     &                type(k).eq.typek) then
                     if (j .ne. k) then
                        molk = molcule(k)
                        if (intramol .or. molj.ne.molk) then
                           dx = x(k) - xj
                           dy = y(k) - yj
                           dz = z(k) - zj
                           call image (dx,dy,dz)             
                           rjk = sqrt(dx*dx + dy*dy + dz*dz)
                           bin = int(rjk/width) + 1
                           if (bin .le. nbin)
     &                         hist(bin,ctn,f) = hist(bin,ctn,f) + 1
                        end if
                     end if
                  end if
               end do
            end if
         end do
c
c     count the number of ourrences of each atom type
c
         numj = 0
         numk = 0
         do i = 1, n
            if (name(i).eq.namej .or. 
     &       type(i).eq.typej)  numj = numj + 1
            if (name(i).eq.namek .or. 
     &       type(i).eq.typek)  numk = numk + 1
         end do
c
c     normalize the distance bins to give radial distribution
c
         nframe = 1
         if (numj.ne.0 .and. numk.ne.0) then
            factor = (4.0d0/3.0d0) * pi * dble(nframe)
            if (use_bounds) then
               pairs = dble(numj) * dble(numk)
               volume = (gamma_sin*gamma_term) * xbox * ybox * zbox
               if (octahedron)  volume = 0.5d0 * volume
               if (dodecadron)  volume = volume / root2
               factor = factor * pairs / volume
            end if
            do i = 1, nbin
               rupper = dble(i) * width
               rlower = rupper - width
               expect = factor * (rupper**3 - rlower**3)
               gr(i,ctn,f) = dble(hist(i,ctn,f)) / expect
               hist(i,ctn,f) = 0
            end do
         end if   
c
c     find the 5th degree polynomial smoothed distribution function
c
      if (nbin .ge. 5) then
         gs(1,ctn,f) = (69.0d0*gr(1,ctn,f) + 4.0d0*gr(2,ctn,f) 
     &             - 6.0d0*gr(3,ctn,f)   + 4.0d0*gr(4,ctn,f)
     &             - gr(5,ctn,f)) / 70.0d0
         gs(2,ctn,f) = (2.0d0*gr(1,ctn,f) + 27.0d0*gr(2,ctn,f)
     &             + 12.0d0*gr(3,ctn,f) - 8.0d0*gr(4,ctn,f)
     &             + 2.0d0*gr(5,ctn,f)) / 35.0d0
         do i = 3, nbin-2
            gs(i,ctn,f) = (-3.0d0*gr((i-2),ctn,f) 
     &                + 12.0d0*gr((i-1),ctn,f)
     &                + 17.0d0*gr(i,ctn,f)
     &                + 12.0d0*gr((i+1),ctn,f)
     &                - 3.0d0*gr((i+2),ctn,f)) / 35.0d0
         end do
         gs((nbin-1),ctn,f) = (2.0d0*gr((nbin-4),ctn,f) 
     &                    - 8.0d0*gr((nbin-3),ctn,f)
     &                    + 12.0d0*gr((nbin-2),ctn,f) 
     &                    + 27.0d0*gr((nbin-1),ctn,f)
     &                    + 2.0d0*gr(nbin,ctn,f)) / 35.0d0
         gs(nbin,ctn,f) = (-gr((nbin-4),ctn,f) 
     &                + 4.0d0*gr((nbin-3),ctn,f) 
     &                - 6.0d0*gr((nbin-2),ctn,f) 
     &                + 4.0d0*gr((nbin-1),ctn,f)
     &                + 69.0d0*gr(nbin,ctn,f)) 
     &                / 70.0d0
         do i = 1, nbin
            gs(i,ctn,f) = max(0.0d0,gs(i,ctn,f))
         end do
      end if
      

      if (ctn .eq. rdf_mean) then
         savelock = 1
      endif
c
c     open the output file only to write the mean
c
         if (savelock .eq. 1) then
               gr_mean = 0.0d0
               gs_mean = 0.0d0
            do i = 1, size(gr, 1)
               do j = 1, rdf_mean
                  gr_mean(i,f) = gr_mean(i,f) + gr(i,j,f)
                  gs_mean(i,f) = gs_mean(i,f) + gs(i,j,f)
               end do
               gr_mean(i,f) = gr_mean(i,f) / dble(rdf_mean)
               gs_mean(i,f) = gs_mean(i,f) / dble(rdf_mean)
            end do

            unit = freeunit ()
            open(unit,file='rdf.txt',status='unknown',position='append')
            write (unit,220)  labelj,labelk
  220    format (/,' Pairwise Radial Distribution Function :'
     &        //,7x,'First Name or Type :  ',a6,
     &           5x,'Second Name or Type :  ',a6)
            write (unit,230)
  230    format (/,5x,'Bin',9x,'Counts',7x,'Distance',7x,'Raw g(r)',
     &           4x,'Smooth g(r)',/)
            do i = 1, nbin
               write (unit,240)  i,hist(i,ctn,f),
     &         (dble(i)-0.5d0)*width,gr_mean(i,f),gs_mean(i,f)
  240    format (i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)
            end do 
            close(unit)
         end if 
      end do     

      ctn = ctn + 1
c
c     reset the cutoffs
c      
         use_vdw = use_vdw_temp
         use_charge = use_charge_temp
         use_dipole = use_dipole_temp
         use_mpole = use_mpole_temp
         use_ewald = use_ewald_temp
         vdwcut = vdwcut_temp 
         
      end   