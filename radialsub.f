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
      subroutine radialsub
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
      integer i,j,k
      integer nframe,iframe
      integer iarc,next
      integer molj,molk
      integer numj,numk
      integer typej,typek
      integer start,stop
      integer step,skip
      integer nbin,bin
      integer freeunit
      integer unit
c      integer, allocatable :: hist(:)
      real*8 xj,yj,zj
      real*8 dx,dy,dz
      real*8 rjk,rmax,width
      real*8 rlower,rupper
      real*8 factor,pairs
      real*8 volume,expect
c      real*8, allocatable :: gr(:)
c      real*8, allocatable :: gs(:)
      logical first,intramol
      character*3 namej,namek
      character*6 labelj,labelk
c
c     open the trajectory archive and read the initial frame
c
      call getarcmodified (iarc)
c
c     get the unitcell parameters and number of molecules
c
      call unitcell
      call molecule
c
c     set cutoffs small to enforce use of minimum images
c
c!!      use_vdw = .true.
c!!      use_charge = .false.
c!!      use_dipole = .false.
c!!      use_mpole = .false.
c!!     use_ewald = .false.
c!!     vdwcut = 0.01d0
c!!     call lattice
c
c     get numbers of the coordinate frames to be processed
c
      start = rdf_start
      stop = rdf_stop
      step = rdf_step   
c
c     get the names of the atoms to be used in rdf computation
c
      labelj = rdf_labelj
      labelk = rdf_labelk
c
c     convert the labels to either atom names or type numbers
c
      namej = '   '
      typej = -1
      read (labelj,*,err=70,end=70)  typej
   70 continue
      if (typej .le. 0) then
         next = 1
         call gettext (labelj,namej,next)
      end if
      namek = '   '
      typek = -1
      read (labelk,*,err=80,end=80)  typek
   80 continue
      if (typek .le. 0) then
         next = 1
         call gettext (labelk,namek,next)
      end if
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
c     count the number of coordinate frames in the archive file
c
      abort = .false.
      rewind (unit=iarc)
      first = .true.
      nframe = 0
      do while (.not. abort)
         call readcart (iarc,first)
         nframe = nframe + 1
      end do
      nframe = nframe - 1
      stop = min(nframe,stop)
      nframe = (stop-start)/step + 1
      write (iout,170)  nframe
  170 format (/,' Number of Coordinate Frames :',i14)

c
c     if ctn is greater than rdf_mean then ctn = 1
c
      if (ctn > rdf_mean) then
         ctn = 1
      end if    
c
c     get the archived coordinates for each frame in turn
c
      write (iout,190)
  190 format (/,' Reading the Coordinates Archive File :',/)
      rewind (unit=iarc)
      first = .true.
      nframe = 0
      iframe = start
      skip = start
      do while (iframe.ge.start .and. iframe.le.stop)
         do j = 1, skip-1
            call readcart (iarc,first)
         end do
         iframe = iframe + step
         skip = step
         call readcart (iarc,first)
         if (.not. abort) then
            nframe = nframe + 1
            if (mod(nframe,100) .eq. 0) then
               write (iout,200)  nframe
  200          format (4x,'Processing Coordinate Frame',i13)
            end if
            do j = 1, n
               if (name(j).eq.namej .or. type(j).eq.typej) then
                  xj = x(j)
                  yj = y(j)
                  zj = z(j)
                  molj = molcule(j)
                  do k = 1, n
                     if (name(k).eq.namek .or. type(k).eq.typek) then
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
     &                           hist(bin,ctn) = hist(bin,ctn) + 1
                           end if
                        end if
                     end if
                  end do
               end if
            end do
         end if
      end do
c
c     ensure a valid frame is loaded and report total frames
c
      if (abort) then
         rewind (unit=iarc)
         first = .true.
         call readcart (iarc,first)
      end if
      close (unit=iarc)
      if (mod(nframe,100) .ne. 0) then
         write (iout,210)  nframe
  210    format (4x,'Processing Coordinate Frame',i13)
      end if
c
c     count the number of ourrences of each atom type
c
      numj = 0
      numk = 0
      do i = 1, n
         if (name(i).eq.namej .or. type(i).eq.typej)  numj = numj + 1
         if (name(i).eq.namek .or. type(i).eq.typek)  numk = numk + 1
      end do
c
c     normalize the distance bins to give radial distribution
c
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
            gr(i,ctn) = dble(hist(i,ctn)) / expect
            hist(i,ctn) = 0
         end do
      end if
c
c     find the 5th degree polynomial smoothed distribution function
c
      if (nbin .ge. 5) then
         gs(1,ctn) = (69.0d0*gr(1,ctn) + 4.0d0*gr(2,ctn) 
     &             - 6.0d0*gr(3,ctn)   + 4.0d0*gr(4,ctn)
     &             - gr(5,ctn)) / 70.0d0
         gs(2,ctn) = (2.0d0*gr(1,ctn) + 27.0d0*gr(2,ctn)
     &             + 12.0d0*gr(3,ctn) - 8.0d0*gr(4,ctn)
     &             + 2.0d0*gr(5,ctn)) / 35.0d0
         do i = 3, nbin-2
            gs(i,ctn) = (-3.0d0*gr((i-2),ctn) + 12.0d0*gr((i-1),ctn)
     &                + 17.0d0*gr(i,ctn) + 12.0d0*gr((i+1),ctn)
     &                - 3.0d0*gr((i+2),ctn)) / 35.0d0
         end do
         gs((nbin-1),ctn) = (2.0d0*gr((nbin-4),ctn) 
     &                    - 8.0d0*gr((nbin-3),ctn)
     &                    + 12.0d0*gr((nbin-2),ctn) 
     &                    + 27.0d0*gr((nbin-1),ctn)
     &                    + 2.0d0*gr(nbin,ctn)) / 35.0d0
         gs(nbin,ctn) = (-gr((nbin-4),ctn) + 4.0d0*gr((nbin-3),ctn) 
     6                - 6.0d0*gr((nbin-2),ctn) + 4.0d0*gr((nbin-1),ctn)
     &                + 69.0d0*gr(nbin,ctn)) / 70.0d0
         do i = 1, nbin
            gs(i,ctn) = max(0.0d0,gs(i,ctn))
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
               gr_mean(i) = gr_mean(i) + gr(i, j)
               gs_mean(i) = gs_mean(i) + gs(i, j)
            end do
            gr_mean(i) = gr_mean(i) / dble(rdf_mean)
            gs_mean(i) = gs_mean(i) / dble(rdf_mean)
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
            write (unit,240)  i,hist(i,ctn),
     &         (dble(i)-0.5d0)*width,gr_mean(i),gs_mean(i)
  240    format (i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)
         end do 
         close(unit)
      end if   

cc         unit = freeunit ()
cc         open(unit,file='rdf.txt',status='unknown',position='append')
cc         write (unit,250)  labelj,labelk
cc  250    format (/,' Pairwise Radial Distribution Function :'
cc     &        //,7x,'First Name or Type :  ',a6,
cc     &           5x,'Second Name or Type :  ',a6)
cc         write (unit,260)
cc  260    format (/,5x,'Bin',9x,'Counts',7x,'Distance',7x,'Raw g(r)',
cc     &           4x,'Smooth g(r)',/)
cc         do i = 1, nbin
cc            write (unit,270)  i,hist(i,ctn),
cc     &         (dble(i)-0.5d0)*width,gr(i,ctn),gs(i,ctn)
cc  270    format (i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)
cc         hist(i,ctn) = 0
cc         gr(i,ctn) = 0.0d0
cc         gs(i,ctn) = 0.0d0
cc         end do 
cc         close(unit)
     
c
c     increment start and stop for next call to radialsub
c

      rdf_start = rdf_start + 1
      rdf_stop = rdf_stop + 1
      ctn = ctn + 1
      print *, 'ctn = ', ctn
      print *, 'rdf_mean = ', rdf_mean

      end   