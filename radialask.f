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
c     "radial" finds the radial distribution function for a specified
c     pair of atom types via analysis of a set of coordinate frames
c
c
      subroutine radialask
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
      integer next
      integer start,stop
      integer step
      integer counter
      integer mean
      integer nbin
      integer i
      real*8 rmax,width
      logical exist,query
      logical intramol
      logical read_file
      character*1 answer
      integer freeunit
      integer unit
      character*240 record
      character*240 string
c
c     get numbers of the coordinate frames to be processed
c
      start = 1
      stop = 1
      step = 1
      counter = 1
      savelock = 0
      rdf_num = 1
      read_file = .true.
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  stop
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  step
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter mean and num   ',
     &              'of rdfs :  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  mean, rdf_num
   40    continue
      end if      
c
c     get the names of the atoms to be used in rdf computation
c
      allocate (rdf_labelj(rdf_num))
      allocate (rdf_labelk(rdf_num))
      allocate (rdf_namej(rdf_num))
      allocate (rdf_namek(rdf_num))
      allocate (rdf_typej(rdf_num))
      allocate (rdf_typek(rdf_num))  
c
c     set 0 the arrays
c
      rdf_labelj = '      '
      rdf_labelk = '      '
      rdf_namej  = '   '
      rdf_namek  = '   '
      rdf_typej  = -1
      rdf_typek  = -1

      do i=1, rdf_num
         write (iout,50)
   50    format (/,' Enter 1st & 2nd Atom Names or Type Numbers :  ',$)
         read (input,60)  record
   60    format (a240)
         next = 1
         call gettext (record,rdf_labelj(i),next)
         call gettext (record,rdf_labelk(i),next)
      end do
c
c     convert the labels to either atom names or type numbers
c
      do i=1, rdf_num
         read (rdf_labelj(i),
     &        *,err=70,end=70)  rdf_typej(i)
   70    continue
         if (rdf_typej(i) .le. 0) then
         next = 1
         call gettext (rdf_labelj(i),
     &         rdf_namej(i),next)
         end if
         read (rdf_labelk(i),
     &        *,err=80,end=80)  rdf_typek(i)
   80    continue
         if (rdf_typek(i) .le. 0) then
         next = 1
         call gettext (rdf_labelk(i),
     &        rdf_namek(i),next)
         end if
      end do      
c
c     get maximum distance from input or minimum image convention
c
      if (.not. use_bounds) then
         rmax = -1.0d0
        query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=90,end=90)  rmax
            query = .false.
         end if
   90    continue
         if (query) then
            write (iout,100)
  100       format (/,' Enter Maximum Distance to Accumulate',
     &                 ' [10.0 Ang] :  ',$)
            read (input,110)  rmax
  110       format (f20.0)
         end if
         if (rmax .le. 0.0d0)  rmax = 10.0d0
      else if (octahedron) then
         rmax = (sqrt(3.0d0)/4.0d0) * xbox
         rmax = 0.95d0 * rmax
      else
         rmax = min(xbox2*beta_sin*gamma_sin,ybox2*gamma_sin,
     &                         zbox2*beta_sin)
         rmax = 0.95d0 * rmax
      end if
c
c     get the desired width of the radial distance bins
c
      width = -1.0d0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=120,end=120)  width
         query = .false.
      end if
  120 continue
      if (query) then
         write (iout,130)
  130    format (/,' Enter Width of Distance Bins [0.01 Ang] :  ',$)
         read (input,140)  width
  140    format (f20.0)
      end if
      if (width .le. 0.0d0)  width = 0.01d0
c
c     decide whether to restrict to intermolecular atom pairs
c
      intramol = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,150)
  150    format (/,' Include Intramolecular Pairs in Distribution',
     &              ' [N] :  ',$)
         read (input,160)  record
  160    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  intramol = .true.
c
c     overwrite rdf.txt file
c
      unit = freeunit ()
      open(unit,file='rdf.txt',status='unknown')
      write(unit,170)
  170 format ('on the fly rdf calculation') 
      close(unit)     
c
c     set the number of distance bins to be accumulated
c
      nbin = int(rmax/width)   
c
c     store values in the module for later use by radialsub
c
      rdf_rmax = rmax
      rdf_width = width
      rdf_intramol = intramol
      ctn = counter
      rdf_mean = mean
      rdf_nbin = nbin

c
c     allocate hist, gr, and gs arrays
c
        allocate(hist(nbin,rdf_mean,rdf_num))
        allocate(gr(nbin,rdf_mean,rdf_num))
        allocate(gs(nbin,rdf_mean,rdf_num))
        allocate(gr_mean(nbin,rdf_num))
        allocate(gs_mean(nbin,rdf_num))
      
      hist = 0
      gr = 0.0d0
      gs = 0.0d0
      gr_mean = 0.0d0
      gs_mean = 0.0d0
      
      end
