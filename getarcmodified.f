c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getarc  --  get a coordinates archive file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getarc" asks for a coordinate archive or trajectory file
c     name, then reads the formatted or binary archive file
c
c
      subroutine getarcmodified (iarc)
      use files
      use inform
      use iounit
      use output
      implicit none
      integer iarc,iaux
      integer nask
      integer freeunit
      integer i
      logical exist,first
      character*1 letter
      character*240 arcfile
      character*240 auxfile
c
c
c     replace the extension of the archive file name with ".arc"
c
      arcfile = filename
      write (iout,10) arcfile
   10    format (/,' Coordinate Archive File Name :  ',a240)
      i = index(arcfile, '.', back=.true.)
      if (i > 0) then
         arcfile = arcfile(:i) // "arc"
      endif
      write (iout,20) arcfile
   20    format (/,' Coordinate Archive File Name :  ',a240)
      

c     open file the file and read the first set of coordinates
c
      filename = arcfile
      coordtype = 'CARTESIAN'
      iarc = freeunit ()
      open (unit=iarc,file=arcfile,status='old')
      rewind (unit=iarc)
      read (iarc,30)  letter
   30 format (a1)
      archive = .false.
      if (letter .eq. ' ')  archive = .true.
      if (letter.ge.'0' .and. letter.le.'9')  archive = .true.
      binary = (.not. archive)
c
c     read initial Cartesian coordinates from formatted file
c
      if (archive) then
         rewind (unit=iarc)
         call readxyz (iarc)
      end if
c
c     get atom types and connectivity from formatted file
c
      if (binary) then
         call nextarg (auxfile,exist)
         if (exist) then
            call basefile (auxfile)
            call suffix (auxfile,'xyz','old')
            inquire (file=auxfile,exist=exist)
         end if
         nask = 0
         do while (.not.exist .and. nask.lt.maxask)
            nask = nask + 1
            write (iout,40)
   40       format (/,' Enter Formatted Coordinate File Name :  ',$)
            read (input,50)  auxfile
   50       format (a240)
            call basefile (auxfile)
            call suffix (auxfile,'xyz','old')
            inquire (file=auxfile,exist=exist)
         end do
         if (.not. exist)  call fatal
         iaux = freeunit ()
         open (unit=iaux,file=auxfile,status='old')
         rewind (unit=iaux)
         call readxyz (iaux)
         close (unit=iaux)
c
c     read initial Cartesian coordinates from binary file
c
         filename = arcfile
         close (unit=iarc)
         iarc = freeunit ()
         open (unit=iarc,file=arcfile,form='unformatted',status='old')
         rewind (unit=iarc)
         first = .true.
         call readdcd (iarc,first)
      end if
c
c     quit if the coordinates archive file contains no atoms
c
      if (abort) then
         write (iout,60)
   60    format (/,' GETARC  --  Coordinate Archive File',
     &              ' was not Read Correctly')
         close (unit=iarc)
         call fatal
      end if
      return
      end
