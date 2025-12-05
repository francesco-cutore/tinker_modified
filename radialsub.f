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
c     pair of atom types via analysis of a set of coordinate frames and
c     then calculates the structure factor for the pair of atom types.
c
      subroutine radialsub (istep)
      use atomid
      use atoms
      use bound
      use boxes
      use limits
      use math
      use molcul
      use potent
      use rdfparams
      use factors
      implicit none
      integer i,j,k,f,nq
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
c
c     for strfactor calculation
c
c      real(8) Fa, Fb, a, b
      real(8) kkk, gr_local, r ,s_ab
      real(8) kk,d,cut,dens
      real(8) fi,fj
      real(8) xii,xjj
      character*2 temporary(2)
      real(8) temporarynumber(2)
      real(8) somaquad(600)
c      real(8) kernel
      real(8) delta_ab
      real(8) :: dr, q, sinc, w, integrand
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
      if(savelock .eq. 1000) then !impossible
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
      end if

      if (ctn .eq. rdf_mean) then
         savelock = 1
      endif
c
c     open the output file only to write the mean
c
         if (savelock .eq. 1) then
               dens = n / xbox**3
               S = 0.0d0
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
c
c         write the mean to the output file
c
c            if(savelock .eq. 1000) then !impossible
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
c            end if
         

         cut = rdf_rmax
         d = rdf_width
c
c     calculate the atomic fractions please allocate 1 1 1 9 bf29bq 
c
         xii = REAL(numj) / REAL(n)
         xjj = REAL(numk) / REAL(n)

      Sij(:,f) = 0.0d0
      print *, 'xii', xii
      print *, 'xjj', xjj
      print *, int(cut/d)

      ! Kronecker delta
      if (namek.eq.namej .or. 
     &          typek.eq.typej) then
      delta_ab = 1.0d0
      else
      delta_ab = 0.0d0
      end if

      ! Step size in r
      dr = ( xbox / 2.0d0) / int(cut/d)

      ! Loop over q values from 0.01 to 6.00
      do nq = 1, NN
      q = 0.01d0 * real(nq, 8)
      s_ab = xii * delta_ab

      ! Integration over r
      do j = 1, int(cut/d)
         r = (j - 0.5d0) * dr
         gr_local = gr_mean(j, f)

         ! sinc(qr)
         if (q * r == 0.0d0) then
            sinc = 1.0d0
         else
            sinc = sin(q * r) / (q * r)
         end if

         ! Lorch window
         if (r == 0.0d0) then
            w = 1.0d0
         else
c            w = 1.0d0
            w = sin(2.0d0 * pi * r / cut) / (2.0d0 * pi * r / cut)
         end if

         ! Integrand
         integrand = 4.0d0 * pi * r**2 * (gr_local - 1.0d0) * sinc * w
         s_ab = s_ab + xii * xjj * dens * integrand * dr
      end do

      Sij(nq,f) = s_ab
c      print *, 'q = ', q, '  S_ab(q) = ', s_ab
      end do
   
          do i=1, NN
          kkk = 0.01 * i
           if (kkk.gt.6) EXIT
c        compute structure factor factors as funct of q
            fi =  atomF(namej,kkk)
c            print *, 'fi' , fi 
            fj =  atomF(namek,kkk)
c            print *, 'fj' , fj
            S(i) = S(i) + Sij(i,f)
c            S(i) = S(i) + Sij(i,f)*fj*fi*xii*xjj
c            print *, 'sTRFACT', S(i)
          end do 
      end if         !finish mean savelock     
      end do

      if (savelock .eq. 1) then
          temporary(1:2) = ['H', 'O']
          temporarynumber = [0.666667,0.333333]
          somaquad = 0

          unit = freeunit ()
          open(unit,file='st.txt',status='unknown',position='append')
          do i=1,NN
            kk = .01 * i
            do j = 1, 2
               somaquad(i) = somaquad(i) + 
     &            atomF(temporary(j),kk) * temporarynumber(j)
c               print *, 'atomf=', atomF(temporary(j),kk)
             end do
             somaquad(i) = somaquad(i)**2
c             print *, somaquad(i)
           write(unit,250) 0.01*i, S(i)
c            write(unit,250) 0.01*i, dens*S(i)/somaquad(i)
  250    format (3x,f12.4,3x,f12.4)
          end do
         close(unit)
      endif   
      ctn = ctn + 1
      end   