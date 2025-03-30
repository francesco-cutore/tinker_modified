c
c     RDF parameters module
c
      module rdfparams
      implicit none

      integer rdf_start, rdf_stop, rdf_step, ctn
      integer rdf_mean, rdf_nbin
      integer savelock
      real*8  rdf_rmax, rdf_width
      logical rdf_intramol
      character*6 rdf_labelj, rdf_labelk   
      
      integer, allocatable :: hist(:,:)
      real*8, allocatable :: gr(:,:)
      real*8, allocatable :: gs(:,:)
      real*8, allocatable :: gr_mean(:)
      real*8, allocatable :: gs_mean(:)

      end module rdfparams