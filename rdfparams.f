c
c     RDF parameters module
c
      module rdfparams
      implicit none

      integer rdf_start, rdf_stop, rdf_step, ctn
      integer rdf_mean, rdf_nbin
      integer savelock
      integer rdf_num
      real*8  rdf_rmax, rdf_width
      logical rdf_intramol
      logical rdf_read_file
c      character*6 rdf_labelj, rdf_labelk 
      character*6, allocatable :: rdf_labelj(:)
      character*6, allocatable :: rdf_labelk(:)
      character*3, allocatable :: rdf_namej(:)
      character*3, allocatable :: rdf_namek(:) 
      integer, allocatable :: rdf_typej(:)
      integer, allocatable :: rdf_typek(:)      
      integer, allocatable :: hist(:,:,:)
      real*8, allocatable :: gr(:,:,:)
      real*8, allocatable :: gs(:,:,:)
      real*8, allocatable :: gr_mean(:,:)
      real*8, allocatable :: gs_mean(:,:)

      end module rdfparams