c
c     RDF parameters module
c
      module rdfparams
      implicit none

      integer ctn
      integer NN
      integer sfac_nbins
      integer max_vect
      integer rdf_mean, rdf_nbin
      integer savelock
      integer rdf_num
      real*8  rdf_rmax, rdf_width
      real*8  q_min_temp, q_step 
      real*8  q_max_up, q_step_up
      logical rdf_intramol
      logical rdf_read_file, rdf_smooth 
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
      real*8, allocatable :: Sij(:,:)
      real*8, allocatable :: S(:)
      real*8, allocatable :: S_hist(:,:,:)
      real*8, allocatable :: q_vec(:,:)
      end module rdfparams