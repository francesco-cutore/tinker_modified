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
      integer n_species
      real*8  rdf_rmax, rdf_width
      real*8  q_min_temp, q_step 
      logical rdf_intramol
      logical rdf_read_file, rdf_smooth
      logical rdf_sort

      real*8, allocatable :: S_hist(:,:,:)
      real*8, allocatable :: q_vec(:,:)
      real*8, allocatable :: q_magnitude(:)
c      real*8, allocatable :: atom_factors(:,:)
      integer, allocatable :: atype(:)
      real*8, allocatable :: mole_fractions(:)
      end module rdfparams




c      character*6, allocatable :: rdf_labelj(:)
c      character*6, allocatable :: rdf_labelk(:)
c      character*3, allocatable :: rdf_namej(:)
c      character*3, allocatable :: rdf_namek(:) 
c      integer, allocatable :: rdf_typej(:)
c      integer, allocatable :: rdf_typek(:)      
c      integer, allocatable :: hist(:,:,:)
c      real*8, allocatable :: gr(:,:,:)
c      real*8, allocatable :: gs(:,:,:)
c      real*8, allocatable :: gr_mean(:,:)
c      real*8, allocatable :: gs_mean(:,:)
c      real*8, allocatable :: Sij(:,:)
c      real*8, allocatable :: S(:)