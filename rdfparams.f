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
      logical exp_present
      real*8 rdf_kappa

      logical debug_global
      logical verbose_global

      integer current_md_step

      real*8, allocatable :: S_hist(:,:,:)
      real*8, allocatable :: q_vec(:,:)
      real*8, allocatable :: q_magnitude(:)
      integer, allocatable :: atype(:)
      real*8, allocatable :: mole_fractions(:)
      real*8, allocatable :: S_exp_binned(:)
      end module rdfparams
