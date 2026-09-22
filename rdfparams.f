c
c     RDF parameters module
c
      module rdfparams
      implicit none

      integer NN
      integer sfac_nbins
      integer max_vect
      integer rdf_mean, rdf_nbin
      integer n_species
      real*8  rdf_rmax, rdf_width
      real*8  q_min_temp, q_step
      logical rdf_intramol
      logical rdf_read_file, rdf_smooth
      logical exp_present
      real*8 rdf_kappa

      logical debug_global
      logical verbose_global
      logical :: use_scatter = .false.
      logical :: scatter_expminusone = .true.
      logical :: scatter_autotune = .false.

      integer current_md_step

      character*250 :: scatter_tag
c
c     base name without the "_k<value>" suffix, used for diagnostic
c     files that don't depend on the restraint force constant k
c     (qmatrix_imp, qmatrix_red, st_init, debug_experimental_sfac);
c     scatter_tag keeps the _k suffix for the k-dependent outputs
c     (st.csv and energy.csv)
c
      character*250 :: scatter_basetag

c     smoothing factor and per-bin running average for the restraint's
c     exponentially-decaying memory function (see structfactor_forces
c     and the SCATTER-RESTRAIN memory-function comment in radialask.f)
      real*8 :: scatter_alpha
      logical :: ewma_initialized = .false.
      real*8, allocatable :: S_bin_ewma(:)

      real*8, allocatable :: q_vec(:,:)
      real*8, allocatable :: q_magnitude(:)
      integer, allocatable :: atype(:)
      real*8, allocatable :: mole_fractions(:)
      real*8, allocatable :: S_exp_binned(:)

c     restraint normalization and optional peak-region emphasis
c     (see structfactor_forces, and SCATTER-NORMALIZE / SCATTER-PEAKEMPH
c     in radialask.f). When scatter_normalize is on, the restraint energy
c     is divided by W = sum(w_i) (= number of contributing bins when all
c     weights are 1), so k becomes a penalty per unit mean-square
c     deviation, independent of bin count / q-range. Both default off, in
c     which case S_exp_weight is all 1.0 and the energy is the original
c     0.5*k*sum((S-S_exp)^2).
      logical :: scatter_normalize = .false.
      logical :: scatter_peakemph = .false.
      real*8  :: peak_qlo = 0.0d0
      real*8  :: peak_qhi = 0.0d0
      real*8  :: peak_factor = 1.0d0
      real*8, allocatable :: S_exp_weight(:)

c     weighted RMS deviation of the restrained running-average S(q) from
c     experiment, RMS = sqrt(sum_i w_i*(S-S_exp)^2 / sum_i w_i), updated
c     every step in structfactor_forces and logged to *_energy.csv by
c     gradient.f as the restraint convergence metric
      real*8 :: scatter_rms = 0.0d0
      end module rdfparams
