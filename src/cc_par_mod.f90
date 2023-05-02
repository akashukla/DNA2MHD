!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                               cc_par_mod.f90                              !!
!!                                                                           !!
!!  par_mod                                                                  !!
!!  -- get_io_number                                                         !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                  par_mod                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Module containing the definition of DNA's parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE par_mod
  IMPLICIT NONE

  !! Input parameters
  !!!!!!!!!!!!!!!!!!!!!
  CHARACTER(len=6) :: data_precision
  REAL :: nu  !collisionality
  REAL :: hyp_nu=-1.0  !nu for hyper collisions (by default set to nu)
  LOGICAL :: em_conserve=.false.         !implements an energy and momentum conserving operator (Parker)--hasn't been tested
  LOGICAL :: force_kz0eq0=.false.        !Deletes kz=0 modes at each time step
  LOGICAL :: force_ky0eq0=.false.        !Deletes ky=0 modes at each time step
  LOGICAL :: force_kx0eq0=.false.        !Deletes kx=0 modes at each time step
  LOGICAL :: forceb=.false. ! Whether or not b is forced
  REAL :: kzmin=0.1,kymin=0.1,kxmin=0.1  !minimum k's  (i.e. sets box size)
  REAL :: init_amp_bx=0.01,init_amp_by=0.01,init_amp_bz=0.01, init_amp_vx=0.01, init_amp_vy=0.01, init_amp_vz=0.01
  REAL :: init_kolm = 1.66 !initial scaling exponent of the energy spectrum i.e. vk^2 propto k^(-init_kolm)
  REAL :: phdf = 0.25
  REAL :: phdfxy = 0.0
  INTEGER :: kxinit_min=0,kyinit_min=0,kzinit_min=0  !minimum k's  (i.e. sets box size)
  INTEGER :: kxinit_max=8,kyinit_max=16,kzinit_max=18  !maximum k's  (i.e. sets box size)
  INTEGER :: nkxforce=2,nkyforce=2,nkzforce=2  !ks to force
  REAL :: force_amp= 0.5
  LOGICAL :: kmin_eq_0=.false.           
  INTEGER :: hypx_order=16                
  INTEGER :: hypy_order=16
  INTEGER :: hypz_order=16
  INTEGER :: hypv_order=6 !Parker suggestion
  !REAL :: hyp_z_factor=1.0          
  REAL :: hyp_zonal=0.0                  !coefficient for hyp on zonals (kz=0,ky=0)
  REAL :: hyp_conv=0.0                   !coefficient for hyp on kz=0,n=0
  LOGICAL :: hyp_conv_ky=.false.                   !coefficient for hyp on kz=0,n=0
  INTEGER :: num_k_hyp_conv=0               !number of k>0 modes to apply hyp_conv to (ky and kz)
  REAL :: hyp_x=0.02                     !coefficient for hyp on kx
  REAL :: hyp_y=0.02                     ! hyp on ky
  REAL :: hyp_z=0.0                      !hyp on v
  REAL :: hyp_v=0.0 !Parker suggestion: 10.0
  LOGICAL :: nuno_closure=.true.
  !Number of modes for each coord.
  INTEGER :: nkx0=1,nky0=1,nkz0=1,nv0=1 
  INTEGER :: nh0=1,nspec=1
  !Inverse gradient scale lengths
  REAL :: omt=5.0,omn=5.0
  REAL :: Ti0Te=1.0 !Ratio of i to e temperature
  CHARACTER(len=200) :: diagdir='./'  !Output directory
  REAL :: etg_factor=0.0  !Factor for flux surface averaged phi ==> 0.0 = ETG, 1.0 = ITG
  !left/right_vec:  Output of eignevectors for lapack calculation
  LOGICAL :: left_vec=.false.,right_vec=.false.
  !left_ev: SLEPc uses transpose of linear operator==> finds left eigenvectors
  LOGICAL :: left_ev=.false.
  CHARACTER(len=2) :: comp_type='IV'
  REAL :: dt_max=0.00    !initial maximum time step
  REAL :: courant=0.3   !courant factor times 2 pi for dt calculation
  LOGICAL :: ev_slepc=.true.
  !test_nl:  compares nonlinearity for pseudo spectral vs. convolution
  LOGICAL :: test_nl=.false.
  !k*max: for nonlinear grid--calculated from k*min and nk*0
  REAL :: kzmax
  REAL :: kxmax
  REAL :: kymax
  !k*max: for eigenvalue scans (DEFAULT is set high in CASE of misuse in EV
  !calcs
  REAL :: kzmax0=1000.0
  REAL :: kxmax0=1000.0
  REAL :: kymax0=1000.0

  !FLR effects
  LOGICAL :: mu_integrated=.true.
  LOGICAL :: flr_on=.true.
  INTEGER :: flr_version=1
  LOGICAL :: flr_extra=.true.
  LOGICAL :: flr_nonlinear=.false.

  !HANKEL effects
  LOGICAL :: hankel=.false.
  REAL :: hkmax
  !When hankel is off, then the coordinated used is mu (but it is written a v)
  REAL :: vmax = 9   

  !Variables that dynamic procedure needs
  REAL ::  fracx = 0.5
  REAL ::  fracy = 0.5
  INTEGER :: istep_GyroLES = 50 
  LOGICAL :: GyroLES = .false.
  LOGICAL :: Gyroherm = .false.
  LOGICAL :: Gyroz = .false.
  LOGICAL :: Corr = .false.
  REAL, ALLOCATABLE, DIMENSION(:) :: hyp_x_herm1    ! hyp on kx with hermit dependence
  REAL, ALLOCATABLE, DIMENSION(:) :: hyp_y_herm1    ! hyp on ky with hermit dependence
  REAL :: hyp_x_herm =0.001    ! hyp on kx with hermit dependence
  REAL:: hyp_y_herm =0.001   ! hyp on ky with hermit dependence
  !hyp_x_herm(:) = 0.001
  !hyp_y_herm(:) = 0.001
 
  !performance
  LOGICAL :: performance_tests=.true. !! Main performance switch    
  LOGICAL :: perf_test_lin=.false.
  LOGICAL :: perf_test_nl=.false.
  LOGICAL :: perf_test_rhs=.false.
  LOGICAL :: perf_test_par=.false.
  INTEGER :: perf_monitor(2)
  !Note:rhs_lin_version=2 is broken with hyp_x/y (need to debug)!!!
  INTEGER :: rhs_lin_version
  !INTEGER :: rhs_nl_version=2
  INTEGER :: rhs_nl_version !Akash changed nl_version to 1
  INTEGER :: dealias_type = 3
  ! IO Numbers for Diss/NL term debugging
  LOGICAL :: plot_nls = .false.
  LOGICAL :: enone = .true.
  LOGICAL :: nv = .false.
  INTEGER :: dbio = 200 
  INTEGER :: dvio = 201
  INTEGER :: bdvio = 202
  INTEGER :: vdbio = 203
  INTEGER :: bdcbio = 204
  INTEGER :: cbdbio = 205
  INTEGER :: vdvio = 206
  INTEGER :: bdbio = 207
  INTEGER :: db2io = 208
  REAL :: hall = 1.0
  LOGICAL :: guide  = .true.

  LOGICAL :: calc_dt=.false.        !Automatic initial time step calculation
  LOGICAL :: dt_slepc=.false.        !Use slepc or lapack
!  LOGICAL :: adapt_dt_nl=.true.     !Adapt time step in nonlinear sims
  LOGICAL :: adapt_dt_nl=.false.     !FALSE for dna2mhd initially
  LOGICAL :: verbose=.false.        !Extra output
  LOGICAL :: checkpoint_read=.false. !READ checkpoint for restart
  LOGICAL :: checkpoint_write=.true.
!  LOGICAL :: get_chpt_from_gout=.false.

  CHARACTER(len=40) :: init_cond='DEFAULT'
  REAL :: init_prefactor=0.001
  INTEGER :: version_flag

  LOGICAL :: kscan=.false.
  !np_herm:  number of v (total) processors
  INTEGER :: np_kz=1   
  INTEGER :: np_herm=1   
  INTEGER :: np_hank=1   
  INTEGER :: np_spec=1   
  INTEGER :: np_total=1   

  !More
  INTEGER :: lxyzvhs0
  LOGICAL :: evenyz
  LOGICAL :: spatial2D
  !SLEPc
  INTEGER :: ev_size_loc,ev_size,n_ev!,ev_n_test
  REAL :: ev_prec=1.0e-14
  INTEGER :: ev_max_it=10000
  CHARACTER(len=100) :: which_ev='jd'
  COMPLEX :: ev_shift=(10.0,0.0)
  !INTEGER :: ksp_max_it=0
  !CHARACTER(len=100) :: pc_type='DEFAULT'
  !INTEGER :: it_ev=-1
    
  !Diagnostics
  LOGICAL :: continue_run=.true.
  INTEGER :: istep_gamma=0
  INTEGER :: istep_gout=0
  LOGICAL :: gout_nl = .false.  ! output the nonlinear as well as the distribution function
  LOGICAL :: gout_2xt = .false. !output two adjecent time steps
  INTEGER :: istep_real=20
  INTEGER :: istep_hermite=0
  INTEGER :: istep_energy=0
  INTEGER :: istep_energyspec=0
  INTEGER :: istep_ffm=0       !istep for fields, fluxes, and moments
  INTEGER :: istep_energy3d=0  !for REAL quantities, i.e. flux, energy, etc.
  INTEGER :: istep_eshells=0
  INTEGER :: istep_fmom3d=0    !istep for fields and moments
  INTEGER :: istep_nlt_triple=0    !istep for Bogdan's triple transfer
  LOGICAL :: nlt_symmetrize = .true.  !symmetrize the triple transfer function for q,p
  
  INTEGER :: istep_schpt=1000  !istep for security checkpoint
  INTEGER :: istep_nltest=100
  INTEGER :: istep_gk=10
  INTEGER :: istep_gknl=10
  LOGICAL :: dt_output=.false.
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: phi_last
  COMPLEX :: last_frequency
  REAL :: eps_converge=1.0e-4
  INTEGER :: gk_ky_index=5
  INTEGER :: gk_kz_index=1
  INTEGER :: gk_kx_index=0

  !Permanent
  REAL, PARAMETER :: pi=3.141592653589793238
  REAL, PARAMETER :: e=2.71828183
  COMPLEX, PARAMETER :: i_complex=(0.0,1.0)

  !MPI 
  INTEGER :: mype=0
  INTEGER :: n_mpi_procs
  INTEGER :: mpi_comm_cart_4d
  INTEGER :: mpi_comm_kz
  INTEGER :: mpi_comm_herm
  INTEGER :: mpi_comm_hank
  INTEGER :: mpi_comm_spec
  INTEGER :: mype_kz=0
  INTEGER :: mype_herm=0
  INTEGER :: mype_hank=0
  INTEGER :: mype_spec=0


  !Arrays and Matrices

  !distribution function
  !COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: g_1
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: b_1
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: v_1
  REAL, ALLOCATABLE, DIMENSION(:,:) :: kperp2
  REAL, ALLOCATABLE, DIMENSION(:) :: kxgrid,kygrid,kzgrid,herm_grid,hgrid_loc,&
                            hkgrid,vgrid,delta_hk, delta_v
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: phi
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: pre ! Pressure Fourier transform
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: gpsi ! Magnetic helicity gauge transformation added to Coulomb gauge vector potential

  REAL, ALLOCATABLE, DIMENSION(:,:) :: phi_denom
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kmags

  !for parallelization
  !Hermites
  INTEGER :: lv0
  INTEGER :: lv1,lv2,lbv,ubv
  !kz
  INTEGER :: lkz0
  INTEGER :: lkz1,lkz2,lbkz,ubkz
  !Hankel
  INTEGER :: lh0
  INTEGER :: lh1,lh2,lbh,ubh,nhb
  !Species
  INTEGER :: ls0
  INTEGER :: ls1,ls2,lbs,ubs
  
  !For initial_value
  REAL :: time=0.0  
  REAL :: dt
  REAL :: max_time=1.0e15
  INTEGER :: max_itime=1000000,itime=0 
  INTEGER :: itime_start
  REAL :: max_walltime=83000.0
  LOGICAL :: nonlinear
  LOGICAL :: actual_nonlinear
  LOGICAL :: force_turbulence
  LOGICAL :: set_forcing = .true.
  !INTEGER :: which_nonlinear=1
  LOGICAL :: linear_nlbox !Grid, initialization, etc. is the same, but don't USE nonlinearity

  !For eigensolve
  COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: M_loc, loc_lvec, loc_rvec 
  COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: loc_q
  COMPLEX, ALLOCATABLE, DIMENSION(:) :: eig


  !For Hankel transformation
  REAL, ALLOCATABLE, DIMENSION(:,:) :: T_hkv 
  REAL, ALLOCATABLE, DIMENSION(:) ::  m1_v
  REAL, ALLOCATABLE, DIMENSION(:) :: m2_hk
  REAL, ALLOCATABLE, DIMENSION(:) :: f_in

  !n=global
  !l=local
  !r=row
  !c=column
  INTEGER :: nr_M, nc_M
  INTEGER :: lr_M, lc_M
  !Processor block sizes
  INTEGER :: pr_M,pc_M
  !For time step adaptation
  LOGICAL :: first_stage

  !useful indices
  INTEGER :: hkx_ind !index of highest kx value (for dealiasing)
  INTEGER :: hky_ind !index of highest ky value
  INTEGER :: lky_ind !index of lowest (most negative) ky value
  INTEGER :: hkz_ind !index of highest kx value
  INTEGER :: lkz_ind !index of lowest (most negative) kx value
  INTEGER :: lky_big !index of lowest (most negative) ky value in big (for dealiasing) arrays
  INTEGER :: lkz_big !index of lowest (most negative) kx value in big (for dealiasing) arrays

  INTEGER :: io_number=100

  !NLT diagnostics
  !input parameters
  INTEGER :: istep_nlt = 0
  !INTEGER :: istep_nlt_full = 0
  INTEGER :: nlt_version = 1
  !min_shell_width is the factor that multiplies delta k for the last shell
  REAL :: min_shell_width=1.0
  LOGICAL :: output_nlt_n=.false.

  !rk tests
  REAL :: rc0,ic0
  COMPLEX :: c0
  REAL :: dt_rktest=1.0
  REAL :: rmax,imax
  REAL :: delta_lambda
  LOGICAL :: test_rk4=.false.

  ! for checkpoints in io
  INTEGER :: b_out_handle, v_out_handle
  REAL :: vnu ! Viscosity in units ion skin depth * v_A is vnu/(maxk**(2*hyp))
  REAL :: eta ! Magnetic Prandtl number, mag diffusion constant eta * vnu
  INTEGER :: hyp
  

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_io_number                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_io_number
    IMPLICIT NONE

    io_number=io_number+1

  END SUBROUTINE get_io_number


END MODULE par_mod
