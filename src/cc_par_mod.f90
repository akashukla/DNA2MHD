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

  use, intrinsic :: iso_c_binding
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
  REAL :: init_energy=0.01
  REAL :: init_kolm = 1.66 !initial scaling exponent of the energy spectrum i.e. vk^2 propto k^(-init_kolm)
  REAL :: phdf = 0.25
  REAL :: phdfxy = 0.0
  INTEGER :: kxinit_min=0,kyinit_min=0,kzinit_min=0  !minimum k's  (i.e. sets box size)
  INTEGER :: kxinit_max=8,kyinit_max=16,kzinit_max=18  !maximum k's  (i.e. sets box size)
  INTEGER :: nkxforce=2,nkyforce=2,nkzforce=2  !ks to force
  REAL :: force_frac
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
  INTEGER(C_INTPTR_T) :: nx0_big=1, ny0_big=1, nz0_big=1
  INTEGER :: nh0=1,nspec=1
  !Inverse gradient scale lengths
  REAL :: omt=5.0,omn=5.0
  REAL :: Ti0Te=1.0 !Ratio of i to e temperature
  CHARACTER(len=200) :: diagdir='./'  !Output directory
  CHARACTER(len=200) :: loaddir='./'  !Output directory 
  REAL :: etg_factor=0.0  !Factor for flux surface averaged phi ==> 0.0 = ETG, 1.0 = ITG
  !left/right_vec:  Output of eignevectors for lapack calculation
  LOGICAL :: left_vec=.false.,right_vec=.false.
  !left_ev: SLEPc uses transpose of linear operator==> finds left eigenvectors
  LOGICAL :: left_ev=.false.
  CHARACTER(len=2) :: comp_type='IV'
  REAL(C_DOUBLE) :: dt_max=0.01    !initial maximum time step
  REAL :: turnover
  LOGICAL :: fix_dt 
  REAL :: courant=0.3   !courant factor times 2 pi for dt calculation
  LOGICAL :: ev_slepc=.true.
  !test_nl:  compares nonlinearity for pseudo spectral vs. convolution
  LOGICAL :: test_nl=.false.
  !k*max: for nonlinear grid--calculated from k*min and nk*0
  REAL :: kzmax
  REAL :: kxmax
  REAL :: kymax
  REAL :: kmax

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
  INTEGER :: rhs_lin_version = 12
  INTEGER :: rhs_nl_version !Akash changed nl_version to 1
  INTEGER :: intorder = 4
  LOGICAL :: linen = .false.
  LOGICAL :: keepzero = .true. ! Retain nl calculation for modes at start at zero energy
  INTEGER :: dealias_type = 3
  LOGICAL :: shifted = .true.
  LOGICAL :: splitx = .true.
  ! IO Numbers for Diss/NL term debugging
  LOGICAL :: plot_nls = .false.
  LOGICAL :: enone = .true.
  LOGICAL :: nv = .false.
  LOGICAL :: test_ho = .false.
  REAL :: hall = 1.0
  LOGICAL :: guide  = .true.
  LOGICAL :: uni = .true.
  LOGICAL :: shear = .false.
  LOGICAL :: beltrami = .false.
  LOGICAL :: helical = .false.
  LOGICAL :: walenp = .false.
  LOGICAL :: walenn = .false.
  LOGICAL :: mhc = .false.
  LOGICAL :: init_wave = .false.
  LOGICAL :: init_null = .false.
  LOGICAL :: force_trunc = .false.
  LOGICAL :: bc_norm = .false.
  LOGICAL :: track_divs = .true.
  LOGICAL :: debug_energy = .true.
  LOGICAL :: taylorgreen = .false.
  INTEGER :: random_state = 0
  
  REAL :: en_leftwhist = 1.0
  REAL :: en_leftcyclo = 0.0
  REAL :: en_rightwhist = 0.0
  REAL :: en_rightcyclo = 0.0

  REAL :: force_lw = 1.0
  REAL :: force_lc = 0.0
  REAL :: force_rw = 0.0
  REAL :: force_rc = 0.0
  
  LOGICAL :: calc_dt=.false.        !Automatic initial time step calculation
  LOGICAL :: dt_slepc=.false.        !Use slepc or lapack
!  LOGICAL :: adapt_dt_nl=.true.     !Adapt time step in nonlinear sims
  LOGICAL :: adapt_dt_nl=.false.     !FALSE for dna2mhd initially
  LOGICAL :: verbose=.false.        !Extra output
  LOGICAL :: timer = .false.        ! Timing specific output only
  LOGICAL :: checkpoint_read=.false. !READ checkpoint for restart
  LOGICAL :: checkpoint_write=.true.
!  LOGICAL :: get_chpt_from_gout=.false.

  INTEGER(4) :: init_cond = 0
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
  REAL, PARAMETER :: e=exp(1.0)
  COMPLEX, PARAMETER :: i_complex=(0.0,1.0)

  !MPI 
  INTEGER(kind=4) :: mype=0
  INTEGER(kind=4) :: n_mpi_procs
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
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: b_1
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: v_1
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: reader
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: reader2
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:)  :: fullsmallarray,fullbigarray
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: scatter_big(:),scatter_small(:),gather_big(:),gather_small(:)
  
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: LW,LC,RW,RC
  INTEGER(4), ALLOCATABLE, DIMENSION(:,:,:) :: paddingmask
  ! Forcing masks
  LOGICAL, ALLOCATABLE :: mask(:,:,:)
  REAL, ALLOCATABLE :: mask1(:,:,:)
  
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: LWp,LCp,RWp,RCp
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: LWp2,LCp2,RWp2,RCp2
  REAL :: last_reset = -1.0! time of last forcing phase change

  INTEGER :: rkstage
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: bdv,vdb,cbdb,bdcb,vdv,bdb,db2
  REAL, ALLOCATABLE, DIMENSION(:,:) :: kperp2
  REAL, ALLOCATABLE, DIMENSION(:) :: kxgrid,kygrid,kzgrid,herm_grid,hgrid_loc,&
                            hkgrid,vgrid,delta_hk, delta_v
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: phi
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: pre ! Pressure Fourier transform
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: gpsi ! Magnetic helicity gauge transformation added to Coulomb gauge vector potential

  REAL, ALLOCATABLE, DIMENSION(:,:) :: phi_denom
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kmags
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kperps
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kzs
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: alpha_leftwhist
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: alpha_leftcyclo
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: pcurleig
  
  REAL(C_DOUBLE) :: mhelcorr
  REAL :: precorr,preendiv
  REAL :: vdvcorr,vdvendiv

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
  REAL(C_DOUBLE) :: dt
  REAL :: max_time=1.0e15
  INTEGER(4) :: max_itime=1000,itime=0
  INTEGER(4) :: itime_start
  REAL :: max_walltime=83000.0
  REAL :: time_start = 0.0
  REAL :: time_tot = 10000.0
  LOGICAL :: nonlinear
  LOGICAL :: actual_nonlinear
  LOGICAL :: force_turbulence
  INTEGER :: forcetype
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
  INTEGER :: lkx_ind ! index of lowest kx value
  INTEGER :: lkz_ind !index of lowest (most negative) kx value
  INTEGER :: lkx_big !index of most negative kx value in big array
  INTEGER :: lky_big !index of lowest (most negative) ky value in big (for dealiasing) arrays
  INTEGER :: lkz_big !index of lowest (most negative) kx value in big (for dealiasing) arrays

  INTEGER :: io_number=200

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
  INTEGER :: s_handle,f_handle
  
  REAL :: vnu ! Viscosity in units ion skin depth * v_A; scaled to vnu/(maxk**(2*hyp))
  REAL :: eta ! Magnetic Prandtl number, mag diffusion constant eta * vnu
  REAL :: rey = 0 ! Reynolds Number; if rey > 0 set vnu = 1/rey
  INTEGER :: hyp
  

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_io_number                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_io_number
    IMPLICIT NONE

    io_number=io_number+1

  END SUBROUTINE get_io_number


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                clear_padding                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE clear_padding(temp_big,temp_small)

    IMPLICIT NONE

    complex(C_DOUBLE_COMPLEX), intent(in)  :: temp_big(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: temp_small(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2)
    integer(C_INTPTR_T) :: low_index, high_index,i

    temp_small = cmplx(0.0,0.0)

    if (n_mpi_procs.eq.1) then
      ! In the one process case hkz_ind and lkz_ind both lie between lkz1 and lkz2                                
      DO i = 0,nkx0-1

         temp_small(i,0:hky_ind,0:hkz_ind) = temp_big(i+1,1:1+hky_ind,1:hkz_ind+1)
         temp_small(i,0:hky_ind,lkz_big:nz0_big-1) = temp_big(i+1,1:1+hky_ind,1+lkz_big:nz0_big)
         temp_small(i,lky_big:ny0_big-1,lkz_big:nz0_big-1) = temp_big(i+1,1+lky_big:ny0_big,1+lkz_big:nz0_big)
         temp_small(i,lky_big:ny0_big-1,0:hkz_ind) = temp_big(i+1,1+lky_big:ny0_big,1:hkz_ind+1)

      ENDDO
   else
      ! Two cases, lkz1 in first region (and maybe lkz2 is), or lkz2 in last (and maybe lkz1 is)
      ! Both accounted for by choice of lbkz and ubkz
      DO i = 0,nkx0-1
         temp_small(i,0:hky_ind,lbkz:ubkz) = temp_big(i+1,1:1+hky_ind,lbkz+1:ubkz+1)
         temp_small(i,lky_big:ny0_big-1,lbkz:ubkz) = temp_big(i+1,1+lky_big:ny0_big,lbkz+1:ubkz+1)
      ENDDO
      
   endif

 end subroutine clear_padding
 
 SUBROUTINE clear_4d(array)

   IMPLICIT NONE

   complex(C_DOUBLE_COMPLEX), intent(inout) :: array(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)

   reader2 = array(:,:,:,0)
   CALL clear_padding(reader2,reader)
   array(:,:,:,0) = reader

   reader2 = array(:,:,:,1)
   CALL clear_padding(reader2,reader)
   array(:,:,:,1) = reader

   reader2 = array(:,:,:,2)
   CALL clear_padding(reader2,reader)
   array(:,:,:,2) = reader

 END SUBROUTINE clear_4d

 END MODULE par_mod
