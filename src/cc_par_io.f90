!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                                cc_par_io.f90                              !!
!!                                                                           !!
!! read_parameters                                                           !!
!! output_parameters                                                         !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                read_par                                   !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_parameters
  USE par_mod
  !USE GaussQuadrature, only: mu_grid_type
  USE mpi
  IMPLICIT NONE

  INTEGER :: ierr, par_handle

  NAMELIST /diagnostics/ &
      diagdir, istep_ffm, istep_energy3d,istep_fmom3d, istep_gamma,istep_nltest,&
      istep_schpt,dt_output,istep_energy,istep_hermite,istep_gout,istep_real,&
      istep_nlt,min_shell_width,istep_eshells,output_nlt_n,&
      istep_gk,istep_gknl,gk_ky_index,gk_kz_index,gk_kx_index,istep_GyroLES, &
      istep_nlt_triple, nlt_symmetrize, gout_nl, gout_2xt

  NAMELIST /numerical_parameters/ &
      kzmin,kymin,kxmin,nkx0,nky0,nkz0,nv0,nh0,nspec,&
      hyp_v,hyp_x,hyp_y,hyp_z,hypv_order,np_herm,&
      np_kz,np_spec,np_hank,&
      courant,hyp_conv,num_k_hyp_conv,hyp_conv_ky,hypx_order,hypy_order,hypz_order,&
      kxinit_min, kxinit_max, kyinit_min, kyinit_max, kzinit_min, kzinit_max,&
      init_amp_bx,init_amp_by,init_amp_bz, init_amp_vx, init_amp_vy, init_amp_vz,hyp
  !,&
      !mu_grid_type, vmax,hyp_nu,fracx, fracy

  NAMELIST /physical_parameters/ &
      nu,omt,omn,Ti0Te,eta,vnu

  NAMELIST /flags/ &
      nonlinear, actual_nonlinear, test_nl, calc_dt, comp_type,adapt_dt_nl,&
      linear_nlbox,verbose,checkpoint_read,checkpoint_write,&
      em_conserve,flr_on,force_kz0eq0,force_ky0eq0,force_kx0eq0,flr_version,&
      flr_extra,flr_nonlinear,etg_factor, &!, which_nonlinear,etg_factor
      perf_test_lin,perf_test_nl,perf_test_rhs,rhs_lin_version,rhs_nl_version,&
      perf_test_par, version_flag, hankel, dt_slepc, nuno_closure,mu_integrated,&
      GyroLES, Gyroherm, Gyroz, Corr, &
      plot_nls,dbio,dvio,bdvio,vdbio,bdcbio,cbdbio,vdvio,bdbio,db2io,&
      hall
 
  NAMELIST /eigensolve/ &
      left_vec,right_vec,ev_slepc, kxmax0, kymax0, kzmax0, kscan,n_ev,&
      left_ev, ev_prec, ev_max_it, ev_shift

  NAMELIST /initial_value/ &
      dt_max,  max_itime, max_time ,init_cond,init_prefactor,max_walltime

  NAMELIST /rk_test/ &
      rc0,ic0,dt_rktest,rmax,imax,delta_lambda,test_rk4

  IF (mype==0) WRITE(*,*) "Reading parameters."
  IF (mype==0) WRITE(*,*)

  IF(mype==0) THEN
     WRITE(*,*) "Reading input parameters."
  END IF
  
  CALL get_io_number
  par_handle=io_number

  OPEN(unit = par_handle, file = 'parameters', status = 'unknown')

  READ(par_handle, nml = diagnostics, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect diagnostics NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = numerical_parameters, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect numerical_parameters NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = physical_parameters, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect physical_parameters NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = flags, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect flags NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = eigensolve, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect eigensolve NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = initial_value, iostat = ierr)
  IF (ierr.ne.0) STOP 'on i/o error: incorrect initial_value NAMELIST'

  REWIND(par_handle)
  READ(par_handle, nml = rk_test, iostat = ierr)
  IF (ierr.gt.0) STOP 'on i/o error: incorrect rk_test NAMELIST'

  CLOSE(par_handle)


  !! Initialization and checks
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(mype==0) THEN
     WRITE(*,*) "Initialization and checks for input parameters."
  END IF

  !Should update according to svn revision number when significant changes are
  !made
  version_flag=53

  np_total=np_kz*np_herm*np_hank*np_spec

  !Tau=1.0/Ti0Te
  IF(hyp_nu==-1.0) hyp_nu=nu

  IF(calc_dt) THEN
    !IF(mype==0) WRITE(*,*) "!!!!!!!!calc_dt not working at the moment!!!!!!!!!!"
    !STOP
  END IF

  IF(np_total.ne.n_mpi_procs) STOP "wrong number of processors!"

  IF(nh0==1) THEN
     mu_integrated=.true.
  ELSE 
     mu_integrated=.false.
     rhs_nl_version = 1
  END IF

  IF (.not.GyroLES) Gyroherm = .false. 
  IF (Gyroz==.true.) THEN
    Gyroherm = .false. 
    GyroLES = .false.
  ENDIF

  IF(MOD(nv0,np_herm).ne.0) STOP "nv0 must be divisible by np_herm!"
  IF(nv0/np_herm.lt.2) STOP "nv0 must be greater than 2*np_herm!"

  IF(MOD(nkz0,np_kz).ne.0) STOP "nkz0 must be divisible by np_kz!"

  IF(MOD(nh0,np_hank).ne.0) STOP "nh0 must be divisible by np_hank!"

  IF(MOD(nspec,np_spec).ne.0) STOP "nspec must be divisible by np_spec!"
  IF((nspec.ne.1).and.(nspec.ne.2)) STOP "nspec must by one or two!"

  IF(MOD(nkx0,2).ne.0.and.nonlinear) STOP "nkx0 must be even!"
  IF(MOD(nky0,2).ne.0.and.nonlinear) STOP "nky0 must be even!"
  IF(MOD(nkz0,2).ne.0.and.nonlinear.and.nkz0.ne.1) STOP "nkz0 must be even!"
  IF(nkz0==1.and.nonlinear) THEN
      rhs_nl_version=4
      init_cond='SINGLE_KZ'
      spatial2d= .true.
  END IF

  IF((hyp_x.ne.0.0).or. &
    (hyp_y.ne.0.0).or.  &
    (hyp_z.ne.0.0)) THEN
    IF(comp_type=='EV'.and.mype==0) WRITE(*,*) "Warning, k*max0's must be set &
    for accurate hyp_xyz EV calculations."
  END IF

  IF(nonlinear) THEN
    IF(comp_type=='EV') STOP "Cannot DO nonlinear EV calculation."
    evenyz=.true.
    kmin_eq_0=.true.
    IF(nkx0.le.2.or.nky0.le.2) STOP "Must have more than one Fourier mode &
           & for nonlinear."
  ELSE
    kmin_eq_0=.false.
    evenyz=.false.
    !IF(nkx0.ne.1.or.nky0.ne.1.or.nkz0.ne.1) THEN
    !  IF(mype==0) WRITE(*,*) "Linear run ==> setting nk's to 1"
    !  nkx0=1
    !  nky0=1
    !  nkz0=1
    !END IF
  END IF

  IF(nonlinear) THEN
    kxmax0=(nkx0-1)*kxmin
    kymax0=(nky0/2-1)*kymin !max not nyquist
    kzmax0=(nkz0/2-1)*kzmin !max not nyquist
  END IF

  IF(test_nl.and..not.nonlinear) STOP "Error! Must USE nonlinear=T for test_nl=T."

  IF(dt_max==0.0.and..not.calc_dt) STOP "Must define dt_max or set calc_dt=T."

  IF(linear_nlbox) adapt_dt_nl=.false.

  IF(etg_factor==0) THEN
    IF(mype==0) WRITE(*,*) "ETG run."
  ELSE IF (etg_factor==1) THEN
    IF(mype==0) WRITE(*,*) "ITG run."
  ELSE
    IF(mype==0) WRITE(*,*) "Using fraction of ITG zonal flow response:",etg_factor
  END IF

  IF(nonlinear.and.left_ev) STOP "Cannot DO nonlinear with left_ev."
  IF(left_ev.and.comp_type=='IV') STOP "Cannot DO comp_type='IV' with left_ev."
  IF(nonlinear.and..not.linear_nlbox.and.istep_gamma.ne.0) THEN
    istep_gamma=0
  END IF

  IF(flr_nonlinear) STOP "flr_nonlinear not allowed at this time!"

  IF(rhs_lin_version==2) THEN
      rhs_lin_version=1
      IF(mype==0) WRITE(*,*) "rhs_lin_version changed to 1!!!!"
      IF(mype==0) WRITE(*,*) "rhs_lin_version=2 is broken for hyp_xy.  Need to debug!"
  END IF

  !IF(istep_energy.gt.1.and.etg_factor.ne.0) STOP "Must have etg_factor=0.0 &
          !for energy diagnostics at this time."

  IF(mype==0) WRITE(*,*) "Done with read_par."

  !CALL output_parameters

END SUBROUTINE read_parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                             output_parameters                             !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Outputting parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE output_parameters
  USE par_mod
  INTEGER :: out_handle

  IF (verbose) WRITE(*,*) "In Output Parameters"


  CALL get_io_number
  out_handle=io_number
  IF (verbose) WRITE(*,*) "Got io number, outhandle is: ", out_handle

  IF(mype==0)  THEN
    IF (verbose) WRITE(*,*) "Mype was 0"
  OPEN(unit=out_handle,file=trim(diagdir)//'/parameters.dat',status='unknown')
    IF (verbose) WRITE(*,*) "Opened Parameters"
  
    IF (verbose) WRITE(*,*) "Physical Parameters"

    !****** physical_parameters NAMELIST ********
    WRITE(out_handle,"(A)")    "&physical_parameters"
    WRITE(out_handle,"(A,G12.4)") "nu = ",nu
    WRITE(out_handle,"(A,G12.4)") "omt = ",omt
    WRITE(out_handle,"(A,G12.4)") "omn = ",omn
    WRITE(out_handle,"(A,G12.4)") "Ti0Te = ",Ti0Te
    WRITE(out_handle,"(A,G12.4)") "vnu = ",vnu
    WRITE(out_handle,"(A,G12.4)") "eta = ",eta
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""

    IF (verbose) WRITE(*,*) "Numerical Parameters"

    !****** numerical_parameters NAMELIST ********
    WRITE(out_handle,"(A)")    "&numerical_parameters"
    WRITE(out_handle,"(A,G12.4)") "kxmin = ",kxmin
    WRITE(out_handle,"(A,G12.4)") "kymin = ",kymin
    WRITE(out_handle,"(A,G12.4)") "kzmin = ",kzmin
    WRITE(out_handle,"(A,I4)") "nkx0 = ",nkx0    
    WRITE(out_handle,"(A,I4)") "nky0 = ",nky0    
    WRITE(out_handle,"(A,I4)") "nkz0 = ",nkz0    
    WRITE(out_handle,"(A,I4)") "nv0 = ",nv0    
    WRITE(out_handle,"(A,I4)") "nh0 = ",nh0    
    WRITE(out_handle,"(A,I4)") "nspec = ",nspec    
    WRITE(out_handle,"(A,G12.4)") "hyp_x = ",hyp_x
    WRITE(out_handle,"(A,G12.4)") "hyp_y = ",hyp_y
    WRITE(out_handle,"(A,G12.4)") "hyp_z = ",hyp_z
    WRITE(out_handle,"(A,I4)") "hypx_order = ",hypx_order
    WRITE(out_handle,"(A,I4)") "hypy_order = ",hypy_order
    WRITE(out_handle,"(A,I4)") "hypz_order = ",hypz_order
    WRITE(out_handle,"(A,G12.4)") "hyp_v = ",hyp_v
    WRITE(out_handle,"(A,I4)") "hypv_order = ",hypv_order
    !WRITE(out_handle,"(A,G12.4)") "hyp_z_factor = ",hyp_z_factor
    !WRITE(out_handle,"(A,G12.4)") "hyp_zonal = ",hyp_zonal
    WRITE(out_handle,"(A,G12.4)") "hyp_conv = ",hyp_conv
    IF(hyp_conv.gt.0.0) WRITE(out_handle,"(A,I4)") "num_k_hyp_conv = ",num_k_hyp_conv
    IF(hyp_conv.gt.0.0) WRITE(out_handle,"(A,L1)") "hyp_conv_ky = ",hyp_conv_ky
    WRITE(out_handle,"(A,I4)") "np_herm = ",np_herm    
    WRITE(out_handle,"(A,I4)") "np_kz = ",np_kz    
    WRITE(out_handle,"(A,I4)") "np_hank = ",np_hank    
    WRITE(out_handle,"(A,I4)") "np_spec = ",np_spec    
    !WRITE(out_handle,"(A,L1)") "kmin_eq_0  = ", kmin_eq_0
    WRITE(out_handle,"(A,G12.4)") "courant = ",courant
    WRITE(out_handle,"(A,G12.4)") "hyp_nu = ",hyp_nu
    IF (.not.mu_integrated.and..not.hankel) WRITE(out_handle,"(A,G12.4)") "vmax = ",vmax
    IF (GyroLES.or.Gyroz.or.Corr) WRITE(out_handle,"(A,G12.4)") "fracx = ",fracx
    IF (GyroLES.or.Corr) WRITE(out_handle,"(A,G12.4)") "fracy = ",fracy
    WRITE(out_handle,"(A,I4)") "hyp = ",hyp
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""

    IF (verbose) WRITE(*,*) "Namelist Parameters"
    !****** diagnostics NAMELIST ********
    WRITE(out_handle,"(A)")    "&diagnostics"
    WRITE(out_handle,"(3A)")   "diagdir = '", TRIM(diagdir),"'"
    WRITE(out_handle,"(A,I4)") "istep_ffm = ",istep_ffm    
    WRITE(out_handle,"(A,I4)") "istep_energy3d = ",istep_energy3d    
    WRITE(out_handle,"(A,I4)") "istep_energy = ",istep_energy    
    WRITE(out_handle,"(A,I4)") "istep_hermite = ",istep_hermite    
    WRITE(out_handle,"(A,I4)") "istep_gout = ",istep_gout    
    if(gout_nl) then
       WRITE(out_handle,"(A,L1)") "gout_nl  = ", gout_nl
    end if
    if( gout_2xt) then
       WRITE(out_handle,"(A,L1)") "gout_2xt  = ", gout_2xt
    end if
    WRITE(out_handle,"(A,I4)") "istep_gk = ",istep_gk    
    WRITE(out_handle,"(A,I4)") "istep_gknl = ",istep_gknl    
    WRITE(out_handle,"(A,I4)") "istep_nlt = ",istep_nlt    
    WRITE(out_handle,"(A,I4)") "istep_nlt_triple = ",istep_nlt_triple    
    IF(istep_nlt_triple.gt.0) WRITE(out_handle,"(A,L1)") "nlt_symmetrize = ",nlt_symmetrize    
    WRITE(out_handle,"(A,I4)") "istep_eshells = ",istep_eshells    
    IF(istep_real.ne.20) WRITE(out_handle,"(A,I4)") "istep_real = ",istep_real    
    WRITE(out_handle,"(A,I4)") "istep_fmom3d = ",istep_fmom3d    
    WRITE(out_handle,"(A,I4)") "istep_gamma = ",istep_gamma    
    IF(istep_nltest.gt.0) WRITE(out_handle,"(A,I4)") "istep_nltest = ",istep_nltest    
    WRITE(out_handle,"(A,I4)") "istep_schpt = ",istep_schpt    
    IF (GyroLES==.true..or.Gyroz) WRITE(out_handle,"(A,I4)") "istep_GyroLES = ",istep_GyroLES    
    IF(dt_output) WRITE(out_handle,"(A,L1)") "dt_output  = ", dt_output
    IF(output_nlt_n) WRITE(out_handle,"(A,L1)") "output_nlt_n  = ", output_nlt_n
    IF(istep_eshells.gt.0) WRITE(out_handle,"(A,G12.4)") "min_shell_width = ", min_shell_width
    WRITE(out_handle,"(A,I4)") "gk_ky_index = ",gk_ky_index    
    WRITE(out_handle,"(A,I4)") "gk_kz_index = ",gk_kz_index    
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""

    IF (verbose) WRITE(*,*) "Flags Parameters"

    !****** flags NAMELIST ********
    WRITE(out_handle,"(A)")    "&flags"
    WRITE(out_handle,"(A,L1)") "nonlinear  = ", nonlinear
    WRITE(out_handle,"(A,L1)") "actual_nonlinear  = ", actual_nonlinear
    !WRITE(out_handle,"(A,I4)") "which_nonlinear = ",which_nonlinear
    IF(test_nl) WRITE(out_handle,"(A,L1)") "test_nl  = ", test_nl
    WRITE(out_handle,"(A,L1)") "calc_dt  = ", calc_dt
    IF(dt_slepc) WRITE(out_handle,"(A,L1)") "dt_slepc  = ", dt_slepc
    IF(nuno_closure) WRITE(out_handle,"(A,L1)") "nuno_closure  = ", nuno_closure
    !WRITE(out_handle,"(A,L1)") "gamma1  = ", gamma1
    WRITE(out_handle,"(3A)")   "comp_type = '",trim(comp_type),"'"
    IF(.not.adapt_dt_nl) WRITE(out_handle,"(A,L1)") "adapt_dt_nl  = ", adapt_dt_nl
    IF(linear_nlbox) WRITE(out_handle,"(A,L1)") "linear_nlbox  = ", linear_nlbox
    WRITE(out_handle,"(A,L1)") "checkpoint_read  = ", checkpoint_read
    IF(.not.checkpoint_write) WRITE(out_handle,"(A,L1)") "checkpoint_write  = ", checkpoint_write
    WRITE(out_handle,"(A,G12.4)") "etg_factor  = ", etg_factor
    IF(verbose) WRITE(out_handle,"(A,L1)") "verbose  = ", verbose
    IF(.not.flr_on) WRITE(out_handle,"(A,L1)") "flr_on = ",flr_on
    WRITE(out_handle,"(A,L1)") "mu_integrated = ",mu_integrated
    WRITE(out_handle,"(A,L1)") "hankel = ",hankel
    WRITE(out_handle,"(A,L1)") "GyroLES = ",GyroLES
    IF (GyroLES) WRITE(out_handle,"(A,L1)") "Gyroherm = ",Gyroherm
    IF (Gyroz==.true.) WRITE(out_handle,"(A,L1)") "Gyroz = ",Gyroz
    IF (Corr==.true.) WRITE(out_handle,"(A,L1)") "Corr = ",Corr
    IF(flr_version.ne.1) WRITE(out_handle,"(A,I4)") "flr_version = ",flr_version
    IF(.not.flr_extra) WRITE(out_handle,"(A,L1)") "flr_extra = ",flr_extra
    !WRITE(out_handle,"(A,L1)") "flr_nonlinear = ",flr_nonlinear
    IF(perf_test_lin) WRITE(out_handle,"(A,L1)") "perf_test_lin = ",perf_test_lin
    IF(perf_test_par) WRITE(out_handle,"(A,L1)") "perf_test_par = ",perf_test_par
    IF(perf_test_nl) WRITE(out_handle,"(A,L1)") "perf_test_nl = ",perf_test_nl
    IF(perf_test_rhs) WRITE(out_handle,"(A,L1)") "perf_test_rhs = ",perf_test_rhs
    WRITE(out_handle,"(A,I4)") "version_flag = ",version_flag
    !IF(istep_nlt.ne.0) WRITE(out_handle,"(A,I4)") "nlt_version = ",nlt_version
    WRITE(out_handle,"(A,I4)") "rhs_lin_version = ",rhs_lin_version
    WRITE(out_handle,"(A,I4)") "rhs_nl_version = ",rhs_nl_version
    WRITE(out_handle,"(A,L1)") "em_conserve = ",em_conserve
    IF(force_kz0eq0) WRITE(out_handle,"(A,L1)") "force_kz0eq0 = ",force_kz0eq0
    IF(force_ky0eq0) WRITE(out_handle,"(A,L1)") "force_ky0eq0 = ",force_ky0eq0
    IF(force_kx0eq0) WRITE(out_handle,"(A,L1)") "force_kx0eq0 = ",force_kx0eq0
    WRITE(out_handle,"(A,L1)") "plot_nls = ",plot_nls
    IF (plot_nls) THEN 
    WRITE(out_handle,"(A,I4)") "dbio = ",dbio
    WRITE(out_handle,"(A,I4)") "dvio = ",dvio
    WRITE(out_handle,"(A,I4)") "bdvio = ",bdvio
    WRITE(out_handle,"(A,I4)") "vdbio = ",vdbio
    WRITE(out_handle,"(A,I4)") "bdcbio = ",bdcbio
    WRITE(out_handle,"(A,I4)") "cbdbio = ",cbdbio
    WRITE(out_handle,"(A,I4)") "vdvio = ",vdvio
    WRITE(out_handle,"(A,I4)") "bdbio = ",bdbio
    WRITE(out_handle,"(A,I4)") "db2io = ",db2io
    ENDIF
    WRITE(out_handle,"(A,I4)") "hall = ",hall
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""


    IF (verbose) WRITE(*,*) "Eigensolve "

    !****** eigensolve NAMELIST ********
    WRITE(out_handle,"(A)")    "&eigensolve"
    WRITE(out_handle,"(A,I4)") "n_ev = ",n_ev    
    IF(left_ev) WRITE(out_handle,"(A,L1)") "left_ev  = ", left_ev
    IF(kscan) WRITE(out_handle,"(A,L1)") "kscan  = ", kscan
    IF(ev_slepc) WRITE(out_handle,"(A,L1)") "ev_slepc  = ", ev_slepc
    IF(left_vec) WRITE(out_handle,"(A,L1)") "left_vec  = ", left_vec
    IF(right_vec) WRITE(out_handle,"(A,L1)") "right_vec  = ", right_vec
    WRITE(out_handle,"(A,G12.4)") "kxmax0 = ",kxmax0
    WRITE(out_handle,"(A,G12.4)") "kymax0 = ",kymax0
    WRITE(out_handle,"(A,G12.4)") "kzmax0 = ",kzmax0
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""

    IF (verbose) WRITE(*,*) "Initial Value"

    !****** initial_value NAMELIST ********
    WRITE(out_handle,"(A)")    "&initial_value"
    WRITE(out_handle,"(A,I10)") "max_itime = ",max_itime    
    WRITE(out_handle,"(A,G12.4)") "max_walltime = ",max_walltime    
    WRITE(out_handle,"(A,G12.4)") "dt_max = ",dt_max
    WRITE(out_handle,"(A,G12.4)") "max_time = ",max_time
    !WRITE(out_handle,"(A)") "init_cond = '",init_cond,"'"
    IF(trim(init_cond).ne.'DEFAULT') WRITE(out_handle,"(3A)") "init_cond = '", TRIM(init_cond),"'"
    IF(init_prefactor.ne.0.001) WRITE(out_handle,"(A,G12.4)") "init_prefactor =", init_prefactor
    WRITE(out_handle,"(A)")    "/"
    WRITE(out_handle,"(A)")    ""


    IF(test_rk4) THEN
        IF (verbose) WRITE(*,*) "In test rk4"
      !****** rk_test NAMELIST ********
      WRITE(out_handle,"(A)")    "&rk_test"
      WRITE(out_handle,"(A,L1)") "test_rk4  = ", test_rk4
      WRITE(out_handle,"(A,G12.4)") "rc0 = ",rc0
      WRITE(out_handle,"(A,G12.4)") "ic0 = ",ic0
      WRITE(out_handle,"(A,G12.4)") "delta_lambda = ",delta_lambda
      WRITE(out_handle,"(A,G12.4)") "rmax = ",rmax
      WRITE(out_handle,"(A,G12.4)") "imax = ",imax
      WRITE(out_handle,"(A,G12.4)") "dt_rktest = ",dt_rktest
      WRITE(out_handle,"(A)")    "/"
      WRITE(out_handle,"(A)")    ""
    END IF

   CLOSE(out_handle)
   IF (verbose) WRITE(*,*) "closed out handle"

   END IF !mype==0

END SUBROUTINE output_parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                              checkpoint_out                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Writing checkpoints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE checkpoint_out(purpose)
  USE par_mod
  !USE mpi
  !USE field_solver, only: get_phi
  USE nonlinearity, only: get_rhs_nl
  IMPLICIT NONE

  INTEGER, INTENT(in) :: purpose
  CHARACTER(len=100) :: chp_name, chp_name_b, chp_name_v
  INTEGER :: chp_handle_b, chp_handle_v, chp_handle
  INTEGER :: l,p
  !COMPLEX :: g_out(0:nkx0-1,0:nky0-1,0:nkz0-1,0:lv0-1)
  COMPLEX :: b_out(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2)
  COMPLEX :: v_out(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2)
  !INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: send_proc,recv_proc,ierr
  !LOGICAL :: g_output,not_first
  LOGICAL :: not_first,b_output,v_output

   IF(np_hank.gt.1) STOP "checkpoint_out not yet implemented for np_hank.gt.1"
   IF(np_spec.gt.1) STOP "checkpoint_out not yet implemented for np_spec.gt.1"
   IF(np_kz.gt.1) STOP "checkpoint_out not yet implemented for np_kz.gt.1"

  not_first=.false.
  IF(purpose==1) THEN !WRITE security checkpoint
    chp_name='/s_checkpoint'
    !g_output=.false.
    b_output=.false.
    v_output=.false.
  ELSE IF(purpose==2) THEN !WRITE final checkpoint
    chp_name='/checkpoint'
      !g_output=.false.
      b_output=.false.
      v_output=.false.
  ELSE IF(purpose==3) THEN !WRITE out distribution function-start new file
      chp_name_b='/b_out.dat'
      chp_name_v='/v_out.dat'
      !g_output=.true.
      b_output=.true.
      v_output=.true.
  ELSE IF(purpose==4) THEN !WRITE out distribution function to possibly existing g_out file
      chp_name_b='/b_out.dat'
      chp_name_v='/v_out.dat'
      !g_output=.true.
      b_output=.true.
      v_output=.true.
      INQUIRE(file=trim(diagdir)//trim(chp_name_b),exist=not_first)
      INQUIRE(file=trim(diagdir)//trim(chp_name_v),exist=not_first)
  !ELSE IF(purpose==5) THEN !WRITE out nonlinearity-start new file
  !    chp_name='/gnl_out.dat'
  !    !g_output=.true.
  !    b_output=.true.
  !    v_output=.true.
  !ELSE IF(purpose==6) THEN !WRITE out nonlinearity to possibly existing g_out file
  !    chp_name='/gnl_out.dat'
  !    !g_output=.true.
  !    b_output=.true.
  !    v_output=.true.
  !    INQUIRE(file=trim(diagdir)//trim(chp_name),exist=not_first)
  END IF

  IF(v_output.and.b_output.and.not_first) THEN
    chp_handle_b=b_out_handle
  ELSE
    CALL get_io_number
    chp_handle=io_number
    IF(b_output) b_out_handle=io_number
    IF(v_output) v_out_handle=io_number
  END IF

  IF(not_first) THEN
    IF(mype==0) THEN
      IF(b_output.and.v_output) THEN
        OPEN(unit=chp_handle_b,file=trim(diagdir)//trim(chp_name_b), form='unformatted', status='unknown',access='stream',position='append')
        OPEN(unit=chp_handle_v,file=trim(diagdir)//trim(chp_name_v), form='unformatted', status='unknown',access='stream',position='append')
      ELSE 
        WRITE(*,*) 'In this Section,chp_handle is: ', chp_handle, 'chp_name is: ', chp_name
        OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name), form='unformatted', status='REPLACE',access='stream')
      END IF
    END If
  ELSE
    IF(mype==0) THEN
      IF(b_output.and.v_output) THEN
        OPEN(unit=chp_handle_b,file=trim(diagdir)//trim(chp_name_b), form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=chp_handle_v,file=trim(diagdir)//trim(chp_name_v), form='unformatted', status='REPLACE',access='stream')
      ELSE 
        OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name), form='unformatted', status='REPLACE',access='stream')
      END IF
    END IF
  END IF

  IF(mype==0) THEN
    IF(.not.b_output) THEN
        WRITE(*,*) 'Trying to write to chp_handle'
	  !Output info necessary for restarts using the checkpoint
	  WRITE(chp_handle) itime 
  	  WRITE(chp_handle) dt 
	  WRITE(chp_handle) nkx0 
	  WRITE(chp_handle) nky0 
	  WRITE(chp_handle) nkz0 
      WRITE(*,*) 'Wrote partially'
  	  !WRITE(chp_handle) nv0 
      WRITE(chp_handle) time 
      WRITE(*,*) 'Wrote fully'
    END IF
    IF(b_output.and.v_output) THEN
      WRITE(*,*) 'Trying to write to chp_handle_b and v'
      WRITE(chp_handle_b) time 
      WRITE(chp_handle_v) time 
      WRITE(*,*) 'Wrote to chp_handle_b and v'
    END IF
  END IF

  IF(purpose==5.or.purpose==6) THEN
  !   IF(verbose) WRITE(*,*) "gnlout:p,mype",p,mype
  !   CALL get_phi(g_1) 
  !   g_out=cmplx(0.0,0.0)
  !   CALL get_rhs_nl(g_1,phi,g_out)
  ELSE
     !g_out = g_1(:,:,:,:,0,0)
     b_out = b_1(:,:,:,:)
     v_out = v_1(:,:,:,:)
  ENDIF

  IF(mype==0) THEN
    !DO l=lv1,lv2
    DO l=0,2
      WRITE(*,*) 'l=0,2 loop'
	  !WRITE(chp_handle) g_out(:,:,:,l)
      IF(b_output.and.v_output) THEN
	    WRITE(chp_handle_b) b_out(:,:,:,l)
	    WRITE(chp_handle_v) v_out(:,:,:,l)
      !ELSE
	  !  WRITE(chp_handle) b_out(:,:,:,l)
      END IF
	  !WRITE(chp_handle) v_out(:,:,:,l)
      WRITE(*,*) 'Done with l=0,2'
    END DO
  END IF 

!  DO p=1,np_herm-1
!    send_proc=p
!    recv_proc=0
!    IF(verbose) WRITE(*,*) "p,mype",p,mype
!
!    !IF(mype==send_proc) CALL MPI_Send(g_1(0,0,0,lv1,0,0), nkx0*nky0*nkz0*lv0, &
!    !IF(mype==send_proc) CALL MPI_Send(g_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!	!  	   	      MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!    IF(mype==send_proc) CALL MPI_Send(b_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!	   	   	      MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!    IF(mype==send_proc) CALL MPI_Send(v_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!	   	   	      MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!
!    !IF(mype==recv_proc) CALL MPI_Recv(g_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!	!			  MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!    IF(mype==recv_proc) CALL MPI_Recv(b_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!				  MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!    IF(mype==recv_proc) CALL MPI_Recv(v_out(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!				  MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!    IF(mype==0) THEN
!      DO l=0,lv0-1
!	    WRITE(chp_handle) g_out(:,:,:,l)
!      END DO
!    END IF
!  END DO

  IF(b_output.and.v_output) THEN
      WRITe(*,*) 'Closing v and b handles'
    IF(mype==0) CLOSE(chp_handle_b)
    IF(mype==0) CLOSE(chp_handle_v)
  ELSE
      WRITe(*,*) 'Closing chp handle'
    IF(mype==0) CLOSE(chp_handle)
  END IF
  IF(verbose) WRITE(*,*) "checkpoint_out,mype",mype

  !CALL mpi_barrier(mpi_comm_world,ierr)

END SUBROUTINE checkpoint_out



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                             checkpoint_in                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  For reading checkpoints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SUBROUTINE checkpoint_in
!!      USE par_mod
!!      USE mpi
!!      IMPLICIT NONE
!!        
!!      CHARACTER(len=100) :: chp_name
!!      INTEGER :: chp_handle
!!      INTEGER :: l,p
!!      INTEGER :: stat(MPI_STATUS_SIZE)
!!      INTEGER :: send_proc,recv_proc,ierr
!!      INTEGER :: nkx0_in,nky0_in,nkz0_in,nv0_in
!!      COMPLEX :: g_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0:lv0-1)
!!    
!!       IF(np_hank.gt.1) STOP "checkpoint_in not yet implemented for np_hank.gt.1"
!!       IF(np_spec.gt.1) STOP "checkpoint_in not yet implemented for np_spec.gt.1"
!!       IF(np_kz.gt.1) STOP "checkpoint_in not yet implemented for np_kz.gt.1"
!!    
!!      chp_name='/checkpoint'
!!      
!!      CALL get_io_number
!!      chp_handle=io_number
!!    
!!      ! me
!!      IF(mype==0) WRITE(*,*) "Set the chp handle."
!!    
!!      IF(mype==0) OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name),form='unformatted',status='unknown',access='stream')
!!    
!!      IF(mype==0) THEN
!!        READ(chp_handle) itime 
!!        READ(chp_handle) dt 
!!        READ(chp_handle) nkx0_in 
!!        READ(chp_handle) nky0_in
!!        READ(chp_handle) nkz0_in 
!!        READ(chp_handle) nv0_in 
!!        READ(chp_handle) time 
!!        IF(nkx0_in.ne.nkx0) THEN
!!      	  IF(mype==0) WRITE(*,*) "Error in checkpoint_in: incorrect nkx0",nkx0,nkx0_in
!!    	  STOP
!!        END IF
!!        IF(nky0_in.ne.nky0) STOP "Error in checkpoint_in: incorrect nky0"
!!        IF(nkz0_in.ne.nkz0) STOP "Error in checkpoint_in: incorrect nkz0"
!!        IF(nv0_in.ne.nv0) STOP "Error in checkpoint_in: incorrect nv0"
!!        WRITE(*,*) "Starting from checkpoint with:"
!!        WRITE(*,*) "itime=",itime
!!        WRITE(*,*) "time=",time
!!        WRITE(*,*) "dt=",dt
!!        WRITE(*,*) "dt_max=",dt_max
!!      END IF
!!    
!!      !Send time info to other processors
!!      CALL MPI_BCAST(itime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
!!      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!      CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!    
!!      itime_start=itime
!!    
!!      IF(mype==0) THEN
!!        DO l=lv1,lv2
!!      	!READ(chp_handle) g_1(:,:,:,l,0,0)
!!      	READ(chp_handle) b_1(:,:,:,:)
!!      	READ(chp_handle) v_1(:,:,:,:)
!!        END DO
!!      END IF 
!!    
!!      DO p=1,np_herm-1
!!        IF(mype==0) THEN
!!          DO l=0,lv0-1
!!    	    READ(chp_handle) g_in(:,:,:,l)
!!          END DO
!!        END IF
!!    
!!        send_proc=0
!!        recv_proc=p
!!    
!!        !IF(mype==send_proc) CALL MPI_Send(g_in(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!!        !             MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!        IF(mype==send_proc) CALL MPI_Send(b_in(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!!                     MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!        IF(mype==send_proc) CALL MPI_Send(v_in(0,0,0,0), nkx0*nky0*nkz0*lv0, &
!!                     MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!        !IF(mype==recv_proc) CALL MPI_Recv(g_1(0,0,0,lv1,0,0), nkx0*nky0*nkz0*lv0, &
!!        !		     MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!        IF(mype==recv_proc) CALL MPI_Recv(b_1(0,0,0,lv1,0,0), nkx0*nky0*nkz0*lv0, &
!!        		     MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!    
!!        IF(mype==recv_proc) CALL MPI_Recv(v_1(0,0,0,lv1,0,0), nkx0*nky0*nkz0*lv0, &
!!        		     MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!      END DO
!!    
!!      CLOSE(chp_handle)
!!    
!!      CALL mpi_barrier(mpi_comm_world,ierr)
!!      WRITE(*,*) "Done reading checkpoint.",mype
!!    
!!    END SUBROUTINE checkpoint_in
