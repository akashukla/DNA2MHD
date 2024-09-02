!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                                cc_init.f90                                !!
!!                                                                           !!
!!  init_run                                                                 !!
!!  arrays                                                                   !!
!!  finalize_arrays                                                          !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                    init_run                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Initiating the current run of the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_run
  !USE calculate_time_step, ONLY: calc_initial_dt, initialize_adapt_dt
  USE diagnostics, ONLY: initialize_diagnostics
  USE nonlinearity, ONLY: initialize_fourier  
  USE par_mod
  !USE hk_effects
  !USE GaussQuadrature
  IMPLICIT NONE
  
  CALL arrays
  
  !! Initial value computation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) "Starting initial value computation.",mype
  !IF(calc_dt.and.nonlinear) THEN
  !    CALL calc_initial_dt
  !    CALL output_parameters
  !END IF
  !IF(calc_dt.and.dt_slepc.and.nonlinear) THEN 
  !  !calculates time step by solving for ev's of every kx,ky,kz
  !  IF(mype==0) WRITE(*,*) "Calculation dt_max with SLEPc"
  !  CALL calc_dt_slepc
  !  CALL output_parameters
  !END IF

  !IF(calc_dt.and..not.dt_slepc.and.nonlinear) THEN 
  !  !calculates time step by solving for ev's of every kx,ky,kz
  !  !CALL calc_dt_lapack
  !  IF(mype==0) WRITE(*,*) "Calculation dt_max with LAPACK"
  !  CALL calc_dt_lapack
  !  CALL output_parameters
  !END IF

  itime=0
  time=0.0
  
  !! Initialised run
  !!!!!!!!!!!!!!!!!!!!
  CALL initial_condition(init_cond)   !checkpoint must be READ before others
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) "Called initial condition.",mype
  CALL initialize_diagnostics
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) "Called initial diagnostics.",mype
  !IF(nonlinear) CALL initialize_adapt_dt
  IF(nonlinear) CALL initialize_fourier
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) "Called initial fourier.",mype
  
END SUBROUTINE init_run



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                   arrays                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Sets up indices and allocates important arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE arrays
  USE par_mod
  !USE flr_effects
  !USE hk_effects
  !USE Gaussquadrature
  IMPLICIT NONE
  
  INTEGER :: i,j,l,hypv_handle

  !Info for k grids
  !Note the format for the ky/kx indices:
  !Index:	0,    1,      2,. . . .       hky_ind,  nky0/2,  lky_ind, . . . nky0-1
  !ky val:	0.0,  kymin,  2*kymin . . . . kymax,    *,       -kymax, . . .  -kymin
  if (splitx) then
    hkx_ind = nkx0-1
    hkz_ind = nkz0/2-1
    lkz_ind = nkz0-hkz_ind
  else
    hkx_ind = nkx0/2-1
    lkx_ind = nkx0-hkx_ind
    hkz_ind = nkz0-1
  endif

  hky_ind=nky0/2-1     !Index of maximum (used,i.e. not dummy) ky value 
  lky_ind=nky0-hky_ind !Index of minimum (most negative) ky value
!  IF(nkz0==1) THEN
!    hkz_ind=0
!    lkz_ind=0
!  ELSE
!    hkz_ind=nkz0-1 !Index of maximum (used,i.e. not dummy) kz value
!    lkz_ind=0 !Index of minimum (most negative) kz value
!  END IF
  

  IF(verbose.and.(mype.eq.0)) THEN
    WRITE(*,*) "hkx_ind",hkx_ind
    WRITE(*,*) "hky_ind",hky_ind
    WRITE(*,*) "lky_ind",lky_ind
    WRITE(*,*) "hkz_ind",hkz_ind
    WRITE(*,*) "lkz_ind",lkz_ind
  END IF

  !Info for v grid and parallelization
  lv0=nv0/np_herm   !Number of hermites per processor
  lv1=mype_herm*lv0          !Lower index on processor mype
  lv2=(mype_herm+1)*lv0-1    !Upper index on processor mype
  lbv=lv1-1             !Lower boundary index 
  IF(mype_herm==0) lbv=0     
  ubv=lv2+1             !Upper boundary index
  IF(mype_herm==np_herm-1) ubv=lv2  

  !Info for kz grid
  lkz0=nkz0
  lkz1=mype/n_mpi_procs*nkz0
  lkz2=(mype+1)/n_mpi_procs*nkz0-1
  !Verify the following when implementing kz parallelization
  lbkz=lkz1
  ubkz=lkz2

  !Info for Hankel grid
  lh0=nh0/np_hank
  lh1=mype_hank*lh0
  lh2=(mype_hank+1)*lh0-1
  !As a first implementation, no collision operator, so nhb = 0
  nhb  = 0
  lbh=lh1 - nhb
  ubh=lh2 + nhb


  !Info for species grid
  ls0=nspec/np_spec
  ls1=mype_spec*ls0
  ls2=(mype_spec+1)*ls0-1
  !Verify the following when implementing s parallelization
  lbs=ls1
  ubs=ls2

  lxyzvhs0=nkx0*nky0*lkz0*lv0*lh0*ls0   
  ev_size_loc=lxyzvhs0
  ev_size=lxyzvhs0*np_herm*np_hank*np_kz*np_spec
  
!  IF(.not.allocated(g_1))&
!      ALLOCATE(g_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)) 
  IF(.not.allocated(b_1))&
      ALLOCATE(b_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(v_1))&
      ALLOCATE(v_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(kxgrid)) ALLOCATE(kxgrid(0:nkx0-1))
  IF(.not.allocated(kygrid)) ALLOCATE(kygrid(0:nky0-1))
  IF(spatial2d) THEN
    !kzgrid proportional to ky grid (see Watanabe and Sugama '04)
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nky0-1))
  ELSE
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nkz0-1))
 END IF
 
  IF(.not.allocated(kmags)) ALLOCATE(kmags(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  IF(.not.allocated(kperps)) ALLOCATE(kperps(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  IF(.not.allocated(kzs)) ALLOCATE(kzs(0:nkx0-1,0:nky0-1,lkz1:lkz2))
 ! IF(.not.allocated(gpsi)) ALLOCATE(gpsi(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ! IF(.not.allocated(pre)) ALLOCATE(pre(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  IF(.not.allocated(pcurleig)) ALLOCATE(pcurleig(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  
  !IF(.not.allocated(herm_grid)) ALLOCATE(herm_grid(0:nv0-1))
  !IF(.not.allocated(hgrid_loc)) ALLOCATE(hgrid_loc(lv1:lv2))
  !IF(.not.allocated(hkgrid)) ALLOCATE(hkgrid(0:nh0-1))
  !IF(.not.allocated(delta_hk)) ALLOCATE(delta_hk(0:nh0-1))
  !IF(.not.allocated(vgrid)) ALLOCATE(vgrid(0:nh0-1))
  !IF(.not.allocated(delta_v)) ALLOCATE(delta_v(0:nh0-1))
  !IF(.not.allocated(kperp2)) ALLOCATE(kperp2(0:nkx0-1,0:nky0-1))
  !IF(.not.allocated(hyp_x_herm1)) ALLOCATE(hyp_x_herm1(lv1:lv2)) 
  !IF(.not.allocated(hyp_y_herm1)) ALLOCATE(hyp_y_herm1(lv1:lv2)) 

  !hyp_x_herm1(:) = 0.001 
  !hyp_y_herm1(:) = 0.001 

  IF(kmin_eq_0) THEN
    if (splitx) then

    DO i=0,nkx0-1
      kxgrid(i)=i*kxmin
    END DO
    kxmax=(nkx0-1)*kxmin

    kzgrid(0)=0.0
    DO i=1,nkz0/2-1
      kzgrid(i)=i*kzmin
      kzgrid(nkz0-i)=-i*kzmin
    END DO
    kzmax=(nkz0/2-1)*kzmin
    kzgrid(nkz0/2)=kzmax+kzmin

    else

    DO i=0,nkz0-1
      kzgrid(i)=i*kzmin
    END DO 
    kzmax=(nkz0-1)*kzmin

    kxgrid(0)=0.0
    DO i=1,nkx0/2-1
      kxgrid(i)=i*kxmin
      kxgrid(nkx0-i)=-i*kxmin
    END DO
    kxmax=(nkx0/2-1)*kxmin
    kxgrid(nkx0/2)=kxmax+kxmin 

    endif

    kygrid(0)=0.0
    DO i=1,nky0/2-1
      kygrid(i)=i*kymin
      kygrid(nky0-i)=-i*kymin
    END DO 
    kymax=(nky0/2-1)*kymin
    kygrid(nky0/2)=kymax+kymin

  ELSE !.not.kmin_eq_0, i.e. linear
    kxgrid(0)=kxmin
    kxmax=kxmin
    kygrid(0)=kymin
    kymax=kymin
    kzgrid(0)=kzmin
    kzmax=kzmin
 END IF

 kmax = sqrt(kxmax**2.0 + kymax**2.0 + kzmax**2.0)

  !DO i=0,nv0-1
  !  herm_grid(i)=REAL(i)
  !END DO 

  !DO i=lv1,lv2
  !  hgrid_loc(i)=REAL(i)
  !END DO 

  !DO i=0,nkx0-1
  !  DO j=0,nky0-1
  !    kperp2(i,j) = kxgrid(i)**2 + kygrid(j)**2
  !  ENDDO
  !ENDDO

  !Hankel grid

!  IF (.not.mu_integrated) THEN
!     IF (hankel) THEN
!        IF(.not.allocated(T_hkv)) ALLOCATE(T_hkv(0:nh0-1,0:nh0-1))
!        IF(.not.allocated(m1_v)) ALLOCATE(m1_v(0:nh0-1))
!        IF(.not.allocated(m2_hk)) ALLOCATE(m2_hk(0:nh0-1))
!        !IF(.not.allocated(f_in)) ALLOCATE(f_in(lh1:lh2))
!        hkmax = sqrt(2*maxval(kperp2))
!        call Gethankelgrid(hkmax,vmax,delta_hk,delta_v,hkgrid,vgrid,m1_v,m2_hk,T_hkv)
!        !f_in(lh1:lh2) = e**(-hkgrid(lh1:lh2)**2./2)
!        !CALL hankel_transform(f_in,.true.)
!        !CALL hankel_transform(f_in,.false.)
!     ELSE 
!        CALL GetMuWeightsAndKnots(delta_v,vgrid,vmax,nh0)
!     END IF
!  END IF

  !Set up FLR terms
!  CALL get_J0
!  
!  !Set up Hankel terms
!  CALL get_hk
!
!  IF(hyp_v.ne.0.0) THEN
!    CALL get_io_number
!    hypv_handle=io_number
!    IF(mype==0) THEN
!      OPEN(unit=hypv_handle,file=trim(diagdir)//'/hypv_vs_coll.dat',status='unknown')
!      WRITE(hypv_handle,*) "#coll=",nu
!      WRITE(hypv_handle,*) "#hypv=",hyp_v
!      WRITE(hypv_handle,*) "#hypv_order=",hypv_order
!      DO l=0,nv0-1
!        WRITE(hypv_handle,*) l,nu*herm_grid(l),hyp_v*(REAL(herm_grid(l))/REAL(nv0))**hypv_order 
!      END DO
!      CLOSE(hypv_handle)
!    END IF
!  END IF

  itime_start=itime

END SUBROUTINE arrays



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                   arrays_temp                              !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Sets up indices and allocates important arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE arrays_temp
  USE par_mod
  !USE flr_effects
  !USE hk_effects
  !USE Gaussquadrature
  IMPLICIT NONE
  
  INTEGER :: i,l,j,hypv_handle

  !Info for k grids
  !Note the format for the ky/kz indices:
  !Index:	0,    1,      2,. . . .       hky_ind,  nky0/2,  lky_ind, . . . nky0-1
  !ky val:	0.0,  kymin,  2*kymin . . . . kymax,    *,       -kymax, . . .  -kymin
  hkx_ind=nkx0-1     !Index of maximum kx value
  hky_ind=nky0/2-1     !Index of maximum (used,i.e. not dummy) ky value 
  lky_ind=nky0-hky_ind !Index of minimum (most negative) ky value
  hkz_ind=nkz0/2-1 !Index of maximum (used,i.e. not dummy) kz value
  lkz_ind=nkz0-hkz_ind !Index of minimum (most negative) kz value

  IF(verbose.and.(mype.eq.0)) THEN
    WRITE(*,*) "hkx_ind",hkx_ind
    WRITE(*,*) "hky_ind",hky_ind
    WRITE(*,*) "lky_ind",lky_ind
    WRITE(*,*) "hkz_ind",hkz_ind
    WRITE(*,*) "lkz_ind",lkz_ind
  END IF

!  !Info for v grid and parallelization
!  lv0=nv0  !Number of hermites per processor
!  lv1=0
!  lv2=lv0-1    !Upper index on processor mype
!  lbv=lv1-1             !Lower boundary index 
!  IF(mype_herm==0) lbv=0     
!  ubv=lv2+1             !Upper boundary index
!  IF(mype_herm==np_herm-1) ubv=lv2  

  !Info for kz grid
  lkz0=nkz0
  lkz1=0
  lkz2=lkz0-1
  !Verify the following when implementing kz parallelization
  lbkz=lkz1
  ubkz=lkz2

!  !Info for Hankel grid
!  lh0=nh0
!  lh1=0
!  lh2=lh0-1
!  !Verify the following when implementing h parallelization
!  !As a first implementation, no collision operator, so nhb = 0
!  nhb  = 0
!  lbh=lh1 - nhb
!  ubh=lh2 + nhb


  !Info for species grid
  ls0=nspec
  ls1=0
  ls2=ls0-1
  !Verify the following when implementing s parallelization
  lbs=ls1
  ubs=ls2

!   lxyzvhs0=nkx0*nky0*lkz0*lv0*lh0*ls0   
!   ev_size_loc=lxyzvhs0
!   ev_size=lxyzvhs0*np_herm*np_hank*np_kz*np_spec
  
!  IF(.not.allocated(g_1))&
!      ALLOCATE(g_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)) 
  IF(.not.allocated(b_1))&
      ALLOCATE(b_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(v_1))&
      ALLOCATE(v_1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(kxgrid)) ALLOCATE(kxgrid(0:nkx0-1))
  IF(.not.allocated(kygrid)) ALLOCATE(kygrid(0:nky0-1))
  IF(spatial2d) THEN
    !See Watanabe and Sugama '04
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nky0-1))
  ELSE
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nkz0-1))
  END IF
!  IF(.not.allocated(herm_grid)) ALLOCATE(herm_grid(0:nv0-1))
!  IF(.not.allocated(hgrid_loc)) ALLOCATE(hgrid_loc(lv1:lv2))
!  IF(.not.allocated(hkgrid)) ALLOCATE(hkgrid(0:nh0-1))
!  IF(.not.allocated(vgrid)) ALLOCATE(vgrid(0:nh0-1))
!  IF(.not.allocated(delta_hk)) ALLOCATE(delta_hk(0:nh0-1))
!  IF(.not.allocated(delta_v)) ALLOCATE(delta_v(0:nh0-1))
!  IF(.not.allocated(kperp2)) ALLOCATE(kperp2(0:nkx0-1,0:nky0-1))
!  IF(.not.allocated(hyp_x_herm1)) ALLOCATE(hyp_x_herm1(lv1:lv2)) 
!  IF(.not.allocated(hyp_y_herm1)) ALLOCATE(hyp_y_herm1(lv1:lv2)) 

!  hyp_x_herm1(:) = 0.001 
!  hyp_y_herm1(:) = 0.001 


  IF(kmin_eq_0) THEN
    DO i=0,nkx0-1
      kxgrid(i)=i*kxmin
    END DO 
    kxmax=(nkx0-1)*kxmin
    !WRITE(*,*) "kx grid:", kxgrid

    kygrid(0)=0.0
    DO i=1,nky0/2-1
      kygrid(i)=i*kymin
      kygrid(nky0-i)=-i*kymin
    END DO 
    kymax=(nky0/2-1)*kymin
    kygrid(nky0/2)=kymax+kymin  !dummy index
    !WRITE(*,*) "ky grid:", kygrid

    IF(nkz0.ge.2) THEN
      kzgrid(0)=0.0
      DO i=1,nkz0/2-1
        kzgrid(i)=i*kzmin
        kzgrid(nkz0-i)=-i*kzmin
      END DO 
      kzmax=(nkz0/2-1)*kzmin
      kzgrid(nkz0/2)=kzmax+kzmin  !dummy index
    ELSE
      !See Watanabe and Sugama '04
      kzgrid=kygrid*kzmin/kymin
      kzmax=kymax*kzmin/kymin
    END IF

  ELSE !.not.kmin_eq_0, i.e. linear
    kxgrid(0)=kxmin
    kxmax=kxmin
    kygrid(0)=kymin
    kymax=kymin
    kzgrid(0)=kzmin
    kzmax=kzmin
  END IF

!  DO i=0,nv0-1
!    herm_grid(i)=REAL(i)
!  END DO 
!
!  DO i=lv1,lv2
!    hgrid_loc(i)=REAL(i)
!  END DO 
!
!
!   DO i=0,nkx0-1
!    DO j=0,nky0-1
!      kperp2(i,j) = kxgrid(i)**2 + kygrid(j)**2
!    ENDDO
!   ENDDO
!
!  !Hankel grid
!  IF (.not.mu_integrated) THEN
!     IF (hankel) THEN
!        IF(.not.allocated(T_hkv)) ALLOCATE(T_hkv(0:nh0-1,0:nh0-1))
!        IF(.not.allocated(m1_v)) ALLOCATE(m1_v(0:nh0-1))
!        IF(.not.allocated(m2_hk)) ALLOCATE(m2_hk(0:nh0-1))
!        !IF(.not.allocated(f_in)) ALLOCATE(f_in(lh1:lh2))
!        hkmax = sqrt(2*maxval(kperp2))
!        call Gethankelgrid(hkmax,vmax,delta_hk,delta_v,hkgrid,vgrid,m1_v,m2_hk,T_hkv)
!        !f_in(lh1:lh2) = e**(-hkgrid(lh1:lh2)**2./2)
!        !CALL hankel_transform(f_in,.true.)
!        !CALL hankel_transform(f_in,.false.)
!     ELSE 
!        CALL GetMuWeightsAndKnots(delta_v,vgrid,vmax,nh0)
!     END IF
!  END IF

  !Set up FLR terms
!  CALL get_J0
!  
!  !Set up Hankel terms
!  CALL get_hk
!
!  IF(hyp_v.ne.0.0) THEN
!    CALL get_io_number
!    hypv_handle=io_number
!    IF(mype==0) THEN
!      OPEN(unit=hypv_handle,file=trim(diagdir)//'/hypv_vs_coll.dat',status='unknown')
!      WRITE(hypv_handle,*) "#coll=",nu
!      WRITE(hypv_handle,*) "#hypv=",hyp_v
!      WRITE(hypv_handle,*) "#hypv_order=",hypv_order
!      DO l=0,nv0-1
!        WRITE(hypv_handle,*) l,nu*herm_grid(l),hyp_v*(REAL(herm_grid(l))/REAL(nv0))**hypv_order 
!      END DO
!      CLOSE(hypv_handle)
!    END IF
!  END IF

  itime_start=itime

END SUBROUTINE arrays_temp




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                              finalize_arrays                              !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE finalize_arrays
  USE par_mod

!  IF(allocated(g_1)) DEALLOCATE(g_1)
  IF(allocated(b_1)) DEALLOCATE(b_1)
  IF(allocated(v_1)) DEALLOCATE(v_1)
!  IF(allocated(gpsi)) DEALLOCATE(gpsi)
!  IF(allocated(pre)) DEALLOCATE(pre)
  IF(allocated(pcurleig)) DEALLOCATE(pcurleig)
  
  if (verbose.and.(mype.eq.0)) print *, 'Deallocated b,v'

  IF(allocated(kxgrid)) DEALLOCATE(kxgrid)
  IF(allocated(kygrid)) DEALLOCATE(kygrid)
  IF(allocated(kzgrid)) DEALLOCATE(kzgrid)
  IF(allocated(kmags)) DEALLOCATE(kmags)
  IF(allocated(kzs)) DEALLOCATE(kzs)
  IF(allocated(alpha_leftwhist)) DEALLOCATE(alpha_leftwhist)
  IF(allocated(alpha_leftcyclo)) DEALLOCATE(alpha_leftcyclo)
  
!  IF(allocated(herm_grid)) DEALLOCATE(herm_grid)
!  IF(allocated(hgrid_loc)) DEALLOCATE(hgrid_loc)
!  IF(allocated(delta_hk)) DEALLOCATE(delta_hk)
!  IF(allocated(delta_v)) DEALLOCATE(delta_v)
!  IF(allocated(hkgrid)) DEALLOCATE(hkgrid)
!  IF(allocated(vgrid)) DEALLOCATE(vgrid)
!  IF(allocated(kperp2)) DEALLOCATE(kperp2)
!  IF(allocated(T_hkv)) DEALLOCATE(T_hkv)
!  IF(allocated(m1_v)) DEALLOCATE(m1_v)
!  IF(allocated(m2_hk)) DEALLOCATE(m2_hk)
!  IF(allocated(hyp_x_herm1)) DEALLOCATE(hyp_x_herm1)
!  IF(allocated(hyp_y_herm1)) DEALLOCATE(hyp_y_herm1)
  !IF(allocated(f_in)) DEALLOCATE(f_in)
  IF(verbose.and.(mype.eq.0)) print *, 'Deallocated kgrids'

END SUBROUTINE finalize_arrays

