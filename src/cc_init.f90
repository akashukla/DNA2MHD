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
  
  INTEGER :: i,j,k,l,hypv_handle
  INTEGER(C_INTPTR_T) :: zpad
  REAL :: b1r,b1i,b2r,b2i,b3r,b3i

  !Info for k grids
  !Note the format for the ky/kx indices:
  !Index:	0,    1,      2,. . . .       hky_ind,  nky0/2,  lky_ind, . . . nky0-1
  !ky val:	0.0,  kymin,  2*kymin . . . . kymax,    *,       -kymax, . . .  -kymin

  hkx_ind = nkx0-1
  hky_ind = nky0/2-1
  hkz_ind = nkz0/2-1
  lky_ind = nky0-hky_ind
  lkz_ind = nkz0-hkz_ind
  

  if (dealias_type.eq.1) zpad = 2
  if (((dealias_type.eq.3).or.(dealias_type.eq.4)).or.((dealias_type.eq.6).or.(dealias_type.eq.8))) zpad = dealias_type

  ny0_big=  zpad*nky0/2
  nx0_big = zpad*nkx0
  nz0_big = zpad*nkz0/2

  lky_big = ny0_big - hky_ind
  lkz_big = nz0_big - hkz_ind

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
  lkz0=nz0_big
  lkz1=mype*nz0_big/n_mpi_procs
  lkz2=(mype+1)*nz0_big/n_mpi_procs-1

  if (verbose) print *, mype,"Mype Limits",lkz1,lkz2

  ! Define upper and lower bounds for nonzero indices in z
  if (n_mpi_procs.eq.1) then
     lbkz = lkz1
     ubkz = lkz2
  else if (lkz1.le.hkz_ind) then
     lbkz = lkz1
     ubkz = min(hkz_ind,lkz2)
  else if (lkz2.ge.lkz_big) then
     lbkz = max(lkz1,lkz_big)
     ubkz = lkz2
  else
     ! No nonzero indices for lkz1 and lkz2 in the padding region                                                                                                                       
     lbkz = lkz2
     ubkz = lkz1
  endif

  ! Define mask for zero padding
  IF (.not.allocated(paddingmask)) ALLOCATE(paddingmask(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  paddingmask = 0
  if (n_mpi_procs.eq.1) then
     ! In the one process case hkz_ind and lkz_ind both lie between lkz1 and lkz2
     paddingmask(0:nkx0-1,0:hky_ind,0:hkz_ind) = 1
     paddingmask(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1) = 1
     paddingmask(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1) = 1
     paddingmask(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind) = 1
  else
     ! Use choice of lbkz and ubkz to account for zero padding
     paddingmask(0:nkx0-1,0:hky_ind,lbkz:ubkz) = 1
     paddingmask(0:nkx0-1,lky_big:ny0_big-1,lbkz:ubkz) = 1
  endif

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
  IF(.not.allocated(reader)) ALLOCATE(reader(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  IF(.not.allocated(reader2)) ALLOCATE(reader2(1:1+nx0_big/2,1:ny0_big,1+lkz1:1+lkz2))
  IF (.not.allocated(fullsmallarray)) ALLOCATE(fullsmallarray(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  IF (.not.allocated(fullbigarray)) ALLOCATE(fullbigarray(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))

  IF (.not.allocated(scatter_big)) ALLOCATE(scatter_big((1+nx0_big/2)*ny0_big*nz0_big/n_mpi_procs))
  IF (.not.allocated(scatter_small)) ALLOCATE(scatter_small((1+nx0_big/2)*ny0_big*nz0_big/n_mpi_procs))

  IF (.not.allocated(gather_big)) ALLOCATE(gather_big((1+nx0_big/2)*ny0_big*nz0_big))
  IF (.not.allocated(gather_small)) ALLOCATE(gather_small((1+nx0_big/2)*ny0_big*nz0_big))
  
  IF(.not.allocated(b_1))&
      ALLOCATE(b_1(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(v_1))&
      ALLOCATE(v_1(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)) 
  IF(.not.allocated(kxgrid)) ALLOCATE(kxgrid(0:nx0_big/2))
  IF(.not.allocated(kygrid)) ALLOCATE(kygrid(0:ny0_big-1))
  IF(spatial2d) THEN
    !kzgrid proportional to ky grid (see Watanabe and Sugama '04)
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nz0_big-1))
  ELSE
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(0:nz0_big-1))
 END IF
 
  IF(.not.allocated(kmags)) ALLOCATE(kmags(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  IF(.not.allocated(kperps)) ALLOCATE(kperps(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  IF(.not.allocated(kzs)) ALLOCATE(kzs(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
 ! IF(.not.allocated(gpsi)) ALLOCATE(gpsi(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ! IF(.not.allocated(pre)) ALLOCATE(pre(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  IF(.not.allocated(pcurleig)) ALLOCATE(pcurleig(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2))
  
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

    DO i=0,nx0_big/2
      kxgrid(i)=i*kxmin
    END DO
    kxmax=(nkx0-1)*kxmin

    kzgrid(0)=0.0
    DO i=1,nz0_big/2-1
      kzgrid(i)=i*kzmin
      kzgrid(nz0_big-i)=-i*kzmin
    END DO
    kzmax=(nkz0/2-1)*kzmin
    kzgrid(nz0_big/2)=nz0_big/2 * kzmin

    kygrid(0)=0.0
    DO i=1,ny0_big/2-1
      kygrid(i)=i*kymin
      kygrid(ny0_big-i)=-i*kymin
    END DO 
    kymax=(nky0/2-1)*kymin
    kygrid(ny0_big/2)=ny0_big/2 * kymin

  ELSE !.not.kmin_eq_0, i.e. linear
    kxgrid(0)=kxmin
    kxmax=kxmin
    kygrid(0)=kymin
    kymax=kymin
    kzgrid(0)=kzmin
    kzmax=kzmin
 END IF

 kmax = sqrt(kxmax**2.0 + kymax**2.0 + kzmax**2.0)

 DO i = 0,nx0_big/2
    DO j = 0,ny0_big-1
       DO k = lkz1,lkz2
          kmags(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2)
          kperps(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2)
          kzs(i,j,k) = kzgrid(k)
       END DO
    END DO
 END DO

 DO i = 0,nx0_big/2
    DO j = 0,ny0_big-1
       DO k = lkz1,lkz2
          IF ((kmags(i,j,k).ne.0).and.(kperps(i,j,k).ne.0)) THEN
             b1r = -1.0 / (kperps(i,j,k)*sqrt(2.0)) * kygrid(j)
             b2r = 1.0 / (kperps(i,j,k)*sqrt(2.0)) * kxgrid(i)
             b3r = 0.0
             b1i = (kygrid(j) * b3r - kzgrid(k) * b2r)/kmags(i,j,k)
             b2i = (kzgrid(k) * b1r - kxgrid(i) * b3r)/kmags(i,j,k)
             b3i = (kxgrid(i) * b2r - kygrid(j) * b1r)/kmags(i,j,k)
             pcurleig(i,j,k,0) = cmplx(b1r,b1i)
             pcurleig(i,j,k,1) = cmplx(b2r,b2i)
             pcurleig(i,j,k,2) = cmplx(b3r,b3i)
          ELSE
             pcurleig(i,j,k,:) = cmplx(0.0,0.0)
          ENDIF
          ! if (verbose) print *, i,j,k,sum(abs(pcurleig(i,j,k,:))**2)                                                                                                  
       ENDDO
    ENDDO
 ENDDO

 ALLOCATE(alpha_leftwhist(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
 ALLOCATE(alpha_leftcyclo(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

 alpha_leftwhist = -kmags/2.0 - sqrt(1.0 + ((kmags**2.0) / 4.0))
 alpha_leftcyclo = -kmags/2.0 + sqrt(1.0 + ((kmags**2.0) / 4.0))  
 
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
!!                              finalize_arrays                              !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE finalize_arrays
  USE par_mod

  !  IF(allocated(g_1)) DEALLOCATE(g_1)
  IF (allocated(reader)) DEALLOCATE(reader)
  IF (allocated(reader2)) DEALLOCATE(reader2)
  IF (allocated(fullbigarray)) DEALLOCATE(fullbigarray)
  IF (allocated(fullsmallarray)) DEALLOCATE(fullsmallarray)

  IF (allocated(gather_big)) DEALLOCATE(gather_big)
  IF (allocated(gather_small)) DEALLOCATE(gather_small)

  IF (allocated(scatter_big)) DEALLOCATE(scatter_big)
  IF (allocated(scatter_small)) DEALLOCATE(scatter_small)
  
  IF(allocated(b_1)) DEALLOCATE(b_1)
  IF(allocated(v_1)) DEALLOCATE(v_1)
  IF(allocated(paddingmask)) DEALLOCATE(paddingmask)
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

