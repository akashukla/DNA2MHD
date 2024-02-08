!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                           cc_time_advance.f90                             !!
!!                                                                           !!
!!  time_advance                                                             !!
!!  -- iv_solver                                                             !!
!!  -- get_g_next                                                            !!
!!  -- get_rhs                                                               !!
!!  -- get_next_rk_test                                                      !!
!!  -- get_rhs_rktest                                                        !!
!!  -- rk4_stability                                                         !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                             time_advance                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE time_advance
  USE par_mod
  USE linear_rhs
  USE nonlinearity
  USE diagnostics
  !USE field_solver
  !USE calculate_time_step, ONLY: adapt_dt

  PUBLIC :: iv_solver,rk4_stability
  
  PRIVATE
  
!  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: g_2,k1,k2
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: b_2, bk1, bk2, v_2, vk1, vk2
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: bk3,bk4,bk5,bk6,bk7,bk8,bk9,bk10,bk11,bk12,bk13,vk3,vk4,vk5,vk6,vk7,vk8,vk9,vk10,vk11,vk12,vk13,b_3,v_3

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                iv_solver                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE iv_solver

  REAL :: dt_next
  INTEGER :: ionums(9)
  INTEGER :: q
  
 IF(.not.checkpoint_read) dt=dt_max
 !itime=0
 !time=0.0

 IF(mype==0) WRITE(*,*) "max_itime=",max_itime
 IF(mype==0) WRITE(*,*) "max_time=",max_time

 CALL remove_div(b_1,v_1)
 CALL ALLOCATE_STEPS
 
 DO WHILE(time.lt.max_time.and.itime.lt.max_itime.and.continue_run)
 
   !IF(verbose) WRITE(*,*) "Calling diagnostics",time,itime,mype
   CALL diag
   !IF(verbose) WRITE(*,*) "Done with diagnostics",time,itime,mype

   IF(verbose.and.(mype.eq.0)) WRITE(*,*) "iv_solver: before get_g_next",time,itime,mype
   IF(verbose.and.(mype.eq.0)) WRITE(*,*) "iv_solver: before get_g_next dt=",dt
   !CALL save_b(b_1)
   !CALL save_time(itime)
   IF ((mype.eq.0).and.(plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN
      ionums = [dbio,dvio,bdvio,vdbio,bdcbio,cbdbio,vdvio,bdbio,db2io]
      DO q = 1,9
         WRITE(ionums(q)) time
      ENDDO
   ENDIF

   IF (dp547) THEN
      CALL dorpi547M(b_1,v_1)
   ELSE
      CALL get_g_next(b_1, v_1,dt_next)
      dt = minval([dt_next,dt_max])
   ENDIF
   
   itime=itime+1 
   IF(mype==0.and.verbose) WRITE(*,*) "itime: ",itime
   IF(mype==0.and.verbose) WRITE(*,*) "dt: ",dt
   time=time+dt
   !IF(mype==0.and.dt_output) WRITE(*,*) "Before adapt: dt_max,dt",dt_max,dt
   !IF(adapt_dt_nl) CALL adapt_dt 
   !IF(mype==0.and.dt_output) WRITE(*,*) "After adapt: dt_max,dt",dt_max,dt
   CALL check_wallclock
   IF(current_wallclock.gt.REAL(max_walltime)) THEN
   !  IF(mype==0) WRITE(*,*) "Maximum wall time exceeded.", current_wallclock, max_walltime
     continue_run=.false. 
   ENDIF
   IF (dt < 1.0E-5) then 
     WRITE(*,*) "dt too small to proceed" 
     continue_run=.false.
   ENDIF
  !END IF
END DO

CALL DEALLOCATE_STEPS

 write(*,*) 'Simulation Time: ',current_wallclock
 IF(verbose.and.(mype.eq.0)) WRITE(*,*) "time,itime,mype",time,itime,mype

END SUBROUTINE iv_solver

!!!!!!!!!!!!!!!!!!!!!!!!
!! save_b             !!
!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE save_b(b_out)
!COMPLEX, INTENT(in) :: b_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
!CHARACTER(len=100) :: chp_name
!INTEGER :: chp_handle
!LOGICAL :: not_first,b_output,v_output
!
!  b_output=.true.
!  chp_name='/b_out.dat'
!  INQUIRE(file=trim(diagdir)//trim(chp_name),exist=not_first)
!
!  IF(not_first) THEN
!    IF(mype==0) OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name), &
!         form='unformatted', status='unknown',access='stream',position='append')
!  ELSE
!    IF(mype==0) OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name),&
!                           form='unformatted', status='REPLACE',access='stream')
!  END IF
!
!  WRITE(chp_handle) b_out(1,1,1,0)
!
!END SUBROUTINE save_b
!
!SUBROUTINE save_time(itime)
!INTEGER, INTENT(in) :: itime
!CHARACTER(len=100) :: chp_name
!INTEGER :: chp_handle
!LOGICAL :: not_first,t_output
!
!  t_output=.true.
!  chp_name='/time_out.dat'
!  INQUIRE(file=trim(diagdir)//trim(chp_name),exist=not_first)
!
!  IF(not_first) THEN
!    IF(mype==0) OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name), &
!         form='unformatted', status='unknown',access='stream',position='append')
!  ELSE
!    IF(mype==0) OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name),&
!                           form='unformatted', status='REPLACE',access='stream')
!  END IF
!
!  WRITE(chp_handle) itime 
!
!END SUBROUTINE save_time



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_g_next                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_g_next(b_in, v_in,dt_new)

! COMPLEX, INTENT(inout) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
 COMPLEX, INTENT(inout) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(inout) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 REAL :: dt_new
 REAL :: dt_new1,dt_new2,dt_new3,dt_new4
 REAL :: nmhc

!  !4th order Runge-Kutta
!  first_stage=.true.
!  CALL get_rhs(g_in,k1)
!  g_2=g_in+(1.0/6.0)*dt*k1
!  first_stage=.false.
!  !CALL get_rhs(g_in+0.5*dt*k1,k2)
!  k1=g_in+0.5*dt*k1

!  CALL get_rhs(k1,k2)
!  g_2=g_2+(1.0/3.0)*dt*k2
!  k2=g_in+0.5*dt*k2

!  CALL get_rhs(k2,k1)
!  g_2=g_2+(1.0/3.0)*dt*k1
!  k1=g_in+dt*k1

!  CALL get_rhs(k1,k2)
!  g_in=g_2+(1.0/6.0)*dt*k2
!  !g_in=g_2 

 !4th order Runge-Kutta

 first_stage=.true.
 CALL get_rhs(b_in, v_in, bk1, vk1,nmhc,dt_new1)
 b_2=b_in+(1.0/6.0)*dt*bk1
 v_2=v_in+(1.0/6.0)*dt*vk1
 IF (mhc) mhelcorr = mhelcorr + (1.0/6.0)*dt*nmhc
 first_stage=.false.
 !CALL get_rhs(b_in+0.5*dt*bk1,bk2)
 bk1=b_in+0.5*dt*bk1
 vk1=v_in+0.5*dt*vk1
 
 CALL get_rhs(bk1,vk1,bk2,vk2,nmhc,dt_new2)
 b_2=b_2+(1.0/3.0)*dt*bk2
 bk2=b_in+0.5*dt*bk2
 v_2=v_2+(1.0/3.0)*dt*vk2
 vk2=v_in+0.5*dt*vk2
 IF (mhc) mhelcorr = mhelcorr + (1.0/3.0) *dt * nmhc

 CALL get_rhs(bk2,vk2,bk1,vk1,nmhc,dt_new3)
 b_2=b_2+(1.0/3.0)*dt*bk1
 bk1=b_in+dt*bk1
 v_2=v_2+(1.0/3.0)*dt*vk1
 vk1=v_in+dt*vk1
 IF (mhc) mhelcorr = mhelcorr + (1.0/3.0)*dt*nmhc

 CALL get_rhs(bk1,vk1,bk2,vk2,nmhc,dt_new4)
 b_in=b_2+(1.0/6.0)*dt*bk2
 v_in=v_2+(1.0/6.0)*dt*vk2
 IF (mhc) mhelcorr = mhelcorr + (1.0/6.0)*dt*nmhc

 !!4th order Runge-Kutta
 !first_stage=.true.
 !CALL get_rhs(v_in,vk1)
 !v_2=v_in+(1.0/6.0)*dt*vk1
 !first_stage=.false.
 !!CALL get_rhs(v_in+0.5*dt*vk1,vk2)
 !vk1=v_in+0.5*dt*vk1
 !CALL get_rhs(vk1,vk2)
 !v_2=v_2+(1.0/3.0)*dt*vk2
 !vk2=v_in+0.5*dt*vk2
 !CALL get_rhs(vk2,vk1)
 !v_2=v_2+(1.0/3.0)*dt*vk1
 !vk1=v_in+dt*vk1
 !CALL get_rhs(vk1,vk2)
 !v_in=v_2+(1.0/6.0)*dt*vk2
 !!v_in=v_2 

 !Remove below blocks for now for linear run test
!  IF(force_kz0eq0) b_in(:,:,0,:)=cmplx(0.0,0.0)
!  IF(force_ky0eq0) b_in(:,0,:,:)=cmplx(0.0,0.0)
!  IF(force_kx0eq0) b_in(0,:,:,:)=cmplx(0.0,0.0)
!  IF(nkz0.ge.2) b_in(:,:,hkz_ind+1,:)=cmplx(0.0,0.0)
!  b_in(:,hky_ind+1,:,:)=cmplx(0.0,0.0)
!  b_in(0,0,:,:)=cmplx(0.0,0.0)
! 
!  IF(force_kz0eq0) v_in(:,:,0,:)=cmplx(0.0,0.0)
!  IF(force_ky0eq0) v_in(:,0,:,:)=cmplx(0.0,0.0)
!  IF(force_kx0eq0) v_in(0,:,:,:)=cmplx(0.0,0.0)
!  IF(nkz0.ge.2) v_in(:,:,hkz_ind+1,:)=cmplx(0.0,0.0)
!  v_in(:,hky_ind+1,:,:)=cmplx(0.0,0.0)
!  v_in(0,0,:,:)=cmplx(0.0,0.0)

 dt_new = minval([dt_new1,dt_new2,dt_new3,dt_new4])

 if ((walenp).or.(walenn)) then
  if (walenp) print *, maxval(abs(b_in-v_in)),maxloc(abs(b_in-v_in))
  if (walenn) print *, maxval(abs(b_in+v_in)),maxloc(abs(b_in+v_in))
 endif

END SUBROUTINE get_g_next

SUBROUTINE remove_div(b_in,v_in)

 COMPLEX, INTENT(inout) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(inout) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 INTEGER :: i,j,k,l,h,ierr
 COMPLEX :: div_v, div_b
 REAL :: k2, zero

 div_v = 0.0 +i_complex*0.0
 div_b = 0.0 +i_complex*0.0
 k2=0.0
 zero=0.0

 DO i=0,nkx0-1
   DO j=0,nky0-1
     DO k=lkz1,lkz2
        k2 = kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2
        div_v = kxgrid(i)*v_in(i,j,k,0) + kygrid(j)*v_in(i,j,k,1) + kzgrid(k)*v_in(i,j,k,2)
        div_b = kxgrid(i)*b_in(i,j,k,0) + kygrid(j)*b_in(i,j,k,1) + kzgrid(k)*b_in(i,j,k,2)
        IF(k2.ne.zero) THEN 
                v_in(i,j,k,0) = v_in(i,j,k,0) - div_v*kxgrid(i)/k2
                v_in(i,j,k,1) = v_in(i,j,k,1) - div_v*kygrid(j)/k2
                v_in(i,j,k,2) = v_in(i,j,k,2) - div_v*kzgrid(k)/k2

                ! Constructs the pressure since the above term is - grad p
                pre(i,j,k) = i_complex * div_v / k2

              ! The b equation is a curl, so we don't need to remove div b
                b_in(i,j,k,0) = b_in(i,j,k,0) - div_b*kxgrid(i)/k2
                b_in(i,j,k,1) = b_in(i,j,k,1) - div_b*kygrid(j)/k2
                b_in(i,j,k,2) = b_in(i,j,k,2) - div_b*kzgrid(k)/k2

              ! Could be used to effect a gauge change such that dA/dt = v x B - curl B x B
              ! Keeps track of the integral of the pressure

              gpsi(i,j,k,0) = gpsi(i,j,k,0) - div_v*kxgrid(i)/k2
              gpsi(i,j,k,1) = gpsi(i,j,k,1) - div_v*kygrid(j)/k2
              gpsi(i,j,k,2) = gpsi(i,j,k,2) - div_v*kzgrid(k)/k2
         ENDIF 
     ENDDO
   ENDDO
 ENDDO 

if (verbose.and.(mype.eq.0)) print *,'Divergence Removed'

END SUBROUTINE remove_div

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                    get_rhs                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs(b_in,v_in, rhs_out_b,rhs_out_v,nmhc,ndt)

 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 REAL, INTENT(out) :: nmhc
 REAL :: ndt

!  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!  COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  INTEGER :: k

  IF (test_ho) THEN
     CALL get_rhs_test(b_in,v_in,rhs_out_b,rhs_out_v)
     ndt = dt_max
  ELSE
     
     IF (guide) CALL get_rhs_lin(b_in,v_in,rhs_out_b, rhs_out_v,0)
     IF ((guide.and.verbose).and.(mype.eq.0)) WRITE(*,*) 'Lin Max Abs V',maxval(abs(rhs_out_v)),maxloc(abs(rhs_out_v))
     IF ((guide.and.verbose).and.(mype.eq.0)) WRITE(*,*) 'Lin Max Abs B',maxval(abs(rhs_out_b)),maxloc(abs(rhs_out_b))
     !IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v)
     IF(actual_nonlinear) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v,ndt)
     IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'NL Max Abs V',maxval(abs(rhs_out_v)),maxloc(abs(rhs_out_v))
     IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'NL Max Abs B',maxval(abs(rhs_out_b)),maxloc(abs(rhs_out_b))
     IF (.not.(actual_nonlinear)) ndt = dt_max
     if (verbose.and.(mype.eq.0)) print *, ndt

     CALL remove_div(rhs_out_b,rhs_out_v)

     nmhc = 0.0
     IF (mhc) CALL getmhcrk(b_in,v_in,nmhc)

  ! Add forcing
     IF(force_turbulence) CALL get_rhs_force(rhs_out_b, rhs_out_v,ndt)

     if (verbose.and.(mype.eq.0)) print *,'RHS found'
  ENDIF
  

END SUBROUTINE get_rhs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               get_next_rk_test                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_next_rk_test(g_in,g_out,lambda,dt_test)

 COMPLEX, INTENT(in) :: g_in
 COMPLEX, INTENT(in) :: lambda
 COMPLEX, INTENT(out) :: g_out
 REAL, INTENT(in) :: dt_test
 COMPLEX :: k1
 COMPLEX :: k2

 !4th order Runge-Kutta
 CALL get_rhs_rktest(g_in,k1,lambda)
 g_out=g_in+(1.0/6.0)*dt_test*k1
 CALL get_rhs_rktest(g_in+0.5*dt_test*k1,k2,lambda)
 g_out=g_out+(1.0/3.0)*dt_test*k2
 CALL get_rhs_rktest(g_in+0.5*dt_test*k2,k1,lambda)
 g_out=g_out+(1.0/3.0)*dt_test*k1
 CALL get_rhs_rktest(g_in+dt_test*k1,k2,lambda)
 g_out=g_out+(1.0/6.0)*dt_test*k2

END SUBROUTINE get_next_rk_test


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               get_rhs_rktest                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_rktest(g_,k_,lambda)

  IMPLICIT NONE
  COMPLEX, INTENT(in) :: g_
  COMPLEX, INTENT(in) :: lambda
  COMPLEX, INTENT(out) :: k_

  k_=lambda*g_

END SUBROUTINE get_rhs_rktest


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              rk4_stability                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rk4_stability

  IMPLICIT NONE

  COMPLEX :: g1,g2
  COMPLEX :: lambda
  INTEGER :: i,j,numr,numi
  REAL :: rlam,ilam
  !COMPLEX :: k0,kn
  COMPLEX :: omega
  REAL, ALLOCATABLE, DIMENSION(:) :: rout,iout
  REAL, ALLOCATABLE, DIMENSION(:) :: rlamout,ilamout
  CHARACTER(len=5) :: chnumr,chnumi
  !INTEGER :: ierr

  IF(np_herm .gt.1) STOP "rk test only for serial execution."
  
  c0 = cmplx(rc0,ic0)
IF (mype.eq.0) THEN
  OPEN(unit=100,file=trim(diagdir)//'/rout.dat',status='unknown')
  OPEN(unit=200,file=trim(diagdir)//'/iout.dat',status='unknown')
  OPEN(unit=300,file=trim(diagdir)//'/rgrid.dat',status='unknown')
  OPEN(unit=400,file=trim(diagdir)//'/igrid.dat',status='unknown')

  WRITE(*,*) "rmax,imax",rmax,imax
  numr=nint(rmax*2.0/delta_lambda)
  numi=nint(imax*2.0/delta_lambda)
  WRITE(chnumr,'(i5.5)') numr
  WRITE(chnumi,'(i5.5)') numi
  WRITE(*,*) "rmax_index",numr
  WRITE(*,*) "imax_index",numi
  WRITE(*,*) chnumr
  WRITE(*,*) chnumi
  ALLOCATE(rout(numr))
  ALLOCATE(iout(numr))
  ALLOCATE(rlamout(numr))
  ALLOCATE(ilamout(numi))
  DO j=1,numi
    DO i=1,numr
      !IF(j==1) WRITE(100,*) ""
      !IF(j==1) WRITE(200,*) ""
      rlam=-1.0*rmax+delta_lambda*(i-1) 
      ilam=-1.0*imax+delta_lambda*(j-1) 
      IF(j==1) THEN
          rlamout(i)=rlam
          WRITE(300,*) rlam
      END IF
      IF(i==1) THEN
           ilamout(j)=ilam
          WRITE(400,*) ilam
      END IF
      lambda=cmplx(rlam,ilam)
      g1=cmplx(1.0,0.0)

      CALL get_next_rk_test(g1,g2,lambda,dt_rktest)

      omega=(g2-g1)/(0.5*dt_rktest*(g2+g1)) 
      !omega=(g2-g1)/(dt*g1) 
      rout(i)=REAL(omega)
      iout(i)=aimag(omega)

!!!!!!!!!!!!!!Test
      rout(i)=(REAL(sqrt(conjg(g2)*g2))-REAL(sqrt(conjg(g1)*g1)))/dt_rktest
!!!!!!!!!!!!!!Test
      rlamout(i)=rlam 
      ilamout(i)=ilam 
      

    END DO
    WRITE(100,'('//chnumr//'es16.8)') rout
    WRITE(200,'('//chnumr//'es16.8)') iout
  END DO

  CLOSE(100)
  CLOSE(200)
  CLOSE(300)
  CLOSE(400)



  OPEN(unit=100,file=trim(diagdir)//'/SR_test.dat',status='unknown')
  DO j=1,numi
      !IF(j==1) WRITE(100,*) ""
      !IF(j==1) WRITE(200,*) ""
      rlam=0.0
      ilam=-1.0*imax+delta_lambda*(j-1) 
      lambda=cmplx(rlam,ilam)
      g1=cmplx(1.0,0.0)

      CALL get_next_rk_test(g1,g2,lambda,dt_rktest)
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!      CALL get_rhs(g1,lambda,k0)  !get k1
!      g2=g1+1.0/6.0*k0*dt     
!      CALL get_rhs(g1+0.5*dt*k0,lambda,kn) !get k2
!      g2=g2+1.0/3.0*kn*dt
!      CALL get_rhs(g1+0.5*dt*kn,lambda,k0) !get k3
!      g2=g2+1.0/3.0*k0*dt
!      CALL get_rhs(g1+dt*k0,lambda,kn)     !get k4
!      g2=g2+1.0/6.0*kn*dt
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!

      omega=(g2-g1)/(0.5*dt_rktest*(g2+g1)) 
      !omega=(g2-g1)/(dt*g1) 
!!!!!!!!!!!!!!Test
!      rout(i)=(REAL(sqrt(conjg(g2)*g2))-REAL(sqrt(conjg(g1)*g1)))/dt
!!!!!!!!!!!!!!Test
      WRITE(100,*) ilam,(REAL(sqrt(conjg(g2)*g2))-REAL(sqrt(conjg(g1)*g1)))/dt_rktest
  END DO
  CLOSE(100)

  OPEN(unit=100,file=trim(diagdir)//'/SI_test.dat',status='unknown')
  DO i=1,numr
      !IF(j==1) WRITE(100,*) ""
      !IF(j==1) WRITE(200,*) ""
      ilam=0.0
      rlam=-1.0*rmax+delta_lambda*(i-1) 
      lambda=cmplx(rlam,ilam)
      g1=cmplx(1.0,0.0)

      CALL get_next_rk_test(g1,g2,lambda,dt_rktest)
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!      CALL get_rhs(g1,lambda,k0)  !get k1
!      g2=g1+1.0/6.0*k0*dt     
!      CALL get_rhs(g1+0.5*dt*k0,lambda,kn) !get k2
!      g2=g2+1.0/3.0*kn*dt
!      CALL get_rhs(g1+0.5*dt*kn,lambda,k0) !get k3
!      g2=g2+1.0/3.0*k0*dt
!      CALL get_rhs(g1+dt*k0,lambda,kn)     !get k4
!      g2=g2+1.0/6.0*kn*dt
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!RK4!!!!!!!!!!!!!!!!

      !omega=(g2-g1)/(dt*g1) 
      omega=(g2-g1)/(0.5*dt_rktest*(g2+g1)) 

      WRITE(100,*) rlam,REAL(omega)
  END DO
  CLOSE(100)

ENDIF

END SUBROUTINE rk4_stability

SUBROUTINE getmhcrk(b_in,v_in,nmhc)

 IMPLICIT NONE
 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

 REAL, intent(out) :: nmhc

nmhc = - 32.0*(pi**3) * sum(real(b_in(:,:,:,0)*conjg(v_in(:,:,:,1))-b_in(:,:,:,1)*conjg(v_in(:,:,:,0))))

END SUBROUTINE getmhcrk

SUBROUTINE dorpi547M(b_in,v_in)

  !! Integrating using adaptive DORPI 547M algorithm outlined in
  !! Dormand Prince A family of embedded Runge-Kutta formulae 1980
  !! https://doi.org/10.1016/0771-050X(80)90013-3
  !! Butcher tableau for reference
  !! 0    | 
  !! 1/5  | 1/5
  !! 3/10 | 3/40       9/40
  !! 4/5  | 44/45      -56/15      32/9
  !! 8/9  | 19372/6561 -25360/2187 64448/6561  -212/729
  !! 1    | 9017/3168  -355/33     46732/5247  49/176     -5103/18656
  !! 1    | 35/384     0           500/1113    125/192    -2187/6784      11/84
  !! _____________
  !!      | 35/384     0           500/1113    125/192    -2187/6784      11/84
  !!      | 5179/57600 0           7571/16695  393/640     -92097/339200  187/2100  1/40
  
  IMPLICIT NONE

  COMPLEX, INTENT(INOUT) :: b_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2)
  COMPLEX, INTENT(INOUT) :: v_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2)
  REAL :: dt_new1,dt_new2,dt_new3,dt_new4,dt_new5,dt_new6,dt_new7
  REAL :: nmhc1,nmhc2,nmhc3,nmhc4,nmhc5,nmhc6,nmhc7,next_corr
  LOGICAL :: trap = .true.
  
  dt = dt_max

  ! print *, "dt start", dt
  ! print *, "dt_max", dt_max
  !! First step

  trap = .true.
  
  do while (trap)


     ! Use First Same As Last property to save time for higher steps
     if (itime.eq.0) CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,dt_new1)
     IF (mhc.and.(itime.ne.0)) CALL getmhcrk(b_in,v_in,nmhc1)
     ! CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,dt_new1) 
     
     !! Second step
     CALL get_rhs(b_in+dt*bk1*1.0/5.0,&
          v_in+dt*vk1*1.0/5.0,bk2,vk2,nmhc2,dt_new2)
  
     !! Third step
     CALL get_rhs(b_in+dt*(3.0*bk1/40.0+9.0*bk2/40.0),&
          v_in+dt*(3.0*vk1/40.0+9.0*vk2/40.0),bk3,vk3,nmhc3,dt_new3)

     !! Fourth step
     CALL get_rhs(b_in+dt*(44.0/45.0*bk1-56.0/15.0*bk2+32.0/9.0*bk3),&
          v_in+dt*(44.0/45.0*vk1-56.0/15.0*vk2+32.0/9.0*vk3),bk4,vk4,nmhc4,dt_new4)

     !! Fifth step
     CALL get_rhs(b_in+dt*(19372.0/6561.0*bk1-25360.0/2187.0*bk2+64448.0/6561.0*bk3-212.0/729.0*bk4),&
          v_in+dt*(19372.0/6561.0*vk1-25360.0/2187.0*vk2+64448.0/6561.0*vk3-212.0/729.0*vk4),&
          bk5,vk5,nmhc5,dt_new5)

     !! Sixth step
     CALL get_rhs(b_in+dt*(9017.0/3168.0*bk1 - 355.0/33.0*bk2 + 46732.0/5247.0*bk3 &
          + 49.0/176.0*bk4-5103.0/18656.0*bk5),&
          v_in+dt*(9017.0/3168.0*vk1 - 355.0/33.0*vk2 + 46732.0/5247.0*vk3	&
       + 49.0/176.0*vk4-5103.0/18656.0*vk5),bk6,vk6,nmhc6,dt_new6)

     b_2 = b_in + dt*(35.0/384.0*bk1+500.0/1113.0*bk3+125.0/192.0*bk4&
          -2187.0/6784.0*bk5+11.0/84.0*bk6)
     v_2 = v_in + dt*(35.0/384.0*vk1+500.0/1113.0*vk3+125.0/192.0*vk4&
          -2187.0/6784.0*vk5+11.0/84.0*vk6)
     next_corr = dt*(35.0/384.0*nmhc1+500.0/1113.0*nmhc3+125.0/192.0*nmhc4&
          -2187.0/6784.0*nmhc5+11.0/84.0*nmhc6)
     
     !! Seventh step
     CALL get_rhs(b_2,v_2,bk7,vk7,nmhc7,dt_new7)

     b_3 = b_in + dt*(5179.0/57600.0*bk1+7571.0/16695.0*bk3+393.0/640.0*bk4&
          -92097.0/339200.0*bk5+187.0/2100.0*bk6+bk7/40.0)
     v_3 = v_in + dt*(5179.0/57600.0*vk1+7571.0/16695.0*vk3+393.0/640.0*vk4&
          -92097.0/339200.0*vk5+187.0/2100.0*vk6+vk7/40.0)

     trap = ((maxval(abs((b_3-b_2)/b_2),abs(b_2).gt.10**(-10.0)).gt.10**(-3.0)).or.&
          (maxval(abs((v_3-v_2)/v_2),abs(v_2).gt.10**(-10.0)).gt.10**(-3.0)))
     if ((itime.lt.1)) dt = minval([dt_new1,dt_new2,dt_new3,dt_new4,dt_new5,dt_new6,dt_new7,dt])
     if ((itime.gt.1)) dt = minval([dt_new2,dt_new3,dt_new4,dt_new5,dt_new6,dt_new7,dt])
     if (trap) dt = 0.9 * dt
     
     if (verbose) then
        print *, "Trap",trap
        print *, "RATIOS",(maxval(abs((b_3-b_2)/b_2),abs(b_2).gt.10**(-10.0))),maxval(abs((v_3-v_2)/v_2),abs(v_2).gt.10**(-10.0))
        print *, "dt",dt
        print *, "itime",itime
     endif
     
  enddo

  b_in = b_2
  v_in = v_2
  mhelcorr = mhelcorr + next_corr

  bk1 = bk7
  vk1 = vk7

END SUBROUTINE dorpi547M

SUBROUTINE ALLOCATE_STEPS

  IMPLICIT NONE

  ALLOCATE(b_2(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(v_2(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(b_3(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(v_3(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  
  ALLOCATE(bk1(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk1(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk2(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk2(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk3(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk3(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk4(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk4(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk5(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk5(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk6(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk6(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk7(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk7(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk8(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk8(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk9(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk9(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk10(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk10(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk11(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk11(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk12(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk12(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

  ALLOCATE(bk13(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))
  ALLOCATE(vk13(0:nkx0-1,0:nky0-1,0:nkz0-1,0:2))

END SUBROUTINE ALLOCATE_STEPS

SUBROUTINE DEALLOCATE_STEPS

  IMPLICIT NONE

  DEALLOCATE(b_2)
  DEALLOCATE(v_2)

  DEALLOCATE(b_3)
  DEALLOCATE(v_3)

  DEALLOCATE(bk1)
  DEALLOCATE(vk1)

  DEALLOCATE(bk2)
  DEALLOCATE(vk2)

  DEALLOCATE(bk3)
  DEALLOCATE(vk3)

  DEALLOCATE(bk4)
  DEALLOCATE(vk4)

  DEALLOCATE(bk5)
  DEALLOCATE(vk5)

  DEALLOCATE(bk6)
  DEALLOCATE(vk6)

  DEALLOCATE(bk7)
  DEALLOCATE(vk7)

  DEALLOCATE(bk8)
  DEALLOCATE(vk8)

  DEALLOCATE(bk9)
  DEALLOCATE(vk9)

  DEALLOCATE(bk10)
  DEALLOCATE(vk10)

  DEALLOCATE(bk11)
  DEALLOCATE(vk11)

  DEALLOCATE(bk12)
  DEALLOCATE(vk12)

  DEALLOCATE(bk13)
  DEALLOCATE(vk13)
  
END SUBROUTINE DEALLOCATE_STEPS


END MODULE time_advance
