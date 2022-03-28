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
  USE field_solver
  USE calculate_time_step, ONLY: adapt_dt

  PUBLIC :: iv_solver,rk4_stability
  
  PRIVATE

!  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: g_2,k1,k2
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: b_2, bk1, bk2, v_2, vk1, vk2

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                iv_solver                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE iv_solver

 IF(.not.checkpoint_read) dt=dt_max
 !itime=0
 !time=0.0

 IF(mype==0) WRITE(*,*) "max_itime=",max_itime
 IF(mype==0) WRITE(*,*) "max_time=",max_time
 DO WHILE(time.lt.max_time.and.itime.lt.max_itime.and.continue_run)
 
   IF(verbose) WRITE(*,*) "Calling diagnostics",time,itime,mype
   CALL diag
   IF(verbose) WRITE(*,*) "Done with diagnostics",time,itime,mype

   IF(verbose) WRITE(*,*) "iv_solver: before get_g_next",time,itime,mype
   CALL get_g_next(b_1, v_1)
   itime=itime+1 
   IF(mype==0.and.verbose) WRITE(*,*) "itime:",itime
   time=time+dt
   IF(mype==0.and.dt_output) WRITE(*,*) "Before adapt: dt_max,dt",dt_max,dt
   IF(adapt_dt_nl) CALL adapt_dt 
   IF(mype==0.and.dt_output) WRITE(*,*) "After adapt: dt_max,dt",dt_max,dt
   CALL check_wallclock
   IF(current_wallclock.gt.REAL(max_walltime)) THEN
     IF(mype==0) WRITE(*,*) "Maximum wall time exceeded.", current_wallclock, max_walltime
     continue_run=.false. 
   END IF

 END DO
 IF(verbose) WRITE(*,*) "time,itime,mype",time,itime,mype

END SUBROUTINE iv_solver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_g_next                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_g_next(b_in, v_in)

! COMPLEX, INTENT(inout) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 
 
 ALLOCATE(b_2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ALLOCATE(bk1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ALLOCATE(bk2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))

 ALLOCATE(v_2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ALLOCATE(vk1(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ALLOCATE(vk2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))

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
 CALL get_rhs(b_in,bk1)
 b_2=b_in+(1.0/6.0)*dt*bk1
 first_stage=.false.
 !CALL get_rhs(b_in+0.5*dt*bk1,bk2)
 bk1=b_in+0.5*dt*bk1
 CALL get_rhs(bk1,bk2)
 b_2=b_2+(1.0/3.0)*dt*bk2
 bk2=b_in+0.5*dt*bk2
 CALL get_rhs(bk2,bk1)
 b_2=b_2+(1.0/3.0)*dt*bk1
 bk1=b_in+dt*bk1
 CALL get_rhs(bk1,bk2)
 b_in=b_2+(1.0/6.0)*dt*bk2
 !b_in=b_2 

 !4th order Runge-Kutta
 first_stage=.true.
 CALL get_rhs(v_in,vk1)
 v_2=v_in+(1.0/6.0)*dt*vk1
 first_stage=.false.
 !CALL get_rhs(v_in+0.5*dt*vk1,vk2)
 vk1=v_in+0.5*dt*vk1
 CALL get_rhs(vk1,vk2)
 v_2=v_2+(1.0/3.0)*dt*vk2
 vk2=v_in+0.5*dt*vk2
 CALL get_rhs(vk2,vk1)
 v_2=v_2+(1.0/3.0)*dt*vk1
 vk1=v_in+dt*vk1
 CALL get_rhs(vk1,vk2)
 v_in=v_2+(1.0/6.0)*dt*vk2
 !v_in=v_2 


 IF(force_kz0eq0) b_in(:,:,0,:,:,:)=cmplx(0.0,0.0)
 IF(force_ky0eq0) b_in(:,0,:,:,:,:)=cmplx(0.0,0.0)
 IF(force_kx0eq0) b_in(0,:,:,:,:,:)=cmplx(0.0,0.0)
 IF(nkz0.ge.2) b_in(:,:,hkz_ind+1,:,:,:)=cmplx(0.0,0.0)
 b_in(:,hky_ind+1,:,:,:,:)=cmplx(0.0,0.0)
 b_in(0,0,:,:,:,:)=cmplx(0.0,0.0)

 IF(force_kz0eq0) v_in(:,:,0,:,:,:)=cmplx(0.0,0.0)
 IF(force_ky0eq0) v_in(:,0,:,:,:,:)=cmplx(0.0,0.0)
 IF(force_kx0eq0) v_in(0,:,:,:,:,:)=cmplx(0.0,0.0)
 IF(nkz0.ge.2) v_in(:,:,hkz_ind+1,:,:,:)=cmplx(0.0,0.0)
 v_in(:,hky_ind+1,:,:,:,:)=cmplx(0.0,0.0)
 v_in(0,0,:,:,:,:)=cmplx(0.0,0.0)



 DEALLOCATE(b_2)
 DEALLOCATE(bk1)
 DEALLOCATE(bk2)
 DEALLOCATE(v_2)
 DEALLOCATE(vk1)
 DEALLOCATE(vk2)

END SUBROUTINE get_g_next


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                    get_rhs                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs(b_in, rhs_out_b, 0)

 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 !COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 !COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

!  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!  COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  INTEGER :: k

  CALL get_rhs_lin(b_in,v_in,rhs_out_b, rhs_out_v,0)

  IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v,0)

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

END SUBROUTINE rk4_stability


END MODULE time_advance
