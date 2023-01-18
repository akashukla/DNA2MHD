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
  INTEGER(kind=4) :: rkstage


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                iv_solver                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE iv_solver

 REAL :: dt_next
 IF(.not.checkpoint_read) dt=dt_max
 !itime=0
 !time=0.0

 IF(mype==0) WRITE(*,*) "max_itime=",max_itime
 IF(mype==0) WRITE(*,*) "max_time=",max_time
 rkstage = 1
 DO WHILE(time.lt.max_time.and.itime.lt.max_itime.and.continue_run)
 
   !IF(verbose) WRITE(*,*) "Calling diagnostics",time,itime,mype
   CALL diag
   !IF(verbose) WRITE(*,*) "Done with diagnostics",time,itime,mype

   IF(verbose.and.(mype.eq.0)) WRITE(*,*) "iv_solver: before get_g_next",time,itime,mype
   IF(verbose.and.(mype.eq.0)) WRITE(*,*) "iv_solver: before get_g_next dt=",dt
   !CALL save_b(b_1)
   !CALL save_time(itime)
   CALL get_g_next(b_1, v_1,dt_next)
   dt = minval([dt_next,dt_max])
   itime=itime+1 
   IF(mype==0.and.verbose) WRITE(*,*) "itime:",itime
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
 INTEGER :: ionums(9)
 INTEGER :: q

 IF ((mype.eq.0).and.(plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN
 ionums = [dbio,dvio,bdvio,vdbio,bdcbio,cbdbio,vdvio,bdbio,db2io]
 DO q = 1,9
    WRITE(ionums(q)) rkstage
 ENDDO 
 ENDIF

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
 CALL get_rhs(b_in, v_in, bk1, vk1,dt_new1)
 b_2=b_in+(1.0/6.0)*dt*bk1
 v_2=v_in+(1.0/6.0)*dt*vk1
 first_stage=.false.
 !CALL get_rhs(b_in+0.5*dt*bk1,bk2)
 bk1=b_in+0.5*dt*bk1
 vk1=v_in+0.5*dt*vk1
 
! IF (plot_nls.and.(mod(itime,istep_energy).eq.0)) THEN
! DO q = 1,9
!    WRITE(ionums(q)) rkstage
! ENDDO
! ENDIF
 
 CALL get_rhs(bk1,vk1,bk2,vk2,dt_new2)
 b_2=b_2+(1.0/3.0)*dt*bk2
 bk2=b_in+0.5*dt*bk2
 v_2=v_2+(1.0/3.0)*dt*vk2
 vk2=v_in+0.5*dt*vk2

 !IF (plot_nls.and.(mod(itime,istep_energy).eq.0)) THEN 
 !DO q = 1,9
 !   WRITE(ionums(q)) rkstage
 !ENDDO
 !ENDIF
 
 CALL get_rhs(bk2,vk2,bk1,vk1,dt_new3)
 b_2=b_2+(1.0/3.0)*dt*bk1
 bk1=b_in+dt*bk1
 v_2=v_2+(1.0/3.0)*dt*vk1
 vk1=v_in+dt*vk1

 !IF (plot_nls.and.(mod(itime,istep_energy).eq.0)) THEN
 !DO q = 1,9
 !   WRITE(ionums(q)) rkstage
 !ENDDO
 !ENDIF

 CALL get_rhs(bk1,vk1,bk2,vk2,dt_new4)
 b_in=b_2+(1.0/6.0)*dt*bk2
 v_in=v_2+(1.0/6.0)*dt*vk2

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
 CALL remove_div(b_in,v_in)

 DEALLOCATE(b_2)
 DEALLOCATE(bk1)
 DEALLOCATE(bk2)
 DEALLOCATE(v_2)
 DEALLOCATE(vk1)
 DEALLOCATE(vk2)

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

                b_in(i,j,k,0) = b_in(i,j,k,0) - div_b*kxgrid(i)/k2
                b_in(i,j,k,1) = b_in(i,j,k,1) - div_b*kygrid(j)/k2
                b_in(i,j,k,2) = b_in(i,j,k,2) - div_b*kzgrid(k)/k2

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
SUBROUTINE get_rhs(b_in,v_in, rhs_out_b,rhs_out_v,ndt)

 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
 REAL :: ndt

!  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!  COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  INTEGER :: k

  CALL get_rhs_lin(b_in,v_in,rhs_out_b, rhs_out_v,0)
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'Lin Abs Bx kmax',abs(rhs_out_b(nkx0-1,nky0/2,nkz0/2,0))
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'Lin Abs Vx kmax',abs(rhs_out_v(nkx0-1,nky0/2,nkz0/2,0))
  !IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v)
  IF(actual_nonlinear) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v,ndt)
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'NL Max Abs B',maxval(abs(rhs_out_v)),maxloc(abs(rhs_out_v))
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'NL Max Abs V',maxval(abs(rhs_out_b)),maxloc(abs(rhs_out_b))
  IF (.not.(actual_nonlinear)) ndt = dt_max

  ! Add forcing
  IF(force_turbulence) CALL get_rhs_force(rhs_out_b, rhs_out_v,ndt)
if (verbose.and.(mype.eq.0)) print *,'RHS found'
rkstage = rkstage + 1

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


END MODULE time_advance
