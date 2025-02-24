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
  USE mpi
  USE linear_rhs
  USE nonlinearity
  USE diagnostics
  
  !USE field_solver
  !USE calculate_time_step, ONLY: adapt_dt

  PUBLIC :: iv_solver,rk4_stability
  
  PRIVATE
  
!  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: g_2,k1,k2
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: b_2, bk1, bk2, v_2, vk1, vk2
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: bk1s,vk1s,bk2s,vk2s
  
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) :: bk3,bk4,bk5,bk6,bk7,bk8,bk9,bk10,bk11,bk12,bk13,vk3,vk4,vk5,vk6,vk7,vk8,vk9,vk10,vk11,vk12,vk13,b_3,v_3
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: bsph0,vsph0,bsph1,vsph1,bsph3,vsph3,bsphk1,bsphk2,bsphk3,bsphk4,bsphk5,bsphk6,bsphk7,&
       vsphk1,vsphk2,vsphk3,vsphk4,vsphk5,vsphk6,vsphk7 ! z axis spherical chart
  LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: breakz ! if a chart fails; considered unlikely but possible
  INTEGER :: ierr

  !COMPLEX :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  !COMPLEX :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  !COMPLEX :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  !COMPLEX :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                iv_solver                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE iv_solver

  REAL(C_DOUBLE) :: dt_next
  INTEGER :: ionums(9)
  INTEGER :: q,ind
  REAL :: sttime,diagtime

 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
 IF(.not.checkpoint_read) dt=dt_max
 !itime=0
 !time=0.0

 WRITE(*,*) "max_itime=",max_itime,mype
 WRITE(*,*) "max_time=",max_time,mype

 CALL remove_div(b_1,v_1)

 if (verbose) then
    DO ind = 0,2
       print *, "Max bind",maxval(abs(b_1(:,:,:,ind)))
       print *, "Max vind",maxval(abs(v_1(:,:,:,ind)))
    ENDDO
 endif
 
 CALL ALLOCATE_STEPS

 rkstage = 0
 if (verbose) print *, mype,"Starting Main Loop"
 ! print *, "Time Step Start",dt
 
 DO WHILE((time.lt.max_time).and.(itime.lt.max_itime).and.(continue_run))
 
   IF(verbose) WRITE(*,*) "Calling diagnostics",time,itime,mype

   ! Don't call diagnostics on first iteration when doing a warm restart - no doubles in data
   if (timer.and.(mype.eq.0)) sttime = MPI_WTIME()
   if ((.not.checkpoint_read).or.(itime.gt.itime_start)) CALL diag
   if (timer.and.(mype.eq.0)) diagtime = MPI_WTIME()
   if (timer.and.(mype.eq.0)) print *, "Diag Time",diagtime - sttime,mype
   
   !IF(verbose) WRITE(*,*) "Done with diagnostics",time,itime,mype

   IF(verbose) WRITE(*,*) "iv_solver: before get_g_next",time,itime,mype
   IF(verbose) WRITE(*,*) "iv_solver: before get_g_next dt=",dt
   !CALL save_b(b_1)
   !CALL save_time(itime)

   !! Rules for Time Integrators
   !! intorder < 10 are explicit, nonsplitting methods
   !! intorder 20 is a Strang splitting 
   !! intorder 40 is a Suzuki fourth order splitting
   !! intorder 41 is Gauss 4th order collocation; may have trouble converging for nonideal
   !! intorder 81 is Stormer-Verlet 8th order - do not use for now
   
   
   ! Implicit + Strucuture Conserving Integrators with Splitting for Nonideal fluid
   IF (intorder.eq.20) THEN
      CALL split2(b_1,v_1)
   ELSE IF (intorder.eq.40) THEN
      CALL split4(b_1,v_1)
   ! Implicit + Structure Conserving Integrators (Use on Ideal Fluid only until better understood)
   ELSE IF (intorder.eq.41) THEN
      CALL gauss4(b_1,v_1,dt_next)
      dt = minval([dt_max,dt_next])
   ELSE IF (intorder.eq.21) THEN
      CALL gauss2(b_1,v_1,dt_next)
      dt = minval([dt_max,dt_next])
      ! Dormand Prince Adaptive/Embedded Schemes
   ELSE IF (intorder.eq.5) THEN
      CALL dorpi547M(b_1,v_1)
   ! Other Explicit RK Integrators
   ELSE IF (intorder.eq.4) THEN
      CALL get_g_next(b_1, v_1,dt_next)
      dt = minval([dt_next,dt_max])
   ELSE IF (intorder.eq.3) THEN
      CALL ralston3(b_1,v_1,dt_next)
      dt = minval([dt_next,dt_max])
   ELSE
      CALL ralston2(b_1,v_1,dt_next)
      dt = minval([dt_next,dt_max])
   ENDIF
   
   itime=itime+1 
   IF(mype==0.and.verbose) WRITE(*,*) "itime: ",itime
   IF(mype==0.and.verbose) WRITE(*,*) "dt: ",dt
   time=time+dt
   !IF(mype==0.and.dt_output) WRITE(*,*) "Before adapt: dt_max,dt",dt_max,dt
   !IF(adapt_dt_nl) CALL adapt_dt 
   !IF(mype==0.and.dt_output) WRITE(*,*) "After adapt: dt_max,dt",dt_max,dt
   time_tot = MPI_WTIME()
   IF((time_tot-time_start).gt.REAL(max_walltime)) THEN
   !  IF(mype==0) WRITE(*,*) "Maximum wall time exceeded.", current_wallclock, max_walltime
     continue_run=.false. 
   ENDIF
   IF (dt < 1.0E-5) then 
     WRITE(*,*) "dt too small to proceed" ,dt
     continue_run=.false.
   ENDIF
   !END IF

   ! print *, "Time Step itime",itime,dt
END DO

if (mype.eq.0.and.verbose) print *, "Run stopped"
if (force_turbulence) CALL finalize_force
if (mype.eq.0.and.verbose) print *, "Force deallocated"
CALL diag

  CALL DEALLOCATE_STEPS

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

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
    
 REAL(C_DOUBLE) :: dt_new
 REAL(C_DOUBLE) :: dt_new1,dt_new2,dt_new3,dt_new4
 REAL(C_DOUBLE) :: nmhc

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

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  INTEGER :: i,j,k,l,h
  COMPLEX(C_DOUBLE_COMPLEX) :: div_v, div_b
  REAL(C_DOUBLE) :: k2
  COMPLEX(C_DOUBLE_COMPLEX) :: exb(0:2),exv(0:2)

 div_v = 0.0 +i_complex*0.0
 div_b = 0.0 +i_complex*0.0
 k2=0.0
 zero=0.0

 if (mype.eq.0) then
    exb = b_in(1,1,1,:)
    exv = v_in(1,1,1,:)
 endif 

 DO i=cstart(1),cend(1)
   DO j=cstart(2),cend(2)
     DO k=cstart(3),cend(3)
        k2 = kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2
        div_v = kxgrid(i)*v_in(i,j,k,0) + kygrid(j)*v_in(i,j,k,1) + kzgrid(k)*v_in(i,j,k,2)
        div_b = kxgrid(i)*b_in(i,j,k,0) + kygrid(j)*b_in(i,j,k,1) + kzgrid(k)*b_in(i,j,k,2)

        v_in(i,j,k,0) = v_in(i,j,k,0) - div_v*kxgrid(i)/k2
        v_in(i,j,k,1) = v_in(i,j,k,1) - div_v*kygrid(j)/k2
        v_in(i,j,k,2) = v_in(i,j,k,2) - div_v*kzgrid(k)/k2
        
        ! The b equation is a curl, so we don't need to remove div b (except at start)
        b_in(i,j,k,0) = b_in(i,j,k,0) - div_b*kxgrid(i)/k2
        b_in(i,j,k,1) = b_in(i,j,k,1) - div_b*kygrid(j)/k2
        b_in(i,j,k,2) = b_in(i,j,k,2) - div_b*kzgrid(k)/k2
     
     ENDDO
   ENDDO
ENDDO

if (mype.eq.0) then
   b_in(1,1,1,:) = exb
   v_in(1,1,1,:) = exv
endif

if (verbose.and.(mype.eq.0)) print *,'Divergence Removed'

END SUBROUTINE remove_div

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                    get_rhs                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs(b_in,v_in, rhs_out_b,rhs_out_v,nmhc,ndt)
  
  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(out) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(out) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  REAL(C_DOUBLE), INTENT(out) :: nmhc
  REAL(C_DOUBLE) :: ndt
  REAL :: sttime,lintime,nltime,disstime,forcetime
  
  INTEGER :: i,j,k

  IF (test_ho) THEN
     CALL get_rhs_test(b_in,v_in,rhs_out_b,rhs_out_v)
     ndt = dt_max
  ELSE

     if (verbose) print *, "Pre Linear","Mype",mype,"MaxVal v_1",maxval(abs(v_in))
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if (timer.and.(mype.eq.0)) sttime = MPI_WTIME()
     CALL get_rhs_lin(b_in,v_in,rhs_out_b, rhs_out_v,0)
     if (timer.and.(mype.eq.0)) lintime = MPI_WTIME()
     if (timer.and.(mype.eq.0)) print *, "Linear Time",lintime-sttime
     
     !IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v)
     if (timer.and.(mype.eq.0)) sttime = MPI_WTIME()
     if (verbose) print *, "Pre NL","Mype",mype,"MaxVal k1",maxval(abs(rhs_out_v))
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
     IF(nonlinear) CALL get_rhs_nl(b_in, v_in,rhs_out_b,rhs_out_v,ndt)
     if (timer.and.(mype.eq.0)) nltime = MPI_WTIME()
     if (timer.and.(mype.eq.0)) print *, "Nonlinear Time",nltime-sttime
     
      IF (.not.(actual_nonlinear)) ndt = dt_max
     if (verbose.and.(mype.eq.0)) print *, ndt

     nmhc = 0.0
     if (verbose) print *, "Pre MHC","Mype",mype,"MaxVal k1",maxval(abs(rhs_out_v))
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     IF (mhc) CALL getmhcrk(b_in,v_in,nmhc)

     ! Add forcing
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     IF (mod(intorder,20).ne.0) THEN
        if (timer.and.(mype.eq.0)) sttime = MPI_WTIME()
        if (verbose) print *, "Pre Diss","Mype",mype,"MaxVal k1",maxval(abs(rhs_out_v))
        CALL get_rhs_diss(b_in,v_in,rhs_out_b,rhs_out_v)
        if (timer.and.(mype.eq.0)) disstime =MPI_WTIME()
        if (timer.and.(mype.eq.0)) print *, "Dissipation Time",disstime-sttime
        if (timer.and.(mype.eq.0)) sttime = MPI_WTIME()
        if (verbose) print *, "Pre Force","Mype",mype,"MaxVal k1",maxval(abs(rhs_out_v))
        IF (force_turbulence) CALL get_rhs_force(rhs_out_b, rhs_out_v)
        if (timer.and.(mype.eq.0)) forcetime = MPI_WTIME()
        if (timer.and.(mype.eq.0)) print *, "Forcing Time",forcetime-sttime
     ENDIF
     
     if (verbose) print *, "Pre Div","Mype",mype,"MaxVal k1",maxval(abs(rhs_out_v))
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL remove_div(rhs_out_b,rhs_out_v)
     
     if (verbose.and.(mype.eq.0)) print *,'RHS found'

  ENDIF
  rkstage = rkstage + 1
  

END SUBROUTINE get_rhs

SUBROUTINE getmhcrk(b_in,v_in,nmhc)
  
  IMPLICIT NONE
  
  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  
  REAL(C_DOUBLE), intent(out) :: nmhc
  REAL(C_DOUBLE) :: nmhcsm(1:2),nmhcs(1:2)
  INTEGER :: i,j,k
  
 nmhcsm = 0.0

 DO j = cstart(2),cend(2)
    DO k = cstart(3),cend(3)
       nmhcsm(2) = nmhcsm(2) - 1.0 * real(b_in(1,j,k,0)*conjg(v_in(1,j,k,1))-b_in(1,j,k,1)*conjg(v_in(1,j,k,0)))
       DO i = xst,cend(1)
          nmhcsm(1) = nmhcsm(1) -2.0 * real(b_in(i,j,k,0)*conjg(v_in(i,j,k,1))-b_in(i,j,k,1)*conjg(v_in(i,j,k,0)))
       ENDDO
    ENDDO
 ENDDO

 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
 CALL MPI_ALLREDUCE(nmhcsm,nmhcs,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
 
 nmhc = sum(nmhcs)*(16.0*pi**3)
 
 if (verbose.and.(mype.eq.0)) print *, "Through MHC"
   
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

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  REAL(C_DOUBLE) :: dt_new1,dt_new2,dt_new3,dt_new4,dt_new5,dt_new6,dt_new7
  REAL(C_DOUBLE) :: nmhc1,nmhc2,nmhc3,nmhc4,nmhc5,nmhc6,nmhc7,next_corr
  REAL(C_DOUBLE) :: errm,err
  
  dt = dt_max

  ! print *, "dt start", dt
  ! print *, "dt_max", dt_max
  !! First step

  ! Use First Same As Last property to save time for higher steps
  if (itime.eq.0) CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,dt_new1)
     
  !! Second step
  bk1s = b_in+dt*bk1*1.0/5.0
  vk1s = v_in+dt*bk1*1.0/5.0
  CALL get_rhs(bk1s,vk1s,bk2,vk2,nmhc2,dt_new2)
     
  !! Third step
  bk1s = b_in+dt*(3.0*bk1/40.0+9.0*bk2/40.0)
  vk1s = v_in+dt*(3.0*vk1/40.0+9.0*vk2/40.0)
  CALL get_rhs(bk1s,vk1s,bk3,vk3,nmhc3,dt_new3)
     
  !! Fourth step
  bk1s = b_in+dt*(44.0/45.0*bk1-56.0/15.0*bk2+32.0/9.0*bk3)
  vk1s = v_in+dt*(44.0/45.0*vk1-56.0/15.0*vk2+32.0/9.0*vk3)
  CALL get_rhs(bk1s,vk1s,bk4,vk4,nmhc4,dt_new4)
     
  !! Fifth step
  bk1s = b_in+dt*(19372.0/6561.0*bk1-25360.0/2187.0*bk2+64448.0/6561.0*bk3-212.0/729.0*bk4)
  vk1s = v_in+dt*(19372.0/6561.0*vk1-25360.0/2187.0*vk2+64448.0/6561.0*vk3-212.0/729.0*vk4)
  CALL get_rhs(bk1s,vk1s,bk5,vk5,nmhc5,dt_new5)
     
  !! Sixth step
  bk1s = b_in+dt*(9017.0/3168.0*bk1 - 355.0/33.0*bk2 + 46732.0/5247.0*bk3 + 49.0/176.0*bk4-5103.0/18656.0*bk5)
  vk1s = v_in+dt*(9017.0/3168.0*vk1 - 355.0/33.0*vk2 + 46732.0/5247.0*vk3 + 49.0/176.0*vk4-5103.0/18656.0*vk5)
  CALL get_rhs(bk1s,vk1s,bk6,vk6,nmhc6,dt_new6)

  
  !print *, "6",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  b_2 = b_in + dt*(35.0/384.0*bk1+500.0/1113.0*bk3+125.0/192.0*bk4&
       -2187.0/6784.0*bk5+11.0/84.0*bk6)

  !print *, "b2",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
  v_2 = v_in + dt*(35.0/384.0*vk1+500.0/1113.0*vk3+125.0/192.0*vk4&
       -2187.0/6784.0*vk5+11.0/84.0*vk6)

  !print *, "v2",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  next_corr = dt*(35.0/384.0*nmhc1+500.0/1113.0*nmhc3+125.0/192.0*nmhc4&
       -2187.0/6784.0*nmhc5+11.0/84.0*nmhc6)

  !print *, mype, "6th stage"
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
  !! Seventh step
  CALL get_rhs(b_2,v_2,bk2,vk2,nmhc7,dt_new7)

  !print *, "last stage",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  b_3 = b_in + dt*(5179.0/57600.0*bk1+7571.0/16695.0*bk3+393.0/640.0*bk4&
       -92097.0/339200.0*bk5+187.0/2100.0*bk6+bk2/40.0)

  !print *, "b3",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  v_3 = v_in + dt*(5179.0/57600.0*vk1+7571.0/16695.0*vk3+393.0/640.0*vk4&
       -92097.0/339200.0*vk5+187.0/2100.0*vk6+vk2/40.0)

  !print *, "v3",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !print *, "Through b3",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  errm = max(maxval(abs(b_2-b_3)),maxval(abs(v_2-v_3)))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(errm,err,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

  !dt = min(dt_new1,dt_new2,dt_new3,dt_new4,dt_new5,dt_new6,dt * (10.0**(-5.0)/err)**(1.0/6.0))
  
  b_in = b_2
  v_in = v_2
  mhelcorr = mhelcorr + next_corr
  !print *, "Through b_in",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  bk1 = bk2
  vk1 = vk2

  !print *, "Through DORPI",mype
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE dorpi547M

SUBROUTINE ralston2(b_in,v_in,dt_new)

  ! Uses Ralston's method to integrate system
  ! Butcher tableau from Wikipedia
  ! 0  |
  ! 2/3| 2/3
  ! ___|____
  !      1/4 3/4 
  
  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  REAL(C_DOUBLE), INTENT(OUT) :: dt_new
  REAL(C_DOUBLE) :: dt_new1 = 0,dt_new2 = 0
  REAL(C_DOUBLE) :: nmhc1,nmhc2
  
  CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,dt_new1)

  if (verbose) print *, "Post First Stage","Mype",mype,"Maxval k1",maxval(abs(vk1))
  b_2 = b_in + (1.0/4.0)*bk1*dt
  v_2 = v_in + (1.0/4.0)*vk1*dt
  if (verbose) print *, "Pre Second Stage","Mype",mype,"MaxVal v_1",maxval(abs(v_2))
  if (mhc) mhelcorr = mhelcorr + (1.0/4.0)*nmhc1*dt
  if (verbose) print *, "Through MHC"

  bk1s = b_in+(2.0/3.0)*bk1*dt
  vk1s = v_in+(2.0/3.0)*vk1*dt

  if (verbose) print *, "Pre Second Stage","Mype",mype,"MaxVal v_1",maxval(abs(vk1s))
  
  CALL get_rhs(bk1s,vk1s,bk2,vk2,nmhc2,dt_new2)
  
  b_in = b_2 + (3.0/4.0) * bk2 * dt
  v_in = v_2 + (3.0/4.0) * vk2 * dt
  
  if (mhc) mhelcorr = mhelcorr + (3.0/4.0) * nmhc2 * dt

  dt_new = minval([dt_new1,dt_new2])

END SUBROUTINE ralston2

SUBROUTINE ralston3(b_in,v_in,dt_new)

  ! Uses Ralston's third order minimum error bound method
  ! Butcher tableau from Ralston 1962 Runge-Kutta Methods with Minimum Error Bounds

  !     |
  ! 1/2 | 1/2
  ! 3/4 |  0   3/4
  !____________________
  !     | 2/9  1/3  4/9

  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  
  REAL(C_DOUBLE), INTENT(OUT) :: dt_new
  REAL(C_DOUBLE) :: dt_new1,dt_new2,dt_new3
  REAL(C_DOUBLE) :: nmhc1,nmhc2,nmhc3

  CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,dt_new1)
  b_2 = b_in + (2.0/9.0)*bk1*dt
  v_2 = v_in + (2.0/9.0)*vk1*dt
  if (verbose)   print *, "First Stage"
  if (mhc) mhelcorr = mhelcorr + (2.0/9.0)*nmhc1*dt

  bk1s = b_in+(1.0/2.0)*bk1*dt
  vk1s = v_in+(1.0/2.0)*vk1*dt

  CALL get_rhs(bk1s,vk1s,bk2,vk2,nmhc2,dt_new2)
  b_2 = b_2 + (1.0/3.0) * bk2 * dt
  v_2 = v_2 + (1.0/3.0) * vk2 * dt
  if (mhc) mhelcorr = mhelcorr + (1.0/3.0) * nmhc2 * dt

  bk1s = b_in+(3.0/4.0)*bk2*dt
  vk1s = v_in+(3.0/4.0)*vk2*dt

  CALL get_rhs(bk1s,vk1s,bk1,vk1,nmhc3,dt_new3)
  b_in = b_2 + (4.0/9.0) * bk1 * dt
  v_in = v_2 + (4.0/9.0) * vk1 * dt
  if (mhc) mhelcorr = mhelcorr + (4.0/9.0) * nmhc3 * dt
  
  dt_new = minval([dt_new1,dt_new2,dt_new3])

END SUBROUTINE ralston3

SUBROUTINE GAUSS4(b_in,v_in,dt_new)

  ! Use the order 4 Gauss Legendre collocation method for time integration
  ! k1 = f(y + a11 * dt * k1 + a12 * dt * k2)
  ! k2 = f(y + a21 * dt * k1 + a22 * dt * k2)
  ! y1 = y + 0.5 * dt * k1 + 0.5 * dt * k2
  ! This method relies on a fixed point iteration solution for the implicit RK stages
  
  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  
  REAL(C_DOUBLE), INTENT(OUT) :: dt_new
  
  REAL(C_DOUBLE) :: a11,a12,a21,a22
  REAL(C_DOUBLE) :: nmhc1,nmhc2,nmhc1s,nmhc2s
  REAL(C_DOUBLE) :: ndt,maxdev,maxdevm
  INTEGER :: solvloop

  a11 = 1.0/4.0
  a22 = 1.0/4.0
  a12 = 1.0/4.0 - sqrt(3.0)/6.0
  a21 = 1.0/4.0 + sqrt(3.0)/6.0

  ! Initial derivative guesses
  CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,ndt)
  CALL get_rhs(b_in,v_in,bk2,vk2,nmhc2,ndt)

  ! Start iterating

  maxdev = 1.0
  solvloop = 0

  DO WHILE ((solvloop.lt.1000).and.(maxdev.gt.10.0**(-16.0)))

     if (solvloop.eq.0) then

        b_2 = b_in+a11*dt*bk1+a12*dt*bk2
        v_2 = v_in+a11*dt*vk1+a12*dt*vk2

        CALL get_rhs(b_2,v_2,bk1s,vk1s,nmhc1s,ndt)

        b_2 = b_in+a21*dt*bk1+a22*dt*bk2
        v_2 = v_in+a21*dt*vk1+a22*dt*vk2

        CALL get_rhs(b_2,v_2,bk2s,vk2s,nmhc2s,ndt)

     endif

     bk1 = bk1s
     vk1 = vk1s

     bk2 = bk2s
     vk2 = vk2s

     b_2 = b_in+a11*dt*bk1+a12*dt*bk2
     v_2 = v_in+a11*dt*vk1+a12*dt*vk2

     CALL get_rhs(b_2,v_2,bk1s,vk1s,nmhc1s,ndt)
     
     b_2 = b_in+a21*dt*bk1+a22*dt*bk2
     v_2 = v_in+a21*dt*vk1+a22*dt*vk2
     
     CALL get_rhs(b_2,v_2,bk2s,vk2s,nmhc2s,ndt)
     
     ! Fixed point solution
     
     maxdevm = max(maxval(abs(bk1-bk1s)),&
          maxval(abs(vk1-vk1s)),&
          maxval(abs(bk2-bk2s)),&
          maxval(abs(vk2-vk2s)))
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_ALLREDUCE(maxdevm,maxdev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
     
     solvloop = solvloop + 1
     
     if (solvloop.eq.999.and.(mype.eq.0)) print *, "Failed to Converge After 1000 iterations Maxdev ",maxdev

  ENDDO

  if (itime.eq.10.and.(mype.eq.0)) print *, "Number of Iterations Needed itime = 10 ",solvloop

  ! Update the fields
  
  b_in = b_in + (1.0/2.0) * dt * bk1s + (1.0/2.0) * dt * bk2s
  v_in = v_in + (1.0/2.0) * dt * vk1s + (1.0/2.0) * dt * vk2s
  mhelcorr = mhelcorr + (1.0/2.0) * dt * nmhc1s + (1.0/2.0) * dt * nmhc2s
  dt_new = ndt

END SUBROUTINE GAUSS4


SUBROUTINE GAUSS2(b_in,v_in,dt_new)

  ! Uses the implicit midpoint method for time integration
  ! Like other Gauss methods this conserves quadratic invariants
  ! This method relies on a fixed point iteration solution for the implicit RK stages                                                                                                                                                                                    

  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  REAL(C_DOUBLE), INTENT(OUT) :: dt_new

  REAL(C_DOUBLE) :: a11
  REAL(C_DOUBLE) :: nmhc1,nmhc1s
  REAL(C_DOUBLE) :: ndt,maxdev,maxdevm
  INTEGER :: solvloop

  a11 = 1.0/2.0

  ! Initial guess
  CALL get_rhs(b_in,v_in,bk1,vk1,nmhc1,ndt)
  ! Iterate
  bk2 = b_in+a11*dt*bk1
  vk2 = v_in+a11*dt*vk1
  CALL get_rhs(bk2,vk2,bk1s,vk1s,nmhc1s,ndt)

  maxdev = 1.0
  solvloop = 0

  DO WHILE ((solvloop.lt.1000).and.(maxdev.gt.10.0**(-16.0)))
     ! Check for now to see if the fixed point iteration converges

     if (verbose) print *, mype,"Iteration ",solvloop,"Discrepancy ",maxdev

     if (verbose) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     bk1 = bk1s 
     vk1 = vk1s

     bk2 = b_in+a11*dt*bk1
     vk2 = v_in+a11*dt*vk1
     if (verbose) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if (verbose) print *, mype,"Through Iteration ",solvloop, "Assignment"

     CALL get_rhs(bk2,vk2,bk1s,vk1s,nmhc1s,ndt)

     if (verbose) print *, mype,"Through Iteration ",solvloop,"Iteration"
     solvloop = solvloop + 1

     if (verbose) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     maxdevm = max(maxval(abs(bk1-bk1s)),maxval(abs(vk1-vk1s)))
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_ALLREDUCE(maxdevm,maxdev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)

     if (verbose) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     if (solvloop.eq.999.and.(mype.eq.0)) print *, "Failed to Converge After 1000 iterations, Maxdev ",maxdev

  ENDDO

  if (itime.eq.10.and.(mype.eq.0)) print *, "Number of Iterations Needed itime = 10 ",solvloop

  ! Update the fields
  
  b_in = b_in + dt * bk1s
  v_in = v_in + dt * vk1s
  mhelcorr = mhelcorr + dt * nmhc1s 
  dt_new = ndt

END SUBROUTINE GAUSS2

SUBROUTINE SPLIT2(b_in,v_in)

  ! Uses second order Strang splitting
  
  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  REAL(C_DOUBLE) :: ndt
  REAL :: dt2

  dt2 = dt
!   if (.not.present(dtm)) dtm = dt_max
  
  dt = 0.5 * dt2
  if (verbose) print *, "Set time step",mype
  CALL get_rhs_diss2(b_in,v_in)
  if (verbose) print *, "Through dissipation 1",mype
  if (force_turbulence) CALL get_rhs_force(bk1,vk1)
  b_in = b_in + dt * bk1
  v_in = v_in + dt * vk1
  if (verbose) print *, "Through force 1",mype

  dt = dt2
  if (intorder.eq.20) CALL GAUSS2(b_in,v_in,ndt)
  if (intorder.eq.40) CALL GAUSS4(b_in,v_in,ndt)
  if (verbose) print *, "Through ideal",mype

  dt = 0.5 * dt2
  bk1 = 0.0
  vk1 = 0.0
  if (force_turbulence) CALL get_rhs_force(bk1,vk1)

  b_in = b_in + dt * bk1
  v_in = v_in + dt * vk1
  if (verbose.and.mype.eq.0) print *, "Through force 2"
  
  CALL get_rhs_diss2(b_in,v_in)
  if (verbose.and.mype.eq.0) print *, "Through dissipation 2"

  dt = dt2

END SUBROUTINE SPLIT2

SUBROUTINE SPLIT4(b_in,v_in)

  ! Uses 4th order splitting methods as described by H. Yoshida 1990, Hairer et al 2006, Suzuki, etc.

  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
    
  REAL(C_DOUBLE) :: dt2

  dt2 = dt

  dt = - (2.0**(1.0/3.0))/(2.0-2.0**(1.0/3.0)) * dt2
  CALL SPLIT2(b_1,v_1)

  dt = 1.0/(2.0-2.0**(1.0/3.0)) * dt2
  CALL SPLIT2(b_1,v_1)

  dt = - (2.0**(1.0/3.0))/(2.0-2.0**(1.0/3.0)) * dt2
  CALL SPLIT2(b_1,v_1)

  dt = dt2

END SUBROUTINE SPLIT4

SUBROUTINE ALLOCATE_STEPS

  IMPLICIT NONE

  
  ALLOCATE(b_2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(v_2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  ALLOCATE(bk1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(vk1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  ALLOCATE(bk2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(vk2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  ALLOCATE(bk1s(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(vk1s(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  ALLOCATE(bk2s(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(vk2s(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  if (intorder.eq.5) then

     ALLOCATE(b_3(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     ALLOCATE(v_3(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     
     ALLOCATE(bk3(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     ALLOCATE(vk3(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     
     ALLOCATE(bk4(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     ALLOCATE(vk4(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     
     ALLOCATE(bk5(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     ALLOCATE(vk5(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     
     ALLOCATE(bk6(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     ALLOCATE(vk6(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
     
  endif
  
END SUBROUTINE ALLOCATE_STEPS

SUBROUTINE DEALLOCATE_STEPS

  IMPLICIT NONE

  if (allocated(b_2)) DEALLOCATE(b_2)
  if (allocated(v_2)) DEALLOCATE(v_2)

  if (allocated(bk1)) DEALLOCATE(bk1)
  if (allocated(vk1)) DEALLOCATE(vk1)

  if (allocated(bk2)) DEALLOCATE(bk2)
  if (allocated(vk2)) DEALLOCATE(vk2)

  if (allocated(bk1s)) DEALLOCATE(bk1s)
  if (allocated(vk1s)) DEALLOCATE(vk1s)
  if (allocated(bk2s)) DEALLOCATE(bk2s)
  if (allocated(vk2s)) DEALLOCATE(vk2s)

  if (allocated(b_3)) DEALLOCATE(b_3)
  if (allocated(v_3)) DEALLOCATE(v_3)

  if (allocated(bk3)) DEALLOCATE(bk3)
  if (allocated(vk3)) DEALLOCATE(vk3)

  if (allocated(bk4)) DEALLOCATE(bk4)
  if (allocated(vk4)) DEALLOCATE(vk4)

  if (allocated(bk5)) DEALLOCATE(bk5)
  if (allocated(vk5)) DEALLOCATE(vk5)

  if (allocated(bk6)) DEALLOCATE(bk6)
  if (allocated(vk6)) DEALLOCATE(vk6)

END SUBROUTINE DEALLOCATE_STEPS

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
