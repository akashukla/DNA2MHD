!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                         cc_initial_condition.f90                          !!
!!                                                                           !!
!!  initial_condition                                                        !!
!!  get_real                                                                 !!
!!  check_real                                                               !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                            initial_condition                              !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Provides initial condition
!!  To DO: implement a better i.c.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initial_condition(which_init0)
  USE par_mod
  USE mtrandom
  IMPLICIT NONE

  INTEGER :: i,j,k, l
  REAL :: kfactor,err
  REAL :: kmags(0:nkx0-1,0:nky0-1,1:nkz0-1),divratio(0:nkx0-1,0:nky0-1,1:nkz0-1)
  REAL :: s1, s2
  !REAL :: init_prefactor
  COMPLEX :: phase
  REAL :: phase1,phase2,kspect
  CHARACTER(len=40), INTENT(in) :: which_init0
  CHARACTER(len=40) :: which_init
  INTEGER, DIMENSION(:), ALLOCATABLE :: rseed
  INTEGER :: rseed_size
  REAL :: zerocmplx

  zerocmplx=cmplx(0.0,0.0)
   
  !g_1(:,:,:,:,:,:)=cmplx(0.0,0.0)
  b_1(:,:,:,:)=cmplx(0.0,0.0)
  v_1(:,:,:,:)=cmplx(0.0,0.0)

  DO i = 0,nkx0-1
    DO j = 0,nky0-1
      DO k = 1,nkz0-1
        kmags(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2)
        divratio(i,j,k) = (kxgrid(i)+kygrid(j))/kzgrid(k)
      END DO
    END DO
  END DO
  s1 = sum(kmags(:,:,1:nkz0-1) ** (-1.0*init_kolm))
  s2 = sum((divratio(:,:,1:nkz0-1) ** 2) * kmags(:,:,1:nkz0-1) ** (-1.0 * init_kolm))

  !init_prefactor=0.001
!Default Initialization
!      CALL RANDOM_SEED
      CALL RANDOM_SEED(SIZE=rseed_size)
      ALLOCATE(rseed(rseed_size))
      rseed(:) = 1
      CALL RANDOM_SEED(PUT=rseed)
      DEALLOCATE(rseed)
 
      DO i=kxinit_min,kxinit_max-1
       DO j=kyinit_min,kyinit_max-1
        DO k=kzinit_min,kzinit_max-1
         !DO l=0,2
         IF(kzgrid(k).eq.zerocmplx) THEN
             b_1(i,j,k,0)=cmplx(0.0,0.0)
             b_1(i,j,k,1)=cmplx(0.0,0.0)
             b_1(i,j,k,2)=cmplx(0.0,0.0)
             v_1(i,j,k,0)=cmplx(0.0,0.0)
             v_1(i,j,k,1)=cmplx(0.0,0.0)
             v_1(i,j,k,2)=cmplx(0.0,0.0)
         ELSE IF (.not.(enone)) THEN
             b_1(i,j,k,0)=init_amp_bx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))* 1/(kmags(i,j,k)**(init_kolm/2.0))
             b_1(i,j,k,1)=init_amp_by*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kmags(i,j,k)**(init_kolm/2.0))
             b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
             !b_1(i,j,k,2)=init_amp_bz
             v_1(i,j,k,0)=init_amp_vx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kmags(i,j,k)**(init_kolm/2.0))
             v_1(i,j,k,1)=init_amp_vy*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kmags(i,j,k)**(init_kolm/2.0))
             !v_1(i,j,k,2)=init_amp_vz
             v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
         ELSE
             b_1(i,j,k,0)= sqrt(init_amp_bx/(8 * pi**3 * (2*s1 + s2)))/(kmags(i,j,k)**(init_kolm/2.0)) 
             !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
             b_1(i,j,k,1)= sqrt(init_amp_bx/(8 * pi**3 * (2*s1 + s2)))/(kmags(i,j,k)**(init_kolm/2.0))
             !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
             b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
             !b_1(i,j,k,2)=init_amp_bz              
             v_1(i,j,k,0)= sqrt(init_amp_bx/(8 * pi**3 * (2*s1 + s2)))/(kmags(i,j,k)**(init_kolm/2.0))
             !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
             v_1(i,j,k,1)= sqrt(init_amp_bx/(8 * pi**3 * (2*s1 + s2)))/(kmags(i,j,k)**(init_kolm/2.0))
             !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
             !v_1(i,j,k,2)=init_amp_vz 
             v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
        END IF
        END DO
       END DO
      END DO
      if (nv) b_1(:,:,:,:) = cmplx(0.0,0.0)
      gpsi(:,:,:,:) = cmplx(0.0,0.0)

! Only use default for now
!  which_init=which_init0 
!  IF(checkpoint_read) which_init='checkpoint'
!
!  IF(nonlinear) THEN
!    IF(trim(which_init)=='low_k') THEN
!      IF(mype==0) WRITE(*,*) "low_k initial condition"
!      DO i=0,nkx0-1
!       DO j=0,nky0-1
!        DO k=0,nkz0-1
!         IF((lv1==0).and.(i.le.2).and.(j.le.2).and.(k.le.2)) g_1(i,j,k,0,:,:)=cmplx(0.001,0.001)
!         IF(lv1.le.1.and.lv2.ge.1.and.(i.le.2).and.(j.le.2).and.(k.le.2)) g_1(i,j,k,1,:,:)=cmplx(0.001,0.001)
!         IF(lv1.le.2.and.lv2.ge.2.and.(i.le.2).and.(j.le.2).and.(k.le.2)) g_1(i,j,k,2,:,:)=cmplx(0.001,0.001)
!         IF(lv1.le.3.and.lv2.ge.3.and.(i.le.2).and.(j.le.2).and.(k.le.2)) g_1(i,j,k,3,:,:)=cmplx(0.001,0.001)
!        END DO
!       END DO
!      END DO
!    ELSE IF(trim(which_init)=='DEFAULT') THEN
!      IF(mype==0.and.verbose) WRITE(*,*) "DEFAULT initial condition"
!!      CALL RANDOM_SEED
!      CALL RANDOM_SEED(SIZE=rseed_size)
!      ALLOCATE(rseed(rseed_size))
!      rseed(:) = 1
!      CALL RANDOM_SEED(PUT=rseed)
!      DEALLOCATE(rseed)
! 
!      DO i=0,nkx0-1
!       DO j=0,nky0-1
!        IF((abs(kxgrid(i)+kygrid(j)).gt.epsilon(1.0))) THEN  
!           kfactor=0.0
!           kfactor=abs(kxgrid(i)/kxmin)
!           kfactor=kfactor+abs(kygrid(j)/kymin)
!           kfactor=1.0/kfactor
!           !Some stuff to mix the phases:
!           CALL RANDOM_NUMBER(phase1)
!           CALL RANDOM_NUMBER(phase2)
!           phase=cmplx(phase1,phase2)
!           !IF(lv1==0) g_1(i,j,0,0)=                 init_prefactor*phase*kfactor
!           !IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,0,1)=  init_prefactor*phase*kfactor
!           !IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,0,2)=  init_prefactor*phase*kfactor
!           !IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,0,3)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,1,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,1,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,1,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,1,3,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,2,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,2,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,2,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,2,3,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,3,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,3,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,3,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,3,3,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,nkz0-1,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,nkz0-1,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,nkz0-1,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,nkz0-1,3,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,nkz0-2,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,nkz0-2,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,nkz0-2,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,nkz0-2,3,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1==0) g_1(i,j,nkz0-3,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,nkz0-3,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,nkz0-3,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,nkz0-3,3,:,:)=  init_prefactor*phase*kfactor
!        END IF
!       END DO
!      END DO
!
!      CALL get_real(g_1)
!      CALL check_real(g_1,err)
!      IF(mype==0.and.verbose) WRITE(*,*) "Done with initial condition."
!      !WRITE(*,*) "mype,check REAL error",mype,err
!
!
!    ELSE IF(trim(which_init)=='SINGLE_KZ') THEN
!      IF(mype==0.and.verbose) WRITE(*,*) "SINGLE_KZ initial condition"
!!      CALL RANDOM_SEED
!      CALL RANDOM_SEED(SIZE=rseed_size)
!      ALLOCATE(rseed(rseed_size))
!      rseed(:) = 1
!      CALL RANDOM_SEED(PUT=rseed)
!      DEALLOCATE(rseed)
! 
!      DO i=0,nkx0-1
!       DO j=0,nky0-1
!        IF((abs(kxgrid(i)+kygrid(j)).gt.epsilon(1.0))) THEN  
!           kfactor=0.0
!           kfactor=abs(kxgrid(i)/kxmin)
!           kfactor=kfactor+abs(kygrid(j)/kymin)
!           kfactor=1.0/kfactor
!           !Some stuff to mix the phases:
!           CALL RANDOM_NUMBER(phase1)
!           CALL RANDOM_NUMBER(phase2)
!           phase=cmplx(phase1,phase2)
!           IF(lv1==0) g_1(i,j,0,0,:,:)=                 init_prefactor*phase*kfactor
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,0,1,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,0,2,:,:)=  init_prefactor*phase*kfactor
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,0,3,:,:)=  init_prefactor*phase*kfactor
!        END IF
!       END DO
!      END DO
!
!      CALL get_real(g_1)
!      CALL check_real(g_1,err)
!      IF(mype==0.and.verbose) WRITE(*,*) "Done with initial condition."
!      !WRITE(*,*) "mype,check REAL error",mype,err
!
!
!    ELSE IF(trim(which_init)=='old_default') THEN
!
!     IF(mype==0) WRITE(*,*) "Initial condition:",trim(which_init)
!      CALL RANDOM_SEED(SIZE=rseed_size)
!      ALLOCATE(rseed(rseed_size))
!      rseed(:) = 1
!      CALL RANDOM_SEED(PUT=rseed)
!      DEALLOCATE(rseed)
!      !DO i=1,nkx0-1
!      ! DO j=0,nky0-1
!      !  DO k=0,nkz0-1
!      DO i=0,nkx0-1
!       DO j=0,nky0-1
!        DO k=0,nkz0-1
!
!         IF(abs(kxgrid(i)+kygrid(j)+kygrid(k)).gt.epsilon(1.0)) THEN  
!
!           kfactor=kxgrid(i)
!           kfactor=kfactor+kygrid(j)
!           kfactor=kfactor+kzgrid(k)
!           IF(kfactor.gt.epsilon(1.0)) kfactor=1.0/kfactor
!           !Some stuff to mix the phases:
!           CALL RANDOM_NUMBER(phase1)
!           CALL RANDOM_NUMBER(phase2)
!           phase=cmplx(phase1,phase2)
!           IF(lv1==0) THEN
!               !g_1(i,j,k,0)= REAL(nv0-0)/REAL(nv0)*init_prefactor*phase*kfactor
!               g_1(i,j,k,0,0,:)= cmplx(REAL(i),REAL(j))*phase
!               !WRITE(*,*) i,j,k,0,g_1(i,j,k,0)
!           END IF
!           IF(lv1.le.1.and.lv2.ge.1) THEN 
!               g_1(i,j,k,1,:,:)= (cmplx(1.0,1.0)+cmplx(REAL(i),REAL(j)))*phase
!               !g_1(i,j,k,1)= 5.0*REAL(nv0-1)/REAL(nv0)*init_prefactor*phase*kfactor
!               !WRITE(*,*) i,j,k,1,g_1(i,j,k,1)
!           END IF
!                 
!
!         END IF
!
!        END DO
!       END DO
!      END DO
!
!      CALL get_real(g_1)
!      CALL check_real(g_1,err)
!      !WRITE(*,*) "mype,check REAL error",mype,err
!
!    ELSE IF(trim(which_init)=='ic_test') THEN
!         i=1;j=1;k=0
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=1;j=1;k=1
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=nkx0-1;j=nky0-1;k=0
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=nkx0-1;j=nky0-1;k=nkz0-1
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=nkx0/4;j=nky0/4;k=0
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=nkx0/4;j=nky0/4;k=nkz0/4
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=3*nkx0/4;j=3*nky0/4;k=0
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!         i=3*nkx0/4;j=3*nky0/4;k=3*nkz0/4
!           IF(lv1==0) g_1(i,j,k,0,:,:)=                 cmplx(1.0,1.0)
!           IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=  cmplx(1.0,1.0)
!           IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=  cmplx(1.0,1.0)
!    ELSE IF(trim(which_init)=='high_amp') THEN
!      g_1=cmplx(0.0,0.0)
!      DO i=1,nkx0-1
!        DO j=1,nky0-1
!          DO k=1,nkz0-1
!          kspect=0.01*(kxgrid(i)**(-2.0)+kygrid(j)**(-2.0)+kzgrid(k)**(-2.0))
!            IF(lv1==0) g_1(i,j,k,0,:,:)=cmplx(kspect,kspect)
!            IF(lv1.le.1.and.lv2.ge.1) g_1(i,j,k,1,:,:)=cmplx(kspect,kspect)
!            IF(lv1.le.2.and.lv2.ge.2) g_1(i,j,k,2,:,:)=cmplx(kspect,kspect)
!            IF(lv1.le.3.and.lv2.ge.3) g_1(i,j,k,3,:,:)=cmplx(kspect,kspect)
!          END DO
!        END DO
!      END DO
!    ELSE IF(trim(which_init)=='simple') THEN
!      !Note: designed for kxmin=kymin=kzmin=0.1
!      g_1=cmplx(0.0,0.0)
!     !WRITE(*,*) "nl_test initial condition!!!!!"
!     IF(mype==0) THEN
!      !CALL RANDOM_SEED
!      CALL RANDOM_NUMBER(phase1)
!      CALL RANDOM_NUMBER(phase2)
!      phase=cmplx(phase1,phase2)
!
!!      i=0;j=1;k=2
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=0;j=1;k=1
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=0;j=2;k=1
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=0;j=2;k=2
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=1;j=1;k=2
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=1;j=1;k=1
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=1;j=2;k=1
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!!
!!      i=1;j=2;k=2
!!      g_1(i,j,k,0)=cmplx(REAL(i),REAL(j))*phase
!!      g_1(i,j,k,1)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!
!      i=nkx0-1;j=1;k=0
!      g_1(i,j,k,0,:,:)=cmplx(REAL(i),REAL(j))*phase
!      WRITE(*,*) i,j,k,0,g_1(i,j,k,0,:,:)
!      g_1(i,j,k,1,:,:)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!      WRITE(*,*) i,j,k,1,g_1(i,j,k,1,:,:)
!
!      i=nkx0-3;j=1;k=0
!      g_1(i,j,k,0,0,:)=cmplx(REAL(i),REAL(j))*phase
!      WRITE(*,*) i,j,k,0,g_1(i,j,k,0,:,:)
!      g_1(i,j,k,1,:,:)=(cmplx(REAL(i),REAL(j))+(1.0,1.0))*phase
!      WRITE(*,*) i,j,k,1,g_1(i,j,k,1,:,:)
!
!      CALL get_real(g_1)
!
!     END IF
!
!    ELSE IF(trim(which_init)=='cosn') THEN
!        CALL cosnoise(g_1)
!        CALL get_real(g_1)
!        CALL check_real(g_1,err)
!      IF(mype==0.and.verbose) WRITE(*,*) "Done with initial condition."
! 
!    ELSE IF(trim(which_init)=='checkpoint') THEN
!
!        IF(mype==0) WRITE(*,*) "Reading checkpoint."
!        CALL checkpoint_in
!        WRITE(*,*) "Done reading checkpoint.",mype
!        
!    ELSE
!      STOP "Invalid parameter: init_cond!"
!    END IF
!  ELSE
!    IF(lv1==0) g_1(:,:,:,0,:,:)=cmplx(0.01,0.01)
!    IF(lv1.le.1.and.lv2.ge.1) g_1(:,:,:,1,:,:)=cmplx(0.01,0.01)
!    IF(lv1.le.2.and.lv2.ge.2) g_1(:,:,:,2,:,:)=cmplx(0.01,0.01)
!    IF(lv1.le.3.and.lv2.ge.3) g_1(:,:,:,3,:,:)=cmplx(0.01,0.01)
!  END IF


END SUBROUTINE initial_condition



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                  get_real                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE get_real(g_in)
!    USE par_mod
!  
!    COMPLEX, INTENT(inout) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) 
!  
!    INTEGER :: j,k
!  
!    if(np_kz.ne.1) STOP "get_real not yet modified for kz parallelization."
!  
!    !kz=0
!    DO j=1,hky_ind
!      g_in(0,nky0-j,0,:,:,:) = conjg(g_in(0,j,0,:,:,:))
!    END DO
!  
!    !ky=0
!    DO k=1,hkz_ind
!      g_in(0,0,nkz0-k,:,:,:) = conjg(g_in(0,0,k,:,:,:))
!    END DO
!  
!    !the rest
!    DO j=1,nky0-1
!      DO k=1,hkz_ind
!        g_in(0,nky0-j,nkz0-k,:,:,:) = conjg(g_in(0,j,k,:,:,:))
!      END DO
!    END DO
!  
!    g_in(0,hky_ind+1,:,:,:,:)=cmplx(0.0,0.0)
!    IF(nkz0.ge.2) g_in(0,:,hkz_ind+1,:,:,:)=cmplx(0.0,0.0)
!    g_in(0,0,0,:,:,:)=cmplx(0.0,0.0)
!  
!  END SUBROUTINE get_real
!  
!  !SUBROUTINE cosnoise(aux_x, aux_y, aux_z, aux_amp)
!  !  SUBROUTINE cosnoise(g_in)
!  !    USE par_mod
!  !    USE mtrandom
!  !      COMPLEX, INTENT(inout) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) 
!  !      COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2) :: dens
!  !  !    REAL, INTENT(IN) :: aux_x, aux_y, aux_z, aux_amp
!  !      REAL :: AMPLITUDE
!  !      INTEGER :: i, j, k, n, h, ex, ey, es
!  !      COMPLEX, DIMENSION(0:nkx0-1, 0:nky0-1) :: noise
!  !      COMPLEX :: imag  = (0.0,1.0)
!  !  
!  !      g_in = CMPLX(0.0,0.0)
!  !      ex = 2
!  !      ey = 2
!  !      es = 1
!  !      amplitude = 1.0e-2
!  !  
!  !  !    IF (aux_x.ne.-100) ex = Abs(aux_x)
!  !  !    IF (aux_y.ne.-100) ey = Abs(aux_y)
!  !  !    IF (aux_z.ne.-100) es = Abs(aux_z)
!  !  !    IF (aux_amp.ne.-100) amplitude = Abs(aux_amp)
!  !  
!  !      call sgrnd(1)
!  !      noise=0.
!  !    
!  !    !reality condition for noise
!  !    !  DO j =0,hky_ind
!  !    !      noise(0,j)= amplitude*grnd() + amplitude*imag*grnd()
!  !    !      noise(0,nky0-j) =  conjg(noise(0,j))
!  !    !  END DO
!  !  
!  !      DO i = 0,nkx0-1
!  !          DO j = 0, nky0-1
!  !              noise(i,j)= amplitude*grnd() + amplitude*imag*grnd()
!  !          END DO
!  !      END DO
!  !  
!  !      dens=0.
!  !      DO k = lkz1, lkz2
!  !          dens(ex,0,k) = 1.0
!  !          DO j=0,nky0-1
!  !              IF (j==ey.or.j==(nky0-ey)) dens(0,j,k) = 1.0
!  !              DO i=0,nkx0-1
!  !                  dens(i,j,k)=dens(i,j,k)+noise(i,j)
!  !              END DO
!  !          END DO
!  !      END DO
!  !  
!  !      IF (.not.mu_integrated) THEN
!  !       IF (mype_herm == 0) THEN
!  !          DO i=0,nkx0-1
!  !              DO j=0,nky0-1
!  !                  DO k = lkz1,lkz2
!  !                      DO h = lh1,lh2
!  !                      g_in(i,j,k,0,h,:) = dens(i,j,k)*pi**(-5/4.)*e**(-vgrid(h))
!  !                      END DO
!  !                  END DO
!  !              END DO
!  !          END DO
!  !       END IF
!  !      ELSE 
!  !        IF (mype_herm == 0) THEN
!  !          DO i=0,nkx0-1
!  !              DO j=0,nky0-1
!  !                  DO k = lkz1,lkz2
!  !                      g_in(i,j,k,0,:,:) = dens(i,j,k)*pi**(-0.25)
!  !                  END DO
!  !              END DO
!  !          END DO
!  !       END IF
!  !      END IF
!  !  
!  !  END SUBROUTINE cosnoise 
!  
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !!                                                                           !!
!  !!                                check_real                                 !!
!  !!                                                                           !!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE check_real(g_in,sum_tot)
!    USE par_mod
!  
!    COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) 
!    REAL, INTENT(out) :: sum_tot
!    REAL :: gsum
!  
!    INTEGER :: j,k
!  
!    if(np_kz.ne.1) STOP "check_real not yet modified for kz parallelization."
!  
!    sum_tot=0.0
!  
!    !kz=0
!    DO j=1,hky_ind
!      sum_tot=sum_tot+sum( abs(g_in(0,nky0-j,0,:,:,:) - conjg(g_in(0,j,0,:,:,:))) )
!    END DO
!  
!    !ky=0
!    DO k=1,hkz_ind
!      sum_tot=sum_tot+sum( abs( g_in(0,0,nkz0-k,:,:,:) - conjg(g_in(0,0,k,:,:,:))))
!    END DO
!  
!    !the rest
!    DO j=1,nky0-1
!      DO k=1,hkz_ind
!        sum_tot=sum_tot+sum( abs( g_in(0,nky0-j,nkz0-k,:,:,:) - conjg(g_in(0,j,k,:,:,:)) ))
!      END DO
!    END DO
!  
!    !g_in(0,hky_ind+1,:,:)=cmplx(0.0,0.0)
!    !g_in(0,:,hkz_ind+1,:)=cmplx(0.0,0.0)
!    !g_in(0,0,0,:)=cmplx(0.0,0.0)
!  
!    gsum= sum(sum(sum(sum(sum(sum(abs(g_in),1),1),1),1),1),1)
!  
!    IF(gsum.gt.epsilon(1.0)) sum_tot=sum_tot/gsum
!  
!  END SUBROUTINE check_real


