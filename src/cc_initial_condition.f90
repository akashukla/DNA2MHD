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
SUBROUTINE initial_condition
  USE par_mod
  USE mtrandom
  USE mpi
  USE linear_rhs, only: init_force

  use iso_fortran_env, only: int64
  
  IMPLICIT NONE

  INTEGER :: i,j,k, l,ind,ri
  REAL :: kfactor,err
  REAL :: s1, s11,s12,s13,s4
  !REAL :: init_prefactor
  COMPLEX :: phase,phaseb,phasev,phaseby,phasevy
  COMPLEX :: LW1, LC1, RW1, RC1
  REAL :: phase1,phase2,phase1y,phase2y,kspect,myphase1,myphase1y,myphase2,myphase2y,showphase,thb,thv,mythb,mythv,bt1,bt2,bt3,btmag,b1w,b2w,v1w,v2w
  REAL :: b1r,b1i,b2r,b2i,b3r,b3i
  
  REAL :: zerocmplx
  REAL :: kxzeroen=0.0,knzeroen=0.0,kxzeroenm=0.0,knzeroenm=0.0

  INTEGER, DIMENSION(:), ALLOCATABLE :: rseed
  INTEGER :: rseed_size,ierr
  REAL :: testrandoms(10)
  REAL :: truncx,truncy,truncz,nlt
  REAL :: fperpm,fperp
  REAL(C_DOUBLE) :: dt_criticalm,dtc

  REAL :: turnoverm,nltmin,nltmax ! MPI ALL REDUCE variables
  INTEGER(int64) :: t
  REAL(8) :: waveamps1(4) = 0.0,waveamps2(4)

  zerocmplx=0.0
   
  !g_1(:,:,:,:,:,:)=cmplx(0.0,0.0)
  b_1(:,:,:,:)=cmplx(0.0,0.0)
  v_1(:,:,:,:)=cmplx(0.0,0.0)

  xst = max(cstart(1),2)
  yst = max(cstart(2),2)
  zst = max(cstart(3),2)

  if (force_turbulence) CALL init_force

 !init_prefactor=0.001
 !Default Initialization
 !      CALL RANDOM_SEED

 ! This seeding procedure is inspired by the fortran random_seed example
 ! Is probably a bit weaker depending on the clock 
 ! But after 4000 attempts of printing 30 random numbers appeared uniform
 
 CALL RANDOM_SEED(SIZE=rseed_size)
 ALLOCATE(rseed(rseed_size))
 
 CALL system_clock(t)
 t = mod(t,4294967297_int64)

 DO ri = 1,rseed_size
    rseed(ri) = mod(t,4294967297_int64)
    if (random_state.gt.0) rseed(ri) = random_state
 enddo
 
 CALL RANDOM_SEED(PUT=rseed+mype) ! Add mype to ensure different processes have different phases
 DEALLOCATE(rseed)
 
 ! Not truncating initial conditions (for now)
 truncx = 300.0
 truncy = 300.0
 truncz = 300.0
 
 if (enone) s1 = 0.0
 
 DO i=xst,cend(1)
    DO j=yst,cend(2)
       DO k=zst,cend(3)
          ! Write the if statements as subroutines sometime
          
          ! if ((k.ne.nkz0/2).and.(j.ne.nky0/2)) then - this gets purged anyways
          ! Only initalize selected mode numbers
          
          ! Sometime write this loop as a sequence of subroutines
          
!!! Uniform distribution
          
          if ((random_state.ge.0).and.uni) then
             CALL RANDOM_NUMBER(phase1)
             if (phdf.le.1.0) phase2 = phase1 + phdf
             if (phdf.gt.1.0) CALL RANDOM_NUMBER(phase2)
             if (phdfxy.le.1.0) then
                phase1y = phase1+phdfxy
                phase2y = phase2+phdfxy
             else
                CALL RANDOM_NUMBER(phase1y)
                CALL RANDOM_NUMBER(phase2y)
             endif
             CALL RANDOM_NUMBER(thb)
             CALL RANDOM_NUMBER(thv)
             
             
!!! Triangular distribution 
             
          else if (random_state.ge.0) then
             
             CALL RANDOM_NUMBER(myphase1)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_ALLREDUCE(myphase1,phase1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
             if (phdf.le.1.0) phase2 = phase1 + phdf
             if (phdf.gt.1.0) then 
                CALL RANDOM_NUMBER(myphase2) 
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_ALLREDUCE(myphase2,phase2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
             endif
             
             if (phdfxy.le.1.0) then 
                phase1y = phase1 + phdfxy
                phase2y = phase2 + phdfxy
             else
                CALL RANDOM_NUMBER(myphase1y)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_ALLREDUCE(myphase1y,phase1y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
                
                CALL RANDOM_NUMBER(myphase2y)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_ALLREDUCE(myphase2y,phase2y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
             endif
             CALL RANDOM_NUMBER(mythb)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_ALLREDUCE(mythb,thb,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
             
             CALL RANDOM_NUMBER(mythv)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_ALLREDUCE(mythv,thv,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
             
          else ! Make all modes pure real (-1,0) or pure imaginary (0,1)
             if (init_wave.or.beltrami) then
                phase1 = -random_state/4.0
                phase2 = -random_state/4.0
                phase1y = -random_state/4.0
                phase2y = -random_state/4.0
             else
                phase1 = -random_state/4.0
                phase2 = 0.0
                phase1y = 0.0
                phase2y = 0.0
             endif
          endif
          
          ! phase2 = phase1 - 1.0/4.0
          phaseb = cmplx(cos(2*pi*phase1),sin(2*pi*phase1))
          phasev = cmplx(cos(2*pi*phase2),sin(2*pi*phase2))
          phaseby = cmplx(cos(2*pi*phase1y),sin(2*pi*phase1y))
          phasevy = cmplx(cos(2*pi*phase2y),sin(2*pi*phase2y))
               
          CALL RANDOM_NUMBER(showphase)
          showphase = 2.0
          if ((max_itime.lt.100).and.(showphase.lt.10.0/real(nkx0*nky0*nkz0))) then
             print *, i,j,k
             print *, phase1y
             print *, phase2
             print *, phase2y
          endif
          
          IF(kzgrid(k).eq.zerocmplx) THEN
             b_1(i,j,k,0)=cmplx(0.0,0.0)
             b_1(i,j,k,1)=cmplx(0.0,0.0)
             b_1(i,j,k,2)=cmplx(0.0,0.0)
             v_1(i,j,k,0)=cmplx(0.0,0.0)
             v_1(i,j,k,1)=cmplx(0.0,0.0)
             v_1(i,j,k,2)=cmplx(0.0,0.0)
          ENDIF

          ! Rules for initial conditions
          ! All but 0,1 use energy normalization; all but 0 assume init_kolm initial scaling with kperp
          ! 0 - start with zero (force up)
          ! 1 - guess bx,by,vx,vy, start incompressible with init_kolm spectrum; no energy normalization; make divergence free by solving for bz
          ! 2 - use phase and make divergence free by solving for bz
          ! 3 - Taylor Green Vortex, significant for Navier Stokes
          ! <= 10 - use random phase k x z /(sqrt(2) kperp) + i k x (k x z)/(sqrt(2) kperp k) combos
          ! 4: pure helical
          ! 5: minus helical
          ! 6: shear Alfven wave k x z
          ! 7: pseudo Alfven wave k x (k x z)
          ! 8: b = v
          ! 9: b = -v
          ! 10: combination of up/down shear pseudo Alfven waves, require energy fractions
          
          ! >10 - Hall MHD wave decompositon, require energy fraction in each

          IF (init_cond.eq.1) THEN ! Requires guesses of bx,by,vx,vy
             b_1(i,j,k,0)=init_amp_bx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))&
                  *1/(kperps(i,j,k)**(init_kolm/2.0)) * phaseb*cos(2*pi*thb)
             b_1(i,j,k,1)=init_amp_by*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))&
                  *1/(kperps(i,j,k)**(init_kolm/2.0)) * phaseb*phaseby*sin(2*pi*thb)
             b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
             !b_1(i,j,k,2)=init_amp_bz
             v_1(i,j,k,0)=init_amp_vx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))&
                  *1/(kperps(i,j,k)**(init_kolm/2.0)) * phasev*phaseb*cos(2*pi*thv)
             v_1(i,j,k,1)=init_amp_vy*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))&
                  *1/(kperps(i,j,k)**(init_kolm/2.0)) * phasev*phasevy*phaseb*sin(2*pi*thv)
             !v_1(i,j,k,2)=init_amp_vz
             v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
          ELSE IF (init_cond.eq.2) THEN ! Remove divergence by solving divergence free condition
             b_1(i,j,k,0)= phaseb*cos(2*pi*thb)*1/(kperps(i,j,k)**(init_kolm/2.0))
             b_1(i,j,k,1)= phaseb*phaseby*sin(2*pi*thb)*1/(kperps(i,j,k)**(init_kolm/2.0))
             b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
             v_1(i,j,k,0)= phasev*phaseb*cos(2*pi*thv)*1/(kperps(i,j,k)**(init_kolm/2.0))
             v_1(i,j,k,1)= phasev*phaseb*phaseby*sin(2*pi*thv)*1/(kperps(i,j,k)**(init_kolm/2.0))
             v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
          ELSE IF (init_cond.le.10) THEN
             ! These initial conditions write a Beltrami decomposition for b,v in terms of curl eigenstates (b1r,\pm b1i) as given below
             
             CALL RANDOM_NUMBER(b1w)
             CALL RANDOM_NUMBER(b2w)
             CALL RANDOM_NUMBER(v1w)
             CALL RANDOM_NUMBER(v2w)
             
             CALL MPI_BCAST(b1w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(b2w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(v1w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(v2w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
             
             IF (init_cond.eq.4) THEN
                b2w = 0.0
                v2w = 0.0
                phasev = phaseb
                v1w = b1w
             ELSE IF (init_cond.eq.5) THEN
                b1w = 0.0
                v1w = 0.0
                phasevy = phaseby
                v2w = b2w                
             ELSE IF (init_cond.eq.6) THEN ! k x z Alfven shear wave
                phaseby = phaseb
                phasevy = phasev
                b1w = b2w
                v1w = v2w
             ELSE IF (init_cond.eq.7) THEN ! k x (k x z) pseudo Alfven wave
                phaseby = -phaseb
                phasevy = -phasev
                b1w = b2w
                v1w = v2w
             ELSE IF ((init_cond.eq.8).or.(init_cond.eq.9)) THEN
                if (init_cond.eq.8) then
                   phasev = phaseb
                   phasevy = phaseby
                else 
                   phasev = -phaseb
                   phasevy = -phaseby
                endif
                v1w = b1w
                v2w = b2w
             ENDIF

             b_1(i,j,k,0) = (b1w*phaseb*pcurleig(i,j,k,0) &
                  + b2w*phaseby*conjg(pcurleig(i,j,k,0)))/sqrt(b1w**2+b2w**2)
             b_1(i,j,k,1) = (b1w*phaseb*pcurleig(i,j,k,1) &
                  + b2w*phaseby*conjg(pcurleig(i,j,k,1)))/sqrt(b1w**2+b2w**2)
             b_1(i,j,k,2) = (b1w*phaseb*pcurleig(i,j,k,2) &
                  + b2w*phaseby*conjg(pcurleig(i,j,k,2)))/sqrt(b1w**2+b2w**2)
             
             v_1(i,j,k,0) = (v1w*phasev*pcurleig(i,j,k,0) &
                  + v2w*phasevy*conjg(pcurleig(i,j,k,0)))/sqrt(v1w**2+v2w**2)
             v_1(i,j,k,1) = (v1w*phasev*pcurleig(i,j,k,1) &
                  + v2w*phasevy*conjg(pcurleig(i,j,k,1)))/sqrt(v1w**2+v2w**2)
             v_1(i,j,k,2) = (v1w*phasev*pcurleig(i,j,k,2) &
                  + v2w*phasevy*conjg(pcurleig(i,j,k,2)))/sqrt(v1w**2+v2w**2)
             
             b_1(i,j,k,:) = b_1(i,j,k,:) / (kperps(i,j,k)**(init_kolm/2.0))
             v_1(i,j,k,:) = v_1(i,j,k,:) / (kperps(i,j,k)**(init_kolm/2.0))

             
             IF (init_cond.eq.10) THEN ! Mixture of shear, pseudoAlfven waves
             
                b_1(i,j,k,:) = phaseb * sqrt(en_leftwhist) * (pcurleig(i,j,k,:)+conjg(pcurleig(i,j,k,:))) &
                     + phasev * sqrt(en_leftcyclo) * (pcurleig(i,j,k,:)+conjg(pcurleig(i,j,k,:))) &
                     + phaseby * sqrt(en_rightwhist) * (pcurleig(i,j,k,:)-conjg(pcurleig(i,j,k,:))) &
                     + phasevy * sqrt(en_rightcyclo) * (pcurleig(i,j,k,:)-conjg(pcurleig(i,j,k,:)))

                v_1(i,j,k,:) = phaseb *	sqrt(en_leftwhist) * (pcurleig(i,j,k,:)+conjg(pcurleig(i,j,k,:))) &
                     - phasev *	sqrt(en_leftcyclo) * (pcurleig(i,j,k,:)+conjg(pcurleig(i,j,k,:))) &
                     + phaseby * sqrt(en_rightwhist) * (pcurleig(i,j,k,:)-conjg(pcurleig(i,j,k,:))) &
                     - phasevy * sqrt(en_rightcyclo) * (pcurleig(i,j,k,:)-conjg(pcurleig(i,j,k,:)))

                b_1(i,j,k,:) = b_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))
                v_1(i,j,k,:) = v_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))

             ENDIF
             
          ELSE ! Mixture of Hall MHD waves
             
             LW1 = sqrt(en_leftwhist) * phaseb/sqrt(1 + alpha_leftwhist(i,j,k)**2)
             LC1 = sqrt(en_leftcyclo) * phasev/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
             RW1 = sqrt(en_rightwhist) * phaseby/sqrt(1 + alpha_leftwhist(i,j,k)**2)
             RC1 = sqrt(en_rightcyclo) * phasevy/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
             
             b_1(i,j,k,:) = (LW1 * alpha_leftwhist(i,j,k) + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                  - (RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k))*conjg(pcurleig(i,j,k,:))
             v_1(i,j,k,:) = (LW1 + LC1) * pcurleig(i,j,k,:) + (RW1 + RC1)*conjg(pcurleig(i,j,k,:))
             
             b_1(i,j,k,:) = b_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))
             v_1(i,j,k,:) = v_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))
             
          ENDIF
          
       ENDDO
    ENDDO
 ENDDO
 ! Filter modes in the padding region
 DO ind = 0,2
    b_1(:,:,:,ind) = b_1(:,:,:,ind) * paddingmask
    v_1(:,:,:,ind) = v_1(:,:,:,ind) * paddingmask
 ENDDO
 if (force_trunc) then ! Only initialize forced waves
    DO ind = 0,2 
       b_1(:,:,:,ind) = b_1(:,:,:,ind) * mask1
       v_1(:,:,:,ind) = v_1(:,:,:,ind) * mask1
    ENDDO
 ENDIF

 IF (init_cond.ge.31) THEN ! Three wave systems - excite for resonance condition

    b_1 = cmplx(0.0,0.0)
    v_1 = cmplx(0.0,0.0)

    ! conditions 31-35 derived for perp scale 0.05 and par scale 0.01

    ! +Whistler +Whistler Interactions

    ! Parallel 
    IF (init_cond.eq.31) THEN

       ! (0,35,9) (0,3,-1)
       ! (0,10,15) (0,18,-27)
       ! (0,15,9) (0,8,-5)

       CALL setwaveindices(0,15,9,0,8,-5)

       waveamps1(1) = 1.0
       waveamps2(1) = 1.0

    ENDIF

    ! Antiparallel
    IF (init_cond.eq.32) THEN ! Low k Whistler + Whistler Antipar

       ! (0,35,9) (0,-20,-35)
       ! (0,10,15) (0,-28,12)
       ! (0,15,9) (0,-8,-33)

       CALL setwaveindices(0,15,9,0,-8,-33)

       waveamps1(1) = 1.0
       waveamps2(1) = 1.0

    ENDIF

    ! Perpendicular
    IF (init_cond.eq.33) THEN ! Low k +Whistler +Cyclotron Near Perp

       ! (0,35,9) (27,0,-5)
       ! (0,10,15) (8,0,-8)
       ! (0,15,9) (18,0,-14)

       CALL setwaveindices(0,15,9,18,0,-14)

       waveamps1(1) = 1.0
       waveamps2(1) = 1.0
       
    ENDIF

    IF (init_cond.eq.34) THEN

       ! (0,15,9) (28,27,-26)

       CALL setwaveindices(0,15,9,28,27,-26)

       waveamps1(1) = 1.0
       waveamps2(1) = 1.0

    ENDIF

    IF (init_cond.eq.35) THEN

       ! (0,10,15) (12,2,-61)

       CALL setwaveindices(0,10,15,12,2,-61)

       waveamps1(1) = 1.0
       waveamps2(1) = 1.0

    ENDIF

    ! Optimized 0.025 kperpmin 0.005 kzmin modes from contour plot
    IF (init_cond.eq.200) CALL setwaveindices(wave1x,wave1y,wave1z,wave2x,wave2y,wave2z)
    
    waveamps1(1) = 1.0
    waveamps2(1) = 1.0
    CALL isolatedhmhdwave(wave1x,wave1y,wave1z,waveamps1*2.0,mype1)
    CALL isolatedhmhdwave(wave2x,wave2y,wave2z,waveamps1*1.0,mype2)
    CALL isolatedhmhdwave(wave3x,wave3y,wave3z,waveamps2*0.5,mype3)

 ENDIF
 
 
 if (verbose.and.(mype.eq.0)) print *, "Through initial b_1 and v_1",mype
 
 ! Set energy as fraction of 4 pi^3
 if (init_cond.ge.2) then
    
    knzeroenm = sum(abs(b_1(xst:cend(1),:,:,:))**2+abs(v_1(xst:cend(1),:,:,:))**2)
    ! Account for zeros
    if (cstart(1).eq.1) kxzeroenm = sum(0.5*(abs(b_1(1,:,:,:))**2+abs(v_1(1,:,:,:))**2))
    
    if (verbose.and.(mype.eq.0)) print *, mype,"Unnormalized Sum",knzeroenm+kxzeroenm
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(knzeroenm,knzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(kxzeroenm,kxzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    s1 = knzeroen+kxzeroen
    if (verbose.and.(mype.eq.0)) print *, mype,"s1",s1
    
    b_1 = b_1 * sqrt(init_energy / (2.0*s1))
    v_1 = v_1 * sqrt(init_energy / (2.0*s1))
 endif
 
 if (verbose.and.(mype.eq.0)) print *, "Through energy normalization"

  if (init_cond.eq.3) then
     ! Initializes Fourier components of Taylor Green vortex u = sin x cos y cos z , v = - cos x sin y cos z
     v_1 = cmplx(0.0,0.0)
     
     if (((cstart(1).le.2).and.(cend(1).ge.2)).and.((cstart(2).le.2).and.(cend(2).ge.2))) then
        v_1(2,2,0,0) = 1.0
        v_1(2,2,0,1) = -1.0
     else if (((cstart(1).le.2).and.(cend(1).ge.2)).and.(cend(2).eq.nky0)) then
        v_1(2,nky0,0,0) = 1.0
        v_1(2,nky0,0,1) = 1.0
     endif
     
     v_1 = -v_1 * cmplx(0.0,0.125)
  endif
 
  mhelcorr = 0.0
  if (checkpoint_read) call checkpoint_in
 
  turnoverm = 0.0
  nltmax = 10.0**8.0
  nltmin = 0.0
  
  DO i = cstart(1),cend(1)
     DO j = cstart(2),cend(2)
        DO k = cstart(3),cend(3)
           
           turnoverm = turnoverm + sum(abs(v_1(i,j,k,:))**2.0) * sin(kxgrid(i)/(2.0*kxmin))**2.0
           nltmax = min(nltmax,minval(kmags(i,j,k)*abs(v_1(i,j,k,:)),kmags(i,j,k)*abs(v_1(i,j,k,:)).gt.10.0**(-10.0)))
           nltmin = max(nltmin,kmags(i,j,k)*maxval(abs(v_1(i,j,k,:)),abs(v_1(i,j,k,:)).gt.10.0**(-10.0)))
           
        ENDDO
     ENDDO
  ENDDO
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(turnoverm,turnover,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  turnover = 1.0/(2*kxmin*sqrt(turnover))
  
  if (mype.eq.0) print *, mype,"Turnover Time Estimate",turnover
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(nltmax,nlt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
  
  if (mype.eq.0) nlt = 10.0/nlt
  
  if (mype.eq.0) print *, mype,"Equation-Based Max Nonlinear Time Scale", nlt
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(nltmin,nlt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if (mype.eq.0) nlt = 10.0/nlt
  
  if (mype.eq.0) print *, mype,"Equation-Based Min Nonlinear Time Scale", nlt
  
  ! if (nv) b_1(:,:,:,:) = cmplx(0.0,0.0)
  
  ! Scale forcing to represent energy gained per time step |F| |v| ~ energy gain per step as fraction of initial energy
  ! |v| ~ sqrt(|v|^2) ~ sqrt(energy/2)
 
  ! if (force_turbulence) force_amp = force_amp * sqrt(4.0 * pi**3.0 * init_energy)
  ! if (verbose) print *, mype,"Force amp",force_amp
  
  IF (init_cond.eq.0) THEN
     b_1 = cmplx(0.0,0.0)
     v_1 = cmplx(0.0,0.0)
  ENDIF
  
  ! Linear stability maximum time step
  dt_criticalm = 2.0/(maxval(kzgrid)*(maxval(hall*kmags)/2 + sqrt(1 + 0.25*maxval(hall*kmags)**2.0)))
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(dt_criticalm,dtc,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
  if (mype.eq.0) print *, "Gauss2 Critical Time Step", dtc
  if (calc_dt.and.(.not.(test_ho)).and.(hall.ne.0.0)) dt_max = minval([dt_max,dtc/2.0])
  if (verbose.and.(mype.eq.0)) then
     print *, "kzgrid max", maxval(kzgrid)
     print *, "kmags max", kmax
 endif
 
 
 ! Check on Initial Energy

 if (mype.eq.0) b_1(cstart(1),cstart(2),cstart(3),2) = 1.0
 
 knzeroenm = sum(abs(b_1(xst:cend(1),:,:,:))**2+abs(v_1(xst:cend(1),:,:,:))**2)
 if (cstart(1).eq.1) kxzeroenm = sum(0.5*(abs(b_1(1,:,:,:))**2+abs(v_1(1,:,:,:))**2))
 
 print *, mype, "Mype Energy",8*pi**3 * (knzeroenm+kxzeroenm)
 
 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
 CALL MPI_ALLREDUCE(knzeroenm,knzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
 
 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
 CALL MPI_ALLREDUCE(kxzeroenm,kxzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
 
 s1 = knzeroen + kxzeroen
 
 if (mype.eq.0) print *, "All Mype Initial Energy",s1*8*pi**3

 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 if (rey.eq.0) then
    rey = kxmin/vnu * sqrt(2.0*init_energy) * (nkx0)**(2.0*hyp)
    vnu = vnu / (kmax**(2.0*hyp))
 else
    vnu = kxmin/(rey * kxmin**(2*hyp)) * sqrt(force_amp * 8*pi **3 )
 endif

 eta = eta * vnu
 if (mype.eq.0) print *, "Perp Reynolds Number",rey
 if (mype.eq.0) print *, 'Viscosity',vnu
 if (mype.eq.0) print *, "Force Amp",force_amp
 
 dt = dt_max

 
END SUBROUTINE initial_condition

SUBROUTINE isolatedhmhdwave(ix,iy,iz,waveamps,mypewave)

  use par_mod
  use mpi
  
  IMPLICIT NONE

  integer(4) :: ix,iy,iz ! Zero based indices for array
  real(8) :: waveamps(4)
  integer :: mypewave
  integer(4) :: wavex,wavey,wavez
  real(8) :: normamps(4)
  integer :: mypewavem
  integer(4) :: minusiy,minusiz
  integer(4) :: ierr

  mypewavem = 0

  if (ix.ge.cstart(1).and.ix.le.cend(1)) then

     if (iy.ge.cstart(2).and.iy.le.cend(2)) then

        mypewavem = mype
        
        if (verbose) print *, "Isolated Wave Mype",mype,ix,iy,iz

        normamps(1) = waveamps(1)/sqrt(1 + alpha_leftwhist(ix,iy,iz)**2)
        normamps(2) = waveamps(2)/sqrt(1 + alpha_leftcyclo(ix,iy,iz)**2)
        normamps(3) = waveamps(3)/sqrt(1 + alpha_leftwhist(ix,iy,iz)**2)
        normamps(4) = waveamps(4)/sqrt(1 + alpha_leftcyclo(ix,iy,iz)**2)

        b_1(ix,iy,iz,:) = (normamps(1) * alpha_leftwhist(ix,iy,iz) + normamps(2) * alpha_leftcyclo(ix,iy,iz)) * pcurleig(ix,iy,iz,:)&
             - (normamps(3)* alpha_leftwhist(ix,iy,iz) + normamps(4) * alpha_leftcyclo(ix,iy,iz))*conjg(pcurleig(ix,iy,iz,:))
        v_1(ix,iy,iz,:) = (normamps(1) + normamps(2)) * pcurleig(ix,iy,iz,:) + (normamps(3) + normamps(4))*conjg(pcurleig(ix,iy,iz,:))

     endif
     
     if (ix.eq.1) then ! preserve reality

        minusiy = ny0_big+2-iy
        minusiz = nz0_big+2-iz

        if (minusiy.ge.cstart(2).and.minusiy.le.cend(2)) then

           normamps(1) = waveamps(1)/sqrt(1 + alpha_leftwhist(ix,minusiy,minusiz)**2)
           normamps(2) = waveamps(2)/sqrt(1 + alpha_leftcyclo(ix,minusiy,minusiz)**2)
           normamps(3) = waveamps(3)/sqrt(1 + alpha_leftwhist(ix,minusiy,minusiz)**2)
           normamps(4) = waveamps(4)/sqrt(1 + alpha_leftcyclo(ix,minusiy,minusiz)**2)

           b_1(ix,minusiy,minusiz,:) = -(normamps(1) * alpha_leftwhist(ix,minusiy,minusiz) + normamps(2) * alpha_leftcyclo(ix,minusiy,minusiz)) * pcurleig(ix,minusiy,minusiz,:)&
                +(normamps(3)* alpha_leftwhist(ix,minusiy,minusiz) + normamps(4) * alpha_leftcyclo(ix,minusiy,minusiz))*conjg(pcurleig(ix,minusiy,minusiz,:))
           v_1(ix,minusiy,minusiz,:) = -(normamps(1) + normamps(2)) * pcurleig(ix,minusiy,minusiz,:) - (normamps(3) + normamps(4))*conjg(pcurleig(ix,minusiy,minusiz,:))
           
        endif
     
     endif

  endif
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mypewavem,mypewave,1,MPI_INTEGER4,MPI_MAX,MPI_COMM_WORLD,ierr)
  
END SUBROUTINE isolatedhmhdwave

SUBROUTINE setwaveindices(input1x,input1y,input1z,input2x,input2y,input2z)

  ! Set wave indices from inputs
  ! Make sure x is always positive

  use par_mod

  implicit none
  integer(4), intent(in) :: input1x,input1y,input1z,input2x,input2y,input2z
  
  wave1x =	input1x
  wave1y =	input1y
  wave1z =	input1z
  
  wave2x =	input2x
  wave2y =	input2y
  wave2z =	input2z
  
  wave3x = 1 + (wave1x + wave2x)
  wave1x = wave1x + 1
  wave2x = wave2x + 1
  
  CALL adjustwaveindices(wave1y,wave2y,wave3y,1)
  CALL adjustwaveindices(wave1z,wave2z,wave3z,2)  
  
END SUBROUTINE setwaveindices

SUBROUTINE adjustwaveindices(index1,index2,index3,ind)

  ! Helper routine to adjust the wave index from a zero-centered grid to the one-base

  use par_mod

  implicit none
  integer(4), intent(inout) :: index1,index2,index3
  integer(4), intent(in) :: ind
  integer(4) :: dimsize

  if (ind.eq.1) dimsize = ny0_big
  if (ind.eq.2) dimsize = nz0_big
  
  index3 = index1 + index2
  if (index3.ge.0) index3 = index3 + 1
  if (index3.lt.0) index3 = dimsize + (index3+1)

  if (index2.ge.0) index2 = index2 + 1
  if (index2.lt.0) index2 = dimsize + (index2+1)

  if (index1.ge.0) index1 = index1 + 1
  if (index1.lt.0) index1 = dimsize + (index1+1)

END SUBROUTINE adjustwaveindices
