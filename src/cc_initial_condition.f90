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
  USE mpi
  USE diagnostics, only : bv_first

  use iso_fortran_env, only: int64
  
  IMPLICIT NONE

  INTEGER :: i,j,k, l,zst,xst,yst,ri
  REAL :: kfactor,err
  REAL :: s1, s11,s12,s13,s4
  !REAL :: init_prefactor
  COMPLEX :: phase,phaseb,phasev,phaseby,phasevy
  COMPLEX :: LW1, LC1, RW1, RC1
  REAL :: phase1,phase2,phase1y,phase2y,kspect,myphase1,myphase1y,myphase2,myphase2y,showphase,thb,thv,mythb,mythv,bt1,bt2,bt3,b1r,b2r,b3r,b1i,b2i,b3i,btmag,b1w,b2w,v1w,v2w
  CHARACTER(len=40), INTENT(in) :: which_init0
  CHARACTER(len=40) :: which_init
  REAL :: zerocmplx
  REAL :: kxzeroen,knzeroen,kxzeroenm,knzeroenm

  INTEGER, DIMENSION(:), ALLOCATABLE :: rseed
  INTEGER :: rseed_size,ierr
  REAL :: testrandoms(10)
  REAL :: truncx,truncy,truncz,nlt
  REAL :: fperpm,fperp

  REAL :: turnoverm,nltmin,nltmax ! MPI ALL REDUCE variables
  INTEGER(int64) :: t 

  zerocmplx=0.0
   
  !g_1(:,:,:,:,:,:)=cmplx(0.0,0.0)
  b_1(:,:,:,:)=cmplx(0.0,0.0)
  v_1(:,:,:,:)=cmplx(0.0,0.0)

  xst = max(1,kxinit_min)
  yst = max(1,kyinit_min)
  zst = max(1,lkz1)

  DO i = 0,nkx0-1
    DO j = 0,nky0-1
      DO k = lkz1,lkz2
        kmags(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2)
        kperps(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2)
        kzs(i,j,k) = kzgrid(k)
      END DO
    END DO
 END DO

 if (verbose) print *, mype,"MaxVal ks",maxval(kmags),maxval(kperps)

 ALLOCATE(alpha_leftwhist(0:nkx0-1,0:nky0-1,lkz1:lkz2))
 ALLOCATE(alpha_leftcyclo(0:nkx0-1,0:nky0-1,lkz1:lkz2))

 alpha_leftwhist = -kmags/2.0 - sqrt(1.0 + ((kmags**2.0) / 4.0))
 alpha_leftcyclo = -kmags/2.0 + sqrt(1.0 + ((kmags**2.0) / 4.0))

 if (verbose) print *, mype,"MaxVal alphas",maxval(alpha_leftwhist),maxval(alpha_leftcyclo)
 
 DO i = 0,nkx0-1
    DO j = 0,nky0-1
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

 if (verbose) print *, mype,"MaxVal PC",maxval(abs(pcurleig))

 if (verbose) print *, mype,"b1 Extent",sum(0*b_1 + 1)

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

      if (rey.eq.0) then
         ! Set viscosity to set relative rate of dissipation at high scales 
         rey = kxmin/vnu * sqrt(force_amp * 8*pi **3 ) * (nkx0)**(2.0*hyp)
         vnu = vnu / (kmax**(2.0*hyp))
      else
         ! Set viscosity from Reynolds number
         vnu = kxmin/(rey * kxmin**(2*hyp)) * sqrt(force_amp * 8*pi **3 )
      endif
      
      if ((mype.eq.0)) print *, "Perp Reynolds Number",rey
      
      ! Set resistivity from Magnetic Prandtl number
      eta = eta * vnu
      if(mype.eq.0) print *, 'Viscosity',vnu

      ! Set forcing such that a smallest wavenumber mode
      ! would have energy force_amp times guide field energy divided by number of modes forced
      ! Order of magnitude estimate pending phases and forcing, of course

      fperpm = sum(abs( kperps .lt. (maxval(kperps) * force_frac)))

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fperpm,fperp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      force_amp = abs(cmplx(vnu*(kxmin)**(2.0*hyp),0)) * sqrt(1.0/3.0 * force_amp) * sqrt(8.0 * (pi ** 3))/fperp
      
      print *, "Force Amp",force_amp      

      IF (init_null) THEN
         b_1 = cmplx(0.0,0.0)

         ! Set forcing amplitude based on estimate of steady state v_1 = f/(vnu k^2h)
         ! Then (f^2/vnu^2 sum(abs(v1^2/kmags^(4*hyp)))) = force_amp * (4 pi ** 3)

         ! force_amp = vnu * sqrt(init_energy * (4 * pi**3) / sum((abs(v_1)**2) /spread(kmags**(4*hyp),4,3),spread(kmags**(4*hyp),4,3).gt.10**(-10)))
         ! if ((mype.eq.0)) print *, "Force Amp",force_amp

         ! v_1 = v_1 * force_amp/(vnu * spread(kmags**(2*hyp),4,3))
         ! v_1(0,0,0,:) = 0

         ! Guessed v_1 amps

         v_1 = spread(kperps ** (-init_kolm * 0.5),4,3)
         v_1(0,:,:,:) = zerocmplx
         v_1(:,0,:,:) = zerocmplx
         v_1(:,:,0,:) = zerocmplx

         knzeroenm = sum(abs(b_1(1:nkx0-1,:,:,:))**2+abs(v_1(1:nkx0-1,:,:,:))**2)
         kxzeroenm = sum(0.5*(abs(b_1(0,:,:,:))**2+abs(v_1(0,:,:,:))**2))

         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_ALLREDUCE(knzeroenm,knzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	 CALL MPI_ALLREDUCE(kxzeroenm,kxzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

         s1 = knzeroen + kxzeroen

         v_1 = v_1 * sqrt(0.5 * init_energy / (s1 * 2.0))

         turnover = 1/(2*kxmin * sqrt(sum(sum(abs(v_1)**2.0,4)*sin(spread(spread(kxgrid/(2*kxmin),2,nky0),3,lkz2+1-lkz1))**2)))
         print *, mype,"Turnover Time Estimate",turnover

         nlt = 10.0/(minval(abs(spread(kmags,4,3)*v_1),abs(spread(kmags,4,3)*v_1).gt. 10**(-10)))

         print *, mype,"Equation-Based Nonlinear Time Scale", nlt

         ! Revert back to null IC
         v_1 = cmplx(0.0,0.0)

      ELSE         

      DO i=xst,nkx0-1
         DO j=yst,nky0-1
            DO k=zst,lkz2
              ! Write the if statements as subroutines sometime

               if ((k.ne.nkz0/2).and.(j.ne.nky0/2)) then
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
                 ELSE IF (.not.(enone)) THEN
                    b_1(i,j,k,0)=init_amp_bx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kperps(i,j,k)**(init_kolm/2.0)) * phaseb*cos(2*pi*thb)
                    b_1(i,j,k,1)=init_amp_by*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kperps(i,j,k)**(init_kolm/2.0)) * phaseb*phaseby*sin(2*pi*thb)
                    b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
                    !b_1(i,j,k,2)=init_amp_bz
                    v_1(i,j,k,0)=init_amp_vx*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kperps(i,j,k)**(init_kolm/2.0)) * phasev*phaseb*cos(2*pi*thv)
                    v_1(i,j,k,1)=init_amp_vy*1.0/sqrt(real(nkx0*nky0*(nkz0-1)))*1/(kperps(i,j,k)**(init_kolm/2.0)) * phasev*phasevy*phaseb*sin(2*pi*thv)
                    !v_1(i,j,k,2)=init_amp_vz
                    v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
                 ELSE IF (beltrami) THEN
                    ! These initial conditions write a Beltrami decomposition for b,v in terms of curl eigenstates (b1r,\pm b1i) as given below
                    
                    CALL RANDOM_NUMBER(b1w)
                    CALL RANDOM_NUMBER(b2w)
                    CALL RANDOM_NUMBER(v1w)
                    CALL RANDOM_NUMBER(v2w)
                    
                    CALL MPI_BCAST(b1w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
                    CALL MPI_BCAST(b2w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
                    CALL MPI_BCAST(v1w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
                    CALL MPI_BCAST(v2w,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
                    
                    if (bc_norm) shear=.true.
                    
                    IF (helical) THEN
                       b2w = 0.0
                       v2w = 0.0
                       phasev = phaseb
                       v1w = b1w
                    ELSE IF (shear) THEN
                       phaseby = phaseb
                       phasevy = phasev
                       b1w = b2w
                       v1w = v2w
                    ELSE IF ((walenp).or.(walenn)) THEN
                       if (walenp) then
                          phasev = phaseb
                          phasevy = phaseby
                       else 
                          phasev = -phaseb
                          phasevy = -phaseby
                       endif
                       v1w = b1w
                       v2w = b2w
                       
                    ENDIF
                    
                    b_1(i,j,k,0) = (b1w*phaseb*pcurleig(i,j,k,0) + b2w*phaseby*conjg(pcurleig(i,j,k,0)))/sqrt(b1w**2+b2w**2)
                    b_1(i,j,k,1) = (b1w*phaseb*pcurleig(i,j,k,1) + b2w*phaseby*conjg(pcurleig(i,j,k,1)))/sqrt(b1w**2+b2w**2)
                    b_1(i,j,k,2) = (b1w*phaseb*pcurleig(i,j,k,2) + b2w*phaseby*conjg(pcurleig(i,j,k,2)))/sqrt(b1w**2+b2w**2)
                    
                    v_1(i,j,k,0) = (v1w*phasev*pcurleig(i,j,k,0) + v2w*phasevy*conjg(pcurleig(i,j,k,0)))/sqrt(v1w**2+v2w**2)
                    v_1(i,j,k,1) = (v1w*phasev*pcurleig(i,j,k,1) + v2w*phasevy*conjg(pcurleig(i,j,k,1)))/sqrt(v1w**2+v2w**2)
                    v_1(i,j,k,2) = (v1w*phasev*pcurleig(i,j,k,2) + v2w*phasevy*conjg(pcurleig(i,j,k,2)))/sqrt(v1w**2+v2w**2)
                    
                    b_1(i,j,k,:) = b_1(i,j,k,:) / (kperps(i,j,k)**(init_kolm/2.0))
                    v_1(i,j,k,:) = v_1(i,j,k,:) / (kperps(i,j,k)**(init_kolm/2.0))

                 ELSE IF (init_wave) THEN
                    
                    LW1 = sqrt(en_leftwhist) * phaseb/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                    LC1 = sqrt(en_leftcyclo) * phasev/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                    RW1 = sqrt(en_rightwhist) * phaseby/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                    RC1 = sqrt(en_rightcyclo) * phasevy/sqrt(1 + alpha_leftcyclo(i,j,k)**2)

                    ! Given a normal boundary condition choice v = 1/kp k  x z
                    ! This condition sets left-right branch modes equal amplitude and maintains above condition
                    if (bc_norm) then
                       RW1 = LW1
                       RC1 = LC1
                    endif

                    b_1(i,j,k,:) = (LW1 * alpha_leftwhist(i,j,k) + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                         - (RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k))*conjg(pcurleig(i,j,k,:))
                    
                    v_1(i,j,k,:) = (LW1 + LC1) * pcurleig(i,j,k,:) + (RW1 + RC1)*conjg(pcurleig(i,j,k,:))
                    
                    b_1(i,j,k,:) = b_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))
                    v_1(i,j,k,:) = v_1(i,j,k,:) / (sqrt(en_leftwhist + en_leftcyclo + en_rightwhist + en_rightcyclo) * kperps(i,j,k)**(init_kolm/2.0))
                    
                 ELSE   

                    IF ((walenp).or.(walenn)) THEN
                       if (walenp) phasev = 1.0
                       if (walenn) phasev = -1.0
                       phasevy =phaseby
                       thv = thb
                    ENDIF
                    
                   IF (.not.(shear.or.helical)) THEN
                      s11 = kxgrid(i)**2 * cos(2*pi*thb)**2 + kxgrid(i)**2 * cos(2*pi*thv)**2
                      s12 = kygrid(j)**2 * sin(2*pi*thb)**2 + kygrid(j)**2 * sin(2*pi*thv)**2
                      s13 = kxgrid(i) * kygrid(j) * (sin(4*pi*thb) * cos(2*pi*phase1y) + sin(4*pi*thv)*cos(2*pi*phase2y))
                      s1 = s1 + kperps(i,j,k)**(-1.0*init_kolm) * (2+(s11+s12+s13)/(kzgrid(k)**2))
                   ELSE
                      s1 = s1 + 2.0*kperps(i,j,k) ** (-1.0*init_kolm)
                   ENDIF
                   
                   IF (helical) THEN
                      bt1 = 2.0*phase1y - 1.0
                      bt2 = 2.0*phase2 - 1.0
                      bt3 = 2.0*phase2y - 1.0
                      b1r = bt1 - (kxgrid(i) * bt1 + kygrid(j) * bt2 + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kxgrid(i)
                      b2r = bt2 - (kxgrid(i) * bt1 + kygrid(j) * bt2 + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kygrid(j)
                      b3r = bt3 - (kxgrid(i) * bt1 + kygrid(j) * bt2 + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kzgrid(k)
                      btmag = sqrt(b1r**2 + b2r**2 + b3r**2)
                      b1i = (kygrid(j) * b3r - kzgrid(k) * b2r)/kmags(i,j,k)
                      b2i = (kzgrid(k) * b1r - kxgrid(i) * b3r)/kmags(i,j,k)
                      b3i = (kxgrid(i) * b2r - kygrid(j) * b1r)/kmags(i,j,k)
                      b_1(i,j,k,0) = cmplx(b1r,b1i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                      b_1(i,j,k,1) = cmplx(b2r,b2i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                      b_1(i,j,k,2) = cmplx(b3r,b3i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                      v_1(i,j,k,:) = b_1(i,j,k,:)
                   ELSE IF (shear.and.(max(i,j).ne.0)) THEN
                      b_1(i,j,k,0) = - kygrid(j)/sqrt(kxgrid(i)**2 + kygrid(j)**2) * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                      b_1(i,j,k,1) = kxgrid(i)/sqrt(kxgrid(i)**2 + kygrid(j)**2) * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                      b_1(i,j,k,2) = cmplx(0.0,0.0)
                      v_1(i,j,k,0) = -kygrid(j)/sqrt(kxgrid(i)**2 + kygrid(j)**2) * phasev * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                      v_1(i,j,k,1) = kxgrid(i)/sqrt(kxgrid(i)**2 + kygrid(j)**2) * phasev*phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                      v_1(i,j,k,2) = cmplx(0.0,0.0)
                   ELSE
                      b_1(i,j,k,0)= phaseb*cos(2*pi*thb)*1/(kperps(i,j,k)**(init_kolm/2.0))
                      !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
                      b_1(i,j,k,1)= phaseb*phaseby*sin(2*pi*thb)*1/(kperps(i,j,k)**(init_kolm/2.0))
                      !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
                      b_1(i,j,k,2) = (-kxgrid(i)*b_1(i,j,k,0)-kygrid(j)*b_1(i,j,k,1))/kzgrid(k)
                      !b_1(i,j,k,2)=init_amp_bz              
                      v_1(i,j,k,0)= phasev*phaseb*cos(2*pi*thv)*1/(kperps(i,j,k)**(init_kolm/2.0))
                      !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
                      v_1(i,j,k,1)= phasev*phasevy*phaseb*sin(2*pi*thv)*1/(kperps(i,j,k)**(init_kolm/2.0))
                      !0.32/sqrt((2+pi**2 * real(nkx0+nky0)/144.0) * real(nkx0*nky0*(nkz0-1)) * 8 * pi**3)
                      !v_1(i,j,k,2)=init_amp_vz 
                      v_1(i,j,k,2) = (-kxgrid(i)*v_1(i,j,k,0)-kygrid(j)*v_1(i,j,k,1))/kzgrid(k)
                   END IF
                   
                END IF

             endif

             if (force_trunc) then !Only initalize forced waves
                if (kperps(i,j,k) > (maxval(kperps)*force_frac)) then
                   b_1(i,j,k,:) = cmplx(0.0,0.0)
                   v_1(i,j,k,:) = cmplx(0.0,0.0)
                end if
             endif
             
          ENDDO
       ENDDO
    ENDDO

    if (verbose) print *, "Through initial b_1 and v_1",mype

    ! Set energy as fraction of 4 pi^3
    if (enone) then

          knzeroenm = sum(abs(b_1(1:nkx0-1,:,:,:))**2+abs(v_1(1:nkx0-1,:,:,:))**2)
          kxzeroenm = sum(0.5*(abs(b_1(0,:,:,:))**2+abs(v_1(0,:,:,:))**2))

          if (verbose) print *, mype,"Unnormalized Sum",knzeroenm+kxzeroenm
          
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(knzeroenm,knzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(kxzeroenm,kxzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

          s1 = knzeroen+kxzeroen
          if (verbose) print *, mype,"s1",s1

    
       b_1 = b_1 * sqrt(init_energy / (2.0*s1))
       v_1 = v_1 * sqrt(init_energy / (2.0*s1))
    endif

    if (verbose) print *, "Through energy normalization"

    turnoverm = 0.0
    nltmax = 10.0**8.0
    nltmin = 0.0

    DO i = 0,nkx0-1
       DO j = 0,nky0-1
          DO k = lkz1,lkz2

             turnoverm = turnoverm + sum(abs(v_1(i,j,k,:))**2.0) * sin(kxgrid(i)/(2.0*kxmin))**2.0
             nltmax = min(nltmax,minval(kmags(i,j,k)*abs(v_1(i,j,k,:)),kmags(i,j,k)*abs(v_1(i,j,k,:)).gt.10.0**(-10.0)))
             nltmin = max(nltmin,kmags(i,j,k)*maxval(abs(v_1(i,j,k,:)),abs(v_1(i,j,k,:)).gt.10.0**(-10.0)))

          ENDDO
       ENDDO
    ENDDO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(turnoverm,turnover,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

    turnover = 1.0/(2*kxmin*sqrt(turnover))
    
    print *, mype,"Turnover Time Estimate",turnover

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(nltmax,nlt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)

    nlt = 10.0/nlt

    print *, mype,"Equation-Based Max Nonlinear Time Scale", nlt

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(nltmin,nlt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)

    nlt = 10.0/nlt
    
    print *, mype,"Equation-Based Min Nonlinear Time Scale", nlt
    
    if (nv) b_1(:,:,:,:) = cmplx(0.0,0.0)
    
    if (taylorgreen) then
       ! Initializes Fourier components of Taylor Green vortex u = sin x cos y cos z , v = - cos x sin y cos z
       v_1 = cmplx(0.0,0.0)
       if (mype.eq.0) then 
          v_1(1,1,0,0) = 1.0
          v_1(1,nky0-1,0,0) = 1.0
          v_1(1,1,0,1) = -1.0
          v_1(1,nky0-1,0,1) = 1.0
       endif
    
       v_1 = -v_1 * cmplx(0.0,0.125)
    endif

    ! Scale forcing to represent energy gained per time step |F| |v| ~ energy gain per step as fraction of initial energy
    ! |v| ~ sqrt(|v|^2) ~ sqrt(energy/2)
    
    ! if (force_turbulence) force_amp = force_amp * sqrt(4.0 * pi**3.0 * init_energy)
    ! if (verbose) print *, mype,"Force amp",force_amp
    
 ENDIF
 IF (checkpoint_read) CALL checkpoint_in
   
 ! Linear stability maximum time step
    print *, "Gauss2 Critical Time Step", 2.0/(maxval(kzgrid)*(kmax/2 + sqrt(1 + 0.25*kmax**2.0)))
    if (calc_dt.and.(.not.(test_ho)).and.(hall.ne.0.0)) dt_max = minval([dt_max,1.0/(maxval(kzgrid)*(kmax/2 + sqrt(1 + 0.25*kmax**2.0)))])
    if (verbose.and.(mype.eq.0)) then
       print *, "kzgrid max", maxval(kzgrid)
       print *, "kmags max", kmax
    endif

      ! Magnetic Helicity Correction
    mhelcorr = 0.0


    ! Check on Initial Energy
    knzeroenm = sum(abs(b_1(1:nkx0-1,:,:,:))**2+abs(v_1(1:nkx0-1,:,:,:))**2)
    kxzeroenm = sum(0.5*(abs(b_1(0,:,:,:))**2+abs(v_1(0,:,:,:))**2))

    print *, mype, "Mype Energy",8*pi**3 * (knzeroenm+kxzeroenm)
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(knzeroenm,knzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(kxzeroenm,kxzeroen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    s1 = knzeroen + kxzeroen

    if (mype.eq.0) print *, "All Mype Initial Energy",s1*8*pi**3

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


