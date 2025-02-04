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

  use iso_fortran_env, only: int64
  
  IMPLICIT NONE

  INTEGER :: i,j,k, l,zst,xst,yst,ri,ind
  REAL :: kfactor,err
  REAL :: s1, s11,s12,s13,s4
  !REAL :: init_prefactor
  COMPLEX :: phase,phaseb,phasev,phaseby,phasevy
  COMPLEX :: LW1, LC1, RW1, RC1
  REAL :: phase1,phase2,phase1y,phase2y,kspect,myphase1,myphase1y,myphase2,myphase2y,showphase,thb,thv,mythb,mythv,bt1,bt2,bt3,btmag,b1w,b2w,v1w,v2w
  REAL :: b1r,b1i,b2r,b2i,b3r,b3i
  CHARACTER(len=40), INTENT(in) :: which_init0
  CHARACTER(len=40) :: which_init
  REAL :: zerocmplx
  REAL :: kxzeroen,knzeroen,kxzeroenm,knzeroenm

  INTEGER, DIMENSION(:), ALLOCATABLE :: rseed
  INTEGER :: rseed_size,ierr
  REAL :: testrandoms(10)
  REAL :: truncx,truncy,truncz,nlt
  REAL :: fperpm,fperp
  REAL(C_DOUBLE) :: dt_critical

  REAL :: turnoverm,nltmin,nltmax ! MPI ALL REDUCE variables
  INTEGER(int64) :: t 

  zerocmplx=0.0
   
  !g_1(:,:,:,:,:,:)=cmplx(0.0,0.0)
  b_1(:,:,:,:)=cmplx(0.0,0.0)
  v_1(:,:,:,:)=cmplx(0.0,0.0)

  xst = 1
  yst = 1
  zst = max(1,lkz1)

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

      fperp = (force_frac * nky0)**3.0
      
      ! force_amp = abs(cmplx(vnu*(kxmin)**(2.0*hyp),0)) * sqrt(1.0/3.0 * force_amp) * sqrt(8.0 * (pi ** 3))/fperp
      
      print *, "Force Amp",force_amp      

      DO i=xst,nx0_big/2
         DO j=yst,ny0_big-1
            DO k=zst,lkz2
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
                     s13 = kxgrid(i) * kygrid(j) * (sin(4*pi*thb) * cos(2*pi*phase1y) &
                          + sin(4*pi*thv)*cos(2*pi*phase2y))
                     s1 = s1 + kperps(i,j,k)**(-1.0*init_kolm) * (2+(s11+s12+s13)/(kzgrid(k)**2))
                  ELSE
                     s1 = s1 + 2.0*kperps(i,j,k) ** (-1.0*init_kolm)
                  ENDIF
                  
                  IF (helical.and.(max(i,j,k).ne.0)) THEN
                     bt1 = 2.0*phase1y - 1.0
                     bt2 = 2.0*phase2 - 1.0
                     bt3 = 2.0*phase2y - 1.0
                     b1r = bt1 - (kxgrid(i) * bt1 + kygrid(j) * bt2 &
                          + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kxgrid(i)
                     b2r = bt2 - (kxgrid(i) * bt1 + kygrid(j) * bt2 &
                          + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kygrid(j)
                     b3r = bt3 - (kxgrid(i) * bt1 + kygrid(j) * bt2 &
                          + kzgrid(k) * bt3)/(kmags(i,j,k)**2) * kzgrid(k)
                     btmag = sqrt(b1r**2 + b2r**2 + b3r**2)
                     b1i = (kygrid(j) * b3r - kzgrid(k) * b2r)/kmags(i,j,k)
                     b2i = (kzgrid(k) * b1r - kxgrid(i) * b3r)/kmags(i,j,k)
                     b3i = (kxgrid(i) * b2r - kygrid(j) * b1r)/kmags(i,j,k)
                     b_1(i,j,k,0) = cmplx(b1r,b1i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                     b_1(i,j,k,1) = cmplx(b2r,b2i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                     b_1(i,j,k,2) = cmplx(b3r,b3i)/(sqrt(2.0) * btmag * kperps(i,j,k)**(init_kolm/2.0))
                     v_1(i,j,k,:) = b_1(i,j,k,:)
                  ELSE IF (shear.and.(max(i,j).ne.0)) THEN
                     b_1(i,j,k,0) = - kygrid(j)/sqrt(kxgrid(i)**2 &
                          + kygrid(j)**2) * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                     b_1(i,j,k,1) = kxgrid(i)/sqrt(kxgrid(i)**2 &
                          + kygrid(j)**2) * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                     b_1(i,j,k,2) = cmplx(0.0,0.0)
                     v_1(i,j,k,0) = -kygrid(j)/sqrt(kxgrid(i)**2 &
                          + kygrid(j)**2) * phasev * phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
                     v_1(i,j,k,1) = kxgrid(i)/sqrt(kxgrid(i)**2 &
                          + kygrid(j)**2) * phasev*phaseb * (kperps(i,j,k) **(-init_kolm/2.0))
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
               
             ! endif

             if (force_trunc) then !Only initalize forced waves
                if (kperps(i,j,k) > (sqrt(kxmax**2.0 + kymax**2.0)*force_frac)) then
                   b_1(i,j,k,:) = cmplx(0.0,0.0)
                   v_1(i,j,k,:) = cmplx(0.0,0.0)
                end if
             endif
             
          ENDDO
       ENDDO
    ENDDO
    ! Filter modes in the padding region
    DO ind = 0,2
       b_1(:,:,:,ind) = b_1(:,:,:,ind) * paddingmask
       v_1(:,:,:,ind) =	v_1(:,:,:,ind) * paddingmask
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

    if (verbose) then
       DO ind = 0,2
          print *, "Max Bind",mype,maxval(abs(b_1(:,:,:,ind)))
          print *, "Max Vind",mype,maxval(abs(v_1(:,:,:,ind)))
       ENDDO
    endif
    
    

    if (verbose) print *, "Through energy normalization"

    turnoverm = 0.0
    nltmax = 10.0**8.0
    nltmin = 0.0

    DO i = 0,nkx0-1
       DO j = 0,ny0_big-1
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

    IF (init_null) THEN
       b_1 = cmplx(0.0,0.0)
       v_1 = cmplx(0.0,0.0)
    ENDIF
    
    IF (checkpoint_read) CALL checkpoint_in
   
    ! Linear stability maximum time step
    dt_critical = 2.0/(maxval(kzgrid)*(maxval(kmags)/2 + sqrt(1 + 0.25*maxval(kmags)**2.0)))
    print *, "Gauss2 Critical Time Step", dt_critical
    if (calc_dt.and.(.not.(test_ho)).and.(hall.ne.0.0)) dt_max = minval([dt_max,dt_critical/2.0])
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

    dt = dt_max

END SUBROUTINE initial_condition
