!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                               cc_calc_dt.f90                              !!
!!                                                                           !!
!!  calculate_time_step                                                      !!
!!  -- calc_dt_slepc                                                         !!
!!  -- compute_stability_criterion                                           !!
!!  -- zroots                                                                !!
!!  -- zroots_complex                                                        !!
!!  -- laguer                                                                !!
!!  -- initialize_adapt_dt                                                   !!
!!  -- adapt_dt                                                              !!
!!  -- finalize_adapt_dt                                                     !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                             calculate_time_step                           !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Calculates initial time step
!!  Also includes routines for adaptive time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE calculate_time_step
  USE par_mod
  USE mpi
#ifdef WITHSLEPC
  USE eigen_iterative, ONLY: calc_dt_ev, current_evs
#endif
  USE eigen_direct, ONLY: calc_dt_ev_lapack
  USE communications
  USE nonlinearity, ONLY: ve_max

  PUBLIC :: compute_stability_criterion,initialize_adapt_dt,adapt_dt,&
            finalize_adapt_dt, calc_initial_dt 

  PRIVATE

  !quantities needed for nonlinear time step adaption
  INTEGER:: win_dts=100, win_est=25
  REAL,DIMENSION(:),ALLOCATABLE:: last_dts, last_estimates
  REAL:: ltav
  REAL :: rk_corr(2)
  !INTEGER:: TIMEFILE=78


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             calc_dt!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_initial_dt

  IF(dt_slepc) THEN
    IF(mype==0) WRITE(*,*) "Caculating dt_max with SLEPC."
#ifdef WITHSLEPC
    CALL calc_dt_slepc
#else
    stop "Must link with SLEPC in order to use calc_dt_slepc!"
#endif
  ELSE
    IF(mype==0) WRITE(*,*) "Caculating dt_max with LAPACK."
    CALL calc_dt_lapack
  END IF

END SUBROUTINE calc_initial_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             calc_dt_lapack                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_dt_lapack

  INTEGER :: i,j,k,l
  INTEGER :: nkx0_bak,nky0_bak,nkz0_bak
  INTEGER :: np_herm_bak,np_kz_bak,np_spec_bak,np_hank_bak
  LOGICAL :: nonlinear_bak
  REAL :: kxmin_bak,kymin_bak,kzmin_bak
  REAL :: kxmax_bak,kymax_bak,kzmax_bak
  !REAL :: dt_min
  !REAL :: kx_dtmin,ky_dtmin,kz_dtmin
  LOGICAL :: kmin_eq_0_bak,evenyz_bak
  COMPLEX :: max_omega_final
  REAL :: dt_temp
  REAL :: kx_temp(3)
  REAL :: ky_temp(3)
  REAL :: kz_temp(3)
  COMPLEX :: phase
  REAL :: dt_lim,angle
  CHARACTER(len=1) :: chnum1,chnum2,chnum3
  COMPLEX :: current_evs(nv0*nh0*nspec)
  INTEGER :: ierr

  IF(mype==0) WRITE(*,*) "Calculating time step with LAPACK."

  max_omega_final=cmplx(0.0,0.0)
  nonlinear_bak=nonlinear
  nonlinear=.false.

  kmin_eq_0_bak=kmin_eq_0
  kmin_eq_0=.false.

  evenyz_bak=evenyz
  evenyz=.false.

  kx_temp(1)=0.0
  kx_temp(2)=kxmin
  kx_temp(3)=kxmin*hkx_ind
  ky_temp(1)=0.0
  ky_temp(2)=kymin
  ky_temp(3)=kymin*(hky_ind)
  IF(spatial2d) THEN
    kz_temp(1)=0.0
    kz_temp(2)=kzmin
    kz_temp(3)=kzmin*(hky_ind)
  ELSE
    kz_temp(1)=0.0
    kz_temp(2)=kzmin
    kz_temp(3)=kzmin*(hkz_ind)
  END IF

  nkx0_bak=nkx0 
  nky0_bak=nky0 
  nkz0_bak=nkz0 

  nkx0=1
  nky0=1
  nkz0=1

  kxmin_bak=kxmin 
  kymin_bak=kymin 
  kzmin_bak=kzmin 
  kxmax_bak=kxmax 
  kymax_bak=kymax 
  kzmax_bak=kzmax 

  np_herm_bak=np_herm
  np_herm=1
  np_hank_bak=np_hank
  np_hank=1
  np_kz_bak=np_kz
  np_kz=1
  np_spec_bak=np_spec
  np_spec=1

  n_ev=nv0*nh0*nspec

  dt_max=100.0

  CALL finalize_arrays

  !write(*,*) "Before call comm", mype
  !CALL comm
  !write(*,*) "After call comm", mype

  DO i=1,3
    !IF(mype==0) WRITE(*,*) "i=",i
    DO j=1,3
      DO k=1,3
       IF((i.ne.1).or.(j.ne.1)) THEN
        IF(mype==0) WRITE(chnum1,'(i1.1)') i
        IF(mype==0) WRITE(chnum2,'(i1.1)') j
        IF(mype==0) WRITE(chnum3,'(i1.1)') k
        IF(mype==0) OPEN(unit=100,file=trim(diagdir)//'/dtevl_'//chnum1//chnum2//chnum3,status='unknown')
        kxmin=kx_temp(i)
        kymin=ky_temp(j)
        kzmin=kz_temp(k)
        !IF(mype==0) WRITE(*,*) kxmin,kymin,kzmin 
        CALL arrays_temp
        kxmax=kxmax_bak
        kymax=kymax_bak
        kzmax=kzmax_bak
        !WRITE(*,*) "Before CALL calc_dt_ev",mype
        CALL calc_dt_ev_lapack(current_evs)
        !IF(mype==0) WRITE(*,*) current_evs
        dt_temp=100.0
        !WRITE(*,*) "Before CALL compute_stability_criterion",mype
        DO l=1,nv0*nh0*nspec
         IF(mype==0) WRITE(100,*) REAL(current_evs(l)),aimag(current_evs(l))
         IF(REAL(current_evs(l)).ne.0.0.or.aimag(current_evs(l)).ne.0.0) THEN
           IF(REAL(current_evs(l)).lt.0.0) CALL compute_stability_criterion(current_evs(l),dt_temp)
         END IF
         IF(dt_temp.lt.dt_max) THEN
          dt_max=dt_temp
          IF(mype==0) WRITE(*,*) "New dt_max:",dt_max
          IF(mype==0) WRITE(*,*) "Omega",current_evs(l)
          IF(mype==0) WRITE(*,*) "kx,ky,kz",kxmin,kymin,kzmin
         END IF
        END DO
        CALL finalize_arrays 
        IF(mype==0) CLOSE(100)
       END IF
      END DO
    END DO
  END DO

  dt_max=dt_max*0.95
  CALL MPI_ALLREDUCE(dt_max,dt_temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM , MPI_COMM_WORLD,ierr)
  dt_max=dt_temp/float(n_mpi_procs)
  IF(mype==0) WRITE(*,*) "dt_max=",dt_max

  nkx0=nkx0_bak
  nky0=nky0_bak
  nkz0=nkz0_bak

  kxmin=kxmin_bak
  kymin=kymin_bak
  kzmin=kzmin_bak

  nonlinear=nonlinear_bak
  kmin_eq_0=kmin_eq_0_bak
  evenyz=evenyz_bak

  np_herm=np_herm_bak
  np_hank=np_hank_bak
  np_kz=np_kz_bak
  np_spec=np_spec_bak

  !CALL finalize_comm
  !CALL comm

  CALL arrays

END SUBROUTINE calc_dt_lapack

  

#ifdef WITHSLEPC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             calc_dt_slepc                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_dt_slepc

  INTEGER :: i,j,k,l
  INTEGER :: nkx0_bak,nky0_bak,nkz0_bak
  LOGICAL :: nonlinear_bak
  REAL :: kxmin_bak,kymin_bak,kzmin_bak
  REAL :: kxmax_bak,kymax_bak,kzmax_bak
  !REAL :: dt_min
  !REAL :: kx_dtmin,ky_dtmin,kz_dtmin
  LOGICAL :: kmin_eq_0_bak,evenyz_bak
  COMPLEX :: max_omega_final
  REAL :: dt_temp
  REAL :: kx_temp(3)
  REAL :: ky_temp(3)
  REAL :: kz_temp(3)
  COMPLEX :: phase
  REAL :: dt_lim,angle
  CHARACTER(len=1) :: chnum1,chnum2,chnum3

  IF(mype==0) WRITE(*,*) "Calculating time step with SLEPc."

  max_omega_final=cmplx(0.0,0.0)
  nonlinear_bak=nonlinear
  nonlinear=.false.

  kmin_eq_0_bak=kmin_eq_0
  kmin_eq_0=.false.

  evenyz_bak=evenyz
  evenyz=.false.

  kx_temp(1)=0.0
  kx_temp(2)=kxmin
  kx_temp(3)=kxmin*hkx_ind
  ky_temp(1)=0.0
  ky_temp(2)=kymin
  ky_temp(3)=kymin*(hky_ind)
  kz_temp(1)=0.0
  kz_temp(2)=kzmin
  kz_temp(3)=kzmin*(hkz_ind)

  nkx0_bak=nkx0 
  nky0_bak=nky0 
  nkz0_bak=nkz0 

  nkx0=1
  nky0=1
  nkz0=1

  kxmin_bak=kxmin 
  kymin_bak=kymin 
  kzmin_bak=kzmin 
  kxmax_bak=kxmax 
  kymax_bak=kymax 
  kzmax_bak=kzmax 

  n_ev=nv0

  dt_max=100.0

  CALL finalize_arrays
  ALLOCATE(current_evs(nv0))

  DO i=1,3
    !IF(mype==0) WRITE(*,*) "i=",i
    DO j=1,3
      DO k=1,3
       IF((i.ne.1).or.(j.ne.1)) THEN
        WRITE(chnum1,'(i1.1)') i
        WRITE(chnum2,'(i1.1)') j
        WRITE(chnum3,'(i1.1)') k
        OPEN(unit=100,file=trim(diagdir)//'/dtevs_'//chnum1//chnum2//chnum3,status='unknown')
        kxmin=kx_temp(i)
        kymin=ky_temp(j)
        kzmin=kz_temp(k)
        !IF(mype==0) WRITE(*,*) kxmin,kymin,kzmin 
        CALL arrays
        kxmax=kxmax_bak
        kymax=kymax_bak
        kzmax=kzmax_bak
        !WRITE(*,*) "Before CALL calc_dt_ev",mype
        CALL calc_dt_ev
        dt_temp=100.0
        !WRITE(*,*) "Before CALL compute_stability_criterion",mype
        DO l=1,nv0
         WRITE(100,*) REAL(current_evs(l)),aimag(current_evs(l))
         IF(REAL(current_evs(l)).ne.0.0.or.aimag(current_evs(l)).ne.0.0) THEN
           IF(REAL(current_evs(l)).lt.0.0) CALL compute_stability_criterion(current_evs(l),dt_temp)
         END IF
         IF(dt_temp.lt.dt_max) THEN
          dt_max=dt_temp
          IF(mype==0) WRITE(*,*) "New dt_max:",dt_max
          IF(mype==0) WRITE(*,*) "Omega",current_evs(l)
          IF(mype==0) WRITE(*,*) "kx,ky,kz",kxmin,kymin,kzmin
         END IF
        END DO
        CALL finalize_arrays 
        CLOSE(100)
       END IF
      END DO
    END DO
  END DO

  dt_max=dt_max*0.95
  IF(mype==0) WRITE(*,*) "dt_max=",dt_max

  DEALLOCATE(current_evs)
  nkx0=nkx0_bak
  nky0=nky0_bak
  nkz0=nkz0_bak

  kxmin=kxmin_bak
  kymin=kymin_bak
  kzmin=kzmin_bak

  nonlinear=nonlinear_bak
  kmin_eq_0=kmin_eq_0_bak
  evenyz=evenyz_bak

  CALL arrays

  OPEN(unit=100,file=trim(diagdir)//'/rk_stability.dat',status='unknown')
  DO i=0,500
    angle=2.0*pi*(i-1)/(500.0)
    phase=cmplx(cos(angle),sin(angle))
    phase=phase/abs(phase)
    IF(REAL(phase).lt.0.0) THEN
        CALL compute_stability_criterion(phase,dt_lim) 
        WRITE(100,*) REAL(phase)*dt_lim/dt_max,aimag(phase)*dt_lim/dt_max,dt_lim
    END IF
  END DO
  CLOSE(100)


END SUBROUTINE calc_dt_slepc
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                         compute_stability_criterion                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  From RK_coefficients.F90 in GENE: Computes the criterion for neutral stability !!
!!  The linearized GK equation for an eigenvector \f$|\lambda\rangle\f$ can be written as 
!!  \f$\partial_t |\lambda\rangle=L|\lambda\rangle=\lambda |\lambda\rangle\f$ (L is the linear operator, 
!!  \f$\lambda\f$) the eigenvalue. The time evolution is approximated by a Runge-Kutta scheme, 
!!  \f$|\lambda(t_0+\Delta t)\rangle=\sum_n c_n \Delta t^n (\partial_t^n |\lambda\rangle)|_{t_0}=
!!  \sum_n c_n \Delta t^n \lambda^n |\lambda(t_0)\rangle\f$, where the \f$c_n\f$ are the coefficients of the 
!!  dt expansion corresponding to the RK scheme.
!!  The system is neutrally stable IF the amplitude of \f$|\lambda\rangle\f$ remains constant, i.e. IF
!!  \f$|\sum_n c_n \Delta t^n \lambda^n|^2=1\f$. This is a polynomial equation in \f$(\Delta t \lambda)\f$. 
!!  The 1 on the right cancels with the n=0 contribution on the left, the remaining equation can be devided by
!!  \f$(\Delta t \lambda)\f$ because the trivial solution 0 is not of interest. The remainder is a polynomial of 
!!  order (2*rkstages-1) with three (REAL) degrees of freedom: \f$\Delta t, Re(\lambda), Im(\lambda)\f$. IF 
!!  \f$\lambda\f$ is specified, the (maximal) \f$\Delta t\f$ can be computed, IF \f$\Delta t=1,Re(\lambda)=0\f$ are
!!  specified, the root of the polynomial corresponds to the stability interval on the imaginary axis
!!  (i.e. for advection problems) which is used for (approximate) dt estimates.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compute_stability_criterion(cnum,res)
    COMPLEX,INTENT(in):: cnum
    REAL,INTENT(out):: res
    REAL:: re,im
    INTEGER:: nroots, rnum
    REAL, DIMENSION(:), ALLOCATABLE:: poly_ev
    REAL(8), DIMENSION(:), ALLOCATABLE:: poly_ev_nag
    COMPLEX(8), DIMENSION(:), ALLOCATABLE:: droots
    COMPLEX, DIMENSION(:), ALLOCATABLE:: roots
    REAL,DIMENSION(0:11):: arr
    REAL,DIMENSION(0:6)::b    
    REAL,DIMENSION(:),ALLOCATABLE:: a_rk, b_rk, RK_Taylor
    INTEGER:: rkstages

    !For RK4
    rkstages=4
    ALLOCATE(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
    a_rk=(/0.5, 0.5, 1.0/)
    b_rk=(/1.0, 2.0, 2.0, 1.0/)/6.
    RK_Taylor=(/1., 1., 1./2, 1./6, 1./24/)

    b=0.
    b(0:rkstages)=RK_Taylor    
    nroots=2*rkstages-1
    ALLOCATE(poly_ev(2*rkstages), poly_ev_nag(2*rkstages), droots(2*rkstages-1), roots(2*rkstages-1))
    res=-1.
    re=REAL(cnum)
    im=Aimag(cnum)
    !this is stability polynomial for up to six stages
    arr=(/2*b(1)*re, &
         2*b(2)*(-im**2+re**2)+b(1)**2*(im**2+re**2),& 
         2*re*(b(3)*(-3*im**2+re**2)+b(1)*b(2)*(im**2+re**2)),& 
         b(2)**2*(im**2+re**2)**2+2*(b(1)*b(3)*(-im**4+re**4)+b(4)*(im**4-6*im**2*re**2+re**4)),& 
         2*re*(b(2)*b(3)*(im**2+re**2)**2+b(5)*(5*im**4-10*im**2*re**2+re**4)+& 
         b(1)*b(4)*(-3*im**4-2*im**2*re**2+re**4)), &
         b(3)**2*(im**2+re**2)**3+2*(-(b(2)*b(4)*(im**2-re**2)*(im**2+re**2)**2)+&
         b(6)*(-im**6+15*im**4*re**2-15*im**2*re**4+re**6)+b(1)*b(5)*(im**6-5*im**4*re**2-5*im**2*re**4+re**6)),&
         2*re*(im**2+re**2)*(b(3)*b(4)*(im**2+re**2)**2+b(1)*b(6)*(5*im**4-10*im**2*re**2+re**4)+&
         b(2)*b(5)*(-3*im**4-2*im**2*re**2+re**4)), &
         (im**2+re**2)**2*(b(4)**2*(im**2+re**2)**2+2*(b(3)*b(5)*(-im**4+re**4)+&
         b(2)*b(6)*(im**4-6*im**2*re**2+re**4))),&
         2*re*(im**2+re**2)**3*(b(3)*b(6)*(-3*im**2+re**2)+b(4)*b(5)*(im**2+re**2)), &
         (im**2+re**2)**4*(2*b(4)*b(6)*(-im**2+re**2)+b(5)**2*(im**2+re**2)),& 
         2*b(5)*b(6)*re*(im**2+re**2)**5, &
         b(6)**2*(im**2+re**2)**6/)
    poly_ev=arr(0:nroots)
    
    !find the roots of the polynomial
    CALL zroots(poly_ev,nroots,roots,.true.)
    
    !return only the maximal REAL+positive solutions. For cnum=\lambda, his corresponds to the maximal dt 
    !allowed for the input eigenvector; for cnum=(0.,1.)/(-1.,0.) to the extent of the stabile region along 
    !the imaginary /negative realaxis for dt=1.
    DO rnum=1,nroots
       IF (aimag(roots(rnum)).le.1e-8.and.REAL(roots(rnum)).gt.0.0) res=MAX(res,REAL(roots(rnum)))
    enddo
    DEALLOCATE(poly_ev, poly_ev_nag, roots, droots)

  END SUBROUTINE compute_stability_criterion


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 zroots                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  zroots for REAL input  (from GENE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE zroots(a,m,roots,polish)
    INTEGER :: m  !<degree of the polynomial
    !>the m+1 COMPLEX coefficients of the polynomial in ascending
    !! order. a1+a2*z+a3*z^2+...+a(m+1)*z^(m)
    REAL    :: a(m+1)
    COMPLEX :: roots(m) !<result
    LOGICAL :: polish !<enables polishing (whatsoever that is...)
    COMPLEX :: a_complex(m+1)
    
    a_complex = a
    CALL zroots_complex(a_complex,m,roots,polish)

  END SUBROUTINE zroots


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              zroots_complex                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  zroots for COMPLEX input  (from GENE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE zroots_complex(a,m,roots,polish)
    INTEGER :: m  !<degree of the polynomial
    !>the m+1 COMPLEX coefficients of the polynomial in ascending
    !! order. a1+a2*z+a3*z^2+...+a(m+1)*z^(m)
    COMPLEX :: a(m+1)
    COMPLEX :: roots(m) !<result
    LOGICAL :: polish !<enables polishing (whatsoever that is...)
    REAL,PARAMETER :: EPS=1.0e-6
    INTEGER,parameter :: MAXM=101

    INTEGER :: i,j,jj,its
    COMPLEX :: ad(MAXM),x,b,c

    DO j=1,m+1
       ad(j)=a(j)
    END DO

    DO j=m,1,-1
       x=CMPLX(0.,0.)
       CALL laguer(ad,j,x,its)
       IF(ABS(AIMAG(x)).LE.2.*EPS**2*ABS(REAL(x))) x=CMPLX(REAL(x),0.)
       roots(j)=x
       b=ad(j+1)
       DO jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
       END DO
    END DO
    
    IF (polish) THEN
       DO j=1,m
          CALL laguer(a,m,roots(j),its)
       END DO
    END IF
    
    DO j=2,m
       x=roots(j)
       DO i=j-1,1,-1
          IF(REAL(roots(i)).LE.REAL(x)) GOTO 10
          roots(i+1)=roots(i)
       END DO
       i=0
10     roots(i+1)=x
     END DO
   END SUBROUTINE zroots_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  laguer                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE laguer(a,m,x,its)
     INTEGER :: m,its
     COMPLEX :: a(m+1),x

     REAL,PARAMETER :: EPSS=2.0e-7
     INTEGER,PARAMETER :: MR=8,MT=10,MAXIT=MT*MR

     INTEGER :: iter,j
     REAL :: abx,abp,abm,err,frac(MR)=(/.5,.25,.75,.13,.38,.62,.88,1./)
     COMPLEX :: dx,x1,b,d,f,g,h,sq,gp,gm,g2
     SAVE frac

     DO iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=ABS(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        DO j=m,1,-1
           f=x*f+d
           d=x*d+b
           b=x*b+a(j)
           err=ABS(b)+abx*err
        END DO
        err=EPSS*err
        IF(ABS(b).LE.err) THEN
           RETURN
        ELSE
           g=d/b
           g2=g*g
           h=g2-2.*f/b
           sq=SQRT((m-1)*(m*h-g2))
           gp=g+sq
           gm=g-sq
           abp=ABS(gp)
           abm=ABS(gm)
           IF(abp.LT.abm) gp=gm
           IF (MAX(abp,abm).GT.0.) THEN
              dx=m/gp
           ELSE
              dx=EXP(CMPLX(LOG(1.+abx),float(iter)))
           END IF
        END IF
        x1=x-dx
        IF(x.EQ.x1) RETURN
        IF (MOD(iter,MT).NE.0) THEN
           x=x1
        ELSE
           x=x-dx*frac(iter/MT)
        END IF
     END DO
     STOP 'too many iterations in polyroots/laguer'
   END SUBROUTINE laguer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           initialize_adapt_dt                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  The adapt_dt routines are modified from the GENE versions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initialize_adapt_dt

    ALLOCATE(last_dts(0:win_dts-1),last_estimates(0:win_est-1))
    IF(checkpoint_read) THEN
      last_dts=dt
      last_estimates=dt
    ELSE
      last_dts=dt_max
      last_estimates=dt_max
    END IF
    ltav=0.
    !IF(mype.eq.0) OPEN(TIMEFILE, file=trim(diagdir)//'/timeest', form='formatted',&
    !        status='replace', position='REWIND')
    CALL compute_stability_criterion((0.,1.),rk_corr(1))
    CALL compute_stability_criterion((-1.,0.),rk_corr(2))

  END SUBROUTINE initialize_adapt_dt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                adapt_dt                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE adapt_dt
    REAL:: nl_dt, dt_inst, dt_est, dt_all
    REAL, DIMENSION(win_est):: sortdts
    REAL, DIMENSION(win_est-2):: deri
    REAL:: av, var, normvar
    !LOGICAL :: dt_output

    !linear part is neglected (and compensated by the courant factor),
    !IF it is added here the time step is way too small
    nl_dt=rk_corr(1)/(ve_max(1)*kxmax+ve_max(2)*kymax)
    IF(dt_output.and.mype==0) WRITE(*,*) "nl_dt",nl_dt

    last_estimates(MOD(itime-1,win_est)) = nl_dt
    dt_est = courant*nl_dt
    IF(dt_output.and.mype==0) WRITE(*,*) "dt_est",dt_est

    !some statistical analysis to detect instabilities
    sortdts=cshift(last_estimates,MOD(itime-1,win_est)-win_est+1)
    !2nd derivative
    deri=2*sortdts(2:win_est-1)-sortdts(3:win_est)-sortdts(1:win_est-2)
   
    !variance of the sample
    av=sum(deri)/(win_est-2)
    var=sum((deri-av)**2)/(win_est-3)

    !normalize variance to the old average
    IF (((itime.le.win_est+1).or.(ltav.lt.(1000*epsilon(ltav))))) THEN
       normvar=0. 
    ELSE
       normvar=var/ltav
    END IF

    !limit the reduction factor of the variance term to 0.3
    IF (normvar.gt.50) normvar=50

    !update long term average
    IF (itime.gt.win_est) THEN
       IF(itime.gt.win_est+500) THEN
          ltav=0.998*ltav+0.002*var
       ELSE
          ltav=((itime-win_est-1)*ltav+var)/(itime-win_est)
       END IF
    END IF

    !combine linear, nonlinear and variance constraint
    IF(checkpoint_read.and.(itime-itime_start).lt.win_dts) THEN
      dt_all=MIN(dt_est,dt_max)
    ELSE
      dt_all=MIN(dt_est,dt_max)/cosh(normvar/4)**0.1
    END IF
    IF(dt_output.and.mype==0) WRITE(*,*) "dt_all",dt_all

    last_dts(MOD(itime-1,win_dts))=dt_all

    !IF(dt_output.and.mype==0) WRITE(*,*) "last_dts",last_dts
    dt_inst=minval(last_dts)
    IF(dt_output.and.mype==0) WRITE(*,*) "dt_inst",dt_inst

    IF (dt_inst/dt.lt.0.8) THEN
       dt=dt_inst
    ELSE IF(dt_inst/dt.gt.1.2) THEN
       !ramp up dt only gradually
       dt=1.2*dt
       !modify last_dts to store this value for win_dts timesteps
       last_dts(MOD(itime-1,win_dts))=1.2*dt
    END IF
  
    IF(dt_output.and.mype==0) WRITE(*,*) "dt=",dt
    !IF(dt_output.and.mype==0) WRITE(*,*) "itime=",itime

  END SUBROUTINE adapt_dt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             finalize_adapt_dt                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE finalize_adapt_dt
  
    IF (mype==0.and.verbose) WRITE(*,*) "Finalizing adapt_dt."
    !IF(mype.eq.0) CLOSE(TIMEFILE)
    IF(ALLOCATED(last_dts)) DEALLOCATE(last_dts)
    IF(ALLOCATED(last_estimates)) DEALLOCATE(last_estimates)
    
  END SUBROUTINE finalize_adapt_dt


END MODULE calculate_time_step
