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
  
END SUBROUTINE init_run



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                   arrays                                  !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Sets up indices and allocates important arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE arrays
  USE nonlinearity, ONLY: initialize_fourier
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

  ! Adjust for P3DFFT

  hkx_ind = hkx_ind+1
  hky_ind = hky_ind+1
  hkz_ind = hkz_ind+1
  lky_ind = lky_ind+1
  lkz_ind = lkz_big+1
  lky_big = lky_big+1
  lkz_big = lkz_big+1  

  IF(verbose.and.(mype.eq.0)) THEN
    WRITE(*,*) "hkx_ind",hkx_ind
    WRITE(*,*) "hky_ind",hky_ind
    WRITE(*,*) "lky_ind",lky_ind
    WRITE(*,*) "hkz_ind",hkz_ind
    WRITE(*,*) "lkz_ind",lkz_ind
 END IF
 
 ! Initialize Fourier needed for P3DFFT pencil decomposition
 CALL initialize_fourier
 IF (verbose.and.(mype.eq.0)) WRITE(*,*) "Called initial fourier.",mype 
 
  !Info for kz grid
  lkz0=nz0_big

  if (verbose) print *, mype,"Mype Limits x",cstart(1),cend(1)
  if (verbose) print *, mype,"Mype Limits y",cstart(2),cend(2)
  if (verbose) print *, mype,"Mype Limits z",cstart(3),cend(3)

  ! complex kz is not distributed 

  ! Upper and lower bounds nonzero y
  if (n_mpi_procs.eq.1) then
     lbky = cstart(2)
     ubky = cend(2)
  else if (cstart(2).le.hky_ind) then
     lbky = cstart(2)
     ubky = min(hky_ind,cend(2))
  else if (cend(2).ge.lky_big) then
     lbky = max(cstart(2),lky_ind)
     ubky = cend(2)
  else
     lbky = cend(2)
     ubky = cstart(2)
  endif

  ! Upper and lower bounds nonzero x
  if (n_mpi_procs.eq.1) then
     lbkx = cstart(1)
     ubkx = cend(1)
  else if (cstart(1).le.hkx_ind) then
     lbkx = cstart(1)
     ubkx = min(cend(1),hkx_ind)
  else
     lbkx = cend(1)
     ubkx = cstart(1)
  endif
  

  ! Define mask for zero padding
  IF (.not.allocated(paddingmask)) ALLOCATE(paddingmask(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  paddingmask = 0
  
  if (n_mpi_procs.eq.1) then
     ! In the one process case hkz_ind and lkz_ind both lie between lkz1 and lkz2
     paddingmask(1:nkx0,1:hky_ind,1:hkz_ind) = 1
     paddingmask(1:nkx0,1:hky_ind,lkz_big:nz0_big) = 1
     paddingmask(1:nkx0,lky_big:ny0_big,lkz_big:nz0_big) = 1
     paddingmask(1:nkx0,lky_big:ny0_big,1:hkz_ind) = 1
  else
     ! Use choice of lbs and ubs to select nonzero padding
     paddingmask(lbkx:ubkx,lbky:ubky,1:hkz_ind) = 1
     paddingmask(lbkx:ubkx,lbky:ubky,lkz_big:nz0_big) = 1
  endif

  IF(.not.allocated(reader)) ALLOCATE(reader(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  IF(.not.allocated(reader2)) ALLOCATE(reader2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  IF (.not.allocated(fullsmallarray)) ALLOCATE(fullsmallarray(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  IF (.not.allocated(fullbigarray)) ALLOCATE(fullbigarray(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

  IF (.not.allocated(scatter_big)) ALLOCATE(scatter_big(product(csize)))
  IF (.not.allocated(scatter_small)) ALLOCATE(scatter_small(product(csize)))

  IF (.not.allocated(gather_big)) ALLOCATE(gather_big(product(csize)))
  IF (.not.allocated(gather_small)) ALLOCATE(gather_small(product(csize)))
  
  IF(.not.allocated(b_1))&
      ALLOCATE(b_1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)) 
  IF(.not.allocated(v_1))&
      ALLOCATE(v_1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)) 
  IF(.not.allocated(kxgrid)) ALLOCATE(kxgrid(1:nx0_big/2+1))
  IF(.not.allocated(kygrid)) ALLOCATE(kygrid(1:ny0_big))
  
  IF(spatial2d) THEN
    !kzgrid proportional to ky grid (see Watanabe and Sugama '04)
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(1:nz0_big))
  ELSE
    IF(.not.allocated(kzgrid)) ALLOCATE(kzgrid(1:nz0_big))
 END IF
 
  IF(.not.allocated(kmags)) ALLOCATE(kmags(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  IF(.not.allocated(kperps)) ALLOCATE(kperps(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  IF(.not.allocated(kzs)) ALLOCATE(kzs(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
 ! IF(.not.allocated(gpsi)) ALLOCATE(gpsi(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
 ! IF(.not.allocated(pre)) ALLOCATE(pre(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  IF(.not.allocated(pcurleig)) ALLOCATE(pcurleig(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  
  IF(kmin_eq_0) THEN

    kxgrid(1) = 0.0
    DO i=1,nx0_big/2
      kxgrid(i+1)=i*kxmin
    END DO
    kxmax=(nkx0-1)*kxmin

    kzgrid(1)=0.0
    DO i=1,nz0_big/2-1
      kzgrid(i+1)=i*kzmin
      kzgrid(nz0_big+1-i)=-i*kzmin
    END DO
    kzmax=(nkz0/2-1)*kzmin
    kzgrid(1+nz0_big/2)=nz0_big/2 * kzmin

    kygrid(1)=0.0
    DO i=1,ny0_big/2-1
      kygrid(i+1)=i*kymin
      kygrid(ny0_big+1-i)=-i*kymin
    END DO 
    kymax=(nky0/2-1)*kymin
    kygrid(1+ny0_big/2)=ny0_big/2 * kymin

  ELSE !.not.kmin_eq_0, i.e. linear
    kxgrid(1)=kxmin
    kxmax=kxmin
    kygrid(1)=kymin
    kymax=kymin
    kzgrid(1)=kzmin
    kzmax=kzmin
 END IF

 kmax = sqrt(kxmax**2.0 + kymax**2.0 + kzmax**2.0)

 DO i = cstart(1),cend(1)
    DO j = cstart(2),cend(2)
       DO k = cstart(3),cend(3)
          kmags(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2)
          kperps(i,j,k) = sqrt(kxgrid(i)**2 + kygrid(j)**2)
          kzs(i,j,k) = kzgrid(k)
       END DO
    END DO
 END DO

 DO i = cstart(1),cend(1)
    DO j = cstart(2),cend(2)
       DO k = cstart(3),cend(3)
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

 ALLOCATE(alpha_leftwhist(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
 ALLOCATE(alpha_leftcyclo(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

 alpha_leftwhist = -kmags/2.0 - sqrt(1.0 + ((kmags**2.0) / 4.0))
 alpha_leftcyclo = -kmags/2.0 + sqrt(1.0 + ((kmags**2.0) / 4.0))  
 
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
  
  IF(verbose.and.(mype.eq.0)) print *, 'Deallocated kgrids'

END SUBROUTINE finalize_arrays

