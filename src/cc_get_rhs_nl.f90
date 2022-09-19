!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                            cc_get_rhs_nl.f90                              !!
!!                                                                           !!
!!  nonlinearity                                                             !!
!!  -- initialize_fourier                                                    !!
!!  -- initialize_fourier_ae_nu0                                             !!
!!  -- get_rhs_nl                                                            !!
!!  -- get_rhs_nl1                                                           !!
!!  -- get_rhs_nl2                                                           !!
!!  -- get_rhs_nl3                                                           !!
!!  -- get_rhs_nl_convolution                                                !!
!!  -- get_k_indices                                                         !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               nonlinearity                                !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE nonlinearity
  USE mpi
  USE par_mod
  !USE hk_effects
  !USE flr_effects
  IMPLICIT NONE

  PUBLIC :: initialize_fourier,get_rhs_nl,&
            get_k_indices,get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  REAL, PUBLIC :: ve_max(2)

  PRIVATE

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0, bmag_in, bmagk
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,temp_big, temp_bigx, temp_bigy, temp_bigz
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: ekd

  REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  store_x, store_y, store_z
  REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  bmag, dxbmag, dybmag, dzbmag,bmag_inbig


  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz, dxbx, dybx,dzbx, dxby,dyby,dzby, dxbz,dybz,dzbz
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz, dxvx, dyvx,dzvx, dxvy,dyvy,dzvy, dxvz,dyvz,dzvz

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbx, dxdybx, dxdzbx, dydybx, dydzbx, dzdzbx
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxby, dxdyby, dxdzby, dydyby, dydzby, dzdzby
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbz, dxdybz, dxdzbz, dydybz, dydzbz, dzdzbz
    
  !For fft's

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:):: g_kbig
  REAL, ALLOCATABLE, DIMENSION(:,:,:):: g_rbig
  COMPLEX, ALLOCATABLE, DIMENSION(:,:):: g_kbig_2d
  REAL, ALLOCATABLE, DIMENSION(:,:):: g_rbig_2d
  INTEGER :: nx0_big,ny0_big,nz0_big
  INTEGER(kind=8) :: plan_r2c,plan_c2r
  INTEGER(kind=8) :: plan_kz2z,plan_z2kz
  INTEGER(kind=8) :: plan_ky2y,plan_y2ky
  INTEGER(kind=8) :: plan_kx2x,plan_x2kx
  INTEGER :: zpad
  REAL :: fft_norm  !normalization factor for inverse fft
 
  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             initialize_fourier                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialize_fourier

    CALL initialize_fourier_ae_mu0

END SUBROUTINE initialize_fourier


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             initialize_fourier_ae_mu0                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Note: only for adiabatic electrons (ae) and mu-integrated (mu0) version
SUBROUTINE initialize_fourier_ae_mu0

  include 'fftw3.f'
  
  !for dealiasing
  if (dealias_type.eq.1) zpad = 2
  if (dealias_type.eq.4) zpad = 4
  if (dealias_type.eq.3) zpad = 3
  
  nx0_big = zpad*nkx0
  ny0_big = zpad*nky0/2
  nz0_big = zpad*nkz0/2

  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))

  ALLOCATE(g_rbig(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(g_kbig(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))

  !WRITE(*,*) "making plans"
  CALL dfftw_plan_dft_c2r_3d(plan_c2r,nx0_big,ny0_big,nz0_big,&
                             g_kbig,g_rbig,FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_plan_dft_r2c_3d(plan_r2c,nx0_big,ny0_big,nz0_big,&
                             g_rbig,g_kbig,FFTW_FORWARD,FFTW_ESTIMATE)
  
  DEALLOCATE(g_rbig)
  DEALLOCATE(g_kbig)

  lky_big=ny0_big-hky_ind !Index of minimum (most negative) FILLED ky value for big arrays 
  lkz_big=nz0_big-hkz_ind !Index of minimum (most negative) FILLED kz value for big arrays 

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp(0,0,0),treal(0,0,0))
  !CALL dfftw_execute_dft_r2c(plan_r2c,treal(0,0,0),tcomp(0,0,0))
  !tcomp=tcomp/REAL(n1*n2*n3)

END SUBROUTINE initialize_fourier_ae_mu0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             initialize_fourier_ae_mu0_2d                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Note: only for adiabatic electrons (ae) and mu-integrated (mu0) version
!SUBROUTINE initialize_fourier_ae_mu0_2d

  ! deleted no need

!END SUBROUTINE initialize_fourier_ae_mu0_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             initialize_fourier2                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE initialize_fourier2
!
!  include 'fftw3.f'
!DELETED NO NEED
!
!END SUBROUTINE initialize_fourier2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   get_rhs_nl                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl(b_in, v_in, rhs_out_b, rhs_out_v,ndt)
  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt
  INTEGER :: delta_x,delta_y,delta_z
  
  !IF(mype==0) WRITE(*,*) "In get_rhs_nl"
  !IF(mype==0) WRITE(*,*) "Version is: ",rhs_nl_version
  IF ((rhs_nl_version==1).or.(rhs_nl_version == 12)) THEN
    !IF(mype==0) WRITE(*,*) "version was 1"
    IF (dealias_type.gt.2) CALL get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
    IF (dealias_type.eq.1) THEN
      rhs_out_b = 8.0 * rhs_out_b
      rhs_out_v = 8.0 * rhs_out_v
      DO delta_x = 0,1
        DO delta_y = 0,1
          DO delta_z = 0,1
            CALL get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt,&
              delta_x*pi/nkx0,delta_y*pi/nky0,delta_z*pi/nkz0)
          ENDDO
        ENDDO
      ENDDO
      rhs_out_b = rhs_out_b / 8.0
      rhs_out_v = rhs_out_v / 8.0
    ENDIF
!  ELSE IF(rhs_nl_version==2) THEN
!    CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==3) THEN
!    CALL get_rhs_nl3(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==4) THEN
!    CALL get_rhs_nl4(b_in,v_in,rhs_out_b,rhs_out_v)
  END IF
 
END SUBROUTINE get_rhs_nl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_rhs_nl2                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE get_rhs_nl2(g_in,phi_in,rhs_out)
!
!  DELETED NO NEED
!END SUBROUTINE get_rhs_nl2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl1                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt,dx,dy,dz)

  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt
  REAL, OPTIONAL :: dx,dy,dz
  INTEGER :: i,j,l,k,h,ierr

  IF (PRESENT(dz)) THEN
    ALLOCATE(ekd(0:nkx0-1,0:nky0-1,0:nkz0-1))
    DO i = 0,nx0_big/2
      DO j = 0,ny0_big-1
        DO k = 0,nz0_big-1
          ekd(i,j,k) = exp(i_complex * dx * i) * exp(i_complex * dy * j) * exp(i_complex * dz * k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(np_spec.gt.1) STOP "get_rhs_nl1 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl1 not yet implemented for np_kz.gt.1"
  IF(np_kz.ne.1) STOP "get_rhs_nl1 only suitable for np_kz=1"

  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigx(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigy(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigz(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  
  ALLOCATE(store_x(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_y(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_z(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

! All b arrays 
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

if (rhs_nl_version==1) then
    !! real space term
    ALLOCATE(bmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))  
    !! k space term
    ALLOCATE(bmagk(0:nkx0-1,0:nky0-1,lkz1:lkz2))  
    ALLOCATE(dxbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))  
    ALLOCATE(dybmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))  
    ALLOCATE(dzbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))  
endif

! All v arrays 
  ALLOCATE(vx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
! all first order v arrays
  !vx
  ALLOCATE(dxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !vy
  ALLOCATE(dxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !vz
  ALLOCATE(dxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
! all first order b arrays 
 !bx
  ALLOCATE(dxbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !by
  ALLOCATE(dxby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !bz
  ALLOCATE(dxbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  
! all  second order bx arrays  DXDXBX,   DXDYBX,   DXDZBX,  DYDYBX,   DYDZBX, DZDZBX
  ALLOCATE(dxdxbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order by arrays
  ALLOCATE(dxdxby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order bz arrays
  ALLOCATE(dxdxbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
    
  ALLOCATE(b_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))

  ! I dont want to change g_in, so I copy temporaly to g_in0
  !g_in0 = g_in
  b_inx0 = b_in(:,:,:,0)
  b_iny0 = b_in(:,:,:,1)
  b_inz0 = b_in(:,:,:,2)
  v_inx0 = v_in(:,:,:,0)
  v_iny0 = v_in(:,:,:,1)
  v_inz0 = v_in(:,:,:,2)
  !IF(mype==0) WRITE(*,*) "Actually in nl1"

! START SECOND ORDER b terms  

 ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
    !dxdxbx
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inx0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxbx(0,0,0))
   
   !dxdybx
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inx0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdybx(0,0,0))
    
     !dxdzbx
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inx0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzbx(0,0,0))
    
   ! DYDYbX
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inx0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydybx(0,0,0))
    
     ! DYDzbX
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inx0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzbx(0,0,0))

     ! DzDzbX
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inx0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzbx(0,0,0))
   
   ! FINISHED SECOND ORDER  BX TERMS 
   
! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY 

!dxdxby
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_iny0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxby(0,0,0))
   
   !dxdyvy
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_iny0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:nky0-1,lkz1:lkz2) = temp_small(i,:,:)
    ENDDO
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdyby(0,0,0))
    
     !dxdzby
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_iny0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzby(0,0,0))
   ! DYDYby
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_iny0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyby(0,0,0))
    
     ! DYDzbY
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_iny0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzby(0,0,0))
    
     ! DzDzbY
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_iny0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzby(0,0,0))
    
! FINISHED SECOND ORDER bY TERMS

! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
!dxdxbz
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inz0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxbz(0,0,0))
   
   !dxdyvz
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inz0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdybz(0,0,0))
    
     !dxdzbz
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inz0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzbz(0,0,0))

   ! DYDYbz
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inz0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydybz(0,0,0))
    
     ! DYDzbz
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inz0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzbz(0,0,0))
    
     ! DzDzbz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inz0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzbz(0,0,0))
    
!finished END SECOND ORDER BZ TERMS

!    completed  ALL SECOND ORDER B TERMS  

! TERMS BX BY BZ

    !bx
    DO i=0,nkx0-1
        temp_small(i,:,:)=b_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),bx(0,0,0))

    !by
    DO j=0,nky0-1
        temp_small(:,j,:)=b_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),by(0,0,0))

    !bz
    DO k=0,nkz0-1
        temp_small(:,:,k)=b_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),bz(0,0,0))

!!! bmag terms
if (rhs_nl_version == 1) then
!! Calculate b^2(x) in real space.
  DO i=0,nkx0-1
    DO j=0,nky0-1
        DO k=0,nkz0-1
            bmag(i,j,k) = bx(i,j,k)*bx(i,j,k)+ by(i,j,k)*by(i,j,k)+ bz(i,j,k)*bz(i,j,k)
        END DO
    END DO
  END DO
!!! Now forward transfrom to get b^2(k) stored in temp_big
temp_big=cmplx(0.0,0.0)
CALL dfftw_execute_dft_r2c(plan_r2c,bmag(0,0,0),temp_big(0,0,0))
!!! Now cut down to dealias and store result in bmagk which is b^2(k)
    
    call unpack(temp_big,bmagk)

!!! Now that we have bmagk, need to do IFFT(ikx b^2) IFFT(iky b^2) IFFT(ikz b^2) to get grad b^2/2 in real space
              
    !dxbmag                                                                                                                                                       
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*bmagk(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxbmag(0,0,0))

    !dybmag                                                                                                                                                                      
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*bmagk(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dybmag(0,0,0))

    !dzbmag                                                                                                                                                                      
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*bmagk(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzbmag(0,0,0))

    if (verbose) print *, 'Bmag gradients dealiased'
endif

!!! TERMS  vx,vy,vz 
!vx
    DO i=0,nkx0-1
        temp_small(i,:,:)=v_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vx(0,0,0))

    !vy
    DO j=0,nky0-1
                              !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
        temp_small(:,j,:)=v_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vy(0,0,0))

    !vz
    DO k=0,nkz0-1
        temp_small(:,:,k)=v_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vz(0,0,0))
    if (verbose) print *, 'V dealiased'

!  FIRST ORDER VX TERMS DXVX , DYVX,  DZVX

    ! dxvx
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvx(0,0,0))

    ! dyvx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inx0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvx(0,0,0))

    ! dzvx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inx0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvx(0,0,0))
    if (verbose) print *, 'dVx dealiased'
    
   !  FIRST ORDER VY TERMS  dxvy dyvy dzvz, 

    ! dxvy
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_iny0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvy(0,0,0))
 
    ! dyvy
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvy(0,0,0))

    ! dzvy
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_iny0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvy(0,0,0))
    if (verbose) print *, 'dVy dealiased'
    
    !  FIRST ORDER VZ TERMS  dxvz dyvz dzvz
    ! dxvz
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_inz0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvz(0,0,0))
 
    ! dyvz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvz(0,0,0))

    ! dzvz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvz(0,0,0))
    if (verbose) print *, 'dVz dealiased'
    
    ! DONE ALL FIRST ORDER VX,VY AND VZ TERMS i.e.  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
    
          ! FIRST ORDER BX TERMS ie. dxbx dybx dzbx` 
     ! dxbx
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxbx(0,0,0))

    ! dybx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_inx0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dybx(0,0,0))

    ! dzbx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inx0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzbx(0,0,0))
    if (verbose) print *, 'dBx dealiased'
        
    !  FIRST ORDER BY TERMS ie. dxby dyby dzby
     ! dxby
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_iny0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxby(0,0,0))
 
    ! dyby
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyby(0,0,0))

    ! dzby
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_iny0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzby(0,0,0))
    if (verbose) print *, 'dBy dealiased'
    
        !! FIRST ORDER BZ TERMS ie. dxbz dybz dzbz
! dxbz
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_inz0(i,:,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxbz(0,0,0))
 
    ! dybz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dybz(0,0,0))

    ! dzbz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    call zeropad(temp_small,temp_big)
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzbz(0,0,0))
    if (verbose) print *, 'dBz dealiased'

! DONE ALL FIRST ORDER BX,BY BZ TERMS ie. dxbx dybx dzbx,  dxby, dyby,dzby,  dxbz,dybz,dzbz

! EQUATION (14) 1=xcomp, 2=ycomp 3=zcomp
!     eq14x=P1-Q1-R1+S1
!     eq14y=P2-Q2-R2+S2
!     eq14z=P3-Q3-R3+S3
!P1 = bx*dxvx+by*dyvx+bz*dzvx
!Q1=vx*dxbx+vy*dybx+vz*dzbx
!R1=bx*dxdybz+by*dydybz+bz*dzdybz-bx*dxdzby-by*dydzby-bz*dzdzby
!S1=dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dz*bx 
!
!P2 = bx*dxvy+by*dyvy+bz*dzvy
!Q2=vx*dxby+vy*dyby+vz*dzby
!R2=bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dydxbz-bz*dzdxbz
!S2=dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dz*by 
!
!P3 = bx*dxvz+by*dyvz+bz*dzvz
!Q3=vx*dxbz+vy*dybz+vz*dzbz
!R3=bx*dxdxby+by*dydxby+bz*dzdxby-bx*dxdybx-by*dydybx-bz*dzdybx
!S3=dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dz*bz 

!Use store_i to store each product
!store_x = (bx*dxvx+by*dyvx+bz*dzvx) - (vx*dxbx+vy*dybx+vz*dzbx) - (bx*dxdybz+by*dydybz+bz*dzdybz-bx*dxdzby-by*dydzby-bz*dzdzby) + (dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dz*bx)
!
!store_y = (bx*dxvy+by*dyvy+bz*dzvy) - (vx*dxby+vy*dyby+vz*dzby) - (bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dydxbz-bz*dzdxbz) + (dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dz*by)
!
!store_z = (bx*dxvz+by*dyvz+bz*dzvz)- (vx*dxbz+vy*dybz+vz*dzbz)- (bx*dxdxby+by*dydxby+bz*dzdxby-bx*dxdybx-by*dydybx-bz*dzdybx)+ (dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dz*bz )
!
!Redo with correct order in terms
    store_x = (bx*dxvx+by*dyvx+bz*dzvx) - (vx*dxbx+vy*dybx+vz*dzbx)&
         - hall*((bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)&
         + (dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx))

    store_y = (bx*dxvy+by*dyvy+bz*dzvy) - (vx*dxby+vy*dyby+vz*dzby)&
         - hall*((bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)&
         + (dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby))

    store_z = (bx*dxvz+by*dyvz+bz*dzvz)- (vx*dxbz+vy*dybz+vz*dzbz)&
         - hall*((bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)&
         + (dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz ))

if (verbose) print *, 'B nls calculated'

IF (plot_nls) THEN 
! b.grad v 
CALL fft_write(bx*dxvx+by*dyvx+bz*dzvx,bdvio)
CALL fft_write(bx*dxvy+by*dyvy+bz*dzvy,bdvio)
CALL fft_write(bx*dxvz+by*dyvz+bz*dzvz,bdvio)

! v.grad b
CALL fft_write(vx*dxbx+vy*dybx+vz*dzbx,vdbio)
CALL fft_write(vx*dxby+vy*dyby+vz*dzby,vdbio)
CALL fft_write(vx*dxbz+vy*dybz+vz*dzbz,vdbio)

! b. grad curl b
CALL fft_write(bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby,bdcbio)
CALL fft_write(bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz,bdcbio)
CALL fft_write(bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx,bdcbio)

! curl b . grad b
CALL fft_write(dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx,cbdbio)
CALL fft_write(dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby,cbdbio)
CALL fft_write(dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz,cbdbio)

if (verbose) print *, 'B nl plot info written'

ENDIF

!inverse FFT to get back to Fourier
CALL dfftw_execute_dft_r2c(plan_r2c,store_x(0,0,0),temp_bigx(0,0,0))
CALL dfftw_execute_dft_r2c(plan_r2c,store_y(0,0,0),temp_bigy(0,0,0))
CALL dfftw_execute_dft_r2c(plan_r2c,store_z(0,0,0),temp_bigz(0,0,0))

!Now fill in appropriate rhs elements

CALL unpack(temp_bigx,temp_small)
rhs_out_b(:,:,:,0) = rhs_out_b(:,:,:,0) + temp_small
CALL unpack(temp_bigy,temp_small)
rhs_out_b(:,:,:,1) = rhs_out_b(:,:,:,1) + temp_small
CALL unpack(temp_bigz,temp_small)
rhs_out_b(:,:,:,2) = rhs_out_b(:,:,:,2) + temp_small

if (verbose) print *, 'B eqn stored'

!  DO i=0,nkx0-1
  !First x component
!            rhs_out_b(i,0:hky_ind,0:hkz_ind,0)=rhs_out_b(i,0:hky_ind,0:hkz_ind,0)+&           !kz positive, ky positive
!                                       temp_bigx(i,0:hky_ind,0:hkz_ind)*fft_norm
!            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,0)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,0)+& !kz negative, ky positive
!                                       temp_bigx(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)+& !kz negative, ky negative
!                                       temp_bigx(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,0)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,0)+& !kz positive, ky negative
!                                       temp_bigx(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !y component
!            rhs_out_b(i,0:hky_ind,0:hkz_ind,1)=rhs_out_b(i,0:hky_ind,0:hkz_ind,1)+&           !kz positive, ky positive
!                                       temp_bigy(i,0:hky_ind,0:hkz_ind)*fft_norm
!            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,1)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,1)+& !kz negative, ky positive
!                                       temp_bigy(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)+& !kz negative, ky negative
!                                       temp_bigy(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,1)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,1)+& !kz positive, ky negative
!                                       temp_bigy(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !z component
!            rhs_out_b(i,0:hky_ind,0:hkz_ind,2)=rhs_out_b(i,0:hky_ind,0:hkz_ind,2)+&           !kz positive, ky positive
!                                       temp_bigz(i,0:hky_ind,0:hkz_ind)*fft_norm
!            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,2)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,2)+& !kz negative, ky positive
!                                       temp_bigz(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)+& !kz negative, ky negative
!                                       temp_bigz(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,2)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,2)+& !kz positive, ky negative
!                                       temp_bigz(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
 ! END DO!i loop

!! ! EQUATION (15)
!! U1= vx*dxvx+vy*dyvx+vz*dzvx
!! V2 = bx*dxvx+by*dyvx+bz*dzvx
!! W1= (0.5)*dyb^2  ! TO BE DONE
!! 
!! U2= vx*dxvy+vy*dyvy+vz*dzvy
!! V2 = bx*dxvy+by*dyvy+bz*dzvy
!! W2= (0.5)*dyb^2  ! TO BE DONE
!! 
!! U3= vx*dxvz+vy*dyvz+vz*dzvz
!! V3 = bx*dxvz+by*dyvz+bz*dzvz
!! W3= (0.5)*dzb^2  ! TO BE DONE
!! 
!!       
!!     
!!      eq15x= -U1+V1 -W1
!!      eq15y=-U2+V2-W2
!!      eq15z=-U3+V3-W3

if (rhs_nl_version == 1) then
store_x = -(vx*dxvx+vy*dyvx+vz*dzvx) + (bx*dxbx+by*dybx+bz*dzbx) - 0.5*dxbmag
store_y = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+by*dyby+bz*dzby) - 0.5*dybmag
store_z = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz+bz*dzbz) - 0.5*dzbmag

IF (plot_nls) THEN 
! v . grad v
CALL fft_write(vx*dxvx+vy*dyvx+vz*dzvx,vdvio) 
CALL fft_write(vx*dxvy+vy*dyvy+vz*dzvy,vdvio)
CALL fft_write(vx*dxvz+vy*dyvz+vz*dzvz,vdvio)

! b . grad b
CALL fft_write(bx*dxbx+by*dybx+bz*dzbx,bdbio)
CALL fft_write(bx*dxby+by*dyby+bz*dzby,bdbio)
CALL fft_write(bx*dxbz+by*dybz+bz*dzbz,bdbio)

! 0.5 grad b^2
CALL  fft_write(0.5*dxbmag,db2io)
CALL  fft_write(0.5*dybmag,db2io)
CALL  fft_write(0.5*dzbmag,db2io)

if (verbose) print *, 'v1 nl equation stored'
endif
ENDIF

if (rhs_nl_version == 12) then
store_x = -(vx*dxvx+vy*dyvx+vz*dzvx) + (by*dybx+bz*dzbx) - (by*dxby+bz*dxbz)
store_y = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+bz*dzby) - (bz*dybz+bx*dybx)
store_z = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz) - (by*dzby+bx*dzbx)

IF (plot_nls) THEN 
! v . grad v
CALL fft_write(vx*dxvx+vy*dyvx+vz*dzvx,vdvio)
CALL fft_write(vx*dxvy+vy*dyvy+vz*dzvy,vdvio)
CALL fft_write(vx*dxvz+vy*dyvz+vz*dzvz,vdvio)

! b . grad b
CALL fft_write(bx*dxbx+by*dybx+bz*dzbx,bdbio)
CALL fft_write(bx*dxby+by*dyby+bz*dzby,bdbio)
CALL fft_write(bx*dxbz+by*dybz+bz*dzbz,bdbio)

! 0.5 grad b^2
CALL fft_write(bx*dxbx+by*dxby+bz*dxbz,db2io)
CALL fft_write(bx*dybx+by*dyby+bz*dybz,db2io)
CALL fft_write(bx*dzbx+by*dzby+bz*dzbz,db2io)

if (verbose) print *, 'v12 nl equation stored'
endif
ENDIF

CALL dfftw_execute_dft_r2c(plan_r2c,store_x(0,0,0),temp_bigx(0,0,0))
CALL dfftw_execute_dft_r2c(plan_r2c,store_y(0,0,0),temp_bigy(0,0,0))
CALL dfftw_execute_dft_r2c(plan_r2c,store_z(0,0,0),temp_bigz(0,0,0))

!Now fill in appropriate rhs elements

CALL unpack(temp_bigx,temp_small) 
rhs_out_v(:,:,:,0) = rhs_out_v(:,:,:,0) + temp_small
CALL unpack(temp_bigy,temp_small)
rhs_out_v(:,:,:,1) = rhs_out_v(:,:,:,1) + temp_small
CALL unpack(temp_bigz,temp_small)
rhs_out_v(:,:,:,2) = rhs_out_v(:,:,:,2) + temp_small

!  DO i=0,nkx0-1
  !First x component
!            rhs_out_v(i,0:hky_ind,0:hkz_ind,0)=rhs_out_v(i,0:hky_ind,0:hkz_ind,0)+&           !kz positive, ky positive
!                                       temp_bigx(i,0:hky_ind,0:hkz_ind)*fft_norm
!            rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,0)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,0)+& !kz negative, ky positive
!                                       temp_bigx(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)+& !kz negative, ky negative
!                                       temp_bigx(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
!            rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,0)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,0)+& !kz positive, ky negative
!                                       temp_bigx(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !y component
 !           rhs_out_v(i,0:hky_ind,0:hkz_ind,1)=rhs_out_v(i,0:hky_ind,0:hkz_ind,1)+&           !kz positive, ky positive
 !                                      temp_bigy(i,0:hky_ind,0:hkz_ind)*fft_norm
 !           rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,1)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,1)+& !kz negative, ky positive
 !                                      temp_bigy(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
 !           rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)+& !kz negative, ky negative
 !                                      temp_bigy(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
 !           rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,1)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,1)+& !kz positive, ky negative
 !                                      temp_bigy(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !z component
 !           rhs_out_v(i,0:hky_ind,0:hkz_ind,2)=rhs_out_v(i,0:hky_ind,0:hkz_ind,2)+&           !kz positive, ky positive
 !                                      temp_bigz(i,0:hky_ind,0:hkz_ind)*fft_norm
 !           rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,2)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,2)+& !kz negative, ky positive
 !                                      temp_bigz(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
 !           rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)+& !kz negative, ky negative
 !                                      temp_bigz(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
 !           rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,2)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,2)+& !kz positive, ky negative
 !                                      temp_bigz(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
 ! END DO!i loop

if (verbose) print *, 'rhs out v nl found'

if (calc_dt) CALL next_dt(ndt)
if (.not.(calc_dt)) ndt = dt_max
if (verbose) print *, 'next dt calculated ',ndt

IF (PRESENT(dz)) DEALLOCATE(ekd)

DEALLOCATE(temp_small)
if (verbose) print *, 'ts deallocated'
DEALLOCATE(temp_big)
if (verbose) print *, 'tb deallocated'

DEALLOCATE(temp_bigx)
if (verbose) print *, 'tbx deallocated'
DEALLOCATE(temp_bigy)
if (verbose) print *, 'tby deallocated'
DEALLOCATE(temp_bigz)
if (verbose) print *, 'tbz deallocated'

DEALLOCATE(store_x)
if (verbose) print *, 'stx deallocated'
DEALLOCATE(store_y)
if (verbose) print *, 'sty deallocated'
DEALLOCATE(store_z)
if (verbose) print *, 'stz deallocated'

! All b arrays 
DEALLOCATE(bx)
if (verbose) print *, 'bx deallocated'
DEALLOCATE(by)
if (verbose) print *, 'by deallocated'
DEALLOCATE(bz)
if (verbose) print *, 'first third deallocated'

if (rhs_nl_version==1) then
!DEALLOCATE(bmag_in) 
!if (verbose) print *,'bmagin de'
DEALLOCATE(bmagk)
if (verbose) print *,'bmagk deallocated'
DEALLOCATE(bmag)
if (verbose) print *,'bmag deallocated'
DEALLOCATE(dxbmag)
if (verbose) print *,'xbmag deallocated'
DEALLOCATE(dybmag)
if (verbose) print *,'ybmag deallocated'
DEALLOCATE(dzbmag)
if (verbose) print *,'ybmag deallocated'

!if (verbose) print *, 'bmag terms deallocated'
endif

! All v arrays 
  DEALLOCATE(vx)
  DEALLOCATE(vy)
  DEALLOCATE(vz)
! all first order v arrays
  !vx
  DEALLOCATE(dxvx)
  DEALLOCATE(dyvx)
  DEALLOCATE(dzvx)
  !vy
  DEALLOCATE(dxvy)
  DEALLOCATE(dyvy)
  DEALLOCATE(dzvy)
  !vz
  DEALLOCATE(dxvz)
  DEALLOCATE(dyvz)
  DEALLOCATE(dzvz)
! all first order b arrays 
 !bx
  DEALLOCATE(dxbx)
  DEALLOCATE(dybx)
  DEALLOCATE(dzbx)
  !by
  DEALLOCATE(dxby)
  DEALLOCATE(dyby)
  DEALLOCATE(dzby)
  !bz
  DEALLOCATE(dxbz)
  DEALLOCATE(dybz)
  DEALLOCATE(dzbz)
  
! all  second order bx arrays  DXDXBX,   DXDYBX,   DXDZBX,  DYDYBX,   DYDZBX, DZDZBX
  DEALLOCATE(dxdxbx)
  DEALLOCATE(dxdybx)
  DEALLOCATE(dxdzbx)
  DEALLOCATE(dydybx)
  DEALLOCATE(dydzbx)
  DEALLOCATE(dzdzbx)
  ! all  second o)
  DEALLOCATE(dxdxby)
  DEALLOCATE(dxdyby)
  DEALLOCATE(dxdzby)
  DEALLOCATE(dydyby)
  DEALLOCATE(dydzby)
  DEALLOCATE(dzdzby)
  ! all  second o)
  DEALLOCATE(dxdxbz)
  DEALLOCATE(dxdybz)
  DEALLOCATE(dxdzbz)
  DEALLOCATE(dydybz)
  DEALLOCATE(dydzbz)
  DEALLOCATE(dzdzbz)

  ! all  second o)

  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

if (verbose) print *, 'deallocated nl code'

END SUBROUTINE get_rhs_nl1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the fastest version
!Note: only for np_kz=1
!SUBROUTINE get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
!  DELETED NO NEED

!END SUBROUTINE get_rhs_nl2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl3                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE get_rhs_nl3(g_in,phi_in0,rhs_out)

   ! DELETED get_rhs_nl3

!END SUBROUTINE get_rhs_nl3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl4                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!For nkz0=1

!SUBROUTINE get_rhs_nl4(g_in,phi_in0,rhs_out)
  ! get_rhs_nl4
!END SUBROUTINE get_rhs_nl4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              get_rhs_nl_convolution                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Only for comparison with pseudospectral version.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE get_rhs_nl_convolution(g_in,phi_in0,rhs_out)
!    USE par_mod
!    IMPLICIT NONE
!  
!    COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!    COMPLEX, INTENT(in) :: phi_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!    COMPLEX :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!    COMPLEX, INTENT(inout) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!   
!    INTEGER :: i,j,k,l,ip,jp,kp
!    REAL :: ckxmax,ckymax,ckzmax
!    REAL :: kx,ky,kz,kxp,kyp,kzp,kxpp,kypp,kzpp
!    INTEGER :: xi,yi,zi,xip,yip,zip,xipp,yipp,zipp
!    LOGICAL :: take_conjgp,take_conjgpp
!    COMPLEX :: phi_kp,g_kpp
!  
!    IF(np_hank.gt.1) STOP "get_rhs_convolution not yet implemented for np_hank.gt.1"
!    IF(np_spec.gt.1) STOP "get_rhs_convolution not yet implemented for np_spec.gt.1"
!    IF(np_kz.gt.1) STOP "get_rhs_convolution not yet implemented for np_kz.gt.1"
!  
!    IF(np_kz.ne.1) STOP "get_rhs_convolution only suitable for np_kz=1"
!    DO k=0,nkz0-1
!      phi_in(:,:,k)=J0a(:,:)*phi_in0(:,:,k)
!    END DO
!  
!    ckxmax=REAL(hkx_ind)*kxmin
!    ckymax=REAL(hky_ind)*kymin
!    ckzmax=REAL(hkz_ind)*kzmin
!  
!  !  IF(mype==0) OPEN(unit=200,file='temp.dat',status='unknown')
!  
!    DO l=lv1,lv2
!      !WRITE(*,*) "l=",l
!      DO i=0,hkx_ind      !kx loop
!        DO j=-hky_ind,hky_ind    !ky loop
!          DO k=-hkz_ind,hkz_ind          !kz loop
!  
!             kx=REAL(i)*kxmin
!             ky=REAL(j)*kymin
!             kz=REAL(k)*kzmin
!             CALL get_k_indices(kx,ky,kz,xi,yi,zi,take_conjgp)
!             IF(take_conjgp) STOP "Error!"                   
!  
!             DO ip=-hkx_ind,hkx_ind   !kxp loop
!               DO jp=-hky_ind,hky_ind !kyp loop
!                 DO kp=-hkz_ind,hkz_ind     !kzp loop
!  
!                    kxp=REAL(ip)*kxmin 
!                    kyp=REAL(jp)*kymin 
!                    kzp=REAL(kp)*kzmin 
!  
!                    kxpp=kx-kxp
!                    kypp=ky-kyp
!                    kzpp=kz-kzp
!   
!                    !IF(i==2.and.j==0.and.k==0.and.ip==-13.and.jp==-1.and.kp==0.and.l==1) WRITE(*,*) "Before check.",&
!                    !           kxp,kyp,kzp,kxpp,kypp,kzpp,ckxmax,ckymax,ckzmax
!  
!                    IF((abs(kxpp).le.(ckxmax+0.001)).and.(abs(kypp).le.(ckymax+0.001)).and.(abs(kzpp).le.(ckzmax+0.001))) THEN
!  
!  
!                      CALL get_k_indices(kxp,kyp,kzp,xip,yip,zip,take_conjgp)
!                      IF(take_conjgp) THEN
!                        phi_kp=conjg(phi_in(xip,yip,zip))
!                      ELSE
!                        phi_kp=phi_in(xip,yip,zip)
!                      END IF
!  
!                      CALL get_k_indices(kxpp,kypp,kzpp,xipp,yipp,zipp,take_conjgpp)
!                      IF(take_conjgpp) THEN
!                        g_kpp=conjg(g_in(xipp,yipp,zipp,l,0,0))
!                      ELSE
!                        g_kpp=g_in(xipp,yipp,zipp,l,0,0)
!                      END IF
!          
!                      rhs_out(xi,yi,zi,l,0,0)=rhs_out(xi,yi,zi,l,0,0)+&
!                                  (kxp*ky-kx*kyp)*phi_kp*g_kpp
!  
!                    !IF(i==2.and.j==0.and.k==0.and.ip==-13.and.jp==-1.and.kp==0.and.l==1) WRITE(*,*) "After check.",&
!                    !              (kxp*ky-kx*kyp)*phi_kp*g_kpp
!  
!  !           IF(mype==0) THEN 
!  !             IF(i==2.and.j==0.and.k==0.and.(abs((kxp*ky-kx*kyp)*phi_kp*g_kpp).ne.0.0)) THEN
!  !                  WRITE(200,*) "i,j,k",i,j,k
!  !                  WRITE(200,*) "xip,yip,zip",xip,yip,zip
!  !                  WRITE(200,*) "xipp,yipp,zipp",xipp,yipp,zipp
!  !                  WRITE(200,*) "take_conjgp",take_conjgp
!  !                  WRITE(200,*) "take_conjgpp",take_conjgpp
!  !                  WRITE(200,*) "kxp,ky,kx,kyp",kxp,ky,kx,kyp
!  !                  WRITE(200,*) "C_k,kp",kxp*ky-kx*kyp
!  !                  WRITE(200,*) "phi_kp",phi_kp
!  !                  WRITE(200,*) "g_kpp",g_kpp
!  !                  WRITE(200,*) "(kxp*ky-kx*kyp)*phi_kp*g_kpp",(kxp*ky-kx*kyp)*phi_kp*g_kpp
!  !             END IF
!  !           END IF
!  
!                    END IF
!                    
!                 END DO
!               END DO
!             END DO
!  
!          END DO
!        END DO
!      END DO
!    END DO
!  
!  
!  !  IF(mype==0) CLOSE(200)
!  
!  END SUBROUTINE get_rhs_nl_convolution


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   get_k_indices                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This SUBROUTINE gets the correct indices for 
!!  a give kx,ky,kz (including negative).  IF kz is negative, take_conjg=T
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_k_indices(kx_in,ky_in,kz_in,i,j,k,take_conjg)
  USE par_mod
  IMPLICIT NONE

  REAL, INTENT(in) :: kx_in,ky_in,kz_in
  INTEGER, INTENT(out) :: i,j,k
  LOGICAL, INTENT(out) :: take_conjg
  REAL :: kx0,ky0,kz0

  take_conjg=.false.
  kx0=kx_in
  ky0=ky_in 
  kz0=kz_in 

  i=nint(kx0/kxmin)

  IF(kx0.lt.0.0) THEN
    take_conjg=.true.
    i=-1*i
    ky0=-1.0*ky0
    kz0=-1.0*kz0
  END IF

  IF(ky0.ge.0.0) THEN
    j=nint(ky0/kymin)
  ELSE
    j=nint(ky0/kymin+nky0)
  END IF

  IF(kz0.ge.0.0) THEN
    k=nint(kz0/kzmin)
  ELSE
    k=nint(kz0/kzmin+nkz0)
  END IF

END SUBROUTINE get_k_indices

SUBROUTINE next_dt(dtn)

 real, intent(out) :: dtn

 real :: ndt1xr,ndt1yr,ndt1zr,ndt2xr,ndt2yr,ndt2zr,ndt3xr,ndt3yr,ndt3zr
 real :: ndtr

 ndt1xr = maxval(abs(kxgrid))*maxval(abs(bx))
 ndt1yr = maxval(abs(kygrid))*maxval(abs(by))
 ndt1zr = maxval(abs(kzgrid))*maxval(1+abs(bz))
 ndt2xr = maxval(abs(kxgrid))*maxval(abs(vx))
 ndt2yr = maxval(abs(kygrid))*maxval(abs(vy))
 ndt2zr = maxval(abs(kzgrid))*maxval(abs(vz))
 ndt3xr = maxval(abs(kxgrid))*maxval(abs(dybz)+abs(dzby))
 ndt3yr = maxval(abs(kygrid))*maxval(abs(dzbx)+abs(dxbz))
 ndt3zr = maxval(abs(kzgrid))*maxval(abs(dxby)+abs(dybx))
 ndtr = ndt1xr + ndt1yr + ndt1zr &
   + ndt2xr + ndt2yr + ndt2zr &
   + ndt3xr + ndt3yr + ndt3zr
 dtn = courant/ndtr

END SUBROUTINE next_dt

SUBROUTINE fft_write(arr_real,unit)


real :: arr_real(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
integer :: unit
complex, allocatable, dimension(:,:,:) :: arr_spec
complex, allocatable, dimension(:,:,:) :: temporary
integer :: s

CALL dfftw_execute_dft_r2c(plan_r2c,arr_real(0,0,0),temporary(0,0,0))
CALL unpack(temporary,arr_spec)
WRITE(unit) arr_spec

END SUBROUTINE fft_write

SUBROUTINE zeropad(tempsmall,tempbig)


complex, allocatable, dimension(:,:,:) :: tempsmall
complex, allocatable, dimension(:,:,:) :: tempbig
integer :: i,j,k

tempbig = cmplx(0.0,0.0)
DO i = 0,nkx0-1
  IF ((zpad.eq.4).or.(zpad.eq.2)) THEN
    DO j = 0,nky0-1
      DO k = 0,nkz0-1
        IF (zpad.eq.4) tempbig(i,j,k) = tempsmall(i,j,k)
        IF (zpad.eq.2) tempbig(i,j,k) = tempsmall(i,j,k)*ekd(i,j,k)
      ENDDO
    ENDDO
  ENDIF
  IF(zpad.eq.3) THEN
      tempbig(i,0:hky_ind,0:hkz_ind)= tempsmall(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive
      tempbig(i,0:hky_ind,lkz_big:nz0_big-1)= tempsmall(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
      tempbig(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)= tempsmall(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative 
      tempbig(i,lky_big:ny0_big-1,0:hkz_ind)= tempsmall(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative  
  ENDIF
ENDDO
tempsmall = cmplx(0.0,0.0)

END SUBROUTINE zeropad

SUBROUTINE unpack(tempbig,tempsmall)


complex, allocatable, dimension(:,:,:) :: tempbig
complex, allocatable, dimension(:,:,:) :: tempsmall
integer :: i,j,k

tempsmall = cmplx(0.0,0.0)
DO i = 0,nkx0-1
  IF((zpad.eq.4).or.(zpad.eq.2)) THEN 
    DO j = 0,nky0-1
      DO k = 0,nkz0-1
        IF (zpad.eq.4) tempsmall(i,j,k) = fft_norm * tempbig(i,j,k)
        IF (zpad.eq.2) tempsmall(i,j,k) = fft_norm * tempbig(i,j,k) / ekd(i,j,k)
      ENDDO
    ENDDO
  ENDIF
  IF(zpad.eq.3) THEN 
      tempsmall(i,0:hky_ind,0:hkz_ind)= fft_norm * tempbig(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive
      tempsmall(i,0:hky_ind,lkz_ind:nkz0-1)= fft_norm * tempbig(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
      tempsmall(i,lky_ind:nky0-1,lkz_ind:nz0_big-1)= fft_norm * tempbig(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
      tempsmall(i,lky_ind:nky0-1,0:hkz_ind)= fft_norm * tempbig(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative  
  ENDIF
ENDDO
tempbig = cmplx(0.0,0.0)

END SUBROUTINE unpack

END MODULE nonlinearity

