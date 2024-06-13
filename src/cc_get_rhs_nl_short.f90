!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  REAL, PUBLIC :: ve_max(2)

  PRIVATE

  COMPLEX(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,temp_big,temp_big1

  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) ::  store

  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz,curlbx,curlby,curlbz
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz,curlvx,curlvy,curlvz

  !For fft's

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:):: g_kbig
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:):: g_rbig
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:):: g_kbig_2d
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:):: g_rbig_2d
  INTEGER :: nx0_big,ny0_big,nz0_big
  INTEGER(kind=8), allocatable :: plan_r2c,plan_c2r
  INTEGER(kind=8) :: plan_kz2z,plan_z2kz
  INTEGER(kind=8) :: plan_ky2y,plan_y2ky
  INTEGER(kind=8) :: plan_kx2x,plan_x2kx
  REAL :: fft_norm  !normalization factor for inverse fft
  INTEGER :: zpad
  INTEGER :: i,j,k
 
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
  if (((dealias_type.eq.3).or.(dealias_type.eq.4)).or.((dealias_type.eq.6).or.(dealias_type.eq.8))) zpad = dealias_type

  ny0_big=zpad*nky0/2
  if (splitx) then
  nx0_big = zpad*nkx0
  nz0_big = zpad*nkz0/2
  else
  nx0_big = zpad*nkx0/2
  nz0_big = zpad*nkz0
  endif

  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))
  ALLOCATE(plan_r2c)
  ALLOCATE(plan_c2r)

  ALLOCATE(store(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))

  !WRITE(*,*) "making plans"
  CALL dfftw_plan_dft_c2r_3d(plan_c2r,nx0_big,ny0_big,nz0_big,&
                             temp_big,store,FFTW_ESTIMATE,FFTW_BACKWARD)
  CALL dfftw_plan_dft_r2c_3d(plan_r2c,nx0_big,ny0_big,nz0_big,&
       store,temp_big,FFTW_ESTIMATE,FFTW_FORWARD)

  lky_big=ny0_big-hky_ind !Index of minimum (most negative) FILLED ky value for big arrays
  if (.not.splitx) lkx_big = nx0_big-hkx_ind
  if (splitx) lkz_big=nz0_big-hkz_ind !Index of minimum (most negative) FILLED kz value for big arrays 

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  CALL ALLOCATIONS

  IF (verbose) CALL DEALIASINGTEST
  if (verbose) CALL PRECISIONTEST
  rkstage = 0

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp,treal)
  !CALL dfftw_execute_dft_r2c(plan_r2c,treal,tcomp)
  !tcomp=tcomp/REAL(n1*n2*n3)

END SUBROUTINE initialize_fourier_ae_mu0

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
  
  !IF(mype==0) WRITE(*,*) "In get_rhs_nl"
  !IF(mype==0) WRITE(*,*) "Version is: ",rhs_nl_version
  IF ((rhs_nl_version.eq.1).or.(rhs_nl_version.eq.12)) THEN
    !IF(mype==0) WRITE(*,*) "version was 1"
    CALL get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
!  ELSE IF(rhs_nl_version==2) THEN
!    CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==3) THEN
!    CALL get_rhs_nl3(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==4) THEN
!    CALL get_rhs_nl4(b_in,v_in,rhs_out_b,rhs_out_v)
  END IF
  rkstage = rkstage + 1
 
END SUBROUTINE get_rhs_nl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl1                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt)

  USE par_mod
  include 'fftw3.f'

  
  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt

                                                                !COMPLEX :: temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1)
                                                                !COMPLEX :: temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1)
                                                                !REAL :: dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                                !REAL :: dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                                !REAL :: dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                                !REAL :: dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
   INTEGER :: l,h, ierr

  IF(np_spec.gt.1) STOP "get_rhs_nl1 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl1 not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_nl1 only suitable for np_kz=1"

  ! I dont want to change g_in, so I copy temporaly to g_in0
  !g_in0 = g_in
  b_inx0 = b_in(:,:,:,0)
  b_iny0 = b_in(:,:,:,1)
  b_inz0 = b_in(:,:,:,2)
  v_inx0 = v_in(:,:,:,0)
  v_iny0 = v_in(:,:,:,1)
  v_inz0 = v_in(:,:,:,2)
  !IF(mype==0) WRITE(*,*) "Actually in nl1"

! TERMS BX BY BZ


  if (.not.nv) then ! Skip b if Navier Stokes
    !bx
    DO i=0,nkx0-1
        temp_small(i,:,:)=b_inx0(i,:,:)
    END DO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    bx = store

    !by
    DO j=0,nky0-1
        temp_small(:,j,:)=b_iny0(:,j,:)
    END DO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    by = store

    !bz
    DO k=0,nkz0-1
        temp_small(:,:,k)=b_inz0(:,:,k)
    END DO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    bz = store

    ! curlbx
    DO j = 0,nky0-1
       DO k = 0,nkz0-1
          temp_small(:,j,k) = i_complex * kygrid(j) * b_inz0(:,j,k) - i_complex  * kzgrid(k) * b_iny0(:,j,k)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlbx = store
    
    ! curlby
    DO k = 0,nkz0-1
       DO i = 0,nkx0-1
          temp_small(i,:,k) = i_complex * kzgrid(k) * b_inx0(i,:,k) - i_complex * kxgrid(i) * b_inz0(i,:,k)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlby = store

    ! curlbz
    DO i = 0,nkx0-1
       DO j = 0,nky0-1
          temp_small(i,j,:) = i_complex * kxgrid(i) * b_iny0(i,j,:) - i_complex * kygrid(j) * b_inx0(i,j,:)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlbz = store
    
endif !Skip b FFTs if Navier Stokes 

!!! TERMS  vx,vy,vz 
!vx
    DO i=0,nkx0-1
        temp_small(i,:,:)=v_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    vx = store

    !vy
    DO j=0,nky0-1
                              !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
        temp_small(:,j,:)=v_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    vy = store

    !vz
    DO k=0,nkz0-1
        temp_small(:,:,k)=v_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    vz = store

    ! curlvx
    DO j = 0,nky0-1
       DO k = 0,nkz0-1
          temp_small(:,j,k) = i_complex * kygrid(j) * v_inz0(:,j,k) - i_complex  * kzgrid(k) * v_iny0(:,j,k)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlvx = store

    ! curlvy
    DO k = 0,nkz0-1
       DO i = 0,nkx0-1
          temp_small(i,:,k) = i_complex * kzgrid(k) * v_inx0(i,:,k) - i_complex * kxgrid(i) * v_inz0(i,:,k)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlvy = store

    ! curlvz
    DO i = 0,nkx0-1
       DO j = 0,nky0-1
          temp_small(i,j,:) = i_complex * kxgrid(i) * v_iny0(i,j,:) - i_complex * kygrid(j) * v_inx0(i,j,:)
       ENDDO
    ENDDO
    CALL ZEROPAD
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
    curlvz = store
    
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
!store_y = (bx*dxvy+by*dyvy+bz*dzvy) - (vx*dxby+vy*dyby+vz*dzby) - (bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dydxbz-bz*dzdxbz) + (dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dz*by)
!store_z = (bx*dxvz+by*dyvz+bz*dzvz)- (vx*dxbz+vy*dybz+vz*dzbz)- (bx*dxdxby+by*dydxby+bz*dzdxby-bx*dxdybx-by*dydybx-bz*dzdybx)+ (dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dz*bz )

    ! The RHS of Mahajan 2021 equation 14 is b dot grad v - v dot grad b - b dot grad curl b + curl b dot grad b
    ! = curl (v x b - curl b x b)
    ! We dodge all second order b derivatives and most b derivatives by only finding F = v x b - curl b x b pseudospectrally
    ! and then finding curl F in Fourier space

    ! To save memory, only use one store by adding terms to rhs_out as needed

    if (.not.nv) then ! Skip b ffts if Navier Stokes

       ! x: vy bz - vz by - (curlby bz - curlbz by)
       store = vy * bz - vz * by - hall * (curlby * bz - curlbz * by)
       CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
       CALL UNPACK

       if ((verbose).and.mype.eq.0) print *, "MaxNLBx",&
            maxval(abs(temp_small(1:nkx0-1,1:nky0-1,1:nkz0-1))/abs(rhs_out_b(1:nkx0-1,1:nky0-1,1:nkz0-1,0)),&
            ((abs(rhs_out_b(1:nkx0-1,1:nky0-1,1:nkz0-1,0)).gt.10.0**(-8.0)))),&
    (/1,1,1/) + maxloc(abs(temp_small(1:nkx0-1,1:nky0-1,1:nkz0-1))/abs(rhs_out_b(1:nkx0-1,1:nky0-1,1:nkz0-1,0)),&
            ((abs(rhs_out_b(1:nkx0-1,1:nky0-1,1:nkz0-1,0)).gt.10.0**(-8.0))))

       DO k = 0,nkz0-1
          rhs_out_b(:,:,k,1) = rhs_out_b(:,:,k,1) + i_complex * kzgrid(k) * temp_small(0:nkx0-1,:,k)
       ENDDO
       DO j = 0,nky0-1
          rhs_out_b(:,j,:,2) = rhs_out_b(:,j,:,2) - i_complex * kygrid(j) * temp_small(0:nkx0-1,j,:)
       ENDDO

       ! y: vz bx - vx bz - hall * (curlbz bx - curlbx bz)
       store = vz * bx - vx * bz - hall * (curlbz * bx - curlbx * bz)
       CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
       CALL UNPACK

       DO i = 0,nkx0-1
          rhs_out_b(i,:,:,2) = rhs_out_b(i,:,:,2) + i_complex * kxgrid(i) * temp_small(i,:,:)
       ENDDO
       DO k = 0,nkz0-1
          rhs_out_b(:,:,k,0) = rhs_out_b(:,:,k,0) - i_complex * kzgrid(k) * temp_small(0:nkx0-1,:,k)
       ENDDO

       ! z: vx by - vy bx - hall (curlbx by - curlby bx)
       store = vx * by - vy * bx - hall*(curlbx * by - curlby * bx)
       CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
       CALL UNPACK

       DO j = 0,nky0-1
          rhs_out_b(:,j,:,0) = rhs_out_b(:,j,:,0) + i_complex * kygrid(j) * temp_small(0:nkx0-1,j,:)
       ENDDO
       DO i = 0,nkx0-1
          rhs_out_b(i,:,:,1) = rhs_out_b(i,:,:,1) - i_complex * kxgrid(i) * temp_small(i,:,:)
       ENDDO

endif ! Skip b if Navier Stokes
  
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

! Above equations are - v dot grad v + b dot grad b - 1/2 grad b^2
! Equal to v x curl v - 1/2 grad v^2 + curl b x b
! We're incompressible, so the grad term will be projected out with pressure
! (If compressibility gets included, then would handle v^2 pseudospectrally with three more FFTs (or rewrite whole as grad))

! x: vy curlvz - vz curlvy + curlby bz - curlbz by

store = vy * curlvz - vz * curlvy + curlby * bz - curlbz * by
CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
CALL UNPACK

       if ((verbose).and.mype.eq.0) print *, "MaxNLVx",&
            maxval(abs(temp_small(1:nkx0-1,1:nky0-1,1:nkz0-1))/abs(rhs_out_v(1:nkx0-1,1:nky0-1,1:nkz0-1,0)),&
            ((abs(rhs_out_v(1:nkx0-1,1:nky0-1,1:nkz0-1,0)).gt.10.0**(-8.0)))),&
    (/1,1,1/) + maxloc(abs(temp_small(1:nkx0-1,1:nky0-1,1:nkz0-1))/abs(rhs_out_v(1:nkx0-1,1:nky0-1,1:nkz0-1,0)),&
            ((abs(rhs_out_v(1:nkx0-1,1:nky0-1,1:nkz0-1,0)).gt.10.0**(-8.0))))
       
rhs_out_v(:,:,:,0) = rhs_out_v(:,:,:,0) + temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1)

! y: vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz
store = vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz
CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
CALL UNPACK
if (verbose) print *, "vy fft complete"
rhs_out_v(:,:,:,1) = rhs_out_v(:,:,:,1) + temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1)
if (verbose) print *, "vy done"

! z: vx * curlvy - vy * curlvx + curlbx * by - curlby * bx
store = vx * curlvy - vy * curlvx + curlbx * by - curlby * bx
CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
CALL UNPACK
if (verbose) print *, "vz fft complete"
rhs_out_v(:,:,:,2) = rhs_out_v(:,:,:,2) + temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1)
  if (verbose) print *, "vz done"

  ! to preserve reality of the fields, remove v,b terms at nky0/2,nkz0/2
  rhs_out_b(:,nky0/2,:,:) = 0
  rhs_out_v(:,nky0/2,:,:) = 0
  rhs_out_b(:,:,nkz0/2,:) = 0
  rhs_out_v(:,:,nkz0/2,:) = 0
      
if (verbose) print *, 'rhs out v nl found'

if (calc_dt) CALL next_dt(ndt)
if (.not.(calc_dt)) ndt = dt_max

if (verbose) print *, 'next dt calculated ',ndt

END SUBROUTINE get_rhs_nl1

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
 ndt3xr = maxval(abs(kxgrid))*maxval(abs(curlbx))
 ndt3yr = maxval(abs(kygrid))*maxval(abs(curlby))
 ndt3zr = maxval(abs(kzgrid))*maxval(abs(curlbz))
 ndtr = ndt1xr + ndt1yr + ndt1zr &
   + ndt2xr + ndt2yr + ndt2zr &
   + ndt3xr + ndt3yr + ndt3zr
 dtn = courant/ndtr

END SUBROUTINE next_dt

SUBROUTINE finalize_fourier

implicit none

CALL DEALLOCATIONS
call dfftw_destroy_plan(plan_r2c)
call dfftw_destroy_plan(plan_c2r)
call dfftw_cleanup()

DEALLOCATE(store)
DEALLOCATE(temp_big)

DEALLOCATE(plan_r2c)
DEALLOCATE(plan_c2r)

END SUBROUTINE finalize_fourier

SUBROUTINE ALLOCATIONS

  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))

! All b arrays 
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

! All v arrays 
  ALLOCATE(vx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all first order v arrays

  ALLOCATE(curlvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(curlvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(curlvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  
! all first order b arrays 

  ALLOCATE(curlbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(curlby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(curlbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  
  ALLOCATE(b_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))

END SUBROUTINE ALLOCATIONS

SUBROUTINE DEALLOCATIONS

DEALLOCATE(temp_small)
if (verbose) print *, 'ts deallocated'

! All b arrays 
DEALLOCATE(bx)
if (verbose) print *, 'bx deallocated'
DEALLOCATE(by)
if (verbose) print *, 'by deallocated'
DEALLOCATE(bz)
if (verbose) print *, 'first third deallocated'

DEALLOCATE(curlbx)
DEALLOCATE(curlby)
DEALLOCATE(curlbz)

! All v arrays 
  DEALLOCATE(vx)
  DEALLOCATE(vy)
  DEALLOCATE(vz)

  DEALLOCATE(curlvx)
  DEALLOCATE(curlvy)
  DEALLOCATE(curlvz)
  
  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

END SUBROUTINE DEALLOCATIONS

SUBROUTINE DEALIASINGTEST

  IMPLICIT NONE

  temp_small = 0.0
  temp_small(1,1,1) = 1

  CALL ZEROPAD
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
  vx = store  
  
  temp_small = 0.0
  DO i = 0,nkx0-1
     temp_small(i,:,:) = i
  ENDDO
  temp_small(:,nky0/2,:) = 0
  temp_small(:,:,nkz0/2) = 0
  temp_small(:,0,:) = 0
  temp_small(:,:,0) = 0

  CALL ZEROPAD
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
  vy = store
  
  temp_small = 0.0
  temp_small(0,0,0) = 1
  CALL ZEROPAD
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
  vz = store
  
  store = vx * vy
  CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
  CALL UNPACK
  
  print *, "Dealiasing Test Result",temp_small(1,1,1)

  store = vx * vz
  CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)
  CALL UNPACK
  
  print *, "Control Result", temp_small(1,1,1)

  temp_small = cmplx(0.0,0.0)
  temp_big = cmplx(0.0,0.0)
  store = cmplx(0.0,0.0)

END SUBROUTINE DEALIASINGTEST

SUBROUTINE ZEROPAD

  IMPLICIT NONE

  integer :: i 
  
    temp_big=cmplx(0.0,0.0)

    DO i=0,nkx0-1
       temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive
       temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
       temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
       temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
     END DO!k loop                             

END SUBROUTINE ZEROPAD

SUBROUTINE UNPACK

  IMPLICIT NONE

  integer :: i

  temp_small = cmplx(0.0,0.0)
  DO i = 0,nkx0-1
     temp_small(i,0:hky_ind,0:hkz_ind) = temp_big(i,0:hky_ind,0:hkz_ind)*fft_norm
     temp_small(i,0:hky_ind,lkz_ind:nkz0-1) = temp_big(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
     temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) = temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
     temp_small(i,lky_ind:nky0-1,0:hkz_ind) = temp_big(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  ENDDO

END SUBROUTINE UNPACK

SUBROUTINE PRECISIONTEST

  IMPLICIT NONE

  ! Estimate the precision of the transforms by doing a transform and then an inverse transform

  ALLOCATE(temp_big1(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))

  temp_big1 = 1.0
  temp_big = temp_big1

  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big,store)
  CALL dfftw_execute_dft_r2c(plan_r2c,store,temp_big)

  print *,"Precision Test",maxval(abs(temp_big1-temp_big*fft_norm))

  DEALLOCATE(temp_big1)

END SUBROUTINE PRECISIONTEST


END MODULE nonlinearity

