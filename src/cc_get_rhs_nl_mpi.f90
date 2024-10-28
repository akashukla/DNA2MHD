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
  
  use, intrinsic :: iso_c_binding
  
  IMPLICIT NONE

  include 'fftw3-mpi.f03'

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  REAL, PUBLIC :: ve_max(2)

  PRIVATE

  COMPLEX(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,bigmask,smallmask

  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz,curlbx,curlby,curlbz
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz,curlvx,curlvy,curlvz

  !For fft's

  INTEGER(C_INTPTR_T) :: nx0_big,ny0_big,nz0_big
  INTEGER(C_INTPTR_T) :: local_N,local_k_offset,alloc_local
  
  type(C_PTR) :: plan_r2c,plan_c2r
  type(C_PTR) :: rinout,cinout

  COMPLEX(C_DOUBLE_COMPLEX), pointer :: temp_big(:,:,:)
  REAL(C_DOUBLE), pointer ::  store(:,:,:)
  
  REAL :: fft_norm  !normalization factor for inverse fft
  INTEGER(C_INTPTR_T) :: zpad
  INTEGER(C_INTPTR_T) :: i,j,k
  INTEGER :: ierr

  LOGICAL :: padv1
  LOGICAL :: padv2
  LOGICAL :: printpad,printunpack
  
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
  
  !for dealiasing

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL fftw_mpi_init()
  
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

  alloc_local = fftw_mpi_local_size_3d(nz0_big,ny0_big,nx0_big/2+1,MPI_COMM_WORLD,local_N,local_k_offset)
  rinout = fftw_alloc_real(2*alloc_local)
  cinout = fftw_alloc_complex(alloc_local)
  
  CALL c_f_pointer(rinout,store,[2*(nx0_big/2+1),ny0_big,local_N])
  CALL c_f_pointer(cinout,temp_big,[nx0_big/2+1,ny0_big,local_N])
  
  plan_c2r = fftw_mpi_plan_dft_c2r_3d(nz0_big,ny0_big,nx0_big,temp_big,store,MPI_COMM_WORLD,FFTW_PATIENT)
  plan_r2c = fftw_mpi_plan_dft_r2c_3d(nz0_big,ny0_big,nx0_big,store,temp_big,MPI_COMM_WORLD,FFTW_PATIENT)
     
  !WRITE(*,*) "making plans"

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
  WRITE(*,*) "local_N",local_N,mype
  WRITE(*,*) "local_k_offset",local_k_offset,mype

  CALL ALLOCATIONS

  padv1 = (.false.)
  padv2 = (.true.)
  printpad = (.true.)
  printunpack = (.true.)

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp,treal)
  !CALL dfftw_execute_dft_r2c(plan_r2c,treal,tcomp)
  
  !tcomp=tcomp/REAL(n1*n2*n3)

END SUBROUTINE initialize_fourier_ae_mu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   get_rhs_nl                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl(b_in, v_in, rhs_out_b, rhs_out_v,ndt)
  USE par_mod

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
 
END SUBROUTINE get_rhs_nl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl1                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt)

  USE par_mod
  
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
   temp_small = b_inx0
   CALL ZEROPAD
   bx = store

   !by
   temp_small = b_iny0
   CALL ZEROPAD    
   by = store
   
   !bz
   temp_small = b_inz0
   CALL ZEROPAD
   bz = store

    ! curlbx
    DO j = 0,nky0-1
       DO k = lkz1,lkz2
          temp_small(:,j,k) = i_complex * kygrid(j) * b_inz0(:,j,k) - i_complex  * kzgrid(k) * b_iny0(:,j,k)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlbx = store
    
    ! curlby
    DO k = lkz1,lkz2
       DO i = 0,nkx0-1
          temp_small(i,:,k) = i_complex * kzgrid(k) * b_inx0(i,:,k) - i_complex * kxgrid(i) * b_inz0(i,:,k)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlby = store

    ! curlbz
    DO i = 0,nkx0-1
       DO j = 0,nky0-1
          temp_small(i,j,:) = i_complex * kxgrid(i) * b_iny0(i,j,:) - i_complex * kygrid(j) * b_inx0(i,j,:)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlbz = store
    
 endif !Skip b FFTs if Navier Stokes 

!!! TERMS  vx,vy,vz 
!vx
    temp_small = v_inx0
    !Add padding for dealiasing
    CALL ZEROPAD
    vx = store

    !vy
    temp_small = v_iny0
    !Add padding for dealiasing
    CALL ZEROPAD    
    vy = store

    !vz
    temp_small = v_inz0
    !Add padding for dealiasing    
    CALL ZEROPAD    
    vz = store

    ! curlvx
    DO j = 0,nky0-1
       DO k = lkz1,lkz2
          temp_small(:,j,k) = i_complex * kygrid(j) * v_inz0(:,j,k) - i_complex  * kzgrid(k) * v_iny0(:,j,k)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlvx = store

    ! curlvy
    DO k = lkz1,lkz2
       DO i = 0,nkx0-1
          temp_small(i,:,k) = i_complex * kzgrid(k) * v_inx0(i,:,k) - i_complex * kxgrid(i) * v_inz0(i,:,k)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlvy = store

    ! curlvz
    DO i = 0,nkx0-1
       DO j = 0,nky0-1
          temp_small(i,j,:) = i_complex * kxgrid(i) * v_iny0(i,j,:) - i_complex * kygrid(j) * v_inx0(i,j,:)
       ENDDO
    ENDDO
    CALL ZEROPAD    
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
       
       CALL UNPACK

       DO k = lkz1,lkz2
          rhs_out_b(:,:,k,1) = rhs_out_b(:,:,k,1) + i_complex * kzgrid(k) * temp_small(:,:,k)
       ENDDO
       DO j = 0,nky0-1
          rhs_out_b(:,j,:,2) = rhs_out_b(:,j,:,2) - i_complex * kygrid(j) * temp_small(:,j,:)
       ENDDO

       ! y: vz bx - vx bz - hall * (curlbz bx - curlbx bz)
       store = vz * bx - vx * bz - hall * (curlbz * bx - curlbx * bz)
       
       CALL UNPACK

       DO i = 0,nkx0-1
          rhs_out_b(i,:,:,2) = rhs_out_b(i,:,:,2) + i_complex * kxgrid(i) * temp_small(i,:,:)
       ENDDO
       DO k = lkz1,lkz2
          rhs_out_b(:,:,k,0) = rhs_out_b(:,:,k,0) - i_complex * kzgrid(k) * temp_small(:,:,k)
       ENDDO

       ! z: vx by - vy bx - hall (curlbx by - curlby bx)
       store = vx * by - vy * bx - hall*(curlbx * by - curlby * bx)
       
       CALL UNPACK

       DO j = 0,nky0-1
          rhs_out_b(:,j,:,0) = rhs_out_b(:,j,:,0) + i_complex * kygrid(j) * temp_small(:,j,:)
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

CALL UNPACK

rhs_out_v(:,:,:,0) = rhs_out_v(:,:,:,0) + temp_small

! y: vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz
store = vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz

CALL UNPACK
if (verbose) print *, "vy fft complete"
rhs_out_v(:,:,:,1) = rhs_out_v(:,:,:,1) + temp_small
if (verbose) print *, "vy done"

! z: vx * curlvy - vy * curlvx + curlbx * by - curlby * bx
store = vx * curlvy - vy * curlvx + curlbx * by - curlby * bx

CALL UNPACK
if (verbose) print *, "vz fft complete"
rhs_out_v(:,:,:,2) = rhs_out_v(:,:,:,2) + temp_small
if (verbose) print *, "vz done"

! to preserve reality of the fields, remove v,b terms at nky0/2,nkz0/2
rhs_out_b(:,nky0/2,:,:) = 0
rhs_out_v(:,nky0/2,:,:) = 0

! Condition to only do this for when nkz0/2 in mype needed
if ((nkz0/2 .ge. lkz1).and.(nkz0/2 .le. lkz2)) then
   rhs_out_b(:,:,nkz0/2,:) = 0
   rhs_out_v(:,:,nkz0/2,:) = 0
endif
  
if (verbose) print *, 'rhs out v nl found'

if (calc_dt) CALL next_dt(ndt)
if (.not.(calc_dt)) ndt = dt_max

if (calc_dt) print *, 'next dt calculated ',ndt

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

call fftw_destroy_plan(plan_r2c)
call fftw_destroy_plan(plan_c2r)
 
! if (verbose) print *, "plans destroyed"

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL fftw_free(cinout)
if (verbose) print *, "temp_big freed"
CALL fftw_free(rinout)

if (verbose) print *, "store freed"

CALL fftw_mpi_cleanup()

END SUBROUTINE finalize_fourier

SUBROUTINE ALLOCATIONS

  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  

  ! All b arrays
  ALLOCATE(bx(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(by(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(bz(1:nx0_big+2,1:ny0_big,1:local_N))
  
  ! All v arrays
  ALLOCATE(vx(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(vy(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(vz(1:nx0_big+2,1:ny0_big,1:local_N))

  ! all first order v arrays
  ALLOCATE(curlvx(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(curlvy(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(curlvz(1:nx0_big+2,1:ny0_big,1:local_N))
  
  ! all first order b arrays
  ALLOCATE(curlbx(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(curlby(1:nx0_big+2,1:ny0_big,1:local_N))
  ALLOCATE(curlbz(1:nx0_big+2,1:ny0_big,1:local_N))

  ALLOCATE(bigmask(1:nx0_big/2+1,1:ny0_big,1:local_N))
  ALLOCATE(smallmask(0:nkx0-1,0:nky0-1,0:nkz0/n_mpi_procs-1))
  
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

  if (verbose) print *, "all derivatives deallocated"

  DEALLOCATE(bigmask)

  if (verbose) print *, "big mask deallocated"
  DEALLOCATE(smallmask)

  if (verbose) print *, "small mask deallocated"
  
  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

END SUBROUTINE DEALLOCATIONS

SUBROUTINE ZEROPAD

  IMPLICIT NONE

  integer :: lkz1_rank,lkz2_rank
  integer :: rank,k,kp,i
  integer :: ind
  logical :: zmask

  temp_big=cmplx(0.0,0.0)
  if (verbose) print *, "Entering Zeropad"

  if (padv1) then
  
  DO i=0,nkx0-1
     temp_big(i+1,1:1+hky_ind,1:hkz_ind+1)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive
!     !if (verbose) print *, i,"++ done"
     temp_big(i+1,1:1+hky_ind,1+lkz_big:nz0_big)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive             
!     ! if (verbose) print *, i,"+- done"
     temp_big(i+1,1+lky_big:ny0_big,1+lkz_big:nz0_big)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative  
!     ! if (verbose) print *, i,"-- done"
     temp_big(i+1,1+lky_big:ny0_big,1:hkz_ind+1)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative             
!     ! if (verbose) print *, i,"-+ done"
  END DO

endif

if (padv2) then

  DO rank = 0,n_mpi_procs-1
     
     lkz1_rank = rank*(3*nkz0/2)/n_mpi_procs
     lkz2_rank = (rank+1)*(3*nkz0/2)/n_mpi_procs-1

     ! print *, "MaxVal temp_small",maxval(abs(temp_small))

     DO k = lkz1,lkz2
     
        ! Shift up values of k greater than nkz0/2 - match the larger array
        if (k.gt.nkz0/2) kp = k + nkz0/2
        if (k.le.nkz0/2) kp = k

        ! Mask away values that either are outside the rank bounds
        zmask = ((kp.ge.lkz1_rank).and.(kp.le.lkz2_rank))
        ! Or that are taken out because of dealiasing
        zmask = (zmask.and.((kp.ge.lkz_big).or.(kp.le.hkz_ind)))

        ! Find position of ind in bigmask - lkz1_rank is min in rank, kp has been shifted up
        ! If the minimum below is lkz2_rank and not both equal, then zmask will be zero so works
        ! But also make sure that the index doesn't fall below one - if it would have, zmask takes out
        ind = max(min(1+kp - lkz1_rank,1+lkz2_rank),1)
        
        if (zmask) bigmask(1:nkx0,1:hky_ind+1,ind) = temp_small(0:nkx0-1,0:hky_ind,k)
        if (zmask) bigmask(1:nkx0,lky_big+1:ny0_big,ind) = temp_small(0:nkx0-1,lky_ind:nky0-1,k)

     ENDDO

     ! print *, "MaxVal bigmask",maxval(abs(bigmask))
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_REDUCE(bigmask,temp_big,nx0_big*ny0_big*local_N,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,rank,MPI_COMM_WORLD,ierr)

     ! print *, "MaxVal temp_big",maxval(abs(temp_big))
     
  ENDDO
endif

if (printpad) CALL WRITE_PADDING(1)
  if (verbose) print *, "Padded"
  CALL fftw_mpi_execute_dft_c2r(plan_c2r,temp_big,store)
  if (verbose) print *, "Through IRFFT"

END SUBROUTINE ZEROPAD

SUBROUTINE UNPACK

  IMPLICIT NONE

  integer :: i
  integer :: lkz1_rank,lkz2_rank,ind,rank,kp,k
  logical :: zmask

  if (verbose) print *, "Entering Unpack"

  CALL fftw_mpi_execute_dft_r2c(plan_r2c,store,temp_big)
  if (verbose) print *, "Through RFFT"
  temp_small = cmplx(0.0,0.0)

 if (padv1) then
  
  DO i = 0,nkx0-1
     temp_small(i,0:hky_ind,0:hkz_ind) = temp_big(i+1,1:1+hky_ind,1:hkz_ind+1)*fft_norm
     !   if (verbose) print *, i,"++ done"
    temp_small(i,0:hky_ind,lkz_ind:nkz0-1) = temp_big(i+1,1:1+hky_ind,1+lkz_big:nz0_big)*fft_norm
    !   if (verbose) print *, i,"+- done"
    temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) = temp_big(i+1,1+lky_big:ny0_big,1+lkz_big:nz0_big)*fft_norm
    !   if (verbose) print *, i,"-- done"
    temp_small(i,lky_ind:nky0-1,0:hkz_ind) = temp_big(i+1,1+lky_big:ny0_big,1:hkz_ind+1)*fft_norm
    !   if (verbose) print *, i,"-+ done"
 ENDDO

endif
  
if (padv2) then

    DO rank = 0,n_mpi_procs-1

       lkz1_rank = rank*(nkz0)/n_mpi_procs
       lkz2_rank = (rank+1)*(nkz0)/n_mpi_procs-1

       DO k = 1,local_N

          ! Find full position in large array
          kp = local_k_offset + k-1
          
          ! Mask away values from dealiasing
          zmask = (((kp.ge.lkz_big).or.(kp.le.hkz_ind)))

          ! Pull kp down to original array size if over half way
          if (kp.ge.lkz_big) kp = kp - nkz0/2
          ! Now take out values that aren't compatible with rank indices
          zmask = (zmask.and.((kp.ge.lkz1_rank).and.(kp.le.lkz2_rank)))
          
          ! Find position of ind in smallmask - lkz1_rank is min in rank
          ! If the minimum below is lkz2_rank and not both equal, then zmask will be zero so works
          ind = max(min(kp-lkz1_rank,lkz2_rank),lkz1_rank)
          
          if (zmask) smallmask(0:nkx0-1,0:hky_ind,ind) = temp_big(1:nkx0,1:1+hky_ind,k)*fft_norm
          if (zmask) smallmask(0:nkx0-1,lky_ind:nky0-1,ind) = temp_big(1:nkx0,1+lky_big:ny0_big,k)*fft_norm
        
       ENDDO

     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_REDUCE(smallmask,temp_small,nkx0*nky0*nkz0/n_mpi_procs,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,rank,MPI_COMM_WORLD,ierr)

  ENDDO

endif

if (printunpack) CALL WRITE_PADDING(2)
  
  if (verbose) print *, "All Done"

END SUBROUTINE UNPACK

SUBROUTINE WRITE_PADDING(purpose)

  IMPLICIT NONE

  INTEGER, intent(in) :: purpose
  INTEGER :: chp_handle
  CHARACTER(len=100) :: chp_name

  CALL get_io_number
  chp_handle = io_number

  if (purpose.eq.1) chp_name = "/padding.dat"
  if (purpose.eq.2) chp_name = "/unpack.dat"

  OPEN(unit=chp_handle,file=trim(diagdir)//trim(chp_name),&
       form='unformatted', status='replace',access='stream')

  if (purpose.eq.1) WRITE(chp_handle) temp_big
  if (purpose.eq.2) WRITE(chp_handle) temp_small

  CLOSE(chp_handle)

  if (purpose.eq.1) printpad = (.false.)
  if (purpose.eq.2) printunpack = (.false.)

END SUBROUTINE WRITE_PADDING


END MODULE nonlinearity

