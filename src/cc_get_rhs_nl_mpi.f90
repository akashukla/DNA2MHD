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
  !USE mpi
  USE par_mod
  !USE hk_effects
  !USE flr_effects
  
  use, intrinsic :: iso_c_binding
  
  include 'fftw3-mpi.f03'
  include 'mpif.h'

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  PRIVATE

  COMPLEX(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,bigmask,smallmask

  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz,curlbx,curlby,curlbz
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz,curlvx,curlvy,curlvz

  !For fft's

  INTEGER(C_INTPTR_T) :: local_N,local_k_offset,alloc_local
  
  type(C_PTR) :: plan_r2c,plan_c2r
  type(C_PTR) :: rinout,cinout

  type(C_PTR) :: plan_1c2r,plan_c2c
  type(C_PTR) :: rout1,cin1x,cin1y
  COMPLEX(C_DOUBLE_COMPLEX), pointer :: tempbig1(:),testyz(:)
  REAL(C_DOUBLE), pointer :: storex(:)

  COMPLEX(C_DOUBLE_COMPLEX), pointer :: temp_big(:,:,:)
  REAL(C_DOUBLE), pointer ::  store(:,:,:)
  
  REAL :: fft_norm  !normalization factor for inverse fft
  INTEGER(C_INTPTR_T) :: i,j,k
  INTEGER :: ierr
  REAL :: t1, t2,t3,t4

  LOGICAL :: padv1,padv2,padv3,padv4
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
  
  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))

  alloc_local = fftw_mpi_local_size_3d(nz0_big,ny0_big,nx0_big/2+1,MPI_COMM_WORLD,local_N,local_k_offset)
  rinout = fftw_alloc_real(2*alloc_local)
  cinout = fftw_alloc_complex(alloc_local)
  
  CALL c_f_pointer(rinout,store,[2*(nx0_big/2+1),ny0_big,local_N])
  CALL c_f_pointer(cinout,temp_big,[nx0_big/2+1,ny0_big,local_N])

  plan_c2r = fftw_mpi_plan_dft_c2r_3d(nz0_big,ny0_big,nx0_big,temp_big,store,MPI_COMM_WORLD,FFTW_PATIENT)
  plan_r2c = fftw_mpi_plan_dft_r2c_3d(nz0_big,ny0_big,nx0_big,store,temp_big,MPI_COMM_WORLD,FFTW_PATIENT)
     
  !WRITE(*,*) "making plans"

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big
  WRITE(*,*) "local_N, lkz2+1-lkz1",local_N,lkz2+1-lkz1,mype
  WRITE(*,*) "local_k_offset",local_k_offset,mype

  CALL ALLOCATIONS

  printpad = (.false.)
  printunpack = (.false.)

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp,treal)
  !CALL dfftw_execute_dft_r2c(plan_r2c,treal,tcomp)
  
  !tcomp=tcomp/REAL(n1*n2*n3)

END SUBROUTINE initialize_fourier_ae_mu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   get_rhs_nl                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl(b_in, v_in, rhs_out_b, rhs_out_v,ndt)
  USE par_mod

  COMPLEX, INTENT(in) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
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
  
  COMPLEX, INTENT(in) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
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


  if (timer.and.(mype.eq.0)) t1 = MPI_WTIME()
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
    DO j = 0,ny0_big-1
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
       DO j = 0,ny0_big-1
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
    DO j = 0,ny0_big-1
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
       DO j = 0,ny0_big-1
          temp_small(i,j,:) = i_complex * kxgrid(i) * v_iny0(i,j,:) - i_complex * kygrid(j) * v_inx0(i,j,:)
       ENDDO
    ENDDO
    CALL ZEROPAD    
    curlvz = store

    if (timer.and.(mype.eq.0)) t2 = MPI_WTIME()
    if (timer.and.(mype.eq.0)) print *, "Time for Derivs and IRFFTs",t2-t1 
    
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

    if (timer.and.(mype.eq.0)) t1 = MPI_WTIME()
    if (.not.nv) then ! Skip b ffts if Navier Stokes

       ! x: vy bz - vz by - (curlby bz - curlbz by)
       store = vy * bz - vz * by - hall * (curlby * bz - curlbz * by)
       
       CALL UNPACK

       DO k = lkz1,lkz2
          rhs_out_b(:,:,k,1) = rhs_out_b(:,:,k,1) + i_complex * kzgrid(k) * temp_small(:,:,k)
       ENDDO
       DO j = 0,ny0_big-1
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

       DO j = 0,ny0_big-1
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
! (If compressibility gets included, then would handle v^2 pseudospectrally with one more FFT (or rewrite whole as grad))

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

if (timer.and.(mype.eq.0)) t2 = MPI_WTIME()
if (timer.and.(mype.eq.0)) print *, "equations and RFFTs",t2-t1

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

! if (calc_dt) print *, 'next dt calculated ',ndt

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

  ALLOCATE(temp_small(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

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

  if (padv2) then
  ALLOCATE(bigmask(1:nx0_big/2+1,1:ny0_big,1:local_N))
  ALLOCATE(smallmask(0:nkx0-1,0:nky0-1,0:nkz0/n_mpi_procs-1))
endif

  
  ALLOCATE(b_inx0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  ALLOCATE(b_iny0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  ALLOCATE(b_inz0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  ALLOCATE(v_inx0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  ALLOCATE(v_iny0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
  ALLOCATE(v_inz0(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

END SUBROUTINE ALLOCATIONS

SUBROUTINE DEALLOCATIONS


  if (allocated(temp_small)) DEALLOCATE(temp_small)
if (verbose) print *, 'ts deallocated'

! All b arrays

if (allocated(bx)) DEALLOCATE(bx)
if (verbose) print *, 'bx deallocated'
if (allocated(by)) DEALLOCATE(by)
if (verbose) print *, 'by deallocated'
if (allocated(bz)) DEALLOCATE(bz)
if (verbose) print *, 'first third deallocated'

if (allocated(curlbx)) DEALLOCATE(curlbx)
if (allocated(curlby)) DEALLOCATE(curlby)
if (allocated(curlbz)) DEALLOCATE(curlbz)

! All v arrays 
if (allocated(vx))  DEALLOCATE(vx)
if (allocated(vy))  DEALLOCATE(vy)
if (allocated(vz))  DEALLOCATE(vz)

if (allocated(curlvx))  DEALLOCATE(curlvx)
if (allocated(curlvy))  DEALLOCATE(curlvy)
if (allocated(curlvz))  DEALLOCATE(curlvz)

  if (verbose) print *, "all derivatives deallocated"

if (allocated(bigmask))  DEALLOCATE(bigmask)

  if (verbose) print *, "big mask deallocated"
if (allocated(smallmask))  DEALLOCATE(smallmask)

  if (verbose) print *, "small mask deallocated"
  
if (allocated(b_inx0))  DEALLOCATE(b_inx0)
if (allocated(b_iny0))  DEALLOCATE(b_iny0)
if (allocated(b_inz0))  DEALLOCATE(b_inz0)
if (allocated(v_inx0))  DEALLOCATE(v_inx0)
if (allocated(v_iny0))  DEALLOCATE(v_iny0)
if (allocated(v_inz0))  DEALLOCATE(v_inz0)

END SUBROUTINE DEALLOCATIONS

SUBROUTINE ZEROPAD

  IMPLICIT NONE

  integer :: lkz1_rank,lkz2_rank
  integer :: rank,k,kp,i
  integer :: ind
  logical :: zmask
  integer :: llkz1,llkz2

  if (verbose) print *, "In Zeropad","Mype",mype,"MaxVal temp_small",maxval(abs(temp_small)),maxloc(abs(temp_small))
  temp_big=cmplx(0.0,0.0)
  
  temp_big = temp_small
  if (verbose) print *, "In Zeropad","Mype",mype,"MaxVal temp_big",maxval(abs(temp_big)),maxloc(abs(temp_big))

  if (printpad) CALL WRITE_PADDING(1)
  !if (verbose) print *, "Padded"
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL fftw_mpi_execute_dft_c2r(plan_c2r,temp_big,store)
  if (verbose) print *, "In Zeropad","Mype",mype,"MaxVal store",maxval(abs(store)),maxloc(abs(store))
  ! if (verbose) print *, "Through IRFFT"

END SUBROUTINE ZEROPAD

SUBROUTINE UNPACK

  IMPLICIT NONE

  integer :: i
  integer :: lkz1_rank,lkz2_rank,ind,rank,kp,k
  logical :: zmask

  !if (verbose) print *, "Entering Unpack"

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL fftw_mpi_execute_dft_r2c(plan_r2c,store,temp_big)
  if (verbose) print *, "In Unpack","Mype",mype,"MaxVal tempbig",maxval(abs(temp_big))

  ! print *, "Post RFFT",maxval(abs(temp_big))
  !if (verbose) print *, "Through RFFT"
  temp_small = cmplx(0.0,0.0)
  fullsmallarray = cmplx(0.0,0.0)

   CALL clear_padding(temp_big,temp_small)

   temp_small = temp_small * fft_norm
   if (verbose) print *, "In Unpack","Mype",mype,"MaxVal tempsmall",maxval(abs(temp_small))
   
   if (printunpack) CALL WRITE_PADDING(2)
  
  ! if (verbose) print *, "All Done"

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
  
   if (purpose.eq.1) then 
      if (padv1) WRITE(chp_handle) temp_big
      if (padv3.and.(mype.eq.0)) WRITE(chp_handle) fullbigarray
   endif

   if (purpose.eq.2) then
      if (padv1) WRITE(chp_handle) temp_small
      if (padv3.and.(mype.eq.0)) WRITE(chp_handle) fullsmallarray
   endif
  

  CLOSE(chp_handle)

  if (purpose.eq.1) printpad = (.false.)
  if (purpose.eq.2) printunpack = (.false.)

END SUBROUTINE WRITE_PADDING

SUBROUTINE TESTING_1D

  IMPLICIT NONE

  INTEGER :: i,j,k

  type(C_PTR) :: planx,plan_c2c
  type(C_PTR) :: rout1,cin1x,cin1y
  COMPLEX(C_DOUBLE_COMPLEX), pointer :: tempbig1(:),testyz(:)
  REAL(C_DOUBLE), pointer :: storex(:)
  INTEGER(C_INTPTR_T) :: allx,allyz
  
  
  
  
  DO i = 0,nx0_big/2
     DO j = 0,ny0_big-1
        DO k = 0,nz0_big-1

        ENDDO
     ENDDO
  ENDDO
  

  

END SUBROUTINE TESTING_1D


END MODULE nonlinearity

