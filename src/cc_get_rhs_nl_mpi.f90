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
  use p3dfft
  
  use, intrinsic :: iso_c_binding

  include 'mpif.h'

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  PRIVATE

  COMPLEX(p3dfft_type), ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0
  COMPLEX(p3dfft_type), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small

  REAL(p3dfft_type), ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz,curlbx,curlby,curlbz
  REAL(p3dfft_type), ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz,curlvx,curlvy,curlvz

  REAL(p3dfft_type), ALLOCATABLE, DIMENSION(:,:,:,:) :: realarrays

  !For fft's

  COMPLEX(p3dfft_type), allocatable :: temp_big(:,:,:)
  REAL(p3dfft_type), allocatable ::  store(:,:,:)
  
  REAL :: fft_norm  !normalization factor for inverse fft
  INTEGER :: i,j,k
  INTEGER :: ierr
  REAL :: t1, t2,t3,t4

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
SUBROUTINE initialize_fourier_ae_mu0
  
  ! Initialize P3DFFT
  implicit none

  integer(4) :: dims(2)
  
  dims(1) = 2
  dims(2) = n_mpi_procs/2
  
  if (mype.eq.0) print *, "Processor Grid",dims

  t1 = MPI_WTIME()
  CALL p3dfft_setup(dims,nx0_big,ny0_big,nz0_big,MPI_COMM_WORLD,overwrite=.true.)

  ! Get Dimensions for Complex
  CALL p3dfft_get_dims(cstart,cend,csize,2) 

  ! Get Dimensions for Real
  CALL p3dfft_get_dims(rstart,rend,rsize,1)

  if (verbose) print *, mype,cstart(1),cend(1),cstart(2),cend(2),cstart(3),cend(3)

  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  CALL ALLOCATIONS

  t2 = MPI_WTIME()

  if (mype.eq.0) print *, "Time for FFT Plan",t2-t1

END SUBROUTINE initialize_fourier_ae_mu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   get_rhs_nl                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl(b_in, v_in, rhs_out_b, rhs_out_v,ndt)
  USE par_mod

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(in) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(in) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  REAL(C_DOUBLE) :: ndt
  
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
  
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(in) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(in) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  REAL(C_DOUBLE) :: ndt

  INTEGER :: l,h, ierr
  REAL(8) :: looprvec(0:11),loopkx,loopky,loopkz
  COMPLEX(8) :: loopcvec(0:2)
  

  ! I dont want to change g_in, so I copy temporaly to g_in0
  !g_in0 = g_in
  
  !IF(mype==0) WRITE(*,*) "Actually in nl1"

! TERMS BX BY BZ


  if (timer.and.(mype.eq.0)) t1 = MPI_WTIME()
  if (.not.nv) then ! Skip b if Navier Stokes
     !bx
     temp_big = b_in(:,:,:,0)
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,0) = store
     if (verbose.and.(mype.eq.0)) print *, "Through bx"
     
     !by
     temp_big = b_in(:,:,:,1)
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,4) = store
     if (verbose.and.(mype.eq.0)) print *, "Through by"
     
     !bz
     temp_big = b_in(:,:,:,2)
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,8) = store     
     if (verbose.and.(mype.eq.0)) print *, "Through bz"
     
     ! curlbx
     temp_big = cmplx(0.0,0.0)
     DO j = cstart(2),cend(2)
        DO k = cstart(3),cend(3)
           temp_big(:,j,k) = i_complex * kygrid(j) * b_in(:,j,k,2) &
                - i_complex  * kzgrid(k) * b_in(:,j,k,1)
        ENDDO
     ENDDO
     if (verbose.and.(mype.eq.0)) print *, "Through assignment"
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,2) = store
     
     if (verbose.and.(mype.eq.0)) print *, "Through curl bx" 
     
     ! curlby
     temp_big = cmplx(0.0,0.0)
     DO k = cstart(3),cend(3)
        DO i = cstart(1),cend(1)
           temp_big(i,:,k) = i_complex * kzgrid(k) * b_in(i,:,k,0) &
                - i_complex * kxgrid(i) * b_in(i,:,k,2)
        ENDDO
     ENDDO
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,6) = store
     
     ! curlbz
     temp_big = cmplx(0.0,0.0)
     DO i = cstart(1),cend(1)
        DO j = cstart(2),cend(2)
           temp_big(i,j,:) = i_complex * kxgrid(i) * b_in(i,j,:,1) - i_complex * kygrid(j) * b_in(i,j,:,0)
        ENDDO
     ENDDO
     CALL p3dfft_btran_c2r(temp_big,store,"fff")
     realarrays(:,:,:,10) = store

     if (verbose.and.(mype.eq.0)) print *, "Through b derivatives"
     
  endif !Skip b FFTs if Navier Stokes 
  
!!! TERMS  vx,vy,vz 
  !vx
  temp_big = v_in(:,:,:,0)
  !Add padding for dealiasing
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,1) = store
  
  !vy
  temp_big = v_in(:,:,:,1)
  !Add padding for dealiasing
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,5) = store
  
  !vz
  temp_big = v_in(:,:,:,2)
  !Add padding for dealiasing
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,9) = store

  if (verbose.and.(mype.eq.0)) print *, "Through v FFTs"
  
  temp_big  = cmplx(0.0,0.0)
  DO j = cstart(2),cend(2)
     DO k = cstart(3),cend(3)
        temp_big(:,j,k) = i_complex * kygrid(j) * v_in(:,j,k,2) &
             - i_complex  * kzgrid(k) * v_in(:,j,k,1)
     ENDDO
  ENDDO
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,3) = store
  
  temp_big = cmplx(0.0,0.0)
  DO k = cstart(3),cend(3)
     DO i = cstart(1),cend(1)
        temp_big(i,:,k) = i_complex * kzgrid(k) * v_in(i,:,k,0) &
             - i_complex * kxgrid(i) * v_in(i,:,k,2)
     ENDDO
  ENDDO
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,7) = store
  
  temp_big = cmplx(0.0,0.0)
  DO i = cstart(1),cend(1)
     DO j = cstart(2),cend(2)
        temp_big(i,j,:) = i_complex * kxgrid(i) * v_in(i,j,:,1) - i_complex * kygrid(j) * v_in(i,j,:,0)
     ENDDO
  ENDDO
  CALL p3dfft_btran_c2r(temp_big,store,"fff")
  realarrays(:,:,:,11) = store
 
  if (verbose.and.(mype.eq.0)) print *, "Through derivatives"
  
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

     store = (realarrays(:,:,:,5) - hall * realarrays(:,:,:,6)) * realarrays(:,:,:,8) - (realarrays(:,:,:,9) - hall * realarrays(:,:,:,10)) * realarrays(:,:,:,4)

     CALL UNPACK
     
     DO k = cstart(3),cend(3)
        rhs_out_b(:,:,k,1) = rhs_out_b(:,:,k,1) + i_complex * kzgrid(k) * temp_small(:,:,k)
     ENDDO
     DO j = cstart(2),cend(2)
        rhs_out_b(:,j,:,2) = rhs_out_b(:,j,:,2) - i_complex * kygrid(j) * temp_small(:,j,:)
     ENDDO
     
     ! y: vz bx - vx bz - hall * (curlbz bx - curlbx bz)
     store = (realarrays(:,:,:,9) - hall * realarrays(:,:,:,10)) * realarrays(:,:,:,0) - (realarrays(:,:,:,1) - hall * realarrays(:,:,:,2)) * realarrays(:,:,:,8) 
     
     CALL UNPACK
     
     DO i = cstart(1),cend(1)
        rhs_out_b(i,:,:,2) = rhs_out_b(i,:,:,2) + i_complex * kxgrid(i) * temp_small(i,:,:)
     ENDDO
     DO k = cstart(3),cend(3)
        rhs_out_b(:,:,k,0) = rhs_out_b(:,:,k,0) - i_complex * kzgrid(k) * temp_small(:,:,k)
     ENDDO
     
     ! z: vx by - vy bx - hall (curlbx by - curlby bx)
     store = (realarrays(:,:,:,1) - hall * realarrays(:,:,:,2)) * realarrays(:,:,:,4) - (realarrays(:,:,:,5) - hall * realarrays(:,:,:,6)) * realarrays(:,:,:,0)
     
     CALL UNPACK
     
     DO j = cstart(2),cend(2)
        rhs_out_b(:,j,:,0) = rhs_out_b(:,j,:,0) + i_complex * kygrid(j) * temp_small(:,j,:)
     ENDDO
     DO i = cstart(1),cend(1)
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
  
  store = realarrays(:,:,:,5) * realarrays(:,:,:,11) - realarrays(:,:,:,9) * realarrays(:,:,:,7) + realarrays(:,:,:,6) * realarrays(:,:,:,8) - realarrays(:,:,:,10) * realarrays(:,:,:,4)  
  
  CALL UNPACK
  
  rhs_out_v(:,:,:,0) = rhs_out_v(:,:,:,0) + temp_small
  
  ! y: vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz
  store = realarrays(:,:,:,9) * realarrays(:,:,:,3) - realarrays(:,:,:,1) * realarrays(:,:,:,11) + realarrays(:,:,:,10) * realarrays(:,:,:,0) - realarrays(:,:,:,2) * realarrays(:,:,:,8)
  
  CALL UNPACK
  if (verbose.and.(mype.eq.0)) print *, "vy fft complete"
  rhs_out_v(:,:,:,1) = rhs_out_v(:,:,:,1) + temp_small
  if (verbose.and.(mype.eq.0)) print *, "vy done"
  
  ! z: vx * curlvy - vy * curlvx + curlbx * by - curlby * bx
  store = realarrays(:,:,:,1) * realarrays(:,:,:,7) - realarrays(:,:,:,5) * realarrays(:,:,:,3) + realarrays(:,:,:,2) * realarrays(:,:,:,4) - realarrays(:,:,:,6) * realarrays(:,:,:,0)
  
  CALL UNPACK
  if (verbose.and.(mype.eq.0)) print *, "vz fft complete"
  rhs_out_v(:,:,:,2) = rhs_out_v(:,:,:,2) + temp_small
  if (verbose.and.(mype.eq.0)) print *, "vz done"
  
  if (timer.and.(mype.eq.0)) t2 = MPI_WTIME()
  if (timer.and.(mype.eq.0)) print *, "equations and RFFTs",t2-t1
  
  ! to preserve reality of the fields, remove v,b terms at nky0/2,nkz0/2 - padding mask should take care of this
  
  if (verbose.and.(mype.eq.0)) print *, 'rhs out v nl found'
  
  CALL next_dt(ndt)
  
  if ((mod(itime,100).eq.0).and.mype.eq.0) print *, 'next dt calculated ',ndt
  if (.not.(calc_dt)) ndt = dt_max
  
END SUBROUTINE get_rhs_nl1

SUBROUTINE next_dt(dtn)
  
  real(C_DOUBLE), intent(out) :: dtn
  
  real :: ndt1xr,ndt1yr,ndt1zr,ndt2xr,ndt2yr,ndt2zr,ndt3xr,ndt3yr,ndt3zr
  real :: ndtr
  
  ndt1xr = maxval(abs(kxgrid))*maxval(abs(realarrays(:,:,:,0)))
  ndt1yr = maxval(abs(kygrid))*maxval(abs(realarrays(:,:,:,4)))
  ndt1zr = maxval(abs(kzgrid))*maxval(abs(realarrays(:,:,:,8)))
  ndt2xr = maxval(abs(kxgrid))*maxval(abs(realarrays(:,:,:,1)))
  ndt2yr = maxval(abs(kygrid))*maxval(abs(realarrays(:,:,:,5)))
  ndt2zr = maxval(abs(kzgrid))*maxval(abs(realarrays(:,:,:,9)))
  ndt3xr = maxval(abs(kxgrid))*maxval(abs(realarrays(:,:,:,2)))*hall
  ndt3yr = maxval(abs(kygrid))*maxval(abs(realarrays(:,:,:,6)))*hall
  ndt3zr = maxval(abs(kzgrid))*maxval(abs(realarrays(:,:,:,10)))*hall
  ndtr = ndt1xr + ndt1yr + ndt1zr &
       + ndt2xr + ndt2yr + ndt2zr &
       + ndt3xr + ndt3yr + ndt3zr
  dtn = courant/ndtr
  
END SUBROUTINE next_dt

SUBROUTINE finalize_fourier
  
  implicit none
  
  CALL DEALLOCATIONS
  
  CALL p3dfft_clean
  
END SUBROUTINE finalize_fourier

SUBROUTINE ALLOCATIONS
  
  ALLOCATE(temp_big(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  ALLOCATE(store(rstart(1):rend(1),rstart(2):rend(2),rstart(3):rend(3)))
  
  ALLOCATE(temp_small(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

  ! bx, vx, curlbx, curlvx, by, vy, curlby, curlvy, bz, vz, curlbz, curlvz 
  ALLOCATE(realarrays(rstart(1):rend(1),rstart(2):rend(2),rstart(3):rend(3),0:11))
  
END SUBROUTINE ALLOCATIONS

SUBROUTINE DEALLOCATIONS
  
  
  if (allocated(temp_big)) DEALLOCATE(temp_big)
  if (allocated(store)) DEALLOCATE(store)
  
  if (allocated(temp_small)) DEALLOCATE(temp_small)
  if (verbose.and.(mype.eq.0)) print *, 'ts deallocated'
  
  ! All b arrays
  
  if (verbose.and.(mype.eq.0)) print *, "all derivatives deallocated"

  if (allocated(realarrays)) DEALLOCATE(realarrays)
  
END SUBROUTINE DEALLOCATIONS

SUBROUTINE UNPACK
  
  IMPLICIT NONE
  
  integer :: i,j,k
  integer :: lkz1_rank,lkz2_rank,ind,rank,kp
  logical :: zmask

  !if (verbose) print *, "Entering Unpack"
  
  CALL p3dfft_ftran_r2c(store,temp_big,"fff")
  
  ! print *, "Post RFFT",maxval(abs(temp_big))
  if (verbose.and.(mype.eq.0)) print *, "Through RFFT"
  
  temp_small = temp_big * paddingmask * fft_norm
  
  if (verbose.and.(mype.eq.0)) print *, "All Done"
  
END SUBROUTINE UNPACK

END MODULE nonlinearity

