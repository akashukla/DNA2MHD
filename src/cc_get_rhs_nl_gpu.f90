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
  include 'fftw3.f03'
  include 'mpif.h'

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  PRIVATE

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small
  REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:,:,:,:) :: realarrays

  !For fft's


  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: temp_big(:,:,:)
  REAL(C_DOUBLE), allocatable :: store(:,:,:)
  
  !COMPLEX(C_DOUBLE_COMPLEX), pointer :: temp_big(:,:,:),temp_big1(:,:,:)
  !REAL(C_DOUBLE), pointer ::  store(:,:,:),store1(:,:,:)

  type(C_PTR) :: plan_r2c,plan_c2r,rdata,cdata
  
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
  
  ! Initialize FFT plans
  implicit none

  t1 = MPI_WTIME()

  ! This was convenient from P3DFFT, we'll keep the cstart/cend/rstart/rend arrays for indexing
  cstart = [1,1,1]
  cend = [1+nx0_big/2,ny0_big,nz0_big]
  rstart = [1,1,1]
  rend = [nx0_big,ny0_big,nz0_big]
  if (verbose) print *, mype,cstart(1),cend(1),cstart(2),cend(2),cstart(3),cend(3)

  ALLOCATE(temp_big(1:1+nx0_big/2,1:ny0_big,1:nz0_big))
  ALLOCATE(store(1:nx0_big,1:ny0_big,1:nz0_big))
  
!  cdata = fftw_alloc_complex(int((nx0_big/2 + 1)*ny0_big*nz0_big,C_SIZE_T))
!  call c_f_pointer(cdata,store,[2*(nx0_big/2+1),ny0_big,nz0_big])
!  call c_f_pointer(cdata,temp_big,[nx0_big/2+1,ny0_big,nz0_big])

!  rdata = fftw_alloc_complex(int((nx0_big/2 + 1)*ny0_big*nz0_big,C_SIZE_T))
!  call c_f_pointer(rdata,store1,[2*(nx0_big/2+1),ny0_big,nz0_big])
!  call c_f_pointer(rdata,temp_big1,[nx0_big/2+1,ny0_big,nz0_big])

  plan_c2r = fftw_plan_dft_c2r_3d(nz0_big,ny0_big,nx0_big,temp_big,store,FFTW_ESTIMATE)
  plan_r2c = fftw_plan_dft_r2c_3d(nz0_big,ny0_big,nx0_big,store,temp_big,FFTW_ESTIMATE)
  
  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  ! This was convenient from P3DFFT, we'll keep the cstart/cend/rstart/rend arrays for indexing
  cstart = [1,1,1]
  cend = [1+nx0_big/2,ny0_big,nz0_big]
  rstart = [1,1,1]
  rend = [nx0_big,ny0_big,nz0_big]
  if (verbose) print *, mype,cstart(1),cend(1),cstart(2),cend(2),cstart(3),cend(3)
    
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
     CALL ZEROPAD
     realarrays(:,:,:,0) = store
     if (verbose.and.(mype.eq.0)) print *, "Through bx"

     !by
     temp_big = b_in(:,:,:,1)
     CALL ZEROPAD
     realarrays(:,:,:,4) = store
     if (verbose.and.(mype.eq.0)) print *, "Through by"
     
     !bz
     temp_big =	b_in(:,:,:,2)
     CALL ZEROPAD
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
     CALL ZEROPAD
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
     CALL ZEROPAD
     realarrays(:,:,:,6) = store
     
     ! curlbz
     temp_big = cmplx(0.0,0.0)
     DO i = cstart(1),cend(1)
        DO j = cstart(2),cend(2)
           temp_big(i,j,:) = i_complex * kxgrid(i) * b_in(i,j,:,1) - i_complex * kygrid(j) * b_in(i,j,:,0)
        ENDDO
     ENDDO
     CALL ZEROPAD
     realarrays(:,:,:,10) = store

     if (verbose.and.(mype.eq.0)) print *, "Through b derivatives"
     
  endif !Skip b FFTs if Navier Stokes 
  
!!! TERMS  vx,vy,vz 
  !vx
  temp_big = v_in(:,:,:,0)
  !Add padding for dealiasing
  CALL ZEROPAD
  realarrays(:,:,:,1) = store
  
  !vy
  temp_big = v_in(:,:,:,1)
  CALL ZEROPAD
  !Add padding for dealiasing
  realarrays(:,:,:,5) = store
  
  !vz
  temp_big = v_in(:,:,:,2)
  CALL ZEROPAD
  realarrays(:,:,:,9) = store
  
  !Add padding for dealiasing

  if (verbose.and.(mype.eq.0)) print *, "Through v FFTs"

  temp_big = cmplx(0.0,0.0)
  DO j = cstart(2),cend(2)
     DO k = cstart(3),cend(3)
        temp_big(:,j,k) = i_complex * kygrid(j) * v_in(:,j,k,2) &
             - i_complex  * kzgrid(k) * v_in(:,j,k,1)
     ENDDO
  ENDDO
  CALL ZEROPAD
  realarrays(:,:,:,3) = store

  temp_big = cmplx(0.0,0.0)
  DO k = cstart(3),cend(3)
     DO i = cstart(1),cend(1)
        temp_big(i,:,k) = i_complex * kzgrid(k) * v_in(i,:,k,0) &
             - i_complex * kxgrid(i) * v_in(i,:,k,2)
     ENDDO
  ENDDO
  CALL ZEROPAD
  realarrays(:,:,:,7) = store

  temp_big = cmplx(0.0,0.0)
  DO i = cstart(1),cend(1)
     DO j = cstart(2),cend(2)
        temp_big(i,j,:) = i_complex * kxgrid(i) * v_in(i,j,:,1) - i_complex * kygrid(j) * v_in(i,j,:,0)
     ENDDO
  ENDDO
  CALL ZEROPAD
  realarrays(:,:,:,11) = store

  if (verbose.and.(mype.eq.0)) print *, "Through derivatives"
  
  if (timer.and.(mype.eq.0)) t2 = MPI_WTIME()
  if (timer.and.(mype.eq.0)) print *, "Time for Derivs and IRFFTs",t2-t1 
  
  if (timer.and.(mype.eq.0)) t1 = MPI_WTIME()
  if (.not.nv) then ! Skip b ffts if Navier Stokes
     
     ! x: vy bz - vz by - (curlby bz - curlbz by)

     store = (realarrays(:,:,:,5) - hall * realarrays(:,:,:,6)) * realarrays(:,:,:,8) &
          - (realarrays(:,:,:,9) - hall * realarrays(:,:,:,10)) * realarrays(:,:,:,4)

     CALL UNPACK
     DO k = cstart(3),cend(3)
        rhs_out_b(:,:,k,1) = rhs_out_b(:,:,k,1) + i_complex * kzgrid(k) * temp_big(:,:,k)
     ENDDO
     DO j = cstart(2),cend(2)
        rhs_out_b(:,j,:,2) = rhs_out_b(:,j,:,2) - i_complex * kygrid(j) * temp_big(:,j,:)
     ENDDO
     
     ! y: vz bx - vx bz - hall * (curlbz bx - curlbx bz) 
     store = (realarrays(:,:,:,9) - hall * realarrays(:,:,:,10)) * realarrays(:,:,:,0) &
          - (realarrays(:,:,:,1) - hall * realarrays(:,:,:,2)) * realarrays(:,:,:,8)

     CALL UNPACK

     DO i = cstart(1),cend(1)
        rhs_out_b(i,:,:,2) = rhs_out_b(i,:,:,2) + i_complex * kxgrid(i) * temp_big(i,:,:)
     ENDDO
     DO k = cstart(3),cend(3)
        rhs_out_b(:,:,k,0) = rhs_out_b(:,:,k,0) - i_complex * kzgrid(k) * temp_big(:,:,k)
     ENDDO

     ! z: vx by - vy bx - hall (curlbx by - curlby bx) 
     store = (realarrays(:,:,:,1) - hall * realarrays(:,:,:,2)) * realarrays(:,:,:,4) &
          - (realarrays(:,:,:,5) - hall * realarrays(:,:,:,6)) * realarrays(:,:,:,0)

     CALL UNPACK
     
     DO j = cstart(2),cend(2)
        rhs_out_b(:,j,:,0) = rhs_out_b(:,j,:,0) + i_complex * kygrid(j) * temp_big(:,j,:)
     ENDDO
     DO i = cstart(1),cend(1)
        rhs_out_b(i,:,:,1) = rhs_out_b(i,:,:,1) - i_complex * kxgrid(i) * temp_big(i,:,:)
     ENDDO
     

  endif

  ! x: vy curlvz - vz curlvy + curlby bz - curlbz by
  store = realarrays(:,:,:,5) * realarrays(:,:,:,11) - realarrays(:,:,:,9) * realarrays(:,:,:,7) &
       + realarrays(:,:,:,6) * realarrays(:,:,:,8) - realarrays(:,:,:,10) * realarrays(:,:,:,4)
  CALL UNPACK
  rhs_out_v(:,:,:,0) = rhs_out_v(:,:,:,0) + temp_big

  ! y: vz * curlvx - vx * curlvz + curlbz * bx - curlbx * bz 
  store = realarrays(:,:,:,9) * realarrays(:,:,:,3) - realarrays(:,:,:,1) * realarrays(:,:,:,11) &
       + realarrays(:,:,:,10) * realarrays(:,:,:,0) - realarrays(:,:,:,2) * realarrays(:,:,:,8)
  CALL UNPACK
  rhs_out_v(:,:,:,1) = rhs_out_v(:,:,:,1) + temp_big

  !z : vx * curlvy - vy * curlvx + curlbx * by - curlby * bx
  store = realarrays(:,:,:,1) * realarrays(:,:,:,7) - realarrays(:,:,:,5) * realarrays(:,:,:,3) &
       + realarrays(:,:,:,2) * realarrays(:,:,:,4) - realarrays(:,:,:,6) * realarrays(:,:,:,0)
  CALL UNPACK
  rhs_out_v(:,:,:,2) = rhs_out_v(:,:,:,2) + temp_big
  
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
  
  CALL fftw_destroy_plan(plan_c2r)
  call fftw_destroy_plan(plan_r2c)
  !CALL fftw_free(rdata)
  !CALL fftw_free(cdata)
  
END SUBROUTINE finalize_fourier

SUBROUTINE ALLOCATIONS
  
  ! ALLOCATE(temp_small(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  
  ! we might be able to get out of needing the real arrays with in place transforms
  ALLOCATE(realarrays(rstart(1):rend(1),rstart(2):rend(2),rstart(3):rend(3),0:11))
  
END SUBROUTINE ALLOCATIONS

SUBROUTINE DEALLOCATIONS
  
  
  ! if (allocated(temp_small)) DEALLOCATE(temp_small)
  if (verbose.and.(mype.eq.0)) print *, 'ts deallocated'
  
  ! All b arrays

  DEALLOCATE(temp_big)
  DEALLOCATE(store)
  
  if (verbose.and.(mype.eq.0)) print *, "all derivatives deallocated"

  if (allocated(realarrays)) DEALLOCATE(realarrays)
  
END SUBROUTINE DEALLOCATIONS

SUBROUTINE ZEROPAD

  CALL fftw_execute_dft_c2r(plan_c2r,temp_big,store)

END SUBROUTINE ZEROPAD


SUBROUTINE UNPACK
  
  IMPLICIT NONE
  
  integer :: i,j,k
  integer :: lkz1_rank,lkz2_rank,ind,rank,kp
  logical :: zmask

  !if (verbose) print *, "Entering Unpack"
  
  CALL fftw_execute_dft_r2c(plan_r2c,store,temp_big)
  
  ! print *, "Post RFFT",maxval(abs(temp_big))
  if (verbose.and.(mype.eq.0)) print *, "Through RFFT"

  temp_big = temp_big * fft_norm * paddingmask
    
  if (verbose.and.(mype.eq.0)) print *, "All Done"
  
END SUBROUTINE UNPACK

END MODULE nonlinearity

