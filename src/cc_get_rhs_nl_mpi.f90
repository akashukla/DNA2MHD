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
            get_k_indices,get_rhs_nl2,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl_convolution
  
  REAL, PUBLIC :: ve_max(2)

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: temp_small, temp_bigx, temp_bigy, temp_bigz,temp_bigbx, temp_bigby, temp_bigbz,temp_bigvx, temp_bigvy, temp_bigvz,temp_smallm,temp_bigm

  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), pointer :: data,datai
  TYPE(C_PTR) :: plan_c2r,plan_r2c
  TYPE(C_PTR) :: cdatai,cdata
  INTEGER(C_INTPTR_T) :: alloc_local,local_N,local_j_offset,alloc_locali,local_Ni,local_j_offseti,i,j,k,l,h

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) ::  store_bx, store_by, store_bz,store_bxm,store_bym,store_bzm,store_vx, store_vy, store_vz,store_vxm,store_vym,store_vzm
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:,:) ::  bdv,vdb,bdcb,cbdb,vdv,bdb,db2

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz, dxbx, dybx,dzbx, dxby,dyby,dzby, dxbz,dybz,dzbz
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz, dxvx, dyvx,dzvx, dxvy,dyvy,dzvy, dxvz,dyvz,dzvz

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbx, dxdybx, dxdzbx, dydybx, dydzbx, dzdzbx
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxby, dxdyby, dxdzby, dydyby, dydzby, dzdzby
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbz, dxdybz, dxdzbz, dydybz, dydzbz, dzdzbz
    
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvx, dxdyvx, dxdzvx, dydyvx, dydzvx, dzdzvx
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvy, dxdyvy, dxdzvy, dydyvy, dydzvy, dzdzvy
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvz, dxdyvz, dxdzvz, dydyvz, dydzvz, dzdzvz

  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: output,outputm

  !For fft's

  INTEGER(C_INTPTR_T) :: nx0_big,ny0_big,nz0_big
  REAL(C_DOUBLE) :: fft_norm  !normalization factor for inverse fft
  INTEGER :: zpad,counter
  INTEGER :: ierr
 
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

  IMPLICIT NONE
  

  ! determine amount of padding for dealiasing
  if (dealias_type.eq.3) then
  nx0_big = 2*nkx0
  ny0_big = 2*nky0
  nz0_big = 4*nkz0
  endif
  if (dealias_type.eq.4) then 
  nx0_big=2*nkx0
  ny0_big=2*nky0
  nz0_big=2*nkz0
  endif

  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))
  counter = 0

  ! Set Up Modern FFTW Plan

  alloc_local = fftw_mpi_local_size_3d(nz0_big, ny0_big, nx0_big, MPI_COMM_WORLD, &
    local_N, local_j_offset)
  cdata = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cdata, data, [nx0_big,ny0_big,local_N])
  ALLOCATE(data(1:nx0_big,1:ny0_big,1:local_N))
  ! create MPI plan for in-place forward DFT (note dimension reversal)
  plan_r2c = fftw_mpi_plan_dft_3d(nz0_big,ny0_big,nx0_big, data, data, MPI_COMM_WORLD, &
    FFTW_FORWARD, FFTW_ESTIMATE)
  ! backward plan
  alloc_locali = fftw_mpi_local_size_3d(nz0_big,ny0_big,nx0_big,MPI_COMM_WORLD, &
       local_Ni,local_j_offseti)
  cdatai = fftw_alloc_complex(alloc_locali)
  call c_f_pointer(cdatai,datai,[nx0_big,ny0_big,local_Ni])
  ALLOCATE(datai(1:nx0_big,1:ny0_big,1:local_Ni))
  plan_c2r = fftw_mpi_plan_dft_3d(nz0_big,ny0_big,nx0_big,datai,datai,MPI_COMM_WORLD, &
       FFTW_BACKWARD,FFTW_ESTIMATE)

  lky_big=ny0_big-hky_ind !Index of minimum (most negative) FILLED ky value for big arrays  
  lkz_big=nz0_big-hkz_ind !Index of minimum (most negative) FILLED kz value for big arrays 

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  ! Practice FFTs
!  do k = 1, local_N
!    do j = 1,ny0_big
!    do i = 1, nx0_big
!       data(i, j,k) = 1.0/real(i) + 1.0/(real(k+local_j_offset))+1.0/real(j)
!   end do
!   end do 
!  end do
!

  print *,'Is local_N local_Ni?',(local_N.eq.local_Ni),local_N,local_Ni
  print *,'Are offsets 0?',local_j_offset,local_j_offseti

!  print *, maxval(abs(data))
!
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)                                                               
!  call fftw_mpi_execute_dft(plan_r2c, data, data)
!  print *, 'Through FFFT'
!  print *, 'Max FFTd data',maxval(abs(data))
!
!  do k = 1, local_Ni
!    do j = 1,ny0_big
!    do i = 1, nx0_big
!       datai(i, j,k) = data(i,j,k)**2
!   end do
!   end do
!  end do
!
!  print *,maxval(abs(datai))

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fftw_mpi_execute_dft(plan_c2r,datai,datai)
  print *, 'Through IFFT'

  ! ALLOCATIONS
  ALLOCATE(temp_small(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(temp_smallm(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(temp_bigm(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(temp_bigx(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(temp_bigy(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(temp_bigz(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(temp_bigbx(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigby(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigbz(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigvx(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigvy(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigvz(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(store_bxm(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_bym(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_bzm(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vxm(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vym(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vzm(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vx(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vy(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_vz(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_bx(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_by(1:nx0_big,1:ny0_big,1:local_N))
  ALLOCATE(store_bz(1:nx0_big,1:ny0_big,1:local_N))

  ALLOCATE(output(1:nx0_big,1:ny0_big,1:nz0_big))
 
! All b arrays  
  ALLOCATE(bx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(by(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(bz(1:nx0_big,1:ny0_big,1:local_Ni))

! All v arrays 
  ALLOCATE(vx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(vy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(vz(1:nx0_big,1:ny0_big,1:local_Ni))
! all first order v arrays  
  !vx                                                                                                                                                                                   
  ALLOCATE(dxvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dyvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzvx(1:nx0_big,1:ny0_big,1:local_Ni))
  !vy 
  ALLOCATE(dxvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dyvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzvy(1:nx0_big,1:ny0_big,1:local_Ni))
  !vz  
  ALLOCATE(dxvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dyvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzvz(1:nx0_big,1:ny0_big,1:local_Ni))
! all first order b arrays 
 !bx
  ALLOCATE(dxbx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dybx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzbx(1:nx0_big,1:ny0_big,1:local_Ni))
  !by
  ALLOCATE(dxby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dyby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzby(1:nx0_big,1:ny0_big,1:local_Ni))
  !bz
  ALLOCATE(dxbz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dybz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzbz(1:nx0_big,1:ny0_big,1:local_Ni))
  
! all  second order bx arrays  DXDXBX,   DXDYBX,   DXDZBX,  DYDYBX,   DYDZBX, DZDZBX
  ALLOCATE(dxdxbx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdybx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzbx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydybx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzbx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzbx(1:nx0_big,1:ny0_big,1:local_Ni))
  ! all  second order by arrays
  ALLOCATE(dxdxby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdyby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydyby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzby(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzby(1:nx0_big,1:ny0_big,1:local_Ni))
  ! all  second order bz arrays
  ALLOCATE(dxdxbz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdybz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzbz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydybz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzbz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzbz(1:nx0_big,1:ny0_big,1:local_Ni))
    
  ! all  second order vx arrays i.e. DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  ALLOCATE(dxdxvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdyvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydyvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzvx(1:nx0_big,1:ny0_big,1:local_Ni))
  ! all  second order vy arrays
  ALLOCATE(dxdxvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdyvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydyvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzvy(1:nx0_big,1:ny0_big,1:local_Ni))
  ! all  second order bz arrays
  ALLOCATE(dxdxvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdyvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dxdzvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydyvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dydzvz(1:nx0_big,1:ny0_big,1:local_Ni))
  ALLOCATE(dzdzvz(1:nx0_big,1:ny0_big,1:local_Ni))

  ! nonlinearities
  IF(mod(counter,4).eq.0) THEN
  ALLOCATE(bdv(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(vdb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(bdcb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(cbdb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(vdv(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(bdb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(db2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ENDIF

  ! initial arrays
  ALLOCATE(b_inx0(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(b_iny0(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(b_inz0(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(v_inx0(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(v_iny0(1:nkx0,1:nky0,1:nkz0))
  ALLOCATE(v_inz0(1:nkx0,1:nky0,1:nkz0))
 
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
  
    CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
    counter = counter + 1
 
END SUBROUTINE get_rhs_nl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v,ndt)

 USE par_mod
  

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX :: rhs_out_nlb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX :: rhs_out_nlv(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2) 
  COMPLEX :: rhs_out_nlbm(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2),rhs_out_nlvm(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt,myndt
  REAL :: sttime,dum
  REAL :: rsnls

  IF (verbose.and.(mype.eq.0)) CALL cpu_time(sttime)
  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL Start Time:',dum - sttime
  ! I dont want to change g_in, so I copy temporaly to g_in0
  !g_in0 = g_in
  b_inx0(1:nkx0,1:nky0,1:nkz0) = b_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0)
  b_iny0(1:nkx0,1:nky0,1:nkz0) = b_in(0:nkx0-1,0:nky0-1,0:nkz0-1,1)
  b_inz0(1:nkx0,1:nky0,1:nkz0) = b_in(0:nkx0-1,0:nky0-1,0:nkz0-1,2)
  v_inx0(1:nkx0,1:nky0,1:nkz0) = v_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0)
  v_iny0(1:nkx0,1:nky0,1:nkz0) = v_in(0:nkx0-1,0:nky0-1,0:nkz0-1,1)
  v_inz0(1:nkx0,1:nky0,1:nkz0) = v_in(0:nkx0-1,0:nky0-1,0:nkz0-1,2)
  if (verbose.and.(mype.eq.0)) print *,'Through Initialization'
  
  if (nv.eq..false.) then

  ! START SECOND ORDER b terms

  ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
  !dxdxbx
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*i_complex*kxgrid(i-1)*b_inx0(i,:,:) ! there is  two i's in the
  END DO
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&                                        
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (verbose.and.(mype.eq.0)) print *,'Through temp small'

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdxbx)

  !dxdybx
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO j=0,nky0-1
        temp_smallm(i,j,:)=i_complex*kxgrid(i-1)*i_complex*kygrid(j-1)*b_inx0(i,j,:)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdybx)

  !dxdzbx
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(i,:,k)=i_complex*kxgrid(i-1)*i_complex*kzgrid(k-1)*b_inx0(i,:,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdzbx)

  ! DYDYbX
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*i_complex*kygrid(j-1)*b_inx0(:,j,:)  ! towo y grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydybx)
 
  ! DYDzbX
  DO j=mype+1,nky0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(:,j,k)=i_complex*kygrid(j-1)*i_complex*kzgrid(k-1)*b_inx0(:,j,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydzbx)

  ! DzDzbX
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*i_complex*kzgrid(k-1)*b_inx0(:,:,k)  !two kz grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzdzbx)

  ! FINISHED SECOND ORDER  BX TERMS

  ! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY

  !dxdxby
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*i_complex*kxgrid(i-1)*b_iny0(i,:,:) ! there is  two i's in the
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdxby)

  !dxdyby
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO j=0,nky0-1
        temp_smallm(i,j,:)=i_complex*kxgrid(i-1)*i_complex*kygrid(j-1)*b_iny0(i,j,:)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdyby)

  !dxdzby
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(i,:,k)=i_complex*kxgrid(i-1)*i_complex*kzgrid(k-1)*b_iny0(i,:,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdzby)

  ! DYDYby
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*i_complex*kygrid(j-1)*b_iny0(:,j,:)  ! two y grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydyby)

  ! DYDzbY
  DO j=mype+1,nky0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(:,j,k)=i_complex*kygrid(j-1)*i_complex*kzgrid(k-1)*b_iny0(:,j,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydzby)
   
  ! DzDzbY
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*i_complex*kzgrid(k-1)*b_iny0(:,:,k)  !two kz grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzdzby)

  ! FINISHED SECOND ORDER bY TERMS

  ! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
  !dxdxbz

  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*i_complex*kxgrid(i-1)*b_inz0(i,:,:) ! there is  two i's in the
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdxbz)

  !dxdybz
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO j=0,nky0-1
        temp_smallm(i,j,:)=i_complex*kxgrid(i-1)*i_complex*kygrid(j-1)*b_inz0(i,j,:)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdybz)

  !dxdzbz
  DO i=mype+1,nkx0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(i,:,k)=i_complex*kxgrid(i-1)*i_complex*kzgrid(k-1)*b_inz0(i,:,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxdzbz)

  ! DYDYbz
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*i_complex*kygrid(j-1)*b_inz0(:,j,:)  ! towo y grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydybz)

  ! DYDzbz
  DO j=mype+1,nky0-1,n_mpi_procs
     DO k=0,nkz0-1
        temp_smallm(:,j,k)=i_complex*kygrid(j-1)*i_complex*kzgrid(k-1)*b_inz0(:,j,k)
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dydzbz)

  ! DzDzbz
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*i_complex*kzgrid(k-1)*b_inz0(:,:,k)  !two kz grid
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzdzbz)

  !finished END SECOND ORDER BZ TERMS

  !    completed  ALL SECOND ORDER B TERMS

  ! TERMS BX BY BZ

  !bx
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=b_inx0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,bx)

  !by
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=b_iny0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,by)

  !bz
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=b_inz0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,bz)

  endif

!!! TERMS  vx,vy,vz
  !vx
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=v_inx0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,vx)

  !vy
  DO j=mype+1,nky0-1,n_mpi_procs
     !temp_small(i,:,:)=i_complex*kxgrid(i-1)*phi_in(i,:,:)
     temp_smallm(:,j,:)=v_iny0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,vy)

  !vz
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=v_inz0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,vz)

  !  FIRST ORDER VX TERMS DXVX , DYVX,  DZVX

  ! dxvx
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*v_inx0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxvx)

  ! dyvx
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*v_inx0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dyvx)

  ! dzvx
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*v_inx0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzvx)

  !  FIRST ORDER VY TERMS  dxvy dyvy dzvz,

  ! dxvy
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*v_iny0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxvy)

  ! dyvy
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*v_iny0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dyvy)

  ! dzvy
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*v_iny0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzvy)

  !  FIRST ORDER VZ TERMS  dxvz dyvz dzvz
  ! dxvz
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*v_inz0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxvz)

  ! dyvz
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*v_inz0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dyvz)

  ! dzvz
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*v_inz0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzvz)

  ! DONE ALL FIRST ORDER VX,VY AND VZ TERMS i.e.  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
  if (nv.eq..false.) then
  ! FIRST ORDER BX TERMS ie. dxbx dybx dzbx`
  ! dxbx
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*b_inx0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxbx)

  ! dybx
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*b_inx0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dybx)

  ! dzbx
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*b_inx0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzbx)

  !  FIRST ORDER BY TERMS ie. dxby dyby dzby
  ! dxby
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*b_iny0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxby)

  ! dyby
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*b_iny0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dyby)

  ! dzby
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*b_iny0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzby)

  !! FIRST ORDER BZ TERMS ie. dxbz dybz dzbz
  ! dxbz
  DO i=mype+1,nkx0-1,n_mpi_procs
     temp_smallm(i,:,:)=i_complex*kxgrid(i-1)*b_inz0(i,:,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dxbz)

  ! dybz
  DO j=mype+1,nky0-1,n_mpi_procs
     temp_smallm(:,j,:)=i_complex*kygrid(j-1)*b_inz0(:,j,:)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dybz)

  ! dzbz
  DO k=mype+1,nkz0-1,n_mpi_procs
     temp_smallm(:,:,k)=i_complex*kzgrid(k-1)*b_inz0(:,:,k)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(temp_smallm,temp_small,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  !Add padding for dealiasing
  CALL torealspace(temp_small,dzbz)

  endif

  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL PostIFFTs Time:',dum-sttime

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
  if (nv.eq..false.) then

  store_bx = (bx*dxvx+by*dyvx+bz*dzvx)&
       - (vx*dxbx+vy*dybx+vz*dzbx)&
       - hall*((bx*dxdybz+by*dydybz+bz*dydzbz&
          -bx*dxdzby-by*dydzby-bz*dzdzby)&
       + (dybz*dxbx-dzby*dxbx-dxbz*dybx&
          +dzbx*dybx+dxby*dzbx-dybx*dzbx))

  store_by = (bx*dxvy+by*dyvy+bz*dzvy)&
       - (vx*dxby+vy*dyby+vz*dzby)&
       - hall*((bx*dxdzbx+by*dydzbx+bz*dzdzbx&
          -bx*dxdxbz-by*dxdybz-bz*dxdzbz)&
       + (dybz*dxby-dzby*dxby-dxbz*dyby&
          +dzbx*dyby+dxby*dzby-dybx*dzby))

  store_bz = (bx*dxvz+by*dyvz+bz*dzvz)&
       - (vx*dxbz+vy*dybz+vz*dzbz)&
       - hall*((bx*dxdxby+by*dxdyby+bz*dxdzby&
          -bx*dxdybx-by*dydybx-bz*dydzbx)&
       + (dybz*dxbz-dzby*dxbz-dxbz*dybz&
          +dzbx*dybz+dxby*dzbz-dybx*dzbz))

  endif

  !! ! EQUATION (15)  
  !! U1= vx*dxvx+vy*dyvx+vz*dzvx 
  !! V2 = bx*dxvx+by*dyvx+bz*dzvx 
  !! W1= (0.5)*dyb^2  ! TO BE DONE 

  !! U2= vx*dxvy+vy*dyvy+vz*dzvy   
  !! V2 = bx*dxvy+by*dyvy+bz*dzvy  
  !! W2= (0.5)*dyb^2  ! TO BE DONE

  !! U3= vx*dxvz+vy*dyvz+vz*dzvz  
  !! V3 = bx*dxvz+by*dyvz+bz*dzvz 
  !! W3= (0.5)*dzb^2  ! TO BE DONE

  !!      eq15x= -U1+V1 -W1       
  !!      eq15y=-U2+V2-W2         
  !!      eq15z=-U3+V3-W3         

   store_vx = -(vx*dxvx+vy*dyvx+vz*dzvx)&
         + (by*dybx+bz*dzbx) - (by*dxby+bz*dxbz)
   store_vy = -(vx*dxvy+vy*dyvy+vz*dzvy)&
         + (bx*dxby+bz*dzby) - (bz*dybz+bx*dybx)
   store_vz = -(vx*dxvz+vy*dyvz+vz*dzvz)&
         + (bx*dxby+by*dybz) - (by*dzby+bx*dzbx)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_bxm,store_bx,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_bym,store_by,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_bzm,store_bz,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_vxm,store_vx,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_vym,store_vy,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  CALL MPI_ALLREDUCE(store_vzm,store_vz,nx0_big*ny0_big*local_N&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)

  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL PostNLEq Time:',dum-sttime
  IF(verbose.and.(mype.eq.0)) print *, 'NL BVx Store Values', maxval(abs(store_bx)),maxval(abs(store_vx))

if (nv.eq..false.) then
  temp_bigbx = fft_spec2(store_bx)
  temp_bigby = fft_spec2(store_by) 
  temp_bigbz = fft_spec2(store_bz) 
endif

  temp_bigvx = fft_spec2(store_vx)
  temp_bigvy = fft_spec2(store_vy)
  temp_bigvz = fft_spec2(store_vz)

  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL PostFFFTs Time:',dum-sttime

  !Now fill in appropriate rhs elements

  if (nv.eq..false.) then
  DO i=mype,nkx0-1,n_mpi_procs
        rhs_out_nlbm(i,:,:,0)=temp_bigbx(i,0:nky0-1,0:nkz0-1)
        rhs_out_nlbm(i,:,:,1)=temp_bigby(i,0:nky0-1,0:nkz0-1)
        rhs_out_nlbm(i,:,:,2)=temp_bigbz(i,0:nky0-1,0:nkz0-1)
  END DO!i loop                                                 
  endif

  DO i=mype,nkx0-1,n_mpi_procs
        rhs_out_nlvm(i,:,:,0)=temp_bigvx(i,0:nky0-1,0:nkz0-1)
        rhs_out_nlvm(i,:,:,1)=temp_bigvy(i,0:nky0-1,0:nkz0-1)
        rhs_out_nlvm(i,:,:,2)=temp_bigvz(i,0:nky0-1,0:nkz0-1)
  END DO!i loop

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(rhs_out_nlvm,rhs_out_nlv,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (nv.eq..false.) CALL MPI_ALLREDUCE(rhs_out_nlbm,rhs_out_nlb,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  rhs_out_v = rhs_out_v + rhs_out_nlv
  if (nv.eq..false.) rhs_out_b = rhs_out_b + rhs_out_nlb

  if (verbose.and.(mype.eq.0)) print *, 'rhs out v nl found'

  if (calc_dt) CALL next_dt(myndt)
  if (calc_dt) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (calc_dt) CALL MPI_ALLREDUCE(myndt,ndt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
  if (.not.(calc_dt)) ndt = dt_max
  if (verbose.and.(mype.eq.0)) print *, 'next dt calculated ',ndt

  IF ((plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN
     rsnls = 1
     if (mod(counter,4).eq.0) rsnls = 0

     ! b.grad v
     bdv(:,:,:,0) = rsnls*bdv(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxvx+by*dyvx+bz*dzvx)
     bdv(:,:,:,1) = rsnls*bdv(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxvy+by*dyvy+bz*dzvy)
     bdv(:,:,:,2) = rsnls*bdv(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxvz+by*dyvz+bz*dzvz)

     ! v.grad b
     vdb(:,:,:,0) = rsnls*vdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxbx+vy*dybx+vz*dzbx)
     vdb(:,:,:,1) = rsnls*vdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxby+vy*dyby+vz*dzby)
     vdb(:,:,:,2) = rsnls*vdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxbz+vy*dybz+vz*dzbz)

     ! b. grad curl b
     bdcb(:,:,:,0) = rsnls*bdcb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)
     bdcb(:,:,:,1) = rsnls*bdcb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)
     bdcb(:,:,:,2) = rsnls*bdcb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)

     ! curl b . grad b
     cbdb(:,:,:,0) = rsnls*cbdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx)
     cbdb(:,:,:,1) = rsnls*cbdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby)
     cbdb(:,:,:,2) = rsnls*cbdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz)

     ! v . grad v
     vdv(:,:,:,0) = rsnls*vdv(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxvx+vy*dyvx+vz*dzvx)
     vdv(:,:,:,1) = rsnls*vdv(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxvy+vy*dyvy+vz*dzvy)
     vdv(:,:,:,2) = rsnls*vdv(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(vx*dxvz+vy*dyvz+vz*dzvz)

     ! b . grad b
     bdb(:,:,:,0) = rsnls*bdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxbx+by*dybx+bz*dzbx)
     bdb(:,:,:,1) = rsnls*bdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxby+by*dyby+bz*dzby)
     bdb(:,:,:,2) = rsnls*bdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxbz+by*dybz+bz*dzbz)
        
     ! 0.5 grad b^2
     db2(:,:,:,0) = rsnls*db2(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dxbx+by*dxby+bz*dxbz)
     db2(:,:,:,1) = rsnls*db2(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dybx+by*dyby+bz*dybz)
     db2(:,:,:,2) = rsnls*db2(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 1.5)) * fft_spec2(bx*dzbx+by*dzby+bz*dzbz)

     IF(mod(counter,4).eq.3) THEN
        ! b . grad v
        if(mod(mype-1,n_mpi_procs).eq.0) WRITE(bdvio) bdv(:,:,:,0)
        if(mod(mype-1,n_mpi_procs).eq.0) WRITE(bdvio) bdv(:,:,:,1)
        if(mod(mype-1,n_mpi_procs).eq.0) WRITE(bdvio) bdv(:,:,:,2)
        ! v . grad b
        if(mod(mype-2,n_mpi_procs).eq.0) WRITE(vdbio) vdb(:,:,:,0)
        if(mod(mype-2,n_mpi_procs).eq.0) WRITE(vdbio) vdb(:,:,:,1)
        if(mod(mype-2,n_mpi_procs).eq.0) WRITE(vdbio) vdb(:,:,:,2)
        ! b . grad curl b
        if(mod(mype-3,n_mpi_procs).eq.0) WRITE(bdcbio) bdcb(:,:,:,0)
        if(mod(mype-3,n_mpi_procs).eq.0) WRITE(bdcbio) bdcb(:,:,:,1)
        if(mod(mype-3,n_mpi_procs).eq.0) WRITE(bdcbio) bdcb(:,:,:,2)
        ! curl b . grad b
        if(mod(mype-4,n_mpi_procs).eq.0) WRITE(cbdbio) cbdb(:,:,:,0)
        if(mod(mype-4,n_mpi_procs).eq.0) WRITE(cbdbio) cbdb(:,:,:,1)
        if(mod(mype-4,n_mpi_procs).eq.0) WRITE(cbdbio) cbdb(:,:,:,2)
        ! v . grad v
        if(mod(mype-5,n_mpi_procs).eq.0) WRITE(vdvio) vdv(:,:,:,0)
        if(mod(mype-5,n_mpi_procs).eq.0) WRITE(vdvio) vdv(:,:,:,1)
        if(mod(mype-5,n_mpi_procs).eq.0) WRITE(vdvio) vdv(:,:,:,2)
        ! b . grad b
        if(mod(mype-6,n_mpi_procs).eq.0) WRITE(bdbio) bdb(:,:,:,0)
        if(mod(mype-6,n_mpi_procs).eq.0) WRITE(bdbio) bdb(:,:,:,1)
        if(mod(mype-6,n_mpi_procs).eq.0) WRITE(bdbio) bdb(:,:,:,2)
        ! 0.5 grad b^2
        if(mod(mype-7,n_mpi_procs).eq.0) WRITE(db2io) db2(:,:,:,0)
        if(mod(mype-7,n_mpi_procs).eq.0) WRITE(db2io) db2(:,:,:,1)
        if(mod(mype-7,n_mpi_procs).eq.0) WRITE(db2io) db2(:,:,:,2)
    ENDIF
  ENDIF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL PostWrite Time:',dum-sttime

  if (calc_dt) CALL next_dt(ndt)
  if (.not.(calc_dt)) ndt = dt_max
  if (verbose.and.(mype.eq.0)) print *, 'next dt calculated ',ndt
  IF (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  IF(verbose.and.(mype.eq.0)) print *, 'RHS NL Total Time:',dum-sttime

END SUBROUTINE get_rhs_nl2

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

 if (.not.nv) then
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
 else
 dtn = courant/(maxval(abs(kxgrid))*maxval(abs(vx)) &
    +maxval(abs(kygrid))*maxval(abs(vy)) &
    + maxval(abs(kzgrid))*maxval(abs(vz)))
 endif

END SUBROUTINE next_dt

FUNCTION fft_spec2(arr_real) result(arr_spec)

implicit none

complex :: arr_real(1:nx0_big,1:ny0_big,1:local_N)
complex :: arr_spec(0:nkx0-1,0:nky0-1,0:nkz0-1)
complex :: arr_specm(0:nkx0-1,0:nky0-1,0:nkz0-1)
integer :: s

  data = cmplx(0.0,0.0)
  arr_specm = cmplx(0.0,0.0)

  data = arr_real
  CALL fftw_mpi_execute_dft(plan_r2c,data,data)
  if (verbose.and.(mype.eq.0)) print *, 'Through ffft'
  if (local_j_offset.lt.nkz0) arr_specm(0:nkx0-1,0:nky0-1,local_j_offset:min(nkz0-1,local_j_offset+local_N-1)) = data(1:nkx0,1:nky0,local_j_offset+1:min(nkz0,local_j_offset+local_N))
  if (verbose.and.(mype.eq.0)) print *, 'Through post assignment'

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(arr_specm,arr_spec,nkx0*nky0*nkz0,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  arr_spec = arr_spec * fft_norm
  if (verbose.and.(mype.eq.0)) print *, 'Through fft_spec2'
  IF(verbose.and.(mype.eq.0)) print *, 'NL FFFT Maxes', maxval(abs(arr_spec)),maxloc(abs(arr_spec))

END FUNCTION fft_spec2

SUBROUTINE finalize_fourier

implicit none

DEALLOCATE(temp_small)
DEALLOCATE(temp_smallm)
if (verbose.and.(mype.eq.0)) print *, 'ts deallocated'

DEALLOCATE(temp_bigx)
DEALLOCATE(temp_bigbx)
DEALLOCATE(temp_bigvx)
DEALLOCATE(temp_bigm)
if (verbose.and.(mype.eq.0)) print *, 'tbx deallocated'
DEALLOCATE(temp_bigy)
DEALLOCATE(temp_bigby)
DEALLOCATE(temp_bigvy)
if (verbose.and.(mype.eq.0)) print *, 'tby deallocated'
DEALLOCATE(temp_bigz)
DEALLOCATE(temp_bigbz)
DEALLOCATE(temp_bigvz)
if (verbose.and.(mype.eq.0)) print *, 'tbz deallocated'

DEALLOCATE(store_bxm)
DEALLOCATE(store_bym)
DEALLOCATE(store_bzm)
DEALLOCATE(store_vxm)
DEALLOCATE(store_vym)
DEALLOCATE(store_vzm)
DEALLOCATE(store_bx)
DEALLOCATE(store_vx)
if (verbose.and.(mype.eq.0)) print *, 'stx deallocated'
DEALLOCATE(store_by)
DEALLOCATE(store_vy)
if (verbose.and.(mype.eq.0)) print *, 'sty deallocated'
DEALLOCATE(store_bz)
DEALLOCATE(store_vz)
if (verbose.and.(mype.eq.0)) print *, 'stz deallocated'

! All b arrays 
DEALLOCATE(bx)
if (verbose.and.(mype.eq.0)) print *, 'bx deallocated'
DEALLOCATE(by)
if (verbose.and.(mype.eq.0)) print *, 'by deallocated'
DEALLOCATE(bz)
if (verbose.and.(mype.eq.0)) print *, 'first third deallocated'

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
  ! all  second o) DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  DEALLOCATE(dxdxvx)
  DEALLOCATE(dxdyvx)
  DEALLOCATE(dxdzvx)
  DEALLOCATE(dydyvx)
  DEALLOCATE(dydzvx)
  DEALLOCATE(dzdzvx)
  ! all  second o)
  DEALLOCATE(dxdxvy)
  DEALLOCATE(dxdyvy)
  DEALLOCATE(dxdzvy)
  DEALLOCATE(dydyvy)
  DEALLOCATE(dydzvy)
  DEALLOCATE(dzdzvy)
  ! all  second o)
  DEALLOCATE(dxdxvz)
  DEALLOCATE(dxdyvz)
  DEALLOCATE(dxdzvz)
  DEALLOCATE(dydyvz)
  DEALLOCATE(dydzvz)
  DEALLOCATE(dzdzvz)
  
  if (plot_nls) then
  DEALLOCATE(bdv)
  DEALLOCATE(vdb)
  DEALLOCATE(bdcb)
  DEALLOCATE(cbdb)
  DEALLOCATE(vdv)
  DEALLOCATE(bdb)
  DEALLOCATE(db2)
  endif
  
  DEALLOCATE(output)

  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

  DEALLOCATE(data)
  DEALLOCATE(datai)

  if (verbose.and.(mype.eq.0)) print *, 'deallocated nl code'
  call fftw_destroy_plan(plan_r2c)
  call fftw_destroy_plan(plan_c2r)

if (verbose.and.(mype.eq.0)) print *, "Destroyed plans"
call fftw_cleanup()
if (verbose.and.(mype.eq.0)) print *, 'cleaned up'
call fftw_free(cdata)
call fftw_free(cdatai)

END SUBROUTINE finalize_fourier

SUBROUTINE torealspace(temp_small,output)

  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: temp_small(1:nkx0,1:nky0,1:nkz0)
  COMPLEX(C_DOUBLE_COMPLEX), intent(out) :: output(1:nx0_big,1:ny0_big,1:nz0_big)
  INTEGER(C_INTPTR_T) :: xcind,ycind,zcind

 IF (verbose.and.(mype.eq.0)) print *, maxval(abs(temp_small))

  datai = cmplx(0.0,0.0)
  output = cmplx(0.0,0.0)

  DO k = 1,min(local_Ni,nkz0),n_mpi_procs
    DO j = 1,nky0
      DO i = 1,nkx0
        if (k.le.nkz0) datai(i,j,k) = temp_small(i,j,k+local_j_offseti)
        if (dealias_type.eq.3) then 
          if ((k.le.2*nkz0).and.(k.gt.nkz0+1)) then 
            datai(i,j,k) = conjg(temp_small(1+mod(nkx0+1-i,nkx0),1+mod(nky0+1-j,nky0),1+mod(2*nkz0+1-(k+local_j_offseti),2*nkz0)))
          endif
        endif
      ENDDO
    ENDDO
  ENDDO

 IF (verbose.and.(mype.eq.0)) print *, maxval(abs(datai))

!  if (verbose.and.(mype.eq.0)) print *,'Through \"temp big\"'

 CALL fftw_mpi_execute_dft(plan_c2r,datai,datai)
 output = datai
 temp_small = cmplx(0.0,0.0)

 IF (verbose.and.(mype.eq.0)) print *, maxval(abs(datai)),maxval(abs(output))

END SUBROUTINE torealspace

END MODULE nonlinearity

