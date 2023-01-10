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

  PUBLIC :: initialize_fourier,finalize_fourier,get_rhs_nl,&
            get_k_indices,get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2, get_rhs_nl2, get_rhs_nl_convolution
  
  REAL, PUBLIC :: ve_max(2)

  PRIVATE

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: b_inx0, b_iny0,  b_inz0, v_inx0, v_iny0,  v_inz0

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) ::  bmag_in, bmagk, ekd
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,temp_big, temp_bigx, temp_bigy, temp_bigz,mytemp_small

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) ::  store_x, store_y, store_z,store
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) ::  bmag, dxbmag, dybmag, dzbmag,bmag_inbig


  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz, dxbx, dybx,dzbx, dxby,dyby,dzby, dxbz,dybz,dzbz
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz, dxvx, dyvx,dzvx, dxvy,dyvy,dzvy, dxvz,dyvz,dzvz

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbx, dxdybx, dxdzbx, dydybx, dydzbx, dzdzbx
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxby, dxdyby, dxdzby, dydyby, dydzby, dzdzby
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbz, dxdybz, dxdzbz, dydybz, dydzbz, dzdzbz
    
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvx, dxdyvx, dxdzvx, dydyvx, dydzvx, dzdzvx
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvy, dxdyvy, dxdzvy, dydyvy, dydzvy, dzdzvy
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvz, dxdyvz, dxdzvz, dydyvz, dydzvz, dzdzvz

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mybx,myby,mybz, mydxbx, mydybx,mydzbx, mydxby,mydyby,mydzby, mydxbz,mydybz,mydzbz
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: myvx,myvy,myvz, mydxvx, mydyvx,mydzvx, mydxvy,mydyvy,mydzvy, mydxvz,mydyvz,mydzvz

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mydxdxbx, mydxdybx, mydxdzbx, mydydybx, mydydzbx, mydzdzbx
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mydxdxby, mydxdyby, mydxdzby, mydydyby, mydydzby, mydzdzby
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mydxdxbz, mydxdybz, mydxdzbz, mydydybz, mydydzbz, mydzdzbz

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: rhs_out_nlb,rhs_out_nlv,myrhs_out_nlb,myrhs_out_nlv
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: bdv,vdb,bdcb,cbdb,vdv,bdb,db2
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: mybdv,myvdb,mybdcb,mycbdb,myvdv,mybdb,mydb2

  !For fft's

  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:):: g_kbig
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:):: g_rbig
  COMPLEX, ALLOCATABLE, DIMENSION(:,:):: g_kbig_2d
  COMPLEX, ALLOCATABLE, DIMENSION(:,:):: g_rbig_2d
  INTEGER :: nx0_big,ny0_big,nz0_big
  INTEGER(kind=8), allocatable :: plan_r2c,plan_c2r
  INTEGER(kind=8) :: plan_kz2z,plan_z2kz
  INTEGER(kind=8) :: plan_ky2y,plan_y2ky
  INTEGER(kind=8) :: plan_kx2x,plan_x2kx
  REAL :: fft_norm  !normalization factor for inverse fft
  INTEGER :: zpad,counter
 
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
  if ((dealias_type.eq.3).or.(dealias_type.eq.4).or.(dealias_type.eq.20)) zpad = dealias_type
  nx0_big=zpad*nkx0/2
  ny0_big=zpad*nky0/2
  nz0_big=zpad*nkz0/2
  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))
  counter = 0
  ALLOCATE(plan_r2c)
  ALLOCATE(plan_c2r)

  ALLOCATE(store(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_big(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
 
  !WRITE(*,*) "making plans"
  CALL dfftw_plan_dft_3d(plan_c2r,nx0_big,ny0_big,nz0_big,&
                             temp_big,store,FFTW_ESTIMATE,FFTW_BACKWARD)
  CALL dfftw_plan_dft_3d(plan_r2c,nx0_big,ny0_big,nz0_big,&
                             store,temp_big,FFTW_ESTIMATE,FFTW_FORWARD)
  
  lky_big=ny0_big-hky_ind !Index of minimum (most negative) FILLED ky value for big arrays
  lkz_big=nz0_big-hkz_ind !Index of minimum (most negative) FILLED kz value for big arrays 

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big,nz0_big",nx0_big,ny0_big,nz0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big
  IF(mype==0) WRITE(*,*) "hkz_ind,lkz_ind",hkz_ind,lkz_ind
  IF(mype==0) WRITE(*,*) "lkz_big",lkz_big

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp,treal)
  !CALL dfftw_execute_dft(plan_r2c,treal,tcomp)
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
  COMPLEX :: rhs_out_b2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX :: rhs_out_v2(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt
  
  !IF(mype==0) WRITE(*,*) "In get_rhs_nl"
  !IF(mype==0) WRITE(*,*) "Version is: ",rhs_nl_version
  IF ((rhs_nl_version.eq.1).or.(rhs_nl_version.eq.12)) THEN
    !IF(mype==0) WRITE(*,*) "version was 1"
    IF (dealias_type.eq.3) CALL get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
    IF (dealias_type.eq.4) CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
    IF (dealias_type.eq.1) CALL get_rhs_nlps(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
!    if ((verbose).and.(mod(itime,100).eq.1)) CALL wherenezero(rhs_out_b(0:nkx0-1,0:nky0-1,0:nkz0-1,0)-b_in(0:nkx0-1,0:nky0-1,0:nkz0-1,0))
!    if (verbose) print *, 'Checked nonzero convolution terms'
!    if ((verbose).and.(mod(itime,100).eq.1)) CALL wherebvnotreal(rhs_out_b,rhs_out_v)
!    if (verbose) print *, 'Checked b,v reality condition'
    if (dealias_type.eq.20) THEN
      CALL finalize_fourier
      dealias_type = 4
      CALL initialize_fourier_ae_mu0
      rhs_out_b2 = rhs_out_b
      rhs_out_v2 = rhs_out_v
      CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
      rhs_out_b = -1.0 * rhs_out_b
      rhs_out_v = -1.0 * rhs_out_v
      CALL finalize_fourier
      dealias_type = 1
      CALL initialize_fourier_ae_mu0
      CALL get_rhs_nlps(b_in,v_in,rhs_out_b,rhs_out_v,ndt)
      print *, 'B Phase Shift - Padding Comparison'
      CALL wherenezero(rhs_out_b-rhs_out_b2)
      print *, 'V Phase Shift - Padding Comparison'
      CALL wherenezero(rhs_out_v-rhs_out_v2)
      rhs_out_v = rhs_out_v2
      rhs_out_b = rhs_out_b2
      dealias_type = 20
    endif
    counter = counter + 1
!  ELSE IF(rhs_nl_version==2) THEN
!    CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==3) THEN
!    CALL get_rhs_nl3(b_in,v_in,rhs_out_b,rhs_out_v)
!  ELSE IF(rhs_nl_version==4) THEN
!    CALL get_rhs_nl4(b_in,v_in,rhs_out_b,rhs_out_v)

  END IF
 
END SUBROUTINE get_rhs_nl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  get_rhs_nlps                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_rhs_nlps(b_in, v_in, rhs_out_b, rhs_out_v,ndt)
  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

  REAL :: ndt
  INTEGER :: delta_x, delta_y, delta_z

  rhs_out_b = 8.0 * rhs_out_b
  rhs_out_v = 8.0 * rhs_out_v
  DO delta_x = 0,1
    DO delta_y = 0,1
      DO delta_z = 0,1
        CALL get_rhs_nl3(b_in,v_in,delta_x * pi/(nx0_big),delta_y * pi/(ny0_big), delta_z * pi/(nz0_big),rhs_out_b,rhs_out_v,ndt)
      ENDDO
    ENDDO
  ENDDO
  rhs_out_b = rhs_out_b / 8.0
  rhs_out_v = rhs_out_v / 8.0

END SUBROUTINE get_rhs_nlps

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
   INTEGER :: i,j,l,k,h, ierr

  IF(np_spec.gt.1) STOP "get_rhs_nl1 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl1 not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_nl1 only suitable for np_kz=1"
  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  
  ALLOCATE(store_x(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_y(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_z(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

! All b arrays 
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

!  if (rhs_nl_version==1) then
!    ALLOCATE(bmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))  
!    ALLOCATE(bmag_in(0:nkx0-1,0:nky0-1,0:nkz0-1))
!    ALLOCATE(bmag_inbig(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
!    ALLOCATE(dxbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
!    ALLOCATE(dybmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
!    ALLOCATE(dzbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
!  endif
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
    
  ! all  second order vx arrays i.e. DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  ALLOCATE(dxdxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order vy arrays
  ALLOCATE(dxdxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order bz arrays
  ALLOCATE(dxdxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

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


! bmag should be calculated in real space
!  DO i=0,nkx0-1
!    DO j=0,nky0-1
!        DO k=0,nkz0-1
!            bmag_in(i,j,k) = b_inx0(i,j,k)*b_inx0(i,j,k)+ b_iny0(i,j,k)*b_iny0(i,j,k)+ b_inz0(i,j,k)*b_inz0(i,j,k)
!        END DO
!    END DO
!  END DO

  
   !   SECOND ORDER  VX TERMS DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  ! second order v terms are not needed 
  if (.false.) then 

    !dxdxvx
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inx0(i,:,:) ! there is  two i's in the 
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxvx = store
   
   !dxdyvx
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inx0(i,j,:)
    END DO
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdyvx = store    

  !dxdzvx
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inx0(i,:,k)
    END DO
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzvx = store

  ! DYDYVX
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inx0(:,j,:)  ! towo y grid
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydyvx = store    

     ! DYDzVX
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inx0(:,j,k)  
    END DO
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzvx = store  
  
     ! DzDzVX
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kygrid(k)*i_complex*kzgrid(k)*v_inx0(:,:,k)  !two kz grid
  END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzvx = store
   
   ! finished SECOND ORDER  VX TERMS
   
! SECOND ORDER  VY TERMS DXDXVY,DXDYVY,DXDZVY, DYDYVY,DYDZVY, DZDZVY 

  !dxdxvy
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_iny0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxvy = store
   
   !dxdyvy
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_iny0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdyvy = store
    
     !dxdzvy
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_iny0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzvy = store

   ! DYDYVy
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_iny0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydyvy = store
    
     ! DYDzVY
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_iny0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzvy = store
    
     ! DzDzVY
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kygrid(k)*i_complex*kzgrid(k)*v_iny0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzvy = store
    
! finished SECOND ORDER VY TERMS

! SECOND ORDER  VZ TERMS DXDXVZ,DXDYVZ,DXDZVZ, DYDYVz,DYDZVz, DZDZVz
!dxdxvz
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inz0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxvz = store
   
   !dxdyvz
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inz0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdyvz = store
    
     !dxdzvz
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inz0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzvz = store

   ! DYDYVz
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inz0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydyvz = store
    
     ! DYDzVz
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inz0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzvz = store
    
     ! DzDzVz
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_inz0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzvz = store
    
! finished SECOND ORDER VZ TERMS
  endif

! START SECOND ORDER b terms  

 ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
    !dxdxbx
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inx0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxbz = store
   
   !dxdybx
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inx0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdybx = store
    
     !dxdzbx
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inx0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzbx = store
    
   ! DYDYbX
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inx0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydybx = store
    
     ! DYDzbX
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inx0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzbx = store

     ! DzDzbX
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inx0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzbx = store
   
   ! FINISHED SECOND ORDER  BX TERMS 
   
! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY 

!dxdxby
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_iny0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxby = store
   
   !dxdyvy
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_iny0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdyby = store
    
     !dxdzby
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_iny0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzby = store

   ! DYDYby
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_iny0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydyby = store
    
     ! DYDzbY
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_iny0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzby = store
    
     ! DzDzbY
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_iny0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzby = store
    
! FINISHED SECOND ORDER bY TERMS

! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
!dxdxbz
  DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inz0(i,:,:) ! there is  two i's in the 
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdxbz = store
   
   !dxdybz
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inz0(i,j,:)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdybz = store
    
     !dxdzbz
  DO i=0,nkx0-1
    DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inz0(i,:,k)
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxdzbz = store

   ! DYDYbz
  DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inz0(:,j,:)  ! towo y grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydybz = store
    
     ! DYDzbz
  DO j=0,nky0-1
    DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inz0(:,j,k)  
    END DO
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dydzbz = store
    
     ! DzDzbz
  DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inz0(:,:,k)  !two kz grid
    END DO
        !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzdzbz = store
    
!finished END SECOND ORDER BZ TERMS

!    completed  ALL SECOND ORDER B TERMS  

! TERMS BX BY BZ

    !bx
    DO i=0,nkx0-1
        temp_small(i,:,:)=b_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    bx = store

    !by
    DO j=0,nky0-1
        temp_small(:,j,:)=b_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    by = store

    !bz
    DO k=0,nkz0-1
        temp_small(:,:,k)=b_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    bz = store

!!! bmag terms
if (rhs_nl_version == 1) then
!! Calculate b^2(x) in real space.
  DO i=0,nx0_big-1
    DO j=0,ny0_big-1
        DO k=0,nz0_big-1
            bmag(i,j,k) = bx(i,j,k)*bx(i,j,k)+ by(i,j,k)*by(i,j,k)+ bz(i,j,k)*bz(i,j,k)
        END DO
    END DO
  END DO
!!! Now forward transfrom to get b^2(k) stored in temp_big
temp_big=cmplx(0.0,0.0)
store = bmag
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
!!! Now cut down to dealias and store result in bmagk which is b^2(k)
  DO i=0,nkx0-1
  !First x component
            bmagk(i,0:hky_ind,0:hkz_ind)= temp_big(i,0:hky_ind,0:hkz_ind)*fft_norm !kz positive, ky positive
            bmagk(i,0:hky_ind,lkz_ind:nkz0-1)=temp_big(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm !kz negative, ky positive
            bmagk(i,lky_ind:nky0-1,lkz_ind:nkz0-1)= temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm !kz negative, ky negative
            bmagk(i,lky_ind:nky0-1,0:hkz_ind)= temp_big(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm !kz positive, ky negative
  END DO!i loop
!!! Now that we have bmagk, need to do IFFT(ikx b^2) IFFT(iky b^2) IFFT(ikz b^2) to get grad b^2/2 in real space

    !bmag_inbig = bx*bx + by*by + bz*bz
    !CALL dfftw_execute_dft(plan_r2c,bmag_inbig,bmag)
    !do i = 0,nkx0-1
    !        bmag_in(i,0:hky_ind,0:hkz_ind)=bmag(i,0:hky_ind,0:hkz_ind)*fft_norm   !kz positive, ky positive 
    !        bmag_in(i,0:hky_ind,lkz_ind:nkz0-1)=bmag(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm !kz negative, ky positive 
    !        bmag_in(i,lky_ind:nky0-1,lkz_ind:nkz0-1) = bmag(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm ! kz negative, ky negative
    !        bmag_in(i,lky_ind:nky0-1,0:hkz_ind)=bmag(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm !kz positive, ky negative
    !end do
              
    !dxbmag                                                                                                                                                       
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*bmagk(i,:,:)
    END DO
    !Add padding for dealiasing                                                                                                                                                 
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive                                                          
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive   
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative                                   
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative                         
    END DO!k loop                                                                      
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxbmag = store

    !dybmag                                                                                                                                                                      
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*bmagk(:,j,:)
    END DO
    !Add padding for dealiasing                                                                                                                                                  
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive                                                            
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive                                                
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative                                     
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative                                                 
    END DO!k loop                                       
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dybmag = store

    !dzbmag                                                                                                                                                                      
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*bmagk(:,:,k)
    END DO
    !Add padding for dealiasing                                                                                                                                                 
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive                                   
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive                         
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative            
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative                                                 
    END DO!k loop                                                                                                                                                
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzbmag = store

    if (verbose.and.(mype.eq.0)) print *, 'Bmag gradients dealiased'
endif





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
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
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
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
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
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    vz = store

!  FIRST ORDER VX TERMS DXVX , DYVX,  DZVX

    ! dxvx
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxvx = store

    ! dyvx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inx0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dyvx = store

    ! dzvx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inx0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzvx = store
    
   !  FIRST ORDER VY TERMS  dxvy dyvy dzvz, 

    ! dxvy
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_iny0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxvy = store
 
    ! dyvy
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dyvy = store

    ! dzvy
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_iny0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzvy = store
    
    !  FIRST ORDER VZ TERMS  dxvz dyvz dzvz
    ! dxvz
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*v_inz0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxvz = store
 
    ! dyvz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dyvz = store

    ! dzvz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzvz = store
    
    ! DONE ALL FIRST ORDER VX,VY AND VZ TERMS i.e.  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
    
          ! FIRST ORDER BX TERMS ie. dxbx dybx dzbx` 
     ! dxbx
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_inx0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxbx = store

    ! dybx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_inx0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dybx = store

    ! dzbx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inx0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzbx = store
        
    !  FIRST ORDER BY TERMS ie. dxby dyby dzby
     ! dxby
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_iny0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxby = store
 
    ! dyby
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_iny0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dyby = store

    ! dzby
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_iny0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzby = store
    
        !! FIRST ORDER BZ TERMS ie. dxbz dybz dzbz
! dxbz
    DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*b_inz0(i,:,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dxbz = store
 
    ! dybz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dybz = store

    ! dzbz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inz0(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft(plan_c2r,temp_big,store)
    dzbz = store

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

IF ((plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN 
! b.grad v 
WRITE(bdvio) fft_spec(bx*dxvx+by*dyvx+bz*dzvx)
WRITE(bdvio) fft_spec(bx*dxvy+by*dyvy+bz*dzvy)
WRITE(bdvio) fft_spec(bx*dxvz+by*dyvz+bz*dzvz)

! v.grad b
WRITE(vdbio) fft_spec(vx*dxbx+vy*dybx+vz*dzbx)
WRITE(vdbio) fft_spec(vx*dxby+vy*dyby+vz*dzby)
WRITE(vdbio) fft_spec(vx*dxbz+vy*dybz+vz*dzbz)

! b. grad curl b
WRITE(bdcbio) fft_spec(bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)
WRITE(bdcbio) fft_spec(bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)
WRITE(bdcbio) fft_spec(bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)

! curl b . grad b
WRITE(cbdbio) fft_spec(dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx)
WRITE(cbdbio) fft_spec(dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby)
WRITE(cbdbio) fft_spec(dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz)
ENDIF

!inverse FFT to get back to Fourier
store = store_x
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigx = temp_big

store = store_y
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigy = temp_big

store = store_z
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigz = temp_big

!Now fill in appropriate rhs elements
  DO i=0,nkx0-1
  !First x component
            rhs_out_b(i,0:hky_ind,0:hkz_ind,0)=rhs_out_b(i,0:hky_ind,0:hkz_ind,0)+&           !kz positive, ky positive
                                       temp_bigx(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,0)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,0)+& !kz negative, ky positive
                                       temp_bigx(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)+& !kz negative, ky negative
                                       temp_bigx(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,0)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,0)+& !kz positive, ky negative
                                       temp_bigx(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !y component
            rhs_out_b(i,0:hky_ind,0:hkz_ind,1)=rhs_out_b(i,0:hky_ind,0:hkz_ind,1)+&           !kz positive, ky positive
                                       temp_bigy(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,1)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,1)+& !kz negative, ky positive
                                       temp_bigy(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)+& !kz negative, ky negative
                                       temp_bigy(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,1)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,1)+& !kz positive, ky negative
                                       temp_bigy(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !z component
            rhs_out_b(i,0:hky_ind,0:hkz_ind,2)=rhs_out_b(i,0:hky_ind,0:hkz_ind,2)+&           !kz positive, ky positive
                                       temp_bigz(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,2)=rhs_out_b(i,0:hky_ind,lkz_ind:nkz0-1,2)+& !kz negative, ky positive
                                       temp_bigz(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)=rhs_out_b(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)+& !kz negative, ky negative
                                       temp_bigz(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,2)=rhs_out_b(i,lky_ind:nky0-1,0:hkz_ind,2)+& !kz positive, ky negative
                                       temp_bigz(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  END DO!i loop

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

IF ((plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN 
! v . grad v
WRITE(vdvio) fft_spec(vx*dxvx+vy*dyvx+vz*dzvx) 
WRITE(vdvio) fft_spec(vx*dxvy+vy*dyvy+vz*dzvy)
WRITE(vdvio) fft_spec(vx*dxvz+vy*dyvz+vz*dzvz)

! b . grad b
WRITE(bdbio) fft_spec(bx*dxbx+by*dybx+bz*dzbx)
WRITE(bdbio) fft_spec(bx*dxby+by*dyby+bz*dzby)
WRITE(bdbio) fft_spec(bx*dxbz+by*dybz+bz*dzbz)

! 0.5 grad b^2
WRITE(db2io) fft_spec(0.5*dxbmag)
WRITE(db2io) fft_spec(0.5*dybmag)
WRITE(db2io) fft_spec(0.5*dzbmag)

if (verbose.and.(mype.eq.0)) print *, 'v1 nl equation stored'
endif
ENDIF

if (rhs_nl_version == 12) then
store_x = -(vx*dxvx+vy*dyvx+vz*dzvx) + (by*dybx+bz*dzbx) - (by*dxby+bz*dxbz)
store_y = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+bz*dzby) - (bz*dybz+bx*dybx)
store_z = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz) - (by*dzby+bx*dzbx)

IF ((plot_nls.and.(mod(itime,istep_energy).eq.0))) THEN 
! v . grad v
WRITE(vdvio) fft_spec(vx*dxvx+vy*dyvx+vz*dzvx)
WRITE(vdvio) fft_spec(vx*dxvy+vy*dyvy+vz*dzvy)
WRITE(vdvio) fft_spec(vx*dxvz+vy*dyvz+vz*dzvz)

! b . grad b
WRITE(bdbio) fft_spec(bx*dxbx+by*dybx+bz*dzbx)
WRITE(bdbio) fft_spec(bx*dxby+by*dyby+bz*dzby)
WRITE(bdbio) fft_spec(bx*dxbz+by*dybz+bz*dzbz)

! 0.5 grad b^2
WRITE(db2io) fft_spec(bx*dxbx+by*dxby+bz*dxbz)
WRITE(db2io) fft_spec(bx*dybx+by*dyby+bz*dybz)
WRITE(db2io) fft_spec(bx*dzbx+by*dzby+bz*dzbz)

if (verbose.and.(mype.eq.0)) print *, 'v12 nl equation stored'
endif
ENDIF

store = store_x
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigx = temp_big

store = store_y
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigy = temp_big

store = store_z
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temp_bigz = temp_big

!Now fill in appropriate rhs elements
  DO i=0,nkx0-1
  !First x component
            rhs_out_v(i,0:hky_ind,0:hkz_ind,0)=rhs_out_v(i,0:hky_ind,0:hkz_ind,0)+&           !kz positive, ky positive
                                       temp_bigx(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,0)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,0)+& !kz negative, ky positive
                                       temp_bigx(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,0)+& !kz negative, ky negative
                                       temp_bigx(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,0)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,0)+& !kz positive, ky negative
                                       temp_bigx(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !y component
            rhs_out_v(i,0:hky_ind,0:hkz_ind,1)=rhs_out_v(i,0:hky_ind,0:hkz_ind,1)+&           !kz positive, ky positive
                                       temp_bigy(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,1)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,1)+& !kz negative, ky positive
                                       temp_bigy(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,1)+& !kz negative, ky negative
                                       temp_bigy(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,1)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,1)+& !kz positive, ky negative
                                       temp_bigy(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  !z component
            rhs_out_v(i,0:hky_ind,0:hkz_ind,2)=rhs_out_v(i,0:hky_ind,0:hkz_ind,2)+&           !kz positive, ky positive
                                       temp_bigz(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,2)=rhs_out_v(i,0:hky_ind,lkz_ind:nkz0-1,2)+& !kz negative, ky positive
                                       temp_bigz(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)=rhs_out_v(i,lky_ind:nky0-1,lkz_ind:nkz0-1,2)+& !kz negative, ky negative
                                       temp_bigz(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,2)=rhs_out_v(i,lky_ind:nky0-1,0:hkz_ind,2)+& !kz positive, ky negative
                                       temp_bigz(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
  END DO!i loop

if (verbose.and.(mype.eq.0)) print *, 'rhs out v nl found'

if (calc_dt) CALL next_dt(ndt)
if (.not.(calc_dt)) ndt = dt_max
if (verbose.and.(mype.eq.0)) print *, 'next dt calculated ',ndt

DEALLOCATE(temp_small)
if (verbose.and.(mype.eq.0)) print *, 'ts deallocated'

DEALLOCATE(temp_bigx)
if (verbose.and.(mype.eq.0)) print *, 'tbx deallocated'
DEALLOCATE(temp_bigy)
if (verbose.and.(mype.eq.0)) print *, 'tby deallocated'
DEALLOCATE(temp_bigz)
if (verbose.and.(mype.eq.0)) print *, 'tbz deallocated'

DEALLOCATE(store_x)
if (verbose.and.(mype.eq.0)) print *, 'stx deallocated'
DEALLOCATE(store_y)
if (verbose.and.(mype.eq.0)) print *, 'sty deallocated'
DEALLOCATE(store_z)
if (verbose.and.(mype.eq.0)) print *, 'stz deallocated'

! All b arrays 
DEALLOCATE(bx)
if (verbose.and.(mype.eq.0)) print *, 'bx deallocated'
DEALLOCATE(by)
if (verbose.and.(mype.eq.0)) print *, 'by deallocated'
DEALLOCATE(bz)
if (verbose.and.(mype.eq.0)) print *, 'first third deallocated'

if (rhs_nl_version==1) then
!DEALLOCATE(bmag_in) 
!if (verbose) print *,'bmagin deallocated'
DEALLOCATE(bmagk)
if (verbose) print *,'bmagk deallocated'
DEALLOCATE(bmag)
if (verbose) print *,'bmag deallocated'
DEALLOCATE(dxbmag)
if (verbose) print *,'dxbmag deallocated'
DEALLOCATE(dybmag)
if (verbose) print *,'dybmag deallocated'
DEALLOCATE(dzbmag)
if (verbose) print *,'dzbmag deallocated'
!DEALLOCATE(bmag_inbig)
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
  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

if (verbose.and.(mype.eq.0)) print *, 'deallocated nl code'

END SUBROUTINE get_rhs_nl1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v,ndt)

 USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: sttime,dum
  REAL :: ndt

  INTEGER :: i,j,l,k,h, ierr

  if (verbose.and.(mype.eq.0)) CALL cpu_time(sttime)

  IF(np_spec.gt.1) STOP "get_rhs_nl1 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl1 not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_nl1 only suitable for np_kz=1"

  ALLOCATE(rhs_out_nlb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(rhs_out_nlv(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(myrhs_out_nlb(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  ALLOCATE(myrhs_out_nlv(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))

  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(store_x(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_y(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_z(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))


  ! All b arrays
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  !  if (rhs_nl_version==1) then
  !    ALLOCATE(bmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(bmag_in(0:nkx0-1,0:nky0-1,0:nkz0-1))
  !    ALLOCATE(bmag_inbig(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dxbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dybmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dzbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !  endif
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

  ! all  second order vx arrays i.e. DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  ALLOCATE(dxdxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order vy arrays
  ALLOCATE(dxdxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order bz arrays
  ALLOCATE(dxdxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ! Ind processors

  ALLOCATE(mytemp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))

  ALLOCATE(mybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(myby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(myvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(myvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(myvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
 
  ALLOCATE(mydxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
 
  ALLOCATE(mydxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
 
  ALLOCATE(mydxbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxdxbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzdzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxdxby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzdzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(mydxdxbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydxdzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydydzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(mydzdzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ! nonlinearities
  IF(mod(counter,4).eq.0) THEN
  ALLOCATE(bdv(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(vdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(bdcb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(cbdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(vdv(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(bdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(db2(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(mybdv(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(myvdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(mybdcb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(mycbdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(myvdv(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(mybdb(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ALLOCATE(mydb2(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1,0:2))
  ENDIF


  ! current b and v arrays

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


  ! bmag should be calculated in real space
  !  DO i=0,nkx0-1
  !    DO j=0,nky0-1
  !        DO k=0,nkz0-1
  !            bmag_in(i,j,k) = b_inx0(i,j,k)*b_inx0(i,j,k)+ b_iny0(i,j,k)*b_iny0(i,j,k)+ b_inz0(i,j,k)*b_inz0(i,j,k)
  !        END DO
  !    END DO
  !  END DO

  if (nv.eq..false.) then
  
  ! second order v terms are not needed
  if (.false.) then

  !   SECOND ORDER  VX TERMS DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  !dxdxvx
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inx0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvx = store

  !dxdyvx
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inx0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvx = store

  !dxdzvx
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inx0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvx = store

  ! DYDYVX
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inx0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvx = store

  ! DYDzVX
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inx0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvx = store

  ! DzDzVX
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_inx0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvx = store

    ! finished SECOND ORDER  VX TERMS

  ! SECOND ORDER  VY TERMS DXDXVY,DXDYVY,DXDZVY, DYDYVY,DYDZVY, DZDZVY

  !dxdxvy
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_iny0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvy = store

  !dxdyvy
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_iny0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvy = store

  !dxdzvy
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_iny0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvy = store

  ! DYDYVy
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_iny0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvy = store

  ! DYDzVY
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_iny0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvy = store

  ! DzDzVY
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_iny0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvy = store

  ! finished SECOND ORDER VY TERMS

  ! SECOND ORDER  VZ TERMS DXDXVZ,DXDYVZ,DXDZVZ, DYDYVz,DYDZVz, DZDZVz
  !dxdxvz
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inz0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvz = store

  !dxdyvz
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inz0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvz = store

  !dxdzvz
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inz0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvz = store

  ! DYDYVz
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inz0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvz = store

  ! DYDzVz
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inz0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvz = store

  ! DzDzVz
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_inz0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvz = store

  ! finished SECOND ORDER VZ TERMS   ! second order v terms are not needed  
  endif

  ! START SECOND ORDER b terms

  ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
  !dxdxbx
  IF (mod(mype-1,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inx0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdxbx = store
  ENDIF

  !dxdybx
  IF (mod(mype-2,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO j=0,nky0-1
        mytemp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inx0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdybx = store
  ENDIF

  !dxdzbx
  IF (mod(mype-3,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        mytemp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inx0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdzbx = store
  ENDIF

  ! DYDYbX
  IF (mod(mype-4,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inx0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydybx = store
  ENDIF

  ! DYDzbX
  IF (mod(mype-5,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     DO k=0,nkz0-1
        mytemp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inx0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydzbx = store
  ENDIF

  ! DzDzbX
  IF (mod(mype-6,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inx0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzdzbx = store
  ENDIF
  ! FINISHED SECOND ORDER  BX TERMS

  ! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY

  !dxdxby
  IF (mod(mype-7,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_iny0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdxby = store
  ENDIF

  !dxdyby
  IF (mod(mype-8,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO j=0,nky0-1
        mytemp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_iny0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdyby = store
  ENDIF

  !dxdzby
  IF (mod(mype-9,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        mytemp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_iny0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdzby = store
  ENDIF

  ! DYDYby
  IF (mod(mype-10,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_iny0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydyby = store
  ENDIF

  ! DYDzbY
  IF (mod(mype-11,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     DO k=0,nkz0-1
        mytemp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_iny0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydzby = store
  ENDIF

  ! DzDzbY
  IF (mod(mype-12,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_iny0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzdzby = store
  ENDIF
  ! FINISHED SECOND ORDER bY TERMS

  ! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
  !dxdxbz

  IF (mod(mype-13,n_mpi_procs).eq.0) THEN 
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inz0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdxbz = store
  ENDIF


  !dxdybz
  IF (mod(mype-14,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO j=0,nky0-1
        mytemp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inz0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdybz = store
  ENDIF

  !dxdzbz
  IF (mod(mype-15,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        mytemp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inz0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxdzbz = store
  ENDIF

  ! DYDYbz
  IF (mod(mype-16,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inz0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydybz = store
  ENDIF

  ! DYDzbz
  IF (mod(mype-17,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     DO k=0,nkz0-1
        mytemp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inz0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydydzbz = store
  ENDIF

  ! DzDzbz
  IF (mod(mype-18,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inz0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzdzbz = store
  ENDIF

  !finished END SECOND ORDER BZ TERMS

  !    completed  ALL SECOND ORDER B TERMS

  ! TERMS BX BY BZ

  !bx
  IF (mod(mype-19,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=b_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mybx = store
  ENDIF

  !by
  IF (mod(mype-20,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=b_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  myby = store
  ENDIF

  !bz
  IF (mod(mype-21,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=b_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mybz = store
  ENDIF

!!! bmag terms
  if (rhs_nl_version == 1) then
     !! Calculate b^2(x) in real space.
     DO i=0,nx0_big-1
        DO j=0,ny0_big-1
           DO k=0,nz0_big-1
              bmag(i,j,k) = bx(i,j,k)*conjg(bx(i,j,k))+ by(i,j,k)*conjg(by(i,j,k))+ bz(i,j,k)*conjg(bz(i,j,k))
           END DO
        END DO
     END DO
!!! Now forward transfrom to get b^2(k) stored in temp_big
     temp_big=cmplx(0.0,0.0)
     store = bmag
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
!!! Now cut down to dealias and store result in bmagk which is b^2(k)
     
     DO i=0,nkx0-1
        bmagk(i,0:nky0-1,0:nkz0-1)= temp_big(i,0:nky0-1,0:nkz0-1)*fft_norm
     END DO!i loop
     
!!! Now that we have bmagk, need to do IFFT(ikx b^2) IFFT(iky b^2) IFFT(ikz b^2) to get grad b^2/2 in real space

     !bmag_inbig = bx*bx + by*by + bz*bz
     !CALL dfftw_execute_dft(plan_r2c,bmag_inbig,bmag)
     !do i = 0,nkx0-1
     !        bmag_in(i,0:hky_ind,0:hkz_ind)=bmag(i,0:hky_ind,0:hkz_ind)*fft_norm   !kz positive, ky positive
     !        bmag_in(i,0:hky_ind,lkz_ind:nkz0-1)=bmag(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm !kz negative, ky positive
     !        bmag_in(i,lky_ind:nky0-1,lkz_ind:nkz0-1) = bmag(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm ! kz negative, ky negative
     !        bmag_in(i,lky_ind:nky0-1,0:hkz_ind)=bmag(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm !kz positive, ky negative
     !end do

     !dxbmag
     DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*bmagk(i,:,:)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dxbmag = store

     !dybmag
     DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*bmagk(:,j,:)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dybmag = store

     !dzbmag
     DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*bmagk(:,:,k)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dzbmag = store

     if (verbose.and.(mype.eq.0)) print *, 'Bmag gradients dealiased'
  endif
  endif
!!! TERMS  vx,vy,vz
  !vx
  IF (mod(mype-22,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=v_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  myvx = store
  ENDIF

  !vy
  IF (mod(mype-23,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
     mytemp_small(:,j,:)=v_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  myvy = store
  ENDIF

  !vz
  IF (mod(mype-24,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=v_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  myvz = store
  ENDIF

  !  FIRST ORDER VX TERMS DXVX , DYVX,  DZVX

  ! dxvx
  IF (mod(mype-25,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*v_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxvx = store
  ENDIF

  ! dyvx
  IF (mod(mype-26,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*v_inx0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydyvx = store
  ENDIF

  ! dzvx
  IF (mod(mype-27,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*v_inx0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzvx = store
  ENDIF

  !  FIRST ORDER VY TERMS  dxvy dyvy dzvz,

  ! dxvy
  IF (mod(mype-28,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*v_iny0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxvy = store
  ENDIF

  ! dyvy
  IF (mod(mype-29,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*v_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydyvy = store
  ENDIF

  ! dzvy
  IF (mod(mype-30,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*v_iny0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzvy = store
  ENDIF

  !  FIRST ORDER VZ TERMS  dxvz dyvz dzvz
  ! dxvz
  IF (mod(mype-31,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*v_inz0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxvz = store
  ENDIF

  ! dyvz
  IF (mod(mype-32,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydyvz = store
  ENDIF

  ! dzvz
  IF (mod(mype-33,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*v_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzvz = store
  ENDIF

  ! DONE ALL FIRST ORDER VX,VY AND VZ TERMS i.e.  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
  if (nv.eq..false.) then
  ! FIRST ORDER BX TERMS ie. dxbx dybx dzbx`
  ! dxbx
  IF (mod(mype-34,n_mpi_procs).eq.0) THEN  
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*b_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxbx = store
  ENDIF

  ! dybx
  IF (mod(mype-35,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*b_inx0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydybx = store
  ENDIF

  ! dzbx
  IF (mod(mype-36,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*b_inx0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzbx = store
  ENDIF

  !  FIRST ORDER BY TERMS ie. dxby dyby dzby
  ! dxby
  IF (mod(mype-37,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*b_iny0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxby = store
  ENDIF

  ! dyby
  IF (mod(mype-38,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*b_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydyby = store
  ENDIF

  ! dzby
  IF (mod(mype-39,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*b_iny0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzby = store
  ENDIF

  !! FIRST ORDER BZ TERMS ie. dxbz dybz dzbz
  ! dxbz
  IF (mod(mype-40,n_mpi_procs).eq.0) THEN
  DO i=0,nkx0-1
     mytemp_small(i,:,:)=i_complex*kxgrid(i)*b_inz0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydxbz = store
  ENDIF

  ! dybz
  IF (mod(mype-41,n_mpi_procs).eq.0) THEN
  DO j=0,nky0-1
     mytemp_small(:,j,:)=i_complex*kygrid(j)*b_inz0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydybz = store
  ENDIF

  ! dzbz
  IF (mod(mype-42,n_mpi_procs).eq.0) THEN
  DO k=0,nkz0-1
     mytemp_small(:,:,k)=i_complex*kzgrid(k)*b_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=mytemp_small(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  mydzbz = store
  ENDIF
  endif

  ! Synchronize and merge all parallel arrays

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (nv.eq..false.) then 
  CALL MPI_ALLREDUCE(mybx,bx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myby,by,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mybz,bz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(mydxbx,dxbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydybx,dybx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzbx,dzbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxby,dxby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydyby,dyby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzby,dzby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxbz,dxbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydybz,dybz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzbz,dzbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(mydxdxbx,dxdxbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdybx,dxdybx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdzbx,dxdzbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdxby,dxdxby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdyby,dxdyby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdzby,dxdzby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdxbz,dxdxbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdybz,dxdybz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxdzbz,dxdzbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(mydydybx,dydybx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydydzbx,dydzbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydydyby,dydyby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydydzby,dydzby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydydybz,dydybz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydydzbz,dydzbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(mydzdzbx,dzdzbx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzdzby,dzdzby,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzdzbz,dzdzbz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  endif
  CALL MPI_ALLREDUCE(myvx,vx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myvy,vy,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myvz,vz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 
  CALL MPI_ALLREDUCE(mydxvx,dxvx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydyvx,dyvx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzvx,dzvx,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxvy,dxvy,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydyvy,dyvy,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzvy,dzvy,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydxvz,dxvz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydyvz,dyvz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydzvz,dzvz,nkx0*nky0*nkz0&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)


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
  if (nv.eq..false.) then

  IF (mod(mype-1,n_mpi_procs).eq.0) THEN 
  store = (bx*dxvx+by*dyvx+bz*dzvx) - (vx*dxbx+vy*dybx+vz*dzbx)&
       - hall*((bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)&
       + (dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx))
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  myrhs_out_nlb(:,:,:,0) = temp_big*fft_norm
  ENDIF

  IF (mod(mype-2,n_mpi_procs).eq.0) THEN
  store = (bx*dxvy+by*dyvy+bz*dzvy) - (vx*dxby+vy*dyby+vz*dzby)&
       - hall*((bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)&
       + (dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby))
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  myrhs_out_nlb(:,:,:,1) = temp_big*fft_norm
  ENDIF

  IF (mod(mype-3,n_mpi_procs).eq.0) THEN
  store = (bx*dxvz+by*dyvz+bz*dzvz)- (vx*dxbz+vy*dybz+vz*dzbz)&
       - hall*((bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)&
       + (dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz ))
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  myrhs_out_nlb(:,:,:,2) = temp_big*fft_norm
  ENDIF

  endif

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


  if (rhs_nl_version.eq.12) then

     IF (mod(mype-4,n_mpi_procs).eq.0) THEN 
     store = -(vx*dxvx+vy*dyvx+vz*dzvx) + (by*dybx+bz*dzbx) - (by*dxby+bz*dxbz)
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,0) = temp_big*fft_norm
     IF(verbose) WRITE(*,*), "Max V NL ",mype,maxval(abs(myrhs_out_nlv))
     ENDIF 

     IF (mod(mype-5,n_mpi_procs).eq.0) THEN
     store = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+bz*dzby) - (bz*dybz+bx*dybx)
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,1) = temp_big*fft_norm
     IF(verbose) WRITE(*,*), "Max V NL ",mype,maxval(abs(myrhs_out_nlv))
     ENDIF

     IF (mod(mype-6,n_mpi_procs).eq.0) THEN
     store = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz) - (by*dzby+bx*dzbx)
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,2) = temp_big*fft_norm
     IF(verbose) WRITE(*,*), "Max V NL ",mype,maxval(abs(myrhs_out_nlv))
     ENDIF

  endif
  if (rhs_nl_version.eq.1) then 
     IF (mod(mype-4,n_mpi_procs).eq.0) THEN
     store = -(vx*dxvx+vy*dyvx+vz*dzvx) + (bx*dxbx+by*dybx+bz*dzbx) - 0.5*dxbmag
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,0) = temp_big*fft_norm
     ENDIF
     IF (mod(mype-5,n_mpi_procs).eq.0) THEN
     store = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+by*dyby+bz*dzby) - 0.5*dybmag
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,1) = temp_big*fft_norm
     ENDIF
     IF (mod(mype-6,n_mpi_procs).eq.0) THEN
     store = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz+bz*dzbz) - 0.5*dzbmag     
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
     myrhs_out_nlv(:,:,:,2) = temp_big*fft_norm
     ENDIF
  endif

! Find nonlinearities at energy timestep

     IF (plot_nls.and.(mod(itime,istep_energy).eq.0)) THEN
        ! b.grad v
        IF (mod(mype-7,n_mpi_procs).eq.0) THEN
        mybdv(:,:,:,0) = mybdv(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxvx+by*dyvx+bz*dzvz)
        ENDIF
        IF (mod(mype-8,n_mpi_procs).eq.0) THEN
        mybdv(:,:,:,1) = mybdv(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxvy+by*dyvy+bz*dzvy)
        ENDIF
        IF (mod(mype-9,n_mpi_procs).eq.0) THEN
        mybdv(:,:,:,2) = mybdv(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxvz+by*dyvz+bz*dzvz)
        ENDIF

        ! v.grad b
        IF (mod(mype-10,n_mpi_procs).eq.0) THEN
        myvdb(:,:,:,0) = myvdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxbx+vy*dybx+vz*dzbx)
        ENDIF
        IF (mod(mype-11,n_mpi_procs).eq.0) THEN
        myvdb(:,:,:,1) = myvdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxby+vy*dyby+vz*dzby)
        ENDIF
        IF (mod(mype-12,n_mpi_procs).eq.0) THEN
        myvdb(:,:,:,2) = myvdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxbz+vy*dybz+vz*dzbz)
        ENDIF

        ! b. grad curl b
        IF (mod(mype-13,n_mpi_procs).eq.0) THEN
        mybdcb(:,:,:,0) = mybdcb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)
        ENDIF
        IF (mod(mype-14,n_mpi_procs).eq.0) THEN
        mybdcb(:,:,:,1) = mybdcb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)
        ENDIF
        IF (mod(mype-15,n_mpi_procs).eq.0) THEN
        mybdcb(:,:,:,2) = mybdcb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)
        ENDIF

        ! curl b . grad b                       
        IF (mod(mype-16,n_mpi_procs).eq.0) THEN
        mycbdb(:,:,:,0) = mycbdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx)
        ENDIF
        IF (mod(mype-17,n_mpi_procs).eq.0) THEN
        mycbdb(:,:,:,1) = mycbdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby)
        ENDIF
        IF (mod(mype-18,n_mpi_procs).eq.0) THEN
        mycbdb(:,:,:,2) = mycbdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz)
        ENDIF

        ! v . grad v
        IF (mod(mype-19,n_mpi_procs).eq.0) THEN
        myvdv(:,:,:,0) = myvdv(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxvx+vy*dyvx+vz*dzvz)
        ENDIF
        IF (mod(mype-20,n_mpi_procs).eq.0) THEN
        myvdv(:,:,:,1) = myvdv(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxvy+vy*dyvy+vz*dzvy)
        ENDIF
        IF (mod(mype-21,n_mpi_procs).eq.0) THEN
        myvdv(:,:,:,2) = myvdv(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(vx*dxvz+vy*dyvz+vz*dzvz)
        ENDIF

        ! b . grad b
        IF (mod(mype-22,n_mpi_procs).eq.0) THEN
        mybdb(:,:,:,0) = mybdb(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxbx+by*dybx+bz*dzbz)
        ENDIF
        IF (mod(mype-23,n_mpi_procs).eq.0) THEN
        mybdb(:,:,:,1) = mybdb(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxby+by*dyby+bz*dzby)
        ENDIF
        IF (mod(mype-24,n_mpi_procs).eq.0) THEN
        mybdb(:,:,:,2) = mybdb(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxbz+by*dybz+bz*dzbz)
        ENDIF

        ! 0.5 grad b^2
        if (rhs_nl_version.eq.12) then
        IF (mod(mype-25,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,0) = mydb2(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dxbx+by*dxby+bz*dxbz)
        ENDIF
        IF (mod(mype-26,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,1) = mydb2(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dybx+by*dyby+bz*dybz)
        ENDIF
        IF (mod(mype-27,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,2) = mydb2(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(bx*dzbx+by*dzby+bz*dzbz)
        ENDIF
        endif 
        if (rhs_nl_version.eq.1) then
        IF (mod(mype-25,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,0) = mydb2(:,:,:,0) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(0.5*dxbmag)
        ENDIF
        IF (mod(mype-26,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,1) = mydb2(:,:,:,1) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(0.5*dybmag)
        ENDIF
        IF (mod(mype-27,n_mpi_procs).eq.0) THEN
        mydb2(:,:,:,2) = mydb2(:,:,:,2) + (5.0/12.0 - 1.0/6.0 * abs(real(counter) - 3.0/2.0)) * fft_spec2(0.5*dzbmag)
        ENDIF
        endif
    ENDIF

! Get nonlinear right hand sides and record nonlinearities

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mybdv,bdv,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myvdb,vdb,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mybdcb,bdcb,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mycbdb,cbdb,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myvdv,vdv,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mybdb,bdb,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mydb2,db2,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (nv.eq..false.) CALL MPI_ALLREDUCE(myrhs_out_nlb,rhs_out_nlb,nkx0*nky0*nkz0*3&
       ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(myrhs_out_nlv,rhs_out_nlv,nkx0*nky0*nkz0*3&
      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

  IF (verbose.and.(mype.eq.0)) WRITE(*,*) ierr
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'Abs NL Term Vx kmax',abs(rhs_out_nlv(nkx0-1,nky0/2,nkz0/2,0))
  IF (verbose.and.(mype.eq.0)) WRITE(*,*) 'Abs NL Term Bx kmax',abs(rhs_out_nlb(nkx0-1,nky0/2,nkz0/2,0))

  if (nv.eq..false.) rhs_out_b = rhs_out_b + rhs_out_nlb
  if (nv) rhs_out_b = cmplx(0.0,0.0)
  rhs_out_v = rhs_out_v + rhs_out_nlv

  if (verbose.and.(mype.eq.0)) print *, 'rhs out v nl found'

    IF ((mod(counter,4).eq.3).and.(mod(itime,istep_energy).eq.0)) THEN
       IF (mod(mype-1,n_mpi_procs).eq.0) THEN
          WRITE(bdvio) bdv(:,:,:,0)
          WRITE(bdvio) bdv(:,:,:,1)
          WRITE(bdvio) bdv(:,:,:,2)
       ENDIF
       IF (mod(mype-2,n_mpi_procs).eq.0) THEN
          WRITE(vdbio) vdb(:,:,:,0)
          WRITE(vdbio) vdb(:,:,:,1)
          WRITE(vdbio) vdb(:,:,:,2)
       ENDIF
       IF (mod(mype-3,n_mpi_procs).eq.0) THEN
          WRITE(bdcbio) bdcb(:,:,:,0)
          WRITE(bdcbio) bdcb(:,:,:,1)
          WRITE(bdcbio) bdcb(:,:,:,2)
       ENDIF
       IF (mod(mype-4,n_mpi_procs).eq.0) THEN
          WRITE(cbdbio) cbdb(:,:,:,0)
          WRITE(cbdbio) cbdb(:,:,:,1)
          WRITE(cbdbio) cbdb(:,:,:,2)
       ENDIF
       IF (mod(mype-5,n_mpi_procs).eq.0) THEN
          WRITE(vdvio) vdv(:,:,:,0)
          WRITE(vdvio) vdv(:,:,:,1)
          WRITE(vdvio) vdv(:,:,:,2)
       ENDIF
       IF (mod(mype-6,n_mpi_procs).eq.0) THEN
          WRITE(bdbio) bdb(:,:,:,0)
          WRITE(bdbio) bdb(:,:,:,1)
          WRITE(bdbio) bdb(:,:,:,2)
       ENDIF
       IF (mod(mype-7,n_mpi_procs).eq.0) THEN
          WRITE(db2io) db2(:,:,:,0)
          WRITE(db2io) db2(:,:,:,1)
          WRITE(db2io) db2(:,:,:,2)
       ENDIF
     
        if (verbose.and.(mype.eq.0)) print *, 'nonlinearities written'

    ENDIF



  if (calc_dt) CALL next_dt(ndt)
  if (.not.(calc_dt)) ndt = dt_max
  if (verbose.and.(mype.eq.0)) print *, 'next dt calculated ',ndt

  DEALLOCATE(temp_small)
  if (verbose.and.(mype.eq.0)) print *, 'ts deallocated'
 
  DEALLOCATE(temp_bigx)
  if (verbose.and.(mype.eq.0)) print *, 'tbx deallocated'
  DEALLOCATE(temp_bigy)
  if (verbose.and.(mype.eq.0)) print *, 'tby deallocated'
  DEALLOCATE(temp_bigz)
  if (verbose.and.(mype.eq.0)) print *, 'tbz deallocated'

  DEALLOCATE(store_x)
  if (verbose.and.(mype.eq.0)) print *, 'stx deallocated'
  DEALLOCATE(store_y)
  if (verbose.and.(mype.eq.0)) print *, 'sty deallocated'
  DEALLOCATE(store_z)
  if (verbose.and.(mype.eq.0)) print *, 'stz deallocated'

  ! All b arrays
  DEALLOCATE(bx)
  if (verbose.and.(mype.eq.0)) print *, 'bx deallocated'
  DEALLOCATE(by)
  if (verbose.and.(mype.eq.0)) print *, 'by deallocated'
  DEALLOCATE(bz)
  if (verbose.and.(mype.eq.0)) print *, 'first third deallocated'

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
     !DEALLOCATE(bmag_inbig)
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
  DEALLOCATE(dzvz) ! all first order b arrays
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
  ! all ind processors
  DEALLOCATE(mytemp_small)

  DEALLOCATE(mybx)
  DEALLOCATE(myby)
  DEALLOCATE(mybz)

  DEALLOCATE(myvx)
  DEALLOCATE(myvy)
  DEALLOCATE(myvz)

  DEALLOCATE(mydxvx)
  DEALLOCATE(mydyvx)
  DEALLOCATE(mydzvx)

  DEALLOCATE(mydxvy)
  DEALLOCATE(mydyvy)
  DEALLOCATE(mydzvy)

  DEALLOCATE(mydxvz)
  DEALLOCATE(mydyvz)
  DEALLOCATE(mydzvz)

  DEALLOCATE(mydxbx)
  DEALLOCATE(mydybx)
  DEALLOCATE(mydzbx)

  DEALLOCATE(mydxby)
  DEALLOCATE(mydyby)
  DEALLOCATE(mydzby)

  DEALLOCATE(mydxbz)
  DEALLOCATE(mydybz)
  DEALLOCATE(mydzbz)
  
  DEALLOCATE(mydxdxbx)
  DEALLOCATE(mydxdybx)
  DEALLOCATE(mydxdzbx)
  DEALLOCATE(mydydybx)
  DEALLOCATE(mydydzbx)
  DEALLOCATE(mydzdzbx)

  DEALLOCATE(mydxdxby)
  DEALLOCATE(mydxdyby)
  DEALLOCATE(mydxdzby)
  DEALLOCATE(mydydyby)
  DEALLOCATE(mydydzby)
  DEALLOCATE(mydzdzby)

  DEALLOCATE(mydxdxbz)
  DEALLOCATE(mydxdybz)
  DEALLOCATE(mydxdzbz)
  DEALLOCATE(mydydybz)
  DEALLOCATE(mydydzbz)
  DEALLOCATE(mydzdzbz)

  DEALLOCATE(rhs_out_nlb)
  DEALLOCATE(rhs_out_nlv)
  DEALLOCATE(myrhs_out_nlb)
  DEALLOCATE(myrhs_out_nlv)

  ! current bv arrays
  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

  ! nonlinearities
  IF(mod(counter,4).eq.3) THEN
  DEALLOCATE(bdv)
  DEALLOCATE(vdb)
  DEALLOCATE(bdcb)
  DEALLOCATE(cbdb)
  DEALLOCATE(vdv)
  DEALLOCATE(bdb)
  DEALLOCATE(db2)
  DEALLOCATE(mybdv)
  DEALLOCATE(myvdb)
  DEALLOCATE(mybdcb)
  DEALLOCATE(mycbdb)
  DEALLOCATE(myvdv)
  DEALLOCATE(mybdb)
  DEALLOCATE(mydb2)
  ENDIF

  if (verbose.and.(mype.eq.0)) CALL cpu_time(dum)
  if (verbose.and.(mype.eq.0)) WRITE(*,*) 'time for nl code', dum-sttime

  if (verbose.and.(mype.eq.0)) print *, 'deallocated nl code'

END SUBROUTINE get_rhs_nl2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl3                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_rhs_nl3(b_in,v_in,dx,dy,dz,rhs_out_b,rhs_out_v,ndt)

USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: dx, dy, dz
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  REAL :: ndt

  INTEGER :: i,j,l,k,h, ierr

  ALLOCATE(ekd(0:nkx0-1,0:nky0-1,0:nkz0-1))

  ekd = cmplx(0.0,0.0)
  DO i = 0,nkx0-1
     DO j = 0,nky0-1
        DO k = 0,nkz0-1
           ekd(i,j,k) = exp(i_complex * i * dx) * exp(i_complex * j * dy) * exp(i_complex * k * dz)
        ENDDO
     ENDDO
  ENDDO
  
  IF(np_spec.gt.1) STOP "get_rhs_nl1 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl1 not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_nl1 only suitable for np_kz=1"
  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_bigx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(temp_bigz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(store_x(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_y(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(store_z(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ! All b arrays
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  !  if (rhs_nl_version==1) then
  !    ALLOCATE(bmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(bmag_in(0:nkx0-1,0:nky0-1,0:nkz0-1))
  !    ALLOCATE(bmag_inbig(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dxbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dybmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !    ALLOCATE(dzbmag(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  !  endif
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

  ! all  second order vx arrays i.e. DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  ALLOCATE(dxdxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order vy arrays
  ALLOCATE(dxdxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ! all  second order bz arrays
  ALLOCATE(dxdxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

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

  ! bmag should be calculated in real space
  !  DO i=0,nkx0-1
  !    DO j=0,nky0-1
  !        DO k=0,nkz0-1
  !            bmag_in(i,j,k) = b_inx0(i,j,k)*b_inx0(i,j,k)+ b_iny0(i,j,k)*b_iny0(i,j,k)+ b_inz0(i,j,k)*b_inz0(i,j,k)
  !        END DO
  !    END DO
  !  END DO

  !   SECOND ORDER  VX TERMS DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
  if (nv.eq..false.) then
  if (.false.) then
  !dxdxvx
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inx0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvx = store

  !dxdyvx
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inx0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvx = store

  !dxdzvx
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inx0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvx = store

  ! DYDYVX
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inx0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvx = store

  ! DYDzVX
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inx0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvx = store

  ! DzDzVX
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_inx0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvx = store

  ! finished SECOND ORDER  VX TERMS

  ! SECOND ORDER  VY TERMS DXDXVY,DXDYVY,DXDZVY, DYDYVY,DYDZVY, DZDZVY

  !dxdxvy
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_iny0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvy = store

  !dxdyvy
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_iny0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvy = store

  !dxdzvy
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_iny0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvy = store

  ! DYDYVy
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_iny0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvy = store

  ! DYDzVY
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_iny0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvy = store

  ! DzDzVY
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_iny0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvy = store

  ! finished SECOND ORDER VY TERMS

  ! SECOND ORDER  VZ TERMS DXDXVZ,DXDYVZ,DXDZVZ, DYDYVz,DYDZVz, DZDZVz
  !dxdxvz
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*v_inz0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxvz = store

  !dxdyvz
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*v_inz0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyvz = store

  !dxdzvz
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*v_inz0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzvz = store

  ! DYDYVz
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*v_inz0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyvz = store

  ! DYDzVz
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*v_inz0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzvz = store

  ! DzDzVz
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*v_inz0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzvz = store

  ! finished SECOND ORDER VZ TERMS
  endif

  ! START SECOND ORDER b terms

  ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
  !dxdxbx
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inx0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxbx = store

  !dxdybx
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inx0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdybx = store

  !dxdzbx
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inx0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzbx = store

  ! DYDYbX
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inx0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydybx = store

  ! DYDzbX
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inx0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzbx = store

  ! DzDzbX
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inx0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzbx = store

  ! FINISHED SECOND ORDER  BX TERMS

  ! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY

  !dxdxby
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_iny0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxby = store

  !dxdyvy
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_iny0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdyby = store

  !dxdzby
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_iny0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzby = store

  ! DYDYby
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_iny0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydyby = store

  ! DYDzbY
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_iny0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzby = store

  ! DzDzbY
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_iny0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzby = store

  ! FINISHED SECOND ORDER bY TERMS

  ! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
  !dxdxbz
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*i_complex*kxgrid(i)*b_inz0(i,:,:) ! there is  two i's in the
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdxbz = store

  !dxdybz
  DO i=0,nkx0-1
     DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*i_complex*kygrid(j)*b_inz0(i,j,:)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdybz = store

  !dxdzbz
  DO i=0,nkx0-1
     DO k=0,nkz0-1
        temp_small(i,:,k)=i_complex*kxgrid(i)*i_complex*kzgrid(k)*b_inz0(i,:,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxdzbz = store

  ! DYDYbz
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*i_complex*kygrid(j)*b_inz0(:,j,:)  ! towo y grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydybz = store

  ! DYDzbz
  DO j=0,nky0-1
     DO k=0,nkz0-1
        temp_small(:,j,k)=i_complex*kygrid(j)*i_complex*kzgrid(k)*b_inz0(:,j,k)
     END DO
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dydzbz = store

  ! DzDzbz
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*i_complex*kzgrid(k)*b_inz0(:,:,k)  !two kz grid
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzdzbz = store

  !finished END SECOND ORDER BZ TERMS

  !    completed  ALL SECOND ORDER B TERMS

  ! TERMS BX BY BZ

  !bx
  DO i=0,nkx0-1
     temp_small(i,:,:)=b_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  bx = store

  !by
  DO j=0,nky0-1
     temp_small(:,j,:)=b_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  by = store

  !bz
  DO k=0,nkz0-1
     temp_small(:,:,k)=b_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  bz = store

!!! bmag terms
  if (rhs_nl_version == 1) then
     !! Calculate b^2(x) in real space.
     DO i=0,nx0_big-1
        DO j=0,ny0_big-1
           DO k=0,nz0_big-1
              bmag(i,j,k) = bx(i,j,k)*conjg(bx(i,j,k))+ by(i,j,k)*conjg(by(i,j,k))+ bz(i,j,k)*conjg(bz(i,j,k))
           END DO
        END DO
     END DO
!!! Now forward transfrom to get b^2(k) stored in temp_big
     temp_big=cmplx(0.0,0.0)
     store = bmag
     CALL dfftw_execute_dft(plan_r2c,store,temp_big)
!!! Now cut down to dealias and store result in bmagk which is b^2(k)
     
     DO i=0,nkx0-1
        bmagk(i,0:nky0-1,0:nkz0-1)= temp_big(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
     END DO!i loop
     
!!! Now that we have bmagk, need to do IFFT(ikx b^2) IFFT(iky b^2) IFFT(ikz b^2) to get grad b^2/2 in real space

     !bmag_inbig = bx*bx + by*by + bz*bz
     !CALL dfftw_execute_dft(plan_r2c,bmag_inbig,bmag)
     !do i = 0,nkx0-1
     !        bmag_in(i,0:hky_ind,0:hkz_ind)=bmag(i,0:hky_ind,0:hkz_ind)*fft_norm   !kz positive, ky positive
     !        bmag_in(i,0:hky_ind,lkz_ind:nkz0-1)=bmag(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm !kz negative, ky positive
     !        bmag_in(i,lky_ind:nky0-1,lkz_ind:nkz0-1) = bmag(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm ! kz negative, ky negative
     !        bmag_in(i,lky_ind:nky0-1,0:hkz_ind)=bmag(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm !kz positive, ky negative
     !end do

     !dxbmag
     DO i=0,nkx0-1
        temp_small(i,:,:)=i_complex*kxgrid(i)*bmagk(i,:,:)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dxbmag = store

     !dybmag
     DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*bmagk(:,j,:)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dybmag = store

     !dzbmag
     DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*bmagk(:,:,k)
     END DO
     !Add padding for dealiasing
     temp_big=cmplx(0.0,0.0)
     DO i=0,nkx0-1
        temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
         !kz negative, ky positive
         !kz negative, ky negative
         !kz positive, ky negative
     END DO!k loop
     CALL dfftw_execute_dft(plan_c2r,temp_big,store)
     dzbmag = store
     if (verbose) print *, 'Bmag gradients dealiased'
  endif
  endif 
!!! TERMS  vx,vy,vz
  !vx
  DO i=0,nkx0-1
     temp_small(i,:,:)=v_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  vx = store

  !vy
  DO j=0,nky0-1
     !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
     temp_small(:,j,:)=v_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  vy = store

  !vz
  DO k=0,nkz0-1
     temp_small(:,:,k)=v_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  vz = store

  !  FIRST ORDER VX TERMS DXVX , DYVX,  DZVX

  ! dxvx
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*v_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxvx = store

  ! dyvx
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*v_inx0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dyvx = store

  ! dzvx
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*v_inx0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzvx = store

  !  FIRST ORDER VY TERMS  dxvy dyvy dzvz,

  ! dxvy
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*v_iny0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxvy = store

  ! dyvy
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*v_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dyvy = store 

  ! dzvy
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*v_iny0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzvy = store

  !  FIRST ORDER VZ TERMS  dxvz dyvz dzvz
  ! dxvz
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*v_inz0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxvz = store

  ! dyvz
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*v_inz0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dyvz = store

  ! dzvz
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*v_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzvz = store

  ! DONE ALL FIRST ORDER VX,VY AND VZ TERMS i.e.  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
  if (nv.eq..false.) then
  ! FIRST ORDER BX TERMS ie. dxbx dybx dzbx`
  ! dxbx
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*b_inx0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxbx = store

  ! dybx
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*b_inx0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dybx = store

  ! dzbx
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*b_inx0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzbx = store

  !  FIRST ORDER BY TERMS ie. dxby dyby dzby
  ! dxby
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*b_iny0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxby = store

  ! dyby
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*b_iny0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dyby = store

  ! dzby
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*b_iny0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzby = store

  !! FIRST ORDER BZ TERMS ie. dxbz dybz dzbz
  ! dxbz
  DO i=0,nkx0-1
     temp_small(i,:,:)=i_complex*kxgrid(i)*b_inz0(i,:,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dxbz = store

  ! dybz
  DO j=0,nky0-1
     temp_small(:,j,:)=i_complex*kygrid(j)*b_inz0(:,j,:)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dybz = store

  ! dzbz
  DO k=0,nkz0-1
     temp_small(:,:,k)=i_complex*kzgrid(k)*b_inz0(:,:,k)
  END DO
  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
     temp_big(i,0:nky0-1,0:nkz0-1)=temp_small(i,0:nky0-1,0:nkz0-1)*ekd(i,0:nky0-1,0:nkz0-1)    !kz positive, ky positive
      !kz negative, ky positive
      !kz negative, ky negative
      !kz positive, ky negative
  END DO!k loop
  CALL dfftw_execute_dft(plan_c2r,temp_big,store)
  dzbz = store

  ! DONE ALL FIRST ORDER BX,BY BZ TERMS ie. dxbx dybx dzbx,  dxby, dyby,dzby,  dxbz,dybz,dzbz
  endif
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
  store_x = (bx*dxvx+by*dyvx+bz*dzvx) - (vx*dxbx+vy*dybx+vz*dzbx)&
       - hall*((bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)&
       + (dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx))

  store_y = (bx*dxvy+by*dyvy+bz*dzvy) - (vx*dxby+vy*dyby+vz*dzby)&
       - hall*((bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)&
       + (dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby))

  store_z = (bx*dxvz+by*dyvz+bz*dzvz)- (vx*dxbz+vy*dybz+vz*dzbz)&
       - hall*((bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)&
       + (dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz ))

  IF ((mod(counter,4).eq.0).and.((plot_nls.and.(max(dx,dy,dz).lt.10.0**(-10.0))).and.(mod(itime,istep_energy).eq.0)).and.(max(dx,dy,dz).le.10.0**(-10.0))) THEN
     ! b.grad v 
     WRITE(bdvio) fft_spec3(bx*dxvx+by*dyvx+bz*dzvx)
     WRITE(bdvio) fft_spec3(bx*dxvy+by*dyvy+bz*dzvy)
     WRITE(bdvio) fft_spec3(bx*dxvz+by*dyvz+bz*dzvz)

     ! v.grad b
     WRITE(vdbio) fft_spec3(vx*dxbx+vy*dybx+vz*dzbx)
     WRITE(vdbio) fft_spec3(vx*dxby+vy*dyby+vz*dzby)
     WRITE(vdbio) fft_spec3(vx*dxbz+vy*dybz+vz*dzbz)

     ! b. grad curl b
     WRITE(bdcbio) fft_spec3(bx*dxdybz+by*dydybz+bz*dydzbz-bx*dxdzby-by*dydzby-bz*dzdzby)
     WRITE(bdcbio) fft_spec3(bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dxdybz-bz*dxdzbz)
     WRITE(bdcbio) fft_spec3(bx*dxdxby+by*dxdyby+bz*dxdzby-bx*dxdybx-by*dydybx-bz*dydzbx)

     ! curl b . grad b
     WRITE(cbdbio) fft_spec3(dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dzbx)
     WRITE(cbdbio) fft_spec3(dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dzby)
     WRITE(cbdbio) fft_spec3(dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dzbz)
  ENDIF

  !inverse FFT to get back to Fourier

  store = store_x
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigx = temp_big

  store = store_y
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigy = temp_big

  store = store_z
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigz = temp_big

  !Now fill in appropriate rhs elements                                                                                                                                                                
  DO i=0,nkx0-1
        rhs_out_b(i,:,:,0)=rhs_out_b(i,:,:,0)+temp_bigx(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
        rhs_out_b(i,:,:,1)=rhs_out_b(i,:,:,1)+temp_bigy(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
        rhs_out_b(i,:,:,2)=rhs_out_b(i,:,:,2)+temp_bigz(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
  END DO!i loop  
 
  endif

  if (nv) rhs_out_b = cmplx(0.0,0.0)
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

     IF ((mod(counter,4).eq.0).and.((plot_nls.and.(max(dx,dy,dz).lt.10.0**(-10.0))).and.(mod(itime,istep_energy).eq.0)).and.(max(dx,dy,dz).le.10.0**(-10.0))) THEN
        ! v . grad v
        WRITE(vdvio) fft_spec3(vx*dxvx+vy*dyvx+vz*dzvx)
        WRITE(vdvio) fft_spec3(vx*dxvy+vy*dyvy+vz*dzvy)
        WRITE(vdvio) fft_spec3(vx*dxvz+vy*dyvz+vz*dzvz)

        ! b . grad b
        WRITE(bdbio) fft_spec3(bx*dxbx+by*dybx+bz*dzbx)
        WRITE(bdbio) fft_spec3(bx*dxby+by*dyby+bz*dzby)
        WRITE(bdbio) fft_spec3(bx*dxbz+by*dybz+bz*dzbz)

        ! 0.5 grad b^2
        WRITE(db2io) fft_spec3(0.5*dxbmag)
        WRITE(db2io) fft_spec3(0.5*dybmag)
        WRITE(db2io) fft_spec3(0.5*dzbmag)

        if (verbose) print *, 'v1 nl equation stored'
     endif
  ENDIF

  if (rhs_nl_version == 12) then
     store_x = -(vx*dxvx+vy*dyvx+vz*dzvx) + (by*dybx+bz*dzbx) - (by*dxby+bz*dxbz)
     store_y = -(vx*dxvy+vy*dyvy+vz*dzvy) + (bx*dxby+bz*dzby) - (bz*dybz+bx*dybx)
     store_z = -(vx*dxvz+vy*dyvz+vz*dzvz) + (bx*dxbz+by*dybz) - (by*dzby+bx*dzbx)

     IF ((mod(counter,4).eq.0).and.((plot_nls.and.(max(dx,dy,dz).lt.10.0**(-10.0))).and.(mod(itime,istep_energy).eq.0)).and.(max(dx,dy,dz).le.10.0**(-10.0))) THEN
        ! v . grad v
        WRITE(vdvio) fft_spec3(vx*dxvx+vy*dyvx+vz*dzvx)
        WRITE(vdvio) fft_spec3(vx*dxvy+vy*dyvy+vz*dzvy)
        WRITE(vdvio) fft_spec3(vx*dxvz+vy*dyvz+vz*dzvz)

        ! b . grad b
        WRITE(bdbio) fft_spec3(bx*dxbx+by*dybx+bz*dzbx)
        WRITE(bdbio) fft_spec3(bx*dxby+by*dyby+bz*dzby)
        WRITE(bdbio) fft_spec3(bx*dxbz+by*dybz+bz*dzbz)

        ! 0.5 grad b^2
        WRITE(db2io) fft_spec3(bx*dxbx+by*dxby+bz*dxbz)
        WRITE(db2io) fft_spec3(bx*dybx+by*dyby+bz*dybz)
        WRITE(db2io) fft_spec3(bx*dzbx+by*dzby+bz*dzbz)

        if (verbose) print *, 'v12 nl equation stored'
     endif
  ENDIF

  store = store_x
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigx = temp_big

  store = store_y
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigy = temp_big

  store = store_z
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temp_bigz = temp_big

  !Now fill in appropriate rhs elements
  DO i=0,nkx0-1
        rhs_out_v(i,:,:,0)=rhs_out_v(i,:,:,0)+temp_bigx(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
        rhs_out_v(i,:,:,1)=rhs_out_v(i,:,:,1)+temp_bigy(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
        rhs_out_v(i,:,:,2)=rhs_out_v(i,:,:,2)+temp_bigz(i,0:nky0-1,0:nkz0-1)/ekd(i,0:nky0-1,0:nkz0-1)*fft_norm
  END DO!i loop

  if (verbose) print *, 'rhs out v nl found'

  if (calc_dt) CALL next_dt(ndt)
  if (.not.(calc_dt)) ndt = dt_max
  if (verbose) print *, 'next dt calculated ',ndt

  DEALLOCATE(ekd)
  if (verbose) print *, 'ekd deallocated'
  DEALLOCATE(temp_small)
  if (verbose) print *, 'ts deallocated'
 
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
     !DEALLOCATE(bmag_inbig)
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
  DEALLOCATE(b_inx0)
  DEALLOCATE(b_iny0)
  DEALLOCATE(b_inz0)
  DEALLOCATE(v_inx0)
  DEALLOCATE(v_iny0)
  DEALLOCATE(v_inz0)

  if (verbose) print *, 'deallocated nl code'

END SUBROUTINE get_rhs_nl3

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

FUNCTION fft_spec(arr_real) result(arr_spec)

implicit none
complex :: arr_real(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
complex :: arr_spec(0:nkx0-1,0:nky0-1,lkz1:lkz2)
complex :: temporary(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
integer :: s

store = arr_real
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temporary = temp_big

!Now fill in appropriate rhs elements                                          

  DO s=0,nkx0-1
            arr_spec(s,0:hky_ind,0:hkz_ind)=temporary(s,0:hky_ind,0:hkz_ind)*fft_norm             !kz positive, ky positive 
            arr_spec(s,0:hky_ind,lkz_ind:nkz0-1)=temporary(s,0:hky_ind,lkz_big:nz0_big-1)*fft_norm     !kz negative, ky positive
            arr_spec(s,lky_ind:nky0-1,lkz_ind:nkz0-1)=temporary(s,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm       !kz negative, ky negative  
            arr_spec(s,lky_ind:nky0-1,0:hkz_ind)=temporary(s,lky_big:ny0_big-1,0:hkz_ind)*fft_norm     !kz positive, ky negative  
  END DO

END FUNCTION fft_spec

FUNCTION fft_spec2(arr_real) result(arr_spec)

implicit none
complex :: arr_real(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
complex :: arr_spec(0:nkx0-1,0:nky0-1,lkz1:lkz2)
complex :: temporary(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
integer :: s

  store = arr_real
  CALL dfftw_execute_dft(plan_r2c,store,temp_big)
  temporary = temp_big
!Now fill in appropriate rhs elements 

  DO s=0,nkx0-1
     arr_spec(s,:,:)=temporary(s,0:nky0-1,0:nkz0-1)*fft_norm
  END DO!i loop

END FUNCTION fft_spec2

FUNCTION fft_spec3(arr_real) result(arr_spec)

implicit none
complex :: arr_real(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
complex :: arr_spec(0:nkx0-1,0:nky0-1,lkz1:lkz2)
complex :: temporary(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
integer :: s

store = arr_real
CALL dfftw_execute_dft(plan_r2c,store,temp_big)
temporary = temp_big

!Now fill in appropriate rhs elements 

  DO s=0,nkx0-1
        arr_spec(s,:,:)=temporary(s,0:nky0-1,0:nkz0-1)/ekd(s,0:nky0-1,0:nkz0-1)*fft_norm
  END DO!i loop  

END FUNCTION fft_spec3

SUBROUTINE finalize_fourier

implicit none

call dfftw_destroy_plan(plan_r2c)
call dfftw_destroy_plan(plan_c2r)
if (verbose) print *, "Destroyed plans"
call dfftw_cleanup()
if (verbose.and.(mype.eq.0)) print *, "Destroyed plans"
DEALLOCATE(store)
if (verbose.and.(mype.eq.0)) print *, "Deallocated store"
DEALLOCATE(temp_big)
if (verbose.and.(mype.eq.0)) print *, "Deallocated temp_big"
DEALLOCATE(plan_r2c)
DEALLOCATE(plan_c2r)
if (verbose.and.(mype.eq.0)) print *, "Deallocated plans"

END SUBROUTINE finalize_fourier

SUBROUTINE wherenezero(fftop)

implicit none
complex :: fftop(0:nkx0-1,0:nky0-1,0:nkz0-1)
integer :: ix,iy,iz

DO ix = 0,nkx0-1
  DO iy = 0,nky0-1
    DO iz = 0,nkz0-1
      IF (abs(fftop(ix,iy,iz)).gt.10.0**(-10.0)) print *, ix,iy,iz,fftop(ix,iy,iz)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE wherenezero

SUBROUTINE wherebvnotreal(fftop1,fftop2)

implicit none
complex :: fftop1(0:nkx0-1,0:nky0-1,0:nkz0-1)
complex :: fftop2(0:nkx0-1,0:nky0-1,0:nkz0-1)
integer :: ix,iy,iz

DO ix = 0,nkx0/2
  DO iy = 0,nky0/2
    DO iz = 0,nkz0/2
      IF (abs(fftop1(ix,iy,iz)-conjg(fftop1(mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0)))).gt.10.0**(-10.0)) &
        print *, 'B',ix,iy,iz,mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0),fftop1(ix,iy,iz),fftop1(mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0))
      IF (abs(fftop1(ix,iy,iz)-conjg(fftop1(mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0)))).gt.10.0**(-10.0)) &
        print *, 'V',ix,iy,iz,mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0),fftop2(ix,iy,iz),fftop2(mod(nkx0-ix,nkx0),mod(nky0-iy,nky0),mod(nkz0-iz,nkz0))
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE wherebvnotreal


!SUBROUTINE convcomp(x,y)
!
!implicit none
!
!complex :: x(0:nkx0-1,0:nky0-1,0:nkz0-1)
!complex :: y(0:nkx0-1,0:nky0-1,0:nkz0-1)
!complex :: convo(0:nkx0-1,0:nky0-1,0:nkz0-1)
!integer :: i,j,k
!integer :: l,m,n

! compute convolution

!DO i = 0,nkx0-1
!  DO j = 0,nky0-1
!    DO k = 0,nkz0-1
!      DO l = 0,i+1
!        DO m = 0,j+1
!          DO n = 0,k+1
!            if ((i.ge.l).and.(j.ge.m).and.(k.ge.n)) convo(i,j,k) = convo(i,j,k) + x(l,m,n) * y(i-l,j-m,k-n)
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO

! compute fft convolution

!SELECT CASE (dealias_type)

!CASE(1)

!temp_big = cmplx(0.0,0.0)
!temp_big(0:nkx0-1,0:nky0-1,0:nkz0-1) = x * ekd
!temp_big(0:nkx0-1,0:nky0-1,0:nkz0-1) 

!CASE(3)

!CASE(4)

!END SELECT

!END SUBROUTINE convcomp

END MODULE nonlinearity

