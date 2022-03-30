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
  USE hk_effects
  USE flr_effects
  IMPLICIT NONE

  PUBLIC :: initialize_fourier,get_rhs_nl,&
            get_rhs_nl_convolution,get_k_indices,get_rhs_nl2,get_rhs_nl1,&
            initialize_fourier_ae_mu0 !,initialize_fourier2
  
  REAL, PUBLIC :: ve_max(2)

  PRIVATE

                                                          !COMPLEX :: temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1)
                                                          !COMPLEX :: temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1)
                                                          !REAL :: dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                          !REAL :: dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                          !REAL :: dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
                                                          !REAL :: dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: g_in0
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: temp_small,temp_big
                                                          !REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxphi,dyphi,dxg,dyg
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: bx,by,bz, dxbx, dybx,dzbx, dxby,dyby,dzbz, dxbz,dybz,dzbz
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbx, dxdybx, dxdzbx, dxdybx, dydybx, dzdybx, dxdzbx, dydzbx, dzdzbx
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxby, dxdyby, dxdzby, dxdyby, dydyby, dzdyby, dxdzby, dydzby, dzdzby
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxbz, dxdybz, dxdzbz, dxdybz, dydybz, dzdybz, dxdzbz, dydzbz, dzdzbz
    
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: vx,vy,vz, dxvx, dyvx,dzvx, dxvy,dyvy,dzvz, dxvz,dyvz,dzvz
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvx, dxdyvx, dxdzvx, dxdyvx, dydyvx, dzdyvx, dxdzvx, dydzvx, dzdzvx
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvy, dxdyvy, dxdzvy, dxdyvy, dydyvy, dzdyvy, dxdzvy, dydzvy, dzdzvy
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dxdxvz, dxdyvz, dxdzvz, dxdyvz, dydyvz, dzdyvz, dxdzvz, dydzvz, dzdzvz
  
                                                        ! COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: temp_small_2d,temp_big_2d
                                                        ! REAL, ALLOCATABLE, DIMENSION(:,:) :: dxphi_2d,dyphi_2d,dxg_2d,dyg_2d


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
  nx0_big=3*nkx0
  ny0_big=3*nky0/2
  nz0_big=3*nkz0/2
  fft_norm=1.0/(REAL(nx0_big*ny0_big*nz0_big))

  ALLOCATE(g_rbig(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(g_kbig(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))

  !WRITE(*,*) "making plans"
  CALL dfftw_plan_dft_c2r_3d(plan_c2r,nx0_big,ny0_big,nz0_big,&
                             g_kbig,g_rbig,FFTW_ESTIMATE)
  CALL dfftw_plan_dft_r2c_3d(plan_r2c,nx0_big,ny0_big,nz0_big,&
                             g_rbig,g_kbig,FFTW_ESTIMATE)
  
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
SUBROUTINE get_rhs_nl(b_in, v_in, rhs_out_b, rhs_out_v)
  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

  IF(rhs_nl_version==1) THEN
    CALL get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v)
  ELSE IF(rhs_nl_version==2) THEN
    CALL get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
  ELSE IF(rhs_nl_version==3) THEN
    CALL get_rhs_nl3(b_in,v_in,rhs_out_b,rhs_out_v)
  ELSE IF(rhs_nl_version==4) THEN
    CALL get_rhs_nl4(b_in,v_in,rhs_out_b,rhs_out_v)
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
SUBROUTINE get_rhs_nl1(b_in,v_in,rhs_out_b,rhs_out_v)

  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX, INTENT(inout) :: rhs_out_v(0:nkx0-1,0:nky0-1,0:2)
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
  ALLOCATE(temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
! All b arrays 
  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
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
                                                        !ALLOCATE(dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
                                                        ! ALLOCATE(dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
                                                        !ALLOCATE(dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
                                                        ! ALLOCATE(dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
                                                        !ALLOCATE(g_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
  ALLOCATE(b_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(b_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inx0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_iny0(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  ALLOCATE(v_inz0(0:nkx0-1,0:nky0-1,lkz1:lkz2))

                                                                  ! Some test for the inverse hankel transform
                                                                  !  IF (mype==0) THEN
                                                                  !  DO i=0,nkx0-1
                                                                  !DELETED NO NEED 
                                                                  !    END DO
                                                                  !  ENDIF
 
  ! I dont want to change g_in, so I copy temporaly to g_in0
  !g_in0 = g_in
  b_inx0 = b_in(:,:,:,0)
  b_iny0 = b_in(:,:,:,1)
  b_inz0 = b_in(:,:,:,2)
  v_inx0 = v_in(:,:,:,0)
  v_iny0 = v_in(:,:,:,1)
  v_inz0 = v_in(:,:,:,2)

                                              !IF (hankel) THEN
                                              !!First I need the Hankel transform of g_in
                                              !call hankel_transform(g_in0,.true.)    
                                              !ENDIF

                                             ! Now a loop versus the hankel index
                                             !This check for hankel, but for v_on should work with lh1:lh2
                                              !DO h=lh1,lh2
                                                !IF(mype==0.and.first_stage) WRITE(*,*) "mype,itime,rhs_lin",mype,itime,abs(sum(sum(sum(sum(rhs_out,1),1),1),1))
                                               !  DO k=0,nkz0-1
                                               !     phi_in(:,:,k)=J0a(:,:)*J0_fac(:,:,h)*phi_in0(:,:,k)
                                               ! END DO
   
   !   SECOND ORDER  VX TERMS DXDXVX,   DXDYVX,   DXDZVX,  DYDYVX,   DYDZVX, DZDZVX
    !dxdxvx
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*v_inx0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxvx(0,0,0))
   
   !dxdyvx
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*v_inx0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdyvx(0,0,0))
    
     !dxdzvx
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*v_inx0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzvx(0,0,0))
   ! DYDYVX
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*v_inx0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyvx(0,0,0))
    
     ! DYDzVX
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*v_inx0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzvx(0,0,0))
    
     ! DzDzVX
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*v_inx0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzvx(0,0,0))
   
   ! finished SECOND ORDER  VX TERMS
   
! SECOND ORDER  VY TERMS DXDXVY,DXDYVY,DXDZVY, DYDYVY,DYDZVY, DZDZVY 

  !dxdxvy
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*v_iny0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxvy(0,0,0))
   
   !dxdyvy
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*v_iny0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdyvy(0,0,0))
    
     !dxdzvy
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*v_iny0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzvy(0,0,0))
   ! DYDYVy
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*v_iny0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyvy(0,0,0))
    
     ! DYDzVY
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*v_iny0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzvy(0,0,0))
    
     ! DzDzVY
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*v_iny0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzvy(0,0,0))
    
! finished SECOND ORDER VY TERMS

! SECOND ORDER  VZ TERMS DXDXVZ,DXDYVZ,DXDZVZ, DYDYVz,DYDZVz, DZDZVz
!dxdxvz
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*v_inz0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxvz(0,0,0))
   
   !dxdyvz
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*v_inz0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdyvz(0,0,0))
    
     !dxdzvz
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*v_inz0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzvz(0,0,0))
   ! DYDYVz
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*v_inz0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyvz(0,0,0))
    
     ! DYDzVz
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*v_inz0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzvz(0,0,0))
    
     ! DzDzVz
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*v_inz0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzvz(0,0,0))
    
! finished SECOND ORDER VZ TERMS

! START SECOND ORDER b terms  

 ! SECOND ORDER  BX TERMS DXDXBX,DXDYBX,DXDZBX,  DYDYBX,DYDZBX, DZDZBX
    !dxdxbx
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*b_inx0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxbx(0,0,0))
   
   !dxdybx
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*b_inx0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdybx(0,0,0))
    
     !dxdzbx
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*b_inx0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzbx(0,0,0))
    
   ! DYDYbX
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*b_inx0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyvb(0,0,0))
    
     ! DYDzbX
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*b_inx0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzbx(0,0,0))
    
     ! DzDzbX
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*b_inx0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzvbx(0,0,0))
   
   ! FINISHED SECOND ORDER  BX TERMS 
   
! SECOND ORDER  by TERMS DXDXBY,DXDYBY,DXDZBY, DYDYBY,DYDZBY, DZDZBY 

!dxdxby
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*b_iny0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxby(0,0,0))
   
   !dxdyvy
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*b_iny0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdyby(0,0,0))
    
     !dxdzby
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*b_iny0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzby(0,0,0))
   ! DYDYby
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*b_iny0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydyby(0,0,0))
    
     ! DYDzbY
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*b_iny0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzby(0,0,0))
    
     ! DzDzbY
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*b_iny0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzby(0,0,0))
    
! FINISHED SECOND ORDER bY TERMS

! SECOND ORDER  BZ TERMS DXDXBZ,DXDYBZ,DXDZBZ, DYDYBz,DYDZBz, DZDZBz
!dxdxbz
  DO i=0,nkx0-1
    DO i=0,nkx0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kxgrid(i)*b_inz0(i,:,:) ! there is  two i's in the 
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdxbz(0,0,0))
   
   !dxdyvz
  DO i=0,nkx0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kygrid(j)*b_inz0(i,j,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdybz(0,0,0))
    
     !dxdzbz
  DO i=0,nkx0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kxgrid(i)*kzgrid(k)*b_inz0(i,:,k)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxdzbz(0,0,0))
   ! DYDYbz
  DO i=0,nky0-1
    DO j=0,nky0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kygrid(j)*b_inz0(:,j,:)  ! towo y grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydybz(0,0,0))
    
     ! DYDzbz
  DO i=0,nky0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(i)*kzgrid(k)*b_inz0(:,j,k)  
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dydzbz(0,0,0))
    
     ! DzDzbz
  DO i=0,nkz0-1
    DO j=0,nkz0-1
        temp_small(i,j,:)=i_complex*kygrid(k)*kzgrid(k)*b_inz0(:,:,k)  !two kz grid
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzdzvbz(0,0,0))
    
!finished END SECOND ORDER BZ TERMS

!    completed  ALL SECOND ORDER B TERMS  

! TERMS BX BY BZ

    !dx phi
    !bx
    DO i=0,nkx0-1
                       !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),bx(0,0,0))
    !by
    DO j=0,nky0-1
                               !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),by(0,0,0))
    !bz
    DO i=0,nkz0-1
                                             !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),bz(0,0,0))

!!! TERMS  vx,vy,vz 
!vx
    DO i=0,nkx0-1
                               !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vx(0,0,0))
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vy(0,0,0))
    !vz
    DO i=0,nkz0-1
                      !temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),vz(0,0,0))

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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvx(0,0,0))

    ! dyvx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inx(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvx(0,0,0))

    ! dzvx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inx(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvx(0,0,0))
    
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvy(0,0,0))
 
    ! dyvy
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_iny(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvy(0,0,0))

    ! dzvy
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_iny(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvy(0,0,0))
    
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxvz(0,0,0))
 
    ! dyvz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyvz(0,0,0))

    ! dzvz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*v_inz(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzvz(0,0,0))
    
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxbx(0,0,0))

    ! dybx
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_inx(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dybx(0,0,0))

    ! dzbx
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inx(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzbx(0,0,0))
        
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxby(0,0,0))
 
    ! dyby
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*b_iny(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyby(0,0,0))

    ! dzby
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_iny(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzby(0,0,0))
    
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
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxbz(0,0,0))
 
    ! dybz
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*v_inz(:,j,:)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dybz(0,0,0))

    ! dzbz
    DO k=0,nkz0-1
        temp_small(:,:,k)=i_complex*kzgrid(k)*b_inz(:,:,k)
    END DO
    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!k loop
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dzbz(0,0,0))

! DONE ALL FIRST ORDER BX,BY BZ TERMS ie. dxbx dybx dzbx,  dxby, dyby,dzby,  dxbz,dybz,dzbz

! EQUATION (14) 1=xcomp, 2=ycomp 3=zcomp

P1 = bx*dxvx+by*dyvx+bz*dzvx
Q1=vx*dxbx+vy*dybx+vz*dzbx
R1=bx*dxdybz+by*dydybz+bz*dzdybz-bx*dxdzby-by*dydzby-bz*dzdzby
S1=dybz*dxbx-dzby*dxbx-dxbz*dybx+dzbx*dybx+dxby*dzbx-dybx*dz*bx 

P2 = bx*dxvy+by*dyvy+bz*dzvy
Q2=vx*dxby+vy*dyby+vz*dzby
R2=bx*dxdzbx+by*dydzbx+bz*dzdzbx-bx*dxdxbz-by*dydxbz-bz*dzdxbz
S2=dybz*dxby-dzby*dxby-dxbz*dyby+dzbx*dyby+dxby*dzby-dybx*dz*by 

P3 = bx*dxvz+by*dyvz+bz*dzvz
Q3=vx*dxbz+vy*dybz+vz*dzbz
R3=bx*dxdxby+by*dydxby+bz*dzdxby-bx*dxdybx-by*dydybx-bz*dzdybx
S3=dybz*dxbz-dzby*dxbz-dxbz*dybz+dzbx*dybz+dxby*dzbz-dybx*dz*bz 

! EQUATION (15)
U1= vx*dxvx+vy*dyvx+vz*dzvx
V2 = bx*dxvx+by*dyvx+bz*dzvx
W1= (0.5)*dyb^2  ! TO BE DONE

U2= vx*dxvy+vy*dyvy+vz*dzvy
V2 = bx*dxvy+by*dyvy+bz*dzvy
W2= (0.5)*dyb^2  ! TO BE DONE

U3= vx*dxvz+vy*dyvz+vz*dzvz
V3 = bx*dxvz+by*dyvz+bz*dzvz
W3= (0.5)*dzb^2  ! TO BE DONE

    IF(first_stage) ve_max(1)=maxval(abs(dxphi)) 

                                          !Now dy phi
                                         ! DO j=0,nky0-1
                                         !     temp_small(:,j,:)=i_complex*kygrid(j)*phi_in(:,j,:)
                                         ! END DO

                                          !Add padding for dealiasing
                                        !  temp_big=cmplx(0.0,0.0)
                                        !  DO i=0,nkx0-1
                                           !   temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
                                           !   temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
                                            !  temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
                                            !  temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
                                        !  END DO!i loop

                                         ! CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyphi(0,0,0))

    IF(first_stage) ve_max(2)=maxval(abs(dyphi)) 

                                 ! DO l=lv1,lv2

                                      !dx g
                                     ! DO i=0,nkx0-1
                                     !     temp_small(i,:,:)=i_complex*kxgrid(i)*g_in0(i,:,:,l,h,0)
                                     ! END DO
                                    !  temp_big=cmplx(0.0,0.0)
                                     ! DO i=0,nkx0-1
                                      !    temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
                                       !   temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
                                       !   temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
                                        !  temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
                                     ! END DO!i loop

                                      !CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxg(0,0,0))

                                          !dy g
                                     ! DO j=0,nky0-1
                                       !   temp_small(:,j,:)=i_complex*kygrid(j)*g_in0(:,j,:,l,h,0)
                                      !END DO
                                     ! temp_big=cmplx(0.0,0.0)
                                    !  DO i=0,nkx0-1
                                       !   temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
                                       !   temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
                                        !  temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
                                        !  temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
                                      ! END DO!i loop
                                     ! CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyg(0,0,0))

      
     
        !re-USE dxg to store the product
      ! dxg=dyphi*dxg-dxphi*dyg
      
     eq14x=P1-Q1-R1+S1
     eq14y=P2-Q2-R2+S2
     eq14z=P3-Q3-R3+S3
     
     eq15x= -U1+V1 -W1
     eq15y=-U2+V2-W2
     eq15z=-U3+V3-W3
                                              !inverse FFT to get back to Fourier
                                              !CALL dfftw_execute_dft_r2c(plan_r2c,dxg(0,0,0),temp_big(0,0,0))
    
    CALL dfftw_execute_dft_r2c(plan_r2c,eq14x(0,0,0),temp_big(0,0,0))
    CALL dfftw_execute_dft_r2c(plan_r2c,eq14y(0,0,0),temp_big(0,0,0))
    CALL dfftw_execute_dft_r2c(plan_r2c,eq14z(0,0,0),temp_big(0,0,0))
    
    CALL dfftw_execute_dft_r2c(plan_r2c,eq15x(0,0,0),temp_big(0,0,0))
    CALL dfftw_execute_dft_r2c(plan_r2c,eq15y(0,0,0),temp_big(0,0,0))
    CALL dfftw_execute_dft_r2c(plan_r2c,eq15z(0,0,0),temp_big(0,0,0))
               
                                                 !Now fill in appropriate rhs elements
        DO i=0,nkx0-1
            rhs_out(i,0:hky_ind,0:hkz_ind,l,h,0)=rhs_out(i,0:hky_ind,0:hkz_ind,l,h,0)+&           !kz positive, ky positive
                                       temp_big(i,0:hky_ind,0:hkz_ind)*fft_norm
            rhs_out(i,0:hky_ind,lkz_ind:nkz0-1,l,h,0)=rhs_out(i,0:hky_ind,lkz_ind:nkz0-1,l,h,0)+& !kz negative, ky positive
                                       temp_big(i,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
            rhs_out(i,lky_ind:nky0-1,lkz_ind:nkz0-1,l,h,0)=rhs_out(i,lky_ind:nky0-1,lkz_ind:nkz0-1,l,h,0)+& !kz negative, ky negative
                                       temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
            rhs_out(i,lky_ind:nky0-1,0:hkz_ind,l,h,0)=rhs_out(i,lky_ind:nky0-1,0:hkz_ind,l,h,0)+& !kz positive, ky negative
                                       temp_big(i,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
        END DO!i loop

                                                                          !END DO !v loop
                                                                        !END DO ! h loop

                                                                        ! The results should be hankel transformed back
                                                                      ! IF (hk_on) THEN 
                                                                      !    DO i=0,nkx0-1
                                                                      !        DO j=0,nky0-1
                                                                      !            DO k=lkz1,lkz2
                                                                      !                DO l=lv1,lv2
                                                                      !                    call hankel_transform_c(rhs_out(i,j,k,l,:,0),.false.)    
                                                                      !                END DO
                                                                      !            END DO
                                                                      !        END DO
                                                                      !    END DO
                                                                      ! ENDIF

                                                                         !IF (hankel) call hankel_transform(rhs_out,.false.) 
     
     ! NEED DEALLOCATE ALL                                                                   

                                           ! DEALLOCATE(g_in0)
  DEALLOCATE(temp_small)
  DEALLOCATE(temp_big)
                                            !DEALLOCATE(dxphi)
                                             ! DEALLOCATE(dyphi)
                                             ! DEALLOCATE(dxg)
                                             ! DEALLOCATE(dyg)

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
SUBROUTINE get_rhs_nl_convolution(g_in,phi_in0,rhs_out)
  USE par_mod
  IMPLICIT NONE

  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  COMPLEX, INTENT(in) :: phi_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  COMPLEX :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  COMPLEX, INTENT(inout) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
 
  INTEGER :: i,j,k,l,ip,jp,kp
  REAL :: ckxmax,ckymax,ckzmax
  REAL :: kx,ky,kz,kxp,kyp,kzp,kxpp,kypp,kzpp
  INTEGER :: xi,yi,zi,xip,yip,zip,xipp,yipp,zipp
  LOGICAL :: take_conjgp,take_conjgpp
  COMPLEX :: phi_kp,g_kpp

  IF(np_hank.gt.1) STOP "get_rhs_convolution not yet implemented for np_hank.gt.1"
  IF(np_spec.gt.1) STOP "get_rhs_convolution not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_convolution not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_convolution only suitable for np_kz=1"
  DO k=0,nkz0-1
    phi_in(:,:,k)=J0a(:,:)*phi_in0(:,:,k)
  END DO

  ckxmax=REAL(hkx_ind)*kxmin
  ckymax=REAL(hky_ind)*kymin
  ckzmax=REAL(hkz_ind)*kzmin

!  IF(mype==0) OPEN(unit=200,file='temp.dat',status='unknown')

  DO l=lv1,lv2
    !WRITE(*,*) "l=",l
    DO i=0,hkx_ind      !kx loop
      DO j=-hky_ind,hky_ind    !ky loop
        DO k=-hkz_ind,hkz_ind          !kz loop

           kx=REAL(i)*kxmin
           ky=REAL(j)*kymin
           kz=REAL(k)*kzmin
           CALL get_k_indices(kx,ky,kz,xi,yi,zi,take_conjgp)
           IF(take_conjgp) STOP "Error!"                   

           DO ip=-hkx_ind,hkx_ind   !kxp loop
             DO jp=-hky_ind,hky_ind !kyp loop
               DO kp=-hkz_ind,hkz_ind     !kzp loop

                  kxp=REAL(ip)*kxmin 
                  kyp=REAL(jp)*kymin 
                  kzp=REAL(kp)*kzmin 

                  kxpp=kx-kxp
                  kypp=ky-kyp
                  kzpp=kz-kzp
 
                  !IF(i==2.and.j==0.and.k==0.and.ip==-13.and.jp==-1.and.kp==0.and.l==1) WRITE(*,*) "Before check.",&
                  !           kxp,kyp,kzp,kxpp,kypp,kzpp,ckxmax,ckymax,ckzmax

                  IF((abs(kxpp).le.(ckxmax+0.001)).and.(abs(kypp).le.(ckymax+0.001)).and.(abs(kzpp).le.(ckzmax+0.001))) THEN


                    CALL get_k_indices(kxp,kyp,kzp,xip,yip,zip,take_conjgp)
                    IF(take_conjgp) THEN
                      phi_kp=conjg(phi_in(xip,yip,zip))
                    ELSE
                      phi_kp=phi_in(xip,yip,zip)
                    END IF

                    CALL get_k_indices(kxpp,kypp,kzpp,xipp,yipp,zipp,take_conjgpp)
                    IF(take_conjgpp) THEN
                      g_kpp=conjg(g_in(xipp,yipp,zipp,l,0,0))
                    ELSE
                      g_kpp=g_in(xipp,yipp,zipp,l,0,0)
                    END IF
        
                    rhs_out(xi,yi,zi,l,0,0)=rhs_out(xi,yi,zi,l,0,0)+&
                                (kxp*ky-kx*kyp)*phi_kp*g_kpp

                  !IF(i==2.and.j==0.and.k==0.and.ip==-13.and.jp==-1.and.kp==0.and.l==1) WRITE(*,*) "After check.",&
                  !              (kxp*ky-kx*kyp)*phi_kp*g_kpp

!           IF(mype==0) THEN 
!             IF(i==2.and.j==0.and.k==0.and.(abs((kxp*ky-kx*kyp)*phi_kp*g_kpp).ne.0.0)) THEN
!                  WRITE(200,*) "i,j,k",i,j,k
!                  WRITE(200,*) "xip,yip,zip",xip,yip,zip
!                  WRITE(200,*) "xipp,yipp,zipp",xipp,yipp,zipp
!                  WRITE(200,*) "take_conjgp",take_conjgp
!                  WRITE(200,*) "take_conjgpp",take_conjgpp
!                  WRITE(200,*) "kxp,ky,kx,kyp",kxp,ky,kx,kyp
!                  WRITE(200,*) "C_k,kp",kxp*ky-kx*kyp
!                  WRITE(200,*) "phi_kp",phi_kp
!                  WRITE(200,*) "g_kpp",g_kpp
!                  WRITE(200,*) "(kxp*ky-kx*kyp)*phi_kp*g_kpp",(kxp*ky-kx*kyp)*phi_kp*g_kpp
!             END IF
!           END IF

                  END IF
                  
               END DO
             END DO
           END DO

        END DO
      END DO
    END DO
  END DO


!  IF(mype==0) CLOSE(200)

END SUBROUTINE get_rhs_nl_convolution


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

END MODULE nonlinearity

