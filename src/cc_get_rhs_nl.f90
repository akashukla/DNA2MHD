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
SUBROUTINE initialize_fourier_ae_mu0_2d

  include 'fftw3.f'
  
  !for dealiasing
  nx0_big=3*nkx0
  ny0_big=3*nky0/2
  fft_norm=1.0/(REAL(nx0_big*ny0_big))

  ALLOCATE(g_rbig_2d(0:nx0_big-1,0:ny0_big-1))
  ALLOCATE(g_kbig_2d(0:nx0_big/2,0:ny0_big-1))

  !WRITE(*,*) "making plans"
  CALL dfftw_plan_dft_c2r_2d(plan_c2r,nx0_big,ny0_big,nz0_big,&
                             g_kbig_2d,g_rbig_2d,FFTW_ESTIMATE)
  CALL dfftw_plan_dft_r2c_2d(plan_r2c,nx0_big,ny0_big,nz0_big,&
                             g_rbig_2d,g_kbig_2d,FFTW_ESTIMATE)
  
  DEALLOCATE(g_rbig_2d)
  DEALLOCATE(g_kbig_2d)

  lky_big=ny0_big-hky_ind !Index of minimum (most negative) FILLED ky value for big arrays

  IF(mype==0) WRITE(*,*) "Initializing FFT"
  IF(mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
  IF(mype==0) WRITE(*,*) "nx0_big,ny0_big",nx0_big,ny0_big
  IF(mype==0) WRITE(*,*) "hky_ind,lky_ind",hky_ind,lky_ind
  IF(mype==0) WRITE(*,*) "lky_big",lky_big

  !CALL dfftw_execute_dft_c2r(plan_c2r,tcomp(0,0,0),treal(0,0,0))
  !CALL dfftw_execute_dft_r2c(plan_r2c,treal(0,0,0),tcomp(0,0,0))
  !tcomp=tcomp/REAL(n1*n2*n3)

END SUBROUTINE initialize_fourier_ae_mu0_2d

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

  ALLOCATE(bx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(by(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(bz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(vx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(vz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(dxvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvy(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzvz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(dxbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(dxdybx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydyby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdybz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxdzbx(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dydzby(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dzdzbz(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  ALLOCATE(dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  
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
   
   !$$$$$$$%%%%%%%&&&&&&&&&&& ! SECOND ORDER  VX TERMS DXDXVX,DXDYVX,DXDZVX,  DYDYVX,DYDZVX, DZDZVX
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
   
   !$$$$%%%%%%%%%%%&&&&&&&&&&&&&&& SECOND ORDER  VX TERMS END
   
! SECOND ORDER  VY TERMS DXDXVY,DXDYVY,DXDZVY, DYDYVY,DYDZVY, DZDZVY 
!$$$$$$$$$$#########%%%%%%%%%$$$$$$$$$$$$
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
    
!!!!!!!!!!!!!!$$$$$$$$$%%%%%%%%%%%%%%&&&&&&&&&&&&& END SECOND ORDER VY TERMS

!@@@@@@@@@@%%%%%%%%%%%%%$$$$$$$$$$ SECOND ORDER  VZ TERMS DXDXVZ,DXDYVZ,DXDZVZ, DYDYVz,DYDZVz, DZDZVz
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
    
!$$$$$$$$$$$$###############@@@@@@&&&&&&&&&&&&&&& VZ THINGS END

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

!!! Prerana do vx,vy,vz terms just like bx,by,bz above  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!! PRERANA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    
    !Prerana do
!!  ! &&&&&&&&&&&&&&&&&&&&&&&  Prerana do dxvy dyvy dzvz, 
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
    
    ! &&&&&&&&&&&&&&&&&&&&&&& dxvz dyvz dzvz
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
    
    !&&&&&&&&&&&&&&& DONE ALL  dxvx dyvx dzvx,dxvy dyvy,dzvy, dxvz,dyvz,dzvz
    
    ! NOW FOR bx
     !! Prerana do dxbx dybx dzbx
    !! dxby dyby dzby
    !! dxbz dybz dzbz
    
    !&&&&&&&&&&& dxbx dybx dzbx` 
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
    
    
    ! !&&&&&&&&&&&&&&&&&&&&&& dxby dyby dzby
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
    
    
    !! &&&&&&&&&&&&&&&& dxbz dybz dzbz

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

!&&&&&&&&&&&&&&& DONE ALL  dxbx dybx dzbx,dxbydyby,dzby, dxbz,dybz,dzbz

    IF(first_stage) ve_max(1)=maxval(abs(dxphi)) 

    !Now dy phi
    DO j=0,nky0-1
        temp_small(:,j,:)=i_complex*kygrid(j)*phi_in(:,j,:)
    END DO

    !Add padding for dealiasing
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
        temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
        temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
        temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
        temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    END DO!i loop
 
    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyphi(0,0,0))

    IF(first_stage) ve_max(2)=maxval(abs(dyphi)) 

    DO l=lv1,lv2

        !dx g
        DO i=0,nkx0-1
            temp_small(i,:,:)=i_complex*kxgrid(i)*g_in0(i,:,:,l,h,0)
        END DO
        temp_big=cmplx(0.0,0.0)
        DO i=0,nkx0-1
            temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
            temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
            temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
            temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
        END DO!i loop

        CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxg(0,0,0))
    
        !dy g
        DO j=0,nky0-1
            temp_small(:,j,:)=i_complex*kygrid(j)*g_in0(:,j,:,l,h,0)
        END DO
        temp_big=cmplx(0.0,0.0)
        DO i=0,nkx0-1
            temp_big(i,0:hky_ind,0:hkz_ind)=temp_small(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
            temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=temp_small(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
            temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
            temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=temp_small(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
        END DO!i loop

        CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyg(0,0,0))

        !re-USE dxg to store the product
        dxg=dyphi*dxg-dxphi*dyg
        !inverse FFT to get back to Fourier
        CALL dfftw_execute_dft_r2c(plan_r2c,dxg(0,0,0),temp_big(0,0,0)) 

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

    END DO !v loop
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
  
  IF (hankel) call hankel_transform(rhs_out,.false.)    

  DEALLOCATE(g_in0)
  DEALLOCATE(temp_small)
  DEALLOCATE(temp_big)
  DEALLOCATE(dxphi)
  DEALLOCATE(dyphi)
  DEALLOCATE(dxg)
  DEALLOCATE(dyg)

END SUBROUTINE get_rhs_nl1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the fastest version
!Note: only for np_kz=1
SUBROUTINE get_rhs_nl2(b_in,v_in,rhs_out_b,rhs_out_v)
  USE par_mod
  IMPLICIT NONE
  include 'fftw3.f'


  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  COMPLEX, INTENT(inout) :: phi_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  !COMPLEX :: phi_in(0:nkx0-1,0:nky0-1,0:nkz0-1)
  COMPLEX, INTENT(inout) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
 
  INTEGER :: i,j,l,k

  IF(np_hank.gt.1) STOP "get_rhs_nl2 not yet implemented for np_hank.gt.1"
  IF(np_spec.gt.1) STOP "get_rhs_nl2 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl2 not yet implemented for np_kz.gt.1"

  ALLOCATE(temp_small(0:nkx0-1,0:nky0-1,0:nkz0-1))
  ALLOCATE(temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))
  ALLOCATE(dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1))

  !Gyroaverage
  DO k=0,nkz0-1
    phi_in0(0:nkx0-1,0:nky0-1,k)=J0a(0:nkx0-1,0:nky0-1)*phi_in0(0:nkx0-1,0:nky0-1,k)
  END DO

  !dx phi
!  DO i=0,nkx0-1
!    temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
!  END DO

  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  DO i=0,nkx0-1
    temp_big(i,0:hky_ind,0:hkz_ind)=i_complex*kxgrid(i)*phi_in0(i,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
    temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=i_complex*kxgrid(i)*phi_in0(i,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
    temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=i_complex*kxgrid(i)*phi_in0(i,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
    temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=i_complex*kxgrid(i)*phi_in0(i,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
  END DO!k loop
  
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxphi(0,0,0))

  IF(first_stage) ve_max(1)=maxval(abs(dxphi)) 

  !Now dy phi
  DO j=0,nky0-1
    temp_small(:,j,:)=i_complex*kygrid(j)*phi_in0(:,j,:)
  END DO

  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
  !DO i=0,nkx0-1
    temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
    temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
    temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
    temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
  !END DO!i loop
 
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyphi(0,0,0))

  IF(first_stage) ve_max(2)=maxval(abs(dyphi)) 

  DO l=lv1,lv2

    !dx g
    !DO i=0,nkx0-1
    !  temp_small(i,:,:)=i_complex*kxgrid(i)*g_in(i,:,:,l)
    !END DO
    temp_big=cmplx(0.0,0.0)
    DO i=0,nkx0-1
      temp_big(i,0:hky_ind,0:hkz_ind)=i_complex*kxgrid(i)*g_in(i,0:hky_ind,0:hkz_ind,l,0,0)    !kz positive, ky positive    
      temp_big(i,0:hky_ind,lkz_big:nz0_big-1)=i_complex*kxgrid(i)*g_in(i,0:hky_ind,lkz_ind:nkz0-1,l,0,0) !kz negative, ky positive
      temp_big(i,lky_big:ny0_big-1,lkz_big:nz0_big-1)=i_complex*kxgrid(i)*g_in(i,lky_ind:nky0-1,lkz_ind:nkz0-1,l,0,0) !kz negative, ky negative
      temp_big(i,lky_big:ny0_big-1,0:hkz_ind)=i_complex*kxgrid(i)*g_in(i,lky_ind:nky0-1,0:hkz_ind,l,0,0) !kz positive, ky negative
    END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxg(0,0,0))
    
    !dy g
    DO j=0,nky0-1
      temp_small(:,j,:)=i_complex*kygrid(j)*g_in(:,j,:,l,0,0)
    END DO
    temp_big=cmplx(0.0,0.0)
    !DO i=0,nkx0-1
      temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
      temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
      temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
      temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    !END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyg(0,0,0))

    !re-USE dxg to store the product
    dxg=dyphi*dxg-dxphi*dyg
    !inverse FFT to get back to Fourier
    CALL dfftw_execute_dft_r2c(plan_r2c,dxg(0,0,0),temp_big(0,0,0)) 

    !Now fill in appropriate rhs elements
    !DO i=0,nkx0-1
      rhs_out(0:nkx0-1,0:hky_ind,0:hkz_ind,l,0,0)=rhs_out(0:nkx0-1,0:hky_ind,0:hkz_ind,l,0,0)+&           !kz positive, ky positive
                                       temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)*fft_norm
      rhs_out(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1,l,0,0)=rhs_out(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1,l,0,0)+& !kz negative, ky positive
                                       temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
      rhs_out(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1,l,0,0)=rhs_out(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1,l,0,0)+& !kz negative, ky negative
                                       temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
      rhs_out(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind,l,0,0)=rhs_out(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind,l,0,0)+& !kz positive, ky negative
                                       temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
    !END DO!i loop

  END DO

  !de-gyroaverage
  DO k=0,nkz0-1
    phi_in0(0:nkx0-1,0:nky0-1,k)=phi_in0(0:nkx0-1,0:nky0-1,k)/J0a(0:nkx0-1,0:nky0-1)
  END DO

  DEALLOCATE(temp_small)
  DEALLOCATE(temp_big)
  DEALLOCATE(dxphi)
  DEALLOCATE(dyphi)
  DEALLOCATE(dxg)
  DEALLOCATE(dyg)

END SUBROUTINE get_rhs_nl2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl3                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_nl3(g_in,phi_in0,rhs_out)
  USE par_mod
  include 'fftw3.f'

  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  COMPLEX, INTENT(in) :: phi_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  COMPLEX :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  COMPLEX, INTENT(inout) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  COMPLEX :: temp_small(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  COMPLEX :: temp_big(0:nx0_big/2,0:ny0_big-1,0:nz0_big-1)
  REAL :: dxphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
  REAL :: dyphi(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
  REAL :: dxg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
  REAL :: dyg(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
!  REAL :: nl_tot_real(0:nx0_big-1,0:ny0_big-1,0:nz0_big-1)
 
  INTEGER :: i,j,l,k

  IF(np_hank.gt.1) STOP "get_rhs_nl3 not yet implemented for np_hank.gt.1"
  IF(np_spec.gt.1) STOP "get_rhs_nl3 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl3 not yet implemented for np_kz.gt.1"

  IF(np_kz.ne.1) STOP "get_rhs_nl3 only suitable for np_kz=1"
  !IF(mype==0.and.first_stage) WRITE(*,*) "mype,itime,rhs_lin",mype,itime,abs(sum(sum(sum(sum(rhs_out,1),1),1),1))
  DO k=0,nkz0-1
    phi_in(:,:,k)=J0a(:,:)*phi_in0(:,:,k)
  END DO

  !dx phi
  DO i=0,nkx0-1
    temp_small(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
  END DO

  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
!  DO i=0,nkx0-1
    temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
    temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
    temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
    temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
!  END DO!k loop
  
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxphi(0,0,0))

  IF(first_stage) ve_max(1)=maxval(abs(dxphi)) 

  !Now dy phi
  DO j=0,nky0-1
    temp_small(:,j,:)=i_complex*kygrid(j)*phi_in(:,j,:)
  END DO

  !Add padding for dealiasing
  temp_big=cmplx(0.0,0.0)
!  DO i=0,nkx0-1
    temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
    temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
    temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
    temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
!  END DO!i loop
 
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyphi(0,0,0))

  IF(first_stage) ve_max(2)=maxval(abs(dyphi)) 

  DO l=lv1,lv2

    !dx g
    DO i=0,nkx0-1
      temp_small(i,:,:)=i_complex*kxgrid(i)*g_in(i,:,:,l,0,0)
    END DO
    temp_big=cmplx(0.0,0.0)
    !DO i=0,nkx0-1
      temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
      temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
      temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
      temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    !END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dxg(0,0,0))
    
    !dy g
    DO j=0,nky0-1
      temp_small(:,j,:)=i_complex*kygrid(j)*g_in(:,j,:,l,0,0)
    END DO
    temp_big=cmplx(0.0,0.0)
    !DO i=0,nkx0-1
      temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)=temp_small(0:nkx0-1,0:hky_ind,0:hkz_ind)    !kz positive, ky positive    
      temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1) !kz negative, ky positive
      temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)=temp_small(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1) !kz negative, ky negative
      temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)=temp_small(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind) !kz positive, ky negative
    !END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big(0,0,0),dyg(0,0,0))

    !re-USE dxg to store the product
    dxg=dyphi*dxg-dxphi*dyg
    !inverse FFT to get back to Fourier
    CALL dfftw_execute_dft_r2c(plan_r2c,dxg(0,0,0),temp_big(0,0,0)) 

    !Now fill in appropriate rhs elements
    !DO i=0,nkx0-1
      rhs_out(0:nkx0-1,0:hky_ind,0:hkz_ind,l,0,0)=rhs_out(0:nkx0-1,0:hky_ind,0:hkz_ind,l,0,0)+&           !kz positive, ky positive
                                       temp_big(0:nkx0-1,0:hky_ind,0:hkz_ind)*fft_norm
      rhs_out(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1,l,0,0)=rhs_out(0:nkx0-1,0:hky_ind,lkz_ind:nkz0-1,l,0,0)+& !kz negative, ky positive
                                       temp_big(0:nkx0-1,0:hky_ind,lkz_big:nz0_big-1)*fft_norm
      rhs_out(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1,l,0,0)=rhs_out(0:nkx0-1,lky_ind:nky0-1,lkz_ind:nkz0-1,l,0,0)+& !kz negative, ky negative
                                       temp_big(0:nkx0-1,lky_big:ny0_big-1,lkz_big:nz0_big-1)*fft_norm
      rhs_out(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind,l,0,0)=rhs_out(0:nkx0-1,lky_ind:nky0-1,0:hkz_ind,l,0,0)+& !kz positive, ky negative
                                       temp_big(0:nkx0-1,lky_big:ny0_big-1,0:hkz_ind)*fft_norm
    !END DO!i loop

  END DO

END SUBROUTINE get_rhs_nl3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 get_rhs_nl4                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!For nkz0=1
SUBROUTINE get_rhs_nl4(g_in,phi_in0,rhs_out)
  USE par_mod
  IMPLICIT NONE
  include 'fftw3.f'


  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
  COMPLEX, INTENT(inout) :: phi_in0(0:nkx0-1,0:nky0-1,lkz1:lkz2)
  !COMPLEX :: phi_in(0:nkx0-1,0:nky0-1,0:nkz0-1)
  COMPLEX, INTENT(inout) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
 
  INTEGER :: i,j,l,k

  IF(np_hank.gt.1) STOP "get_rhs_nl2 not yet implemented for np_hank.gt.1"
  IF(np_spec.gt.1) STOP "get_rhs_nl2 not yet implemented for np_spec.gt.1"
  IF(np_kz.gt.1) STOP "get_rhs_nl2 not yet implemented for np_kz.gt.1"

  ALLOCATE(temp_small_2d(0:nkx0-1,0:nky0-1))
  ALLOCATE(temp_big_2d(0:nx0_big/2,0:ny0_big-1))
  ALLOCATE(dxphi_2d(0:nx0_big-1,0:ny0_big-1))
  ALLOCATE(dyphi_2d(0:nx0_big-1,0:ny0_big-1))
  ALLOCATE(dxg_2d(0:nx0_big-1,0:ny0_big-1))
  ALLOCATE(dyg_2d(0:nx0_big-1,0:ny0_big-1))

  !Gyroaverage
  DO k=0,nkz0-1
    phi_in0(0:nkx0-1,0:nky0-1,k)=J0a(0:nkx0-1,0:nky0-1)*phi_in0(0:nkx0-1,0:nky0-1,k)
  END DO

  !dx phi
!  DO i=0,nkx0-1
!    temp_small_2d(i,:,:)=i_complex*kxgrid(i)*phi_in(i,:,:)
!  END DO

  !Add padding for dealiasing
  temp_big_2d=cmplx(0.0,0.0)
  DO i=0,nkx0-1
    temp_big_2d(i,0:hky_ind)=i_complex*kxgrid(i)*phi_in0(i,0:hky_ind,0)    !ky positive    
    temp_big_2d(i,lky_big:ny0_big-1)=i_complex*kxgrid(i)*phi_in0(i,lky_ind:nky0-1,0) !ky negative
  END DO!k loop
  
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big_2d(0,0),dxphi_2d(0,0))

  IF(first_stage) ve_max(1)=maxval(abs(dxphi_2d)) 

  !Now dy phi
  DO j=0,nky0-1
    temp_small_2d(:,j)=i_complex*kygrid(j)*phi_in0(:,j,0)
  END DO

  !Add padding for dealiasing
  temp_big_2d=cmplx(0.0,0.0)
  !DO i=0,nkx0-1
    temp_big_2d(0:nkx0-1,0:hky_ind)=temp_small_2d(0:nkx0-1,0:hky_ind)    !kz positive, ky positive    
    temp_big_2d(0:nkx0-1,lky_big:ny0_big-1)=temp_small_2d(0:nkx0-1,lky_ind:nky0-1) !kz negative, ky negative
  !END DO!i loop
 
  CALL dfftw_execute_dft_c2r(plan_c2r,temp_big_2d(0,0),dyphi_2d(0,0))

  IF(first_stage) ve_max(2)=maxval(abs(dyphi_2d)) 

  DO l=lv1,lv2

    !dx g
    !DO i=0,nkx0-1
    !  temp_small(i,:,:)=i_complex*kxgrid(i)*g_in(i,:,:,l)
    !END DO
    temp_big_2d=cmplx(0.0,0.0)
    DO i=0,nkx0-1
      temp_big_2d(i,0:hky_ind)=i_complex*kxgrid(i)*g_in(i,0:hky_ind,0,l,0,0)    ! ky positive    
      temp_big_2d(i,lky_big:ny0_big-1)=i_complex*kxgrid(i)*g_in(i,lky_ind:nky0-1,0,l,0,0) ! ky negative
    END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big_2d(0,0),dxg_2d(0,0))
    
    !dy g
    DO j=0,nky0-1
      temp_small_2d(:,j)=i_complex*kygrid(j)*g_in(:,j,0,l,0,0)
    END DO
    temp_big_2d=cmplx(0.0,0.0)
    !DO i=0,nkx0-1
      temp_big_2d(0:nkx0-1,0:hky_ind)=temp_small_2d(0:nkx0-1,0:hky_ind)    !ky positive    
      temp_big_2d(0:nkx0-1,lky_big:ny0_big-1)=temp_small_2d(0:nkx0-1,lky_ind:nky0-1) ! ky negative
    !END DO!i loop

    CALL dfftw_execute_dft_c2r(plan_c2r,temp_big_2d(0,0),dyg_2d(0,0))

    !re-USE dxg_2d to store the product
    dxg_2d=dyphi_2d*dxg_2d-dxphi_2d*dyg_2d
    !inverse FFT to get back to Fourier
    CALL dfftw_execute_dft_r2c(plan_r2c,dxg_2d(0,0),temp_big_2d(0,0)) 

    !Now fill in appropriate rhs elements
    !DO i=0,nkx0-1
      rhs_out(0:nkx0-1,0:hky_ind,0,l,0,0)=rhs_out(0:nkx0-1,0:hky_ind,0,l,0,0)+&           ! ky positive
                                       temp_big_2d(0:nkx0-1,0:hky_ind)*fft_norm
      rhs_out(0:nkx0-1,lky_ind:nky0-1,0,l,0,0)=rhs_out(0:nkx0-1,lky_ind:nky0-1,0,l,0,0)+& !kz negative, ky negative
                                       temp_big_2d(0:nkx0-1,lky_big:ny0_big-1)*fft_norm
    !END DO!i loop

  END DO

  !de-gyroaverage
  DO k=0,nkz0-1
    phi_in0(0:nkx0-1,0:nky0-1,k)=phi_in0(0:nkx0-1,0:nky0-1,k)/J0a(0:nkx0-1,0:nky0-1)
  END DO

  DEALLOCATE(temp_small_2d)
  DEALLOCATE(temp_big_2d)
  DEALLOCATE(dxphi_2d)
  DEALLOCATE(dyphi_2d)
  DEALLOCATE(dxg_2d)
  DEALLOCATE(dyg_2d)

END SUBROUTINE get_rhs_nl4



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

