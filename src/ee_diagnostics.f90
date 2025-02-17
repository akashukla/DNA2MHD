!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                            ee_diagnostics.f90                             !!
!!                                                                           !!
!!  diagnostics                                                              !!
!!  -- initialize_diagnostics                                                !!
!!  -- finalize_diagnostics                                                  !!
!!  -- diag                                                                  !!
!!  --                                                                       !!
!!  -- ...                                                                   !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               diagnostics                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE diagnostics
  USE par_mod
  USE mpi
  USE nonlinearity
  USE linear_rhs
  !USE flr_effects, only: J0a
  !USE hk_effects
  !USE field_solver, only: get_phi
  !USE Gyro_LES
  IMPLICIT NONE

  PUBLIC :: initialize_diagnostics, finalize_diagnostics, diag
     !output_data, nl_test,

  PRIVATE

  INTEGER :: ffm_handle,omega_handle,glast_handle,ev_handle,en_handle,enspec_handle,xi_handle,mode_handle,&
    herm_handle,eshells_handle,fmom3d_handle,energy3d_handle, gk_nkyout,&
    gknl_nkyout
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gk_indices
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gk_handle
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gknl_indices
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gknl_handle
  !!!!BENCHMARKING!!!!
  !!!!BENCHMARKING!!!!
  !INTEGER, ALLOCATABLE, DIMENSION(:) :: gknl_temphandle
  !!!!BENCHMARKING!!!!
  !!!!BENCHMARKING!!!!
  REAL :: energy_last,time_last
  REAL :: real_err
  CHARACTER(len=4) :: char_nv0
  INTEGER :: initial_wallclock
  REAL :: current_wallclock
  LOGICAL :: file_exists

  !NLT diagnostics
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: NLT
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: NLT_n
  REAL, ALLOCATABLE, DIMENSION(:) :: shell_bounds
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: nltk

  REAL :: magbound,canbound
  
  INTEGER :: num_shells
  INTEGER :: ikx_minmax
  INTEGER :: iky_minmax
  INTEGER :: ikz_minmax
  INTEGER :: shell_handle
  INTEGER :: shell_handle_n(11)
  INTEGER :: nlt_handle
  INTEGER :: nlt_status
  LOGICAL :: shells_initialized
  !nlt_triple diagnostics
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: AVP
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: WVORT
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: OSPEC
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WRITESPEC
  REAL :: ham0,ham1
  INTEGER :: nlt3_handle 
  !NLT Testing
  !REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: NLT_bak
  INTEGER :: test_handle1,test_handle2
  INTEGER :: ierr
  INTEGER(4) :: i,j,k

  !eshells


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            initialize_diagnostics                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_diagnostics
  IMPLICIT NONE
  !LOGICAL :: file_exists

  shells_initialized=.false.
  CALL get_io_number
  ffm_handle=io_number
  CALL get_io_number
  omega_handle=io_number
  CALL get_io_number
  glast_handle=io_number
  CALL get_io_number
  ev_handle=io_number

  if(np_kz.gt.1) STOP "Must implement diagnostics for kz parallelization."


  IF(istep_energy.gt.0) THEN
    !OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',status='unknown')
    IF (mype.eq.0) THEN
       CALL get_io_number
       en_handle = io_number
       print *, "Energy Handle",en_handle
    IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/energy_out.dat',exist=file_exists)
      IF(file_exists) THEN
        OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',form='unformatted', status='unknown',access='stream',position='append')
      ELSE
        OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',form='unformatted', status='REPLACE',access='stream')
        ! WRITE(en_handle,*) "#time,entropy,phi^2 energy,dE/dt total,flux,coll,hcoll,hyps,N.L.,hyp_conv,dE/dt"
      END IF
    ELSE
      OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',form='unformatted', status='REPLACE',access='stream')
      ! WRITE(en_handle,*) "#time,entropy,phi^2 energy,dE/dt total,flux,coll,hcoll,hyps,N.L.,hyp_conv,dE/dt"
     ! WRITE(en_handle,*) &
      !   "#time,Hamiltonian,restvty.,viscsty."
   END IF
   ENDIF
    energy_last=0.0
    time_last=time
    ALLOCATE(AVP(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2))
    ALLOCATE(WVORT(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2))
 END IF
  IF(istep_energyspec.gt.0.and.mype==0) THEN
    CALL get_io_number
    enspec_handle=io_number
    print *, "Enspec Handle",enspec_handle
    CALL get_io_number
    mode_handle = io_number
    print *, "Mode Handle",io_number
    !OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',status='unknown') 
    IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/energyspec_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/mode_out.dat',exist=file_exists)
      IF(file_exists) THEN
        OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='unknown',access='stream',position='append')
        OPEN(unit=mode_handle,file=trim(diagdir)//'/mode_out.dat',form='unformatted', status='unknown',access='stream',position='append')
      ELSE
        OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=mode_handle,file=trim(diagdir)//'/mode_out.dat',form='unformatted', status='REPLACE',access='stream')
      END IF
    ELSE
      OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='REPLACE',access='stream')
      OPEN(unit=mode_handle,file=trim(diagdir)//'/mode_out.dat',form='unformatted', status='REPLACE',access='stream')
    END IF
    ALLOCATE(LW(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(LC(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(RW(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(RC(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(OSPEC(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(WRITESPEC(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
 END IF

END SUBROUTINE initialize_diagnostics


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            finalize_diagnostics                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE finalize_diagnostics
  IMPLICIT NONE

  !IF(verbose) WRITE(*,*) "Finalizing checkpoint.",mype
  !WRITE(*,*) "Finalizing diagnostics.",mype
  IF(checkpoint_write) CALL checkpoint_out(2)
  IF(mype==0.and.verbose) WRITE(*,*) "Done writing checkpoint."
  IF(mype==0) CLOSE(ffm_handle)

  IF(istep_energy.gt.0) THEN
    if (mype.eq.0) CLOSE(en_handle)
    DEALLOCATE(AVP)
    DEALLOCATE(WVORT)
  END IF
  IF(istep_energyspec.gt.0.and.mype==0) THEN
    CLOSE(enspec_handle)
    CLOSE(mode_handle)
    DEALLOCATE(LW)
    DEALLOCATE(LC)
    DEALLOCATE(RW)
    DEALLOCATE(RC)
    DEALLOCATE(OSPEC)
    DEALLOCATE(WRITESPEC)
  END IF
  IF(verbose.and.(mype==0)) print *, 'Closed energy handles'
 
  IF((verbose.and.plot_nls).and.(mype.eq.0)) print *, 'Closed nl handles'

END SUBROUTINE finalize_diagnostics


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 diag                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Main routine for diagnostics.  Note: IF you change something here, you must
!!  make the corresponding modification in the python diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE diag
  IMPLICIT NONE

  !COMPLEX :: phi_kz(0:nkz0-1)
  REAL :: phi_tot,flux_tot
  COMPLEX :: gam
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: energy3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: energy3d_temp
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mom3d

  if (itime.eq.0) CALL bound_hels

     IF((istep_energy.ne.0).and.(MOD(itime,istep_energy)==0))THEN
         IF(verbose) WRITE(*,*) "Starting energy diag.",mype
         IF (mype.eq.0) WRITE(en_handle) time
         CALL hmhdhmtn(0)
         if (verbose) write(*,*) "Found Hamiltonian",mype
         CALL mag_helicity()
         CALL cross_helicity()
         CALL hmhdhmtn(1)
         CALL hmhdhmtn(2)
         CALL mode_energy
         if(mype.eq.0) WRITE(en_handle) magbound
         if (mype.eq.0) WRITE(en_handle) canbound
         if (mype.eq.0) WRITE(en_handle) mhelcorr
         if (verbose) write(*,*) "Found Helicities",mype
     END IF

   IF (.false.) THEN
     IF((istep_energyspec.ne.0).and.(mype.eq.0).and.(MOD(itime,istep_energyspec)==0)) THEN
         ! IF(verbose) WRITE(*,*) "Starting energyspec diag.",mype
         WRITE(enspec_handle) time
         CALL en_spec()
         ! IF(verbose) WRITE(*,*) "Done with energyspec diag.",mype
      END IF
   ENDIF
   
  IF((istep_schpt.ne.0).and.(MOD(itime,istep_schpt)==0)) THEN
      ! IF(verbose) WRITE(*,*) "Writing s_checkpoint.",mype
      CALL checkpoint_out(1)
      ! IF(verbose) WRITE(*,*) "Done writing s_checkpoint.",mype
   END IF

  IF((istep_gout.ne.0).and.(MOD(itime,istep_gout)==0)) THEN

      ! IF(verbose) WRITE(*,*) "Starting vbout diag.",mype

      CALL checkpoint_out(4)

      ! IF(itime==0.and.gout_nl) THEN
      !   CALL checkpoint_out(5)
      ! ELSE IF (gout_nl) THEN
      !   CALL checkpoint_out(6)
      ! END IF

      ! IF(verbose) WRITE(*,*) "Done with vbout diag.",mype

   END IF

END SUBROUTINE diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             get_indices_from_ks                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_indices_from_ks(kx_in,ky_in,kz_in,ikx,iky,ikz,take_conjg)
    IMPLICIT NONE

    REAL, INTENT(in) :: kx_in,ky_in,kz_in
    INTEGER, INTENT(out) :: ikx,iky,ikz
    LOGICAL, INTENT(out) :: take_conjg
    REAL :: kx0,ky0,kz0

    take_conjg=.false.
    kx0=kx_in
    ky0=ky_in 
    kz0=kz_in 

    ikx=nint(kx0/kxmin)
    IF(kx0.lt.0.0) THEN
      take_conjg=.true.
      ikx=-1*ikx
      ky0=-1.0*ky0
      kz0=-1.0*kz0
    END IF

    IF(ky0.ge.0.0) THEN
      iky=nint(ky0/kymin)
    ELSE
      iky=nint(ky0/kymin)+nky0
    END IF

    IF(kz0.ge.0.0) THEN
      ikz=nint(kz0/kzmin)
    ELSE
      ikz=nint(kz0/kzmin)+nkz0
    END IF

  END SUBROUTINE get_indices_from_ks
  
subroutine hmhdhmtn(opt) 

implicit none
real :: ham
integer :: opt ! 0 Total Energy, 1 Kinetic Energy, 2 Magnetic Energy
real :: hamsm(2),hams(2)

!! Computes the Hall MHD Hamiltonian 8 pi^3 (sum(v_k^2 + b_k^2) + 2b_z0)/2 

hamsm = 0.0

if ((opt.eq.0).or.(opt.eq.1)) then
   !if (splitx) hamsm(1) = hamsm(1) + sum(abs(v_1(1:nkx0-1,:,:,:))**2)
   !if (splitx) hamsm(2) = hamsm(2) + 0.5 * sum(abs(v_1(0,:,:,:))**2)
   hamsm(1) = hamsm(1) + sum(abs(v_1(1:nkx0-1,:,:,:))**2)
   hamsm(2) = hamsm(2) + 0.5 * sum(abs(v_1(0,:,:,:))**2)

endif
if ((opt.eq.0).or.(opt.eq.2)) then
   !if (splitx) hamsm(1) = hamsm(1) + sum(abs(b_1(1:nkx0-1,:,:,:))**2)
   !if (splitx) hamsm(2) = hamsm(2) + 0.5 * sum(abs(b_1(0,:,:,:))**2)

   hamsm(1) = hamsm(1) + sum(abs(b_1(1:nkx0-1,:,:,:))**2)
   hamsm(2) = hamsm(2) + 0.5 * sum(abs(b_1(0,:,:,:))**2)
   if (mype.eq.0) hamsm(2) = hamsm(2) + real(b_1(0,0,0,2))
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(hamsm,hams,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

ham = (hams(1)+hams(2)) * (8*(pi**3))
if (verbose) print *, mype,"Hamiltonian",(hamsm(1)+hamsm(2)) * (8*(pi**3))
if (verbose) print *, mype,"Total Ham",ham

if (mype.eq.0) WRITE(en_handle) ham

if (opt.eq.0) then
   if (itime.eq.0) ham0 = ham

!   if (verbose.and.(mype.eq.0)) print *, "Initial Hamiltonian",ham

if (mod(itime,istep_energyspec).eq.0.and.itime.gt.0) then
   ham1 = ham
   write(*,"(A20,E10.3)") "Energy Fractional Change",(ham1-ham0)/(ham0)
   ham0 = ham1
endif
endif


end subroutine hmhdhmtn

subroutine vec_potential

implicit none
integer :: i,j,k,ind

! Ak = (i k x Bk)/k^2
AVP = cmplx(0.0,0.0)

do i = 0,nkx0-1
  do j = 0,ny0_big-1
     do k = lkz1,lkz2

          AVP(i,j,k,0) = cmplx(0.0,1.0) * (kygrid(j)*b_1(i,j,k,2)-kzgrid(k)*b_1(i,j,k,1))/(kmags(i,j,k)**2)
          AVP(i,j,k,1) = cmplx(0.0,1.0) * (kzgrid(k)*b_1(i,j,k,0)-kxgrid(i)*b_1(i,j,k,2))/(kmags(i,j,k)**2)
          AVP(i,j,k,2) = cmplx(0.0,1.0) * (kxgrid(i)*b_1(i,j,k,1)-kygrid(j)*b_1(i,j,k,0))/(kmags(i,j,k)**2)

    enddo
   enddo
enddo

if (mype.eq.0) AVP(0,0,0,:) = cmplx(0.0,0.0)


end subroutine vec_potential

SUBROUTINE mag_helicity 

implicit none
real :: maghel
integer :: i,j,k,ind
real :: mhs(2),mhsm(2)

mhsm = 0.0

CALL vec_potential()

! if (verbose) print *, mype,"Vector Potential"

mhsm(1) = mhsm(1) + 2.0*sum(real(AVP(1:nkx0-1,:,:,:)*conjg(b_1(1:nkx0-1,:,:,:))))
mhsm(2) = mhsm(2) + sum(AVP(0,:,:,:)*conjg(b_1(0,:,:,:)))

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(mhsm,mhs,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

! if (verbose) print *, mype,"Helicity Allreduce"
maghel = (mhs(1) + mhs(2)) * (2.0 * pi)**3.0

if (mype.eq.0) WRITE(en_handle) maghel

if (verbose) print *,mype,"Magnetic Helicity",(maghel+mhelcorr)

END SUBROUTINE mag_helicity

subroutine vorticity

implicit none
integer :: ind,i,j,k

WVORT = cmplx(0.0,0.0)

  do i = 0,nkx0-1
    do j = 0,ny0_big-1
      do k = lkz1,lkz2
          WVORT(i,j,k,0) = kygrid(j) * v_1(i,j,k,2) - kzgrid(k) * v_1(i,j,k,1)
          WVORT(i,j,k,1) = kzgrid(k) * v_1(i,j,k,0) - kxgrid(i) * v_1(i,j,k,2)
          WVORT(i,j,k,2) = kxgrid(i) * v_1(i,j,k,1) - kygrid(j) * v_1(i,j,k,0)
      enddo
    enddo
   enddo

WVORT = cmplx(0.0,1.0) * WVORT

! if ((verbose).and.(itime.lt.200)) write(*,*) 'w',w

end subroutine vorticity

subroutine cross_helicity

  implicit none
  real :: crosshel,crosshel1,crosshel2
  real :: mh,vb,vw
  integer :: i,j,k,ind
  real :: chs(2),chsm(2)
  
  chsm = 0.0
  CALL vec_potential()
  CALL vorticity()
  
  chsm(1) = 2.0 * sum(real((AVP(1:nkx0-1,:,:,:)+hall*v_1(1:nkx0-1,:,:,:))*conjg(b_1(1:nkx0-1,:,:,:)+hall*WVORT(1:nkx0-1,:,:,:))))
  chsm(2) = real(sum((AVP(0,:,:,:)+hall*v_1(0,:,:,:))*conjg(b_1(0,:,:,:)+hall*WVORT(0,:,:,:))))
  if (mype.eq.0) chsm(2) = chsm(2) + hall*v_1(0,0,0,2)
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(chsm,chs,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   
  crosshel = (chs(1)+chs(2)) * (2.0*pi)**3.0
   
  if (mype.eq.0) WRITE(en_handle) crosshel
  if (verbose) print *, mype,"Cross Helicity",(crosshel+mhelcorr)

end subroutine cross_helicity

subroutine en_spec() 

  implicit none

  integer :: test1, test2

! This doesn't need to be adjusted for 2x only non ksplit = 0 modes

! Use mode tools to get steady state estimates
LW = (8*pi**3)* (abs(v_1(:,:,:,0))**2 + abs(v_1(:,:,:,1))**2+abs(v_1(:,:,:,2))**2)+ &
     (8*pi**3)* (abs(b_1(:,:,:,0))**2+ abs(b_1(:,:,:,1))**2+abs(b_1(:,:,:,2))**2)

! Rate of Energy Flow Constant in Inertial Range

test1 = nint(3*kmax*force_frac/kxmin)+2
test2 = nint(7*kmax*force_frac/kxmin)+2

! CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
print *, mype,"Fractional IR Energy Change Diff",abs((LW(test1,test1,lkz1+1)-OSPEC(test1,test1,lkz1+1)) &
     -(LW(test2,test2,lkz1+1)-OSPEC(test2,test2,lkz1+1))) &
     /abs(LW(test2,test2,lkz1+1)-OSPEC(test2,test2,lkz1+1))

! Fractional Change in Spectrum

print *, mype,"Maximum spectrum fractional change",maxval((abs(LW)-OSPEC)/OSPEC,OSPEC.gt.10.0**(-10.0))

! Writing

if (mype.eq.0) WRITE(enspec_handle) time
WRITESPEC = (8*pi**3)* (abs(v_1(:,:,:,0))**2 + abs(v_1(:,:,:,1))**2+abs(v_1(:,:,:,2))**2)

!CALL GATHER_WRITE(enspec_handle,WRITESPEC)

WRITESPEC = (8*pi**3)* (abs(b_1(:,:,:,0))**2 + abs(b_1(:,:,:,1))**2+abs(b_1(:,:,:,2))**2)

!CALL GATHER_WRITE(enspec_handle,WRITESPEC)

OSPEC = abs(LW)

end subroutine en_spec

subroutine bound_hels 

implicit none

LOGICAL :: magnetic

REAL :: bound
REAL :: magboundm,canboundm
REAL :: boundsm(2),bounds(2)


if (verbose) print *, "In Bound Hels"

magboundm = 0.0
canboundm = 0.0

DO j = 0,ny0_big - 1
   
   DO k = lkz1,lkz2
      DO i = 1,nkx0-1
         magboundm = magboundm + 2.0 * sum(abs(b_1(i,j,k,:))**2.0)/kmags(i,j,k)
         canboundm = canboundm + 2.0 * sum(abs(v_1(i,j,k,:))**2.0)*kmags(i,j,k)
         canboundm = canboundm + 4.0 * sqrt(sum(abs(b_1(i,j,k,:))**2.0))*sqrt(sum(abs(v_1(i,j,k,:))**2.0))
      ENDDO
   ENDDO
   !print *, mype,j,"CanBoundx",canboundm
   DO k = lkz1,lkz2
      magboundm = magboundm + sum(abs(b_1(0,j,k,:))**2.0 /kmags(0,j,k),kmags(0,j,k).gt.10**(-10.0))
      canboundm = canboundm + sum(abs(v_1(0,j,k,:))**2.0)*kmags(0,j,k)
      !if (j.eq.0) print *, sum(abs(v_1(0,j,k,:))**2.0),kmags(0,j,k)
      canboundm = canboundm + 2.0 * sqrt(sum(abs(b_1(0,j,k,:))**2.0))*sqrt(sum(abs(v_1(0,j,k,:))**2.0))
   ENDDO

   !print *, mype,j,"CanBound0 ",canboundm

   
ENDDO

if (verbose) print *, "Through Bound Loop"

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(magboundm,magbound,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(canboundm,canbound,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

magbound = 8.0 * pi**3.0 * magbound
canbound = 8.0 * pi**3.0 * canbound

canbound = canbound + magbound

if (mype.eq.0) WRITE(*,*) "Magnetic Helicity Bound",magbound
if (mype.eq.0) WRITE(*,*) "Canonical Helicity Bound",canbound

end subroutine bound_hels

subroutine mode_energy

  implicit none
  integer :: i,j,k
  real :: lwm,lcm,rwm,rcm,lw,lc,rw,rc

  lwm = 0.0
  lcm = 0.0
  rwm = 0.0
  rcm = 0.0
  
  do i = 0, nkx0-1
     do j = 0, ny0_big-1
        do k = lkz1,lkz2
           lwm = lwm + abs(sum(conjg(pcurleig(i,j,k,:))*(alpha_leftwhist(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2+1)))**2.0
           lcm = lcm + abs(sum(conjg(pcurleig(i,j,k,:))*(alpha_leftcyclo(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftcyclo(i,j,k)**2+1)))**2.0
           rwm = rwm + abs(sum(pcurleig(i,j,k,:)*(-alpha_leftwhist(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2+1)))**2.0
           rcm = rcm + abs(sum(pcurleig(i,j,k,:)*(-alpha_leftcyclo(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftcyclo(i,j,k)**2+1)))**2.0
        enddo
     enddo
  enddo

  lwm = lwm * (8.0*pi**3)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(lwm,lw,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if ((verbose)) print *, "Max LW", lw
  
  if (mype.eq.0) write(en_handle) lw

  lcm =	lcm * (8.0*pi**3)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(lcm,lc,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if ((verbose)) print *, "Max LC",lc
  
  if (mype.eq.0) write(en_handle) lc

  rwm =	rwm * (8.0*pi**3)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(rwm,rw,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if ((verbose)) print *, "Max RW",rw
  
  if (mype.eq.0) write(en_handle) rw

  rcm =	rcm * (8.0*pi**3)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(rcm,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if ((verbose)) print *, "Max RC",rc
  
  if (mype.eq.0) write(en_handle) rc
  
end subroutine mode_energy

subroutine divs

  implicit none
  complex :: divv(0:nkx0-1,0:nky0-1,0:nkz0-1),divb(0:nkx0-1,0:nky0-1,0:nkz0-1)
  integer :: i,j,k

  do i = 0,nkx0-1
     do j = 0,nky0-1
        do k = lkz1,lkz2
           divv(i,j,k) = kxgrid(i) * v_1(i,j,k,0) + kygrid(j) * v_1(i,j,k,1) + kzgrid(k) * v_1(i,j,k,2)
           divb(i,j,k) = kxgrid(i) * b_1(i,j,k,0) + kygrid(j) * b_1(i,j,k,1) +	kzgrid(k) * b_1(i,j,k,2)
        enddo
     enddo
  enddo

  if (timer.and.(mype.eq.0)) print *, "Max Abs Div v " ,maxval(abs(divv))
  if (timer.and.(mype.eq.0)) print *, "Max Abs Div b " ,maxval(abs(divb))

end subroutine divs

END MODULE diagnostics


