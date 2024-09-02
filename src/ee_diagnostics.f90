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

  PUBLIC :: initialize_diagnostics, finalize_diagnostics, diag,bv_last,bv_first,&
     !output_data, nl_test, &
     initial_wallclock,start_wallclock,&
            check_wallclock,current_wallclock!, sum3d_real

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
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: NLT3
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: AVP
  COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: WVORT
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: xi
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: MH
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CH
  REAL :: ham0,ham1
  INTEGER :: nlt3_handle 
  !NLT Testing
  !REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: NLT_bak
  INTEGER :: test_handle1,test_handle2
  INTEGER :: ierr

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
    ALLOCATE(AVP(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
    ALLOCATE(WVORT(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2))
  END IF
  IF(istep_energyspec.gt.0.and.mype==0) THEN
    CALL get_io_number
    enspec_handle=io_number
    CALL get_io_number
    mode_handle = io_number
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
    ALLOCATE(LW(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(LC(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(RW(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(RC(0:nkx0-1,0:nky0-1,lkz1:lkz2))
  END IF

  if (plot_nls) CALL initialize_debug

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
  END IF
  IF(verbose.and.(mype==0)) print *, 'Closed energy handles'
 
  if (plot_nls) CALL finalize_debug
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

     IF((istep_energy.ne.0))THEN
       IF(MOD(itime,istep_energy)==0) THEN
         IF(verbose.and.(mype.eq.0)) WRITE(*,*) "Starting energy diag.",mype
         IF (mype.eq.0) WRITE(en_handle) time
         CALL hmhdhmtn(0)
         if (verbose.and.(mype.eq.0)) write(*,*) "Found Hamiltonian",mype
         CALL mag_helicity()
         CALL err_hel1(1.0,.true.,.true.,.true.)
         CALL bound_hels(.true.)
         CALL cross_helicity()
         CALL err_hel1(1.0,.false.,.true.,.true.)
         
         CALL bound_hels(.false.)
         CALL hmhdhmtn(1)
         CALL hmhdhmtn(2)
         if (mype.eq.0) WRITE(en_handle) mhelcorr
         if (verbose.and.(mype.eq.0)) write(*,*) "Found Helicities",mype
         CALL resvischange("b")
         CALL resvischange("v")
         if (track_divs) CALL divs    
         IF(verbose.and.(mype.eq.0)) WRITE(*,*) "Done with energy diag.",mype
       END IF
     END IF

     IF((istep_energyspec.ne.0).and.(mype.eq.0)) THEN
       IF(MOD(itime,istep_energyspec)==0) THEN
         IF(verbose) WRITE(*,*) "Starting energyspec diag.",mype
         WRITE(enspec_handle) time
         CALL en_spec()
         WRITE(mode_handle) time
         CALL mode_spec()
         IF(verbose) WRITE(*,*) "Done with energyspec diag.",mype
       END IF
     END IF

  IF((istep_schpt.ne.0).and.(mype.eq.0)) THEN
    IF(MOD(itime,istep_schpt)==0) THEN
      IF(verbose) WRITE(*,*) "Writing s_checkpoint.",mype
      CALL checkpoint_out(1)
      IF(verbose) WRITE(*,*) "Done writing s_checkpoint.",mype
    END IF
  END IF

  IF((istep_gout.ne.0).and.(mype.eq.0)) THEN
    IF(MOD(itime,istep_gout)==0.or.(MOD(itime,istep_gout)==1.and.gout_2xt)) THEN

      IF(verbose) WRITE(*,*) "Starting vbout diag.",mype

      CALL checkpoint_out(4)

      IF(itime==0.and.gout_nl) THEN
        CALL checkpoint_out(5)
      ELSE IF (gout_nl) THEN
        CALL checkpoint_out(6)
      END IF

      IF(verbose) WRITE(*,*) "Done with vbout diag.",mype

    END IF
  END IF

END SUBROUTINE diag

SUBROUTINE bv_last

  IMPLICIT NONE

  ! Write the last b_1 and v_1 to a file for a warm restart
  ! This is a trial using MPI routines
  ! Code almost all from MPI 3-1 Standard

  integer :: blast_handle,vlast_handle
  integer :: sizes_4d(4),subsizes_4d(4),starts_4d(4)
  integer :: MPI_4D_COMPLEX_ARRAY
  TYPE(MPI_Status) :: status

  INTEGER :: size
  
  ! Subarray type

  print *, "In BV Last"
  sizes_4d = (/nkx0,nky0,nkz0,3/)
  subsizes_4d = (/nkx0,nky0,nkz0/n_mpi_procs,3/)
  starts_4d = (/0,0,mype*nkz0/n_mpi_procs,0/)
  
  CALL MPI_TYPE_CREATE_SUBARRAY(4,sizes_4d,subsizes_4d,starts_4d,&
       MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX,MPI_4D_COMPLEX_ARRAY,ierr)

  ! Open file and set view
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(diagdir)//'/blast_out.dat',&
       MPI_MODE_WRONLY.or.MPI_MODE_CREATE, MPI_INFO_NULL,blast_handle,ierr)
  CALL MPI_FILE_SET_VIEW(blast_handle,0, MPI_DOUBLE_COMPLEX,MPI_4D_COMPLEX_ARRAY,&
       "native",MPI_INFO_NULL,ierr)

  ! Write into file
  CALL MPI_FILE_WRITE_SHARED(blast_handle,b_1,product(subsizes_4d),MPI_DOUBLE_COMPLEX,status,ierr)
  CALL MPI_FILE_GET_SIZE(blast_handle,size,ierr)
  print	*, size	
  ! Close file
  CALL MPI_FILE_CLOSE(blast_handle, ierr)

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(diagdir)//'/vlast_out.dat',&
       MPI_MODE_WRONLY.or.MPI_MODE_CREATE, MPI_INFO_NULL,vlast_handle,ierr)
  CALL MPI_FILE_SET_VIEW(vlast_handle,0, MPI_DOUBLE_COMPLEX, MPI_4D_COMPLEX_ARRAY,&
       "native",MPI_INFO_NULL,ierr)
  CALL MPI_FILE_WRITE_SHARED(vlast_handle,v_1,product(subsizes_4d),MPI_DOUBLE_COMPLEX,status,ierr)
  CALL MPI_FILE_GET_SIZE(vlast_handle,size,ierr)
  print *, size
  CALL MPI_FILE_CLOSE(vlast_handle, ierr)

  print *, "End bv_last"

END SUBROUTINE bv_last

SUBROUTINE bv_first

  IMPLICIT NONE

  ! Read the last b_1 and v_1 to a file for a warm restart
  ! This is a trial using MPI routines
  ! Code almost all from MPI 3-1 Standard
  
  integer :: bfirst_handle,vfirst_handle
  integer :: sizes_4d(4),subsizes_4d(4),starts_4d(4)
  integer :: MPI_4D_COMPLEX_ARRAY
  TYPE(MPI_Status) :: status

  INTEGER :: size

  ! Subarray type 
  print *, "In bv_first"
  sizes_4d = (/nkx0,nky0,nkz0,3/)
  subsizes_4d = (/nkx0,nky0,nkz0/n_mpi_procs,3/)
  starts_4d = (/0,0,mype*nkz0/n_mpi_procs,0/)

  CALL MPI_TYPE_CREATE_SUBARRAY(4,sizes_4d,subsizes_4d,starts_4d,&
       MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX,MPI_4D_COMPLEX_ARRAY,ierr)

  ! Open file and set view
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(loaddir)//'/blast_out.dat',&
       MPI_MODE_RDONLY, MPI_INFO_NULL,bfirst_handle,ierr)
  CALL MPI_FILE_SET_VIEW(bfirst_handle,0, MPI_DOUBLE_COMPLEX,MPI_4D_COMPLEX_ARRAY,&
       "native",MPI_INFO_NULL,ierr)

  ! Write into file
  CALL MPI_FILE_READ_SHARED(bfirst_handle,b_1,product(subsizes_4d),MPI_DOUBLE_COMPLEX,status,ierr)

  if ((mype.eq.0)) print *, "Max b_1 first",maxval(abs(b_1))
  
  ! Close file
  CALL MPI_FILE_CLOSE(bfirst_handle, ierr)

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(loaddir)//'/vlast_out.dat',&
       MPI_MODE_RDONLY, MPI_INFO_NULL,vfirst_handle,ierr)
  CALL MPI_FILE_SET_VIEW(vfirst_handle,0, MPI_DOUBLE_COMPLEX, MPI_4D_COMPLEX_ARRAY,&
       "native",MPI_INFO_NULL,ierr)
  CALL MPI_FILE_READ_SHARED(vfirst_handle,v_1,product(subsizes_4d),MPI_DOUBLE_COMPLEX,status,ierr)
  if ((mype.eq.0)) print *, "Max v_1 first",maxval(abs(v_1))
  CALL MPI_FILE_CLOSE(vfirst_handle, ierr)

  if ((mype.eq.0)) print *, "Max b_1 first",maxval(abs(b_1))
  if ((mype.eq.0)) print *, "Max v_1 first",maxval(abs(v_1))

  print *, "End bv_first"

END SUBROUTINE bv_first
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           start_wallclock                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Some stuff for testing performance
  SUBROUTINE start_wallclock
    IMPLICIT NONE

    INTEGER :: MCLOCK

    initial_wallclock=MCLOCK()

  END SUBROUTINE start_wallclock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            check_wallclock                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE check_wallclock

   IMPLICIT NONE

   INTEGER :: current_time,max_time
   INTEGER :: MCLOCK
   INTEGER :: diff

   current_time=MCLOCK()  
   diff=current_time-initial_wallclock
   

   CALL MPI_ALLREDUCE(diff,max_time,1 &
    ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

   current_wallclock=REAL(max_time)/1000.0
   !IF(mype==0) WRITE(*,*) "initial_wallclock",initial_wallclock
   !IF(mype==0) WRITE(*,*) "current_time", current_time
   !IF(mype==0) WRITE(*,*) "diff", diff
   !IF(mype==0) WRITE(*,*) "current_wallclock", current_wallclock

  END SUBROUTINE check_wallclock

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
   if (splitx) hamsm(1) = hamsm(1) + sum(abs(v_1(1:nkx0-1,:,:,:))**2)
   if (splitx) hamsm(2) = hamsm(2) + 0.5 * sum(abs(v_1(0,:,:,:))**2)

endif
if ((opt.eq.0).or.(opt.eq.2)) then
   if (splitx) hamsm(1) = hamsm(1) + sum(abs(b_1(1:nkx0-1,:,:,:))**2)
   if (splitx) hamsm(2) = hamsm(2) + 0.5 * sum(abs(b_1(0,:,:,:))**2)
   if (mype.eq.0) hamsm(2) = hamsm(2) + real(b_1(0,0,0,2))
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(hamsm,hams,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

ham = sum(hams) * (8*(pi**3))

if (mype.eq.0) WRITE(en_handle) ham

if (opt.eq.0) then
   if (itime.eq.0) ham0 = ham

   if (verbose.and.(mype.eq.0)) print *, "Initial Hamiltonian",ham

if (mod(itime,istep_energyspec).eq.0.and.itime.gt.0) then
   ham1 = ham
   write(*,"(A20,E10.3)") "Energy Fractional Change",(ham1-ham0)/(ham0)
   ham0 = ham1
endif
endif


end subroutine hmhdhmtn

subroutine resvischange(arr_label) 

!! Returns the rate of change of the above Hamiltonian 
!! From including a resistivity term in the induction equation
!! with dimensionless resistivity eta/(B0/(e ne c)) = 1
!! Can also be used to find change due to viscosity nu/(vA^2/omega_i) = 1
!! with arr taken as v instead of b

implicit none
character(1) :: arr_label 
integer :: i,j,k
real :: k2
real :: res
real :: ind(2),all(2)

ind = 0
if (arr_label.eq."v") then
   ind(1) = vnu*sum(spread(kmags(1:nkx0-1,:,:)**(2.0*hyp),4,3)*abs(v_1(1:nkx0-1,:,:,:)**2.0))
   ind(2) = vnu*0.5 * sum(spread(kmags(0,:,:)**(2.0*hyp),3,3)*abs(v_1(0,:,:,:)**2.0))
endif
if (arr_label.eq."b") then
   ind(1) = eta*sum(spread(kmags(1:nkx0-1,:,:)**(2.0*hyp),4,3)*abs(b_1(1:nkx0-1,:,:,:)**2.0))
   ind(2) = eta*0.5* sum(spread(kmags(0,:,:)**(2.0*hyp),3,3)*abs(b_1(0,:,:,:)**2.0))
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(ind,all,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

res = sum(all)*(16 * pi**3)

if (mype.eq.0) WRITE(en_handle) res

end subroutine resvischange

subroutine initialize_debug

    CALL get_io_number
    dbio = io_number
    CALL get_io_number
    dvio = io_number 
    CALL get_io_number
    bdvio = io_number
    CALL get_io_number
    vdbio = io_number
    CALL get_io_number
    bdcbio = io_number
    CALL get_io_number
    cbdbio = io_number
    CALL get_io_number
    vdvio = io_number
    CALL get_io_number
    bdbio = io_number
    CALL get_io_number
    db2io = io_number


    IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/dissb_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/dissv_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/bdv_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/vdb_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/bdcb_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/cbdb_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/vdv_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/bdb_out.dat',exist=file_exists)
      INQUIRE(file=trim(diagdir)//'/db2_out.dat',exist=file_exists)

      IF(file_exists) THEN
        OPEN(unit=dbio,file=trim(diagdir)//'/dissb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=dvio,file=trim(diagdir)//'/dissv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdvio,file=trim(diagdir)//'/bdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdbio,file=trim(diagdir)//'/vdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdcbio,file=trim(diagdir)//'/bdcb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=cbdbio,file=trim(diagdir)//'/cbdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdvio,file=trim(diagdir)//'/vdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdbio,file=trim(diagdir)//'/bdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=db2io,file=trim(diagdir)//'/db2_out.dat',form='unformatted', status='REPLACE',access='stream')
      ELSE
        OPEN(unit=dbio,file=trim(diagdir)//'/dissb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=dvio,file=trim(diagdir)//'/dissv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdvio,file=trim(diagdir)//'/bdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdbio,file=trim(diagdir)//'/vdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdcbio,file=trim(diagdir)//'/bdcb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=cbdbio,file=trim(diagdir)//'/cbdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdvio,file=trim(diagdir)//'/vdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdbio,file=trim(diagdir)//'/bdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=db2io,file=trim(diagdir)//'/db2_out.dat',form='unformatted', status='REPLACE',access='stream')
      END IF
    ELSE
        OPEN(unit=dbio,file=trim(diagdir)//'/dissb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=dvio,file=trim(diagdir)//'/dissv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdvio,file=trim(diagdir)//'/bdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdbio,file=trim(diagdir)//'/vdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdcbio,file=trim(diagdir)//'/bdcb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=cbdbio,file=trim(diagdir)//'/cbdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=vdvio,file=trim(diagdir)//'/vdv_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=bdbio,file=trim(diagdir)//'/bdb_out.dat',form='unformatted', status='REPLACE',access='stream')
        OPEN(unit=db2io,file=trim(diagdir)//'/db2_out.dat',form='unformatted', status='REPLACE',access='stream')
    END IF

if (verbose) print *, 'Debug files opened'

end subroutine initialize_debug

subroutine finalize_debug

CLOSE(dbio)
CLOSE(dvio)
CLOSE(bdvio)
CLOSE(vdbio)
CLOSE(bdcbio)
CLOSE(cbdbio)
CLOSE(vdvio)
CLOSE(bdbio)
CLOSE(db2io)

if (verbose.and.(mype.eq.0)) print *, 'Debug files closed' 

end subroutine finalize_debug

subroutine vec_potential

implicit none
integer :: i,j,k,ind

! Ak = (i k x Bk)/k^2
AVP = cmplx(0.0,0.0)

do i = 0,nkx0-1
  do j = 0,nky0-1
    do k = lkz1,lkz2
      do ind = 0,2
        if (max(i,j,k).gt.0) then
          if (ind.eq.0) AVP(i,j,k,ind) = cmplx(0.0,1.0) * (kygrid(j)*b_1(i,j,k,2)-kzgrid(k)*b_1(i,j,k,1))/(kmags(i,j,k)**2)
          if (ind.eq.1) AVP(i,j,k,ind) = cmplx(0.0,1.0) * (kzgrid(k)*b_1(i,j,k,0)-kxgrid(i)*b_1(i,j,k,2))/(kmags(i,j,k)**2)
          if (ind.eq.2) AVP(i,j,k,ind) = cmplx(0.0,1.0) * (kxgrid(i)*b_1(i,j,k,1)-kygrid(j)*b_1(i,j,k,0))/(kmags(i,j,k)**2)
        endif
      enddo
    enddo
   enddo
enddo

end subroutine vec_potential

SUBROUTINE mag_helicity 

implicit none
real :: maghel
integer :: i,j,k,ind
real :: mhs(2),mhsm(2)

mhsm = 0.0

CALL vec_potential()

if (splitx) then
   mhsm(1) = mhsm(1) + 2.0*sum(real(AVP(1:nkx0-1,:,:,:)*conjg(b_1(1:nkx0-1,:,:,:))))
   mhsm(2) = mhsm(2) + sum(AVP(0,:,:,:)*conjg(b_1(0,:,:,:)))
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(mhsm,mhs,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

maghel = sum(mhs) * (2.0 * pi)**3.0

if (mype.eq.0) WRITE(en_handle) maghel

END SUBROUTINE mag_helicity

subroutine vorticity

implicit none
integer :: ind,i,j,k

WVORT = cmplx(0.0,0.0)

do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
          if (ind.eq.0) WVORT(i,j,k,ind) = kygrid(j) * v_1(i,j,k,2) - kzgrid(k) * v_1(i,j,k,1)
          if (ind.eq.1) WVORT(i,j,k,ind) = kzgrid(k) * v_1(i,j,k,0) - kxgrid(i) * v_1(i,j,k,2)
          if (ind.eq.2) WVORT(i,j,k,ind) = kxgrid(i) * v_1(i,j,k,1) - kygrid(j) * v_1(i,j,k,0)
      enddo
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
  
  if (splitx) chsm(1) = 2.0 * sum(real((AVP(1:nkx0-1,:,:,:)+v_1(1:nkx0-1,:,:,:))*conjg(b_1(1:nkx0-1,:,:,:)+WVORT(1:nkx0-1,:,:,:))))
  if (splitx) chsm(2) = real(sum((AVP(0,:,:,:)+v_1(0,:,:,:))*conjg(b_1(0,:,:,:)+WVORT(0,:,:,:))))
  chsm(2) = chsm(2) + v_1(0,0,0,2)
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(chsm,chs,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
   
  crosshel = sum(chs) * (2.0*pi)**3.0
   
  if (mype.eq.0) WRITE(en_handle) crosshel

end subroutine cross_helicity

subroutine en_spec() 

implicit none

! This doesn't need to be adjusted for 2x only non ksplit = 0 modes

WRITE(enspec_handle) (8*pi**3)* (abs(v_1(:,:,:,0))**2 + abs(v_1(:,:,:,1))**2+abs(v_1(:,:,:,2))**2)
WRITE(enspec_handle) (8*pi**3)* (abs(b_1(:,:,:,0))**2+ abs(b_1(:,:,:,1))**2+abs(b_1(:,:,:,2))**2)

end subroutine en_spec

subroutine err_hel1(sf,m,b,num)

implicit none

REAL :: sf
LOGICAL :: m
LOGICAL :: b
LOGICAL :: num
REAL :: mherr
REAL :: Ck,k2xderiv,kxderiv,k2yderiv,kyderiv,k2zderiv,kzderiv
REAL :: xcont,ycont,zcont
INTEGER :: i,j,k,ind

mherr = 0.0

if (mype.eq.0) WRITE(en_handle) 0.0

end subroutine err_hel1

subroutine bound_hels(m) 

implicit none

LOGICAL :: m
REAL :: bound
REAL :: boundsm(2),bounds(2)

if (splitx) then
   boundsm(1) = 2.0*sum(sum(b_1(1:nkx0-1,:,:,:)*conjg(b_1(1:nkx0-1,:,:,:)),4)/kmags(1:nkx0-1,:,:),mask=(kmags.gt.0))
   boundsm(2) = sum(sum(abs(b_1(0,:,:,:))**2.0,dim=3)/kmags(0,:,:),mask=(kmags(0,:,:).gt.0))
if (m.eq..false.) then 
   boundsm(1) = boundsm(1) + 2.0*sum(sum(v_1(1:nkx0-1,:,:,:)*conjg(v_1(1:nkx0-1,:,:,:)),4)*kmags(1:nkx0-1,:,:),mask=(kmags.gt.0))
   boundsm(2) = boundsm(2) + sum(sum(v_1(0,:,:,:)*conjg(v_1(0,:,:,:)),dim=3)*kmags(0,:,:),mask=(kmags(0,:,:).gt.0))
   boundsm(1) = boundsm(1) + 4.0*sum(sqrt(sum(abs(v_1(1:nkx0-1,:,:,:))**2,4))*sqrt(sum(abs(b_1(1:nkx0-1,:,:,:))**2,4)),mask=(kmags.gt.0))
   boundsm(2) = boundsm(2) + 2.0*sum(sqrt(sum(abs(v_1(0,:,:,:))**2,3))*sqrt(sum(abs(b_1(0,:,:,:))**2,3)),mask=(kmags(0,:,:).gt.0))
endif
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(boundsm,bounds,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

bound = 8.0 * pi**3.0 * sum(bounds)

if (mype.eq.0) WRITE(en_handle) bound

end subroutine bound_hels

subroutine mode_spec

  implicit none
  integer :: i,j,k
  
  do i = 0, nkx0-1
     do j = 0, nky0-1
        do k = lkz1,lkz2
           LW(i,j,k) = sum(conjg(pcurleig(i,j,k,:))*(alpha_leftwhist(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2+1))
           LC(i,j,k) = sum(conjg(pcurleig(i,j,k,:))*(alpha_leftcyclo(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftcyclo(i,j,k)**2+1))
           RW(i,j,k) = sum(pcurleig(i,j,k,:)*(-alpha_leftwhist(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2+1))
           RC(i,j,k) = sum(pcurleig(i,j,k,:)*(-alpha_leftcyclo(i,j,k)*b_1(i,j,k,:)+v_1(i,j,k,:))/sqrt(alpha_leftcyclo(i,j,k)**2+1))
        enddo
     enddo
  enddo

  if ((verbose).and.itime.lt.100) then
     print *, "Max LW ", sum(0.5*abs(LW)**2)*(16.0*pi**3)
     print *, "Max LC ", sum(0.5*abs(LC)**2)*(16.0*pi**3)
     print *, "Max RW ", sum(0.5*abs(RW)**2)*(16.0*pi**3)
     print *, "Max RC ", sum(0.5*abs(RC)**2)*(16.0*pi**3)
  endif
     
  WRITE(mode_handle) 0.5 * abs(LW)**2 * (16.0*pi**3)
  WRITE(mode_handle) 0.5 * abs(LC)**2 * (16.0*pi**3)
  WRITE(mode_handle) 0.5 * abs(RW)**2 * (16.0*pi**3)
  WRITE(mode_handle) 0.5 * abs(RC)**2 * (16.0*pi**3)

end subroutine mode_spec

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

  if (verbose.and.(mype.eq.0)) print *, "Max Abs Div v " ,maxval(abs(divv))
  if (verbose.and.(mype.eq.0)) print *, "Max Abs Div b " ,maxval(abs(divb))

end subroutine divs

END MODULE diagnostics


