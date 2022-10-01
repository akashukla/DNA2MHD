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

  PUBLIC :: initialize_diagnostics, finalize_diagnostics, diag,&
     !output_data, nl_test, &
     initial_wallclock,start_wallclock,&
            check_wallclock,current_wallclock!, sum3d_real

  PRIVATE

  INTEGER :: ffm_handle,omega_handle,glast_handle,ev_handle,en_handle,enspec_handle,&
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
  INTEGER :: nlt3_handle 
  !NLT Testing
  !REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: NLT_bak
  INTEGER :: test_handle1,test_handle2

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

  !IF(checkpoint_read) THEN
  !  INQUIRE(file=trim(diagdir)//'/ffm.dat',exist=file_exists)
  !  IF(file_exists) THEN
  !     IF(mype==0) OPEN(unit=ffm_handle,file=trim(diagdir)//'/ffm.dat',status='unknown',position='append')
  !  ELSE
  !     IF(mype==0) OPEN(unit=ffm_handle,file=trim(diagdir)//'/ffm.dat',status='unknown')
  !  END IF
  !ELSE
  !  IF(mype==0) OPEN(unit=ffm_handle,file=trim(diagdir)//'/ffm.dat',status='unknown')
  !END IF

  !IF(mype==0) WRITE(ffm_handle,*) "#time,phi_tot,flux_tot"
  !IF(istep_gamma.gt.0) ALLOCATE(phi_last(0:nkx0-1,0:nky0-1,0:nkz0-1))
  !!IF((istep_energy3d.gt.0.or.istep_fmom3d.gt.0).and.mype==0)  CALL initialize_netcdf
  !IF(istep_energy3d.gt.0.and.mype==0) THEN
  !  CALL initialize_energy3d
  !END IF
  !IF(istep_fmom3d.gt.0.and.mype==0) THEN
  !  CALL initialize_fmom3d
  !END IF
  IF(istep_energy.gt.0.and.mype==0) THEN
    CALL get_io_number
    en_handle=io_number
    !OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',status='unknown')
    IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/energy_out.dat',exist=file_exists)
      IF(file_exists) THEN
        OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',form='unformatted', status='REPLACE',access='stream')
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
    energy_last=0.0
    time_last=time
  END IF
  IF(istep_energyspec.gt.0.and.mype==0) THEN
    CALL get_io_number
    enspec_handle=io_number
    !OPEN(unit=en_handle,file=trim(diagdir)//'/energy_out.dat',status='unknown')                                                                    
    IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/energyspec_out.dat',exist=file_exists)
      IF(file_exists) THEN
        OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='REPLACE',access='stream')
      ELSE
        OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='REPLACE',access='stream')
      END IF
    ELSE
      OPEN(unit=enspec_handle,file=trim(diagdir)//'/energyspec_out.dat',form='unformatted', status='REPLACE',access='stream')
    END IF
  END IF


  !IF(mype==0) WRITE(*,*) "Done initializing netcdf."
  ! IF(istep_hermite.gt.0.and.mype==0) THEN
  !   CALL get_io_number
  !   herm_handle=io_number
  !   !OPEN(unit=en_handle,file=trim(diagdir)//'/energy.dat',status='unknown')
  !   IF(checkpoint_read) THEN
  !     INQUIRE(file=trim(diagdir)//'/energy_hermite.dat',exist=file_exists)
  !     IF(file_exists) THEN
  !       OPEN(unit=herm_handle,file=trim(diagdir)//'/energy_hermite.dat',status='unknown',position='append')
  !     ELSE
  !       OPEN(unit=herm_handle,file=trim(diagdir)//'/energy_hermite.dat',status='unknown')
  !     END IF
  !   ELSE
  !     OPEN(unit=herm_handle,file=trim(diagdir)//'/energy_hermite.dat',status='unknown')
  !   END IF
  !   WRITE(char_nv0,'(i4.4)') nv0+1
  !   !IF(mype==0) WRITE(*,*) "char_nv0",char_nv0
  ! END IF

  !IF(istep_gout.gt.0) THEN
  !  IF(checkpoint_read) THEN
  !    !CALL checkpoint_out(4)
  !  ELSE
  !    CALL checkpoint_out(3)
  !  END IF
  !END IF

  ! IF(istep_nlt.gt.0) THEN
  !   CALL initialize_diag_nlt
  ! END IF

  ! IF(istep_nlt_triple.gt.0) THEN
  !   CALL initialize_diag_nlt_triple
  ! END IF

  ! IF(istep_eshells.gt.0) THEN
  !   CALL initialize_diag_eshells
  ! END IF

  ! IF(istep_gk.gt.0) THEN
  !   CALL initialize_diag_gk
  ! END IF

  ! IF(istep_gknl.gt.0) THEN
  !   CALL initialize_diag_gknl
  ! END IF

  !Initialize GyroLES
  !GyroLES should be the last diagnostic to execute
  !IF (GyroLES.or.Gyroz) call initialize_GyroLES
  !IF (Corr) call initialize_corr

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
  ! IF(allocated(phi_last)) DEALLOCATE(phi_last)
  ! !IF(mype==0.and.verbose) WRITE(*,*) "finalizing netcdf."
  ! !IF((istep_energy3d.gt.0.or.istep_fmom3d.gt.0).and.mype==0)  CALL finalize_netcdf
  ! IF(istep_energy3d.gt.0.and.mype==0) THEN
  !    IF(mype==0.and.verbose) WRITE(*,*) "finalizing energy3d."
  !    CALL finalize_energy3d
  ! END IF
  ! IF(istep_fmom3d.gt.0.and.mype==0) THEN
  !    IF(mype==0.and.verbose) WRITE(*,*) "finalizing fmom3d."
  !    CALL finalize_fmom3d
  ! END IF
  ! !IF(verbose) WRITE(*,*) "Writing checkpoint.",mype
  ! !WRITE(*,*) "Writing checkpoint.",mype

  IF(istep_energy.gt.0.and.mype==0) THEN
    CLOSE(en_handle)
  END IF
  IF(istep_energyspec.gt.0.and.mype==0) THEN
    CLOSE(enspec_handle)
  END IF
  IF(verbose) print *, 'Closed energy handles'
  ! IF(istep_hermite.gt.0.and.mype==0) CLOSE(herm_handle)

  ! IF(istep_gk.gt.0) THEN
  !   CALL finalize_diag_gk
  ! END IF

  ! IF(istep_nlt.gt.0) THEN
  !   CALL finalize_diag_nlt
  ! END IF

  ! IF(istep_nlt_triple.gt.0) THEN
  !   CALL finalize_diag_nlt_triple
  ! END IF

  ! IF(istep_gknl.gt.0) THEN
  !   CALL finalize_diag_gknl
  ! END IF

  ! !GyroLES should be the last diagnostic to execute
  ! IF (GyroLES.or.Gyroz) call finalize_GyroLES
  ! IF (Corr) call finalize_corr
 
  if (plot_nls) CALL finalize_debug
  IF(verbose.and.plot_nls) print *, 'Closed nl handles'

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

!!   !This calculates average fields, fluxes and moments (ffm) (acutally just phi
!!   !and Q at this time
!!   IF(istep_ffm.ne.0) THEN
!!     IF(MOD(itime,istep_ffm)==0) THEN
!!       IF(verbose) WRITE(*,*) "Starting ffm diag.",mype
!!       CALL get_phi(g_1)
!!       !phi_kz(:)=sum(sum(conjg(phi(:,:,:))*phi(:,:,:),1),1)
!!       IF(nkz0.gt.1) THEN
!!         !phi_tot=REAL(2*sum(conjg(phi_kz(1:nkz0-1))*phi_kz(1:nkz0-1)+conjg(phi_kz(0))*phi_kz(0)))
!!         CALL sum3d_real(REAL(conjg(phi)*phi),phi_tot)
!!       ELSE
!!         phi_tot=REAL(sum(sum(sum(   (conjg(phi(:,:,:))*phi(:,:,:))   ,1),1),1))
!!       END IF
!!       IF(phi_tot.gt.1.0e30) THEN
!!         IF(mype==0) WRITE(ffm_handle,'(3es16.8)') "#!!!!!Overflow!!!!"
!!         IF(mype==0) WRITE(ffm_handle,'(3es16.8)') "#!!!!!Overflow!!!!"
!!         IF(mype==0) WRITE(ffm_handle,'(3es16.8)') "#!!!!!Overflow!!!!"
!!         STOP "Overflow!"
!!       END IF
!!       CALL get_flux_total(flux_tot)
!!       IF(mype==0) WRITE(ffm_handle,'(3es16.8)') time,phi_tot,flux_tot
!!       IF(mype==0) flush(ffm_handle)
!!       IF(verbose) WRITE(*,*) "Done with ffm diag.",mype
!!     END IF
!!   END IF
!! 
!!   !This calculates 3D fields and moments (phi and pressure at this time)
!!   IF(istep_fmom3d.ne.0) THEN
!!     IF(MOD(itime,istep_fmom3d)==0) THEN
!!       IF(verbose) WRITE(*,*) "Starting fmom3d diag.",mype
!!       CALL get_phi(g_1)
!!       IF(mype==0) THEN 
!!         WRITE(fmom3d_handle) time
!!         WRITE(fmom3d_handle) phi
!!       END IF
!!       ALLOCATE(mom3d(0:nkx0-1,0:nky0-1,lkz1:lkz2))
!!       CALL get_press_3d(mom3d)
!!       IF(mype==0) WRITE(fmom3d_handle) mom3d 
!!       IF(mype==0) flush(fmom3d_handle)
!!       DEALLOCATE(mom3d)
!!       IF(verbose) WRITE(*,*) "Done with fmom3d diag.",mype
!!     END IF
!!   END IF
!! 
!!   !This calculates 3D energy quantities
!!   IF(istep_energy3d.ne.0) THEN
!!     IF(MOD(itime,istep_energy3d)==0) THEN
!!       IF(verbose) WRITE(*,*) "Starting energy3d diag.",mype
!!       ALLOCATE(energy3d(0:nkx0-1,0:nky0-1,lkz1:lkz2))
!!       ALLOCATE(energy3d_temp(0:nkx0-1,0:nky0-1,lkz1:lkz2))
!! 
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       CALL check_real(g_1,real_err)
!!       IF(real_err.gt.1.0e-10) THEN
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "Error! reality constraint violated. mype,real_err",mype,real_err
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       END IF
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!! 
!!       IF(mype==0) WRITE(energy3d_handle) time
!! 
!!       !Total energy
!!       CALL get_energy3d(g_1,energy3d,1) 
!!       CALL get_energy3d(g_1,energy3d_temp,2) 
!!       energy3d=energy3d+energy3d_temp
!!       IF(mype==0) WRITE(energy3d_handle) energy3d
!! 
!!       !Drive term
!!       CALL get_energy3d(g_1,energy3d,4) 
!!       IF(mype==0) WRITE(energy3d_handle) energy3d
!! 
!!       !Collision term
!!       CALL get_energy3d(g_1,energy3d,5) 
!!       CALL get_energy3d(g_1,energy3d_temp,6) 
!!       energy3d=energy3d+energy3d_temp
!!       IF(mype==0) WRITE(energy3d_handle) energy3d
!! 
!!       !Hyp_x,y,z,conv
!!       CALL get_energy3d(g_1,energy3d,7) 
!!       CALL get_energy3d(g_1,energy3d_temp,9) 
!!       energy3d=energy3d+energy3d_temp
!!       IF(mype==0) WRITE(energy3d_handle) energy3d
!! 
!!       IF(mype==0) flush(energy3d_handle)
!!       DEALLOCATE(energy3d_temp)
!!       DEALLOCATE(energy3d)
!! 
!!       IF(verbose) WRITE(*,*) "Done with energy3d diag.",mype
!!     END IF
!!   END IF
!! 
     IF(istep_energy.ne.0) THEN
       IF(MOD(itime,istep_energy)==0) THEN
         IF(verbose) WRITE(*,*) "Starting energy diag.",mype
         WRITE(en_handle) time
         WRITE(en_handle) hmhdhmtn(v_1,b_1)
         WRITE(en_handle) mag_helicity(b_1)
         WRITE(en_handle) cross_helicity(v_1,b_1)
         WRITE(en_handle) eta*resvischange(b_1)
         WRITE(en_handle) vnu*resvischange(v_1)
 
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!! !      CALL get_real(g_1)
!!       !CALL check_real(g_1,real_err)
!!       !IF(real_err.gt.1.0e-10) THEN
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !  WRITE(*,*) "Error! reality constraint violated. mype,real_err",mype,real_err
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !  WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       !END IF
!! 
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       !!!!!!!!!!!Temporary!!!!!!!!!!!
!!       CALL get_energy_terms
!! 
         IF(verbose) WRITE(*,*) "Done with energy diag.",mype
       END IF
     END IF

     IF(istep_energyspec.ne.0) THEN
       IF(MOD(itime,istep_energyspec)==0) THEN
         IF(verbose) WRITE(*,*) "Starting energyspec diag.",mype
         WRITE(enspec_handle) time
         WRITE(enspec_handle) en_spec(v_1,b_1)
         IF(verbose) WRITE(*,*) "Done with energyspec diag.",mype
       END IF
     END IF

  IF(istep_schpt.ne.0) THEN
    IF(MOD(itime,istep_schpt)==0) THEN
      IF(verbose) WRITE(*,*) "Writing s_checkpoint.",mype
      CALL checkpoint_out(1)
      IF(verbose) WRITE(*,*) "Done writing s_checkpoint.",mype
    END IF
  END IF

!!   !Outputs hermite spectrum summed over all k's
!!   IF(istep_hermite.ne.0) THEN
!!     IF(MOD(itime,istep_hermite)==0) THEN
!! 
!!       IF(verbose) WRITE(*,*) "Starting hermite diag.",mype
!!       CALL get_energy_hermite
!!       IF(verbose) WRITE(*,*) "Done with hermite diag.",mype
!! 
!!     END IF
!!   END IF

  !Outputs the entire distribution function
  IF(istep_gout.ne.0) THEN
    IF(MOD(itime,istep_gout)==0.or.(MOD(itime,istep_gout)==1.and.gout_2xt)) THEN

      IF(verbose) WRITE(*,*) "Starting gout diag.",mype

      IF(itime==0) THEN
        CALL checkpoint_out(3)
      ELSE
        CALL checkpoint_out(4)
      END IF

      IF(itime==0.and.gout_nl) THEN
        CALL checkpoint_out(5)
      ELSE IF (gout_nl) THEN
        CALL checkpoint_out(6)
      END IF

      IF(verbose) WRITE(*,*) "Done with gout diag.",mype

    END IF
  END IF

!!Maybe include constraints later
!!   !Enforces reality constraint for kx=0 modes
!!   IF(istep_real.ne.0) THEN
!!     IF(MOD(itime,istep_real)==0) THEN
!! 
!!       IF(verbose) WRITE(*,*) "Starting get_real.",mype
!!       CALL get_real(g_1)
!!       IF(verbose) WRITE(*,*) "Done with get_real.",mype
!! 
!!     END IF
!!   END IF
!! 
!!   !Tests nonlinear energy conservation
!!   IF(istep_nltest.ne.0) THEN
!!     IF(MOD(itime,istep_nltest)==0.and.itime.ne.0) THEN
!!       !CALL nl_test 
!!     END IF
!!   END IF
!! 
!!   !Calculates growth rate for linear runs
!!   IF(istep_gamma.ne.0) THEN
!!    IF(MOD(itime,istep_gamma)==0) THEN
!!     IF(itime==0) THEN
!!       last_frequency=cmplx(0.0,0.0)
!!       phi_last=phi
!!     ELSE
!!       IF(np_kz.gt.1) STOP "Must implement istep_gamma diag for kz parallelization!"
!!       gam=sum(sum(sum(  &
!!             (phi-phi_last)/(istep_gamma*dt*phi)  &
!!             ,1),1),1)/REAL(nkx0*nky0*nkz0)
!!       IF(abs((gam-last_frequency)/gam).lt.eps_converge) THEN
!!         IF(mype==0) THEN
!!           WRITE(*,*) "Simulation converged.  gamma, omega", REAL(gam),aimag(gam)
!!           OPEN(unit=omega_handle,file=trim(diagdir)//'/omega.dat',status='unknown')
!!           WRITE(omega_handle,*) "#gamma,omega"
!!           WRITE(omega_handle,*) REAL(gam),aimag(gam)
!!           CLOSE(omega_handle)
!!         END IF
!!         continue_run=.false.
!!       END IF
!!       phi_last=phi
!!       last_frequency=gam
!!     END IF
!!    END IF
!!   END IF
!! 
!!   !Nonlinear transfer diagnostics
!!   IF(istep_nlt.ne.0) THEN
!!     IF(np_kz.gt.1) STOP "Must implement istep_nlt diagnostics for kz parallelization."
!!     IF(MOD(itime,istep_nlt)==0.and.itime.ne.0) THEN
!!       !CALL nlt_shell_test1 
!!       ALLOCATE(NLT(num_shells,num_shells,0:ikz_minmax,0:ikz_minmax))
!!       !n-resolved NLT for selected n:
!!       IF(output_nlt_n) ALLOCATE(NLT_n(num_shells,num_shells,0:ikz_minmax,0:ikz_minmax,11))
!!       !ALLOCATE(nltk(-ikx_minmax:ikx_minmax,-iky_minmax:iky_minmax,-ikz_minmax:ikz_minmax))
!!       !CALL nlt_test 
!!       IF(mype==0.and.verbose) WRITE(*,*) "!!!!!!!!!Calling get_nlt_shells."
!!       CALL get_nlt_shells
!!       DEALLOCATE(NLT)
!!       IF(output_nlt_n) DEALLOCATE(NLT_n)
!!       !DEALLOCATE(nltk)
!!     END IF
!!   END IF
!! 
!!   !Triple shell nonlinear transfer diagnostics
!!   IF(istep_nlt_triple.ne.0) THEN
!!     IF(np_kz.gt.1) STOP "Must implement istep_nlt diagnostics for kz parallelization."
!!     IF(MOD(itime,istep_nlt_triple)==0.and.itime.ne.0) THEN
!!       ALLOCATE(NLT3(num_shells,num_shells,num_shells))
!!       IF(mype==0.and.verbose) WRITE(*,*) "!!!!!!!!!Calling get_nlt_triple."
!!       CALL get_nlt_triple
!!       DEALLOCATE(NLT3)
!!     END IF
!!   END IF
!! 
!! 
!!   !Energy shell diagnostics
!!   IF(istep_eshells.ne.0) THEN
!!     IF(np_kz.gt.1) STOP "Must implement istep_eshells diagnostic for kz parallelization."
!!     IF(MOD(itime,istep_eshells)==0.and.itime.ne.0) THEN
!!         IF(verbose.and.mype==0) write(*,*) "Calling diag_eshells"
!!         CALL diag_eshells
!!     END IF
!!   END IF
!! 
!!   !Higher time-resolution output of distribution function for specific k's
!!   IF(istep_gk.ne.0) THEN
!!     IF(np_kz.gt.1) STOP "Must implement istep_gk diagnostic for kz parallelization."
!!     IF(MOD(itime,istep_gk)==0) THEN
!!       IF(np_kz.gt.1) STOP "Must implement istep_gk diagnostic for kz parallelization."
!! 
!!       IF(verbose) WRITE(*,*) "Starting gk diag.",mype
!!       CALL diag_gk
!!       IF(verbose) WRITE(*,*) "Done with gk diag.",mype
!! 
!!     END IF
!!   END IF
!! 
!!   !Higher time-resolution output of nonlinearity for specific k's
!!   IF(istep_gknl.ne.0) THEN
!!     IF(np_kz.gt.1) STOP "Must implement istep_gknl diagnostic for kz parallelization."
!!     IF(MOD(itime,istep_gknl)==0) THEN
!!       IF(np_kz.gt.1) STOP "Must implement istep_gknl diagnostic for kz parallelization."
!! 
!!       IF(verbose) WRITE(*,*) "Starting gknl diag.",mype
!!       CALL diag_gknl
!!       IF(verbose) WRITE(*,*) "Done with gknl diag.",mype
!! 
!!     END IF
!!   END IF
!! 
!!   !For gyroLES routine  
!!   !GyroLES should be the last diagnostic to execute
!!   IF (GyroLES) THEN 
!!     IF (MODULO(itime,istep_GyroLES) .EQ. 0) THEN 
!!        CALL get_phi(g_1)
!!        IF (Gyroherm) THEN
!!             CALL exec_GyroLES_herm(g_1)
!!        ELSE 
!!             CALL exec_GyroLES 
!!        ENDIF 
!!     ENDIF
!!   ENDIF
!! 
!!   IF (Gyroz) THEN 
!!     IF (MODULO(itime,istep_GyroLES) .EQ. 0) THEN 
!!        CALL get_phi(g_1)
!!        CALL exec_GyroLES_z 
!!     ENDIF
!!   ENDIF
!! 
!!   !Correlation diagnostics
!!   IF (Corr) THEN 
!!     IF (MODULO(itime,istep_GyroLES) .EQ. 0) THEN 
!!        CALL get_phi(g_1)
!!        CALL Correlation 
!!     ENDIF
!!   ENDIF
!! 
!! 
!!   IF(verbose) WRITE(*,*) "END of diag,mype",mype

END SUBROUTINE diag

!SUBROUTINE diagnostics_last
!
!  IMPLICIT NONE
!  INTEGER :: i,ierr
!  REAL :: norm
!  COMPLEX :: g_out(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1)
!  COMPLEX :: g_out_temp(0:nkx0-1,0:nky0-1,0:nkz0-1,0:nv0-1)
!
!  g_out_temp=cmplx(0.0,0.0)
!  g_out_temp(:,:,:,lv1:lv2)=g_1(:,:,:,lv1:lv2)
!  CALL MPI_ALLREDUCE(g_out_temp,g_out,nkx0*nky0*nkz0*nv0&
!      ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!  IF(mype==0) THEN
!
!     OPEN(unit=glast_handle,file=trim(diagdir)//'/g_last.dat',status='unknown')
!  
!     norm=REAL(sqrt(sum(conjg(g_out(0,0,0,:))*g_out(0,0,0,:))))
!     DO i=0,nv0-1
!       WRITE(glast_handle,*) REAL(i),REAL(g_out(0,0,0,i))/norm,aimag(g_out(0,0,0,i))/norm
!     END DO
!
!     CLOSE(glast_handle)
!
!  END IF
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!END SUBROUTINE diagnostics_last


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                output_data                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Output for eigenvalue solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE output_data
!! 
!!   INTEGER :: i
!!   COMPLEX :: vector_temp(nv0)
!!   INTEGER :: j
!!   REAL :: eig_sort(nv0,2)
!!   INTEGER :: sort_ind(nv0)
!! 
!!   IF(mype==0) OPEN(unit=ev_handle,file=trim(diagdir)//'/eigenvalues.dat',status='unknown')
!! 
!!   CALL sort_evs_diag(eig_sort,sort_ind)
!! 
!!   DO i=1,nv0
!!     IF(mype==0) WRITE(ev_handle,'(2es16.6)') aimag(eig(sort_ind(i))),REAL(eig(sort_ind(i)))
!!   END DO
!!   
!!   IF(mype==0) CLOSE(ev_handle)
!! 
!!    IF(right_vec) THEN
!!    IF(mype==0) OPEN(unit=ev_handle,file=trim(diagdir)//'/right_evec.dat',status='unknown')
!!    DO j=1,nv0
!!     !WRITE vector to file
!!     vector_temp=loc_rvec(:,sort_ind(j))
!!     IF(mype==0) THEN
!!       DO i=1,nv0
!!         !WRITE(*,*) i-0
!!         IF(mype==0) WRITE(ev_handle,*) i-1.0,REAL(vector_temp(i)),aimag(vector_temp(i)) 
!!       END DO
!!       IF(mype==0) WRITE(ev_handle,*) ""
!!       IF(mype==0) WRITE(ev_handle,*) ""
!!     END IF
!!    END DO !j loop
!!    IF(mype==0) CLOSE(ev_handle)
!!     
!!    END IF !right_vec
!! 
!!    IF(left_vec) THEN
!!    IF(mype==0) OPEN(unit=ev_handle,file=trim(diagdir)//'/left_evec.dat',status='unknown')
!!    DO j=1,nv0
!!     vector_temp=loc_lvec(:,sort_ind(j))
!!     !WRITE vector to file
!!     IF(mype==0) THEN
!!       DO i=1,nv0
!!         !WRITE(*,*) i-0
!!         WRITE(ev_handle,*) i-1.0,REAL(vector_temp(i)),aimag(vector_temp(i)) 
!!       END DO
!!       WRITE(ev_handle,*) ""
!!       WRITE(ev_handle,*) ""
!!     END IF
!!    END DO !j loop
!!    IF(mype==0) CLOSE(ev_handle)
!!     
!!    END IF !left_vec
!! 
!! END SUBROUTINE output_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            sort_evs_diag                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE sort_evs_diag(eig_sort,sort_ind)
!!   
!!   REAL, INTENT(out) :: eig_sort(nv0,2)
!!   INTEGER, INTENT(out) :: sort_ind(nv0)
!!   REAL :: grate(nv0),freq(nv0)
!!   INTEGER :: i,mloc(1)
!!   REAL :: mingrate
!! 
!!   grate=REAL(eig)
!!   freq=aimag(eig)
!!   mingrate=minval(grate)
!!   
!!   DO i=1,nv0
!! 
!!     mloc=maxloc(grate)
!!     eig_sort(i,1)=grate(mloc(1))
!!     eig_sort(i,2)=freq(mloc(1))
!!     sort_ind(i)=mloc(1) 
!!     grate(mloc(1))=mingrate-100.0
!! 
!!   END DO
!! 
!! END SUBROUTINE sort_evs_diag


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               nl_test                                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Tests energy conservation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE nl_test
!!   IMPLICIT NONE
!! 
!!   COMPLEX  :: rhs_nl1(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   COMPLEX  :: rhs_nl2(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   COMPLEX  :: diff(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   COMPLEX  :: energy_temp1(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   COMPLEX  :: energy_temp2(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   REAL :: nl_error
!!   REAL :: tot_nl1,tot_nl2,norm
!!   REAL :: etot_nl1,etot_nl2
!!   INTEGER :: i,j,k,l,ierr
!! 
!!   IF(np_kz.gt.1) STOP "nl_test not yet implemented for np_kz.g.t1."
!!   IF(np_hank.gt.1) STOP "nl_test not yet implemented for np_hank.g.t1."
!!   IF(np_spec.gt.1) STOP "nl_test not yet implemented for np_spec.g.t1."
!! 
!!   rhs_nl1=cmplx(0.0,0.0)
!!   rhs_nl2=cmplx(0.0,0.0)
!! 
!!   CALL get_phi(g_1)
!!   !WRITE(*,*) "total g_1",REAL(sum(sum(sum(sum(conjg(g_1)*g_1,1),1),1),1))
!!   !WRITE(*,*) "total phi",REAL(sum(sum(sum(conjg(phi)*phi,1),1),1))
!! 
!!   !CALL get_real(g_1)
!!       CALL check_real(g_1,real_err)
!!       IF(real_err.gt.1.0e-10) THEN
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "Error! reality constraint violated. mype,real_err",mype,real_err
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!         WRITE(*,*) "!!!!!!!!!!!!!!!!!!"
!!       END IF
!! 
!!   IF(mype==0) WRITE(*,*) "Calculating nonlinearity (pseudo-spectral)"
!!   CALL get_rhs_nl(g_1,phi,rhs_nl1)
!! 
!!   IF(mype==0) WRITE(*,*) "Calculating nonlinearity (convolution)"
!!   CALL get_rhs_nl_convolution(g_1,phi,rhs_nl2)
!! 
!!   diff=rhs_nl1-rhs_nl2
!! 
!!   tot_nl1=REAL(sum(sum(sum(sum(conjg(rhs_nl1(:,:,:,:,0,0))*rhs_nl1(:,:,:,:,0,0),1),1),1),1))
!!   tot_nl2=REAL(sum(sum(sum(sum(conjg(rhs_nl2(:,:,:,:,0,0))*rhs_nl2(:,:,:,:,0,0),1),1),1),1))
!! 
!!   WRITE(*,*) "mype, total rhs_nl1",mype, tot_nl1
!!   WRITE(*,*) "mype, total rhs_nl2",mype, tot_nl2
!! 
!!   nl_error=REAL(sum(sum(sum(sum(conjg(diff(:,:,:,:,0,0))*diff(:,:,:,:,0,0),1),1),1),1))
!!   norm=REAL(sum(sum(sum(sum(conjg(rhs_nl1(:,:,:,:,0,0))*rhs_nl1(:,:,:,:,0,0),1),1),1),1))
!!   IF(norm.gt.1.0e-10.or.verbose) THEN
!!     !nl_error=nl_error/norm
!!     IF(nl_error/norm.gt.1.0e-10.or.verbose) WRITE(*,*) "!!!!!!!mype, nl_error=",mype, nl_error/norm
!!   ELSE
!!     !WRITE(*,*) "norm.lt.1.0e-10,mype,norm",mype,norm
!!   END IF
!! 
!!   IF(mype==0) THEN
!!     OPEN(unit=100,file=trim(diagdir)//'/nl_test1.dat',status='unknown')
!!     OPEN(unit=101,file=trim(diagdir)//'/nl_test2.dat',status='unknown')
!!     !OPEN(unit=102,file=trim(diagdir)//'/diff.dat',status='unknown')
!!     DO i=0,nkx0-1
!!      DO j=0,nky0-1
!!       DO k=0,nkz0-1
!!        DO l=lv1,lv2
!!          IF(abs(rhs_nl1(i,j,k,l,0,0)).ge.1.0e-13) WRITE(100,'(6es16.8)') REAL(i),REAL(j),REAL(k),&
!!                                    REAL(l),REAL(rhs_nl1(i,j,k,l,0,0)),aimag(rhs_nl1(i,j,k,l,0,0))
!!          IF(abs(rhs_nl2(i,j,k,l,0,0)).ge.1.0e-13) WRITE(101,'(6es16.8)') REAL(i),REAL(j),REAL(k),&
!!                                    REAL(l),REAL(rhs_nl2(i,j,k,l,0,0)),aimag(rhs_nl2(i,j,k,l,0,0))
!!          !IF(abs(rhs_nl2(i,j,k,l)).gt.1.0e-13) THEN
!!             !IF(abs(diff(i,j,k,l))/abs(rhs_nl2(i,j,k,l)).ge.1.0e-15) WRITE(102,'(5es16.8)') REAL(i),REAL(j),REAL(k),REAL(l),&
!!             !      abs(diff(i,j,k,l))/nl_error
!!          !END IF
!!        END DO    
!!       END DO
!!      END DO
!!     END DO
!!    !CLOSE(100)
!!    !CLOSE(101)
!!    !CLOSE(102)
!!   END IF
!! 
!!   !Check energy conservation 
!!   !nl 1
!!   !energy_temp1=cmplx(0.0,0.0)
!!   energy_temp1(:,:,:)=sqrt(pi)*sum(conjg(g_1(:,:,:,:,0,0))*rhs_nl1(:,:,:,:,0,0),4)
!!   IF(mype==0) THEN
!!     DO k=0,nkz0-1
!!       energy_temp1(:,:,k)=energy_temp1(:,:,k)+pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs_nl1(:,:,k,0,0,0)
!!       !energy_temp1(:,:,k)=energy_temp1(:,:,k)+pi**(0.25)*conjg(phi(:,:,k))*rhs_nl1(:,:,k,0)
!!     END DO
!!   END IF
!! 
!!   energy_temp2=cmplx(0.0,0.0)
!!   CALL MPI_ALLREDUCE(energy_temp1,energy_temp2,nkx0*nky0*nkz0 &
!!     ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!!   !energy3d=REAL(energy_temp2)
!!   CALL sum3d_real(REAL(energy_temp2),etot_nl1)
!!   CALL sum3d_real(abs(REAL(energy_temp2)),tot_nl1)
!! 
!!   !Convolution
!!   !energy_temp1=cmplx(0.0,0.0)
!!   energy_temp1(:,:,:)=sqrt(pi)*sum(conjg(g_1(:,:,:,:,0,0))*rhs_nl2(:,:,:,:,0,0),4)
!!   IF(mype==0) THEN
!!     DO k=0,nkz0-1
!!       energy_temp1(:,:,k)=energy_temp1(:,:,k)+pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs_nl2(:,:,k,0,0,0)
!!       !energy_temp1(:,:,k)=energy_temp1(:,:,k)+pi**(0.25)*conjg(phi(:,:,k))*rhs_nl2(:,:,k,0)
!!     END DO
!!   END IF
!!   energy_temp2=cmplx(0.0,0.0)
!!   !Modify
!!   CALL MPI_ALLREDUCE(energy_temp1,energy_temp2,nkx0*nky0*nkz0 &
!!     ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!!   !energy3d=REAL(energy_temp2)
!!   CALL sum3d_real(REAL(energy_temp2),etot_nl2)
!!   CALL sum3d_real(abs(REAL(energy_temp2)),tot_nl2)
!! 
!!   IF(mype==0) THEN
!!     WRITE(100,*) "#NL Energy 1:", etot_nl1/tot_nl1
!!     WRITE(101,*) "#NL Energy 2:", etot_nl2/tot_nl2
!!     CLOSE(100);CLOSE(101)!;CLOSE(102)
!!   END IF
!!   IF(abs(tot_nl1).gt.1.0e-15) THEN
!!     WRITE(*,*) "#(norm) NL Energy 1:", etot_nl1/tot_nl1
!!   ELSE
!!     WRITE(*,*) "# NL Energy 1 (not normalized)",etot_nl1
!!   END IF
!!   IF(abs(tot_nl2).gt.1.0e-15) THEN
!!     WRITE(*,*) "#(norm) NL Energy 2:", etot_nl2/tot_nl2
!!   ELSE
!!     WRITE(*,*) "# NL Energy 2 (not normalized)",etot_nl2
!!   END IF
!!   !WRITE(*,*) "#NL Energy 1:", etot_nl1
!!   !WRITE(*,*) "#NL Energy 2:", etot_nl2
!!   !WRITE(*,*) "#NL ABS 1:", tot_nl1
!!   !WRITE(*,*) "#NL ABS 2:", tot_nl2
!!   
!!   CALL mpi_barrier(mpi_comm_world,ierr)
!!   !WRITE(*,*) "Done with nl_test.",mype
!!   !STOP
!! 
!! END SUBROUTINE nl_test


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           initialize_energy3d                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE initialize_energy3d
!!   IMPLICIT NONE
!! 
!!   CALL get_io_number
!!   energy3d_handle=io_number
!! 
!!   IF(checkpoint_read) THEN
!!     INQUIRE(file=trim(diagdir)//'/energy3d.dat',exist=file_exists)
!!     IF(file_exists) THEN
!!       IF(mype==0) OPEN(unit=energy3d_handle,file=trim(diagdir)//'/energy3d.dat',&
!!                           status='old',form='unformatted',access='stream',position='append')
!!     ELSE
!!       IF(mype==0) OPEN(unit=energy3d_handle,file=trim(diagdir)//'/energy3d.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!     END IF
!!   ELSE
!!     IF(mype==0) OPEN(unit=energy3d_handle,file=trim(diagdir)//'/energy3d.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!   END IF
!! 
!! END SUBROUTINE initialize_energy3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            finalize_energy3d                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE finalize_energy3d
!! 
!!   IF(mype==0) CLOSE(energy3d_handle)
!!   
!! END SUBROUTINE finalize_energy3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              initialize_fmom3d                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE initialize_fmom3d
!!   IMPLICIT NONE
!! 
!!   CALL get_io_number
!!   fmom3d_handle=io_number
!! 
!!   !IF(mype==0) OPEN(unit=fmom3d_handle,file=trim(diagdir)//'/fmom3d.dat',&
!!   !                   status='replace',form='unformatted',access='stream')
!! 
!!   IF(checkpoint_read) THEN
!!     INQUIRE(file=trim(diagdir)//'/fmom3d.dat',exist=file_exists)
!!     IF(file_exists) THEN
!!       IF(mype==0) OPEN(unit=fmom3d_handle,file=trim(diagdir)//'/fmom3d.dat',&
!!                           status='old',form='unformatted',access='stream',position='append')
!!     ELSE
!!       IF(mype==0) OPEN(unit=fmom3d_handle,file=trim(diagdir)//'/fmom3d.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!     END IF
!!   ELSE
!!     IF(mype==0) OPEN(unit=fmom3d_handle,file=trim(diagdir)//'/fmom3d.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!   END IF
!! 
!! 
!! END SUBROUTINE initialize_fmom3d
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                            finalize_fmom3d                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE finalize_fmom3d
!! 
!!   IF(mype==0) CLOSE(fmom3d_handle)
!!   
!! END SUBROUTINE finalize_fmom3d
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                             get_flux_total                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_flux_total(flux_tot)
!!   IMPLICIT NONE
!!   
!!   REAL, INTENT(out) :: flux_tot
!!   REAL :: flux_3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   REAL :: flux_kx0(0:nky0-1,lkz1:lkz2)
!! 
!!   CALL get_flux_3d(flux_3d)
!!   flux_kx0(:,:)=flux_3d(0,:,:)
!!   flux_3d(0,:,:)=0.0
!!   flux_tot=2.0*sum(sum(sum(flux_3d,1),1),1) !factor of 2.0 for -kx modes
!!   flux_tot=flux_tot+sum(sum(flux_kx0,1),1)
!! 
!! END SUBROUTINE get_flux_total
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_flux_3d                                 !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_flux_3d(flux_3d)
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(out) :: flux_3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   COMPLEX :: g2(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   COMPLEX :: g_mu0(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,ls1:ls2)
!!   COMPLEX :: phi_avg(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   REAL    :: pref(0:nkx0-1,0:nky0-1,lh1:lh2)
!!   INTEGER :: n_bcast,ierr
!!   INTEGER :: i,j
!! 
!! 
!!   if(np_spec.gt.1) STOP "get_flux_3d not yet implemented for np_spec.gt.1"
!!   
!!   IF(.not.mu_integrated) THEN
!!    IF (hankel) THEN
!!     n_bcast=2/lv0
!!     pref  = 1.0 
!!     CALL get_phiavg_hk(phi,phi_avg)
!!     CALL integral_hk(conjg(g_1)*phi_avg,pi*pref,g_mu0)
!!     ! Preliminary test  for hankel
!!     IF(lv1.le.2.and.lv2.ge.2) THEN
!!       g2(:,:,:)=g_mu0(:,:,:,2,0)
!!       IF(n_bcast.ne.mype_herm) STOP "Error in get_flux_3d!!!"
!!       DO i=0,nkx0-1
!!          DO j=0,nky0-1
!!             flux_3d(i,j,:)=REAL(-pi**(0.25)*2.0**(-0.5)*i_complex*kygrid(j)*g2(i,j,:))
!!          END DO
!!       END DO
!!     END IF
!!     CALL MPI_BCAST(flux_3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_HERM,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!    ELSE 
!!    n_bcast=2/lv0
!!    pref  = 1.0
!!    CAll get_phiavg(phi,phi_avg)
!!    IF(lv1.le.2.and.lv2.ge.2) THEN
!!       call integral_v(conjg(g_1)*phi_avg,pi*pref,g_mu0)
!!       g2(:,:,:)=g_mu0(:,:,:,2,0)
!!       IF(n_bcast.ne.mype_herm) STOP "Error in get_flux_3d!!!"
!!       DO i=0,nkx0-1
!!          DO j=0,nky0-1
!!             flux_3d(i,j,:)=REAL(-pi**(0.25)*2.0**(-0.5)*i_complex*kygrid(j)*g2(i,j,:))
!!          END DO
!!       END DO
!!     END IF
!!     CALL MPI_BCAST(flux_3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_HERM,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!    ENDIF !hankel condition
!!   ELSE
!!     n_bcast=2/lv0  
!!     IF(lv1.le.2.and.lv2.ge.2) THEN
!!       g2(:,:,:)=g_1(:,:,:,2,0,0)
!!       IF(n_bcast.ne.mype_herm) STOP "Error in get_flux_3d!!!"
!!       DO i=0,nkx0-1
!!        DO j=0,nky0-1
!!         flux_3d(i,j,:)=REAL(-pi**0.25*2.0**(-0.5)*i_complex*kygrid(j)*&
!!              J0a(i,j)*phi(i,j,:)*conjg(g2(i,j,:)))
!!        END DO
!!       END DO
!!     END IF
!!     CALL MPI_BCAST(flux_3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_HERM,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!   ENDIF
!! END SUBROUTINE get_flux_3d
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_press_3d                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_press_3d(press_3d)
!!   IMPLICIT NONE
!! 
!!   COMPLEX, INTENT(out) :: press_3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   COMPLEX :: g_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   COMPLEX :: gn(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   REAL, DIMENSION(0:nkx0-1,0:nky0-1) :: kperp
!!   INTEGER, DIMENSION(1) :: p_kperp
!!   INTEGER :: n_bcast,ierr
!!   INTEGER :: i,j
!! 
!!   kperp = sqrt(2*kperp2)  
!! 
!! 
!!   if(np_hank.gt.1) STOP "get_flux_3d not yet implemented for np_hank.gt.1"
!!   if(np_spec.gt.1) STOP "get_flux_3d not yet implemented for np_spec.gt.1"
!! 
!!   n_bcast=2/lv0  
!!   IF(lv1.le.2.and.lv2.ge.2) THEN
!!     gn(:,:,:)=g_1(:,:,:,2,0,0)
!!     IF(n_bcast.ne.mype) STOP "Error in get_flux_3d!!!"
!!     DO j=0,nky0-1
!!       press_3d(:,j,:)=-pi**0.25*2.0**(-0.5)*gn(:,j,:)
!!     END DO
!!   END IF
!! 
!!   !Modify
!!   CALL MPI_BCAST(press_3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!! 
!!   n_bcast=0
!!   IF(lv1.le.0.and.lv2.ge.0) THEN
!!     gn(:,:,:)=g_1(:,:,:,0,0,0)
!!     IF(n_bcast.ne.mype) STOP "Error in get_flux_3d!!!"
!!     DO j=0,nky0-1
!!       press_3d(:,j,:)=press_3d(:,j,:)-pi**0.25*2.0**(-1.0)*gn(:,j,:)
!!     END DO
!!   END IF
!! 
!!   CALL MPI_BCAST(press_3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!!    
!!   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!! 
!! END SUBROUTINE get_press_3d


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

   INTEGER :: current_time,ierr,max_time
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

  !!!!!!!!!!!!NLT Diagnostics!!!!!!!!!!!!!!!
  !!!!!!!!!!!!NLT Diagnostics!!!!!!!!!!!!!!!
  !!!!!!!!!!!!NLT Diagnostics!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           initialize_diag_nlt                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE initialize_diag_nlt
!!     IMPLICIT NONE
!! 
!!     CHARACTER(len=100), DIMENSION(11) :: file_name
!!     INTEGER :: i
!! 
!!     IF(mype==0.and.verbose) WRITE(*,*) "!!!!!!!!!!!!!!Initializing diag_nlt."
!!     IF(.not.shells_initialized) CALL get_shell_bounds  !Gets the boundaries for the k-shells
!!     
!!     !ALLOCATE NLT, get shell boundaries
!!     ikx_minmax=nkx0-1
!!     iky_minmax=nky0/2-1
!!     ikz_minmax=nkz0/2-1
!!     CALL get_io_number
!!     shell_handle=io_number
!! 
!!     IF(checkpoint_read) THEN
!!       INQUIRE(file=trim(diagdir)//'/nlt_shells.dat',exist=file_exists)
!!       IF(file_exists) THEN
!!         IF(mype==0) OPEN(unit=shell_handle,file=trim(diagdir)//'/nlt_shells.dat',&
!!                             status='old',form='unformatted',access='stream',position='append')
!!       ELSE
!!         IF(mype==0) OPEN(unit=shell_handle,file=trim(diagdir)//'/nlt_shells.dat',&
!!                             status='replace',form='unformatted',access='stream')
!!       END IF
!!     ELSE
!!       IF(mype==0) OPEN(unit=shell_handle,file=trim(diagdir)//'/nlt_shells.dat',&
!!                             status='replace',form='unformatted',access='stream')
!!     END IF
!! 
!!     CALL get_io_number
!!     nlt_status=io_number
!!     IF(mype==0) OPEN(unit=nlt_status,file=trim(diagdir)//'/nlt_status.dat',status='replace')
!!    
!!     IF(output_nlt_n.and.mype==0) THEN
!! 
!!       !These are nlt diagnostics for specifice n's
!!       file_name(1)='/nlt_shells_n0.dat'   !n=0 (i.e. phi)
!!       file_name(2)='/nlt_shells_n1.dat'   !n=1
!!       file_name(3)='/nlt_shells_n2.dat'   !n=2
!!       file_name(4)='/nlt_shells_n3.dat'   !n=3
!!       file_name(5)='/nlt_shells_n1o16.dat' !n=(1/8)nmax
!!       file_name(6)='/nlt_shells_n2o16.dat' !n=(2/8)nmax
!!       file_name(7)='/nlt_shells_n3o16.dat' !n=(3/8)nmax
!!       file_name(8)='/nlt_shells_n4o16.dat' !n=(4/8)nmax
!!       file_name(9)='/nlt_shells_n5o16.dat' !n=(5/8)nmax
!!       file_name(10)='/nlt_shells_n6o16.dat' !n=(6/8)nmax
!!       file_name(11)='/nlt_shells_n7o16.dat' !n=(7/8)nmax
!! 
!!       !n=0
!!       DO i=1,11
!!         CALL get_io_number
!!         shell_handle_n(i)=io_number
!!         IF(checkpoint_read) THEN
!!           INQUIRE(file=trim(diagdir)//trim(file_name(i)),exist=file_exists)
!!           IF(file_exists) THEN
!!             IF(mype==0) OPEN(unit=shell_handle_n(i),file=trim(diagdir)//trim(file_name(i)),&
!!                                 status='old',form='unformatted',access='stream',position='append')
!!           ELSE
!!             IF(mype==0) OPEN(unit=shell_handle_n(i),file=trim(diagdir)//trim(file_name(i)),&
!!                                 status='replace',form='unformatted',access='stream')
!!           END IF
!!         ELSE
!!           IF(mype==0) OPEN(unit=shell_handle_n(i),file=trim(diagdir)//trim(file_name(i)),&
!!                                 status='replace',form='unformatted',access='stream')
!!         END IF
!! 
!!       END DO
!! 
!!     END IF
!! 
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     !CALL get_io_number
!!     !test_handle1=io_number
!!     !IF(mype==0) OPEN(unit=test_handle1,file=trim(diagdir)//'/nlt_test1.dat',status='replace')
!!     !CALL get_io_number
!!     !test_handle2=io_number
!!     !IF(mype==0) OPEN(unit=test_handle2,file=trim(diagdir)//'/nlt_test2.dat',status='replace')
!!     !IF(mype==0) WRITE(test_handle2,*) "# time,dedt,elin,enl,nl from nlt"
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!! 
!!     IF(mype==0.and.verbose) WRITE(*,*) "Done initializing diag_nlt."
!!    
!!   END SUBROUTINE initialize_diag_nlt
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                          finalize_diag_nlt                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE finalize_diag_nlt
!! 
!!     INTEGER :: i
!! 
!!     IF(mype==0) CLOSE(shell_handle)
!!     !IF(mype==0) CLOSE(nlt_handle)
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     IF(mype==0) CLOSE(test_handle1)
!!     IF(mype==0) CLOSE(test_handle2)
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!!     !!!!!!!!Testing!!!!!!!!
!! 
!! 
!!     IF(output_nlt_n.and.mype==0) THEN
!!       DO i=1,11
!!         CLOSE(shell_handle_n(i))
!!       END DO 
!!     END IF
!! 
!!   END SUBROUTINE finalize_diag_nlt
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                           initialize_diag_nlt_triple                      !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE initialize_diag_nlt_triple
!!     IMPLICIT NONE
!! 
!!     CHARACTER(len=100), DIMENSION(11) :: file_name
!!     INTEGER :: i
!! 
!!     IF(mype==0.and.verbose) WRITE(*,*) "!!!!!!!!!!!!!!Initializing diag_nlt_triple."
!!     IF(.not.shells_initialized) CALL get_shell_bounds  !Gets the boundaries for the k-shells
!!     
!!     !ALLOCATE NLT, get shell boundaries
!!     ikx_minmax=nkx0-1
!!     iky_minmax=nky0/2-1
!!     ikz_minmax=nkz0/2-1
!!     CALL get_io_number
!!     nlt3_handle=io_number
!! 
!!     IF(checkpoint_read) THEN
!!       INQUIRE(file=trim(diagdir)//'/nlt_triple.dat',exist=file_exists)
!!       IF(file_exists) THEN
!!         IF(mype==0) OPEN(unit=nlt3_handle,file=trim(diagdir)//'/nlt_triple.dat',&
!!                             status='old',form='unformatted',access='stream',position='append')
!!       ELSE
!!         IF(mype==0) OPEN(unit=nlt3_handle,file=trim(diagdir)//'/nlt_triple.dat',&
!!                             status='replace',form='unformatted',access='stream')
!!       END IF
!!     ELSE
!!       IF(mype==0) OPEN(unit=nlt3_handle,file=trim(diagdir)//'/nlt_triple.dat',&
!!                             status='replace',form='unformatted',access='stream')
!!     END IF
!! 
!!     IF(mype==0.and.verbose) WRITE(*,*) "Done initializing diag_nlt_triple."
!!    
!!   END SUBROUTINE initialize_diag_nlt_triple
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                          finalize_diag_nlt_triple                         !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE finalize_diag_nlt_triple
!! 
!!     IF(mype==0) CLOSE(nlt3_handle)
!! 
!!   END SUBROUTINE finalize_diag_nlt_triple
!! 
!! 
!! 
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              get_shell_bounds                             !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE get_shell_bounds
!! 
!!     IMPLICIT NONE
!!     REAL :: mwid,area,kmax,kperp,kbound0,dkperp
!!     REAL, ALLOCATABLE, DIMENSION(:) :: shell_test
!!     INTEGER :: i,j,s,count0
!!     LOGICAL ::  another_shell
!!     INTEGER :: shell_version
!!     
!! 
!!     shell_version=1
!!     IF(shell_version==1) THEN
!!         another_shell=.true.
!!         kmax=min(kxmax,kymax) 
!!         count0=0
!!         dkperp=0.2              
!!         kperp=dkperp
!!         DO WHILE(another_shell)
!!           kperp=kperp+dkperp
!!           IF(kperp.ge.kmax) another_shell=.false.
!!           count0=count0+1
!!         END DO
!!         ALLOCATE(shell_bounds(count0))
!!         num_shells=count0+1
!!         IF(mype==0) OPEN(unit=900,file=trim(diagdir)//'/shell_info.dat',status='replace')
!!         IF(mype==0) WRITE(900,*) "#k_perp shell boundaries"
!!         kperp=dkperp
!!         DO i=1,count0
!!           shell_bounds(i)=kperp
!!           IF(mype==0) WRITE(900,*) shell_bounds(i)
!!           kperp=kperp+dkperp
!!         END DO
!!         IF(mype==0) CLOSE(900)
!!         shells_initialized=.true.
!! 
!! 
!! 
!!     ELSEIF(shell_version==2) THEN
!!         another_shell=.true.
!!         kmax=min(kxmax,kymax) 
!!         count0=0
!!         kperp=0.8               
!!         DO WHILE(another_shell)
!!           kperp=kperp*2.0**0.2
!!           IF(kperp.ge.kmax) another_shell=.false.
!!           count0=count0+1
!!         END DO
!!         ALLOCATE(shell_bounds(count0))
!!         num_shells=count0+1
!!         IF(mype==0) OPEN(unit=900,file=trim(diagdir)//'/shell_info.dat',status='replace')
!!         IF(mype==0) WRITE(900,*) "#k_perp shell boundaries"
!!         kperp=0.8
!!         DO i=1,count0
!!           shell_bounds(i)=kperp
!!           IF(mype==0) WRITE(900,*) shell_bounds(i)
!!           kperp=kperp*2.0**0.2
!!         END DO
!!         IF(mype==0) CLOSE(900)
!!         shells_initialized=.true.
!! 
!!     ELSEIF(shell_version==3) THEN
!! 
!!       shells_initialized=.true.
!!       !IF(mype==0) WRITE(*,*) "In get_shell_bounds."
!!       mwid=min_shell_width*kxmin
!!       kmax=MIN(kxmax,kymax)
!!       area=2.0*pi*(kmax-mwid/2.0)*mwid
!!       kbound0=(area/pi)**0.5
!!       count0=0
!! 
!!       DO WHILE(kbound0.lt.kmax+0.5*mwid)
!!         kbound0=(kbound0**2+area/pi)**0.5
!!         count0=count0+1
!!       END DO
!!   
!!       ALLOCATE(shell_bounds(count0))
!!       num_shells=count0+1
!!       ALLOCATE(shell_test(num_shells))
!!       shell_test=0.0
!!       IF(mype==0) OPEN(unit=900,file=trim(diagdir)//'/shell_info.dat',status='replace')
!!   
!!       shell_bounds(1)=(area/pi)**0.5
!!       IF(mype==0) WRITE(900,*) "#k_perp shell boundaries"
!!       IF(mype==0) WRITE(900,*) shell_bounds(1)
!!   
!!       DO i=2,count0
!!         shell_bounds(i)=(shell_bounds(i-1)**2+area/pi)**0.5
!!         IF(mype==0) WRITE(900,*) shell_bounds(i)
!!       END DO
!!   
!!       DO i=0,nkx0-1
!!         DO j=0,nky0-1
!!           kperp=sqrt(kxgrid(i)**2+kygrid(j)**2) 
!!           IF(kperp.lt.shell_bounds(1)) THEN
!!             shell_test(1)=shell_test(1)+1
!!           ELSE IF (kperp.ge.shell_bounds(num_shells-1)) THEN
!!             shell_test(num_shells)=shell_test(num_shells)+1
!!           ELSE
!!             DO s=1,num_shells-2
!!               IF(kperp.ge.shell_bounds(s).and.kperp.lt.shell_bounds(s+1)) THEN
!!                 shell_test(s+1)=shell_test(s+1)+1
!!               END IF
!!             END DO
!!           END IF
!!         END DO
!!       END DO
!!       IF(mype==0) WRITE(900,*) "#Number of k's per shell:"
!!       DO i=1,num_shells
!!         IF(mype==0) WRITE(900,*) i,shell_test(i)
!!       END DO
!!   
!!       IF(mype==0) CLOSE(900)
!!     ELSE
!!       STOP "Error in get_shell_bounds!"
!!     END IF
!!   
!!   END SUBROUTINE get_shell_bounds
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_nlt_shells                              !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!  With filter
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE get_nlt_shells
!!     IMPLICIT NONE
!!     
!!     INTEGER :: i,j,k
!!     INTEGER :: ip,jp,kp
!!     !INTEGER :: shell_index
!!     !INTEGER :: shell_index_p
!!     REAL :: kxp,kyp,kzp
!!     INTEGER :: k0_nlt,kg,ierr
!!     COMPLEX :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!     COMPLEX :: phi_f(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     REAL :: sbounds(num_shells+1)
!!     REAL :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     REAL :: NLT_kxn0,NLT_kx0
!!     REAL :: NLT_temp
!!     LOGICAL :: nlt_filled(num_shells,num_shells,0:ikz_minmax,0:ikz_minmax)
!!     INTEGER :: nlt_out_ns(11)
!!     INTEGER :: nout_index
!! 
!!     if(np_hank.gt.1) STOP "get_nlt_shells not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_nlt_shells not yet implemented for np_spec.gt.1"
!!     if(np_kz.gt.1) STOP "get_nlt_shells not yet implemented for np_kz.gt.1"
!! 
!!     IF(output_nlt_n) THEN
!!       nlt_out_ns(1)=0
!!       nlt_out_ns(2)=1
!!       nlt_out_ns(3)=2
!!       nlt_out_ns(4)=3
!!       nlt_out_ns(5)=nv0/16
!!       nlt_out_ns(6)=2*nv0/16
!!       nlt_out_ns(7)=3*nv0/16
!!       nlt_out_ns(8)=4*nv0/16
!!       nlt_out_ns(9)=5*nv0/16
!!       nlt_out_ns(10)=6*nv0/16
!!       nlt_out_ns(11)=7*nv0/16
!!     END IF
!! 
!!     IF(num_shells-1.ne.size(shell_bounds)) THEN
!!         WRITE(*,*) "Error in get_nlt_shells!"
!!         STOP
!!     END IF
!!     sbounds(2:num_shells)=shell_bounds
!!     sbounds(1)=0.0
!!     sbounds(num_shells+1)=100000.0
!!     IF(mype==0.and.verbose) WRITE(*,*) "In get_nlt_shells."
!!     !IF(mype==0) WRITE(nlt_handle) time
!!     NLT=0.0
!!     IF(output_nlt_n) NLT_n=0.0
!! 
!!     IF(mype==0) WRITE(nlt_status,*) time
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     nlt_filled=.false.
!!     DO ip=1,num_shells
!!       IF(mype==0) WRITE(*,*) "get_nlt_shells:g",ip
!!       IF(mype==0) WRITE(nlt_status,*) "get_nlt_shells:g",ip,num_shells
!!       IF(mype==0) flush(nlt_status)
!!       DO kp=0,ikz_minmax
!! 
!!         !First filter the distribution function so that it only has elements in
!!         !the appropriate kperp shell along with +-kp.  Note that we must USE
!!         !+-kp since the RHS calculation implicitly uses the -kp modes for the
!!         !kx<0 modes (i.e., the reality constraint gives kx<0,kp=-kp modes)
!!         CALL g_shell_filter(sbounds(ip),sbounds(ip+1),kp*kzmin,g_f)
!!         CALL get_phi(g_1)
!!         rhs=cmplx(0.0,0.0)
!!         CALL get_rhs_nl(g_f,phi,rhs)
!! 
!!         DO i=1,num_shells
!!           DO k=0,ikz_minmax
!!             !Don't DO more work than necessary: Fill all possible modes by
!!             !applying symmetry considerations
!!             !(Note: I've already tested explicity that nothing changes)
!!             IF(.not.nlt_filled(i,ip,k,kp)) THEN
!!               CALL g_shell_filter(sbounds(i),sbounds(i+1),k*kzmin,g_f)
!! 
!!               NLT_kxn0=REAL(sqrt(pi)*sum(sum(sum(sum(conjg(g_f(1:nkx0-1,:,:,:,0,0))*rhs(1:nkx0-1,:,:,:,0,0),1),1),1),1))
!!               CALL MPI_ALLREDUCE(NLT_kxn0,NLT_temp,1 &
!!                 ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!! 
!!               NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+2.0*NLT_temp
!!               !NLT(i,ip,-k,-kp)=NLT(i,ip,-k,-kp)+2.0*NLT_temp
!! 
!!               NLT_kx0=REAL(sqrt(pi)*sum(sum(sum(conjg(g_f(0,:,:,:,0,0))*rhs(0,:,:,:,0,0),1),1),1))
!!               CALL MPI_ALLREDUCE(NLT_kx0,NLT_temp,1 &
!!                 ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!               NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+NLT_temp
!!               nlt_filled(i,ip,k,kp)=.true.
!!               NLT(ip,i,kp,k)=-NLT(i,ip,k,kp)
!!               nlt_filled(ip,i,kp,k)=.true.
!! 
!!               IF(output_nlt_n) THEN
!!                 DO nout_index=1,11
!!                 !nlt_out_ns(1)=0
!!                   IF((nlt_out_ns(nout_index).ge.lv1).and.(nlt_out_ns(nout_index).le.lv2)) THEN
!!                     NLT_kxn0=REAL(sqrt(pi)*sum(sum(sum(conjg(g_f(1:nkx0-1,:,:,nlt_out_ns(nout_index),0,0))&
!!                                                     *rhs(1:nkx0-1,:,:,nlt_out_ns(nout_index),0,0),1),1),1))
!! 
!!                     NLT_n(i,ip,k,kp,nout_index)=NLT_n(i,ip,k,kp,nout_index)+2.0*NLT_kxn0
!!                     !NLT(i,ip,-k,-kp)=NLT(i,ip,-k,-kp)+2.0*NLT_temp
!! 
!!                     NLT_kx0=REAL(sqrt(pi)*sum(sum(conjg(g_f(0,:,:,nlt_out_ns(nout_index),0,0))&
!!                                                   *rhs(0,:,:,nlt_out_ns(nout_index),0,0),1),1))
!! 
!!                     NLT_n(i,ip,k,kp,nout_index)=NLT_n(i,ip,k,kp,nout_index)+NLT_kx0
!!                     NLT_n(ip,i,kp,k,nout_index)=-NLT_n(i,ip,k,kp,nout_index)
!!                     
!!                   END IF !IF correct processor
!! 
!!                 END DO
!!               END IF
!! 
!!             END IF
!! 
!!           END DO
!!         END DO
!!         
!!       END DO
!!     END DO
!!     
!! 
!!     nlt_filled=.false.
!!     DO ip=1,num_shells
!!       IF(mype==0) WRITE(*,*) "get_nlt_shells:phi",ip
!!       IF(mype==0) WRITE(nlt_status,*) "get_nlt_shells:phi",ip,num_shells
!!       IF(mype==0) flush(nlt_status)
!!       !IF(verbose) WRITE(*,*) "get_nlt_shells:phi",ip,mype
!!       DO kp=0,ikz_minmax
!!         CALL get_phi(g_1)
!!         CALL phi_shell_filter(sbounds(ip),sbounds(ip+1),kp*kzmin,phi,phi_f)
!!         !CALL test_shell_filter(phi_f)
!!         rhs=cmplx(0.0,0.0)
!!         CALL get_rhs_nl(g_1,phi_f,rhs)
!!         DO i=1,num_shells
!!           DO k=0,ikz_minmax
!!              IF(.not.nlt_filled(i,ip,k,kp)) THEN
!!               !CALL get_phi(g_1)
!!               CALL phi_shell_filter(sbounds(i),sbounds(i+1),k*kzmin,phi,phi_f)
!!               IF(mype==0) THEN
!! 
!!                 phi_f(:,:,k)=phi_f(:,:,k)*J0a(:,:)
!!                 IF(k.ne.0) phi_f(:,:,nkz0-k)=phi_f(:,:,nkz0-k)*J0a(:,:)
!!                 NLT_kxn0=REAL((pi)**0.25*sum(sum(sum(conjg(phi_f(1:nkx0-1,:,:))*rhs(1:nkx0-1,:,:,0,0,0),1),1),1))
!!                 NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+2.0*NLT_kxn0
!!                 IF(output_nlt_n) NLT_n(i,ip,k,kp,1)=NLT_n(i,ip,k,kp,1)+2.0*NLT_kxn0
!!                 NLT_kx0=REAL((pi)**0.25*sum(sum(conjg(phi_f(0,:,:))*rhs(0,:,:,0,0,0),1),1))
!!                 NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+NLT_kx0
!!                 IF(output_nlt_n) NLT_n(i,ip,k,kp,1)=NLT_n(i,ip,k,kp,1)+NLT_kx0
!! 
!!                 nlt_filled(i,ip,k,kp)=.true.
!!                 NLT(ip,i,kp,k)=-NLT(i,ip,k,kp)
!!                 IF(output_nlt_n) NLT_n(ip,i,kp,k,1)=-NLT_n(i,ip,k,kp,1)
!!                 nlt_filled(ip,i,kp,k)=.true.
!! 
!!    
!!             END IF
!!             END IF
!!           END DO
!!         END DO
!! 
!!       END DO
!!     END DO
!! 
!!     DEALLOCATE(rhs)
!! 
!!     IF(mype==0) WRITE(shell_handle) time
!!     IF(mype==0) WRITE(shell_handle) NLT
!!     IF(mype==0) flush(shell_handle)
!! 
!!     IF(output_nlt_n) THEN
!!       DO nout_index=1,11 
!!         CALL MPI_ALLREDUCE(NLT_n(1,1,0,0,nout_index),NLT(1,1,0,0),num_shells**2*(ikz_minmax+1)**2 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!         IF(mype==0) WRITE(shell_handle_n(nout_index)) time
!!         IF(mype==0) WRITE(shell_handle_n(nout_index)) NLT
!!         IF(mype==0) flush(shell_handle_n(nout_index))
!!       END DO
!!     END IF
!! 
!!   END SUBROUTINE get_nlt_shells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_shell_index                             !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE get_shell_index(kx,ky,shell_index)
!!     IMPLICIT NONE
!!     REAL, INTENT(in) :: kx,ky
!!     INTEGER, INTENT(out) :: shell_index
!!     REAL :: kperp
!!     INTEGER :: s
!!     
!!     kperp=sqrt(kx**2+ky**2)
!!     IF(kperp.lt.shell_bounds(1)) THEN
!!       shell_index=1
!!     ELSE IF (kperp.ge.shell_bounds(num_shells-1)) THEN
!!       shell_index=num_shells
!!     ELSE
!!       DO s=1,num_shells-2
!!         IF((kperp.ge.shell_bounds(s)).and.(kperp.lt.shell_bounds(s+1))) THEN
!!           shell_index=s+1
!!         END IF
!!       END DO
!!     END IF
!! 
!!   END SUBROUTINE get_shell_index
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                                  get_nlt_k                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE get_nlt_k(ikx_in,iky_in,ikz_in)
!!     IMPLICIT NONE
!! 
!!     INTEGER, INTENT(in) :: ikx_in,iky_in,ikz_in
!!     REAL :: kx0,ky0,kz0
!!     REAL :: kxp,kyp,kzp
!!     REAL :: kxpp,kypp,kzpp
!!     INTEGER :: ikxp,ikyp,ikzp
!!     INTEGER :: ikxpp,ikypp,ikzpp
!!     !REAL, INTENT(out) :: nltk(-ikx_minmax:ikx_minmax,iky_minmax:iky_minmax,&
!!     !                          -ikz_minmax:ikz_minmax)
!!     COMPLEX :: phibk0
!!     COMPLEX :: phibkp
!!     COMPLEX :: phibkpp
!!     COMPLEX :: gk0(lv1:lv2)
!!     COMPLEX :: gkp(lv1:lv2)
!!     COMPLEX :: gkpp(lv1:lv2)
!!     INTEGER :: i,j,k
!!     LOGICAL :: conjg_p,conjg_pp
!!     !INTEGER :: ikx_minmax
!!     !INTEGER :: iky_minmax
!!     !INTEGER :: ikz_minmax
!!     REAL :: ckkp
!!     REAL :: nltv(lv1:lv2),nlt0_temp,nlt0
!!     INTEGER :: ierr
!! 
!!     if(np_hank.gt.1) STOP "get_nlt_k not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_nlt_k not yet implemented for np_spec.gt.1"
!!     if(spatial2d) STOP "get_nlt_k not yet implemented for spatial2d"
!!     nltk=0.0
!! 
!!     kx0=kxgrid(ikx_in)
!!     ky0=kygrid(iky_in)
!!     kz0=kzgrid(ikz_in)
!!     CALL get_phi(g_1)
!!     phibk0=J0a(ikx_in,iky_in)*phi(ikx_in,iky_in,ikz_in)
!!     gk0(:)=g_1(ikx_in,iky_in,ikz_in,:,0,0)
!!     DO i=-ikx_minmax,ikx_minmax
!!       !IF(mype==0) WRITE(*,*) "i",i
!!       DO j=-iky_minmax,iky_minmax
!!         DO k=-ikz_minmax,ikz_minmax
!!            kxp=i*kxmin
!!            kyp=j*kymin
!!            kzp=k*kzmin
!!            kxpp=kx0-kxp
!!            kypp=ky0-kyp
!!            kzpp=kz0-kzp
!! 
!!            IF((abs(kxpp).le.kxmax).and.&
!!               (abs(kypp).le.kymax).and.&
!!               (abs(kzpp).le.kzmax))THEN
!! 
!!              CALL get_indices_from_ks(kxp,kyp,kzp,ikxp,ikyp,ikzp,conjg_p)
!!              CALL get_indices_from_ks(kxpp,kypp,kzpp,ikxpp,ikypp,ikzpp,conjg_pp)
!!              ckkp=kxp*ky0-kx0*kyp
!!              IF(conjg_p) THEN
!!                phibkp=conjg(J0a(ikxp,ikyp)*phi(ikxp,ikyp,ikzp))
!!                gkp(:)=conjg(g_1(ikxp,ikyp,ikzp,:,0,0))
!!              ELSE
!!                phibkp=(J0a(ikxp,ikyp)*phi(ikxp,ikyp,ikzp))
!!                gkp(:)=(g_1(ikxp,ikyp,ikzp,:,0,0))
!!              END IF
!!              IF(conjg_pp) THEN
!!                phibkpp=conjg(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!                gkpp(:)=conjg(g_1(ikxpp,ikypp,ikzpp,:,0,0))
!!              ELSE
!!                phibkpp=(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!                gkpp(:)=(g_1(ikxpp,ikypp,ikzpp,:,0,0))
!!              END IF
!!              nltv(:)= REAL(-pi**0.5*ckkp*conjg(gk0(:))*phibkpp*gkp(:))
!!              IF(mype==0) nltv(0)=nltv(0)+REAL(pi**0.25*ckkp*(conjg(phibk0)*phibkp*gkpp(0)))
!! 
!!              nlt0_temp=sum(nltv)
!!              CALL mpi_allreduce(nlt0_temp,nlt0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
!!              nltk(i,j,k)=nlt0
!!            END IF
!!         END DO
!!       END DO
!!     END DO
!! 
!!   END SUBROUTINE get_nlt_k


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             g_shell_filter                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE g_shell_filter(kperp_min,kperp_max,kz,g_f)
!!     IMPLICIT NONE
!! 
!!     REAL, INTENT(in) :: kperp_min,kperp_max,kz
!!     COMPLEX, INTENT(out) :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!     INTEGER :: k
!!     INTEGER :: i,j
!!     REAL :: kperp
!!     INTEGER :: i0,j0,k0
!!     LOGICAL :: conjg_p
!!     INTEGER :: neg_k
!! 
!!     k=nint(kz/kzmin)
!!     IF(kz.lt.0.0) THEN
!!         k=k+nkz0
!!     END IF
!!     !Fill in both kz and -kz  (note that this satisfies reality constraint for
!!     !kx=0 modes)
!!     IF(spatial2d) THEN
!!       CALL get_indices_from_ks(0.0,0.0,-kzgrid(j),i0,j0,k0,conjg_p)
!!     ELSE
!!       CALL get_indices_from_ks(0.0,0.0,-kzgrid(k),i0,j0,k0,conjg_p)
!!     END IF
!!     neg_k=k0
!! 
!!     g_f=cmplx(0.0,0.0)
!!     DO i=0,nkx0-1
!!     !DO i=1,nkx0-1
!!       DO j=0,nky0-1
!!         IF(j.ne.nky0/2) THEN
!!           kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!           IF((kperp.ge.kperp_min).and.(kperp.lt.kperp_max)) THEN
!!             g_f(i,j,k,:,:,:)=g_1(i,j,k,:,:,:)
!!             g_f(i,j,neg_k,:,:,:)=g_1(i,j,neg_k,:,:,:)
!!           END IF
!!         END IF
!!       END DO
!!     END DO
!! 
!!   END SUBROUTINE g_shell_filter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               g_shell_filter3                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE g_shell_filter3(kperp_min,kperp_max,kz,g_f)
!    IMPLICIT NONE
!
!    REAL, INTENT(in) :: kperp_min,kperp_max,kz
!    COMPLEX, INTENT(out) :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!    INTEGER :: k
!    INTEGER :: i,j
!    REAL :: kperp
!
!    k=nint(kz/kzmin)
!    IF(kz.lt.0.0) THEN
!        k=k+nkz0
!    END IF
!
!    g_f=cmplx(0.0,0.0)
!    DO i=0,nkx0-1
!    !DO i=1,nkx0-1
!      DO j=0,nky0-1
!        IF(j.ne.nky0/2) THEN
!          kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!          IF((kperp.ge.kperp_min).and.(kperp.lt.kperp_max)) THEN
!            g_f(i,j,k,:,:,:)=g_1(i,j,k,:,:,:)
!          END IF
!        END IF
!      END DO
!    END DO
!
!  END SUBROUTINE g_shell_filter3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              phi_shell_filter                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE phi_shell_filter(kperp_min,kperp_max,kz,phi_in,phi_f)
!!     IMPLICIT NONE
!! 
!!     REAL, INTENT(in) :: kperp_min,kperp_max,kz
!!     COMPLEX, INTENT(in) :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!     COMPLEX, INTENT(out) :: phi_f(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!     INTEGER :: k
!!     INTEGER :: i,j
!!     REAL :: kperp
!!     INTEGER :: i0,j0,k0
!!     LOGICAL :: conjg_p
!!     INTEGER :: neg_k
!!     
!!     !CALL get_phi(g_1)
!!     k=nint(kz/kzmin)
!!     IF(kz.lt.0.0) THEN
!!         k=k+nkz0
!!     END IF
!!     IF(spatial2d) THEN
!!       CALL get_indices_from_ks(0.0,0.0,-kzgrid(j),i0,j0,k0,conjg_p)
!!     ELSE
!!       CALL get_indices_from_ks(0.0,0.0,-kzgrid(k),i0,j0,k0,conjg_p)
!!     END IF
!!     neg_k=k0
!! 
!!     phi_f=cmplx(0.0,0.0)
!!     DO i=0,nkx0-1
!!       DO j=0,nky0-1
!!           kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!           IF((kperp.ge.kperp_min.and.kperp.lt.kperp_max)) THEN
!!             phi_f(i,j,k)=phi_in(i,j,k)
!!             phi_f(i,j,neg_k)=phi_in(i,j,neg_k)
!!           END IF
!!       END DO
!!     END DO
!! 
!!   END SUBROUTINE phi_shell_filter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             phi_shell_filter3                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE phi_shell_filter3(kperp_min,kperp_max,kz,phi_f)
!    IMPLICIT NONE
!
!    REAL, INTENT(in) :: kperp_min,kperp_max,kz
!    COMPLEX, INTENT(out) :: phi_f(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!    INTEGER :: k
!    INTEGER :: i,j
!    REAL :: kperp
!    
!    CALL get_phi(g_1)
!
!    k=nint(kz/kzmin)
!    IF(kz.lt.0.0) THEN
!        k=k+nkz0
!    END IF
!
!    phi_f=cmplx(0.0,0.0)
!    DO i=0,nkx0-1
!      DO j=0,nky0-1
!          kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!          IF((kperp.ge.kperp_min.and.kperp.lt.kperp_max)) THEN
!            phi_f(i,j,k)=phi(i,j,k)
!          END IF
!      END DO
!    END DO
!
!  END SUBROUTINE phi_shell_filter3
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           phi_shell_filter_wo0                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE phi_shell_filter_wo0(kperp_min,kperp_max,kz,phi_f)
!!    IMPLICIT NONE
!!
!!    REAL, INTENT(in) :: kperp_min,kperp_max,kz
!!    COMPLEX, INTENT(out) :: phi_f(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!    INTEGER :: k
!!    INTEGER :: i,j
!!    REAL :: kperp
!!    
!!    CALL get_phi(g_1)
!!
!!    k=nint(kz/kzmin)
!!    IF(kz.lt.0.0) THEN
!!        k=k+nkz0
!!    END IF
!!
!!    phi_f=cmplx(0.0,0.0)
!!    DO i=1,nkx0-1
!!      DO j=0,nky0-1
!!          kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!          IF((kperp.ge.kperp_min.and.kperp.lt.kperp_max)) THEN
!!            phi_f(i,j,k)=phi(i,j,k)
!!          END IF
!!      END DO
!!    END DO
!!
!!  END SUBROUTINE phi_shell_filter_wo0



!!!!!!!!!!!Testing!!!!!!!!!!!!!
!!!!!!!!!!!Testing!!!!!!!!!!!!!
!!!!!!!!!!!Testing!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               get_eterms_k                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE get_eterms_k(ikx_in,iky_in,ikz_in,dedt,elin,enl)
!!    IMPLICIT NONE
!!    
!!    INTEGER, INTENT(in) :: ikx_in,iky_in,ikz_in
!!    REAL, INTENT(out) :: dedt,elin,enl
!!    REAL :: energy3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!    
!!    CALL get_energy3d(g_1,energy3d,3) !Total RHS 
!!    dedt=energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,4) !Drive
!!    elin=energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,5) !Coll
!!    elin=elin+energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,6) !HColl
!!    elin=elin+energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,7) !Artificial
!!    elin=elin+energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,9) !hyp_conv
!!    elin=elin+energy3d(ikx_in,iky_in,ikz_in)
!!    CALL get_energy3d(g_1,energy3d,8) !NL
!!    enl=energy3d(ikx_in,iky_in,ikz_in)
!!
!!  END SUBROUTINE get_eterms_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              nlt_test                                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE nlt_test
!!     IMPLICIT NONE
!!     INTEGER :: ikx1,iky1,ikz1
!!     INTEGER :: ikx2,iky2,ikz2
!!     !REAL :: nltk(-1*(nkx0-1):nkx0-1,-1*(nky0/2-1):nky0/2-1,&
!!     !                          -1*(nkz0/2-1):nkz0/2-1)
!!     REAL :: nlt12,nlt21
!!     REAL :: dedt,elin,enl
!!     INTEGER :: i,j,k
!! 
!!     ikx1=0
!!     iky1=5
!!     ikz1=1
!! 
!!     ikx2=3
!!     iky2=3
!!     ikz2=0
!! 
!!     CALL get_nlt_k(ikx2,iky2,ikz2)
!!     nlt21=nltk(ikx1,iky1,ikz1)
!!     CALL get_nlt_k(ikx1,iky1,ikz1)
!!     nlt12=nltk(ikx2,iky2,ikz2)
!! 
!!     !!!!!Temp
!!     !IF(mype==0) THEN
!!     !  OPEN(unit=555,file=trim(diagdir)//'/temp.dat',status='unknown')
!!     !  DO i=-ikx_minmax,ikx_minmax
!!     !    DO j=-iky_minmax,iky_minmax
!!     !      DO k=-ikz_minmax,ikz_minmax
!!     !        WRITE(555,*) i,j,k,nltk(i,j,k)
!!     !      END DO
!!     !    END DO
!!     !  END DO
!!     !  CLOSE(555)
!!     !END IF
!!     !!!!!Temp
!! 
!!     IF(mype==0) WRITE(test_handle1,*) time,nlt12,nlt21
!! 
!!     CALL get_eterms_k(ikx1,iky1,ikz1,dedt,elin,enl)
!! 
!!     IF(mype==0) WRITE(test_handle2,'(5es16.8)') &
!!         time,dedt,elin,enl,sum(sum(sum(nltk,1),1),1)
!! 
!! 
!!   END SUBROUTINE nlt_test
!! 
!! 
!! 
!! !!!!!!!!!!!Energetics!!!!!!!!!!!!!
!! !!!!!!!!!!!Energetics!!!!!!!!!!!!!
!! !!!!!!!!!!!Energetics!!!!!!!!!!!!!
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                            get_energy_terms                               !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_energy_terms
!! 
!!   REAL :: energy3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   REAL :: energy_tot(10)
!!   !REAL :: phi_tot_temp
!! 
!!   !get entropy
!!   !IF(verbose) WRITE(*,*) "Before get_energy3d."
!!   CALL get_energy3d(g_1,energy3d,1) 
!!   !IF(verbose) WRITE(*,*) "After get_energy3d."
!!   !IF(verbose) WRITE(*,*) "Before sum3d_real."
!!   CALL sum3d_real(energy3d,energy_tot(1))
!!   !IF(verbose) WRITE(*,*) "After sum3d_real."
!!   !get electrostatic part
!!   !IF(verbose) WRITE(*,*) "Before 2."
!!   CALL get_energy3d(g_1,energy3d,2) 
!!   CALL sum3d_real(energy3d,energy_tot(2))
!!   !get total rhs
!!   !IF(verbose) WRITE(*,*) "Before 3."
!!   CALL get_energy3d(g_1,energy3d,3) 
!!   CALL sum3d_real(energy3d,energy_tot(3))
!!   !get drive
!!   !IF(verbose) WRITE(*,*) "Before 4."
!!   CALL get_energy3d(g_1,energy3d,4) 
!!   CALL sum3d_real(energy3d,energy_tot(4))
!!   !get collisions
!!   !IF(verbose) WRITE(*,*) "Before 5."
!!   CALL get_energy3d(g_1,energy3d,5) 
!!   CALL sum3d_real(energy3d,energy_tot(5))
!!   !get hyper_collisions
!!   !IF(verbose) WRITE(*,*) "Before 6."
!!   CALL get_energy3d(g_1,energy3d,6) 
!!   CALL sum3d_real(energy3d,energy_tot(6))
!!   !get artificial
!!   !IF(verbose) WRITE(*,*) "Before 7."
!!   CALL get_energy3d(g_1,energy3d,7) 
!!   CALL sum3d_real(energy3d,energy_tot(7))
!!   !get nonlinearity
!!   !IF(verbose) WRITE(*,*) "Before 8."
!!   CALL get_energy3d(g_1,energy3d,8) 
!!   CALL sum3d_real(energy3d,energy_tot(8))
!!   !!!!!!!hyp_conv!!!!!!
!!   CALL get_energy3d(g_1,energy3d,9) 
!!   CALL sum3d_real(energy3d,energy_tot(9))
!!   !!!!!!hyp_conv!!!!!!
!!   !Actual dE/dt
!!   IF(time-time_last.gt.0.0) THEN
!!       energy_tot(10)=( (energy_tot(1)+energy_tot(2)) - energy_last)/(time-time_last)
!!   ELSE
!!       energy_tot(10)=0.0
!!   END IF
!!   energy_last=energy_tot(1)+energy_tot(2)
!!   time_last=time
!!   !IF(verbose) WRITE(*,*) "Done getting energy terms, output energy_tot"
!!   IF(mype==0) WRITE(en_handle,'(11es16.8)') time,energy_tot
!!   IF(mype==0) flush(en_handle)
!! 
!! END SUBROUTINE get_energy_terms
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              get_energy3d                                 !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_energy3d(g_in,energy3d,which_term)
!! 
!! IMPLICIT NONE
!! 
!!   REAL, INTENT(out) :: energy3d(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   REAL :: energy_temp1(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!   !COMPLEX :: energy_temp2(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   INTEGER, INTENT(in) :: which_term
!!   COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) :: g_int, g_out
!!   COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2) :: cfgamma1, cfgamma2
!!   COMPLEX, DIMENSION(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,ls1:ls2) :: g_mu0
!!   REAL, DIMENSION(0:nkx0-1,0:nky0-1,lh1:lh2) :: pref
!!   INTEGER :: i,j,k,ierr
!!   REAL, DIMENSION(0:nkx0-1,0:nky0-1) :: kperp
!!   INTEGER, DIMENSION(1) :: p_kperp
!!   COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!   INTEGER :: n_bcast
!!   INTEGER :: rhs_lin_version_bak
!!  
!!   kperp = sqrt(2*kperp2)  
!! 
!!   if(np_spec.gt.1) STOP "Need to implement get_energy3d for np_spec gt 1."
!! !!!!!!!!!Entropy!!!!!!!!!!!!
!! !!!!!!!!!Entropy!!!!!!!!!!!!
!! !!!!!!!!!Entropy!!!!!!!!!!!!
!! 
!!   IF (.not.mu_integrated) THEN
!!    IF(hankel) THEN
!!      CALL get_phi(g_in)
!!      CALL get_cfgamma_hk(g_in,phi,cfgamma1,cfgamma2)
!!      pref = 1.0
!! 
!!      IF(which_term==1) THEN
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = cfgamma1*g_in/2.
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1=0.0
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       energy3d=0.5*(energy3d)
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!!      ELSE IF(which_term==2) THEN
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = cfgamma2*g_in/2.
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1=0.0
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       energy3d=0.5*(energy3d)
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!!      ELSE IF(which_term==99) THEN
!!       CALL get_phi(g_in)
!!        DO k=0,nkz0-1
!!          energy3d(:,:,k)=REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!        END DO
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!! 
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!!     ELSE IF(which_term==3) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       CALL get_rhs_lin(g_in,phi,rhs,0) 
!!       IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(g_in,phi,rhs)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       !Modify
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!!     ELSE IF(which_term==4) THEN
!!       n_bcast=2/lv0  
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,6) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       !IF(lv1.le.2.and.lv2.ge.2) THEN
!!         g_int = cmplx(0.0,0.0)
!!         g_mu0 = cmplx(0.0,0.0)
!!         g_int = (cfgamma1+cfgamma2)*rhs
!!         CALL integral_hk(g_int,pref,g_mu0)
!!       !  IF(n_bcast.ne.mype) STOP "Error in get_energy3d!!!"
!!       !  energy3d(:,:,:)=REAL(SUM(g_mu0(:,:,:,:,0),4)) 
!!         energy_temp1(:,:,:)=REAL(SUM(g_mu0(:,:,:,:,0),4)) 
!!       !END IF
!!       !CALL MPI_BCAST(energy3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!!     ELSE IF(which_term==5) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,1) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!!     ELSE IF(which_term==6) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,2) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!!     ELSE IF(which_term==7) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,8)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!!     ELSE IF(which_term==8) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       IF(nonlinear.and..not.linear_nlbox) THEN
!!         CALL get_rhs_nl(g_in,phi,rhs)
!!       ELSE
!!         rhs=cmplx(0.0,0.0)
!!       END IF
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       !energy3d=REAL(energy_temp2)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!!     ELSE IF(which_term==9) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_conv
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_hk(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!!       !CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_conv
!! 
!!     ELSE
!!       STOP "Error in get_energy3d!"
!!     END IF
!! 
!!    ELSE !hankel on
!! 
!!      CALL get_phi(g_in)
!!      CALL get_cfgamma(g_in,phi,cfgamma1,cfgamma2)
!!      pref = 1.0
!! 
!!      IF(which_term==1) THEN
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = cfgamma1*g_in/2.
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1=0.0
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       energy3d=0.5*(energy3d)
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!!      ELSE IF(which_term==2) THEN
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = cfgamma2*g_in/2.
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1=0.0
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       energy3d=0.5*(energy3d)
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!!      ELSE IF(which_term==99) THEN
!!       CALL get_phi(g_in)
!!        DO k=0,nkz0-1
!!          energy3d(:,:,k)=REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!        END DO
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!! 
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!!     ELSE IF(which_term==3) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       CALL get_rhs_lin(g_in,phi,rhs,0) 
!!       IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(g_in,phi,rhs)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       !Modify
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!!     ELSE IF(which_term==4) THEN
!!       n_bcast=2/lv0  
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,6) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       !IF(lv1.le.2.and.lv2.ge.2) THEN
!!         g_int = cmplx(0.0,0.0)
!!         g_mu0 = cmplx(0.0,0.0)
!!         g_int = (cfgamma1+cfgamma2)*rhs
!!         CALL integral_v(g_int,pref,g_mu0)
!!       !  IF(n_bcast.ne.mype) STOP "Error in get_energy3d!!!"
!!       !  energy3d(:,:,:)=REAL(SUM(g_mu0(:,:,:,:,0),4)) 
!!         energy_temp1(:,:,:)=REAL(SUM(g_mu0(:,:,:,:,0),4)) 
!!       !END IF
!!       !CALL MPI_BCAST(energy3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!!     ELSE IF(which_term==5) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,1) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!!     ELSE IF(which_term==6) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,2) 
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!!     ELSE IF(which_term==7) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,8)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!! 
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!!     ELSE IF(which_term==8) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       IF(nonlinear.and..not.linear_nlbox) THEN
!!         CALL get_rhs_nl(g_in,phi,rhs)
!!       ELSE
!!         rhs=cmplx(0.0,0.0)
!!       END IF
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       !energy3d=REAL(energy_temp2)
!!       DEALLOCATE(rhs)
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!!     ELSE IF(which_term==9) THEN
!!       ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       rhs=cmplx(0.0,0.0)
!!       rhs_lin_version_bak=rhs_lin_version 
!!       rhs_lin_version=1
!!       CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!       rhs_lin_version=rhs_lin_version_bak
!!       g_int = cmplx(0.0,0.0)
!!       g_mu0 = cmplx(0.0,0.0)
!!       g_int = (cfgamma1+cfgamma2)*rhs
!!       CALL integral_v(g_int,pref,g_mu0)
!!       energy_temp1(:,:,:)=REAL(sum(g_mu0(:,:,:,:,0),4))
!!       CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!       DEALLOCATE(rhs)
!!     ELSE
!!       STOP "Error in get_energy3d!"
!!     END IF
!!  
!!   END IF !Hankel
!! 
!!  ELSE !mu_integrated
!!   IF(which_term==1) THEN
!!     energy_temp1(:,:,:)=REAL(pi**0.5*sum(conjg(g_in(:,:,:,:,0,0))*g_in(:,:,:,:,0,0),4))
!!     !CALL MPI_ALLREDUCE(energy_temp1,energy_temp2,nkx0*nky0*nkz0 &
!!     !  ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     energy3d=0.5*(energy3d)
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!
!!   ELSE IF(which_term==2) THEN
!!     CALL get_phi(g_in)
!!     energy_temp1=0.0
!!     DO k=0,nkz0-1
!!       !energy3d(:,:,k)=0.5*pi**(0.25)*phi_denom(:,:)*REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!       IF(mype_herm==0) energy_temp1(:,:,k)=0.5*pi**(0.25)*REAL(J0a(:,:)*g_in(:,:,k,0,0,0)*conjg(phi(:,:,k)))   
!!     END DO
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!!   ELSE IF(which_term==99) THEN
!!     CALL get_phi(g_in)
!!     DO k=0,nkz0-1
!!       energy3d(:,:,k)=REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!     END DO
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!! 
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!!   ELSE IF(which_term==3) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     CALL get_rhs_lin(g_in,phi,rhs,0) 
!!     IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(g_in,phi,rhs)
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype_herm==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     !Modify
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!!   ELSE IF(which_term==4) THEN
!!     n_bcast=2/lv0  
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,6) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     IF(lv1.le.2.and.lv2.ge.2) THEN
!!       IF(n_bcast.ne.mype_herm) STOP "Error in get_energy3d!!!"
!!       energy3d(:,:,:)=sqrt(pi)*REAL(conjg(g_in(:,:,:,2,0,0))*rhs(:,:,:,2,0,0)) 
!!     END IF
!!     CALL MPI_BCAST(energy3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_HERM,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==5) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,1) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==6) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,2) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!!   ELSE IF(which_term==7) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,8)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype_herm==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     DEALLOCATE(rhs)
!! 
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!!   ELSE IF(which_term==8) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     IF(nonlinear.and..not.linear_nlbox) THEN
!!       CALL get_rhs_nl(g_in,phi,rhs)
!!     ELSE
!!       rhs=cmplx(0.0,0.0)
!!     END IF
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype_herm==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     !energy3d=REAL(energy_temp2)
!!     DEALLOCATE(rhs)
!! 
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!!   ELSE IF(which_term==9) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_conv
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype_herm==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_HERM,ierr)
!!     DEALLOCATE(rhs)
!!   ELSE
!!     STOP "Error in get_energy3d!"
!!   END IF
!! ENDIF !mu_integrated
!! END SUBROUTINE get_energy3d
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_energy3d_lowmem                         !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_energy3d_lowmem(g_in,energy3d,which_term)
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(out) :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   REAL :: energy_temp1(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   !COMPLEX :: energy_temp2(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   INTEGER, INTENT(in) :: which_term
!!   INTEGER :: k,ierr
!!   COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!   INTEGER :: n_bcast
!!   INTEGER :: rhs_lin_version_bak
!! 
!! !!!!!!!!!Entropy!!!!!!!!!!!!1
!! !!!!!!!!!Entropy!!!!!!!!!!!!1
!! !!!!!!!!!Entropy!!!!!!!!!!!!1
!!   IF(which_term==1) THEN
!!     energy_temp1(:,:,:)=REAL(pi**0.5*sum(conjg(g_in(:,:,:,:,0,0))*g_in(:,:,:,:,0,0),4))
!!     !CALL MPI_ALLREDUCE(energy_temp1,energy_temp2,nkx0*nky0*nkz0 &
!!     !  ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     energy3d=0.5*(energy3d)
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!1
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!1
!! !!!!!!!!!Electrostatic Energy!!!!!!!!!!!!1
!!   ELSE IF(which_term==2) THEN
!!     CALL get_phi(g_in)
!!     DO k=0,nkz0-1
!!       energy3d(:,:,k)=0.5*pi**(-0.25)*phi_denom(:,:)*REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!     END DO
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!!   ELSE IF(which_term==99) THEN
!!     CALL get_phi(g_in)
!!     DO k=0,nkz0-1
!!       energy3d(:,:,k)=REAL(conjg(phi(:,:,k))*phi(:,:,k))   
!!     END DO
!! !!!!!!!!!Temporary!!!!!!!!!!!!
!! 
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!!   ELSE IF(which_term==3) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     CALL get_rhs_lin(g_in,phi,rhs,0) 
!!     IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(g_in,phi,rhs)
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     !energy3d=REAL(energy_temp2)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!!   ELSE IF(which_term==4) THEN
!!     n_bcast=2/lv0  
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,6) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     IF(lv1.le.2.and.lv2.ge.2) THEN
!!       IF(n_bcast.ne.mype) STOP "Error in get_energy3d!!!"
!!       energy3d(:,:,:)=sqrt(pi)*REAL(conjg(g_in(:,:,:,2,0,0))*rhs(:,:,:,2,0,0)) 
!!     END IF
!!     CALL MPI_BCAST(energy3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==5) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,1) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==6) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,2) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!!   ELSE IF(which_term==7) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,8)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     DEALLOCATE(rhs)
!! 
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!!   ELSE IF(which_term==8) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     IF(nonlinear.and..not.linear_nlbox) THEN
!!       CALL get_rhs_nl(g_in,phi,rhs)
!!     ELSE
!!       rhs=cmplx(0.0,0.0)
!!     END IF
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     !energy3d=REAL(energy_temp2)
!!     DEALLOCATE(rhs)
!! 
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!! !!!!!!!!!hyp_conv!!!!!!!!!!!!
!!   ELSE IF(which_term==9) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_conv
!!     rhs_lin_version=rhs_lin_version_bak
!!     energy_temp1(:,:,:)=REAL(sqrt(pi)*sum(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0),4))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         energy_temp1(:,:,k)=energy_temp1(:,:,k)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL MPI_ALLREDUCE(energy_temp1,energy3d,nkx0*nky0*nkz0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!     DEALLOCATE(rhs)
!!   ELSE
!!     STOP "Error in get_energy3d!"
!!   END IF
!! 
!! END SUBROUTINE get_energy3d_lowmem
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               sum3d_real                                  !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE sum3d_real(vec3d,scal)
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(in) :: vec3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   REAL, INTENT(out) :: scal
!!   INTEGER :: i
!! 
!!   scal=sum(sum(vec3d(0,:,:),1),1)
!!   DO i=1,nkx0-1
!!     scal=scal+2.0*sum(sum(vec3d(i,:,:),1),1)
!!   END DO
!! 
!! END SUBROUTINE sum3d_real
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              get_energy_hermite                           !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_energy_hermite
!!   IMPLICIT NONE
!! 
!!   COMPLEX :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   REAL :: energy_herm(0:nv0-1)
!!   REAL :: energy_herm2(0:nv0-1)
!!   REAL :: energy_herm_loc(lv1:lv2)
!!   INTEGER :: l,ierr
!! 
!!     if(np_hank.gt.1) STOP "get_energy_hermite not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_energy_hermite not yet implemented for np_spec.gt.1"
!!     if(np_kz.gt.1) STOP "get_energy_hermite not yet implemented for np_kz.gt.1"
!! 
!!   DO l=lv1,lv2 
!!     !energy_herm_loc(l)=REAL(sum(sum(sum(conjg(g_1(:,:,:,l))*g_1(:,:,:,l),1),1),1))
!!     CALL sum3d_real(REAL(conjg(g_1(:,:,:,l,0,0))*g_1(:,:,:,l,0,0)),energy_herm_loc(l))
!!   END DO
!!   energy_herm=0.0
!!   energy_herm(lv1:lv2)=energy_herm_loc(lv1:lv2)
!! 
!!   CALL MPI_ALLREDUCE(energy_herm,energy_herm2,nv0 &
!!       ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!   
!!   !IF(mype==0) WRITE(*,*) '('//char_nv0//'es16.8)'
!!   IF(mype==0) WRITE(herm_handle,'('//char_nv0//'es16.8)') time,energy_herm2
!!   IF(mype==0) flush(herm_handle)
!! 
!! END SUBROUTINE get_energy_hermite
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                            initialize_diag_eshells                        !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE initialize_diag_eshells
!!   IMPLICIT NONE
!! 
!!   IF(.not.shells_initialized) CALL get_shell_bounds
!!   CALL get_io_number
!!   eshells_handle=io_number
!! 
!!   !IF(mype==0) OPEN(unit=eshells_handle,file=trim(diagdir)//'/eshells.dat',status='replace',&
!!   !     form='unformatted',access='stream')
!! 
!!   IF(checkpoint_read) THEN
!!     INQUIRE(file=trim(diagdir)//'/eshells.dat',exist=file_exists)
!!     IF(file_exists) THEN
!!       IF(mype==0) OPEN(unit=eshells_handle,file=trim(diagdir)//'/eshells.dat',&
!!                           status='old',form='unformatted',access='stream',position='append')
!!     ELSE
!!       IF(mype==0) OPEN(unit=eshells_handle,file=trim(diagdir)//'/eshells.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!     END IF
!!   ELSE
!!     IF(mype==0) OPEN(unit=eshells_handle,file=trim(diagdir)//'/eshells.dat',&
!!                           status='replace',form='unformatted',access='stream')
!!   END IF
!! 
!! END SUBROUTINE initialize_diag_eshells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              finalize_diag_eshells                        !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE finalize_diag_eshells
!! 
!!   CLOSE(eshells_handle)
!!   
!! END SUBROUTINE finalize_diag_eshells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               diag_eshells                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE diag_eshells
!!   IMPLICIT NONE
!! 
!!   REAL :: energy_shells(num_shells,0:nkz0-1,lv1:lv2)
!!   REAL :: flux_shells(num_shells,0:nkz0-1)
!!   !!!!!!!!!!!Testing!!!!!!!!!!!!!
!!   !REAL :: shells_temp1(num_shells,0:nkz0-1)
!!   !REAL :: shells_temp2(num_shells,0:nkz0-1)
!!   !INTEGER :: k,kshift
!!   !!!!!!!!!!!Testing!!!!!!!!!!!!!
!!  
!!   !IF(verbose) WRITE(*,*) "in diag_eshells,mype",mype
!!   IF(mype==0) WRITE(eshells_handle) time
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,1) !Total E
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,3) !Total RHS
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,5) !Coll
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,6) !Hcoll
!!   CALL get_energy_shells(g_1,energy_shells,7) !hyps
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,8) !NL
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,9) !PH 1
!!   CALL output_eshells(energy_shells)
!!   energy_shells=0.0
!!   CALL get_energy_shells(g_1,energy_shells,10) !PH 2
!!   CALL output_eshells(energy_shells)
!!   CALL get_flux_shells(g_1,flux_shells)
!!   IF(mype==0) WRITE(eshells_handle) flux_shells
!!   IF(mype==0) flush(eshells_handle)
!! 
!! 
!! !SUBROUTINE get_flux_shells(g_in,flux_shells)
!! 
!! END SUBROUTINE diag_eshells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                             get_energy_shells                             !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_energy_shells(g_in,energy_shells,which_term)
!!   !which_term=1==>total energy
!!   !which_term=3==>RHS
!!   !which_term=5==>collisions
!!   !which_term=6==>hyper_coll
!!   !which_term=7==>artificial
!!   !which_term=8==>nonlinearity
!!   !which_term=9==>phase mixing 1
!!   !which_term=10==>phase mixing 2
!! 
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(inout) :: energy_shells(num_shells,0:nkz0-1,lv1:lv2)
!!   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2,lh1:lh2,ls1:ls2)
!!   INTEGER, INTENT(in) :: which_term
!!   INTEGER :: k,ierr
!!   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: etemp
!!   COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!   INTEGER :: n_bcast
!!   INTEGER :: rhs_lin_version_bak
!!   !REAL :: temp1,temp2
!! 
!! 
!!     if(np_hank.gt.1) STOP "get_energy_shells not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_energy_shells not yet implemented for np_spec.gt.1"
!! 
!!   IF(which_term==2) THEN
!!     IF(mype==0) WRITE(*,*) "Error! which_term=2 not valid for get_energy_shells."
!!     STOP
!!   END IF
!!   IF(which_term==4) THEN
!!     IF(mype==0) WRITE(*,*) "Error! which_term=4 not valid for get_energy_shells."
!!     STOP
!!   END IF
!! 
!!   ALLOCATE(etemp(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2))
!! !!!!!!!!!Total Energy!!!!!!!!!!!!1
!! !!!!!!!!!Total Energy!!!!!!!!!!!!1
!! !!!!!!!!!Total Energy!!!!!!!!!!!!1
!!   IF(which_term==1) THEN
!!     etemp=0.5*REAL(pi**0.5*(conjg(g_in(:,:,:,:,0,0))*g_in(:,:,:,:,0,0)))
!!     CALL get_phi(g_in)
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+0.5*pi**(0.25)*REAL(J0a(:,:)*g_in(:,:,k,0,0,0)*conjg(phi(:,:,k)))   
!!       END DO
!!     END IF
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     !temp1=sum(sum(sum(etemp(0,:,:,:),1),1),1)  
!!     !temp1=temp1+2.0*sum(sum(sum(sum(etemp(1:nkx0-1,:,:,:),1),1),1),1)
!!     !temp2=sum(sum(sum(energy_shells,1),1),1)
!!     !WRITE(*,*) "mype,etemp,eshells",mype,temp1,temp2
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!! !!!!!!!!!Total RHS!!!!!!!!!!!!
!!   ELSE IF(which_term==3) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     CALL get_rhs_lin(g_in,phi,rhs,0) 
!!     IF(nonlinear.and..not.linear_nlbox) CALL get_rhs_nl(g_in,phi,rhs)
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!No need - no n dependenc
!! !!!See get_flux_shells 
!! !!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!! !!!!!!!!!Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==5) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,1) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!! !!!!!!!!!Hyper-Collisions!!!!!!!!!!!!
!!   ELSE IF(which_term==6) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,2) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!! !!!!!!!!!Artificial!!!!!!!!!!!!
!!   ELSE IF(which_term==7) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,8)   !hyp_x,hyp_y,hyp_z,hyp_zonal
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!! 
!!     CALL get_rhs_lin(g_in,phi,rhs,10)   !hyp_conv
!!     etemp(:,:,:,:)=etemp(:,:,:,:)+REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!! 
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!!     rhs_lin_version=rhs_lin_version_bak
!! 
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!! !!!!!!!!!Nonlinearity!!!!!!!!!!!!
!!   ELSE IF(which_term==8) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     rhs=cmplx(0.0,0.0)
!!     CALL get_phi(g_in)
!!     IF(nonlinear.and..not.linear_nlbox) THEN
!!       CALL get_rhs_nl(g_in,phi,rhs)
!!     ELSE
!!       rhs=cmplx(0.0,0.0)
!!     END IF
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0))
!!     IF(mype==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Phase Mixing 1!!!!!!!!!!!!
!! !!!!!!!!!Phase Mixing 1!!!!!!!!!!!!
!! !!!!!!!!!Phase Mixing 1!!!!!!!!!!!!
!!   ELSE IF(which_term==9) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,3)   
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     CALL get_rhs_lin(g_in,phi,rhs,5)   
!!     etemp(:,:,:,:)=etemp(:,:,:,:)+REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     rhs_lin_version=rhs_lin_version_bak
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!! !!!!!!!!!Phase Mixing 2!!!!!!!!!!!!
!! !!!!!!!!!Phase Mixing 2!!!!!!!!!!!!
!! !!!!!!!!!Phase Mixing 2!!!!!!!!!!!!
!!   ELSE IF(which_term==10) THEN
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,9)   
!!     rhs_lin_version=rhs_lin_version_bak
!!     etemp(:,:,:,:)=REAL(sqrt(pi)*(conjg(g_in(:,:,:,:,0,0))*rhs(:,:,:,:,0,0)))
!!     IF(mype_herm==0) THEN
!!       DO k=0,nkz0-1
!!         etemp(:,:,k,0)=etemp(:,:,k,0)+REAL(pi**(0.25)*J0a(:,:)*conjg(phi(:,:,k))*rhs(:,:,k,0,0,0))
!!       END DO
!!     END IF
!!     CALL put_in_shells_real(etemp,energy_shells)
!!     DEALLOCATE(rhs)
!!   ELSE
!!     STOP "Error in get_energy3d!"
!!   END IF
!! 
!!   DEALLOCATE(etemp)
!! 
!! END SUBROUTINE get_energy_shells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               get_flux_shells                             !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE get_flux_shells(g_in,flux_shells)
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(out) :: flux_shells(num_shells,0:nkz0-1)
!!   REAL :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!   INTEGER :: k,ierr,i,j,s
!!   !REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: etemp
!!   COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!   INTEGER :: n_bcast
!!   INTEGER :: rhs_lin_version_bak
!! 
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!! !!!!!!!!!Drive Term!!!!!!!!!!!!
!!     n_bcast=2/lv0  
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     CALL get_phi(g_in)
!!     rhs_lin_version_bak=rhs_lin_version 
!!     rhs_lin_version=1
!!     CALL get_rhs_lin(g_in,phi,rhs,6) 
!!     rhs_lin_version=rhs_lin_version_bak
!!     IF(lv1.le.2.and.lv2.ge.2) THEN
!!       IF(n_bcast.ne.mype) STOP "Error in get_energy3d!!!"
!!       energy3d(:,:,:)=sqrt(pi)*REAL(conjg(g_in(:,:,:,2,0,0))*rhs(:,:,:,2,0,0)) 
!!     END IF
!!     CALL MPI_BCAST(energy3d,nkx0*nky0*nkz0,MPI_DOUBLE_PRECISION,n_bcast,MPI_COMM_WORLD,ierr) 
!!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!!     flux_shells=0.0
!!     DO i=1,nkx0-1
!!       DO j=0,nky0-1
!!         CALL get_shell_index(kxgrid(i),kygrid(j),s)
!!         flux_shells(s,:)=flux_shells(s,:)+2.0*energy3d(i,j,:) 
!!       END DO
!!     END DO
!! 
!!     i=0
!!     DO j=0,nky0-1
!!       CALL get_shell_index(kxgrid(i),kygrid(j),s)
!!       flux_shells(s,:)=flux_shells(s,:)+energy3d(i,j,:) 
!!     END DO
!! 
!!     DEALLOCATE(rhs)
!! 
!! END SUBROUTINE get_flux_shells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              put_in_shells_real                           !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE put_in_shells_real(e_in,e_shells)
!!   IMPLICIT NONE
!!   REAL, INTENT(inout) :: e_shells(num_shells,0:nkz0-1,lv1:lv2)
!!   REAL, INTENT(in) :: e_in(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2)
!!   INTEGER :: i,j,s
!! 
!!   DO i=1,nkx0-1
!!     DO j=0,nky0-1
!!       CALL get_shell_index(kxgrid(i),kygrid(j),s)
!!       e_shells(s,:,:)=e_shells(s,:,:)+2.0*e_in(i,j,:,:) 
!!     END DO
!!   END DO
!! 
!!   i=0
!!   DO j=0,nky0-1
!!     CALL get_shell_index(kxgrid(i),kygrid(j),s)
!!     e_shells(s,:,:)=e_shells(s,:,:)+e_in(i,j,:,:) 
!!   END DO
!! 
!! END SUBROUTINE put_in_shells_real
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               output_eshells                              !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE output_eshells(eshells)
!!   IMPLICIT NONE
!! 
!!   REAL, INTENT(out) :: eshells(num_shells,0:nkz0-1,lv1:lv2)
!!   REAL  :: eout(num_shells,0:nkz0-1,lv1:lv2)
!!   INTEGER :: l,p
!!   INTEGER :: send_proc,recv_proc,ierr
!!   INTEGER :: stat(MPI_STATUS_SIZE)
!! 
!!     IF(verbose) WRITE(*,*) "In output_eshells,mype",mype
!!     IF(mype==0) THEN
!!       WRITE(eshells_handle) eshells
!!     END IF 
!! 
!!     DO p=1,np_herm-1
!! 
!!      send_proc=p
!!      recv_proc=0
!!      IF(verbose) WRITE(*,*) "eshells:p,mype",p,mype
!! 
!!      IF(mype==send_proc) CALL MPI_Send(eshells(1,0,lv1), num_shells*nkz0*lv0, &
!!                          MPI_DOUBLE_PRECISION, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!      IF(mype==recv_proc) CALL MPI_Recv(eout(1,0,0), num_shells*nkz0*lv0, &
!!                          MPI_DOUBLE_PRECISION, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!      IF(mype==0) THEN
!!       WRITE(eshells_handle) eout
!!      END IF
!! 
!!     END DO
!! 
!! END SUBROUTINE output_eshells
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              test_shell_filter                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE test_shell_filter(phi_in)
!!   IMPLICIT NONE
!! 
!!   COMPLEX, INTENT(in) :: phi_in(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!   INTEGER :: i,j,k
!!   REAL :: kperp_min
!!   REAL :: kperp_max
!!   REAL :: kz_min
!!   REAL :: kz_max
!!   REAL :: kperp
!! 
!!   IF(spatial2d) STOP "test_shell_filter not yet implemented for spatial2d."
!! 
!!   kperp_min=10000.0
!!   kperp_max=0.0
!!   kz_min=10000.0
!!   kz_max=0.0
!!   DO i=0,nkx0-1
!!     DO j=0,nky0-1
!!       DO k=0,nkz0-1
!!          kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!         IF((abs(phi_in(i,j,k)).gt.1.0e-15)) THEN
!!           IF(kperp.lt.kperp_min) kperp_min=kperp
!!           IF(kperp.gt.kperp_max) kperp_max=kperp
!!           kz_min=kzgrid(k)
!!           kz_max=kzgrid(k)
!!         END IF
!!       END DO
!!     END DO
!!   END DO
!! 
!!   IF(mype==0) WRITE(*,*) 'kperp_min',kperp_min
!!   IF(mype==0) WRITE(*,*) 'kperp_max',kperp_max
!!   IF(mype==0) WRITE(*,*) 'kz_min',kz_min
!!   IF(mype==0) WRITE(*,*) 'kz_max',kz_max
!! 
!! END SUBROUTINE test_shell_filter
!! 
!! 
!! 
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                                nlt_shell_test1                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE nlt_shell_test1
!!   
!!     INTEGER :: i,j,k
!!     INTEGER :: ip,jp,kp
!!     !INTEGER :: shell_index
!!     !INTEGER :: shell_index_p
!!     REAL :: kxp,kyp,kzp
!!     INTEGER :: k0_nlt,kg,ierr
!!     COMPLEX :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs2
!!     COMPLEX :: phi_f(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     REAL :: sbounds(num_shells+1)
!!     REAL :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     REAL :: NLT_kxn0,NLT_kx0
!!     REAL :: NLT_temp
!!     REAL :: norm1,norm2
!!     INTEGER :: i0,j0,k0
!!     LOGICAL :: conjg_p
!! 
!!     IF(spatial2d) STOP "nlt_shell_test1 not yet implemented for spatial2d."
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     ALLOCATE(rhs2(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!! 
!!     CALL get_real(g_1)
!!     CALL get_phi(g_1)
!!     DO i=0,nkx0-1
!!       DO j=0,nky0-1
!!         DO k=0,nkz0-1
!!           g_f=cmplx(0.0,0.0)
!!           g_f(i,j,k,:,:,:)=g_1(i,j,k,:,:,:)
!!           IF(i.eq.0) THEN !Enforce reality constraint
!!             CALL get_indices_from_ks(0.0,-kygrid(j),-kzgrid(k),i0,j0,k0,conjg_p)
!!             IF(conjg_p) THEN
!!               g_f(i0,j0,k0,:,:,:)=conjg(g_1(i0,j0,k0,:,:,:))
!!             ELSE
!!               g_f(i0,j0,k0,:,:,:)=g_1(i0,j0,k0,:,:,:)
!!             END IF
!!           END IF 
!!           rhs=cmplx(0.0,0.0)
!!           CALL get_rhs_nl(g_f,phi,rhs)
!!           CALL get_rhs_kp(i,j,k,rhs2) 
!!           norm1=sum(sum(sum(sum(abs(rhs(:,:,:,:,0,0)),1),1),1),1)
!!           norm2=sum(sum(sum(sum(abs(rhs2(:,:,:,:,0,0)),1),1),1),1)
!!           !IF(i.eq.0.and.mype==0) THEN
!!           !IF(i.ne.0.and.mype==0) THEN
!!           IF(mype==0) THEN
!!             WRITE(*,*) "rhs1",i,j,k,norm1
!!             WRITE(*,*) "rhs2",i,j,k,norm2
!!             IF(norm1.ne.0.0) WRITE(*,*) "RHS diff:",i,j,k,sum(sum(sum(sum(abs(rhs(:,:,:,:,0,0)&
!!                 -rhs2(:,:,:,:,0,0)),1),1),1),1)/sum(sum(sum(sum(abs(rhs(:,:,:,:,0,0)),1),1),1),1)
!!             !DO i0=0,nkx0-1
!!             !  DO j0=0,nky0-1
!!             !    DO k0=0,nkz0-1
!!             !      norm1=sum(abs(rhs(i0,j0,k0,:)-rhs2(i0,j0,k0,:)))
!!             !      norm2=sum(abs(rhs(i0,j0,k0,:)))
!!             !      IF(norm2.gt.1.0e-15) norm1=norm1/norm2
!!             !      IF(norm2.gt.1.0e-15) THEN
!!             !        WRITE(*,*) "rhs index diff:",i0,j0,k0,norm1
!!             !      ELSE
!!             !        WRITE(*,*) "No diff!!!!!:",i0,j0,k0,norm1
!!             !      END IF
!!             !    END DO
!!             !  END DO
!!             !END DO
!! 
!!           END IF
!!           
!!         END DO
!!       END DO
!!     END DO
!!    
!!     DEALLOCATE(rhs)
!!     DEALLOCATE(rhs2)
!! 
!!   END SUBROUTINE nlt_shell_test1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_kp                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Given ikx_in,iky_in,ikz_in for the k primes, calculate the resulting nl rhs
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE get_rhs_kp(ikx_in,iky_in,ikz_in,rhs)
!!     IMPLICIT NONE
!! 
!!     INTEGER, INTENT(in) :: ikx_in,iky_in,ikz_in
!!     REAL :: kx0,ky0,kz0
!!     REAL :: kxp,kyp,kzp
!!     REAL :: kxpp,kypp,kzpp
!!     INTEGER :: ikxp,ikyp,ikzp
!!     INTEGER :: ikxpp,ikypp,ikzpp
!!     !REAL, INTENT(out) :: nltk(-ikx_minmax:ikx_minmax,iky_minmax:iky_minmax,&
!!     !                          -ikz_minmax:ikz_minmax)
!!     COMPLEX :: phibk0
!!     COMPLEX :: phibkp
!!     COMPLEX :: phibkpp
!!     COMPLEX :: gk0(lv1:lv2)
!!     COMPLEX :: gkp(lv1:lv2)
!!     COMPLEX :: gkpp(lv1:lv2)
!!     INTEGER :: i,j,k
!!     LOGICAL :: conjg_pp,conjg_p
!!     !INTEGER :: ikx_minmax
!!     !INTEGER :: iky_minmax
!!     !INTEGER :: ikz_minmax
!!     REAL :: ckkp
!! 
!!     COMPLEX, INTENT(out) :: rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!! 
!!     if(np_hank.gt.1) STOP "get_rhs_kp not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_rhs_kp not yet implemented for np_spec.gt.1"
!!     if(spatial2d) STOP "get_rhs_kp not yet implemented for spatial2d"
!! 
!!     rhs=cmplx(0.0,0.0)
!! 
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!!     !!!!!!!!!!!!!!!Fix conjg's for reality constraint
!! 
!!     DO i=0,nkx0-1
!!       DO j=0,nky0-1
!!         DO k=0,nkz0-1
!! 
!! 
!!           ikxp=ikx_in
!!           ikyp=iky_in
!!           ikzp=ikz_in
!! 
!!           kx0=kxgrid(i)
!!           ky0=kygrid(j)
!!           kz0=kzgrid(k)
!!           kxp=kxgrid(ikxp)
!!           kyp=kygrid(ikyp)
!!           kzp=kzgrid(ikzp)
!!           kxpp=kx0-kxp
!!           kypp=ky0-kyp
!!           kzpp=kz0-kzp
!! 
!!           IF((abs(kxpp).le.kxmax).and.(abs(kypp).le.kymax).and.(abs(kzpp).le.kzmax)) THEN
!! 
!!             ckkp=kxp*ky0-kx0*kyp
!! 
!!             CALL get_indices_from_ks(kxpp,kypp,kzpp,ikxpp,ikypp,ikzpp,conjg_pp)
!!             phibkp=(J0a(ikxp,ikyp)*phi(ikxp,ikyp,ikzp))
!!             gkp(:)=(g_1(ikxp,ikyp,ikzp,:,0,0))
!! 
!!             IF(conjg_pp) THEN
!!               phibkpp=conjg(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!               gkpp(:)=conjg(g_1(ikxpp,ikypp,ikzpp,:,0,0))
!!             ELSE
!!               phibkpp=(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!               gkpp(:)=(g_1(ikxpp,ikypp,ikzpp,:,0,0))
!!             END IF
!! 
!!             rhs(i,j,k,:,0,0)=rhs(i,j,k,:,0,0)-ckkp*phibkpp*gkp(:)
!! 
!! 
!! 
!!           END IF
!! 
!!           !Now get other part from reality constraint
!!           ikxp=ikx_in
!!           ikyp=iky_in
!!           ikzp=ikz_in
!!           kxp=kxgrid(ikxp)
!!           kyp=kygrid(ikyp)
!!           kzp=kzgrid(ikzp)
!!           !get new indices 
!!           CALL get_indices_from_ks(-kxp,-kyp,-kzp,ikxp,ikyp,ikzp,conjg_p)
!!           !get new kp's
!!           kxp=kxgrid(ikxp)
!!           kyp=kygrid(ikyp)
!!           kzp=kzgrid(ikzp)
!!           !get new kpp's
!!           kxpp=kx0-kxp
!!           kypp=ky0-kyp
!!           kzpp=kz0-kzp
!! 
!!           IF((abs(kxpp).le.(kxmax+0.001)).and.(abs(kypp).le.(kymax+0.001)).and.(abs(kzpp).le.(kzmax+0.001))) THEN
!! 
!!             ckkp=kxp*ky0-kx0*kyp
!! 
!!             CALL get_indices_from_ks(kxpp,kypp,kzpp,ikxpp,ikypp,ikzpp,conjg_pp)
!! 
!!             IF(conjg_p) THEN
!!               !phibkp=conjg(J0a(ikxp,ikyp)*phi(ikxp,ikyp,ikzp))
!!               gkp(:)=conjg(g_1(ikxp,ikyp,ikzp,:,0,0))
!!             ELSE
!!               !phibkp=(J0a(ikxp,ikyp)*phi(ikxp,ikyp,ikzp))
!!               gkp(:)=(g_1(ikxp,ikyp,ikzp,:,0,0))
!!             END IF
!! 
!!             IF(conjg_pp) THEN
!!               phibkpp=conjg(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!               !gkpp(:)=conjg(g_1(ikxpp,ikypp,ikzpp,:))
!!             ELSE
!!               phibkpp=(J0a(ikxpp,ikypp)*phi(ikxpp,ikypp,ikzpp))
!!               !gkpp(:)=(g_1(ikxpp,ikypp,ikzpp,:))
!!             END IF
!! 
!!             rhs(i,j,k,:,0,0)=rhs(i,j,k,:,0,0)-ckkp*phibkpp*gkp(:)
!! 
!!           END IF
!! 
!! 
!! 
!!         END DO
!!       END DO
!!     END DO
!! 
!!   END SUBROUTINE get_rhs_kp
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                                get_nlt_shells2                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE get_nlt_shells2
!!     
!!     INTEGER :: i,j,k
!!     INTEGER :: ip,jp,kp
!!     INTEGER :: shell_index
!!     INTEGER :: shell_index_p
!!     REAL :: kxp,kyp,kzp
!!     INTEGER :: k0_nlt
!!     !!!!!Testing
!!     INTEGER :: a,b
!!     !!!!!Testing
!! 
!!     !IF(mype==0) OPEN(unit=555,file=trim(diagdir)//'/nlttest.dat',status='replace')
!!     IF(mype==0.and.verbose) WRITE(*,*) "In get_nlt_shells2."
!!     IF(mype==0) WRITE(nlt_status,*) time
!!     NLT=0.0
!!     DO i=0,nkx0-1
!!     !DO i=1,nkx0-1
!!       IF(mype==0) WRITE(nlt_status,*) "nlt",i
!!       IF(mype==0) flush(nlt_status)
!!       IF(mype==0) WRITE(*,*) "nlt",i
!!       !DO j=0,nky0-1
!!       DO j=0,nky0-1
!!         IF(j.ne.nky0/2) THEN
!!         !DO k=-ikz_minmax,ikz_minmax
!!         DO k=-ikz_minmax,ikz_minmax
!! 
!!           IF(k.ge.0) THEN
!!             k0_nlt=k
!!           ELSE
!!             k0_nlt=k+nkz0
!!           END IF
!!           CALL get_shell_index(kxgrid(i),kygrid(j),shell_index)
!!           CALL get_nlt_k(i,j,k0_nlt)
!!           !!!!Testing
!!           !IF(k==2) WRITE(555,*) i,j
!!           !DO a=-ikx_minmax,ikx_minmax
!!           !  DO b=-iky_minmax,iky_minmax
!!           !    IF(k==2) WRITE(555,*) a,b,nltk(a,b,2)
!!           !  END DO
!!           !END DO
!!           !IF(k==2) WRITE(555,*) ""
!!           !IF(k==2) WRITE(555,*) ""
!!           !!!!Testing
!!           !IF(mype==0.and.MOD(istep_nlt_full,itime)==0) WRITE(nlt_handle) i
!!           !IF(mype==0.and.MOD(istep_nlt_full,itime)==0) WRITE(nlt_handle) j
!!           !IF(mype==0.and.MOD(istep_nlt_full,itime)==0) WRITE(nlt_handle) k0_nlt
!!           !IF(mype==0.and.MOD(istep_nlt_full,itime)==0) THEN
!!           !    WRITE(nlt_handle) nltk
!!           !END IF
!!           !DO ip=-ikx_minmax,ikx_minmax
!!           !  DO jp=-iky_minmax,iky_minmax
!!           !    DO kp=-ikz_minmax,ikz_minmax
!!           DO ip=-ikx_minmax,ikx_minmax
!!             !IF(ip.ne.0) THEN
!!             DO jp=-iky_minmax,iky_minmax
!!               DO kp=-ikz_minmax,ikz_minmax
!!                 CALL get_shell_index(ip*kxmin,jp*kymin,shell_index_p) 
!!                 IF(mype==0) THEN
!!                   IF(i.eq.0) THEN
!!                     NLT(shell_index,shell_index_p,k,kp)=NLT(shell_index,shell_index_p,k,kp)+nltk(ip,jp,kp)
!!                   ELSE
!!                     NLT(shell_index,shell_index_p,k,kp)=NLT(shell_index,shell_index_p,k,kp)+nltk(ip,jp,kp)
!!                     NLT(shell_index,shell_index_p,-k,-kp)=NLT(shell_index,shell_index_p,-k,-kp)+nltk(ip,jp,kp)
!!                   END IF
!!                   !IF(k==kp.and.abs(nltk(ip,jp,kp)).gt.1.0e-15) WRITE(*,*) "nltk (k==kp)",nltk(ip,jp,kp)
!!                 END IF
!!               END DO
!!             END DO
!!             !END IF
!!           END DO
!! 
!!         END DO
!!         END IF
!!       END DO
!!     END DO
!! 
!!     !IF(mype==0) CLOSE(555)
!! 
!!     IF(mype==0) WRITE(shell_handle) time
!!     IF(mype==0) WRITE(shell_handle) NLT
!! 
!!     IF(mype==0.and.verbose) WRITE(*,*) "Done get_nlt_shells2."
!! 
!!   END SUBROUTINE get_nlt_shells2
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                             initialize_diag_gk                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE initialize_diag_gk
!!     IMPLICIT NONE
!!     INTEGER :: i
!!     CHARACTER(len=2) :: ch_kx,ch_ky!,ch_kz
!!     CHARACTER(len=100) :: file_name
!! 
!! 
!!     gk_nkyout=nky0/4-1
!!     ALLOCATE(gk_handle(gk_nkyout))
!!     ALLOCATE(gk_indices(gk_nkyout))
!! 
!!     DO i=1,gk_nkyout
!!       CALL get_io_number
!!       gk_handle(i)=io_number
!!     END DO
!! 
!!     DO i=1,gk_nkyout
!!       gk_indices(i)=2*i
!!     END DO
!! 
!!     DO i=1,gk_nkyout
!! 
!!       WRITE(ch_kx,'(i2.2)') gk_kx_index
!!       WRITE(ch_ky,'(i2.2)') gk_indices(i)
!!       !WRITE(ch_kz,'(i2.2)') gk_indices(3,i)
!!       file_name='g_kx'//ch_kx//'ky'//ch_ky//'.dat'
!! 
!!       IF((gk_kx_index.le.nkx0-1).and.(gk_indices(i).lt.(nky0)/2-1)) THEN
!! 
!!         INQUIRE(file=trim(diagdir)//'/'//trim(file_name),exist=file_exists)
!!         IF(checkpoint_read.and.file_exists) THEN
!!           IF(mype==0) OPEN(unit=gk_handle(i),file=trim(diagdir)//'/'//trim(file_name),&
!!                               status='old',form='unformatted',access='stream',position='append')
!!         ELSE
!!           IF(mype==0) OPEN(unit=gk_handle(i),file=trim(diagdir)//'/'//trim(file_name),&
!!                                 status='replace',form='unformatted',access='stream')
!!         END IF
!! 
!!       END IF !gk_indices test
!! 
!!     END DO
!!     !IF(mype==0) WRITE(*,*) "gk_handle",gk_handle
!! 
!!   END SUBROUTINE initialize_diag_gk
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                             initialize_diag_gknl                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE initialize_diag_gknl
!!     IMPLICIT NONE
!!     INTEGER :: i
!!     CHARACTER(len=2) :: ch_kx,ch_ky!,ch_kz
!!     CHARACTER(len=100) :: file_name
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !CHARACTER(len=100) :: temp_file_name
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !!!!!!!BENCHMARKING!!!!!!!!
!! 
!! 
!!     gknl_nkyout=nky0/4-1
!!     ALLOCATE(gknl_handle(gknl_nkyout))
!!     ALLOCATE(gknl_indices(gknl_nkyout))
!! 
!!     DO i=1,gknl_nkyout
!!       CALL get_io_number
!!       gknl_handle(i)=io_number
!!     END DO
!! 
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !ALLOCATE(gknl_temphandle(gknl_nkyout))
!!     !DO i=1,gknl_nkyout
!!     !  CALL get_io_number
!!     !  gknl_temphandle(i)=io_number
!!     !END DO
!!     !!!!!!!BENCHMARKING!!!!!!!!
!!     !!!!!!!BENCHMARKING!!!!!!!!
!! 
!!     DO i=1,gknl_nkyout
!!       gknl_indices(i)=2*i
!!     END DO
!! 
!!     DO i=1,gknl_nkyout
!! 
!!       WRITE(ch_kx,'(i2.2)') gk_kx_index
!!       WRITE(ch_ky,'(i2.2)') gknl_indices(i)
!!       !WRITE(ch_kz,'(i2.2)') gk_indices(3,i)
!!       file_name='nl_kx'//ch_kx//'ky'//ch_ky//'.dat'
!!       !!!!!!!BENCHMARKING!!!!!!!!
!!       !!!!!!!BENCHMARKING!!!!!!!!
!!       !temp_file_name='temp_kx'//ch_kx//'ky'//ch_ky//'.dat'
!!       !!!!!!!BENCHMARKING!!!!!!!!
!!       !!!!!!!BENCHMARKING!!!!!!!!
!! 
!!       IF((gk_kx_index.le.nkx0-1).and.(gknl_indices(i).lt.(nky0)/2-1)) THEN
!! 
!!         INQUIRE(file=trim(diagdir)//'/'//trim(file_name),exist=file_exists)
!!         IF(checkpoint_read.and.file_exists) THEN
!!           IF(mype==0) OPEN(unit=gknl_handle(i),file=trim(diagdir)//'/'//trim(file_name),&
!!                               status='old',form='unformatted',access='stream',position='append')
!!         ELSE
!!           IF(mype==0) OPEN(unit=gknl_handle(i),file=trim(diagdir)//'/'//trim(file_name),&
!!                                 status='replace',form='unformatted',access='stream')
!!         END IF
!! 
!!         !!!!!!!BENCHMARKING!!!!!!!!
!!         !!!!!!!BENCHMARKING!!!!!!!!
!!         !INQUIRE(file=trim(diagdir)//'/'//trim(temp_file_name),exist=file_exists)
!!         !IF(checkpoint_read.and.file_exists) THEN
!!         !  IF(mype==0) OPEN(unit=gknl_temphandle(i),file=trim(diagdir)//'/'//trim(temp_file_name),&
!!         !                      status='old',form='unformatted',access='stream',position='append')
!!         !ELSE
!!         !  IF(mype==0) OPEN(unit=gknl_temphandle(i),file=trim(diagdir)//'/'//trim(temp_file_name),&
!!         !                        status='replace',form='unformatted',access='stream')
!!         !END IF
!!         !!!!!!!BENCHMARKING!!!!!!!!
!!         !!!!!!!BENCHMARKING!!!!!!!!
!! 
!! 
!!       END IF !gk_indices test
!! 
!!     END DO
!! 
!!     !IF(mype==0) WRITE(*,*) "gknl_handle",gknl_handle
!! 
!!   END SUBROUTINE initialize_diag_gknl
!! 
!! 
!! 
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               finalize_diag_gk                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE finalize_diag_gk
!!     IMPLICIT NONE
!!     INTEGER :: i
!! 
!!     DO i=1,gk_nkyout
!!       CLOSE(gk_handle(i))
!!     END DO
!!     DEALLOCATE(gk_handle)
!!     DEALLOCATE(gk_indices)
!! 
!!   END SUBROUTINE finalize_diag_gk
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                               finalize_diag_gknl                            !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE finalize_diag_gknl
!!     IMPLICIT NONE
!!     INTEGER :: i
!! 
!!     DO i=1,gknl_nkyout
!!       CLOSE(gknl_handle(i))
!!     END DO
!!     DEALLOCATE(gknl_handle)
!!     DEALLOCATE(gknl_indices)
!! 
!!   END SUBROUTINE finalize_diag_gknl
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                                   diag_gk                                 !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE diag_gk
!!     IMPLICIT NONE
!! 
!!     INTEGER :: i,ikx,iky
!!     COMPLEX :: g_send(0:nkz0-1,lv1:lv2)
!!     COMPLEX :: g_out(0:nkz0-1,lv1:lv2)
!!     INTEGER :: send_proc,recv_proc,p,ierr
!!     INTEGER :: stat(MPI_STATUS_SIZE)
!! 
!! 
!!     if(np_hank.gt.1) STOP "diag_gk not yet implemented for np_hank.gt.1"
!!     if(np_kz.gt.1) STOP "diag_gk not yet implemented for np_kz.gt.1"
!!     if(np_spec.gt.1) STOP "diag_gk not yet implemented for np_spec.gt.1"
!! 
!!     DO i=1,gk_nkyout
!! 
!!       IF((gk_kx_index.le.nkx0-1).and.(gk_indices(i).lt.(nky0)/2-1)) THEN
!! 
!!         IF(mype==0) WRITE(gk_handle(i)) time
!! 
!!         ikx=gk_kx_index
!!         iky=gk_indices(i)
!! 
!!         IF(mype==0) WRITE(gk_handle(i)) g_1(ikx,iky,:,:,0,0)
!! 
!!         DO p=1,np_herm-1
!! 
!!           send_proc=p
!!           recv_proc=0
!!           IF(verbose) WRITE(*,*) "p,mype",p,mype
!!           g_send=g_1(ikx,iky,:,:,0,0)
!!   
!!           !Modify
!!           IF(mype==send_proc) CALL MPI_Send(g_send, lv0*nkz0, &
!!                               MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!           !Modify
!!           IF(mype==recv_proc) CALL MPI_Recv(g_out, lv0*nkz0, &
!!                               MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!           IF(mype==0) WRITE(gk_handle(i)) g_out
!! 
!!         END DO
!! 
!!         IF(mype==0) FLUSH(gk_handle(i)) 
!! 
!!       END IF
!! 
!!     END DO
!! 
!!   END SUBROUTINE diag_gk
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                                   diag_gknl                               !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE diag_gknl
!!     IMPLICIT NONE
!! 
!!     INTEGER :: i,ikx,iky
!!     COMPLEX :: g_send(0:nkz0-1,lv1:lv2)
!!     COMPLEX :: g_out(0:nkz0-1,lv1:lv2)
!!     INTEGER :: send_proc,recv_proc,p,ierr
!!     INTEGER :: stat(MPI_STATUS_SIZE)
!!     COMPLEX :: rhsnl(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2)
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !    REAL :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!! !    REAL :: rhs3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!! !    REAL :: rhsnl3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!! !    REAL :: energy3d_temp(0:nkx0-1,0:nky0-1,0:nkz0-1)
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! 
!!     if(np_hank.gt.1) STOP "diag_gknl not yet implemented for np_hank.gt.1"
!!     if(np_kz.gt.1) STOP "diag_gknl not yet implemented for np_kz.gt.1"
!!     if(np_spec.gt.1) STOP "diag_gknl not yet implemented for np_spec.gt.1"
!! 
!!     CALL get_phi(g_1) 
!!     rhsnl=0.0
!!     CALL get_rhs_nl(g_1,phi,rhsnl)
!! 
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !    !Total energy
!! !    CALL get_energy3d(g_1,energy3d,1) 
!! !    CALL get_energy3d(g_1,energy3d_temp,2) 
!! !    energy3d=energy3d+energy3d_temp
!! !    CALL get_energy3d(g_1,rhsnl3d,8) 
!! !    CALL get_energy3d(g_1,rhs3d,3) 
!! !    DO i=1,gknl_nkyout
!! !
!! !      IF((gk_kx_index.le.nkx0-1).and.(gknl_indices(i).lt.(nky0)/2-1)) THEN
!! !
!! !        IF(mype==0) WRITE(gknl_temphandle(i)) time
!! !        !IF(mype==0) WRITE(*,*) "in diag_gknl",time
!! !
!! !        ikx=gk_kx_index
!! !        iky=gknl_indices(i)
!! !
!! !        IF(mype==0) WRITE(gknl_temphandle(i)) energy3d(ikx,iky,:)
!! !        IF(mype==0) WRITE(gknl_temphandle(i)) rhs3d(ikx,iky,:)
!! !        IF(mype==0) WRITE(gknl_temphandle(i)) rhsnl3d(ikx,iky,:)
!! !
!! !        IF(mype==0) FLUSH(gknl_temphandle(i)) 
!! !
!! !      END IF
!! !
!! !    END DO
!! !
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! !!!!!!!!!BENCHMARKING!!!!!!!!!!!!!!!!
!! 
!! 
!!     DO i=1,gknl_nkyout
!! 
!!       IF((gk_kx_index.le.nkx0-1).and.(gknl_indices(i).lt.(nky0)/2-1)) THEN
!! 
!!         IF(mype==0) WRITE(gknl_handle(i)) time
!!         !IF(mype==0) WRITE(*,*) "in diag_gknl",time
!! 
!!         ikx=gk_kx_index
!!         iky=gknl_indices(i)
!! 
!!         IF(mype==0) WRITE(gknl_handle(i)) rhsnl(ikx,iky,:,:)
!! 
!!         DO p=1,np_herm-1
!! 
!!           send_proc=p
!!           recv_proc=0
!!           IF(verbose) WRITE(*,*) "p,mype",p,mype
!!           g_send=rhsnl(ikx,iky,:,:)
!!   
!!           !Modify
!!           IF(mype==send_proc) CALL MPI_Send(g_send, lv0*nkz0, &
!!                               MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!!           !Modify
!!           IF(mype==recv_proc) CALL MPI_Recv(g_out, lv0*nkz0, &
!!                               MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!           IF(mype==0) WRITE(gknl_handle(i)) g_out
!! 
!!         END DO
!!         IF(mype==0) FLUSH(gknl_handle(i)) 
!! 
!!       END IF
!! 
!!     END DO
!! 
!!   END SUBROUTINE diag_gknl
!! 
!!   SUBROUTINE get_nlt_triple
!! 
!!     IMPLICIT NONE
!!     
!!     INTEGER :: ierr
!!     !COMPLEX :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!     !COMPLEX :: phi_f(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: rhs
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:)  :: g_K
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:)  :: g_P
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:)  :: g_Q
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:)  :: phi_P
!!     COMPLEX, ALLOCATABLE, DIMENSION(:,:,:)  :: phi_Q
!!     REAL :: sbounds(num_shells+1)
!!     REAL :: energy3d(0:nkx0-1,0:nky0-1,0:nkz0-1)
!!     REAL :: NLT_kxn0,NLT_kx0
!!     REAL :: NLT_temp
!!     LOGICAL :: nlt_filled(num_shells,num_shells,0:ikz_minmax,0:ikz_minmax)
!!     INTEGER :: nlt_out_ns(11)
!!     INTEGER :: nout_index
!!     INTEGER :: ik,ip,iq
!! 
!!     if(np_hank.gt.1) STOP "get_nlt_shells not yet implemented for np_hank.gt.1"
!!     if(np_spec.gt.1) STOP "get_nlt_shells not yet implemented for np_spec.gt.1"
!!     if(np_kz.gt.1) STOP "get_nlt_shells not yet implemented for np_kz.gt.1"
!! 
!!     IF(num_shells-1.ne.size(shell_bounds)) THEN
!!         WRITE(*,*) "Error in get_nlt_shells!"
!!         STOP
!!     END IF
!!     sbounds(2:num_shells)=shell_bounds
!!     sbounds(1)=0.0
!!     sbounds(num_shells+1)=100000.0
!!     IF(mype==0.and.verbose) WRITE(*,*) "In get_nlt_shells."
!!     !IF(mype==0) WRITE(nlt_handle) time
!!     NLT3=0.0
!!     CALL get_phi(g_1)
!! 
!!     ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     ALLOCATE(g_K(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     ALLOCATE(g_P(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!     ALLOCATE(phi_Q(0:nkx0-1,0:nky0-1,lkz1:lkz2))
!!     DO ik=1,num_shells
!!       CALL g_shell_filter3(sbounds(ik),sbounds(ik+1),g_K) 
!!       DO ip=1,num_shells
!!         CALL g_shell_filter3(sbounds(ip),sbounds(ip+1),g_P) 
!!         !CALL phi_shell_filter3(sbounds(ip),sbounds(ip+1),phi_P) 
!!         DO iq=1,num_shells
!! 
!!           !CALL g_shell_filter3(sbounds(iq),sbounds(iq+1),g_Q) 
!!           CALL phi_shell_filter3(sbounds(iq),sbounds(iq+1),phi,phi_Q) 
!!           rhs=cmplx(0.0,0.0)
!!           CALL get_rhs_nl(g_P,phi_Q,rhs)
!!           NLT_kxn0=REAL(sqrt(pi)*sum(sum(sum(sum(conjg(g_K(1:nkx0-1,:,:,:,0,0))*rhs(1:nkx0-1,:,:,:,0,0),1),1),1),1))
!!           CALL MPI_ALLREDUCE(NLT_kxn0,NLT_temp,1 &
!!             ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!! 
!!           NLT3(ik,ip,iq)=NLT3(ik,ip,iq)+2.0*NLT_temp
!! 
!!           NLT_kx0=REAL(sqrt(pi)*sum(sum(sum(conjg(g_K(0,:,:,:,0,0))*rhs(0,:,:,:,0,0),1),1),1))
!!           CALL MPI_ALLREDUCE(NLT_kx0,NLT_temp,1 &
!!              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!! 
!!           NLT3(ik,ip,iq)=NLT3(ik,ip,iq)+NLT_temp
!! 
!!         END DO
!!       END DO
!!     END DO
!! 
!!     DEALLOCATE(g_P)
!!     DEALLOCATE(phi_Q)
!! 
!!     IF(nlt_symmetrize) THEN
!!       ALLOCATE(g_Q(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!!       ALLOCATE(phi_P(0:nkx0-1,0:nky0-1,lkz1:lkz2))
!! 
!!       DO ik=1,num_shells
!!         CALL g_shell_filter3(sbounds(ik),sbounds(ik+1),g_K) 
!!         DO ip=1,num_shells
!!           !CALL g_shell_filter3(sbounds(ip),sbounds(ip+1),g_P) 
!!           CALL phi_shell_filter3(sbounds(ip),sbounds(ip+1),phi,phi_P) 
!!           DO iq=1,num_shells
!!   
!!             CALL g_shell_filter3(sbounds(iq),sbounds(iq+1),g_Q) 
!!             !CALL phi_shell_filter3(sbounds(iq),sbounds(iq+1),phi_Q) 
!!             rhs=cmplx(0.0,0.0)
!!             CALL get_rhs_nl(g_Q,phi_P,rhs)
!!             NLT_kxn0=REAL(sqrt(pi)*sum(sum(sum(sum(conjg(g_K(1:nkx0-1,:,:,:,0,0))*rhs(1:nkx0-1,:,:,:,0,0),1),1),1),1))
!!             CALL MPI_ALLREDUCE(NLT_kxn0,NLT_temp,1 &
!!               ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!   
!!             NLT3(ik,ip,iq)=NLT3(ik,ip,iq)+2.0*NLT_temp
!!   
!!             NLT_kx0=REAL(sqrt(pi)*sum(sum(sum(conjg(g_K(0,:,:,:,0,0))*rhs(0,:,:,:,0,0),1),1),1))
!!             CALL MPI_ALLREDUCE(NLT_kx0,NLT_temp,1 &
!!                ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!   
!!             NLT3(ik,ip,iq)=NLT3(ik,ip,iq)+NLT_temp
!!   
!!           END DO
!!         END DO
!!       END DO
!! 
!!       DEALLOCATE(g_Q)
!!       DEALLOCATE(phi_P)
!!     END IF
!! 
!!     DEALLOCATE(g_K)
!!     DEALLOCATE(rhs)
!! 
!! !    ALLOCATE(rhs(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2))
!! 
!! !        !First filter the distribution function so that it only has elements in
!! !        !the appropriate kperp shell along with +-kp.  Note that we must USE
!! !        !+-kp since the RHS calculation implicitly uses the -kp modes for the
!! !        !kx<0 modes (i.e., the reality constraint gives kx<0,kp=-kp modes)
!! !        CALL g_shell_filter(sbounds(ip),sbounds(ip+1),kp*kzmin,g_f)
!! !        CALL get_phi(g_1)
!! !        rhs=cmplx(0.0,0.0)
!! !        CALL get_rhs_nl(g_f,phi,rhs)
!! 
!!               !CALL g_shell_filter(sbounds(i),sbounds(i+1),k*kzmin,g_f)
!! 
!!               !NLT_kxn0=REAL(sqrt(pi)*sum(sum(sum(sum(conjg(g_f(1:nkx0-1,:,:,:,0,0))*rhs(1:nkx0-1,:,:,:,0,0),1),1),1),1))
!!               !CALL MPI_ALLREDUCE(NLT_kxn0,NLT_temp,1 &
!!               !  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!! 
!!               !NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+2.0*NLT_temp
!! 
!! !              NLT_kx0=REAL(sqrt(pi)*sum(sum(sum(conjg(g_f(0,:,:,:,0,0))*rhs(0,:,:,:,0,0),1),1),1))
!! !              CALL MPI_ALLREDUCE(NLT_kx0,NLT_temp,1 &
!! !                ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!! !              NLT(i,ip,k,kp)=NLT(i,ip,k,kp)+NLT_temp
!! !              nlt_filled(i,ip,k,kp)=.true.
!! !              NLT(ip,i,kp,k)=-NLT(i,ip,k,kp)
!! !              nlt_filled(ip,i,kp,k)=.true.
!!              
!! 
!! 
!!   END SUBROUTINE get_nlt_triple
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                             g_shell_filter0                                !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE g_shell_filter3(kperp_min,kperp_max,g_f)
!!     IMPLICIT NONE
!! 
!!     REAL, INTENT(in) :: kperp_min,kperp_max
!!     COMPLEX, INTENT(out) :: g_f(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!!     INTEGER :: k
!!     INTEGER :: i,j
!!     REAL :: kperp
!! 
!!     g_f=cmplx(0.0,0.0)
!!     DO i=0,nkx0-1
!!       DO j=0,nky0-1
!!         IF(j.ne.nky0/2) THEN
!!           kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!           IF((kperp.ge.kperp_min).and.(kperp.lt.kperp_max)) THEN
!!             g_f(i,j,:,:,:,:)=g_1(i,j,:,:,:,:)
!!           END IF
!!         END IF
!!       END DO
!!     END DO
!! 
!!   END SUBROUTINE g_shell_filter3
!! 
!! 
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!                              phi_shell_filter3                             !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE phi_shell_filter3(kperp_min,kperp_max,phi_in,phi_f)
!!     IMPLICIT NONE
!! 
!!     REAL, INTENT(in) :: kperp_min,kperp_max
!!     COMPLEX, INTENT(in) :: phi_in(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!     COMPLEX, INTENT(out) :: phi_f(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!!     INTEGER :: k
!!     INTEGER :: i,j
!!     REAL :: kperp
!!     
!!     phi_f=cmplx(0.0,0.0)
!!     DO i=0,nkx0-1
!!       DO j=0,nky0-1
!!         IF(j.ne.nky0/2) THEN
!!           kperp=sqrt(kxgrid(i)**2+kygrid(j)**2)
!!           IF((kperp.ge.kperp_min.and.kperp.lt.kperp_max)) THEN
!!             phi_f(i,j,:)=phi_in(i,j,:)
!!           END IF
!!         END IF
!!       END DO
!!     END DO
!! 
!!  END SUBROUTINE phi_shell_filter3



function hmhdhmtn(v_in,b_in) result(ham)

implicit none
complex :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
real :: ham

!! Computes the Hall MHD Hamiltonian 8 pi^3 (sum(v_k^2 + b_k^2) + 2b_z0)/2 

ham = sum(abs(v_in)**2 + abs(b_in)**2) + 2*real(b_in(0,0,0,2))
ham = ham * (4*(pi**3))

end function hmhdhmtn

function resvischange(arr) result(res)

!! Returns the rate of change of the above Hamiltonian 
!! From including a resistivity term in the induction equation
!! with dimensionless resistivity eta/(B0/(e ne c)) = 1
!! Can also be used to find change due to viscosity nu/(vA^2/omega_i) = 1
!! with arr taken as v instead of b


implicit none
complex :: arr(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
integer :: i,j,k,ind
real :: k2
real :: res

res = 0.
do i = 0,nkx0-1
    do j = 0,nky0-1
        do k = lkz1,lkz2
            k2 = kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2
            do ind = 0,2
                res = res + (k2**hyp)*abs(arr(i,j,k,ind))**2
            end do
        end do 
    end do 
end do

res = res*(8 * pi**3)

end function resvischange

subroutine initialize_debug

INTEGER :: ionums(9)
! CHARACTER(len=15), DIMENSION(9) :: fnames
CHARACTER(len=15), DIMENSION(9) :: fnames = [CHARACTER(LEN=15) :: &
     'dissb_out.dat','dissv_out.dat', 'bdv_out.dat',&
     'vdb_out.dat','bdcb_out.dat','cbdb_out.dat',&
     'vdv_out.dat','bdb_out.dat','db2_out.dat']

INTEGER :: q

ionums = [dbio,dvio,bdvio,vdbio,bdcbio,cbdbio,vdvio,bdbio,db2io]

DO q = 1,9
   IF(checkpoint_read) THEN
      INQUIRE(file=trim(diagdir)//'/'//trim(fnames(q)),exist=file_exists)      
      IF(file_exists) THEN

         OPEN(unit=ionums(q),file=trim(diagdir)//'/'//trim(fnames(q)),form='unformatted', status='REPLACE',access='stream')
      ELSE
         OPEN(unit=ionums(q),file=trim(diagdir)//'/'//trim(fnames(q)),form='unformatted', status='REPLACE',access='stream')
      END IF
   ELSE
      OPEN(unit=ionums(q),file=trim(diagdir)//'/'//trim(fnames(q)),form='unformatted', status='REPLACE',access='stream')
   END IF

ENDDO

if (verbose) print *, 'Debug files opened'

end subroutine initialize_debug

subroutine finalize_debug

INTEGER :: ionums(9)
INTEGER :: q

ionums = [dbio,dvio,bdvio,vdbio,bdcbio,cbdbio,vdvio,bdbio,db2io]

DO q = 1,9
CLOSE(ionums(q))
ENDDO 

if (verbose) print *, 'Debug files closed' 
end subroutine finalize_debug

function vec_potential(b) result(A)

implicit none
complex :: b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: A(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: bg(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
integer :: i,j,k,ind


bg = b

! A = - (curl B)/k^2
A = cmplx(0.0,0.0)

do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
        if ((i**2+j**2+k**2).ne.0) then 
          if (ind.eq.0) A(i,j,k,ind) = kygrid(j) * bg(i,j,k,2) - kzgrid(k) * bg(i,j,k,1) 
          if (ind.eq.1) A(i,j,k,ind) = kzgrid(k) * bg(i,j,k,0) - kxgrid(i) * bg(i,j,k,2)
          if (ind.eq.2) A(i,j,k,ind) = kxgrid(i) * bg(i,j,k,1) - kygrid(j) * bg(i,j,k,0)
          A(i,j,k,ind) = cmplx(-1.0,0.0) * A(i,j,k,ind) / (kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2)
        endif
      enddo
    enddo
   enddo
enddo

if ((verbose).and.(itime.lt.200)) write(*,*) 'A',A

end function vec_potential

function mag_helicity(b0) result(maghel)

implicit none
complex :: b0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
real :: maghel
complex :: A0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
integer :: i,j,k,ind
real :: intx

maghel = 0.0
A0 = vec_potential(b0)

do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
        maghel = maghel + real(conjg(A0(i,j,k,ind))*b0(i,j,k,ind))
      enddo
    enddo
!    if (i.ne.0) then
!      intx = (2 * sin(kxgrid(i)*pi) - 2 * kxgrid(i) * pi * cos(kxgrid(i)*pi)) &
!        / (kxgrid(i)**2)
!      maghel = maghel - aimag(b0(i,0,0,1)) * intx / (2*pi)
      ! if (verbose) print *, 'Guide\n Field\n Contribution',intx,aimag(b0(i,0,0,1)),- aimag(b0(i,0,0,1)) * intx / (2*pi)
    ! endif
  enddo
enddo

maghel = maghel + real(A0(0,0,0,2))
maghel = real((2*pi)**3 * maghel)

end function mag_helicity

function vorticity(v) result(w)

implicit none
complex :: v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: w(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
integer :: ind,i,j,k

w = cmplx(0.0,0.0)

do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
          if (ind.eq.0) w(i,j,k,ind) = kygrid(j) * v(i,j,k,2) - kzgrid(k) * v(i,j,k,1)
          if (ind.eq.1) w(i,j,k,ind) = kzgrid(k) * v(i,j,k,0) - kxgrid(i) * v(i,j,k,2)
          if (ind.eq.2) w(i,j,k,ind) = kxgrid(i) * v(i,j,k,1) - kygrid(j) * v(i,j,k,0)
      enddo
    enddo
   enddo
enddo
w = cmplx(0.0,1.0) * w

if ((verbose).and.(itime.lt.200)) write(*,*) 'w',w

end function vorticity

function cross_helicity(v0,b0) result(crosshel)

implicit none
complex :: b0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: v0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
real :: crosshel
complex :: A0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: w0(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
integer :: i,j,k,ind

crosshel = mag_helicity(b0) / (8 * pi**3)
w0 = vorticity(v0)

do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
        crosshel = crosshel + conjg(b0(i,j,k,ind)) * v0(i,j,k,ind) + conjg(v0(i,j,k,ind)) * b0(i,j,k,ind)
        crosshel = crosshel + conjg(v0(i,j,k,ind)) * w0(i,j,k,ind)
      enddo
    enddo
  enddo
enddo

crosshel = crosshel + 3*v0(0,0,0,2)
crosshel = real((2*pi)**3 * crosshel)

end function cross_helicity

function en_spec(v_in,b_in) result(esubk)

implicit none
complex :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
complex :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
real :: esubk(0:nkx0-1,0:nky0-1,lkz1:lkz2)
integer :: i,j,k,ind

esubk = 0.0
do ind = 0,2
  do i = 0,nkx0-1
    do j = 0,nky0-1
      do k = lkz1,lkz2
        esubk(i,j,k) = esubk(i,j,k) + abs(v_in(i,j,k,ind))**2 + abs(b_in(i,j,k,ind))**2
      enddo
    enddo
  enddo
enddo

esubk(0,0,0) = esubk(0,0,0) + 2*real(b_in(0,0,0,2))
esubk = (4 * pi**3) * esubk

end function en_spec

END MODULE diagnostics


