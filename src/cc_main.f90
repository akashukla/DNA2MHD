!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                                 cc_main.f90                               !!
!!                                                                           !!
!!  dna                                                                      !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           
!!                        Direct Numerical Analysis                          
!!                of fundamental gyrokinetic turbulence dynamics             
!!
!!  The DNA code solves the reduced gyrokinetic ecuations in slab geometry 
!!  using a pseudo-spectral method.
!!
!!  DNA creator:
!!    D.R. Hatch - davidhat@ipp.mpg.de
!!
!!  DNA developers (AKA Mutants) listed in alphabetical order:
!!    V. Bratanov - Vasil.Bratanov@ipp.mpg.de
!!    A.B. Navarro - alejandro.banon.navarro@ipp.mpg.de
!!    B. Teaca - bteaca@ipp.mpg.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                  dna                                      !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM dna
  !USE calculate_time_step, ONLY:  finalize_adapt_dt!, calc_dt_lapack
  !USE communications, ONLY: finalize_mpi, init_comm, comm
  !USE diagnostics, ONLY: finalize_diagnostics, nl_test, start_wallclock
  !USE eigen_iterative, ONLY: ev_iterative
  !USE flr_effects, ONLY: finalize_flr
  !USE hk_effects, ONLY: finalize_hk
  !USE mpi
  USE nonlinearity, ONLY: initialize_fourier  
  USE par_mod
  !USE performance, ONLY: end_clock, performance_compute, start_clock
  USE time_advance, ONLY: iv_solver, rk4_stability
  !USE triple_transfers, ONLY: get_triple_transfers
  IMPLICIT NONE
  
  CHARACTER(len=40) :: ic_temp
  REAL :: time_tot
  INTEGER :: run_type
  !! Starting the MPI communication enviromen and the wallclock
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !CALL init_comm
  !CALL start_wallclock
  
  !! Reading the input parameters 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL read_parameters

  IF (.not.mu_integrated) THEN
    flr_on = .false.
  END IF
  
  !if(calc_dt.and..not.dt_slepc) call calc_dt_lapack

  !CALL comm
  
  !! Initiating the current run of the simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL init_run
  
  !! Test of rk4--i.e. calculates stability region
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(test_rk4) CALL rk4_stability
  
  !! Runs various routines for performance tests
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !IF(performance_tests) CALL performance_compute

  !! Initial output
  !!!!!!!!!!!!!!!!!!!
  CALL output_parameters


  !! Selects what to do this run 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(comp_type=='IV') run_type=0
  IF(comp_type=='EV') run_type=1
  IF(test_nl)         run_type=2
  IF(comp_type=='TT') run_type=3
  
  SELECT CASE(run_type)
  
  !! Main loop
  !!!!!!!!!!!!!! 
  CASE(0)
    IF (verbose) WRITE(*,*) "Starting time loop.",mype
    !IF (performance_tests) CALL start_clock
    CALL iv_solver
    !IF (performance_tests) THEN 
    !  CALL end_clock(time_tot)
      IF(mype==0) THEN
        WRITE(*,*) "Total wallclock for time loop:",time_tot
        WRITE(*,*) "itime",itime
        OPEN(unit=999,file=trim(diagdir)//'/end_info.dat',status='unknown')
        WRITE(999,*) "Total wall clock:",time_tot
        WRITE(999,*) "itime",itime
        WRITE(999,*) "time",time
        CLOSE(999)
      END IF
    !END IF

  !! Eigenvalue calculation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CASE(1) 
!    IF(ev_slepc) THEN
!      !which_ev='jd'
!      !CALL ev_iterative
!    ELSE
!      !Lapack eigensolve--serial only
!      !Not used at this time--must check
!      !IF(calc_dt.and.n_mpi_procs==1) CALL calc_each_dt_lapack
!      !IF(kscan.and.n_mpi_procs==1) THEN
!      !  CALL scan_kspace
!      !ELSE
!      !  CALL eigensolve
!      !END IF
!      !CALL finalize_ev
!    END IF
!    !!
!    IF(mype==0) WRITE(*,*) "Eigenvalue calculation completed."
! 
!  !! Compares the results of the nonlinearity as calculated by 
!  !! pseudo-spectral method and by the convolution method.
!  !! The latter is slow and only for benchmarking purposes.
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CASE(2)
!    IF(mype==0) WRITE(*,*) "Testing nonlinearity"
!    ic_temp=init_cond                                                        
!    CALL initial_condition(ic_temp)
!    !IF(mype==0) WRITE(*,*) "initialize_fourier"
!    CALL initialize_fourier
!    !IF(mype==0) WRITE(*,*) "nl_test"
!    CALL nl_test     
!    
!  CASE(3) !! Triple shell diagnostics 
!    IF (mype==0) WRITE(*,*) "Triple Transfers"
!    CALL get_triple_transfers
!    
!  CASE DEFAULT
!    IF(mype==0) WRITE(*,*) "Wrong selection for run_type=",run_type
!    STOP
  END SELECT

 
  !! Finalizing the run
  !!!!!!!!!!!!!!!!!!!!!!!
  IF (verbose) WRITE(*,*) "Finalizing diagnostics.",mype
  !CALL finalize_diagnostics
  !CALL finalize_adapt_dt
  !CALL finalize_flr
  !CALL finalize_hk
  CALL finalize_arrays

  !! Shuts down the MPI enviroment 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !CALL finalize_mpi

END PROGRAM dna
