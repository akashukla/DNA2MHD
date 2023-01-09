!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                               cc_comm.f90                                 !!
!!                                                                           !!
!!  communications                                                           !!
!!  -- init_comm                                                             !!
!!  -- comm                                                                  !!
!!  -- finalize_mpi                                                          !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                               communications                              !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE communications
  USE par_mod
  USE mpi
  IMPLICIT NONE

  PUBLIC :: init_comm, comm, finalize_mpi
  PUBLIC :: my_complex_sum_hank


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 init_comm                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_comm
  USE mpi
  USE par_mod  
  IMPLICIT NONE

  INTEGER :: ierr

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,n_mpi_procs,ierr) 
  
  IF(mype==0) WRITE(*,*) "number of processors", n_mpi_procs

END SUBROUTINE init_comm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                   comm                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE comm 
  USE mpi
  USE par_mod
  IMPLICIT NONE
  integer :: arr_dims(4),ierr
  logical :: if_periodic(4),sub_grid(4)
  integer :: mype_cart(4)
  
  !mpi_comm_cart_4d=0
  arr_dims=(/np_kz,np_herm,np_hank,np_spec/)
  if_periodic=(/.false.,.false.,.false.,.false./)
  CALL MPI_CART_CREATE(MPI_COMM_WORLD,4,arr_dims,if_periodic,.false.,mpi_comm_cart_4d,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_CART_COORDS(mpi_comm_cart_4d,mype,4,mype_cart,ierr)
  mype_kz=mype_cart(1)
  mype_herm=mype_cart(2)
  mype_hank=mype_cart(3)
  mype_spec=mype_cart(4)

!  if(mype==0) write(*,*) "Creating kz subcommunicator."
  sub_grid=(/.true.,.false.,.false.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_kz,ierr)

!  if(mype==0) write(*,*) "Creating y subcommunicator."
  sub_grid=(/.false.,.true.,.false.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_herm,ierr)
  IF(verbose.and.(mype.eq.0)) write(*,*) "mpi_comm_herm",mpi_comm_herm

!  if(mype==0) write(*,*) "Creating hank subcommunicator."
  sub_grid=(/.false.,.false.,.true.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_hank,ierr)

!  if(mype==0) write(*,*) "Creating v subcommunicator."
  sub_grid=(/.false.,.false.,.false.,.true./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_spec,ierr)

END SUBROUTINE comm

SUBROUTINE comm_temp 
  USE mpi
  USE par_mod
  IMPLICIT NONE
  integer :: arr_dims(4),ierr
  logical :: if_periodic(4),sub_grid(4)
  integer :: mype_cart(4)
  
  !mpi_comm_cart_4d=0
  arr_dims=(/np_kz,np_herm,np_hank,np_spec/)
  if_periodic=(/.false.,.false.,.false.,.false./)
  !CALL MPI_CART_CREATE(MPI_COMM_WORLD,4,arr_dims,if_periodic,.false.,mpi_comm_cart_4d,ierr)
  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !CALL MPI_CART_COORDS(mpi_comm_cart_4d,mype,4,mype_cart,ierr)
  mype_kz=0
  mype_herm=0
  mype_hank=0
  mype_spec=0

!  if(mype==0) write(*,*) "Creating kz subcommunicator."
  sub_grid=(/.true.,.false.,.false.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_kz,ierr)

!  if(mype==0) write(*,*) "Creating y subcommunicator."
  sub_grid=(/.false.,.true.,.false.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_herm,ierr)
  IF(verbose.and.(mype.eq.0)) write(*,*) "mpi_comm_herm",mpi_comm_herm

!  if(mype==0) write(*,*) "Creating z subcommunicator."
  sub_grid=(/.false.,.false.,.true.,.false./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_hank,ierr)

!  if(mype==0) write(*,*) "Creating v subcommunicator."
  sub_grid=(/.false.,.false.,.false.,.true./)
  call MPI_CART_SUB(mpi_comm_cart_4d,sub_grid,mpi_comm_spec,ierr)

END SUBROUTINE comm_temp


SUBROUTINE my_complex_sum_hank(localsum,len)
    INTEGER, INTENT(IN) :: len
    COMPLEX, INTENT(INOUT) :: localsum(len)
    COMPLEX :: totalsum(len)
    INTEGER :: ierr

    call MPI_ALLREDUCE(localsum, totalsum, len,&
    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_HANK, ierr)

    localsum  = totalsum

END SUBROUTINE my_complex_sum_hank


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               finalize_mpi                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE finalize_mpi
  USE mpi
  
  INTEGER :: ierr

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF (verbose.and.(mype.eq.0)) print *,'MPI Barrier called'
  CALL MPI_FINALIZE(ierr)
  IF (verbose.and.(mype.eq.0)) print *,'MPI Finalize called'

END SUBROUTINE finalize_mpi


END MODULE communications
