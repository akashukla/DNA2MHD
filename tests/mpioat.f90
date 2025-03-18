program mpioat

  use, intrinsic :: iso_c_binding
  use mpi

  implicit none
  
  integer(4) :: ierr,mype,n_mpi_procs
  integer(8) :: offset
  integer(4) :: threehandle
  character(len=100) :: threename="threewave_out.dat"
  integer(4) :: threesize
  logical :: deletefile
  real(8) :: written
  integer(4) :: repeat

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,n_mpi_procs,ierr)

  inquire(file=trim(threename),exist=deletefile)
  if (deletefile) then
     CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(threename),MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,threehandle,ierr)
     CALL MPI_FILE_CLOSE(threehandle,ierr)
  endif

  DO repeat = 1,7
  
     CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(threename),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,threehandle,ierr)
     CALL MPI_FILE_GET_SIZE(threehandle,threesize,ierr)
     
     offset = threesize
     written = -99
     if (mype.eq.0) then
        written = 3.1415
        CALL MPI_FILE_WRITE_AT(threehandle,offset,written,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
     endif
     
     if (mype.eq.1) then
        offset = offset+8
        written = 1.732
        CALL MPI_FILE_WRITE_AT(threehandle,offset,written,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
        CALL MPI_FILE_WRITE_AT(threehandle,offset+8_8,written,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)     
     endif

     if (mype.eq.2) then
        offset = offset + 24
        written = 1.414
        CALL MPI_FILE_WRITE_AT(threehandle,offset,written,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)	
        CALL MPI_FILE_WRITE_AT(threehandle,offset+8_8,written,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
     endif
     
     CALL MPI_FILE_CLOSE(threehandle,ierr)
     
  ENDDO

  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_FINALIZE(ierr)
  

end program mpioat
