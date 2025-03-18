program mpioatdebug

  use, intrinsic :: iso_c_binding
  
  use mpi

  implicit none

  integer(4) :: ierr,mype,n_mpi_procs
  integer(4) :: threehandle,threesize
  character(len=100) :: threename="threewave_out.dat"
  integer(4) :: mype1 = 4,reader,mype2 = 20,mype3 = 17,wavemype(3),i
  real(8) :: time,lwm,lcm,rwm,rcm
  integer(8) :: offset

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,n_mpi_procs,ierr)


  wavemype = [mype1,mype2,mype1]
  
  DO reader = 1,5
     CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(threename),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,threehandle,ierr)
     CALL MPI_FILE_GET_SIZE(threehandle,threesize,ierr)

     offset = threesize
     if (mype.eq.0) then
        time = 0.06*(reader-1)
        CALL MPI_FILE_WRITE_AT(threehandle,offset,time,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
     endif
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     DO i = 1,3
        if (mype.eq.wavemype(i)) then
           offset = threesize+8+32*(i-1)
           
           lwm = dble(mype)*4
           lcm = dble(mype)
           rwm = dble(mype)
           rcm = dble(mype)
           
           CALL MPI_FILE_WRITE_AT(threehandle,offset,lwm,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
           CALL MPI_FILE_WRITE_AT(threehandle,offset+8_8,lcm,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
           CALL MPI_FILE_WRITE_AT(threehandle,offset+16_8,rwm,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
           CALL MPI_FILE_WRITE_AT(threehandle,offset+24_8,rcm,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
        
        endif

     ENDDO
     
     CALL MPI_FILE_CLOSE(threehandle,ierr)
  ENDDO
  
  CALL MPI_FINALIZE(ierr)

end program mpioatdebug



  
