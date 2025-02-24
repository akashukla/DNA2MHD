program mpiop3d

  ! mpif90 -g -I$P3DFFT_INTEL/include -I$FFTW_ROOT/include -o mps mpiosimple.f90 -L$P3DFFT_INTEL/lib -L$FFTW_ROOT/lib -lp3dfft -lfftw3_mpi -lfftw3 
  use mpi
  use p3dfft
  use, intrinsic :: iso_c_binding

  implicit none

  integer(4) :: ierr,n_mpi_procs,mype
  integer(4) :: sizes(4),subsizes(4),starts(4),lsize ! From Gropp Using Advanced MPI
  integer(4) :: carray3
  integer(4) :: fsize,hsize

  character(len=100) :: fname = "/pscratch/sd/e/echansen/testing/mpiotest1.dat"
  
  integer(4) :: cstart(3),cend(3),csize(3)
  complex(8), allocatable :: cin(:,:,:,:),cin2(:,:,:,:),din(:,:,:,:),din2(:,:,:,:)
  integer(4) :: nx = 48, ny = 48, nz = 48
  integer(4) :: npx = 2,npy

  integer(kind=MPI_OFFSET_KIND) :: offset
  
  integer(4) :: i,j,k,l
  integer(4) :: write_handle,read_handle

  integer(4) :: outi1,outi2,outi3,outi4
  real(8) :: outr1,outr2,outr3
  real(8) :: diffm1,diff1,diffm2,diff2

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,n_mpi_procs,ierr)

  npy = n_mpi_procs/npx

  call p3dfft_setup([npx,npy],nx,ny,nz,MPI_COMM_WORLD)
  call p3dfft_get_dims(cstart,cend,csize,2)

  ALLOCATE(cin(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(cin2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(din(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))
  ALLOCATE(din2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2))

  DO i = cstart(1),cend(1)
     DO j  = cstart(2),cend(2)
        DO k = cstart(3),cend(3)
           DO l = 0,2
              cin(i,j,k,l) = cmplx(0.01*i,0.02*j)
              din(i,j,k,l) = cmplx(0.03*j,0.07*i)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  sizes = (/1+nx/2,ny,nz,3/)
  subsizes = [(1+cend(i)-cstart(i), integer :: i = 1,3),3]
  starts = [(cstart(i)-1, integer :: i=1,3),0]
  lsize = product(subsizes)
  
  print *, mype,lsize
  print *, mype,csize
  print *, mype,starts

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_type_create_subarray(4,sizes,subsizes,starts,MPI_ORDER_FORTRAN,&
       MPI_COMPLEX16,carray3,ierr)
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_type_commit(carray3,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_open(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,&
       MPI_INFO_NULL,write_handle,ierr)

  if (mype.eq.0) then
     call mpi_file_write(write_handle,1001,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,42.03_8,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,16,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,32,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,32,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,0.1_8,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
     call mpi_file_write(write_handle,0.000123_8,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
  endif

  offset = 40

  call mpi_file_set_view(write_handle,offset,&
       MPI_COMPLEX16,carray3,"native",MPI_INFO_NULL,ierr)

  if (mype.eq.0) print *, "View set",ierr
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_write_all(write_handle,cin,lsize,MPI_COMPLEX16,MPI_STATUS_IGNORE,ierr)

  offset = 40 + product(sizes)*16
  call mpi_file_set_view(write_handle,offset,MPI_COMPLEX16,carray3,"native",MPI_INFO_NULL,ierr)
  call mpi_file_write_all(write_handle,din,lsize,MPI_COMPLEX16,MPI_STATUS_IGNORE,ierr)
  
  if (mype.eq.0) print *, "Written",ierr
  
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_get_size(write_handle,fsize,ierr)
  call mpi_file_close(write_handle,ierr)

  if (mype.eq.0) print *, "Written size ",fsize

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_open(MPI_COMM_WORLD,trim(fname),&
       MPI_MODE_RDONLY+MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&
       read_handle,ierr)

  call mpi_file_get_size(read_handle,fsize,ierr)
  if (mype.eq.0) print *, "Read size ",fsize

  call mpi_file_read(read_handle,outi1,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outr1,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outi2,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outi3,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outi4,1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outr2,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
  call mpi_file_read(read_handle,outr3,1,MPI_REAL8,MPI_STATUS_IGNORE,ierr)

  offset = 40
  call mpi_file_set_view(read_handle,offset,&
       MPI_COMPLEX16,carray3,"native",MPI_INFO_NULL,ierr)
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_read_all(read_handle,cin2,lsize,MPI_COMPLEX16,MPI_STATUS_IGNORE,ierr)

  print *, "Read",ierr

  offset = 40 + product(sizes)*16
  call mpi_file_set_view(read_handle,offset,MPI_COMPLEX16,carray3,"native",MPI_INFO_NULL,ierr)
  call mpi_file_read_all(read_handle,din2,lsize,MPI_COMPLEX16,MPI_STATUS_IGNORE,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_file_close(read_handle,ierr)
  
  ! Check Headers
  
  print *, mype,"Outi1",outi1-1001
  print *, mype,"Outi2",outi2-16
  print *, mype,"Outi3",outi3-32
  print *, mype,"Outi4",outi4-32
  print *, mype,"Outr1",outr1-42.03_8
  print *, mype,"Outr2",outr2-0.1_8
  print *, mype,"Outr3",outr3-0.000123_8

  diffm1 = maxval(abs(cin-cin2))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(diffm1,diff1,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

  diffm2 = maxval(abs(din-din2))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(diffm2,diff2,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  print *, mype,"Max Difference Array 1",diffm1,maxloc(abs(cin-cin2))
  print *, mype,"Max Difference Array 2",diffm2,maxloc(abs(din-din2))

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_type_free(carray3,ierr)

  DEALLOCATE(cin2)
  print *, "cin2 deallocated"
  DEALLOCATE(cin)
  print *, "cin deallocated"
  
  call p3dfft_clean
  print *, "p3dfft deallocated"

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_finalize(ierr)

end program mpiop3d
