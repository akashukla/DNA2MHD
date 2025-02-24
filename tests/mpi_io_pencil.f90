program mpi_io_pencil
  use mpi
  use p3dfft

  ! Bulk of code courtesy of ChatGPT
  ! I edited to make complex16 which broke the code,
  ! Using complex(8) should work

  ! Build with
  ! mpif90 -g -I$P3DFFT_INTEL/include -I$FFTW_ROOT/include -o mpi_io_pencil mpi_io_pencil -L$P3DFFT_INTEL/lib -L$FFTW_ROOT/lib -lp3dfft -lfftw3_mpi -lfftw3
  ! Requires eight processors
  
  
    implicit none

    integer :: ierr, mype, nprocs
    integer :: dims(2), periods(2), coords(2), cart_comm
    integer :: nx, ny, nz   ! Global dimensions
    integer :: local_nx, local_ny, local_nz
    integer :: sizes(3), subsizes(3), starts(3)
    integer :: cstart(3),cend(3),csize(3)
    integer :: file_handle, subarray_type
    integer(kind=MPI_OFFSET_KIND) :: offset
    complex(8), allocatable :: data(:, :, :), read_data(:, :, :)

    ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

    ! Define 2D Cartesian grid (processors distributed in x & y)
    dims = (/ 0, 0 /)  ! Let MPI determine dimensions
    call MPI_Dims_create(nprocs, 2, dims, ierr)
    periods = (/ 0, 0 /)  ! No periodic boundary conditions
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., cart_comm, ierr)
    call MPI_Cart_coords(cart_comm, mype, 2, coords, ierr)

    ! Global array dimensions (example: nx x ny x nz grid)
    nx = 8
    ny = 8
    nz = 8  ! The last dimension is **not** distributed

    dims = [2,4]
    call p3dfft_setup(dims,nx,ny,nz,MPI_COMM_WORLD)
    call p3dfft_get_dims(cstart,cend,csize,2)    

    ! Define local dimensions (2D decomposition)
    local_nx = nx / dims(1)   ! Divide along x
    local_ny = ny / dims(2)   ! Divide along y
    local_nz = nz             ! The z-dimension remains the same

    ! Allocate local arrays
    allocate(data(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    allocate(read_data(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    ! Fill the array with unique values per process
    data = real(mype + 1)

    ! Define MPI subarray type for writing
    sizes    = (/ 1+nx/2, ny, nz /)                           ! Global array size
    subsizes = (/ (cend(i)-cstart(i)+1, integer :: i=1,3)  /)        ! Local chunk size
    starts   = (/ (cstart(i)-1, integer :: i=1,3)  /)  ! Start position per rank

    call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, subarray_type, ierr)
    call MPI_Type_commit(subarray_type, ierr)

    ! Open file for writing
    call MPI_File_open(MPI_COMM_WORLD, "output_pencil.dat", MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

    ! Set file view
    offset = 0
    call MPI_File_set_view(file_handle, offset, MPI_DOUBLE_COMPLEX, subarray_type, "native", MPI_INFO_NULL, ierr)

    ! Write data collectively
    call MPI_File_write_all(file_handle, data, product(subsizes), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)

    ! Close file
    call MPI_File_close(file_handle, ierr)

    ! Synchronize before reading
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Open file for reading
    call MPI_File_open(MPI_COMM_WORLD, "output_pencil.dat", MPI_MODE_RDONLY+MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, file_handle, ierr)

    ! Set file view for reading
    call MPI_File_set_view(file_handle, offset, MPI_DOUBLE_COMPLEX, subarray_type, "native", MPI_INFO_NULL, ierr)

    ! Read data collectively
    call MPI_File_read_all(file_handle, read_data, product(subsizes), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)

    ! Close file
    call MPI_File_close(file_handle, ierr)

    ! Compare written and read data
    if (maxval(abs(data - read_data)) > 1.0e-10) then
        print *, "Rank", mype, "ERROR: Data mismatch after reading!"
    else
        print *, "Rank", mype, "SUCCESS: Data matches after reading!"
    end if

    ! Free MPI datatype
    call MPI_Type_free(subarray_type, ierr)

    ! Finalize MPI
    call MPI_FINALIZE(ierr)

end program mpi_io_pencil
