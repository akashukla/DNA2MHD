program fftw613

  use, intrinsic :: iso_c_binding
  ! Compile with ftn -Wl,-rpath,$FFTW_ROOT/lib -L$FFTW_ROOT/lib -lfftw3_mpi -lfftw3 -o fftw613ex.o fftw613ex.f90 
  include 'fftw3-mpi.f03'
  include 'mpif.h'
  
  integer(C_INTPTR_T), parameter :: L = 8
  integer(C_INTPTR_T), parameter :: M = 8
  integer(C_INTPTR_T), parameter :: N = 8
  type(C_PTR) :: plan, plani, cdata,rdata
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:)
  real(C_DOUBLE), pointer :: rdataf(:,:,:)
  integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset,k
  integer :: ierr

  !   get local data size and allocate (note dimension reversal)

  CALL MPI_INIT(ierr)
  
  CALL fftw_mpi_init()
  
  alloc_local = fftw_mpi_local_size_3d(N,M, L/2+1, MPI_COMM_WORLD, &
                                       local_M, local_j_offset)
  cdata = fftw_alloc_complex(alloc_local)
  rdata = fftw_alloc_real(2*alloc_local)
  call c_f_pointer(cdata, data, [L/2+1,M,local_M])
  call c_f_pointer(rdata,rdataf,[2*(L/2+1),M,local_M])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  plan = fftw_mpi_plan_dft_c2r_3d(N,M, L, data, rdataf, MPI_COMM_WORLD, &
       FFTW_PATIENT)
  CALL fftw_destroy_plan(plan)
  plan = fftw_mpi_plan_dft_c2r_3d(N,M, L, data, rdataf, MPI_COMM_WORLD, &
       FFTW_PATIENT)
  plani = fftw_mpi_plan_dft_r2c_3d(N,M,L,rdataf,data,MPI_COMM_WORLD,FFTW_PATIENT)

! initialize data to some function my_function(i,j)
  do k = 1, local_M
     do j = 1, M
        do i = 1,L/2
           data(i, j,k) = exp**(-cmplx(1.0,0.0)*(j+local_j_offset))/i**2 
        end do
     end do
  end do
  

  ! compute transform (as many times as desired)
  call fftw_mpi_execute_dft_c2r(plan, data, rdataf)
  call fftw_mpi_execute_dft_r2c(plani,rdataf,data)

  call fftw_destroy_plan(plan)
  call fftw_destroy_plan(plani)

  call fftw_free(cdata)
  call fftw_free(rdata)

  call fftw_mpi_cleanup()

  CALL MPI_FINALIZE(ierr)
  
end program fftw613
