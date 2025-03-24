program exp3d

  use,intrinsic :: iso_c_binding
  
  use mpi
  use p3dfft
  
  implicit none

  ! mpif90 -g -O2 -r8  -L$P3DFFT_STRIDE/lib -L$FFTW_ROOT/lib -I$P3DFFT_STRIDE/include -I$FFTW_ROOT/lib -o exp3ds example_p3dfftstride.f90 -lp3dfft -lfftw3_mpi -lfftw3

  integer(4) :: nx0_big,ny0_big,nz0_big,nx=512
  integer(4) :: ierr
  integer(4) :: cstart(3),cend(3),csize(3),fstart(3),fend(3),fsize(3)
  complex(8), allocatable :: specarray(:,:,:)
  real(8), allocatable :: realarray(:,:,:)
  integer(4) :: i,j,k,repeat
  real(8) :: t1,t2
  integer(4) :: mype,np

  nx0_big = 3*nx/2
  ny0_big = 3*nx/2
  nz0_big = 3*nx/2  

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

  CALL p3dfft_setup([2,np/2],nx0_big,ny0_big,nz0_big,MPI_COMM_WORLD)

  t1 = MPI_WTIME()
  call p3dfft_get_dims(cstart,cend,csize,2)
  call p3dfft_get_dims(fstart,fend,fsize,1)
  t1 = MPI_WTIME()
  allocate(specarray(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
  allocate(realarray(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)))

  DO repeat = 1,500

        DO i = cstart(1),cend(1)
           DO k = cstart(3),cend(3)
              DO j = cstart(2),cend(2)
              
              specarray(i,j,k) = (i-1)/sqrt((2.0*i)**2.0 + (3.0*j)**3.0 + (4.0*k)**4.0)
              
           ENDDO
        ENDDO
     ENDDO

        DO i = cstart(1),cend(1)
           DO k = cstart(3),cend(3)
              DO j = cstart(2),cend(2)
              specarray(i,j,k) = cmplx(0.0,1.0) * k * specarray(i,j,k) &
                   - j * k * specarray(i,j,k)           
              
           ENDDO
        ENDDO
     ENDDO
     
     
     call p3dfft_btran_c2r(specarray,realarray,"fff")


     DO i = fstart(1),fend(1)
     DO k = fstart(3),fend(3)
        DO j = fstart(2),fend(2)
              
              realarray(i,j,k) = realarray(i,j,k)**2.0
              
           ENDDO
        ENDDO
     ENDDO
     
     call p3dfft_ftran_r2c(realarray,specarray,"fff")
DO i = cstart(1),cend(1)
     DO k = cstart(3),cend(3)
        DO j = cstart(2),cend(2)
              
              specarray(i,j,k) = cmplx(0.0,1.0/dble(nx0_big*ny0_big*nz0_big)) * j * specarray(i,j,k) &
                   - i * j * specarray(i,j,k)
              
           ENDDO
        ENDDO
     ENDDO
     
  ENDDO
  
  deallocate(specarray)
  deallocate(realarray)

  t2 = MPI_WTIME()

  if (mype.eq.0) print *, "Time for Operations ",t2-t1

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_FINALIZE(ierr)

end program exp3d
