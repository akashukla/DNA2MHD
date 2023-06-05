!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 29/12/2012                                                                !!
!!                             cc_get_rhs_lin.f90                            !!
!!                                                                           !!
!!  linear_rhs                                                               !!
!!  -- get_rhs_lin                                                           !!
!!  -- get_rhs_lin1                                                          !!
!!  -- get_rhs_lin2                                                          !!
!!  -- get_v_boundaries2                                                     !!
!!                                                                     1.000 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                           !!
!!                                linear_rhs                                 !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE linear_rhs
  USE par_mod
  USE mpi
  USE random
  !USE flr_effects
  !USE hk_effects

  PUBLIC :: get_rhs_lin,get_rhs_lin2 !,get_v_boundaries,get_v_boundaries2


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_lin                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_lin(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
  IMPLICIT NONE
 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 INTEGER, INTENT(in) :: which_term

 IF ((rhs_lin_version==1).or.(rhs_lin_version==12)) THEN
   !If works for mu integrated as well for hankel/vperp version
   CALL get_rhs_lin1_ae(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
 ELSE IF(rhs_lin_version==2) THEN
   !CALL get_rhs_lin2(g_in,g_bounds,phi_in,rhs_out,which_term)
   STOP 'get_rhs_lin2 needs to be benchmarked and updated to 6D g.'
 END IF
 
END SUBROUTINE get_rhs_lin

SUBROUTINE get_rhs_lin1_ae(b_in, v_in, rhs_out_b,rhs_out_v, which_term)

 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 INTEGER, INTENT(in) :: which_term

 INTEGER :: i,j,k,h,ierr
 !for transpose for left ev's
 INTEGER :: grad1_flag
 INTEGER :: grad2_flag
 COMPLEX :: phi_mod1,phi_mod2,g0_bcast
 COMPLEX :: g_closure
 REAL :: L
 REAL :: Larr(0:nkx0-1,0:nky0-1,lkz1:lkz2)

 rhs_out_b=cmplx(0.0,0.0)
 rhs_out_v=cmplx(0.0,0.0)

 !IF(verbose.and.mype==0) WRITE(*,*) "get_rhs_lin1", 68
 DO i=0,nkx0-1
   DO j=0,nky0-1
     DO k=lkz1,lkz2
        L = kxgrid(i)*kxgrid(i) + kygrid(j)*kygrid(j) + kzgrid(k)*kzgrid(k)
        Larr(i,j,k) = L
     !eqn 14
        !rhs_out_b(i,j,k,0) = i_complex*kzgrid(k)*v_in(i,j,k,0) + i_complex*kygrid(j)*b_in(i,j,k,2) -i_complex*kzgrid(k)*b_in(i,j,k,1)
        !rhs_out_b(i,j,k,1) = i_complex*kzgrid(k)*v_in(i,j,k,1) + i_complex*kzgrid(k)*b_in(i,j,k,0) -i_complex*kxgrid(i)*b_in(i,j,k,2)
        !rhs_out_b(i,j,k,2) = i_complex*kzgrid(k)*v_in(i,j,k,2) + i_complex*kxgrid(i)*b_in(i,j,k,1) -i_complex*kygrid(y)*b_in(i,j,k,0)

        rhs_out_b(i,j,k,0) = i_complex*kzgrid(k)*(v_in(i,j,k,0) & 
             - hall*(i_complex*kygrid(j)*b_in(i,j,k,2) - i_complex*kzgrid(k)*b_in(i,j,k,1))) - eta*(L**hyp)*b_in(i,j,k,0)
        rhs_out_b(i,j,k,1) = i_complex*kzgrid(k)*(v_in(i,j,k,1) &
             - hall*(i_complex*kzgrid(k)*b_in(i,j,k,0) - i_complex*kxgrid(i)*b_in(i,j,k,2))) - eta*(L**hyp)*b_in(i,j,k,1)
        rhs_out_b(i,j,k,2) = i_complex*kzgrid(k)*(v_in(i,j,k,2) &
             - hall*(i_complex*kxgrid(i)*b_in(i,j,k,1) - i_complex*kygrid(j)*b_in(i,j,k,0))) - eta*(L**hyp)*b_in(i,j,k,2)
        ! if (verbose) print *, 'b lin equation stored'
        ! if (verbose) print *, 'L, eta',L,eta
      !Eqn 15
if (rhs_lin_version==1) then
        rhs_out_v(i,j,k,0) = i_complex*kzgrid(k)*b_in(i,j,k,0)-i_complex*kxgrid(i)*b_in(i,j,k,2) - vnu*(L**hyp)*v_in(i,j,k,0)
        rhs_out_v(i,j,k,1) = i_complex*kzgrid(k)*b_in(i,j,k,1)-i_complex*kygrid(j)*b_in(i,j,k,2) - vnu*(L**hyp)*v_in(i,j,k,1)
        rhs_out_v(i,j,k,2) = - vnu*(L**hyp)*v_in(i,j,k,2)
        ! if (verbose) print *, 'v1 lin equation stored'
        ! if (verbose) print *, 'L, vnu',L,vnu
endif
    ! Eq 15 v2 from prerana
if (rhs_lin_version==12) then
        rhs_out_v(i,j,k,0) = i_complex*kzgrid(k)*b_in(i,j,k,0) - vnu*(L**hyp)*v_in(i,j,k,0)
        rhs_out_v(i,j,k,1) = i_complex*kzgrid(k)*b_in(i,j,k,1) - vnu*(L**hyp)*v_in(i,j,k,1)
        rhs_out_v(i,j,k,2) = i_complex*kzgrid(k)*b_in(i,j,k,2) - vnu*(L**hyp)*v_in(i,j,k,2)
        ! if (verbose) print *, 'v12 lin equation stored'
        ! if (verbose) print *, 'L, vnu',L,vnu
endif
     END DO
   END DO
 END DO 

if (nv) rhs_out_b = cmplx(0.0,0.0)

if (plot_nls.and.(mod(itime,istep_energy).eq.0)) then
WRITE(dbio) eta*(Larr**hyp)*b_in(:,:,:,0)
WRITE(dbio) eta*(Larr**hyp)*b_in(:,:,:,1)
WRITE(dbio) eta*(Larr**hyp)*b_in(:,:,:,2)
WRITE(dvio) vnu*(Larr**hyp)*v_in(:,:,:,0)
WRITE(dvio) vnu*(Larr**hyp)*v_in(:,:,:,1)
WRITE(dvio) vnu*(Larr**hyp)*v_in(:,:,:,2)

if (verbose.and.(mype.eq.0)) print *, 'Dissipation written'
endif

END SUBROUTINE get_rhs_lin1_ae



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_force                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_force(rhs_out_b, rhs_out_v, dt)
    IMPLICIT NONE
    INTEGER :: i,j,k,h,ierr
    !REAL :: r
    REAL :: myresample,resample
    REAL, INTENT(in) :: dt
    COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
    COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)

    if (mype.eq.0) call random_number(resample)
    CALL MPI_BCAST(resample,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

    IF(resample.lt.dt) THEN
!        print *, "Resample value ",resample
        DO i=1,nkxforce
          DO j=1,nkyforce
            DO k=1,nkzforce
                IF (forceb) THEN
                rhs_out_b(i,j,k,0) = rhs_out_b(i,j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_b(i,j,k,1) = rhs_out_b(i,j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()

                rhs_out_b(nkx0-i,nky0-j,k,0) = rhs_out_b(nkx0-i,nky0-j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_b(nkx0-i,nky0-j,k,1) = rhs_out_b(nkx0-i,nky0-j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()
        
                rhs_out_b(i,nky0-j,k,0) = rhs_out_b(i,nky0-j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_b(i,nky0-j,k,1) = rhs_out_b(i,nky0-j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()

                rhs_out_b(nkx0-i,j,k,0) = rhs_out_b(nkx0-i,j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_b(nkx0-i,j,k,1) = rhs_out_b(nkx0-i,j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                ENDIF

                rhs_out_v(i,j,k,0) = rhs_out_v(i,j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_v(i,j,k,1) = rhs_out_v(i,j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()

                rhs_out_v(nkx0-i,nky0-j,k,0) = rhs_out_v(nkx0-i,nky0-j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_v(nkx0-i,nky0-j,k,1) = rhs_out_v(nkx0-i,nky0-j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()

                rhs_out_v(i,nky0-j,k,0) = rhs_out_v(i,nky0-j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_v(i,nky0-j,k,1) = rhs_out_v(i,nky0-j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()

                rhs_out_v(nkx0-i,j,k,0) = rhs_out_v(nkx0-i,j,k,0) + force_amp*random_normal() + i_complex*force_amp*random_normal()
                rhs_out_v(nkx0-i,j,k,1) = rhs_out_v(nkx0-i,j,k,1) + force_amp*random_normal() + i_complex*force_amp*random_normal()
            END DO
          END DO
        END DO
    ENDIF
END SUBROUTINE get_rhs_force

END MODULE linear_rhs
