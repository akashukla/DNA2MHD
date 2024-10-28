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

  PUBLIC :: get_rhs_lin,get_rhs_lin2,finalize_force,init_force,get_rhs_diss,get_rhs_diss2,get_rhs_force,get_rhs_test !,get_v_boundaries,get_v_boundaries2

  PRIVATE
  
  LOGICAL, ALLOCATABLE :: mask1(:,:,:)

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
    if (verbose) print *, "enter lin rhs"
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

 rhs_out_b=cmplx(0.0,0.0)
 rhs_out_v=cmplx(0.0,0.0)

 !IF(verbose.and.mype==0) WRITE(*,*) "get_rhs_lin1", 68

 DO i=0,nkx0-1
   DO j=0,nky0-1
     DO k=lkz1,lkz2
        ! L = kxgrid(i)*kxgrid(i) + kygrid(j)*kygrid(j) + kzgrid(k)*kzgrid(k)
        L = 0
        
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

END SUBROUTINE get_rhs_lin1_ae



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_force                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_force(rhs_out_b, rhs_out_v)
    IMPLICIT NONE
    INTEGER :: i,j,k,h,ierr
    !REAL :: r
    REAL :: myresample,resample,th1,th2,th3,th4
    COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
    COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
    COMPLEX :: LW1,LC1,RW1,RC1
    REAL :: a,b,c,d,e,f,b1

    if (verbose) a = MPI_WTIME()
    if (mype.eq.0) call random_number(resample)
    if (verbose) b = MPI_WTIME()
    CALL MPI_BCAST(resample,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

    mask1 = (kperps.lt.force_frac*maxval(kperps))
    ! mask2 = (((i.le.nkxforce).and.(j.le.nkyforce)).and.(k.le.nkzforce)).and.(forcetype.eq.12))
    
    IF ((resample.lt.dt).and.(forcetype.eq.11)) THEN

       DO i = 1,nkx0-1
          DO j = 1,nky0-1
             DO k = lkz1,lkz2
                if (verbose) c = MPI_WTIME()
                IF (forceb) THEN
                   rhs_out_b(i,j,k,0) = rhs_out_b(i,j,k,0) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,j,k)
                   rhs_out_b(i,j,k,1) = rhs_out_b(i,j,k,1) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,j,k)
                   
                   rhs_out_b(i,nky0-j,k,0) = rhs_out_b(i,nky0-j,k,0) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,nky0-j,k)
                   rhs_out_b(i,nky0-j,k,1) = rhs_out_b(i,nky0-j,k,1) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,nky0-j,k)

                ENDIF
             
                rhs_out_v(i,j,k,0) = rhs_out_v(i,j,k,0) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,j,k)
                rhs_out_v(i,j,k,1) = rhs_out_v(i,j,k,1) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,j,k)
                
                rhs_out_v(i,nky0-j,k,0) = rhs_out_v(i,nky0-j,k,0) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,nky0-j,k)
                rhs_out_v(i,nky0-j,k,1) = rhs_out_v(i,nky0-j,k,1) + (force_amp*random_normal() + i_complex*force_amp*random_normal())*mask1(i,nky0-j,k)

                if (verbose) d = MPI_WTIME()
                if (verbose) e = MPI_WTIME()
             ENDDO
          ENDDO
       ENDDO
       if (verbose) f = MPI_WTIME()
       
       if (verbose) print *, "Forced"
       if (verbose) print *, "Trial Time ",b-a
       if (verbose) print *, "Rands Time / Mode ",d-c
       if (verbose) print *, "Assignment Time / Mode",e-d
       if (verbose) print *, "Total Time ",f-a
    ENDIF

    IF ((resample.lt.dt).and.(forcetype.eq.21)) THEN

       DO i = 1,nkx0-1
          DO j = 1,nky0-1
             DO k = lkz1,lkz2
                if (verbose) c = MPI_WTIME()
                CALL random_number(th1)
                CALL random_number(th2)
                CALL random_number(th3)
                CALL random_number(th4)
                
                LW1 = exp(20.0*pi*i_complex*th1) !(random_normal()+i_complex*random_normal())
                LC1 = exp(20.0*pi*i_complex*th2) !(random_normal()+i_complex*random_normal())
                RW1 = exp(20.0*pi*i_complex*th3) !(random_normal()+i_complex*random_normal())
                RC1 = exp(20.0*pi*i_complex*th4) !(random_normal()+i_complex*random_normal())
                if (verbose) d = MPI_WTIME()

                LW1 = LW1 * force_amp * sqrt(force_lw)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                LC1 = LC1 * force_amp * sqrt(force_lc)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                RW1 = RW1 * force_amp * sqrt(force_rw)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                RC1 = RC1 * force_amp * sqrt(force_rc)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                
                if (verbose) e = MPI_WTIME()
                rhs_out_b(i,j,k,:) = rhs_out_b(i,j,k,:) + ((LW1 * alpha_leftwhist(i,j,k) + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                     -(RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k)) * conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                rhs_out_v(i,j,k,:) = rhs_out_v(i,j,k,:) + ((LW1+LC1)*pcurleig(i,j,k,:) + (RW1+RC1)*conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                if (verbose) f = MPI_WTIME()
             ENDDO
          ENDDO
       ENDDO

       if (verbose) print *, "Forced"
       if (verbose) print *, "Trial Time ",b-a
       if (verbose) print *, "Rands Time / Mode ",d-c
       if (verbose) print *, "Assignment Time / Mode",e-d
       if (verbose) print *, "Total Time ",f-a

    ENDIF

    IF (forcetype.eq.31) THEN ! continuous mode injection in range

       ! Reset phases every half turnover time and interpolate phases linearly between turnover times

       if (time.eq.0.and.(time > last_reset)) then

          last_reset = time
          CALL random_number(LWp)
          CALL random_number(LWp2)

          CALL random_number(LCp)
          CALL random_number(LCp2)

          CALL random_number(RWp)
          CALL random_number(RWp2)

          CALL random_number(RCp)
          CALL random_number(RCp2)          
          
       else if ((mod(time,turnover/2).lt.dt).and.(time > last_reset)) then

          last_reset = time
          LWp = LWp2
          LCp = LCp2
          RWp = RWp2
          RCp = RCp2

          CALL random_number(LWp2)
          CALL random_number(LCp2)
          CALL random_number(RWp2)
          CALL random_number(RCp2)          
          
       endif

      DO i = 1,nkx0-1
          DO j = 1,nky0-1
             DO k = lkz1,lkz2
                th1 = ((turnover - 2* time) * LWp(i,j,k) + (2*time) * LWp2(i,j,k))/turnover
                th2 = ((turnover - 2* time) * LCp(i,j,k) + (2*time) * LCp2(i,j,k))/turnover
                th3 = ((turnover - 2* time) * RWp(i,j,k) + (2*time) * RWp2(i,j,k))/turnover
                th4 = ((turnover - 2* time) * RCp(i,j,k) + (2*time) * RCp2(i,j,k))/turnover

                LW1 = exp(20.0*pi*i_complex*th1)
                LC1 = exp(20.0*pi*i_complex*th2)
                RW1 = exp(20.0*pi*i_complex*th3)
                RC1 = exp(20.0*pi*i_complex*th4)
                
                LW1 = LW1 * force_amp * sqrt(force_lw)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                LC1 = LC1 * force_amp * sqrt(force_lc)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                RW1 = RW1 * force_amp * sqrt(force_rw)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                RC1 = RC1 * force_amp * sqrt(force_rc)/sqrt(force_lw + force_lc + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)

                rhs_out_b(i,j,k,:) = rhs_out_b(i,j,k,:) + ((LW1 * alpha_leftwhist(i,j,k) + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                     -(RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k)) * conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                rhs_out_v(i,j,k,:) = rhs_out_v(i,j,k,:) + ((LW1+LC1)*pcurleig(i,j,k,:) + (RW1+RC1)*conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)

             ENDDO
          ENDDO
       ENDDO
       
    END IF
    
  END SUBROUTINE get_rhs_force

  SUBROUTINE init_force

    ALLOCATE(mask1(0:nkx0-1,0:nky0-1,lkz1:lkz2))

    ALLOCATE(LWp(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(LWp2(0:nkx0-1,0:nky0-1,lkz1:lkz2))

    ALLOCATE(LCp(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(LCp2(0:nkx0-1,0:nky0-1,lkz1:lkz2))

    ALLOCATE(RWp(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(RWp2(0:nkx0-1,0:nky0-1,lkz1:lkz2))

    ALLOCATE(RCp(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    ALLOCATE(RCp2(0:nkx0-1,0:nky0-1,lkz1:lkz2))
    
  END SUBROUTINE init_force

  SUBROUTINE finalize_force

    if (allocated(mask1)) DEALLOCATE(mask1)

    if (allocated(LWp)) DEALLOCATE(LWp)
    if (allocated(LWp2)) DEALLOCATE(LWp2)

    if (allocated(LCp)) DEALLOCATE(LCp)
    if (allocated(LCp2)) DEALLOCATE(LCp2)

    if (allocated(RWp)) DEALLOCATE(RWp)
    if (allocated(RWp2)) DEALLOCATE(RWp2)

    if (allocated(RCp)) DEALLOCATE(RCp)
    if (allocated(RCp2)) DEALLOCATE(RCp2)
    
  END SUBROUTINE finalize_force
  
SUBROUTINE get_rhs_test(b_in,v_in,rhs_out_b,rhs_out_v)

  IMPLICIT NONE

  COMPLEX,INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)


  rhs_out_b = v_in
  rhs_out_v = - b_in

END SUBROUTINE get_rhs_test

SUBROUTINE get_rhs_diss(b_in,v_in,rhs_out_b,rhs_out_v)

  COMPLEX,INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

  rhs_out_b = rhs_out_b - spread(eta * (kmags ** (2.0*hyp)),4,3) * b_in
  rhs_out_v = rhs_out_v - spread(vnu * (kmags ** (2.0*hyp)),4,3) * v_in
    
END SUBROUTINE get_rhs_diss

SUBROUTINE get_rhs_diss2(b_in,v_in)

  COMPLEX,INTENT(inout) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)
  COMPLEX,INTENT(inout) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,0:2)

  ! b_in = b_in * spread(exp(-eta * (kmags ** (2.0*hyp)) * dt),4,3)
  ! v_in = v_in * spread(exp(-vnu * (kmags ** (2.0*hyp)) * dt),4,3)

  b_in(:,:,:,0) = b_in(:,:,:,0)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,1) = b_in(:,:,:,1)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,2) = b_in(:,:,:,2)*exp(-eta * (kmags ** (2.0*hyp)) * dt)

  v_in(:,:,:,0) = v_in(:,:,:,0)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,1) = v_in(:,:,:,1)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,2) = v_in(:,:,:,2)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)

END SUBROUTINE get_rhs_diss2

END MODULE linear_rhs

