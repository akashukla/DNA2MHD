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
  
  !COMPLEX :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  !COMPLEX :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  !COMPLEX :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  !COMPLEX :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
 

  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_lin                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_lin(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: which_term

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)

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

  INTEGER, INTENT(in) :: which_term

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)

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
       DO j = 0,ny0_big-1

          DO k = lkz1,lkz2

             ! Mahajan equation 14
             rhs_out_b(i,j,k,0) = i_complex*kzgrid(k)*(v_in(i,j,k,0) &
                  - hall*(i_complex*kygrid(j)*b_in(i,j,k,2) - i_complex*kzgrid(k)*b_in(i,j,k,1)))
             rhs_out_b(i,j,k,1) = i_complex*kzgrid(k)*(v_in(i,j,k,1) &
                  - hall*(i_complex*kzgrid(k)*b_in(i,j,k,0) - i_complex*kxgrid(i)*b_in(i,j,k,2)))
             rhs_out_b(i,j,k,2) = i_complex*kzgrid(k)*(v_in(i,j,k,2) &
                  - hall*(i_complex*kxgrid(i)*b_in(i,j,k,1) - i_complex*kygrid(j)*b_in(i,j,k,0)))
             
             ! Mahajan eqn 15b from prerana - no difference between linear versions 
             rhs_out_v(i,j,k,0) = i_complex*kzgrid(k)*b_in(i,j,k,0)
             rhs_out_v(i,j,k,1) = i_complex*kzgrid(k)*b_in(i,j,k,1)
             rhs_out_v(i,j,k,2) = i_complex*kzgrid(k)*b_in(i,j,k,2)
          END DO
          
       END DO

    ENDDO
    

if (nv) rhs_out_b = cmplx(0.0,0.0)

END SUBROUTINE get_rhs_lin1_ae

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_force                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_force(rhs_out_b, rhs_out_v)
  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  
    INTEGER :: i,j,k,h,ierr
    !REAL :: r
    REAL :: myresample,resample,th1,th2,th3,th4
    COMPLEX :: LW1,LC1,RW1,RC1
    REAL :: a,b,c,d,e,f,b1


    if (mype.eq.0) call random_number(resample)
    CALL MPI_BCAST(resample,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

    IF ((forcetype.eq.11).or.(forcetype.eq.12)) THEN
       
       DO i = 1,nkx0-1
          DO j = 1,ny0_big-1
             DO k = lkz1,lkz2
                if (verbose) c = MPI_WTIME()
                IF (forceb) THEN
                   rhs_out_b(i,j,k,0) = rhs_out_b(i,j,k,0) + (force_amp*random_normal() &
                        + i_complex*force_amp*random_normal())*mask1(i,j,k)
                   rhs_out_b(i,j,k,1) = rhs_out_b(i,j,k,1) + (force_amp*random_normal() &
                        + i_complex*force_amp*random_normal())*mask1(i,j,k)
                   rhs_out_b(i,j,k,2) = rhs_out_b(i,j,k,2) + (force_amp*random_normal() &
                        + i_complex*force_amp*random_normal())*mask1(i,j,k)
                   
                ENDIF
             
                rhs_out_v(i,j,k,0) = rhs_out_v(i,j,k,0) + (force_amp*random_normal() &
                     + i_complex*force_amp*random_normal())*mask1(i,j,k)
                rhs_out_v(i,j,k,1) = rhs_out_v(i,j,k,1) + (force_amp*random_normal() &
                     + i_complex*force_amp*random_normal())*mask1(i,j,k)
                rhs_out_v(i,j,k,2) = rhs_out_v(i,j,k,2) + (force_amp*random_normal() &
                     + i_complex*force_amp*random_normal())*mask1(i,j,k)
                
             ENDDO
          ENDDO
       ENDDO
       
    ENDIF

    IF ((resample.lt.dt).and.(forcetype.eq.21)) THEN

       DO i = 1,nkx0-1
          DO j = 1,ny0_big-1
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

                LW1 = LW1 * force_amp * sqrt(force_lw)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                LC1 = LC1 * force_amp * sqrt(force_lc)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                RW1 = RW1 * force_amp * sqrt(force_rw)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                RC1 = RC1 * force_amp * sqrt(force_rc)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                
                rhs_out_b(i,j,k,:) = rhs_out_b(i,j,k,:) &
                     + ((LW1 * alpha_leftwhist(i,j,k) + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                     -(RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k))&
                     * conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                rhs_out_v(i,j,k,:) = rhs_out_v(i,j,k,:) + ((LW1+LC1)*pcurleig(i,j,k,:) &
                     + (RW1+RC1)*conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

    IF (forcetype.gt.30) THEN

       ! Reset phases every half turnover time and interpolate phases linearly between turnover times

       IF (forcetype.eq.31) THEN

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

    ELSE IF (forcetype.eq.32) THEN

       CALL random_number(LWp)
       CALL random_number(LCp)
       CALL random_number(RWp)
       CALL random_number(RCp)
       
       LWp2 = LWp
       LCp2 = LCp
       RWp2 = RWp
       RCp2 = RCp
       
    ENDIF
    
    DO i = 1,nkx0-1
       DO j = 1,ny0_big-1
          DO k = lkz1,lkz2
             ! th1 = ((turnover - 2* time) * LWp(i,j,k) + (2*time) * LWp2(i,j,k))/turnover
             ! th2 = ((turnover - 2* time) * LCp(i,j,k) + (2*time) * LCp2(i,j,k))/turnover
             ! th3 = ((turnover - 2* time) * RWp(i,j,k) + (2*time) * RWp2(i,j,k))/turnover
             ! th4 = ((turnover - 2* time) * RCp(i,j,k) + (2*time) * RCp2(i,j,k))/turnover

             th1 = LWp(i,j,k)
             th2 = LCp(i,j,k)
             th3 = RWp(i,j,k)
             th4 = RCp(i,j,k)
             
             LW1 = exp(20.0*pi*i_complex*th1)
             LC1 = exp(20.0*pi*i_complex*th2)
             RW1 = exp(20.0*pi*i_complex*th3)
             RC1 = exp(20.0*pi*i_complex*th4)
             
             LW1 = LW1 * force_amp * sqrt(force_lw)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
             LC1 = LC1 * force_amp * sqrt(force_lc)/sqrt(force_lw + force_lc &
                  + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
             RW1 = RW1 * force_amp * sqrt(force_rw)/sqrt(force_lw + force_lc &
                  + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftwhist(i,j,k)**2)
                RC1 = RC1 * force_amp * sqrt(force_rc)/sqrt(force_lw + force_lc &
                     + force_rw + force_rc) * 1.0/sqrt(1 + alpha_leftcyclo(i,j,k)**2)
                
                rhs_out_b(i,j,k,:) = rhs_out_b(i,j,k,:) + ((LW1 * alpha_leftwhist(i,j,k) &
                     + LC1 * alpha_leftcyclo(i,j,k)) * pcurleig(i,j,k,:)&
                     -(RW1 * alpha_leftwhist(i,j,k) + RC1 * alpha_leftcyclo(i,j,k)) &
                     * conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                rhs_out_v(i,j,k,:) = rhs_out_v(i,j,k,:) + ((LW1+LC1)*pcurleig(i,j,k,:) &
                     + (RW1+RC1)*conjg(pcurleig(i,j,k,:)))*mask1(i,j,k)
                
             ENDDO
          ENDDO
       ENDDO
       
    END IF
    
  END SUBROUTINE get_rhs_force

  SUBROUTINE init_force

    ALLOCATE(mask(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(mask1(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

    ALLOCATE(LWp(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(LWp2(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

    ALLOCATE(LCp(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(LCp2(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

    ALLOCATE(RWp(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(RWp2(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

    ALLOCATE(RCp(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))
    ALLOCATE(RCp2(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2))

    mask = (kperps.lt.force_frac*maxval(kperps))
    ! masking will remove high k modes
    ! mask2 = (((i.le.nkxforce).and.(j.le.nkyforce)).and.(k.le.nkzforce)).and.(forcetype.eq.12))

    
    mask1 = 0

    if (mod(forcetype,2).eq.1) then
       DO i = 1,nkx0-1
          DO j = 1,ny0_big-1
             DO k = lkz1,lkz2
                
                if (mask(i,j,k)) mask1(i,j,k) = 1
             ENDDO
          ENDDO
       ENDDO
    else
       DO k  = 1,nkzforce
          IF ((k.ge.lkz1).and.(k.le.lkz2)) THEN
             DO i = 1,nkxforce
                DO j = 1,nkyforce
                   mask1(i,j,k) = 1
                   mask1(i,ny0_big-j,k) = 1
                ENDDO
             ENDDO
          ENDIF
          IF (((nz0_big-k).ge.lkz1).and.((nz0_big-k).le.lkz2)) THEN
             DO i = 1,nkxforce
                DO j = 1,nkyforce
                   mask1(i,j,nz0_big-k) = 1
                   mask1(i,ny0_big-j,nz0_big-k) = 1
                ENDDO
             ENDDO
          ENDIF
          
       ENDDO
    endif

    print *, "Force Mask", mype,maxval(mask1)
      
  END SUBROUTINE init_force

  SUBROUTINE finalize_force

    if (allocated(mask)) DEALLOCATE(mask)
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

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2,0:2)

  rhs_out_b = v_in
  rhs_out_v = - b_in

END SUBROUTINE get_rhs_test

SUBROUTINE get_rhs_diss(b_in,v_in,rhs_out_b,rhs_out_v)

  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: rhs_out_b(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: rhs_out_v(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)

  ! Explicit calculation of dissipation
  
  rhs_out_b = rhs_out_b - spread(eta * (kmags ** (2.0*hyp)),4,3) * b_in
  rhs_out_v = rhs_out_v - spread(vnu * (kmags ** (2.0*hyp)),4,3) * v_in
    
END SUBROUTINE get_rhs_diss

SUBROUTINE get_rhs_diss2(b_in,v_in)

  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: b_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: v_in(0:nx0_big/2,0:ny0_big-1,lkz1:lkz2, 0:2)

  ! Exact dissipation through integrating factor method
  
  b_in(:,:,:,0) = b_in(:,:,:,0)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,1) = b_in(:,:,:,1)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,2) = b_in(:,:,:,2)*exp(-eta * (kmags ** (2.0*hyp)) * dt)

  v_in(:,:,:,0) = v_in(:,:,:,0)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,1) = v_in(:,:,:,1)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,2) = v_in(:,:,:,2)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)

END SUBROUTINE get_rhs_diss2

END MODULE linear_rhs

