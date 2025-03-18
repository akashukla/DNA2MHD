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

  PUBLIC :: get_rhs_lin,get_rhs_lin2,finalize_force,init_force,get_rhs_diss,get_rhs_diss2,get_rhs_force,get_rhs_test,remove_div,hmhdnewton !,get_v_boundaries,get_v_boundaries2

  PRIVATE
  
  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_lin                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_lin(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: which_term

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)

  !If works for mu integrated as well for hankel/vperp version
  if (verbose.and.(mype.eq.0)) print *, "enter lin rhs"
  CALL get_rhs_lin1_ae(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
 
END SUBROUTINE get_rhs_lin

SUBROUTINE get_rhs_lin1_ae(b_in, v_in, rhs_out_b,rhs_out_v, which_term)

  INTEGER, INTENT(in) :: which_term

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)

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

 
  DO j = cstart(2),cend(2)
     DO k = cstart(3),cend(3)
        DO i=cstart(1),cend(1)
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

SUBROUTINE remove_div(b_in,v_in)

  !! Project divergence out of magnetic and velocity field perturbations

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  INTEGER :: i,j,k,l,h
  COMPLEX(C_DOUBLE_COMPLEX) :: div_v, div_b
  REAL(C_DOUBLE) :: k2
  COMPLEX(C_DOUBLE_COMPLEX) :: exb(0:2),exv(0:2)

  div_v = 0.0 +i_complex*0.0
  div_b = 0.0 +i_complex*0.0
  k2=0.0
  zero=0.0
  
  if (mype.eq.0) then
     exb = b_in(1,1,1,:)
     exv = v_in(1,1,1,:)
  endif
  
  DO i=cstart(1),cend(1)
     DO j=cstart(2),cend(2)
        DO k=cstart(3),cend(3)
           k2 = kxgrid(i)**2 + kygrid(j)**2 + kzgrid(k)**2
           div_v = kxgrid(i)*v_in(i,j,k,0) + kygrid(j)*v_in(i,j,k,1) + kzgrid(k)*v_in(i,j,k,2)
           div_b = kxgrid(i)*b_in(i,j,k,0) + kygrid(j)*b_in(i,j,k,1) + kzgrid(k)*b_in(i,j,k,2)
           
           v_in(i,j,k,0) = v_in(i,j,k,0) - div_v*kxgrid(i)/k2
           v_in(i,j,k,1) = v_in(i,j,k,1) - div_v*kygrid(j)/k2
           v_in(i,j,k,2) = v_in(i,j,k,2) - div_v*kzgrid(k)/k2
           
           ! The b equation is a curl, so we don't need to remove div b (except at start)
           b_in(i,j,k,0) = b_in(i,j,k,0) - div_b*kxgrid(i)/k2
           b_in(i,j,k,1) = b_in(i,j,k,1) - div_b*kygrid(j)/k2
           b_in(i,j,k,2) = b_in(i,j,k,2) - div_b*kzgrid(k)/k2
           
        ENDDO
     ENDDO
  ENDDO
  
  if (mype.eq.0) then
     b_in(1,1,1,:) = exb
     v_in(1,1,1,:) = exv
  endif
  
  if (verbose.and.(mype.eq.0)) print *,'Divergence Removed'
  
END SUBROUTINE remove_div



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_force                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_force(rhs_out_b, rhs_out_v)
  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), INTENT(inout) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  
    INTEGER :: i,j,k,h,ierr
    !REAL :: r
    REAL :: myresample,resample,th1,th2,th3,th4
    COMPLEX :: LW1,LC1,RW1,RC1
    REAL :: a,b,c,d,e,f,b1


    if (mype.eq.0) call random_number(resample)
    CALL MPI_BCAST(resample,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

    IF ((forcetype.eq.11).or.(forcetype.eq.12)) THEN

       DO i = cstart(1),cend(1)
       DO k = cstart(3),cend(3)
          DO j = cstart(2),cend(2)
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

    IF ((forcetype.ge.20)) THEN

       DO k = cstart(3),cend(3)
          DO i = cstart(1),cend(1)
             DO j = cstart(2),cend(2)
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

    IF (forcetype.ge.30) THEN

       
       CALL random_number(LWp)
       CALL random_number(LWp2)
       
       CALL random_number(LCp)
       CALL random_number(LCp2)
       
       CALL random_number(RWp)
       CALL random_number(RWp2)
       
       CALL random_number(RCp)
       CALL random_number(RCp2)          
          
       DO i = cstart(1),cend(1)
          DO j = cstart(2),cend(2)
             DO k = cstart(3),cend(3)
                !th1 = ((turnover - 2* time) * LWp(i,j,k) + (2*time) * LWp2(i,j,k))/turnover
                !th2 = ((turnover - 2* time) * LCp(i,j,k) + (2*time) * LCp2(i,j,k))/turnover
                !th3 = ((turnover - 2* time) * RWp(i,j,k) + (2*time) * RWp2(i,j,k))/turnover
                !th4 = ((turnover - 2* time) * RCp(i,j,k) + (2*time) * RCp2(i,j,k))/turnover

                
                LW1 = exp(20.0*pi*i_complex*LWp(i,j,k))*1.2533
                LC1 = exp(20.0*pi*i_complex*LCp(i,j,k))*1.2533
                RW1 = exp(20.0*pi*i_complex*RWp(i,j,k))*1.2533
                RC1 = exp(20.0*pi*i_complex*RCp(i,j,k))*1.2533
                
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
       
    ENDIF
    
  END SUBROUTINE get_rhs_force

  SUBROUTINE init_force

    implicit none
    integer(4) :: force_minx,force_maxx,force_miny,force_maxy,force_minz,force_maxz
    integer(4) :: i,j,k

    ALLOCATE(mask(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    ALLOCATE(mask1(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    ALLOCATE(LWp(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    ALLOCATE(LWp2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    ALLOCATE(LCp(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    ALLOCATE(LCp2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    ALLOCATE(RWp(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    ALLOCATE(RWp2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    ALLOCATE(RCp(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))
    ALLOCATE(RCp2(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3)))

    mask = (kperps.lt.force_frac*maxval(kperps))
    ! masking will remove high k modes
    ! mask2 = (((i.le.nkxforce).and.(j.le.nkyforce)).and.(k.le.nkzforce)).and.(forcetype.eq.12))
    
    mask1 = 0

    DO i = xst,cend(1)
       DO j = yst,cend(2)
          DO k = zst,cend(3)
             if ((mod(forcetype,2).eq.1).and.mask(i,j,k)) mask1(i,j,k) = 1
             if ((mod(forcetype,2).eq.0).and.(((i-1.le.nkxforce).and.(j-1.le.nkyforce.or.ny0_big+1-j.le.nkyforce))&
                  .and.(k-1.le.nkzforce.or.nz0_big+1-k.le.nkzforce))) mask1(i,j,k) = 1
          ENDDO
       ENDDO
    ENDDO
   
    ! print *, "Force Mask", mype,maxval(mask1)
      
  END SUBROUTINE init_force

  SUBROUTINE finalize_force

    if (allocated(LWp)) DEALLOCATE(LWp)
    if (allocated(LWp2)) DEALLOCATE(LWp2)

    if (allocated(LCp)) DEALLOCATE(LCp)
    if (allocated(LCp2)) DEALLOCATE(LCp2)

    if (allocated(RWp)) DEALLOCATE(RWp)
    if (allocated(RWp2)) DEALLOCATE(RWp2)

    if (allocated(RCp)) DEALLOCATE(RCp)
    if (allocated(RCp2)) DEALLOCATE(RCp2)

    if (allocated(mask)) DEALLOCATE(mask)
    if (allocated(mask1)) DEALLOCATE(mask1)
    
  END SUBROUTINE finalize_force
  
SUBROUTINE get_rhs_test(b_in,v_in,rhs_out_b,rhs_out_v)

  IMPLICIT NONE

  COMPLEX(C_DOUBLE_COMPLEX) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  COMPLEX(C_DOUBLE_COMPLEX) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  rhs_out_b = v_in
  rhs_out_v = - b_in

END SUBROUTINE get_rhs_test

SUBROUTINE get_rhs_diss(b_in,v_in,rhs_out_b,rhs_out_v)

  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(in) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: rhs_out_b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: rhs_out_v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)

  ! Explicit calculation of dissipation
  
  rhs_out_b = rhs_out_b - spread(eta * (kmags ** (2.0*hyp)),4,3) * b_in
  rhs_out_v = rhs_out_v - spread(vnu * (kmags ** (2.0*hyp)),4,3) * v_in
    
END SUBROUTINE get_rhs_diss

SUBROUTINE get_rhs_diss2(b_in,v_in)

  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: b_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)
  COMPLEX(C_DOUBLE_COMPLEX), intent(inout) :: v_in(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3), 0:2)

  ! Exact dissipation through integrating factor method
  
  b_in(:,:,:,0) = b_in(:,:,:,0)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,1) = b_in(:,:,:,1)*exp(-eta * (kmags ** (2.0*hyp)) * dt)
  b_in(:,:,:,2) = b_in(:,:,:,2)*exp(-eta * (kmags ** (2.0*hyp)) * dt)

  v_in(:,:,:,0) = v_in(:,:,:,0)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,1) = v_in(:,:,:,1)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)
  v_in(:,:,:,2) = v_in(:,:,:,2)*exp(-vnu * (kmags ** (2.0*hyp)) * dt)

END SUBROUTINE get_rhs_diss2

SUBROUTINE hmhdnewton(step0b,step0v,step12b,step12v)

  ! Given the results of fixed point iteration step0b,step0v,step12b,step12v
  ! Solve for the next iteration based on the linear dynamics of the system
  ! i.e. given x0 and F(x0)
  ! We exploit the known normal mode decomposition to perform the matrix inversion
  
  IMPLICIT NONE

  complex(8), intent(inout) :: step0b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  complex(8), intent(inout) :: step0v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  complex(8), intent(in) :: step12b(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)
  complex(8), intent(in) :: step12v(cstart(1):cend(1),cstart(2):cend(2),cstart(3):cend(3),0:2)

  complex(8) :: LWk,LCk,RWk,RCk,stepdifferenceb(0:2),stepdifferencev(0:2)
  integer(4) :: i,j,k
  integer(4) :: zerofoundm = 0,zerofound,ierr,zerodiffm = 0,zerodiff
  real(8) :: maxdevmb = 0,maxdevb,maxdevmv = 0,maxdevv

  maxdevmb = maxval(abs(step12b-step0b))
  maxdevmv = maxval(abs(step12v-step0v))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(maxdevmb,maxdevb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(maxdevmv,maxdevv,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if ((verbose).and.(mype.eq.0)) print *, "Max Deviation",maxdevb,maxdevv
  
  DO i = cstart(1),cend(1)
     DO j = cstart(2),cend(2)
        DO k = cstart(3),cend(3)
           stepdifferenceb = step12b(i,j,k,:)-step0b(i,j,k,:)
           stepdifferencev = step12v(i,j,k,:)-step0v(i,j,k,:)
           if (max(maxval(abs(stepdifferenceb)),maxval(abs(stepdifferencev))).gt.0) zerodiffm = 1

           ! Compute dot products to get normal mode decomposition

           LWk = sum(conjg(pcurleig(i,j,k,:))*(alpha_leftwhist(i,j,k)*stepdifferenceb + stepdifferencev))/(sqrt(1.0+alpha_leftwhist(i,j,k)**2.0))
           LCk = sum(conjg(pcurleig(i,j,k,:))*(-stepdifferenceb/alpha_leftwhist(i,j,k) + stepdifferencev))/(sqrt(1.0+(1.0/alpha_leftwhist(i,j,k))**2.0))
           RWk = sum(pcurleig(i,j,k,:)*(-alpha_leftwhist(i,j,k)*stepdifferenceb + stepdifferencev))/(sqrt(1.0+alpha_leftwhist(i,j,k)**2.0))
           RCk = sum(pcurleig(i,j,k,:)*(stepdifferenceb/alpha_leftwhist(i,j,k) + stepdifferencev))/(sqrt(1.0+(1.0/alpha_leftwhist(i,j,k))**2.0))

           ! Invert for correction normal mode amplitudes
           LWk = LWk/(cmplx(1.0,-0.5*dt*kzgrid(k)*alpha_leftwhist(i,j,k)))
           LCk = LCk/(cmplx(1.0,0.5*dt*kzgrid(k)/alpha_leftwhist(i,j,k)))
           RWk = RWk/(cmplx(1.0,0.5*dt*kzgrid(k)*alpha_leftwhist(i,j,k)))
           RCk = RCk/(cmplx(1.0,-0.5*dt*kzgrid(k)/alpha_leftwhist(i,j,k)))

           if (max(abs(LWk),abs(LCk),abs(RWk),abs(RCk)).gt.0) zerofoundm = 1

           ! Add normal mode corrections
           step0b(i,j,k,:) = step0b(i,j,k,:) + LWk * alpha_leftwhist(i,j,k) * pcurleig(i,j,k,:)/sqrt(alpha_leftwhist(i,j,k)**2.0+1.0)
           step0b(i,j,k,:) = step0b(i,j,k,:) - LCk * 1.0/alpha_leftwhist(i,j,k) * pcurleig(i,j,k,:)/sqrt((1.0/alpha_leftwhist(i,j,k))**2.0+1.0)
           step0b(i,j,k,:) = step0b(i,j,k,:) - RWk * alpha_leftwhist(i,j,k) * conjg(pcurleig(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2.0+1.0)
           step0b(i,j,k,:) = step0b(i,j,k,:) + RCk * 1.0/alpha_leftwhist(i,j,k) * conjg(pcurleig(i,j,k,:))/sqrt((1.0/alpha_leftwhist(i,j,k))**2.0+1.0)

           step0v(i,j,k,:) = step0v(i,j,k,:) + LWk * pcurleig(i,j,k,:)/sqrt(alpha_leftwhist(i,j,k)**2.0+1.0)
           step0v(i,j,k,:) = step0v(i,j,k,:) + LCk * pcurleig(i,j,k,:)/sqrt((1.0/alpha_leftwhist(i,j,k))**2.0+1.0)
           step0v(i,j,k,:) = step0v(i,j,k,:) + RWk * conjg(pcurleig(i,j,k,:))/sqrt(alpha_leftwhist(i,j,k)**2.0+1.0)
           step0v(i,j,k,:) = step0v(i,j,k,:) + RCk * conjg(pcurleig(i,j,k,:))/sqrt((1.0/alpha_leftwhist(i,j,k))**2.0+1.0)

        ENDDO
     ENDDO
  ENDDO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(zerofoundm,zerofound,1,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
  if ((verbose).and.(mype.eq.0).and.zerofound.gt.0) print *, "Only trivial solution found ",n_mpi_procs-zerofound,"nodes"
  if ((verbose).and.zerofoundm.eq.0) print * , "Trivial found", cstart(1),cstart(2)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(zerodiffm,zerodiff,1,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
  if ((verbose).and.(mype.eq.0).and.zerodiff.gt.0) print *, "All differences found equal ",n_mpi_procs-zerodiff,"nodes"
  if ((verbose).and.zerodiffm.eq.0) print * , "Trivial found",	cstart(1),cstart(2)

  maxdevmb = maxval(abs(step12b-step0b))
  maxdevmv = maxval(abs(step12v-step0v))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(maxdevmb,maxdevb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(maxdevmv,maxdevv,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

  if ((verbose).and.(mype.eq.0)) print *, "Max Deviation",maxdevb,maxdevv

  if ((verbose).and.zerofoundm.eq.1) print *, cstart(1),cstart(2),maxval(abs(step12b-step0b)),maxloc(abs(step12b-step0b))
  

END SUBROUTINE hmhdnewton

END MODULE linear_rhs
