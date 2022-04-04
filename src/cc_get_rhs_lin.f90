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
  USE flr_effects
  USE hk_effects

  PUBLIC :: get_rhs_lin,get_rhs_lin2 !,get_v_boundaries,get_v_boundaries2


  CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_lin                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_rhs_lin(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
  IMPLICIT NONE

                                                                          ! COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
                                                                          ! COMPLEX :: g_bounds(0:nkx0-1,0:nky0-1,lkz1:lkz2,lh1:lh2,2)
                                                                          ! !COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2)
                                                                          ! COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
                                                                          ! !COMPLEX :: g_wb(0:nkx0-1,0:nky0-1,0:nkz0-1,lbv:ubv)
                                                                          ! INTEGER, INTENT(in) :: which_term


 !COMPLEX, INTENT(in) :: g_bounds(0:nkx0-1,0:nky0-1,lkz1:lkz2,lh1:lh2,2)
 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 INTEGER, INTENT(in) :: which_term


                                                                                        ! IF(np_herm.gt.1) THEN
                                                                                        !    CALL get_v_boundaries2(g_in,g_bounds)
                                                                                        ! END IF

 IF(rhs_lin_version==1) THEN
   !If works for mu integrated as well for hankel/vperp version
   CALL get_rhs_lin1_ae(b_in, v_in, rhs_out_b, rhs_out_v, which_term)
 ELSE IF(rhs_lin_version==2) THEN
   !CALL get_rhs_lin2(g_in,g_bounds,phi_in,rhs_out,which_term)
   STOP 'get_rhs_lin2 needs to be benchmarked and updated to 6D g.'
 END IF
 
END SUBROUTINE get_rhs_lin

SUBROUTINE get_rhs_lin1_ae(b_in, v_in, rhs_out_b,rhs_out_v, which_term)

! COMPLEX, INTENT(in) :: g_bounds(0:nkx0-1,0:nky0-1,lkz1:lkz2,lh1:lh2,2)
 COMPLEX, INTENT(in) :: b_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(in) :: v_in(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_b(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
 COMPLEX, INTENT(out) :: rhs_out_v(0:nkx0-1,0:nky0-1,lkz1:lkz2, 0:2)
                                                                                    !INTEGER, INTENT(in) :: which_term

 INTEGER :: i,j,k,l,h,ierr
 !for transpose for left ev's
 INTEGER :: grad1_flag
 INTEGER :: grad2_flag
 COMPLEX :: phi_mod1,phi_mod2,g0_bcast
 COMPLEX :: g_closure
                                                                           !which_term=0:  Everything
                                                                           !which_term=1:  Collisions
                                                                           !which_term=2:  Hyper-collisions
                                                                           !which_term=3:  Phase mixing n-1
                                                                           !which_term=4:  omn term
                                                                           !which_term=5:  Landau damping (n=1)
                                                                           !which_term=6:  Drive
                                                                           !which_term=7:  Extra FLR term (n=0)
                                                                           !which_term=8:  hyp_x,y,z
                                                                           !which_term=9:  Phase mixing n+1
                                                                           !which_term=10:  hyp_conv

                                                                          !  !for transpose for left ev's
                                                                          !  grad1_flag=1
                                                                          !  grad2_flag=2
                                                                          !  phi_mod1=cmplx(1.0,0.0)
                                                                          !  phi_mod2=cmplx(1.0,0.0)
                                                                          ! 
                                                                          !  !IF(verbose.and.mype==0) WRITE(*,*) "get_rhs_lin1", 55
                                                                          !  !IF(verbose.and.mype==0) WRITE(*,*) "nkx0,nky0,nkz0",nkx0,nky0,nkz0
                                                                          !  IF(left_ev) THEN
                                                                          !    !The matrix is symmetric with the exception of these two terms
                                                                          !    grad1_flag=0
                                                                          !    grad2_flag=0
                                                                          !    g0_bcast=v_in(0,0,0)
                                                                          !    CALL MPI_BCAST(g0_bcast,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 
                                                                          !    !WRITE(*,*) "!!!g0_bcast",g0_bcast
                                                                          !    IF((1.ge.lv1).and.(1.le.lv2)) phi_mod1=g_in(0,0,0,1,0,0)/g0_bcast
                                                                          !    IF((2.ge.lv1).and.(2.le.lv2)) phi_mod2=g_in(0,0,0,2,0,0)/g0_bcast
                                                                          !  END IF

 rhs_out_b=cmplx(0.0,0.0)
 rhs_out_v=cmplx(0.0,0.0)

 !IF(verbose.and.mype==0) WRITE(*,*) "get_rhs_lin1", 68
 DO i=0,nkx0-1
   DO j=0,nky0-1
     DO k=lkz1,lkz2
     
     !eqn 14
        rhs_out_b(i,j,k,0) = i_complex*kzgrid(k)*v_in(i,j,k,0) + i_complex*kygrid(j)*b_in(i,j,k,2) -i_complex*kzgrid(k)*b_in(i,j,k,1)
        rhs_out_b(i,j,k,1) = i_complex*kzgrid(k)*v_in(i,j,k,1) + i_complex*kzgrid(k)*b_in(i,j,k,0) -i_complex*kxgrid(i)*b_in(i,j,k,2)
        rhs_out_b(i,j,k,2) = i_complex*kzgrid(k)*v_in(i,j,k,2) + i_complex*kxgrid(i)*b_in(i,j,k,1) -i_complex*kygrid(y)*b_in(i,j,k,0)

      !Eqn 15
        rhs_out_v(i,j,k,0) = i_complex*kzgrid(k)*b_in(i,j,k,0)-i_complex*kxgrid(k)*b_in(i,j,k,2)
        rhs_out_v(i,j,k,1) = i_complex*kzgrid(k)*b_in(i,j,k,1)-i_complex*kygrid(k)*b_in(i,j,k,2)
        rhs_out_v(i,j,k,2) = 0

                                            !       !IF(mype==0.and.verbose) WRITE(*,*) "kxgrid(i),kygrid(j),kzgrid(k)",kxgrid(i),kygrid(j),kzgrid(k)
                                            !       !IF(mype==0.and.verbose) WRITE(*,*) "kxmax_hyp,kymax_hyp,kzmax_hyp",kxmax_hyp,kymax_hyp,kzmax_hyp
                                            !       DO l=lv1,lv2
                                            !        DO h = lh1,lh2
                                            !         IF( (which_term==0.or.which_term==1).and. ( (l.gt.2).or.(.not.em_conserve) )) THEN
                                            !           rhs_out(i,j,k,l,h,0)=-nu*herm_grid(l)*g_in(i,j,k,l,h,0) !Collisions
                                            !         ELSE
                                            !           !IF(verbose) WRITE(*,*) "Before.",mype
                                            !           rhs_out(i,j,k,l,h,0)= cmplx(0.0,0.0)
                                            !           !IF(verbose) WRITE(*,*) "After.",mype
                                            !         END IF
                                            !
                                            !         IF( (which_term==0.or.which_term==2).and. ( (l.gt.2).or.(.not.em_conserve) ) ) &
                                            !            rhs_out(i,j,k,l,h,0)= rhs_out(i,j,k,l,h,0) &
                                            !               -hyp_v*(REAL(herm_grid(l))/REAL(nv0-1))**hypv_order*g_in(i,j,k,l,h,0) !Hyper collisions (Parker)
                                            !         IF( (which_term==0.or.which_term==3).and.  l.ne.0) THEN
                                            !           IF(l==lv1) THEN
                                            !             IF(spatial2d) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(j)*sqrt(herm_grid(l))*g_bounds(i,j,k,h,1)  !Phase mixing 1
                                            !             ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(k)*sqrt(herm_grid(l))*g_bounds(i,j,k,h,1)  !Phase mixing 1
                                            !             END IF
                                            !           ELSE
                                            !             IF(spatial2d) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(j)*sqrt(herm_grid(l))*g_in(i,j,k,l-1,h,0)  !Phase mixing 1
                                            !             ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(k)*sqrt(herm_grid(l))*g_in(i,j,k,l-1,h,0)  !Phase mixing 1
                                            !             END IF
                                            !           END IF
                                            !         END IF
                                            !         IF( (which_term==0.or.which_term==9).and.  l.ne.nv0-1) THEN
                                            !           IF(l==lv2) THEN
                                            !             IF(spatial2d) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(j)*sqrt(herm_grid(l+1))*g_bounds(i,j,k,h,2)
                                            !             ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(k)*sqrt(herm_grid(l+1))*g_bounds(i,j,k,h,2)
                                            !             END IF
                                            !           ELSE
                                            !             IF(spatial2d) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(j)*sqrt(herm_grid(l+1))*g_in(i,j,k,l+1,h,0)  !Phase mixing 2
                                            !             ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(k)*sqrt(herm_grid(l+1))*g_in(i,j,k,l+1,h,0)  !Phase mixing 2
                                            !             END IF
                                            !           END IF
                                            !         END IF
                                            !         IF( (which_term==0.or.which_term==9).and.  l.eq.nv0-1 .and. nuno_closure) THEN
                                            !           diss_max= hyp_v*(REAL(nv0)/REAL(nv0-1))**hypv_order + nu*nv0
                                            !           IF(spatial2d) THEN
                                            !             g_closure=-i_complex*kzgrid(j)*sqrt(real(nv0))*g_in(i,j,k,l,h,0)/diss_max
                                            !             rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(j)*SQRT(REAL(nv0))*g_closure
                                            !           ELSE
                                            !             g_closure=-i_complex*kzgrid(k)*sqrt(real(nv0))*g_in(i,j,k,l,h,0)/diss_max
                                            !             rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)-i_complex*kzgrid(k)*SQRT(REAL(nv0))*g_closure
                                            !           END IF
                                            !         END IF
                                            !         IF( (which_term==0.or.which_term==4).and.  l==0) &
                                            !             rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0) &
                                            !              -i_complex*kygrid(j)*omn*pi**(-0.25)*J0a(i,j)*phi_in(i,j,k)*I0a(i,j,h)  !density gradient term
                                            !         IF( (which_term==0.or.which_term==5).and.  l==grad1_flag) THEN
                                            !             IF(spatial2d) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -i_complex*kzgrid(j)*pi**(-0.25)*J0a(i,j)*phi_in(i,j,k)*phi_mod1*I0a(i,j,h)      !delta_n,1 term
                                            !             ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -i_complex*kzgrid(k)*pi**(-0.25)*J0a(i,j)*phi_in(i,j,k)*phi_mod1*I0a(i,j,h)      !delta_n,1 term
                                            !             END IF
                                            !         END IF
                                            !         IF( (which_term==0.or.which_term==6).and.  l==grad2_flag) &
                                            !             rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -2.0**(-0.5)*pi**(-0.25)*omt*i_complex*kygrid(j)*J0a(i,j)*phi_in(i,j,k)*phi_mod2*I0a(i,j,h) !grad T drive
                                            !         !Extra term for FLR effects and Hankel modification 
                                            !         IF( (which_term==0.or.which_term==7).and. ( l==0.and.((flr_on.and.flr_extra).or.(.not.mu_integrated))  )) &
                                            !             rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0) &
                                            !                  +i_complex*kygrid(j)*omt*pi**(-0.25)*( 0.5*(kxgrid(i)**2+kygrid(j)**2) )*J0a(i,j)*phi_in(i,j,k)*I1a(i,j,h)  !Extra flr term
                                            !
                                            !         IF(which_term==0.or.which_term==10) THEN
                                            !            !!!!!hyp_conv
                                            !            !!!!!hyp_conv
                                            !            IF(k==0.and..not.spatial2d) rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -hyp_nu*hyp_conv*g_in(i,j,k,l,h,0)                            !kz=0 crook term
                                            !            IF(num_k_hyp_conv.gt.0) THEN
                                            !              IF((k.le.num_k_hyp_conv.or.k.ge.nkz0-num_k_hyp_conv).and..not.spatial2d.and.k.ne.0) &
                                            !                     rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                    -hyp_nu*hyp_conv*g_in(i,j,k,l,h,0)                            !crook term
                                            !              IF((j.le.num_k_hyp_conv.or.j.ge.nky0-num_k_hyp_conv).and.(j.ne.0).and.hyp_conv_ky) &
                                            !                     rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                    -hyp_nu*hyp_conv*g_in(i,j,k,l,h,0)                            !crook term
                                            !            END IF
                                            !            !!!!!hyp_conv
                                            !            !!!!!hyp_conv
                                            !
                                            !         END IF
                                            !         IF(which_term==0.or.which_term==8) THEN
                                            !            IF(j==0.and.k==0) rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -nu*hyp_zonal*g_in(i,j,k,l,h,0)                            !Zonal crook term
                                            !           !IF(k==0.and.l==0) rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !            IF(.not.spatial2d) rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                 -hyp_z*g_in(i,j,k,l,h,0)*(kzgrid(k)/kzmax0)**hypz_order !hyperviscous z
                                            !            IF(GyroLES.and..not.Gyroherm) THEN 
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&   
                                            !                  -hyp_x_herm*g_in(i,j,k,l,h,0)*(kxgrid(i)/kxmax0)**hypx_order   !hyperviscous x with herm. dep.
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -hyp_y_herm*g_in(i,j,k,l,h,0)*(kygrid(j)/kymax0)**hypy_order   !hyperviscous y with herm. dep.
                                            !            ELSE IF (GyroLES.and.Gyroherm) THEN
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&   
                                            !                  -hyp_x_herm1(l)*g_in(i,j,k,l,h,0)*(kxgrid(i)/kxmax0)**hypx_order   !hyperviscous x with herm. dep.
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -hyp_y_herm1(l)*g_in(i,j,k,l,h,0)*(kygrid(j)/kymax0)**hypy_order   !hyperviscous y with herm. dep.
                                            !            ELSE
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -hyp_x*g_in(i,j,k,l,h,0)*(kxgrid(i)/kxmax0)**hypx_order   !hyperviscous x
                                            !               rhs_out(i,j,k,l,h,0)=rhs_out(i,j,k,l,h,0)&
                                            !                  -hyp_y*g_in(i,j,k,l,h,0)*(kygrid(j)/kymax0)**hypy_order   !hyperviscous y
                                            !            END IF 
                                            !
                                            !         END IF
                                            !        END DO
                                            !       END DO
     END DO
   END DO
 END DO 

END SUBROUTINE get_rhs_lin1_ae


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_rhs_lin2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Note!!!!!!
!!  The following routine is a bit faster, but it needs to be benchmarked!
!!  I think the hyp_x/y/z is messing things up?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE get_rhs_lin2(g_in,g_bounds,phi_in,rhs_out,which_term)
!
! COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2)
! COMPLEX, INTENT(in) :: g_bounds(0:nkx0-1,0:nky0-1,0:nkz0-1,2)
! COMPLEX, INTENT(in) :: phi_in(0:nkx0-1,0:nky0-1,0:nkz0-1)
! COMPLEX, INTENT(out) :: rhs_out(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2)
! INTEGER, INTENT(in) :: which_term
!
! INTEGER :: i,j,k,l,ierr
! !for transpose for left ev's
! INTEGER :: grad1_flag
! INTEGER :: grad2_flag
! COMPLEX :: phi_mod1,phi_mod2,g0_bcast
! !for transpose for left ev's
!
!
! grad1_flag=1
! grad2_flag=2
! phi_mod1=cmplx(1.0,0.0)
! phi_mod2=cmplx(1.0,0.0)
!
! IF(left_ev) THEN
!   !The matrix is symmetric with the exception of these two terms
!   grad1_flag=0
!   grad2_flag=0
!   g0_bcast=g_in(0,0,0,lv1)
!   CALL MPI_BCAST(g0_bcast,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 
!   !WRITE(*,*) "!!!g0_bcast",g0_bcast
!   IF((1.ge.lv1).and.(1.le.lv2)) phi_mod1=g_in(0,0,0,1)/g0_bcast
!   IF((2.ge.lv1).and.(2.le.lv2)) phi_mod2=g_in(0,0,0,2)/g0_bcast
! END IF
!
!   DO l=lv1,lv2
!     IF( (l.gt.2).or.(.not.em_conserve) ) THEN
!       rhs_out(:,:,:,l)= (-nu*herm_grid(l)-hyp_v*(REAL(herm_grid(l))/REAL(nv0-1))**hypv_order)*g_in(:,:,:,l) !Collisions
!     ELSE
!       rhs_out(:,:,:,l)= cmplx(0.0,0.0)
!     END IF
!     DO k=0,nkz0-1
!       IF( l.ne.0) THEN
!         IF(l==lv1) THEN
!          rhs_out(:,:,k,l)=rhs_out(:,:,k,l)-i_complex*kzgrid(k)*sqrt(herm_grid(l))*g_bounds(:,:,k,1)  !Phase mixing 1
!         ELSE
!          rhs_out(:,:,k,l)=rhs_out(:,:,k,l)-i_complex*kzgrid(k)*sqrt(herm_grid(l))*g_in(:,:,k,l-1)  !Phase mixing 1
!         END IF
!       END IF
!       IF( l.ne.nv0-1) THEN
!         IF(l==lv2) THEN
!           rhs_out(:,:,k,l)=rhs_out(:,:,k,l)-i_complex*kzgrid(k)*sqrt(herm_grid(l+1))*g_bounds(:,:,k,2)  !Phase mixing 2
!         ELSE
!           rhs_out(:,:,k,l)=rhs_out(:,:,k,l)-i_complex*kzgrid(k)*sqrt(herm_grid(l+1))*g_in(:,:,k,l+1)  !Phase mixing 2
!         END IF
!       END IF
!     END DO
!   END DO
!
!   DO i=0,nkx0-1
!     IF(hyp_zonal.gt.epsilon(1.0)) rhs_out(i,0,0,:)=rhs_out(i,0,0,:)&
!              -nu*hyp_zonal*g_in(i,0,0,:)                            !Zonal crook term
!     DO j=0,nky0-1
!
!       IF(mype==0)  rhs_out(i,j,:,0)=rhs_out(i,j,:,0) &
!          -i_complex*kygrid(j)*omn*pi**(-0.25)*J0a(i,j)*phi_in(i,j,:)  !density gradient term
!       IF(lv1.le.grad2_flag.and.lv2.ge.grad2_flag) rhs_out(i,j,:,grad2_flag)=rhs_out(i,j,:,grad2_flag)&
!              -2.0**(-0.5)*pi**(-0.25)*omt*i_complex*kygrid(j)*J0a(i,j)*phi_in(i,j,:)*phi_mod2 !grad T drive
!     !Extra term for FLR effects 
!       IF( flr_on.and.flr_extra.and.mype==0 ) &
!         rhs_out(i,j,:,0)=rhs_out(i,j,:,0) &
!              +i_complex*kygrid(j)*omt*pi**(-0.25)*( 0.5*(kxgrid(i)**2+kygrid(j)**2) )*J0a(i,j)*phi_in(i,j,:)  !Extra flr term
!
!       IF(hyp_x.gt.epsilon(1.0)) rhs_out(i,j,:,:)=rhs_out(i,j,:,:)&
!              -hyp_x*g_in(i,j,:,:)*(kxgrid(i)/kxmax0)**hypx_order   !hyperviscous x
!       !IF(mype==0) WRITE(*,'(4es16.8)') -hyp_x*(kxgrid(i)/kxmax)**hypx_order,hyp_x,kxgrid(i), kxmax
!       !IF(mype==0) WRITE(*,*) nkx0
!       IF(hyp_y.gt.epsilon(1.0)) rhs_out(i,j,:,:)=rhs_out(i,j,:,:)&
!              -hyp_y*g_in(i,j,:,:)*(kygrid(j)/kymax0)**hypy_order   !hyperviscous y
!       !IF(mype==0) WRITE(*,'(4es16.8)') -hyp_y*(kygrid(j)/kymax)**hypy_order,hyp_y,kygrid(j),kymax
!       !IF(mype==0) WRITE(*,*) nky0
!
!       DO k=0,nkz0-1
!         IF(lv1.le.grad1_flag.and.lv2.ge.grad1_flag)  rhs_out(i,j,k,grad1_flag)=rhs_out(i,j,k,grad1_flag)&
!              -i_complex*kzgrid(k)*pi**(-0.25)*J0a(i,j)*phi_in(i,j,k)*phi_mod1      !delta_n,1 term
!         IF(hyp_z.gt.epsilon(1.0)) rhs_out(i,j,k,:)=rhs_out(i,j,k,:)&
!              -hyp_z*g_in(i,j,k,:)*(kzgrid(k)/kzmax0)**hypz_order !hyperviscous z 
!
!       END DO
!
!     END DO
!   END DO
!
!   IF(mype==0.and.hyp_conv.ne.0.0) rhs_out(:,:,0,0)=rhs_out(:,:,0,0)-hyp_nu*hyp_conv*g_in(:,:,0,0) !kz=0, l=0 crook term
!
!END SUBROUTINE get_rhs_lin2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                get_v_boundaries                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE get_v_boundaries(g_in,g_wb)
!  USE par_mod
!  USE mpi
!
!  COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,0:nkz0-1,lv1:lv2)
!  COMPLEX, INTENT(out) :: g_wb(0:nkx0-1,0:nky0-1,0:nkz0-1,lbv:ubv)
!  INTEGER :: p,send_proc,recv_proc,ierr
!  !COMPLEX :: send_arr(0:nkx0-1,0:nky0-1,0:nkz0-1)
!  !COMPLEX :: recv_arr(0:nkx0-1,0:nky0-1,0:nkz0-1)
!  INTEGER :: stat(MPI_STATUS_SIZE)
!
!  !WRITE(*,*) "get_v_boundaries,mype",mype
!  g_wb(:,:,:,lv1:lv2)=g_in(:,:,:,lv1:lv2)
!  
!  !fill lower boundaries 
!  DO p=1,np_herm-1
!   !WRITE(*,*) "lower boundaries,mype,p",mype,p
!   send_proc=p-1
!   recv_proc=p
!
!!   IF(mype==send_proc) send_arr(:,:,:)=g_in(:,:,:,lv2)
!   IF(mype==send_proc) CALL MPI_Send(g_wb(0,0,0,lv2), nkx0*nky0*nkz0, &
!                        MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!   IF(mype==recv_proc) CALL MPI_Recv(g_wb(0,0,0,lbv), nkx0*nky0*nkz0, &
!                        MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!!   IF(mype==recv_proc) g_wb(:,:,:,lbv)=recv_arr(:,:,:)
!
!  END DO
!
!  !WRITE(*,*) "Before barrier",mype
!  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  !WRITE(*,*) "After barrier",mype
!
!  !fill upper boundaries 
!  DO p=0,np_herm-2
!    !WRITE(*,*) "upper boundaries,mype,p",mype,p
!   !IF(mype==p+1) send_arr(:,:,:)=g_in(:,:,:,lv1)
!   send_proc=p+1
!   recv_proc=p
!
!   !IF(mype==send_proc) send_arr(:,:,:)=g_in(:,:,:,lv1)
!   IF(mype==send_proc) CALL MPI_Send(g_wb(0,0,0,lv1), nkx0*nky0*nkz0, &
!                       MPI_DOUBLE_COMPLEX, recv_proc, p, MPI_COMM_WORLD, ierr)  
!   IF(mype==recv_proc) CALL MPI_Recv(g_wb(0,0,0,ubv), nkx0*nky0*nkz0, &
!                       MPI_DOUBLE_COMPLEX, send_proc, p, MPI_COMM_WORLD, stat, ierr )  
!   !IF(mype==recv_proc) g_wb(:,:,:,ubv)=recv_arr(:,:,:)
!
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!END SUBROUTINE get_v_boundaries


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               get_v_boundaries2                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE get_v_boundaries2(g_in,g_bounds)
!   USE par_mod
!   USE mpi
! 
!   COMPLEX, INTENT(in) :: g_in(0:nkx0-1,0:nky0-1,lkz1:lkz2,lv1:lv2,lh1:lh2,ls1:ls2)
!   COMPLEX, INTENT(out) :: g_bounds(0:nkx0-1,0:nky0-1,lkz1:lkz2,lh1:lh2,2)
!   INTEGER :: p,send_proc,recv_proc,ierr
!   COMPLEX :: send_arr(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!   COMPLEX :: recv_arr(0:nkx0-1,0:nky0-1,lkz1:lkz2)
!   INTEGER :: stat(MPI_STATUS_SIZE)
!   INTEGER :: h
! 
!   IF(verbose) write(*,*) "In get_v_boundaries2",mype
!   !fill lower boundaries
! 
!   DO h=lh1,lh2
!     IF(mype_herm==0) CALL MPI_Send(g_in(0,0,0,lv2,h,0), nkx0*nky0*lkz0, &
!                         MPI_DOUBLE_COMPLEX, 1, 0, mpi_comm_herm, ierr)  
!     IF(mype_herm.ne.0.and.mype_herm.ne.np_herm-1) THEN
!         CALL MPI_SENDRECV(g_in(0,0,0,lv2,h,0), nkx0*nky0*lkz0, MPI_DOUBLE_COMPLEX,mype_herm+1,0,&
!                     g_bounds(0,0,0,h,1), nkx0*nky0*lkz0, MPI_DOUBLE_COMPLEX,mype_herm-1,0,&
!                     mpi_comm_herm, stat, ierr )
!     END IF
!     IF(mype_herm==np_herm-1) CALL MPI_Recv(g_bounds(0,0,0,h,1), nkx0*nky0*lkz0, &
!                         MPI_DOUBLE_COMPLEX, np_herm-2, 0, mpi_comm_herm, stat, ierr )  
! 
!     !fill upper boundaries
!     IF(mype_herm==np_herm-1) CALL MPI_Send(g_in(0,0,0,lv1,h,0), nkx0*nky0*lkz0, &
!                        MPI_DOUBLE_COMPLEX, np_herm-2, 0, mpi_comm_herm, ierr)  
!     IF(mype_herm.ne.0.and.mype_herm.ne.np_herm-1) THEN
!         CALL MPI_SENDRECV(g_in(0,0,0,lv1,h,0), nkx0*nky0*lkz0, MPI_DOUBLE_COMPLEX,mype_herm-1,0,&
!                     g_bounds(0,0,0,h,2), nkx0*nky0*lkz0, MPI_DOUBLE_COMPLEX,mype_herm+1,0,&
!                     mpi_comm_herm, stat, ierr )
!     END IF
!     IF(mype_herm==0) CALL MPI_Recv(g_bounds(0,0,0,h,2), nkx0*nky0*lkz0, &
!                        MPI_DOUBLE_COMPLEX, 1, 0, mpi_comm_herm, stat, ierr )  
!   END DO
!   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
! END SUBROUTINE get_v_boundaries2

END MODULE linear_rhs
