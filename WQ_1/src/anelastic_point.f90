module anelastic_point

  ! Per-point memory-variable (GSLS) updates for the newer attenuation models
  ! (anelastic-Qf, frequency-Q-4M/8M, constant-Q-4M/8M).
  !
  ! These routines are shared by the interior stencil routines in
  ! JU_xJU_yJU_z6 and the near-boundary/dispatch routines in RHS_Interior
  ! (they used to live in RHS_Interior, but JU_xJU_yJU_z6 cannot use that
  ! module without creating a circular dependency).
  !
  ! Only the non-PML variants live here. The PML variants and the
  ! point-in-PML routing remain in RHS_Interior: the interior loops in
  ! JU_xJU_yJU_z6 exclude PML zones, so the *_any routers below are safe to
  ! call from there.

  use common, only : wp

  implicit none

contains

     subroutine apply_anelastic_Qf_point(M, x, y, z, Ux, Uy, Uz, DU)

          ! anelastic-Qf: frequency-dependent Q with 4 mechanisms
          ! Operates on the _Qf memory variables (eta4Qf..eta9Qf, tau_Qf, weight_Qf,
          ! Qs_inv_Qf, Qp_inv_Qf).
          ! Same kernel structure as apply_anelastic_Q_point.

          use common, only : wp
          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          integer :: i
          real(kind = wp) :: tr

          ! Subtract summed memory variables from stress RHS
          DU(4) = DU(4) - (M%eta4Qf(x,y,z,1) + M%eta4Qf(x,y,z,2) + M%eta4Qf(x,y,z,3) + M%eta4Qf(x,y,z,4))
          DU(5) = DU(5) - (M%eta5Qf(x,y,z,1) + M%eta5Qf(x,y,z,2) + M%eta5Qf(x,y,z,3) + M%eta5Qf(x,y,z,4))
          DU(6) = DU(6) - (M%eta6Qf(x,y,z,1) + M%eta6Qf(x,y,z,2) + M%eta6Qf(x,y,z,3) + M%eta6Qf(x,y,z,4))
          DU(7) = DU(7) - (M%eta7Qf(x,y,z,1) + M%eta7Qf(x,y,z,2) + M%eta7Qf(x,y,z,3) + M%eta7Qf(x,y,z,4))
          DU(8) = DU(8) - (M%eta8Qf(x,y,z,1) + M%eta8Qf(x,y,z,2) + M%eta8Qf(x,y,z,3) + M%eta8Qf(x,y,z,4))
          DU(9) = DU(9) - (M%eta9Qf(x,y,z,1) + M%eta9Qf(x,y,z,2) + M%eta9Qf(x,y,z,3) + M%eta9Qf(x,y,z,4))

          ! Compute volumetric strain rate
          tr = Ux(1) + Uy(2) + Uz(3)

          ! Update memory variables for each mechanism (l=1..4)
          do i = 1, 4
               ! σ_xx (eta4): dη4/dt
               M%Deta4Qf(x,y,z,i) = M%Deta4Qf(x,y,z,i) + ( &
                    ( (M%weight_Qf(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z))*Ux(1) &
                    + (M%weight_Qf(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z) ) * tr ) ) &
                    - M%eta4Qf(x,y,z,i) ) / M%tau_Qf(i)

               ! σ_yy (eta5): dη5/dt
               M%Deta5Qf(x,y,z,i) = M%Deta5Qf(x,y,z,i) + ( &
                    ( (M%weight_Qf(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z))*Uy(2) &
                    + (M%weight_Qf(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z) ) * tr ) ) &
                    - M%eta5Qf(x,y,z,i) ) / M%tau_Qf(i)

               ! σ_zz (eta6): dη6/dt
               M%Deta6Qf(x,y,z,i) = M%Deta6Qf(x,y,z,i) + ( &
                    ( (M%weight_Qf(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z))*Uz(3) &
                    + (M%weight_Qf(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z) ) * tr ) ) &
                    - M%eta6Qf(x,y,z,i) ) / M%tau_Qf(i)

               ! σ_xy (eta7): dη7/dt
               M%Deta7Qf(x,y,z,i) = M%Deta7Qf(x,y,z,i) + ( &
                    (M%weight_Qf(i)*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z)*(Uy(1) + Ux(2)) - M%eta7Qf(x,y,z,i)) / M%tau_Qf(i) )
               
               ! σ_xz (eta8): dη8/dt
               M%Deta8Qf(x,y,z,i) = M%Deta8Qf(x,y,z,i) + ( &
                    (M%weight_Qf(i)*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z)*(Uz(1) + Ux(3)) - M%eta8Qf(x,y,z,i)) / M%tau_Qf(i) )
               
               ! σ_yz (eta9): dη9/dt
               M%Deta9Qf(x,y,z,i) = M%Deta9Qf(x,y,z,i) + ( &
                    (M%weight_Qf(i)*M%M(x,y,z,2)*M%Qs_inv_Qf(x,y,z)*(Uz(2) + Uy(3)) - M%eta9Qf(x,y,z,i)) / M%tau_Qf(i) )
          end do

     end subroutine apply_anelastic_Qf_point

     subroutine apply_const_Q_4M_point(M, x, y, z, Ux, Uy, Uz, DU)

          ! constant-Q-4M: 4-mechanism constant-Q with user-configurable [fmin, fmax]
          ! Operates on the _4M memory variables (eta4_4M..eta9_4M, tau_const_Q_4M, 
          ! weight_const_Q_4M, Qs_inv_const_Q_4M, Qp_inv_const_Q_4M).
          ! Same kernel structure as apply_anelastic_Q_point (N=4).

          use common, only : wp
          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          integer :: i
          real(kind = wp) :: tr

          ! Subtract summed memory variables from stress RHS
          DU(4) = DU(4) - (M%eta4_4M(x,y,z,1) + M%eta4_4M(x,y,z,2) + M%eta4_4M(x,y,z,3) + M%eta4_4M(x,y,z,4))
          DU(5) = DU(5) - (M%eta5_4M(x,y,z,1) + M%eta5_4M(x,y,z,2) + M%eta5_4M(x,y,z,3) + M%eta5_4M(x,y,z,4))
          DU(6) = DU(6) - (M%eta6_4M(x,y,z,1) + M%eta6_4M(x,y,z,2) + M%eta6_4M(x,y,z,3) + M%eta6_4M(x,y,z,4))
          DU(7) = DU(7) - (M%eta7_4M(x,y,z,1) + M%eta7_4M(x,y,z,2) + M%eta7_4M(x,y,z,3) + M%eta7_4M(x,y,z,4))
          DU(8) = DU(8) - (M%eta8_4M(x,y,z,1) + M%eta8_4M(x,y,z,2) + M%eta8_4M(x,y,z,3) + M%eta8_4M(x,y,z,4))
          DU(9) = DU(9) - (M%eta9_4M(x,y,z,1) + M%eta9_4M(x,y,z,2) + M%eta9_4M(x,y,z,3) + M%eta9_4M(x,y,z,4))

          ! Compute volumetric strain rate
          tr = Ux(1) + Uy(2) + Uz(3)

          ! Update memory variables for each mechanism (i=1..4)
          do i = 1, 4
               ! σ_xx (eta4): dη4/dt
               M%Deta4_4M(x,y,z,i) = M%Deta4_4M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_4M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z))*Ux(1) &
                    + (M%weight_const_Q_4M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_4M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z) ) * tr ) ) &
                    - M%eta4_4M(x,y,z,i) ) / M%tau_const_Q_4M(i)

               ! σ_yy (eta5): dη5/dt
               M%Deta5_4M(x,y,z,i) = M%Deta5_4M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_4M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z))*Uy(2) &
                    + (M%weight_const_Q_4M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_4M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z) ) * tr ) ) &
                    - M%eta5_4M(x,y,z,i) ) / M%tau_const_Q_4M(i)

               ! σ_zz (eta6): dη6/dt
               M%Deta6_4M(x,y,z,i) = M%Deta6_4M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_4M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z))*Uz(3) &
                    + (M%weight_const_Q_4M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_4M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z) ) * tr ) ) &
                    - M%eta6_4M(x,y,z,i) ) / M%tau_const_Q_4M(i)

               ! σ_xy (eta7): dη7/dt
               M%Deta7_4M(x,y,z,i) = M%Deta7_4M(x,y,z,i) + ( &
                    (M%weight_const_Q_4M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z)*(Uy(1) + Ux(2)) - M%eta7_4M(x,y,z,i)) / M%tau_const_Q_4M(i) )
               ! σ_xz (eta8): dη8/dt
               M%Deta8_4M(x,y,z,i) = M%Deta8_4M(x,y,z,i) + ( &
                    (M%weight_const_Q_4M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z)*(Uz(1) + Ux(3)) - M%eta8_4M(x,y,z,i)) / M%tau_const_Q_4M(i) )
               ! σ_yz (eta9): dη9/dt
               M%Deta9_4M(x,y,z,i) = M%Deta9_4M(x,y,z,i) + ( &
                    (M%weight_const_Q_4M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_4M(x,y,z)*(Uz(2) + Uy(3)) - M%eta9_4M(x,y,z,i)) / M%tau_const_Q_4M(i) )
          end do

     end subroutine apply_const_Q_4M_point

     subroutine apply_const_Q_8M_point(M, x, y, z, Ux, Uy, Uz, DU)

          use common, only : wp
          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          integer :: i
          real(kind = wp) :: tr

          DU(4) = DU(4) - sum(M%eta4_8M(x,y,z,:))
          DU(5) = DU(5) - sum(M%eta5_8M(x,y,z,:))
          DU(6) = DU(6) - sum(M%eta6_8M(x,y,z,:))
          DU(7) = DU(7) - sum(M%eta7_8M(x,y,z,:))
          DU(8) = DU(8) - sum(M%eta8_8M(x,y,z,:))
          DU(9) = DU(9) - sum(M%eta9_8M(x,y,z,:))

          tr = Ux(1) + Uy(2) + Uz(3)

          do i = 1, 8
               M%Deta4_8M(x,y,z,i) = M%Deta4_8M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_8M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z))*Ux(1) &
                    + (M%weight_const_Q_8M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_8M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z) ) * tr ) ) &
                    - M%eta4_8M(x,y,z,i) ) / M%tau_const_Q_8M(i)

               M%Deta5_8M(x,y,z,i) = M%Deta5_8M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_8M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z))*Uy(2) &
                    + (M%weight_const_Q_8M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_8M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z) ) * tr ) ) &
                    - M%eta5_8M(x,y,z,i) ) / M%tau_const_Q_8M(i)

               M%Deta6_8M(x,y,z,i) = M%Deta6_8M(x,y,z,i) + ( &
                    ( (M%weight_const_Q_8M(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z))*Uz(3) &
                    + (M%weight_const_Q_8M(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_const_Q_8M(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z) ) * tr ) ) &
                    - M%eta6_8M(x,y,z,i) ) / M%tau_const_Q_8M(i)

               M%Deta7_8M(x,y,z,i) = M%Deta7_8M(x,y,z,i) + ( &
                    (M%weight_const_Q_8M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z)*(Uy(1) + Ux(2)) - M%eta7_8M(x,y,z,i)) / M%tau_const_Q_8M(i) )
               M%Deta8_8M(x,y,z,i) = M%Deta8_8M(x,y,z,i) + ( &
                    (M%weight_const_Q_8M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z)*(Uz(1) + Ux(3)) - M%eta8_8M(x,y,z,i)) / M%tau_const_Q_8M(i) )
               M%Deta9_8M(x,y,z,i) = M%Deta9_8M(x,y,z,i) + ( &
                    (M%weight_const_Q_8M(i)*M%M(x,y,z,2)*M%Qs_inv_const_Q_8M(x,y,z)*(Uz(2) + Uy(3)) - M%eta9_8M(x,y,z,i)) / M%tau_const_Q_8M(i) )
          end do

     end subroutine apply_const_Q_8M_point

     subroutine apply_anelastic_Qf8_point(M, x, y, z, Ux, Uy, Uz, DU)

          use common, only : wp
          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          integer :: i
          real(kind = wp) :: tr

          DU(4) = DU(4) - sum(M%eta4Qf8(x,y,z,:))
          DU(5) = DU(5) - sum(M%eta5Qf8(x,y,z,:))
          DU(6) = DU(6) - sum(M%eta6Qf8(x,y,z,:))
          DU(7) = DU(7) - sum(M%eta7Qf8(x,y,z,:))
          DU(8) = DU(8) - sum(M%eta8Qf8(x,y,z,:))
          DU(9) = DU(9) - sum(M%eta9Qf8(x,y,z,:))

          tr = Ux(1) + Uy(2) + Uz(3)

          do i = 1, 8
               M%Deta4Qf8(x,y,z,i) = M%Deta4Qf8(x,y,z,i) + ( &
                    ( (M%weight_Qf8(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z))*Ux(1) &
                    + (M%weight_Qf8(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf8(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z) ) * tr ) ) &
                    - M%eta4Qf8(x,y,z,i) ) / M%tau_Qf8(i)
               M%Deta5Qf8(x,y,z,i) = M%Deta5Qf8(x,y,z,i) + ( &
                    ( (M%weight_Qf8(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z))*Uy(2) &
                    + (M%weight_Qf8(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf8(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z) ) * tr ) ) &
                    - M%eta5Qf8(x,y,z,i) ) / M%tau_Qf8(i)
               M%Deta6Qf8(x,y,z,i) = M%Deta6Qf8(x,y,z,i) + ( &
                    ( (M%weight_Qf8(i)*2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z))*Uz(3) &
                    + (M%weight_Qf8(i) * ( (M%M(x,y,z,1)+2.0_wp*M%M(x,y,z,2))*M%Qp_inv_Qf8(x,y,z) &
                                                     - 2.0_wp*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z) ) * tr ) ) &
                    - M%eta6Qf8(x,y,z,i) ) / M%tau_Qf8(i)
               M%Deta7Qf8(x,y,z,i) = M%Deta7Qf8(x,y,z,i) + ( &
                    (M%weight_Qf8(i)*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z)*(Uy(1) + Ux(2)) - M%eta7Qf8(x,y,z,i)) / M%tau_Qf8(i) )
               M%Deta8Qf8(x,y,z,i) = M%Deta8Qf8(x,y,z,i) + ( &
                    (M%weight_Qf8(i)*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z)*(Uz(1) + Ux(3)) - M%eta8Qf8(x,y,z,i)) / M%tau_Qf8(i) )
               M%Deta9Qf8(x,y,z,i) = M%Deta9Qf8(x,y,z,i) + ( &
                    (M%weight_Qf8(i)*M%M(x,y,z,2)*M%Qs_inv_Qf8(x,y,z)*(Uz(2) + Uy(3)) - M%eta9Qf8(x,y,z,i)) / M%tau_Qf8(i) )
          end do

     end subroutine apply_anelastic_Qf8_point


     subroutine apply_anelastic_Qf_point_any(M, x, y, z, Ux, Uy, Uz, DU)

          ! Frequency-dependent Q at a non-PML point: routes to the
          ! 8-mechanism variant when its memory variables are allocated
          ! (frequency-Q-8M), otherwise the 4-mechanism variant
          ! (anelastic-Qf, frequency-Q-4M).

          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          if (allocated(M%eta4Qf8)) then
               call apply_anelastic_Qf8_point(M, x, y, z, Ux, Uy, Uz, DU)
          else
               call apply_anelastic_Qf_point(M, x, y, z, Ux, Uy, Uz, DU)
          end if

     end subroutine apply_anelastic_Qf_point_any


     subroutine apply_const_Q_point_any(M, x, y, z, Ux, Uy, Uz, DU)

          ! Constant Q at a non-PML point: routes to the 8-mechanism variant
          ! when its memory variables are allocated (constant-Q-8M),
          ! otherwise the 4-mechanism variant (constant-Q-4M).

          use datatypes, only : block_material

          implicit none

          type(block_material), intent(inout) :: M
          integer, intent(in) :: x, y, z
          real(kind = wp), intent(in) :: Ux(:), Uy(:), Uz(:)
          real(kind = wp), intent(inout) :: DU(:)

          if (allocated(M%eta4_8M)) then
               call apply_const_Q_8M_point(M, x, y, z, Ux, Uy, Uz, DU)
          else
               call apply_const_Q_4M_point(M, x, y, z, Ux, Uy, Uz, DU)
          end if

     end subroutine apply_const_Q_point_any

end module anelastic_point
