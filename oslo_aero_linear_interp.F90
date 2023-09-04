module  oslo_aero_linear_interp

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: lininterpol3dim
  public :: lininterpol4dim
  public :: lininterpol5dim

! ==========================================================
contains
! ==========================================================

  subroutine lininterpol3dim (d2mx, dxm1, invd, opt3d, optout1, optout2)

    ! arguments
    real(r8), intent(in)  :: opt3d(2,2,2)
    real(r8), intent(in)  :: d2mx(3)
    real(r8), intent(in)  :: dxm1(3)
    real(r8), intent(in)  :: invd(3)
    real(r8), intent(out) :: optout1
    real(r8), intent(out) :: optout2
    !
    ! local variables
    real(r8) :: opt2d(2,2)
    !------------------------------------

    ! interpolation in the third dimension (except invd(3) factor)
    opt2d(1,1) = d2mx(3)*opt3d(1,1,1) + dxm1(3)*opt3d(1,1,2)
    opt2d(1,2) = d2mx(3)*opt3d(1,2,1) + dxm1(3)*opt3d(1,2,2)
    opt2d(2,1) = d2mx(3)*opt3d(2,1,1) + dxm1(3)*opt3d(2,1,2)
    opt2d(2,2) = d2mx(3)*opt3d(2,2,1) + dxm1(3)*opt3d(2,2,2)

    ! interpolation in the (third and) second dimension
    optout1 = (d2mx(2)*opt2d(1,1) + dxm1(2)*opt2d(1,2))*invd(3)*invd(2)
    optout2 = (d2mx(2)*opt2d(2,1) + dxm1(2)*opt2d(2,2))*invd(3)*invd(2)

  end subroutine lininterpol3dim

  ! ==========================================================
  subroutine lininterpol4dim (d2mx, dxm1, invd, opt4d, optout1, optout2)

    ! arguments
    real(r8), intent(in)  :: opt4d(2,2,2,2)
    real(r8), intent(in)  :: d2mx(4)
    real(r8), intent(in)  :: dxm1(4)
    real(r8), intent(in)  :: invd(4)
    real(r8), intent(out) :: optout1
    real(r8), intent(out) :: optout2
    !
    ! local variables
    real(r8) :: opt3d(2,2,2), opt2d(2,2)
    !------------------------------------

    ! interpolation in the fourth dimension (except invd(4) factor)
    opt3d(1,1,1) = d2mx(4)*opt4d(1,1,1,1) + dxm1(4)*opt4d(1,1,1,2)
    opt3d(1,1,2) = d2mx(4)*opt4d(1,1,2,1) + dxm1(4)*opt4d(1,1,2,2)
    opt3d(1,2,1) = d2mx(4)*opt4d(1,2,1,1) + dxm1(4)*opt4d(1,2,1,2)
    opt3d(1,2,2) = d2mx(4)*opt4d(1,2,2,1) + dxm1(4)*opt4d(1,2,2,2)
    opt3d(2,1,1) = d2mx(4)*opt4d(2,1,1,1) + dxm1(4)*opt4d(2,1,1,2)
    opt3d(2,1,2) = d2mx(4)*opt4d(2,1,2,1) + dxm1(4)*opt4d(2,1,2,2)
    opt3d(2,2,1) = d2mx(4)*opt4d(2,2,1,1) + dxm1(4)*opt4d(2,2,1,2)
    opt3d(2,2,2) = d2mx(4)*opt4d(2,2,2,1) + dxm1(4)*opt4d(2,2,2,2)

    ! interpolation in the third dimension (except invd(3) factor)
    opt2d(1,1) = d2mx(3)*opt3d(1,1,1) + dxm1(3)*opt3d(1,1,2)
    opt2d(1,2) = d2mx(3)*opt3d(1,2,1) + dxm1(3)*opt3d(1,2,2)
    opt2d(2,1) = d2mx(3)*opt3d(2,1,1) + dxm1(3)*opt3d(2,1,2)
    opt2d(2,2) = d2mx(3)*opt3d(2,2,1) + dxm1(3)*opt3d(2,2,2)

    ! interpolation in the (fourth, third and) second dimension
    optout1 = (d2mx(2)*opt2d(1,1) + dxm1(2)*opt2d(1,2))*invd(4)*invd(3)*invd(2)
    optout2 = (d2mx(2)*opt2d(2,1) + dxm1(2)*opt2d(2,2))*invd(4)*invd(3)*invd(2)

  end subroutine lininterpol4dim

  ! ==========================================================
  subroutine lininterpol5dim (d2mx, dxm1, invd, opt5d, optout1, optout2)

    ! arguments
    real(r8), intent(in)  :: opt5d(2,2,2,2,2)
    real(r8), intent(in)  :: d2mx(5)
    real(r8), intent(in)  :: dxm1(5)
    real(r8), intent(in)  :: invd(5)
    real(r8), intent(out) :: optout1
    real(r8), intent(out) :: optout2

    ! local variables
    real(r8) :: opt4d(2,2,2,2), opt3d(2,2,2), opt2d(2,2)
    !------------------------------------

    ! interpolation in the fifth dimension (except invd(5) factor)
    opt4d(1,1,1,1) = d2mx(5)*opt5d(1,1,1,1,1) + dxm1(5)*opt5d(1,1,1,1,2)
    opt4d(1,1,1,2) = d2mx(5)*opt5d(1,1,1,2,1) + dxm1(5)*opt5d(1,1,1,2,2)
    opt4d(1,1,2,1) = d2mx(5)*opt5d(1,1,2,1,1) + dxm1(5)*opt5d(1,1,2,1,2)
    opt4d(1,1,2,2) = d2mx(5)*opt5d(1,1,2,2,1) + dxm1(5)*opt5d(1,1,2,2,2)
    opt4d(1,2,1,1) = d2mx(5)*opt5d(1,2,1,1,1) + dxm1(5)*opt5d(1,2,1,1,2)
    opt4d(1,2,1,2) = d2mx(5)*opt5d(1,2,1,2,1) + dxm1(5)*opt5d(1,2,1,2,2)
    opt4d(1,2,2,1) = d2mx(5)*opt5d(1,2,2,1,1) + dxm1(5)*opt5d(1,2,2,1,2)
    opt4d(1,2,2,2) = d2mx(5)*opt5d(1,2,2,2,1) + dxm1(5)*opt5d(1,2,2,2,2)
    opt4d(2,1,1,1) = d2mx(5)*opt5d(2,1,1,1,1) + dxm1(5)*opt5d(2,1,1,1,2)
    opt4d(2,1,1,2) = d2mx(5)*opt5d(2,1,1,2,1) + dxm1(5)*opt5d(2,1,1,2,2)
    opt4d(2,1,2,1) = d2mx(5)*opt5d(2,1,2,1,1) + dxm1(5)*opt5d(2,1,2,1,2)
    opt4d(2,1,2,2) = d2mx(5)*opt5d(2,1,2,2,1) + dxm1(5)*opt5d(2,1,2,2,2)
    opt4d(2,2,1,1) = d2mx(5)*opt5d(2,2,1,1,1) + dxm1(5)*opt5d(2,2,1,1,2)
    opt4d(2,2,1,2) = d2mx(5)*opt5d(2,2,1,2,1) + dxm1(5)*opt5d(2,2,1,2,2)
    opt4d(2,2,2,1) = d2mx(5)*opt5d(2,2,2,1,1) + dxm1(5)*opt5d(2,2,2,1,2)
    opt4d(2,2,2,2) = d2mx(5)*opt5d(2,2,2,2,1) + dxm1(5)*opt5d(2,2,2,2,2)

    ! interpolation in the fourth dimension (except invd(4) factor)
    opt3d(1,1,1) = d2mx(4)*opt4d(1,1,1,1) + dxm1(4)*opt4d(1,1,1,2)
    opt3d(1,1,2) = d2mx(4)*opt4d(1,1,2,1) + dxm1(4)*opt4d(1,1,2,2)
    opt3d(1,2,1) = d2mx(4)*opt4d(1,2,1,1) + dxm1(4)*opt4d(1,2,1,2)
    opt3d(1,2,2) = d2mx(4)*opt4d(1,2,2,1) + dxm1(4)*opt4d(1,2,2,2)
    opt3d(2,1,1) = d2mx(4)*opt4d(2,1,1,1) + dxm1(4)*opt4d(2,1,1,2)
    opt3d(2,1,2) = d2mx(4)*opt4d(2,1,2,1) + dxm1(4)*opt4d(2,1,2,2)
    opt3d(2,2,1) = d2mx(4)*opt4d(2,2,1,1) + dxm1(4)*opt4d(2,2,1,2)
    opt3d(2,2,2) = d2mx(4)*opt4d(2,2,2,1) + dxm1(4)*opt4d(2,2,2,2)

    ! interpolation in the third dimension (except invd(3) factor)
    opt2d(1,1) = d2mx(3)*opt3d(1,1,1) + dxm1(3)*opt3d(1,1,2)
    opt2d(1,2) = d2mx(3)*opt3d(1,2,1) + dxm1(3)*opt3d(1,2,2)
    opt2d(2,1) = d2mx(3)*opt3d(2,1,1) + dxm1(3)*opt3d(2,1,2)
    opt2d(2,2) = d2mx(3)*opt3d(2,2,1) + dxm1(3)*opt3d(2,2,2)

    ! interpolation in the (fifth, fourth, third and) second dimension
    optout1 = (d2mx(2)*opt2d(1,1) + dxm1(2)*opt2d(1,2))*(invd(5)*invd(4)*invd(3)*invd(2))
    optout2 = (d2mx(2)*opt2d(2,1) + dxm1(2)*opt2d(2,2))*(invd(5)*invd(4)*invd(3)*invd(2))

  end subroutine lininterpol5dim

end module oslo_aero_linear_interp
