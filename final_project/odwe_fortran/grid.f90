module grid_mod
  implicit none
  private
  public :: make_rect_grid, make_wedge_grid

contains

  subroutine make_rect_grid(ni, nj, x_left, x_right, y_bottom, y_top, x, y, dx, dy)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: x_left, x_right, y_bottom, y_top
    real(8), intent(out)   :: x(ni), y(nj), dx, dy
    integer :: i

    dx = (x_right - x_left) / real(ni-1,8)
    dy = (y_top    - y_bottom) / real(nj-1,8)

    do i = 1, ni
       x(i) = x_left + dx * real(i-1,8)
    end do
    do i = 1, nj
       y(i) = y_bottom + dy * real(i-1,8)
    end do
  end subroutine make_rect_grid

  subroutine make_wedge_grid(ni, nj, x_left, x_right, y_bottom, y_top, &
                             x, y, dx, dy, wall_theta, x_ramp, theta_deg)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: x_left, x_right, y_bottom, y_top
    real(8), intent(out)   :: x(ni), y(nj), dx, dy
    real(8), intent(out)   :: wall_theta(ni)
    real(8), intent(in)    :: x_ramp, theta_deg
    integer :: i

    call make_rect_grid(ni, nj, x_left, x_right, y_bottom, y_top, x, y, dx, dy)

    wall_theta(:) = 0.d0
    do i = 1, ni
       if (x(i) >= x_ramp) wall_theta(i) = theta_deg * (acos(-1.d0)/180.d0)
    end do
  end subroutine make_wedge_grid

end module grid_mod
