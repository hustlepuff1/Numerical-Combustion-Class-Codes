module bc_mod
  use state_mod, only: NVAR, prim_to_cons, cons_to_prim
  implicit none
  private
  public :: slip_wall_reflect

contains

  subroutine slip_wall_reflect(U_int, theta, gamma, U_ghost)
    real(8), intent(in)  :: U_int(NVAR), theta, gamma
    real(8), intent(out) :: U_ghost(NVAR)
    real(8) :: rho, u, v, p
    real(8) :: t1(2), n1(2), vvec(2), vn(2), vt(2), vref(2)

    call cons_to_prim(U_int, gamma, rho, u, v, p)

    t1(1) = cos(theta);  t1(2) = sin(theta)
    n1(1) = -sin(theta); n1(2) = cos(theta)

    vvec(1) = u; vvec(2) = v
    vn = dot_product(vvec, n1) * n1
    vt = vvec - vn
    vref = vt - vn   ! reflect normal component

    call prim_to_cons(rho, vref(1), vref(2), p, gamma, U_ghost)
  end subroutine slip_wall_reflect

end module bc_mod
