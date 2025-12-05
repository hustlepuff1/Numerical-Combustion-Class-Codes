module solver_mod
  use state_mod,    only: NVAR, sound_speed
  use flux_mod,     only: residual_2d_odwe
  use chemistry_mod,only: update_chemistry
  implicit none
  private
  public :: compute_dt_cfl_2d, tvd_rk3_step_2d

contains

  !======================================================
  ! Compute dt from CFL condition
  !======================================================
  function compute_dt_cfl_2d(U, ni, nj, dx, dy, gamma, CFL) result(dt)
    real(8), intent(in) :: U(NVAR,ni,nj)
    integer, intent(in) :: ni, nj
    real(8), intent(in) :: dx, dy, gamma, CFL
    real(8) :: dt
    integer :: i, j
    real(8) :: a, rho, u, v, p
    real(8) :: umax, vmax, amax, smax, local_dt

    umax = 0.d0
    vmax = 0.d0
    amax = 0.d0

    do j = 2, nj-1
      do i = 2, ni-1
        call cons_to_prim(U(:,i,j), gamma, rho, u, v, p)
        a = sqrt(gamma * p / rho)
        umax = max(umax, abs(u))
        vmax = max(vmax, abs(v))
        amax = max(amax, a)
      end do
    end do

    smax = max(umax + amax, vmax + amax)
    if (smax > 0.d0) then
      local_dt = CFL * min(dx, dy) / smax
    else
      local_dt = 1.d10
    end if

    dt = local_dt
  end function compute_dt_cfl_2d

  !======================================================
  ! 3rd-order TVD Runge-Kutta step
  !  - Flow part uses residual_2d_odwe (Euler)
  !  - After that, operator-split chemistry update
  !======================================================
  subroutine tvd_rk3_step_2d(Ufield, ni, nj, dt, dx, dy, gamma, Rgas, &
                             inflow_state, wall_theta)

    use state_mod, only: cons_to_prim   ! for debugging

    real(8), intent(inout) :: Ufield(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: dt, dx, dy, gamma, Rgas
    real(8), intent(in)    :: inflow_state(NVAR)
    real(8), intent(in)    :: wall_theta(ni)

    real(8) :: R1(NVAR, ni, nj)
    real(8) :: R2(NVAR, ni, nj)
    real(8) :: R3(NVAR, ni, nj)
    real(8) :: U1(NVAR, ni, nj)
    real(8) :: U2(NVAR, ni, nj)
    integer :: i, j

    ! ---------- Stage 1 ----------
    call residual_2d_odwe(Ufield, ni, nj, dx, dy, gamma, Rgas, inflow_state, &
                          wall_theta, R1)
    U1 = Ufield + dt * R1

    ! ---------- Stage 2 ----------
    call residual_2d_odwe(U1, ni, nj, dx, dy, gamma, Rgas, inflow_state, &
                          wall_theta, R2)
    U2 = 0.75d0 * Ufield + 0.25d0 * (U1 + dt * R2)

    ! ---------- Stage 3 ----------
    call residual_2d_odwe(U2, ni, nj, dx, dy, gamma, Rgas, inflow_state, &
                          wall_theta, R3)
    Ufield = (1.d0/3.d0) * Ufield + (2.d0/3.d0) * (U2 + dt * R3)

    ! ---------- Chemistry operator-split ----------
    call update_chemistry(Ufield, ni, nj, gamma, Rgas, dt)

  end subroutine tvd_rk3_step_2d

end module solver_mod
