module solver_mod
  use state_mod,  only: NVAR, cons_to_prim, sound_speed
  use flux_mod,   only: residual_2d_odwe
  use chemistry_mod, only: update_chemistry
  implicit none
  private
  public :: compute_dt_cfl_2d, tvd_rk3_step_2d

contains

  !======================================================
  ! Compute time step from CFL condition
  !======================================================
  function compute_dt_cfl_2d(U, ni, nj, dx, dy, gamma, CFL) result(dt)
    real(8), intent(in) :: U(NVAR, ni, nj)
    integer, intent(in) :: ni, nj
    real(8), intent(in) :: dx, dy, gamma, CFL
    real(8) :: dt

    integer :: i, j
    real(8) :: rho, uvel, vvel, p, a
    real(8) :: umax, vmax, wavespeed

    umax = 0.d0
    vmax = 0.d0

    do j = 1, nj
      do i = 1, ni
        call cons_to_prim(U(:,i,j), gamma, rho, uvel, vvel, p)
        a = sound_speed(U(:,i,j), gamma)

        umax = max(umax, abs(uvel) + a)
        vmax = max(vmax, abs(vvel) + a)
      end do
    end do

    wavespeed = max(umax/dx, vmax/dy)

    if (wavespeed > 0.d0) then
      dt = CFL / wavespeed
    else
      dt = 1.d30   ! fallback
    end if
  end function compute_dt_cfl_2d

  !======================================================
  ! 2D TVD RK3 step (Euler + chemistry)
  !======================================================
  subroutine tvd_rk3_step_2d(Ufield, ni, nj, dt, dx, dy, gamma, Rgas, &
                             inflow_state, wall_theta, use_chemistry)
    real(8), intent(inout) :: Ufield(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: dt, dx, dy, gamma, Rgas
    real(8), intent(in)    :: inflow_state(NVAR)
    real(8), intent(in)    :: wall_theta(ni)
    logical, intent(in)    :: use_chemistry

    ! --- 1. Declare as ALLOCATABLE arrays (Dynamic Memory) ---
    ! Note: Do NOT specify size here (e.g. ni, nj). Use (:,:,:)
    real(8), allocatable :: U1(:,:,:)
    real(8), allocatable :: U2(:,:,:)
    real(8), allocatable :: R1(:,:,:)
    real(8), allocatable :: R2(:,:,:)
    real(8), allocatable :: R3(:,:,:)

    ! --- 2. Allocate memory explicitly (Heap Memory) ---
    allocate(U1(NVAR, ni, nj), U2(NVAR, ni, nj))
    allocate(R1(NVAR, ni, nj), R2(NVAR, ni, nj), R3(NVAR, ni, nj))

    ! --- Stage 1 ---
    call residual_2d_odwe(Ufield, ni, nj, dx, dy, gamma, Rgas, &
                          inflow_state, wall_theta, R1)
    U1 = Ufield + dt * R1

    ! --- Stage 2 ---
    call residual_2d_odwe(U1, ni, nj, dx, dy, gamma, Rgas, &
                          inflow_state, wall_theta, R2)
    U2 = 0.75d0 * Ufield + 0.25d0 * (U1 + dt * R2)

    ! --- Stage 3 ---
    call residual_2d_odwe(U2, ni, nj, dx, dy, gamma, Rgas, &
                          inflow_state, wall_theta, R3)
    Ufield = (1.d0/3.d0) * Ufield + (2.d0/3.d0) * (U2 + dt * R3)

    ! --- operator-split chemistry ---
    if (use_chemistry) then
      call update_chemistry(Ufield, ni, nj, gamma, Rgas, dt)
    end if

    ! --- 3. Deallocate memory to prevent leaks ---
    deallocate(U1, U2, R1, R2, R3)

  end subroutine tvd_rk3_step_2d

end module solver_mod