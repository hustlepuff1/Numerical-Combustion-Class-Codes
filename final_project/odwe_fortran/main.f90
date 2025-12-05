program odwe_main
  use grid_mod
  use state_mod
  use solver_mod
  use post_mod
  use chemistry_mod
  implicit none

  integer, parameter :: ni = 241, nj = 81
  real(8), parameter :: gamma = 1.4d0, Rgas = 287.d0
  real(8) :: x(ni), y(nj), dx, dy, wall_theta(ni)
  real(8) :: U(NVAR, ni, nj)
  real(8) :: rho_inf, u_inf, v_inf, p_inf
  real(8) :: inflow_state(NVAR)
  real(8) :: t, t_end, dt, CFL
  integer :: it, isnap

  ! 1) grid & wedge
  call make_wedge_grid(ni, nj, 0.d0, 1.d0, 0.d0, 0.2d0, x, y, dx, dy, &
                       wall_theta, 0.3d0, 15.d0)

  ! 2) inflow conditions
  rho_inf = 1.d5 / (Rgas * 800.d0)
  u_inf   = 5.d0 * sqrt(gamma * Rgas * 800.d0)
  v_inf   = 0.d0
  p_inf   = 1.d5
  call prim_to_cons(rho_inf, u_inf, v_inf, p_inf, gamma, inflow_state)

  ! 3) initial flow field = uniform inflow
  call init_uniform(U, ni, nj, inflow_state)

  ! 3b) initialize chemistry field (H2/air mixture)
  call init_chemistry(ni, nj)

  t      = 0.d0
  t_end  = 1.d-4
  CFL    = 0.4d0
  it     = 0
  isnap  = 0

  do while (t < t_end)

     dt = compute_dt_cfl_2d(U, ni, nj, dx, dy, gamma, CFL)
     if (t + dt > t_end) dt = t_end - t

     call tvd_rk3_step_2d(U, ni, nj, dt, dx, dy, gamma, Rgas, inflow_state, &
                          wall_theta)

     t  = t + dt
     it = it + 1

     if (mod(it, 20) == 0 .or. t >= t_end) then
        isnap = isnap + 1
        call write_snapshot(U, ni, nj, x, y, isnap)
        write(*,'(" Wrote snapshot file:snap_",I5.5)') isnap
        write(*,'("it=",I5," t=",ES12.4," dt=",ES10.3)') it, t, dt

        ! debug: print a sample cell on the wall
        call print_sample_cell(U, ni, nj, gamma, Rgas)
     end if

  end do

contains

  subroutine init_uniform(U, ni, nj, U_inflow)
    real(8), intent(out) :: U(NVAR,ni,nj)
    real(8), intent(in)  :: U_inflow(NVAR)
    integer, intent(in)  :: ni, nj
    integer :: i, j
    do j = 1, nj
      do i = 1, ni
        U(:,i,j) = U_inflow
      end do
    end do
  end subroutine init_uniform

  subroutine print_sample_cell(U, ni, nj, gamma, Rgas)
    use state_mod, only: cons_to_prim, temperature_from_primitive
    real(8), intent(in) :: U(NVAR,ni,nj)
    integer, intent(in) :: ni, nj
    real(8), intent(in) :: gamma, Rgas
    integer :: i_sample, j_sample
    real(8) :: rho, u, v, p, T

    i_sample = ni/2
    j_sample = 2   ! near the wall

    call cons_to_prim(U(:,i_sample,j_sample), gamma, rho, u, v, p)
    T = temperature_from_primitive(rho, p, Rgas)
    write(*,'("sample cell: rho=",ES12.4," u=",ES12.4," p=",ES12.4," T=",ES12.4)') &
          rho, u, p, T
  end subroutine print_sample_cell

end program odwe_main
