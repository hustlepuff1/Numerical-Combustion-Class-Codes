program odwe_main
  use grid_mod
  use state_mod          ! now also for allocate_state / deallocate_state
  use solver_mod
  use post_mod
  use chemistry_mod
  implicit none

  integer :: ni, nj
  real(8) :: gamma, Rgas
  real(8), allocatable :: x(:), y(:), wall_theta(:)

  real(8), allocatable :: U(:,:,:)   ! NVAR x ni x nj

  real(8) :: x_left, x_right, y_bottom, y_top
  real(8) :: x_ramp, theta_deg
  real(8) :: rho_inf, u_inf, v_inf, p_inf
  real(8) :: inflow_state(NVAR)
  real(8) :: t, t_end, dt, CFL
  real(8) :: dx, dy
  integer :: it, isnap
  integer :: chem_flag
  logical :: use_chemistry

  integer :: ios
  character(len=*), parameter :: input_name = "odwe_input.dat"

  ! ---- read input file ----
  ios = 0
  open(10, file=input_name, status="old", action="read", iostat=ios)
  if (ios /= 0) then
     ! Try parent directory (for when we run from build/)
     close(10)
     open(10, file="../"//input_name, status="old", action="read", iostat=ios)
     if (ios /= 0) then
        write(*,*) "ERROR: could not open ", trim(input_name), &
                   " in '.' or '../'."
        stop
     end if
  end if

  read(10,*)         ! skip comment
  read(10,*) ni, nj
  read(10,*)         ! skip comment
  read(10,*) x_left, x_right, y_bottom, y_top
  read(10,*)         ! skip comment
  read(10,*) x_ramp, theta_deg
  read(10,*)         ! skip comment
  read(10,*) gamma, Rgas
  read(10,*)         ! skip comment
  read(10,*) rho_inf, u_inf, v_inf
  read(10,*)         ! skip comment
  read(10,*) p_inf
  read(10,*)         ! skip comment
  read(10,*) t_end, CFL
  read(10,*)         ! skip comment
  read(10,*) chem_flag

  close(10)

  use_chemistry = (chem_flag /= 0)


  ! === allocate flow + species fields ===
  allocate(x(ni), y(nj), wall_theta(ni))
  call allocate_state(U, ni, nj)

  ! 1) grid & wedge
  call make_wedge_grid(ni, nj, x_left, x_right, y_bottom, y_top, x, y, dx, dy, &
                       wall_theta, x_ramp, theta_deg)

  ! 2) inflow conditions
  rho_inf = 1.d5 / (Rgas * 800.d0)
  u_inf   = 5.d0 * sqrt(gamma * Rgas * 800.d0)
  v_inf   = 0.d0
  p_inf   = 1.d5
  call prim_to_cons(rho_inf, u_inf, v_inf, p_inf, gamma, inflow_state)

  ! 3) initial flow field = uniform inflow
  call init_uniform(U, ni, nj, inflow_state)

  ! 3b) initialize chemistry field (H2/air mixture)
  if (use_chemistry) then
    call init_chemistry(ni, nj)
  end if

  t      = 0.d0
  it     = 0
  isnap  = 0

  do while (t < t_end)

     dt = compute_dt_cfl_2d(U, ni, nj, dx, dy, gamma, CFL)
     if (t + dt > t_end) dt = t_end - t

     call tvd_rk3_step_2d(U, ni, nj, dt, dx, dy, gamma, Rgas, inflow_state, &
                          wall_theta, use_chemistry)

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

  ! clean up
  call deallocate_state(U)

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
    real(8) :: rho, uvel, vvel, p, T

    ! pick a sample cell near the wall
    i_sample = ni/2
    j_sample = 2

    call cons_to_prim(U(:,i_sample,j_sample), gamma, rho, uvel, vvel, p)
    T = temperature_from_primitive(rho, p, Rgas)

    write(*,'("sample cell: rho=",ES12.4," u=",ES12.4," p=",ES12.4," T=",ES12.4)') &
         rho, uvel, p, T
  end subroutine print_sample_cell

end program odwe_main
