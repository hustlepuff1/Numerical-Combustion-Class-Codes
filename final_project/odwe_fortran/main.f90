program odwe_main
  use grid_mod
  use state_mod          
  use solver_mod
  use post_mod
  use chemistry_mod
  implicit none

  integer :: ni, nj
  integer :: i, j, k
  real(8) :: gamma, Rgas
  real(8), allocatable :: x(:), y(:), wall_theta(:)
  real(8), allocatable :: U(:,:,:)   

  real(8) :: x_left, x_right, y_bottom, y_top
  real(8) :: x_ramp, theta_deg
  real(8) :: rho_inf, u_inf, v_inf, p_inf
  real(8) :: inflow_state(NVAR)
  real(8) :: t, t_end, dt, CFL
  real(8) :: dx, dy
  integer :: it, isnap
  integer :: chem_flag
  logical :: use_chemistry
  real(8) :: max_chem_dt_input
  integer :: isnap_freq 

  integer :: ios
  character(len=*), parameter :: input_name = "odwe_input.dat"

  ! ---- read input file ----
  ios = 0
  open(10, file=input_name, status="old", action="read", iostat=ios)
  if (ios /= 0) then
      close(10)
      open(10, file="../"//input_name, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         write(*,*) "ERROR: could not open input file ", input_name
         stop
      end if
  end if

  read(10,*) 
  read(10,*) ni, nj
  read(10,*)
  read(10,*) x_left, x_right, y_bottom, y_top
  read(10,*)
  read(10,*) x_ramp, theta_deg
  read(10,*)
  read(10,*) gamma, Rgas
  read(10,*)
  read(10,*) rho_inf, u_inf, v_inf
  read(10,*)
  read(10,*) p_inf
  read(10,*)
  read(10,*) t_end, CFL
  read(10,*)
  read(10,*) chem_flag
  read(10,*)
  read(10,*) max_chem_dt_input
  read(10,*)
  read(10,*) isnap_freq
  close(10)

  use_chemistry = (chem_flag == 1)

  ! --- FIX 1: Allocate Grid Arrays BEFORE usage ---
  allocate(x(ni))
  allocate(y(nj))
  allocate(wall_theta(ni))

  call make_wedge_grid(ni, nj, x_left, x_right, y_bottom, y_top, &
                       x, y, dx, dy, wall_theta, x_ramp, theta_deg)

  call allocate_state(U, ni, nj)

  ! --- 1. INITIALIZE NASA COEFFICIENTS ---
  call init_nasa_coeffs()

  ! --- 2. Initialize Chemistry (if needed) ---
  if (use_chemistry) then
     call init_chemistry(ni, nj, max_chem_dt_input)
  end if

  ! --- 3. Initialize Flow ---
  inflow_state = 0.d0
  inflow_state(1) = rho_inf
  inflow_state(2) = rho_inf * u_inf
  inflow_state(3) = rho_inf * v_inf
  inflow_state(4) = p_inf / (gamma - 1.d0) + 0.5d0*rho_inf*(u_inf**2 + v_inf**2)
  if (use_chemistry) then
    inflow_state(5) = rho_inf * 0.028d0 ! H2
    inflow_state(6) = rho_inf * 0.225d0 ! O2
    inflow_state(7) = 0.d0              ! O
    inflow_state(8) = 0.d0              ! H
    inflow_state(9) = 0.d0              ! H2O
    inflow_state(10)= 0.d0              ! HO2
    inflow_state(11)= 0.d0              ! OH
    inflow_state(12)= rho_inf * 0.747d0 ! N2
  end if
  
  call init_uniform(U, ni, nj, inflow_state)

  t = 0.d0
  it = 0
  isnap = 0

  ! Write initial state
  call write_snapshot(U, ni, nj, x, y, isnap)

  ! ---- Time Loop ----
  do while (t < t_end)
      
      dt = compute_dt_cfl_2d(U, ni, nj, dx, dy, gamma, CFL)
      if (t + dt > t_end) dt = t_end - t

      call tvd_rk3_step_2d(U, ni, nj, dt, dx, dy, gamma, Rgas, &
                           inflow_state, &
                           wall_theta, use_chemistry)

      !$omp parallel do private(j, i, k)
      do j = 1, nj
         do i = 1, ni
            do k = 5, NVAR
               if (U(k, i, j) < 0.d0) U(k, i, j) = 0.d0
            end do
         end do
      end do
      !$omp end parallel do
      ! --------------------------------------------

      t  = t + dt
      it = it + 1

      if (mod(it, isnap_freq) == 0 .or. t >= t_end) then
        isnap = isnap + 1
        call write_snapshot(U, ni, nj, x, y, isnap)
        write(*,'(" Wrote snapshot file:snap_",I5.5)') isnap
        write(*,'("it=",I5," t=",ES12.4," dt=",ES10.3)') it, t, dt
        call print_sample_cell(U, ni, nj, gamma, Rgas)
      end if

  end do

  ! Cleanup
  call deallocate_state(U)
  if (allocated(x)) deallocate(x)
  if (allocated(y)) deallocate(y)
  if (allocated(wall_theta)) deallocate(wall_theta)

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
    
    real(8) :: rho_loc, u_vel, v_vel, p_loc, T_loc
    integer :: i_probe, j_probe
    
    i_probe = ni/2
    j_probe = nj/2
    
    call cons_to_prim(U(:,i_probe,j_probe), gamma, rho_loc, u_vel, v_vel, p_loc)
    T_loc = temperature_from_primitive(rho_loc, p_loc, Rgas)
    
    write(*,'(" Probe (",I3,",",I3,"): P=",ES10.3," T=",F7.1," Rho=",ES10.3)') &
          i_probe, j_probe, p_loc, T_loc, rho_loc
  end subroutine print_sample_cell

end program odwe_main