!===========================================================
! File: model_reactive_flow.f90
! Object-like solver for the model reactive flow equation
!===========================================================
module model_reactive_flow_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  !---------------------------------------------------------
  ! "Class": model_reactive_flow
  ! Holds grid and solution arrays + methods
  !---------------------------------------------------------
  type :: model_reactive_flow
     integer :: nx        = 0
     real(dp) :: x_min    = 0.0_dp
     real(dp) :: x_max    = 1.0_dp
     real(dp) :: dx       = 0.0_dp
     real(dp), allocatable :: x(:)      ! grid points
     real(dp), allocatable :: u(:)      ! current solution
     real(dp), allocatable :: rhs(:)    ! RHS for explicit methods
   contains
     procedure :: init                 => mrf_init
     procedure :: set_initial          => mrf_set_initial
     procedure :: compute_rhs          => mrf_compute_rhs
     procedure :: step_euler           => mrf_step_euler
     procedure :: step_rk2             => mrf_step_rk2
     procedure :: step_implicit        => mrf_step_implicit
     procedure :: solve_implicit       => mrf_solve_implicit
     procedure :: solve_explicit_euler => mrf_solve_explicit_euler
     procedure :: solve_explicit_rk2   => mrf_solve_explicit_rk2
  end type model_reactive_flow

contains

  !=========================================================
  ! Reaction term f(u) = -mu * u * (u-1) * (u-0.5)
  !=========================================================
  pure function reaction(u, mu) result(r)
    real(dp), intent(in) :: u, mu
    real(dp) :: r
    r = -mu * u * (u - 1.0_dp) * (u - 0.5_dp)
  end function reaction

  !=========================================================
  ! Initialize grid and allocate arrays
  !=========================================================
  subroutine mrf_init(self, nx, x_min, x_max)
    class(model_reactive_flow), intent(inout) :: self
    integer, intent(in) :: nx
    real(dp), intent(in) :: x_min, x_max
    integer :: i

    self%nx    = nx
    self%x_min = x_min
    self%x_max = x_max

    if (allocated(self%x))   deallocate(self%x)
    if (allocated(self%u))   deallocate(self%u)
    if (allocated(self%rhs)) deallocate(self%rhs)

    allocate(self%x(0:nx-1), self%u(0:nx-1), self%rhs(0:nx-1))

    self%dx = (x_max - x_min) / real(nx-1, dp)

    do i = 0, nx-1
       self%x(i) = x_min + real(i, dp)*self%dx
    end do
  end subroutine mrf_init

  !=========================================================
  ! Set initial condition:
  ! u(x,0) = 1 for x <= 0.3, else 0
  !=========================================================
  subroutine mrf_set_initial(self)
    class(model_reactive_flow), intent(inout) :: self
    integer :: i

    do i = 0, self%nx-1
       if (self%x(i) <= 0.3_dp) then
          self%u(i) = 1.0_dp
       else
          self%u(i) = 0.0_dp
       end if
    end do

    ! Inflow boundary at x=0 kept at 1.0 (consistent with IC extension)
    self%u(0) = 1.0_dp
  end subroutine mrf_set_initial

  !=========================================================
  ! Compute RHS for explicit method-of-lines:
  ! du/dt = -u_x + f(u)
  ! Using upwind for u_x (speed = 1 > 0)
  !=========================================================
  subroutine mrf_compute_rhs(self, mu)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in) :: mu
    integer :: i
    real(dp) :: conv

    ! We do not evolve the boundary point explicitly
    self%rhs(0) = 0.0_dp

    do i = 1, self%nx-1
       conv        = - (self%u(i) - self%u(i-1)) / self%dx
       self%rhs(i) = conv + reaction(self%u(i), mu)
    end do
  end subroutine mrf_compute_rhs

  !=========================================================
  ! One explicit Euler step
  !=========================================================
  subroutine mrf_step_euler(self, mu, dt)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in) :: mu, dt

    call self%compute_rhs(mu)

    ! Update interior points
    self%u(1:self%nx-1) = self%u(1:self%nx-1) + dt * self%rhs(1:self%nx-1)

    ! Inflow boundary at x=0 (Dirichlet)
    self%u(0) = 1.0_dp
  end subroutine mrf_step_euler

  !=========================================================
  ! One explicit RK2 (Heun) step
  !=========================================================
  subroutine mrf_step_rk2(self, mu, dt)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in) :: mu, dt
    real(dp), allocatable :: u_old(:), rhs1(:), rhs2(:)
    integer :: nx_loc

    nx_loc = self%nx

    allocate(u_old(0:nx_loc-1), rhs1(0:nx_loc-1), rhs2(0:nx_loc-1))

    u_old = self%u

    ! k1
    call self%compute_rhs(mu)
    rhs1 = self%rhs

    ! Provisional value
    self%u(1:nx_loc-1) = u_old(1:nx_loc-1) + dt * rhs1(1:nx_loc-1)
    self%u(0)          = 1.0_dp   ! enforce BC

    ! k2
    call self%compute_rhs(mu)
    rhs2 = self%rhs

    ! Final RK2 update
    self%u(1:nx_loc-1) = u_old(1:nx_loc-1) + 0.5_dp * dt * &
                         (rhs1(1:nx_loc-1) + rhs2(1:nx_loc-1))
    self%u(0)          = 1.0_dp   ! BC again

    deallocate(u_old, rhs1, rhs2)
  end subroutine mrf_step_rk2

  !=========================================================
  ! One linearized implicit step:
  !
  ! (u_i^{n+1} - u_i^n)/dt + (u_i^{n+1} - u_{i-1}^{n+1})/dx = f(u_i^n)
  !
  ! Solve for u^{n+1} using forward substitution
  !=========================================================
  subroutine mrf_step_implicit(self, mu, dt)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in) :: mu, dt
    real(dp), allocatable :: u_old(:), u_new(:)
    integer :: i, nx_loc
    real(dp) :: rhs_i, alpha

    nx_loc = self%nx

    allocate(u_old(0:nx_loc-1), u_new(0:nx_loc-1))

    u_old = self%u

    alpha = dt / self%dx

    ! Inflow BC at i=0
    u_new(0) = 1.0_dp

    ! Interior points i = 1..nx-1
    do i = 1, nx_loc-1
       rhs_i   = u_old(i) + dt * reaction(u_old(i), mu)
       u_new(i) = (rhs_i + alpha * u_new(i-1)) / (1.0_dp + alpha)
    end do

    self%u = u_new

    deallocate(u_old, u_new)
  end subroutine mrf_step_implicit

  !=========================================================
  ! Full solve: implicit method to t_final
  !=========================================================
  subroutine mrf_solve_implicit(self, mu, t_final, dt, u_out)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in)  :: mu, t_final, dt
    real(dp), intent(out) :: u_out(:)
    integer :: nsteps, n
    real(dp) :: dt_eff, t

    ! Choose integer number of steps and adjust dt slightly
    nsteps = ceiling(t_final / dt)
    dt_eff = t_final / real(nsteps, dp)

    call self%set_initial()

    t = 0.0_dp
    do n = 1, nsteps
       call self%step_implicit(mu, dt_eff)
       t = t + dt_eff
    end do

    u_out = self%u
  end subroutine mrf_solve_implicit

  !=========================================================
  ! Full solve: explicit Euler to t_final
  !=========================================================
  subroutine mrf_solve_explicit_euler(self, mu, t_final, dt, u_out)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in)  :: mu, t_final, dt
    real(dp), intent(out) :: u_out(:)
    integer :: nsteps, n
    real(dp) :: dt_eff, t

    nsteps = ceiling(t_final / dt)
    dt_eff = t_final / real(nsteps, dp)

    call self%set_initial()

    t = 0.0_dp
    do n = 1, nsteps
       call self%step_euler(mu, dt_eff)
       t = t + dt_eff
    end do

    u_out = self%u
  end subroutine mrf_solve_explicit_euler

  !=========================================================
  ! Full solve: explicit RK2 to t_final
  !=========================================================
  subroutine mrf_solve_explicit_rk2(self, mu, t_final, dt, u_out)
    class(model_reactive_flow), intent(inout) :: self
    real(dp), intent(in)  :: mu, t_final, dt
    real(dp), intent(out) :: u_out(:)
    integer :: nsteps, n
    real(dp) :: dt_eff, t

    nsteps = ceiling(t_final / dt)
    dt_eff = t_final / real(nsteps, dp)

    call self%set_initial()

    t = 0.0_dp
    do n = 1, nsteps
       call self%step_rk2(mu, dt_eff)
       t = t + dt_eff
    end do

    u_out = self%u
  end subroutine mrf_solve_explicit_rk2

end module model_reactive_flow_mod
