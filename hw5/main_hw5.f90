!===========================================================
! File: main_hw5.f90
! Main program for Homework #5 model reactive flow
!===========================================================
program main_hw5
  use model_reactive_flow_mod, only: model_reactive_flow, dp
  implicit none

  integer, parameter :: nx      = 201
  real(dp), parameter :: x_min  = 0.0_dp
  real(dp), parameter :: x_max  = 1.0_dp
  real(dp), parameter :: t_final = 0.4_dp
  integer, parameter :: n_mu    = 5

  integer :: mu_list(n_mu)
  type(model_reactive_flow) :: solver
  real(dp), allocatable :: u_imp(:), u_eul(:), u_rk2(:), u_exact(:)
  real(dp) :: mu, dt_conv, dt_reac, dt, cfl
  real(dp) :: xi
  integer :: i, j, mu_int, unit
  character(len=64) :: filename

  ! Set mu values
  mu_list = (/ 1, 5, 50, 500, 1000 /)

  ! Initialize solver grid
  call solver%init(nx, x_min, x_max)

  allocate(u_imp(0:nx-1), u_eul(0:nx-1), u_rk2(0:nx-1), u_exact(0:nx-1))

  cfl = 0.5_dp   ! CFL number for convection

  do j = 1, n_mu
     mu_int = mu_list(j)
     mu     = real(mu_int, dp)

     !------------------------------------------------------
     ! Choose dt: limited by convection and reaction
     ! dt_conv ~ CFL * dx, dt_reac ~ 0.01 / mu (very safe)
     !------------------------------------------------------
     dt_conv = cfl * solver%dx
     dt_reac = 0.01_dp / mu   ! smaller for large mu
     dt      = min(dt_conv, dt_reac)

     write(*,*) 'Solving for mu =', mu, ' with dt =', dt

     ! Implicit
     call solver%solve_implicit(mu, t_final, dt, u_imp)

     ! Explicit Euler (1st order)
     call solver%solve_explicit_euler(mu, t_final, dt, u_eul)

     ! Explicit RK2 (2nd order)
     call solver%solve_explicit_rk2(mu, t_final, dt, u_rk2)

     ! Exact solution at t_final: u(x,t) = u0(x - t)
     do i = 0, nx-1
        xi = solver%x(i) - t_final
        if (xi <= 0.3_dp) then
           u_exact(i) = 1.0_dp
        else
           u_exact(i) = 0.0_dp
        end if
     end do

     !---------------------------------------------------
     ! Write results to file: solution_mu_XXXX.dat
     ! Columns: x, u_exact, u_implicit, u_explicit1, u_explicit2
     !---------------------------------------------------
     write(filename, '("solution_mu_", I4.4, ".dat")') mu_int
     unit = 10 + j

     open(unit=unit, file=trim(filename), status='replace', action='write')
     write(unit,'(A)') '# x  u_exact  u_implicit  u_explicit1  u_explicit2'

     do i = 0, nx-1
        write(unit,'(5(1X,ES16.8))') solver%x(i), u_exact(i), &
                                     u_imp(i), u_eul(i), u_rk2(i)
     end do

     close(unit)
  end do

  deallocate(u_imp, u_eul, u_rk2, u_exact)

end program main_hw5
