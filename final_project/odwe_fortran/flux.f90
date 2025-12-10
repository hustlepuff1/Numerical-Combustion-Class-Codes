module flux_mod
  use state_mod, only: NVAR, NSPEC, prim_to_cons, cons_to_prim, sound_speed
  implicit none
  private
  public :: euler_flux_x, euler_flux_y
  public :: max_wave_speed_2d
  public :: rusanov_flux_normal
  public :: residual_2d_odwe

contains

  subroutine euler_flux_x(Ucons, gamma, F)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: F(NVAR)
    real(8) :: rho, u, v, p, E
    integer :: k
    rho = Ucons(1)
    if (rho <= 1.d-30) then
       u = 0.d0; v = 0.d0; p = 0.d0
    else
       u = Ucons(2) / rho; v = Ucons(3) / rho
       E = Ucons(4) / rho; p = (gamma - 1.d0) * rho * (E - 0.5d0*(u*u + v*v))
    end if
    F(1) = rho * u
    F(2) = rho * u * u + p
    F(3) = rho * u * v
    F(4) = (Ucons(4) + p) * u
    do k = 5, NVAR
       F(k) = Ucons(k) * u
    end do
  end subroutine euler_flux_x

  subroutine euler_flux_y(Ucons, gamma, G)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: G(NVAR)
    real(8) :: rho, u, v, p, E
    integer :: k
    rho = Ucons(1)
    if (rho <= 1.d-30) then
       u = 0.d0; v = 0.d0; p = 0.d0
    else
       u = Ucons(2) / rho; v = Ucons(3) / rho
       E = Ucons(4) / rho; p = (gamma - 1.d0) * rho * (E - 0.5d0*(u*u + v*v))
    end if
    G(1) = rho * v
    G(2) = rho * u * v
    G(3) = rho * v * v + p
    G(4) = (Ucons(4) + p) * v
    do k = 5, NVAR
       G(k) = Ucons(k) * v
    end do
  end subroutine euler_flux_y

  function max_wave_speed_2d(Ucons, gamma) result(s)
    real(8), intent(in) :: Ucons(NVAR), gamma
    real(8) :: s, rho, u, v, p, a
    call cons_to_prim(Ucons, gamma, rho, u, v, p)
    if (p < 1.d-9) p = 1.d-9
    a = sqrt(gamma*p/rho)
    s = sqrt(u*u + v*v) + a
  end function max_wave_speed_2d

  subroutine rusanov_flux_normal(UL, UR, nx, ny, gamma, F_num)
    real(8), intent(in)  :: UL(NVAR), UR(NVAR)
    real(8), intent(in)  :: nx, ny, gamma
    real(8), intent(out) :: F_num(NVAR)
    real(8) :: FL(NVAR), FR(NVAR), GL(NVAR), GR(NVAR)
    real(8) :: FluxL(NVAR), FluxR(NVAR), sL, sR, s_max
    integer :: k

    call euler_flux_x(UL, gamma, FL); call euler_flux_y(UL, gamma, GL)
    call euler_flux_x(UR, gamma, FR); call euler_flux_y(UR, gamma, GR)

    do k = 1, NVAR
       FluxL(k) = FL(k)*nx + GL(k)*ny
       FluxR(k) = FR(k)*nx + GR(k)*ny
    end do
    sL = max_wave_speed_2d(UL, gamma); sR = max_wave_speed_2d(UR, gamma)
    s_max = max(sL, sR)
    do k = 1, NVAR
       F_num(k) = 0.5d0 * (FluxL(k) + FluxR(k)) - 0.5d0 * s_max * (UR(k) - UL(k))
    end do
  end subroutine rusanov_flux_normal

  subroutine residual_2d_odwe(Ufield, ni, nj, dx, dy, gamma, Rgas, &
                              inflow_state, wall_theta, R)
    use bc_mod, only: slip_wall_reflect
    real(8), intent(in)    :: Ufield(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: dx, dy, gamma, Rgas
    real(8), intent(in)    :: inflow_state(NVAR)
    real(8), intent(in)    :: wall_theta(ni)
    real(8), intent(out)   :: R(NVAR, ni, nj)

    real(8), allocatable   :: Fx(:,:,:), Gy(:,:,:)
    integer :: i, j, k
    real(8) :: UL(NVAR), UR(NVAR), U_ghost(NVAR)
    real(8) :: area_factor, nx, ny, y_frac
    real(8) :: h_ratio_i, h_ratio_ip1, h_ratio_cell

    allocate(Fx(NVAR, ni+1, nj))
    allocate(Gy(NVAR, ni, nj+1))

    ! X-Flux
    !$omp parallel do default(shared) private(j, i, UL, UR)
    do j = 1, nj
       do i = 1, ni+1
          if (i == 1) then
             UL = inflow_state; UR = Ufield(:, 1, j)
          else if (i == ni+1) then
             UL = Ufield(:, ni, j); UR = Ufield(:, ni, j)
          else
             UL = Ufield(:, i-1, j); UR = Ufield(:, i,   j)
          end if
          call rusanov_flux_normal(UL, UR, 1.d0, 0.d0, gamma, Fx(:,i,j))
       end do
    end do
    !$omp end parallel do

    ! Y-Flux
    !$omp parallel do default(shared) private(i, j, UL, UR, U_ghost, area_factor, nx, ny, y_frac)
    do i = 1, ni
       ! Bottom (Wall) - Uses Ghost Cell
       call slip_wall_reflect(Ufield(:,i,1), wall_theta(i), gamma, U_ghost)
       U_ghost(5:NVAR) = Ufield(5:NVAR, i, 1) ! Explicitly copy species
       
       nx = -sin(wall_theta(i)); ny = cos(wall_theta(i))
       call rusanov_flux_normal(U_ghost, Ufield(:,i,1), nx, ny, gamma, Gy(:,i,1))
       Gy(:,i,1) = Gy(:,i,1) / cos(wall_theta(i)) ! Area correction

       ! Interior
       do j = 2, nj
          y_frac = real(j-1, 8) / real(nj-1, 8)
          nx = -(tan(wall_theta(i)) * (1.d0 - y_frac))
          ny = 1.d0
          area_factor = sqrt(nx*nx + ny*ny)
          call rusanov_flux_normal(Ufield(:,i,j-1), Ufield(:,i,j), nx/area_factor, ny/area_factor, gamma, Gy(:,i,j))
          Gy(:,i,j) = Gy(:,i,j) * area_factor
       end do

       ! Top (Outlet)
       call rusanov_flux_normal(Ufield(:,i,nj), Ufield(:,i,nj), 0.d0, 1.d0, gamma, Gy(:,i,nj+1))
    end do
    !$omp end parallel do

    ! Residuals with Area Correction (Compresses the gas as wedge narrows)
    !$omp parallel do default(shared) private(j, i, k, h_ratio_i, h_ratio_ip1, h_ratio_cell)
    do j = 1, nj
       do i = 1, ni
          ! Channel height ratios
          h_ratio_i   = 1.0d0 - (tan(wall_theta(max(1,i-1))) * real(max(0,i-2),8)*dx) / (real(nj-1,8)*dy)
          h_ratio_ip1 = 1.0d0 - (tan(wall_theta(i))          * real(i-1,8)*dx)       / (real(nj-1,8)*dy)
          if(h_ratio_i < 0.1) h_ratio_i = 0.1
          if(h_ratio_ip1 < 0.1) h_ratio_ip1 = 0.1
          h_ratio_cell = 0.5d0 * (h_ratio_i + h_ratio_ip1)

          do k = 1, NVAR
             ! Flux Balance with Area Terms
             R(k,i,j) = -( Fx(k,i+1,j)*h_ratio_ip1 - Fx(k,i,j)*h_ratio_i ) / (dx * h_ratio_cell) &
                        -( Gy(k,i,j+1) - Gy(k,i,j) ) / (dy * h_ratio_cell)
          end do
       end do
    end do
    !$omp end parallel do

    deallocate(Fx, Gy)
  end subroutine residual_2d_odwe
end module flux_mod