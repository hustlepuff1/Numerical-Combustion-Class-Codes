module flux_mod
  use state_mod, only: NVAR, NSPEC, prim_to_cons, cons_to_prim, sound_speed
  implicit none
  private
  public :: euler_flux_x, euler_flux_y
  public :: max_wave_speed_2d
  public :: rusanov_flux_normal
  public :: residual_2d_odwe

contains

  !===========================================================
  ! Euler flux in x-direction: F(U)
  !===========================================================
  subroutine euler_flux_x(Ucons, gamma, F)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: F(NVAR)
    real(8) :: rho, u, v, p, E
    integer :: k

    rho = Ucons(1)
    if (rho <= 1.d-30) then
       u = 0.d0; v = 0.d0; p = 0.d0
    else
       u = Ucons(2) / rho
       v = Ucons(3) / rho
       E = Ucons(4) / rho
       p = (gamma - 1.d0) * rho * (E - 0.5d0*(u*u + v*v))
    end if

    F(1) = rho * u
    F(2) = rho * u * u + p
    F(3) = rho * u * v
    F(4) = (Ucons(4) + p) * u
    
    do k = 5, NVAR
       F(k) = Ucons(k) * u
    end do
  end subroutine euler_flux_x

  !===========================================================
  ! Euler flux in y-direction: G(U)
  !===========================================================
  subroutine euler_flux_y(Ucons, gamma, G)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: G(NVAR)
    real(8) :: rho, u, v, p, E
    integer :: k

    rho = Ucons(1)
    if (rho <= 1.d-30) then
       u = 0.d0; v = 0.d0; p = 0.d0
    else
       u = Ucons(2) / rho
       v = Ucons(3) / rho
       E = Ucons(4) / rho
       p = (gamma - 1.d0) * rho * (E - 0.5d0*(u*u + v*v))
    end if

    G(1) = rho * v
    G(2) = rho * v * u
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

  !===========================================================
  ! Rusanov Flux (Local Lax-Friedrichs)
  ! F_face = 0.5*(F_L + F_R) - 0.5*lambda_max*(U_R - U_L)
  !===========================================================
  subroutine rusanov_flux_normal(UL, UR, nx, ny, gamma, F_num)
    real(8), intent(in)  :: UL(NVAR), UR(NVAR)
    real(8), intent(in)  :: nx, ny, gamma
    real(8), intent(out) :: F_num(NVAR)
    
    real(8) :: FL(NVAR), FR(NVAR), GL(NVAR), GR(NVAR)
    real(8) :: FluxL(NVAR), FluxR(NVAR)
    real(8) :: sL, sR, s_max
    integer :: k

    call euler_flux_x(UL, gamma, FL)
    call euler_flux_y(UL, gamma, GL)
    call euler_flux_x(UR, gamma, FR)
    call euler_flux_y(UR, gamma, GR)

    do k = 1, NVAR
       FluxL(k) = FL(k)*nx + GL(k)*ny
       FluxR(k) = FR(k)*nx + GR(k)*ny
    end do

    sL = max_wave_speed_2d(UL, gamma)
    sR = max_wave_speed_2d(UR, gamma)
    s_max = max(sL, sR)

    do k = 1, NVAR
       F_num(k) = 0.5d0 * (FluxL(k) + FluxR(k)) - &
                  0.5d0 * s_max * (UR(k) - UL(k))
    end do
  end subroutine rusanov_flux_normal

  !===========================================================
  ! Residual Computation (RHS)
  ! R = - [ (F_{i+1/2} - F_{i-1/2})/dx + (G_{j+1/2} - G_{j-1/2})/dy ]
  ! Includes Parallel Directives (OpenMP)
  !===========================================================
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
    real(8) :: sigma, area_factor

    allocate(Fx(NVAR, ni+1, nj))
    allocate(Gy(NVAR, ni, nj+1))

    ! -----------------------------------------------------
    ! 1. Flux in X-direction (Interfaces i=1..ni+1)
    ! -----------------------------------------------------
    !$omp parallel do default(shared) private(j, i, UL, UR)
    do j = 1, nj
       do i = 1, ni+1
          if (i == 1) then
             UL = inflow_state
             UR = Ufield(:, 1, j)
          else if (i == ni+1) then
             UL = Ufield(:, ni, j)
             UR = Ufield(:, ni, j) ! Zero-order extrapolation
          else
             UL = Ufield(:, i-1, j)
             UR = Ufield(:, i,   j)
          end if
          call rusanov_flux_normal(UL, UR, 1.d0, 0.d0, gamma, Fx(:,i,j))
       end do
    end do
    !$omp end parallel do

    ! -----------------------------------------------------
    ! 2. Flux in Y-direction (Interfaces j=1..nj+1)
    ! -----------------------------------------------------
    !$omp parallel do default(shared) private(i, j, UL, UR, U_ghost, area_factor)
    do i = 1, ni
       do j = 1, nj+1
          area_factor = 1.0d0

          if (j == 1) then
             ! Bottom Wall (Slip / Reflective)
             ! Wall angle changes at i, handled by wall_theta(i)
             call slip_wall_reflect(Ufield(:,i,1), wall_theta(i), gamma, U_ghost)

             U_ghost(5:NVAR) = Ufield(5:NVAR, i, 1)
             UL = U_ghost
             UR = Ufield(:, i, 1)
             
             ! Area correction for inclined face (simple length ratio)
             area_factor = 1.0d0 / cos(wall_theta(i))

          else if (j == nj+1) then
             ! Top Wall (Far field / Outlet)
             UL = Ufield(:, i, nj)
             UR = Ufield(:, i, nj)
          else
             UL = Ufield(:, i, j-1)
             UR = Ufield(:, i, j)
          end if

          call rusanov_flux_normal(UL, UR, 0.d0, 1.d0, gamma, Gy(:,i,j))
          Gy(:,i,j) = Gy(:,i,j) * area_factor
       end do
    end do
    !$omp end parallel do

    ! -----------------------------------------------------
    ! 3. Compute Residuals
    ! -----------------------------------------------------
    !$omp parallel do default(shared) private(j, i, k, sigma)
    do j = 1, nj
       do i = 1, ni
          ! Volume correction (approximate for wedge)
          sigma = 1.0d0 - (tan(wall_theta(i)) * real(i-1,8)*dx) / (dy * real(nj-1,8))
          if (sigma < 0.1d0) sigma = 0.1d0 
          
          do k = 1, NVAR
             R(k,i,j) = -( Fx(k,i+1,j) - Fx(k,i,j) ) / dx &
                        -( Gy(k,i,j+1) - Gy(k,i,j) ) / (dy * sigma)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine residual_2d_odwe

end module flux_mod