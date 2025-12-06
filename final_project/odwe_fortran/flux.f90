module flux_mod
  use state_mod, only: NVAR, prim_to_cons, cons_to_prim, sound_speed
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
    real(8), intent(in)  :: Ucons(4), gamma
    real(8), intent(out) :: F(4)
    real(8) :: rho, u, v, p, E

    call cons_to_prim(Ucons, gamma, rho, u, v, p)
    E = Ucons(4) / rho

    F(1) = rho * u
    F(2) = rho * u * u + p
    F(3) = rho * u * v
    F(4) = u * (rho * E + p)
  end subroutine euler_flux_x

  !===========================================================
  ! Euler flux in y-direction: G(U)
  !===========================================================
  subroutine euler_flux_y(Ucons, gamma, G)
    real(8), intent(in)  :: Ucons(4), gamma
    real(8), intent(out) :: G(4)
    real(8) :: rho, u, v, p, E

    call cons_to_prim(Ucons, gamma, rho, u, v, p)
    E = Ucons(4) / rho

    G(1) = rho * v
    G(2) = rho * u * v
    G(3) = rho * v * v + p
    G(4) = v * (rho * E + p)
  end subroutine euler_flux_y

  !===========================================================
  ! Max wave speed in 2D (for CFL)
  !===========================================================
  function max_wave_speed_2d(Ufield, ni, nj, gamma) result(smax)
    real(8), intent(in) :: Ufield(4,ni,nj), gamma
    integer, intent(in) :: ni, nj
    real(8) :: smax
    real(8) :: rho, u, v, p, a
    integer :: i, j

    smax = 0.d0
    do j = 1, nj
      do i = 1, ni
        call cons_to_prim(Ufield(:, i, j), gamma, rho, u, v, p)
        a    = sqrt(gamma * p / rho)
        smax = max(smax, sqrt(u*u + v*v) + a)
      end do
    end do
  end function max_wave_speed_2d

  !===========================================================
  ! Rusanov (local Lax-Friedrichs) flux in direction (nx,ny)
  !===========================================================
  subroutine rusanov_flux_normal(UL, UR, nx, ny, gamma, Fn)
    implicit none
    ! dummies
    real(8), intent(in)  :: UL(4), UR(4)
    real(8), intent(in)  :: nx, ny, gamma
    real(8), intent(out) :: Fn(4)

    ! locals (note: names avoid UL/UR conflict; Fortran is case-insensitive)
    real(8) :: ULloc(4), URloc(4)
    real(8) :: rhoL, uL_s, vL_s, pL, EL
    real(8) :: rhoR, uR_s, vR_s, pR, ER
    real(8) :: FL(4), GL(4), FR(4), GR(4)
    real(8) :: unL, unR, aL, aR, smax
    integer :: k

    ULloc = UL
    URloc = UR

    ! ----- left state -----
    call cons_to_prim(ULloc, gamma, rhoL, uL_s, vL_s, pL)
    EL = ULloc(4) / rhoL

    FL(1) = rhoL * uL_s
    FL(2) = rhoL * uL_s * uL_s + pL
    FL(3) = rhoL * uL_s * vL_s
    FL(4) = uL_s * (rhoL * EL + pL)

    GL(1) = rhoL * vL_s
    GL(2) = rhoL * uL_s * vL_s
    GL(3) = rhoL * vL_s * vL_s + pL
    GL(4) = vL_s * (rhoL * EL + pL)

    ! ----- right state -----
    call cons_to_prim(URloc, gamma, rhoR, uR_s, vR_s, pR)
    ER = URloc(4) / rhoR

    FR(1) = rhoR * uR_s
    FR(2) = rhoR * uR_s * uR_s + pR
    FR(3) = rhoR * uR_s * vR_s
    FR(4) = uR_s * (rhoR * ER + pR)

    GR(1) = rhoR * vR_s
    GR(2) = rhoR * uR_s * vR_s
    GR(3) = rhoR * vR_s * vR_s + pR
    GR(4) = vR_s * (rhoR * ER + pR)

    ! ----- central part & Rusanov dissipation -----
    unL  = uL_s*nx + vL_s*ny
    unR  = uR_s*nx + vR_s*ny
    aL   = sqrt(gamma * pL / rhoL)
    aR   = sqrt(gamma * pR / rhoR)
    smax = max( abs(unL) + aL, abs(unR) + aR )

    do k = 1, 4
      Fn(k) = 0.5d0 * ( nx*FL(k) + ny*GL(k) + nx*FR(k) + ny*GR(k) ) &
              - 0.5d0 * smax * (URloc(k) - ULloc(k))
    end do

  end subroutine rusanov_flux_normal


  !===========================================================
  ! Residual for ODWE (cold Euler): finite-volume Rusanov
  !===========================================================
  subroutine residual_2d_odwe(Ufield, ni, nj, dx, dy, gamma, Rgas, inflow_state, &
                              wall_theta, R)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: Ufield(4,ni,nj)
    real(8), intent(in)    :: dx, dy, gamma, Rgas
    real(8), intent(in)    :: inflow_state(4)
    real(8), intent(in)    :: wall_theta(ni)
    real(8), intent(out)   :: R(4,ni,nj)

    real(8) :: Fx(4,ni+1,nj)     ! fluxes at i±1/2 interfaces
    real(8) :: Gy(4,ni,nj+1)     ! fluxes at j±1/2 interfaces
    real(8) :: UL(4), UR(4), Ughost(4)
    real(8) :: rho, u, v, p, theta, ut, un
    integer :: i, j

      !----------------------------------------
      ! Prevent "unused variable" remark (#7712)
      !----------------------------------------
      if (.false.) then
        print *, "Rgas = ", Rgas
      end if

    !-------------------------
    ! x-direction interfaces
    !-------------------------
    do j = 1, nj
       ! left boundary (inlet): ghost = inflow_state
       UL = inflow_state
       UR = Ufield(:,1,j)
       call rusanov_flux_normal(UL, UR, 1.d0, 0.d0, gamma, Fx(:,1,j))

       ! interior interfaces
       do i = 2, ni
          UL = Ufield(:,i-1,j)
          UR = Ufield(:,i  ,j)
          call rusanov_flux_normal(UL, UR, 1.d0, 0.d0, gamma, Fx(:,i,j))
       end do

       ! right boundary (outlet): zero-gradient (copy last cell)
       UL = Ufield(:,ni,j)
       UR = Ufield(:,ni,j)
       call rusanov_flux_normal(UL, UR, 1.d0, 0.d0, gamma, Fx(:,ni+1,j))
    end do

    !-------------------------
    ! y-direction interfaces
    !-------------------------
    do i = 1, ni
       ! bottom boundary (wedge wall): slip condition
       theta = wall_theta(i)    ! wedge angle at this i (radians)

       UL = Ufield(:,i,1)       ! interior cell
       call cons_to_prim(UL, gamma, rho, u, v, p)

       ! velocity in (t,n) basis
       ut =  u*cos(theta) + v*sin(theta)
       un = -u*sin(theta) + v*cos(theta)

       ! reflect normal component (slip wall)
       un = -un

       ! back to (u,v)
       u =  ut*cos(theta) - un*sin(theta)
       v =  ut*sin(theta) + un*cos(theta)

       call prim_to_cons(rho, u, v, p, gamma, Ughost)
       call rusanov_flux_normal(Ughost, Ufield(:,i,1), 0.d0, 1.d0, gamma, Gy(:,i,1))

       ! interior y-interfaces
       do j = 2, nj
          UL = Ufield(:,i,j-1)
          UR = Ufield(:,i,j  )
          call rusanov_flux_normal(UL, UR, 0.d0, 1.d0, gamma, Gy(:,i,j))
       end do

       ! top boundary: zero-gradient (copy last row)
       UL = Ufield(:,i,nj)
       UR = Ufield(:,i,nj)
       call rusanov_flux_normal(UL, UR, 0.d0, 1.d0, gamma, Gy(:,i,nj+1))
    end do

    !-------------------------
    ! assemble residual
    !-------------------------
    do j = 1, nj
       do i = 1, ni
          R(:,i,j) = -( Fx(:,i+1,j) - Fx(:,i,j) ) / dx  &
                     -( Gy(:,i,j+1) - Gy(:,i,j) ) / dy
       end do
    end do

  end subroutine residual_2d_odwe

end module flux_mod
