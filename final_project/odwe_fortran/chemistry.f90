module chemistry_mod
  use state_mod, only: NSPEC, Wspec, R_univ, species_name
  implicit none
  private

  public :: init_chemistry
  public :: update_chemistry

  ! Species mass fractions on the CFD grid:
  !   Yspec(k, i, j)
  real(8), allocatable :: Yspec(:,:,:)

  ! Number of reactions (framework supports 7-step)
  integer, parameter :: NREAC = 7

  ! Arrhenius parameters for each reaction:
  !   k = A * T^n * exp( -Ea / (R_univ*T) )
  !
  real(8), parameter :: A_reac(NREAC) = (/ &
       3.52d16,  & ! 1: H + O2  -> OH + O
       5.06d04,  & ! 2: H2 + O  -> OH + H
       1.17d09,  & ! 3: H2 + OH -> H2O + H
       5.75d19,  & ! 4: H + O2 + M -> HO2 + M  (using k0 only for now)
       7.08d13,  & ! 5: HO2 + H -> OH + OH
       1.66d13,  & ! 6: HO2 + H -> H2 + O2
       2.89d13   & ! 7: HO2 + OH -> H2O + O2
    /)

  real(8), parameter :: n_reac(NREAC) = (/ &
       -0.7d0,  & ! 1
       2.67d0, & ! 2
       1.3d0,  & ! 3
       -1.4d0, & ! 4 (low-pressure limit)
       0.d0,   & ! 5
       0.d0,   & ! 6
       0.d0    & ! 7
    /)

  real(8), parameter :: Ta_reac(NREAC) = (/ &
       8590.d0, & ! 1
       3166.d0, & ! 2
       1829.d0, & ! 3
       0.d0,    & ! 4
       148.d0,  & ! 5
       414.d0,  & ! 6
       -250.d0  & ! 7
    /)

  ! Heat of reaction for each step [J/mol of reaction]
  ! Again: placeholders – replace with proper ΔH_r.
  real(8), parameter :: dH_reac(NREAC) = (/ &
       -2.5d5, -2.5d5, -2.5d5, -2.5d5, -2.5d5, -2.5d5, -2.5d5 /)

contains

  !======================================================
  ! Initialize species field
  !  - currently: uniform H2/air mixture everywhere
  !======================================================
  subroutine init_chemistry(ni, nj)
    integer, intent(in) :: ni, nj
    integer :: i, j

    if (allocated(Yspec)) deallocate(Yspec)
    allocate(Yspec(NSPEC, ni, nj))

    ! Very simple: stoichiometric H2-air, mass fractions
    ! These numbers are approximate; adjust as you like.
    ! Y_H2 ~ 0.028, Y_O2 ~ 0.225, Y_N2 ~ 0.747
    real(8), parameter :: Y_H2_0 = 0.028d0
    real(8), parameter :: Y_O2_0 = 0.225d0
    real(8), parameter :: Y_N2_0 = 0.747d0

    do j = 1, nj
      do i = 1, ni
        Yspec(1,i,j) = Y_H2_0   ! H2
        Yspec(2,i,j) = Y_O2_0   ! O2
        Yspec(3,i,j) = 0.d0     ! H
        Yspec(4,i,j) = 0.d0     ! O
        Yspec(5,i,j) = 0.d0     ! OH
        Yspec(6,i,j) = 0.d0     ! H2O
        Yspec(7,i,j) = Y_N2_0   ! N2
      end do
    end do

  end subroutine init_chemistry

  !======================================================
  ! Update chemistry in every cell (operator split).
  !
  ! For now:
  !  - Uses a *single* global-step reaction H2 + 0.5 O2 -> H2O
  !    with an Arrhenius rate (placeholder).
  !  - Framework is ready to be extended to full 7-step:
  !    define stoichiometry and loop over NREAC.
  !======================================================
  subroutine update_chemistry(U, ni, nj, gamma, Rgas, dt)
    use state_mod, only: NVAR, cons_to_prim, temperature_from_primitive
    real(8), intent(inout) :: U(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: gamma, Rgas, dt

    integer :: i, j
    real(8) :: rho, u, v, p, T
    real(8) :: Y(NSPEC)
    real(8) :: Cmol(NSPEC)   ! molar conc [mol/m^3]
    real(8) :: kf            ! reaction rate coefficient
    real(8) :: wdot          ! overall progress rate [mol/(m^3 s)]
    real(8) :: qdot          ! heat release rate [J/(m^3 s)]
    real(8) :: dY(NSPEC)
    real(8) :: E, q2

    if (.not. allocated(Yspec)) then
       ! Safety: if user forgot to call init_chemistry
       return
    end if

    do j = 2, nj-1
      do i = 2, ni-1

        ! --- Flow state ---
        call cons_to_prim(U(:,i,j), gamma, rho, u, v, p)
        T = temperature_from_primitive(rho, p, Rgas)

        ! --- Local mass fractions ---
        Y(:) = Yspec(:, i, j)

        ! Convert to molar concentrations: C_k [mol/m^3]
        Cmol(:) = (rho * Y(:)) / Wspec(:)

        ! --------------------------------------------
        ! Example: single global reaction
        !   H2 + 0.5 O2 -> H2O
        ! You will later replace this block with the
        ! full 7-step mechanism (loop over NREAC).
        ! --------------------------------------------

        ! Simple Arrhenius for this global step:
        kf = A_reac(1) * T**n_reac(1) * exp( -Ta_reac(1) / T )

        ! reaction rate limited by H2 and O2
        wdot = kf * Cmol(1) * sqrt(max(Cmol(2), 0.d0))   ! [mol/(m^3 s)]

        ! Species source terms in molar form (stoichiometry)
        ! H2: -1, O2: -0.5, H2O: +1
        dY(:) = 0.d0
        if (rho > 0.d0) then
          dY(1) = dY(1) - (1.d0)   * wdot * Wspec(1) / rho   ! H2
          dY(2) = dY(2) - (0.5d0) * wdot * Wspec(2) / rho   ! O2
          dY(6) = dY(6) + (1.d0)   * wdot * Wspec(6) / rho   ! H2O
        end if

        ! Heat release (per m^3): qdot = -ΔH * wdot
        ! dH_reac is [J/mol of reaction]
        qdot = -dH_reac(1) * wdot     ! [J/(m^3 s)]

        ! --------------------------------------------
        ! Explicit Euler update for this small dt
        ! --------------------------------------------
        Y(:) = Y(:) + dt * dY(:)

        ! Clip Y to [0,1] and renormalize
        call enforce_mass_fraction_limits(Y)

        ! Store back
        Yspec(:, i, j) = Y(:)

        ! Update energy equation in U:
        ! d(rho E)/dt = qdot  =>  Δ(rho E) = qdot * dt
        U(4, i, j) = U(4, i, j) + qdot * dt

        ! Optionally you could recompute p,T here if
        ! you wanted gamma or Rgas to depend on Y.

      end do
    end do

  end subroutine update_chemistry

  !======================================================
  ! Utility: enforce 0 <= Y_k <= 1 and sum(Y_k)=1
  !======================================================
  subroutine enforce_mass_fraction_limits(Y)
    real(8), intent(inout) :: Y(NSPEC)
    integer :: k
    real(8) :: sumY

    do k = 1, NSPEC
      if (Y(k) < 0.d0) Y(k) = 0.d0
    end do
    sumY = sum(Y)
    if (sumY > 0.d0) then
      Y(:) = Y(:) / sumY
    end if
  end subroutine enforce_mass_fraction_limits

end module chemistry_mod
