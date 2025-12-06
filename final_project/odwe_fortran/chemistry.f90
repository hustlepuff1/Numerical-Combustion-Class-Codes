module chemistry_mod
  use state_mod, only: NVAR, NSPEC, Wspec, R_univ, Yspec, cons_to_prim
  implicit none
  private

  public :: init_chemistry
  public :: update_chemistry
  public :: Yspec

  ! Universal gas constant in J/(kmol K)
  real(8), parameter :: Ru_kmol = 8.315d3   ! ~ 8.314e3

  ! Heat of reaction for H2 + 0.5 O2 -> H2O [J/kmol H2]
  real(8), parameter :: DeltaH_H2 = 2.4d8

contains

  !======================================================
  ! Initialize Yspec everywhere (simple premixed H2/air)
  !======================================================
  subroutine init_chemistry(ni, nj)
    integer, intent(in) :: ni, nj
    integer :: i, j

    real(8), parameter :: Y_H2_0 = 0.028d0
    real(8), parameter :: Y_O2_0 = 0.225d0
    real(8), parameter :: Y_N2_0 = 0.747d0

    if (allocated(Yspec)) deallocate(Yspec)
    allocate(Yspec(NSPEC, ni, nj))

    do j = 1, nj
      do i = 1, ni
        Yspec(:,i,j) = 0.0d0
        Yspec(1,i,j) = Y_H2_0   ! H2
        Yspec(2,i,j) = Y_O2_0   ! O2
        Yspec(8,i,j) = Y_N2_0   ! N2

        ! others (O, H, H2O, HO2, OH) start at 0
      end do
    end do
  end subroutine init_chemistry

  !======================================================
  ! Main chemistry update (operator split)
  !  U : [rho, rho*u, rho*v, rho*E]
  !======================================================
  subroutine update_chemistry(U, ni, nj, gamma, Rgas, dt)
    real(8), intent(inout) :: U(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: gamma, Rgas, dt

    integer :: i, j
    real(8) :: rho, uvel, vvel, p, T
    real(8) :: Yloc(NSPEC)
    real(8) :: C(NSPEC)       ! [kmol/m^3]
    real(8) :: w(NSPEC)       ! [kmol/(m^3 s)]
    real(8) :: rhoE
    real(8) :: Rmix

    if (.not. allocated(Yspec)) return

    do j = 2, nj-1
      do i = 2, ni-1

        ! flow state
        call cons_to_prim(U(:,i,j), gamma, rho, uvel, vvel, p)
        if (rho <= 0.d0) cycle

        ! mixture gas constant from composition
        Yloc = Yspec(:,i,j)
        Rmix = Rgas   ! you are using constant Rgas in the rest of the code
        T    = p / (rho * Rmix)
        if (T < 300.d0) T = 300.d0

        call massfrac_to_conc(rho, Yloc, C)
        call reaction_rates_7step(T, C, w)

        rhoE = U(4,i,j)
        call apply_sources(rho, rhoE, Yloc, w, dt)

        U(4,i,j)        = rhoE
        Yspec(:,i,j)    = Yloc

      end do
    end do
  end subroutine update_chemistry

  !------------------------------------------------------
  subroutine massfrac_to_conc(rho, Y, C)
    real(8), intent(in)  :: rho
    real(8), intent(in)  :: Y(NSPEC)
    real(8), intent(out) :: C(NSPEC)
    integer :: k

    do k = 1, NSPEC
      C(k) = (rho * max(Y(k), 0.d0)) / (Wspec(k) * 1.d3) ! [kmol/m^3]
      ! Wspec is kg/mol; multiply by 1e3 to get kg/kmol
    end do
  end subroutine massfrac_to_conc

  !------------------------------------------------------
  ! 7-step mechanism – indices:
  ! 1:H2, 2:O2, 3:O, 4:H, 5:H2O, 6:HO2, 7:OH, 8:N2
  !------------------------------------------------------
  subroutine reaction_rates_7step(T, C, w)
    real(8), intent(in)  :: T, C(NSPEC)
    real(8), intent(out) :: w(NSPEC)

    real(8) :: H2, O2, O, H, H2O, HO2, OH, N2
    real(8) :: k1f, k2f, k3f, k4f, k5f, k6f, k7f
    real(8) :: k1r, k2r, k3r
    real(8) :: k40, k4inf, M, Pr, Fc, c_par, n_par, d_par, F
    real(8) :: gh2, go2, go, gh, goh, gh2o
    real(8) :: g1, g2, g3
    real(8) :: logPr, denom

    H2  = C(1); O2  = C(2); O   = C(3); H   = C(4)
    H2O = C(5); HO2 = C(6); OH  = C(7); N2  = C(8)

    if (T <= 0.d0) then
      w = 0.d0
      return
    end if

    ! Gibbs (J/kmol), rough values from your senior's code
    gh2  = 0.0d0
    go2  = 0.0d0
    go   = 181263.d0 * 1.d3
    gh   = 159830.d0 * 1.d3
    goh  = 21918.d0  * 1.d3
    gh2o = -187100.d0* 1.d3

    g1 = goh + go   - gh  - go2
    g2 = goh + gh   - gh2 - go
    g3 = gh2o + gh  - gh2 - goh

    ! Arrhenius A, n, Ta from the table (A in cm^3, mol, s, K)
    k1f = 3.52d0 * 10.d0**(16.d0 - 3.d0) * T**(-0.7d0) * exp(-8590.d0 / T)
    k2f = 5.06d0 * 10.d0**( 4.d0 - 3.d0) * T**( 2.67d0) * exp(-3166.d0 / T)
    k3f = 1.17d0 * 10.d0**( 9.d0 - 3.d0) * T**( 1.3d0)  * exp(-1829.d0 / T)

    k40   = 5.75d0 * 10.d0**(19.d0 - 3.d0) * T**(-1.4d0)
    k4inf = 4.65d0 * 10.d0**(12.d0 - 3.d0) * T**( 0.44d0)

    k5f = 7.08d0 * 10.d0**(13.d0 - 3.d0) * exp(-148.d0 / T)
    k6f = 1.66d0 * 10.d0**(13.d0 - 3.d0) * exp(-414.d0 / T)
    k7f = 2.89d0 * 10.d0**(13.d0 - 3.d0) * exp( 250.d0 / T)

    k1r = k1f / exp(-g1 / (Ru_kmol * T))
    k2r = k2f / exp(-g2 / (Ru_kmol * T))
    k3r = k3f / exp(-g3 / (Ru_kmol * T))

    M = 2.5d0 * H2 + 16.d0 * H2O + O2 + O + H + HO2 + OH + N2

    if (M <= 0.d0) then
      Pr = 0.d0
    else
      Pr = k40 * M / max(k4inf, 1.d-300)
    end if

    Fc     = 0.5d0
    c_par  = -0.4d0 - 0.67d0 * log(Fc)
    n_par  =  0.75d0 - 1.27d0 * log(Fc)
    d_par  =  0.14d0

    if (Pr <= 0.d0) then
      F = 1.d0
    else
      logPr = log(Pr)
      denom = n_par - d_par * (logPr + c_par)
      F = exp( (1.d0 / (1.d0 + ((logPr + c_par)/denom)**2)) * log(Fc) )
    end if

    k4f = k4inf * Pr / (1.d0 + Pr) * F

    ! Production rates w(k) [kmol/(m^3 s)]
    w(1) = -k2f * H2 * O   + k2r * OH * H &
           -k3f * H2 * OH  + k3r * H2O * H &
           +k6f * HO2 * H

    w(2) = -k1f * H  * O2 + k1r * OH * O &
           -k4f * H  * O2 * M &
           +k6f * HO2 * H  + k7f * HO2 * OH

    w(3) =  k1f * H  * O2 - k1r * OH * O &
           -k2f * H2 * O  + k2r * OH * H

    w(4) = -k1f * H  * O2 + k1r * OH * O &
           +k2f * H2 * O  - k2r * OH * H &
           +k3f * H2 * OH - k3r * H2O * H &
           -k4f * H  * O2 * M &
           -k5f * HO2 * H - k6f * HO2 * H

    w(5) =  k3f * H2 * OH - k3r * H2O * H &
           +k7f * HO2 * OH

    w(6) =  k4f * H  * O2 * M &
           -k5f * HO2 * H - k6f * HO2 * H - k7f * HO2 * OH

    w(7) =  k1f * H  * O2 - k1r * OH * O &
           +k2f * H2 * O  - k2r * OH * H &
           -k3f * H2 * OH + k3r * H2O * H &
           +k5f * HO2 * H + k5f * HO2 * H - k7f * HO2 * OH

    w(8) = 0.d0
  end subroutine reaction_rates_7step

  !------------------------------------------------------
  subroutine apply_sources(rho, rhoE, Y, w, dt)
    real(8), intent(in)    :: rho, dt
    real(8), intent(inout) :: rhoE
    real(8), intent(inout) :: Y(NSPEC)
    real(8), intent(in)    :: w(NSPEC)

    integer :: k
    real(8) :: omega_mass(NSPEC), sumY
    real(8) :: qdot

    if (rho <= 0.d0) return

    do k = 1, NSPEC
      ! w[k] [kmol/m^3/s] * MW[k] [kg/kmol] -> [kg/(m^3 s)]
      omega_mass(k) = w(k) * (Wspec(k) * 1.d3)
    end do

    do k = 1, NSPEC
      Y(k) = Y(k) + (omega_mass(k)/rho) * dt
      if (Y(k) < 0.d0) Y(k) = 0.d0
    end do

    sumY = sum(Y)
    if (sumY > 0.d0) Y = Y / sumY

    ! Heat release from H2 consumption (species 1)
    qdot = -DeltaH_H2 * w(1)         ! [J/(m^3 s)]
    rhoE = rhoE + qdot * dt
  end subroutine apply_sources

end module chemistry_mod
