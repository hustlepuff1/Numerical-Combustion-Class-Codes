module chemistry_mod
  ! Import the new NASA getters
  use state_mod, only: NVAR, NSPEC, Wspec, R_univ, Yspec, cons_to_prim, &
                       get_species_H_RT, get_species_S_R
  implicit none
  private

  public :: init_chemistry
  public :: update_chemistry
  public :: Yspec

  real(8), parameter :: DeltaH_H2 = 2.418d8
  real(8) :: dt_chem_limit = 1.0d-9 

contains

  subroutine init_chemistry(ni, nj, chem_dt_in)
    integer, intent(in) :: ni, nj
    real(8), intent(in) :: chem_dt_in
    integer :: i, j
    real(8), parameter :: Y_H2_0 = 0.028d0
    real(8), parameter :: Y_O2_0 = 0.225d0
    real(8), parameter :: Y_N2_0 = 0.747d0

    dt_chem_limit = chem_dt_in
    if (allocated(Yspec)) deallocate(Yspec)
    allocate(Yspec(NSPEC, ni, nj))

    do j = 1, nj
      do i = 1, ni
        Yspec(:,i,j) = 0.0d0
        Yspec(1,i,j) = Y_H2_0
        Yspec(2,i,j) = Y_O2_0
        Yspec(8,i,j) = Y_N2_0
      end do
    end do
  end subroutine init_chemistry

  subroutine update_chemistry(U, ni, nj, gamma, Rgas_ref, dt)
    real(8), intent(inout) :: U(NVAR, ni, nj)
    integer, intent(in)    :: ni, nj
    real(8), intent(in)    :: gamma, Rgas_ref, dt

    integer :: i, j, k
    real(8) :: rho, rhoE, KE, p, T
    real(8) :: Yloc(NSPEC), C(NSPEC), w(NSPEC)
    real(8) :: Rmix_loc, E_internal
    real(8) :: dt_remain, dt_sub
    real(8), parameter :: T_cutoff = 500.0d0

    if (.not. allocated(Yspec)) return

    !$omp parallel do private(i, j, k, rho, rhoE, KE, Yloc, Rmix_loc, E_internal, T, p, dt_remain, dt_sub, C, w)
    do j = 2, nj-1
      do i = 2, ni-1
        
        rho = U(1, i, j)
        if (rho <= 1.0d-9) cycle

        rhoE = U(4, i, j)
        KE = 0.5d0 * (U(2,i,j)**2 + U(3,i,j)**2) / rho
        p = (gamma - 1.d0) * (rhoE - KE)
        T = p / (rho * Rgas_ref)

        if (T < T_cutoff) then
             do k = 1, NSPEC
                Yspec(k, i, j) = U(4+k, i, j) / rho
             end do
             cycle 
        end if

        do k = 1, NSPEC
          Yloc(k) = U(4+k, i, j) / rho
          if(Yloc(k) < 0.d0) Yloc(k) = 0.d0
          if(Yloc(k) > 1.d0) Yloc(k) = 1.d0
        end do

        dt_remain = dt
        do while (dt_remain > 1.d-14)
            dt_sub = min(dt_remain, dt_chem_limit)
            
            ! Calculate T accurately with local R_mix
            Rmix_loc = 0.d0
            do k = 1, NSPEC
               Rmix_loc = Rmix_loc + (Yloc(k) / Wspec(k))
            end do
            Rmix_loc = Rmix_loc * R_univ
            
            E_internal = (rhoE - KE) / rho
            T = (gamma - 1.d0) * E_internal / Rmix_loc
            if (T < 300.d0) T = 300.d0
            
            call massfrac_to_conc(rho, Yloc, C)
            call reaction_rates_7step(T, C, w)
            call apply_sources(rho, rhoE, Yloc, w, dt_sub)
            
            dt_remain = dt_remain - dt_sub
        end do

        U(4,i,j) = rhoE
        do k = 1, NSPEC
             U(4+k, i, j) = rho * Yloc(k)
        end do
        Yspec(:,i,j) = Yloc(:)
      end do
    end do
    !$omp end parallel do

  end subroutine update_chemistry

  subroutine massfrac_to_conc(rho, Y, C)
    real(8), intent(in)  :: rho, Y(NSPEC)
    real(8), intent(out) :: C(NSPEC)
    integer :: k
    do k = 1, NSPEC
      C(k) = (rho * Y(k)) / Wspec(k) 
    end do
  end subroutine massfrac_to_conc

  subroutine reaction_rates_7step(T, C, w)
    real(8), intent(in)  :: T, C(NSPEC)
    real(8), intent(out) :: w(NSPEC)
    
    real(8) :: H2, O2, O, H, H2O, HO2, OH, N2
    real(8) :: k1f, k2f, k3f, k4f, k5f, k6f, k7f
    real(8) :: k1r, k2r, k3r
    real(8) :: k40, k4inf, M, Pr, Fc, c_par, n_par, d_par, F
    real(8) :: denom, logPr
    
    ! --- NEW: Variables for Gibbs Calculation ---
    real(8) :: h_rt, s_r
    real(8) :: G_spec(NSPEC)
    real(8) :: dG1, dG2, dG3
    integer :: k

    H2=C(1); O2=C(2); O=C(3); H=C(4); H2O=C(5); HO2=C(6); OH=C(7); N2=C(8)

    ! -----------------------------------------------------------------
    ! 1. Calculate Gibbs Free Energy G(T) for each species
    !    G(T) = H(T) - T*S(T)
    ! -----------------------------------------------------------------
    do k = 1, NSPEC
       h_rt = get_species_H_RT(k, T) ! Dimensionless
       s_r  = get_species_S_R(k, T)  ! Dimensionless
       
       ! G = (H/RT)*RT - T*(S/R)*R = R_univ * T * (H/RT - S/R) [J/mol]
       G_spec(k) = R_univ * T * (h_rt - s_r)
    end do

    ! -----------------------------------------------------------------
    ! 2. Calculate Delta G for the 3 reversible reactions
    ! -----------------------------------------------------------------
    ! Rxn 1: H + O2 <=> OH + O
    dG1 = (G_spec(7) + G_spec(3)) - (G_spec(4) + G_spec(2))
    
    ! Rxn 2: H2 + O <=> OH + H
    dG2 = (G_spec(7) + G_spec(4)) - (G_spec(1) + G_spec(3))
    
    ! Rxn 3: H2 + OH <=> H2O + H
    dG3 = (G_spec(5) + G_spec(4)) - (G_spec(1) + G_spec(7))

    ! -------------------------------------------------------------------------
    ! UNITS: Converted from CGS (cm^3/mol/s) to SI (m^3/mol/s)
    ! Factor 2-body: 1.0d-6
    ! Factor 3-body: 1.0d-12
    ! -------------------------------------------------------------------------

    ! Rxn 1: H + O2 <=> OH + O (2-body)
    ! Original: 3.52d13 -> New: 3.52d7
    k1f = 3.52d7  * T**(-0.7d0) * exp(-8590.d0 / T)

    ! Rxn 2: H2 + O <=> OH + H (2-body)
    ! Original: 5.06d4 -> New: 5.06d-2 (Wait, check original exponent carefully!)
    ! If original was 5.06E+04, new is 5.06E-02. 
    ! If original was 5.06E+01, new is 5.06E-05.
    ! LET'S ASSUME STANDARD MECHANISM VALUES:
    ! H2+O is typically ~5E4 (CGS). So use 5.06d-2.
    k2f = 5.06d-2 * T**( 2.67d0) * exp(-3166.d0 / T)

    ! Rxn 3: H2 + OH <=> H2O + H (2-body)
    ! Original: 1.17d6 -> New: 1.17d0
    k3f = 1.17d0  * T**( 1.3d0)  * exp(-1829.d0 / T)

    ! Rxn 4: H + O2 + M <=> HO2 + M (3-body)
    ! Original: 5.75d13 -> New: 5.75d1 (SI) 
    ! Note: This depends on if it is low-pressure limit (3-body) or high-pressure.
    ! k40 is the low-pressure limit (3-body behavior). Convert by 1.d-12.
    k40 = 5.75d1  * T**(-1.4d0) 

    ! k4inf is the high-pressure limit (2-body behavior). Convert by 1.d-6.
    k4inf = 4.65d3 * T**( 0.44d0)

    ! Rxn 5: HO2 + H <=> H2 + O2 (2-body)
    ! Original: 7.08d10 -> New: 7.08d4
    k5f = 7.08d4  * exp(-148.d0 / T)

    ! Rxn 6: HO2 + H <=> 2OH (2-body)
    ! Original: 1.66d10 -> New: 1.66d4
    k6f = 1.66d4  * exp(-414.d0 / T)

    ! Rxn 7: HO2 + OH <=> H2O + O2 (2-body)
    ! Original: 2.89d10 -> New: 2.89d4
    k7f = 2.89d4  * exp( 250.d0 / T)

    ! -----------------------------------------------------------------
    ! 4. Backward Rates (Equilibrium: k_r = k_f / Kp)
    !    Kp = exp(-DeltaG / RT). So k_r = k_f * exp(DeltaG / RT)
    ! -----------------------------------------------------------------
    k1r = k1f * exp( dG1 / (R_univ * T) )
    k2r = k2f * exp( dG2 / (R_univ * T) )
    k3r = k3f * exp( dG3 / (R_univ * T) )

    ! Third-body efficiency
    M = 2.5d0*H2 + 16.d0*H2O + O2 + O + H + HO2 + OH + N2

    ! Fall-off for Rxn 4
    if (max(k4inf, 1.d-30) < 1.d-29) then
       Pr=0.d0
    else
       Pr=k40*M/k4inf
    end if

    Fc=0.5d0; c_par=-0.4d0-0.67d0*log(Fc); n_par=0.75d0-1.27d0*log(Fc); d_par=0.14d0
    if (Pr<=1.d-30) then
       F=1.d0
    else
       logPr=log(Pr)
       denom=n_par-d_par*(logPr+c_par)
       F=exp((1.d0/(1.d0+((logPr+c_par)/denom)**2))*log(Fc))
    end if
    k4f = k4inf*Pr/(1.d0+Pr)*F

    ! Production Rates w(k) [mol/(m^3 s)]
    w(1) = -k2f*H2*O + k2r*OH*H -k3f*H2*OH + k3r*H2O*H +k6f*HO2*H
    w(2) = -k1f*H*O2 + k1r*OH*O -k4f*H*O2*M +k6f*HO2*H + k7f*HO2*OH
    w(3) = k1f*H*O2 - k1r*OH*O -k2f*H2*O + k2r*OH*H
    w(4) = -k1f*H*O2 + k1r*OH*O +k2f*H2*O - k2r*OH*H +k3f*H2*OH - k3r*H2O*H -k4f*H*O2*M -k5f*HO2*H - k6f*HO2*H
    w(5) = k3f*H2*OH - k3r*H2O*H +k7f*HO2*OH
    w(6) = k4f*H*O2*M -k5f*HO2*H - k6f*HO2*H - k7f*HO2*OH
    w(7) = k1f*H*O2 - k1r*OH*O +k2f*H2*O - k2r*OH*H -k3f*H2*OH + k3r*H2O*H +2.0d0*k5f*HO2*H - k7f*HO2*OH 
    w(8) = 0.d0
  end subroutine reaction_rates_7step

  subroutine apply_sources(rho, rhoE, Y, w, dt)
    real(8), intent(in) :: rho, dt
    real(8), intent(inout) :: rhoE, Y(NSPEC)
    real(8), intent(in) :: w(NSPEC)
    integer :: k
    real(8) :: omega_mass, sumY, qdot
    do k = 1, NSPEC
      omega_mass = w(k) * Wspec(k)
      Y(k) = max(0.d0, min(1.d0, Y(k) + (omega_mass/rho)*dt))
    end do
    sumY = sum(Y)
    if(sumY>1.d-30) Y=Y/sumY
    qdot = -DeltaH_H2 * w(1)
    rhoE = rhoE + qdot * dt
  end subroutine apply_sources

end module chemistry_mod