module chemistry_mod
  use state_mod, only: NVAR, NSPEC, Wspec, R_univ, Yspec, cons_to_prim, &
                       get_species_H_RT, get_species_S_R
  implicit none
  private

  public :: init_chemistry
  public :: update_chemistry
  public :: Yspec

  ! Lower Heating Value of H2 ~120 MJ/kg
  real(8), parameter :: DeltaH_H2 = 1.1996d8 
  real(8) :: dt_chem_limit = 1.0d-10

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
    real(8), parameter :: T_cutoff = 300.0d0

    if (.not. allocated(Yspec)) return

    !$omp parallel do private(i, j, k, rho, rhoE, KE, Yloc, Rmix_loc, E_internal, T, p, dt_remain, dt_sub, C, w)
    do j = 1, nj
      do i = 1, ni
        
        rho = U(1, i, j)
        if (rho <= 1.0d-9) cycle

        rhoE = U(4, i, j)
        KE = 0.5d0 * (U(2,i,j)**2 + U(3,i,j)**2) / rho
        p = (gamma - 1.d0) * (rhoE - KE)
        T = p / (rho * Rgas_ref)

        ! Optimization: Skip chemistry if gas is cold
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
            
            ! Approximate T for rate calculation
            E_internal = (rhoE - KE) / rho
            T = (gamma - 1.d0) * E_internal / Rgas_ref
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
    
    real(8) :: h_rt, s_r
    real(8) :: G_spec(NSPEC)
    real(8) :: dG1, dG2, dG3
    integer :: k

    H2=C(1); O2=C(2); O=C(3); H=C(4); H2O=C(5); HO2=C(6); OH=C(7); N2=C(8)

    ! Gibbs Calculation
    do k = 1, NSPEC
       h_rt = get_species_H_RT(k, T)
       s_r  = get_species_S_R(k, T)
       G_spec(k) = R_univ * T * (h_rt - s_r)
    end do

    dG1 = (G_spec(7) + G_spec(3)) - (G_spec(4) + G_spec(2))
    dG2 = (G_spec(7) + G_spec(4)) - (G_spec(1) + G_spec(3))
    dG3 = (G_spec(5) + G_spec(4)) - (G_spec(1) + G_spec(7))

    ! === CRITICAL FIX: MAGNITUDES RESTORED TO C++ LEVELS (10^13) ===
    k1f = 3.52d13 * T**(-0.7d0) * exp(-8590.d0 / T)   
    k2f = 5.06d1  * T**( 2.67d0) * exp(-3166.d0 / T)  
    k3f = 1.17d6  * T**( 1.3d0)  * exp(-1829.d0 / T)  
    
    k40 = 5.75d16 * T**(-1.4d0)   
    k4inf = 4.65d9 * T**( 0.44d0) 
    
    k5f = 7.08d10 * exp(-148.d0 / T) 
    k6f = 1.66d10 * exp(-414.d0 / T)
    k7f = 2.89d10 * exp( 250.d0 / T)

    k1r = k1f * exp( dG1 / (R_univ * T) )
    k2r = k2f * exp( dG2 / (R_univ * T) )
    k3r = k3f * exp( dG3 / (R_univ * T) )

    M = 2.5d0*H2 + 16.d0*H2O + O2 + O + H + HO2 + OH + N2

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
      ! Fix: Convert molar rate (kmol/m3/s) to mass rate (kg/m3/s)
      omega_mass = w(k) * Wspec(k)
      Y(k) = max(0.d0, min(1.d0, Y(k) + (omega_mass/rho)*dt))
    end do
    sumY = sum(Y)
    if(sumY>1.d-30) Y=Y/sumY

    ! --- BUG FIX: Multiply by Molecular Weight of H2 (Wspec(1)) ---
    qdot = -DeltaH_H2 * w(1) * Wspec(1) 
    rhoE = rhoE + qdot * dt
    
  end subroutine apply_sources

end module chemistry_mod