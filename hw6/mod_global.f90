module mod_global
    implicit none
    save

    ! Define double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! --- Simulation Constants ---
    real(dp), parameter :: rho0 = 200.0_dp      ! kg/m^3
    real(dp), parameter :: R_univ = 8.314_dp    ! J/(mol K)

    ! --- Reaction Parameters (High-Speed Case) ---
    ! Pre-exponential Factors (A) in 1/s
    ! Calculated from np.exp(44), np.exp(40), np.exp(30)
    real(dp), parameter :: A1 = 1.285160011435937e19_dp  
    real(dp), parameter :: A2 = 2.353852668370200e17_dp  
    real(dp), parameter :: A3 = 1.068647458152446e13_dp  

    ! Activation Energies (J/mole)
    real(dp), parameter :: E1 = 190000.0_dp
    real(dp), parameter :: E2 = 180000.0_dp
    real(dp), parameter :: E3 = 140000.0_dp

    ! Heat of Reaction (J/kg)
    real(dp), parameter :: dH1 = 270000.0_dp
    real(dp), parameter :: dH2 = -800000.0_dp
    real(dp), parameter :: dH3 = -4000000.0_dp

    ! Specific Heats (J/kg K)
    real(dp), parameter :: Cv_A = 1200.0_dp
    real(dp), parameter :: Cv_B = 1200.0_dp
    real(dp), parameter :: Cv_C = 2000.0_dp
    real(dp), parameter :: Cv_D = 2000.0_dp

end module mod_global