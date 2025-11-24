module mod_global
    implicit none
    save

    ! Define double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! --- Simulation Constants ---
    real(dp), parameter :: rho0 = 200.0_dp      ! kg/m^3
    real(dp), parameter :: R_univ = 8.314_dp    ! J/(mol K)

    ! --- Reaction Parameters (3 Steps) ---
    ! Pre-exponential Factors (A)
    real(dp), parameter :: A1 = 4.4e10_dp         ! 1/s
    real(dp), parameter :: A2 = 4.0e10_dp         ! 1/s
    real(dp), parameter :: A3 = 3.0e10_dp         ! cm^3/s-g (Wait! check units)
    ! Note on A3: 30 cm^3/s-g = 30 * (1e-6 m^3) / s / (1e-3 kg) = 0.03 m^3/kg-s
    ! However, based on the rate law r3 = A3 * exp * rho * Yc, 
    ! standard Arrhenius units for 1st order are 1/s. 
    ! If r3 is mass-based, we'll assume the provided value is scaled for the equation provided.
    ! Let's treat the inputs as consistent with the equations in the PDF.

    ! Activation Energies (Given in KJ/mole, converting to J/mole)
    real(dp), parameter :: E1 = 190000.0_dp
    real(dp), parameter :: E2 = 180000.0_dp
    real(dp), parameter :: E3 = 140000.0_dp

    ! Heat of Reaction (J/g -> Convert to J/kg by * 1000)
    real(dp), parameter :: dH1 = 270000.0_dp    ! 270 J/g = 270,000 J/kg
    real(dp), parameter :: dH2 = -800000.0_dp   ! -800 J/g
    real(dp), parameter :: dH3 = -4000000.0_dp  ! -4000 J/g

    ! Specific Heats (J/kg K)
    real(dp), parameter :: Cv_A = 1200.0_dp
    real(dp), parameter :: Cv_B = 1200.0_dp
    real(dp), parameter :: Cv_C = 2000.0_dp
    real(dp), parameter :: Cv_D = 2000.0_dp

end module mod_global