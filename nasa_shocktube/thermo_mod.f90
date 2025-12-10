module thermo_mod
    use constants_mod, only: dp, R_univ, MW_H2, MW_O2
    implicit none
    
    private
    public :: get_temperature_from_energy, get_species_cp_R, get_species_enthalpy_molar

    ! --- NASA 7-Coefficient Polynomials (HIGH T range: 1000K - 3500K) ---
    ! Source: GRI-Mech 3.0 (thermo30.dat) provided by user
    ! NOTE: We use High-T because our flow is >1000K (Freestream 1388K)
    
    ! H2 Coefficients (High T)
    real(dp), parameter :: H2_coeffs(7) = [ &
        3.33727920E+00_dp, -4.94024731E-05_dp, 4.99456778E-07_dp, &
        -1.79566394E-10_dp, 2.00255376E-14_dp, -9.50158922E+02_dp, &
        -3.20502331E+00_dp ]

    ! O2 Coefficients (High T)
    real(dp), parameter :: O2_coeffs(7) = [ &
        3.28253784E+00_dp, 1.48308754E-03_dp, -7.57966669E-07_dp, &
        2.09470555E-10_dp, -2.16717794E-14_dp, -1.08845772E+03_dp, &
        5.45323129E+00_dp ]

    ! N2 Coefficients (High T) - Optional background
    real(dp), parameter :: N2_coeffs(7) = [ &
        2.92664000E+00_dp, 1.48797680E-03_dp, -5.68476000E-07_dp, &
        1.00970380E-10_dp, -6.75335100E-14_dp, -9.22797700E+02_dp, &
        5.98052800E+00_dp ]
        
    ! OH Coefficients (High T) - Needed for Phase 9
    real(dp), parameter :: OH_coeffs(7) = [ &
        3.09288767E+00_dp, 5.48429716E-04_dp, 1.26505228E-07_dp, &
        -8.79461556E-11_dp, 1.17412376E-14_dp, 3.85865700E+03_dp, &
        4.47669610E+00_dp ]

    ! H2O Coefficients (High T) - Needed for Phase 9
    real(dp), parameter :: H2O_coeffs(7) = [ &
        3.03399249E+00_dp, 2.17691804E-03_dp, -1.64072518E-07_dp, &
        -9.70419870E-11_dp, 1.68200992E-14_dp, -3.00042971E+04_dp, &
        4.96677010E+00_dp ]

contains

    !> Calculate Cp/R (Dimensionless)
    function get_species_cp_R(T, coeffs) result(cp_R)
        real(dp), intent(in) :: T
        real(dp), intent(in) :: coeffs(7)
        real(dp) :: cp_R
        
        ! Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        cp_R = coeffs(1) + coeffs(2)*T + coeffs(3)*T**2 + &
               coeffs(4)*T**3 + coeffs(5)*T**4
    end function get_species_cp_R

    !> Calculate Enthalpy (h) in J/kmol
    !> h/RT = a1 + a2*T/2 + ...
    function get_species_enthalpy_molar(T, coeffs) result(h_molar)
        real(dp), intent(in) :: T
        real(dp), intent(in) :: coeffs(7)
        real(dp) :: h_molar, h_RT
        
        h_RT = coeffs(1) + coeffs(2)*T/2.0_dp + coeffs(3)*T**2/3.0_dp + &
               coeffs(4)*T**3/4.0_dp + coeffs(5)*T**4/5.0_dp + coeffs(6)/T
               
        h_molar = h_RT * R_univ * T
    end function get_species_enthalpy_molar

    !> Calculate Mixture Temperature from Internal Energy (e)
    !> Uses Newton-Raphson Iteration because e(T) is non-linear
    subroutine get_temperature_from_energy(e_internal, Y_H2, Y_O2, T_guess, T_out)
        real(dp), intent(in) :: e_internal  ! J/kg
        real(dp), intent(in) :: Y_H2, Y_O2
        real(dp), intent(in) :: T_guess
        real(dp), intent(out) :: T_out
        
        integer :: iter
        real(dp) :: T_new, T_old, e_calc, Cv_mix
        real(dp) :: h_H2, h_O2, e_H2, e_O2
        real(dp) :: R_spec_H2, R_spec_O2
        real(dp) :: Cp_H2_molar, Cp_O2_molar
        
        R_spec_H2 = R_univ / MW_H2
        R_spec_O2 = R_univ / MW_O2
        
        T_old = T_guess
        
        ! Newton-Raphson
        do iter = 1, 15
            ! Enthalpy (J/kmol) -> convert to J/kg
            h_H2 = get_species_enthalpy_molar(T_old, H2_coeffs) / MW_H2 
            h_O2 = get_species_enthalpy_molar(T_old, O2_coeffs) / MW_O2
            
            ! Internal Energy e = h - RT
            e_H2 = h_H2 - R_spec_H2 * T_old
            e_O2 = h_O2 - R_spec_O2 * T_old
            
            e_calc = Y_H2 * e_H2 + Y_O2 * e_O2
            
            ! Cv = Cp - R
            ! Get Cp in J/kmol-K first
            Cp_H2_molar = get_species_cp_R(T_old, H2_coeffs) * R_univ
            Cp_O2_molar = get_species_cp_R(T_old, O2_coeffs) * R_univ
            
            ! Convert to J/kg-K and subtract R_spec
            Cv_mix = Y_H2 * (Cp_H2_molar/MW_H2 - R_spec_H2) + &
                     Y_O2 * (Cp_O2_molar/MW_O2 - R_spec_O2)
                     
            ! Update T
            if (abs(Cv_mix) < 1.0e-10_dp) Cv_mix = 1.0e-10_dp 
            T_new = T_old - (e_calc - e_internal) / Cv_mix
            
            if (abs(T_new - T_old) < 1.0e-3_dp) exit
            T_old = T_new
        end do
        
        T_out = T_new
    end subroutine get_temperature_from_energy

end module thermo_mod