module mod_chem
    use mod_global
    implicit none

contains

    ! --- Helper: Calculate Mixture Cv ---
    function get_cv_mix(Y) result(cv_mix)
        real(dp), dimension(4), intent(in) :: Y ! Y_A, Y_B, Y_C, Y_D
        real(dp) :: cv_mix
        cv_mix = Y(1)*Cv_A + Y(2)*Cv_B + Y(3)*Cv_C + Y(4)*Cv_D
    end function get_cv_mix

    ! --- Helper: Calculate Reaction Rates (r) ---
    subroutine get_rates(T, Y, r)
        real(dp), intent(in) :: T
        real(dp), dimension(4), intent(in) :: Y
        real(dp), dimension(3), intent(out) :: r

        ! Arrhenius rates: r = A * exp(-E/RT) * rho0 * Y
        r(1) = A1 * exp(-E1 / (R_univ * T)) * rho0 * Y(1)
        r(2) = A2 * exp(-E2 / (R_univ * T)) * rho0 * Y(2)
        r(3) = A3 * exp(-E3 / (R_univ * T)) * rho0 * Y(3)
    end subroutine get_rates

    ! --- Get Source Vector (RHS of ODEs) ---
    ! Returns dU/dt = [dYA/dt, dYB/dt, dYC/dt, dYD/dt, dT/dt]
    subroutine get_source_term(T, Y, S)
        real(dp), intent(in) :: T
        real(dp), dimension(4), intent(in) :: Y
        real(dp), dimension(5), intent(out) :: S
        
        real(dp), dimension(3) :: r
        real(dp) :: cv_mix, heat_release

        call get_rates(T, Y, r)
        cv_mix = get_cv_mix(Y)

        ! Species derivatives: dY/dt = Source / rho0
        S(1) = -r(1) / rho0                       ! A
        S(2) = (r(1) - r(2)) / rho0               ! B
        S(3) = (r(2) - r(3)) / rho0               ! C
        S(4) = r(3) / rho0                        ! D

        ! Energy derivative: dT/dt
        heat_release = r(1)*dH1 + r(2)*dH2 + r(3)*dH3
        S(5) = -heat_release / (rho0 * cv_mix)

    end subroutine get_source_term

    ! --- Calculate Analytical Jacobian Matrix (5x5) ---
    ! FIXED: Calculates 'k' directly to avoid division by zero
    subroutine get_jacobian(T, Y, J)
        real(dp), intent(in) :: T
        real(dp), dimension(4), intent(in) :: Y
        real(dp), dimension(5,5), intent(out) :: J

        real(dp), dimension(3) :: r, dr_dT
        real(dp) :: cv, cv2, dcv_dY(4)
        real(dp) :: sum_rdH, sum_drdH_dT, term1
        
        ! NEW: Rate constants k (1/s)
        real(dp) :: k1, k2, k3

        call get_rates(T, Y, r)
        cv = get_cv_mix(Y)
        cv2 = cv * cv
        
        ! 1. Calculate pure Rate Constants (k)
        ! This avoids calculating (r/Y) which is NaN when Y=0
        k1 = A1 * exp(-E1 / (R_univ * T)) 
        k2 = A2 * exp(-E2 / (R_univ * T)) 
        k3 = A3 * exp(-E3 / (R_univ * T)) 
        
        ! 2. Derivatives of Rates w.r.t Temperature
        ! d(r)/dT = r * (E / R T^2)
        dr_dT(1) = r(1) * (E1 / (R_univ * T**2))
        dr_dT(2) = r(2) * (E2 / (R_univ * T**2))
        dr_dT(3) = r(3) * (E3 / (R_univ * T**2))

        ! 3. Derivatives of Cv w.r.t Species
        dcv_dY(1) = Cv_A
        dcv_dY(2) = Cv_B
        dcv_dY(3) = Cv_C
        dcv_dY(4) = Cv_D

        ! --- FILL JACOBIAN ---
        J = 0.0_dp ! Initialize to zero

        ! Rows 1-4: Species Equations (dYi/dt)
        ! Note: Source S = Rate / rho0. 
        ! So dS/dY = (dRate/dY) / rho0 = (k * rho0) / rho0 = k
        
        ! S_A = -r1/rho
        J(1,1) = -k1                    ! d(S_A)/dYA
        J(1,5) = -dr_dT(1) / rho0       ! d(S_A)/dT

        ! S_B = (r1 - r2)/rho
        J(2,1) = k1                     ! d(S_B)/dYA
        J(2,2) = -k2                    ! d(S_B)/dYB
        J(2,5) = (dr_dT(1) - dr_dT(2)) / rho0

        ! S_C = (r2 - r3)/rho
        J(3,2) = k2                     ! d(S_C)/dYB
        J(3,3) = -k3                    ! d(S_C)/dYC
        J(3,5) = (dr_dT(2) - dr_dT(3)) / rho0

        ! S_D = r3/rho
        J(4,3) = k3                     ! d(S_D)/dYC
        J(4,5) = dr_dT(3) / rho0

        ! Row 5: Energy Equation (dT/dt)
        ! S_T = - (Sum r*dH) / (rho * cv)
        sum_rdH = r(1)*dH1 + r(2)*dH2 + r(3)*dH3

        ! Derivative w.r.t Species Y_k
        ! term1 represents d(Sum_rdH)/dYk
        
        ! For Y_A (k=1): d(Sum)/dYA = d(r1*dH1)/dYA = k1 * rho0 * dH1
        term1 = k1 * rho0 * dH1
        J(5,1) = - (cv * term1 - sum_rdH * dcv_dY(1)) / (rho0 * cv2)

        ! For Y_B (k=2)
        term1 = k2 * rho0 * dH2
        J(5,2) = - (cv * term1 - sum_rdH * dcv_dY(2)) / (rho0 * cv2)

        ! For Y_C (k=3)
        term1 = k3 * rho0 * dH3
        J(5,3) = - (cv * term1 - sum_rdH * dcv_dY(3)) / (rho0 * cv2)

        ! For Y_D (k=4) - Rates don't depend on D
        term1 = 0.0_dp
        J(5,4) = - (cv * term1 - sum_rdH * dcv_dY(4)) / (rho0 * cv2)

        ! Derivative w.r.t Temperature T
        ! d(S_T)/dT = - (1/rho*cv) * d(Sum)/dT
        sum_drdH_dT = dr_dT(1)*dH1 + dr_dT(2)*dH2 + dr_dT(3)*dH3
        J(5,5) = - sum_drdH_dT / (rho0 * cv)

    end subroutine get_jacobian

end module mod_chem