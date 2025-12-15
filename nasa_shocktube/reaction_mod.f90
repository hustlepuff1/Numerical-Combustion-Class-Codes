module reaction_mod
    ! Import ALL MWs from constants
    use constants_mod, only: dp, R_univ, MW_H2, MW_O2, MW_H2O, MW_OH, MW_O, MW_H
    implicit none
    
    private
    public :: reaction_subcycle

contains

    subroutine reaction_subcycle(rho, T, Y, dt_total)
        real(dp), intent(in)    :: rho, dt_total
        real(dp), intent(inout) :: T
        real(dp), intent(inout) :: Y(6) 
        
        integer :: step, n_sub
        real(dp) :: dt_sub
        real(dp) :: w_dot(6)
        
        ! Sub-cycling: Divide global step into smaller chemical steps
        ! to handle stiff Arrhenius source terms (10-50x faster than fluid)
        n_sub = 20 
        dt_sub = dt_total / real(n_sub, dp)
        
        do step = 1, n_sub
            call get_production_rates(rho, T, Y, w_dot)
            
            ! Explicit Euler Update: Y_new = Y_old + (w_dot * MW / rho) * dt
            Y(1) = Y(1) + (w_dot(1) * MW_H2 / rho) * dt_sub
            Y(2) = Y(2) + (w_dot(2) * MW_O2 / rho) * dt_sub
            Y(3) = Y(3) + (w_dot(3) * MW_H2O / rho) * dt_sub
            Y(4) = Y(4) + (w_dot(4) * MW_OH / rho) * dt_sub
            Y(5) = Y(5) + (w_dot(5) * MW_O / rho) * dt_sub
            Y(6) = Y(6) + (w_dot(6) * MW_H / rho) * dt_sub
            
            ! Clip negatives (Mass fractions cannot be < 0)
            if(Y(1)<0.0_dp) Y(1)=0.0_dp
            if(Y(2)<0.0_dp) Y(2)=0.0_dp
            if(Y(3)<0.0_dp) Y(3)=0.0_dp
            if(Y(4)<0.0_dp) Y(4)=0.0_dp
            if(Y(5)<0.0_dp) Y(5)=0.0_dp
            if(Y(6)<0.0_dp) Y(6)=0.0_dp
        end do
    end subroutine reaction_subcycle

    subroutine get_production_rates(rho, T, Y, w_dot)
        real(dp), intent(in) :: rho, T, Y(6)
        real(dp), intent(out) :: w_dot(6)
        
        real(dp) :: C(7), kf(8), Rate(8)
        
        ! 1. Calculate Molar Concentrations [kmol/m^3]
        ! C = (rho * Y) / MW
        C(1) = rho * Y(1) / MW_H2   ! H2
        C(2) = rho * Y(2) / MW_O2   ! O2
        C(3) = rho * Y(3) / MW_H2O  ! H2O
        C(4) = rho * Y(4) / MW_OH   ! OH
        C(5) = rho * Y(5) / MW_O    ! O
        C(6) = rho * Y(6) / MW_H    ! H
        C(7) = sum(C(1:6))          ! M (Third Body Concentration)

        ! 2. Forward Rate Constants [SI Units]
        ! Source: Evans & Schexnayder (1980) / Mani et al. (1991)
        ! Formula: k = A * T^n * exp(-Ta/T)
        ! Conversion:
        !   2-Body: A_si = A_cgs * 1.0e-3 (cm3/mol -> m3/kmol)
        !   3-Body: A_si = A_cgs * 1.0e-6 (cm6/mol2 -> m6/kmol2)

        ! [R1] H + O2 <-> OH + O  (Chain Branching)
        ! A_cgs = 2.24e14 -> A_si = 2.24e11
        kf(1) = 2.24e11_dp * exp(-8455.0_dp/T)         

        ! [R2] O + H2 <-> OH + H  (Chain Branching)
        ! A_cgs = 1.74e13 -> A_si = 1.74e10
        kf(2) = 1.74e10_dp * exp(-4756.0_dp/T)         

        ! [R3] H2 + OH <-> H2O + H (Propagation)
        ! A_cgs = 1.55e13 -> A_si = 1.55e10
        kf(3) = 1.55e10_dp * exp(-2667.0_dp/T)         

        ! [R4] OH + OH <-> H2O + O (Propagation)
        ! A_cgs = 6.02e12 -> A_si = 6.02e9
        kf(4) = 6.02e9_dp  * exp(-550.0_dp/T)          
        
        ! [R5] H2 + M <-> H + H + M (Dissociation)
        ! A_cgs = 5.50e18 -> A_si = 5.50e12
        kf(5) = 5.50e12_dp * (T**(-1.0_dp)) * exp(-51987.0_dp/T) 

        ! [R6] O2 + M <-> O + O + M (Dissociation)
        ! A_cgs = 7.20e18 -> A_si = 7.20e12
        kf(6) = 7.20e12_dp * (T**(-1.0_dp)) * exp(-59340.0_dp/T) 

        ! [R7] H2O + M <-> H + OH + M (Dissociation)
        ! A_cgs = 5.20e21 -> A_si = 5.20e15
        kf(7) = 5.20e15_dp * (T**(-1.5_dp)) * exp(-59386.0_dp/T) 

        ! [R8] H + O + M <-> OH + M (Recombination)
        ! A_cgs = 7.10e18 -> A_si = 7.10e12
        kf(8) = 7.10e12_dp * (T**(-1.0_dp))                      

        ! 3. Calculate Reaction Rates [kmol/(m^3 s)]
        Rate(1) = kf(1)*C(6)*C(2)        ! H + O2
        Rate(2) = kf(2)*C(5)*C(1)        ! O + H2
        Rate(3) = kf(3)*C(4)*C(1)        ! OH + H2
        Rate(4) = kf(4)*C(4)*C(4)        ! OH + OH
        Rate(5) = kf(5)*C(1)*C(7)        ! H2 + M
        Rate(6) = kf(6)*C(2)*C(7)        ! O2 + M
        Rate(7) = kf(7)*C(3)*C(7)        ! H2O + M
        Rate(8) = kf(8)*C(6)*C(5)*C(7)   ! H + O + M

        ! 4. Net Production/Destruction per Species [kmol/(m^3 s)]
        ! H2 (Spec 1): Consumed in 2,3,5
        w_dot(1) = -Rate(2) - Rate(3) - Rate(5)
        
        ! O2 (Spec 2): Consumed in 1,6
        w_dot(2) = -Rate(1) - Rate(6)
        
        ! H2O (Spec 3): Produced in 3,4. Consumed in 7.
        w_dot(3) = Rate(3) + Rate(4) - Rate(7)
        
        ! OH (Spec 4): Produced in 1,2,7,8. Consumed in 3,4.
        w_dot(4) = Rate(1) + Rate(2) - Rate(3) - 2.0_dp*Rate(4) + Rate(7) + Rate(8)
        
        ! O (Spec 5): Produced in 1,4,6. Consumed in 2,8.
        w_dot(5) = Rate(1) - Rate(2) + Rate(4) + 2.0_dp*Rate(6) - Rate(8)
        
        ! H (Spec 6): Produced in 2,3,5,7. Consumed in 1,8.
        w_dot(6) = -Rate(1) + Rate(2) + Rate(3) + 2.0_dp*Rate(5) + Rate(7) - Rate(8)
    end subroutine get_production_rates

end module reaction_mod