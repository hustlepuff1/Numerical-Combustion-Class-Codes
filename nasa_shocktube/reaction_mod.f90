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
        ! REMOVED: unused T_old, S_Y
        
        n_sub = 20 
        dt_sub = dt_total / real(n_sub, dp)
        
        do step = 1, n_sub
            call get_production_rates(rho, T, Y, w_dot)
            
            ! Explicit Euler Update
            Y(1) = Y(1) + (w_dot(1) * MW_H2 / rho) * dt_sub
            Y(2) = Y(2) + (w_dot(2) * MW_O2 / rho) * dt_sub
            Y(3) = Y(3) + (w_dot(3) * MW_H2O / rho) * dt_sub
            Y(4) = Y(4) + (w_dot(4) * MW_OH / rho) * dt_sub
            Y(5) = Y(5) + (w_dot(5) * MW_O / rho) * dt_sub
            Y(6) = Y(6) + (w_dot(6) * MW_H / rho) * dt_sub
            
            ! Clip negatives
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
        
        C(1) = rho * Y(1) / MW_H2
        C(2) = rho * Y(2) / MW_O2
        C(3) = rho * Y(3) / MW_H2O
        C(4) = rho * Y(4) / MW_OH
        C(5) = rho * Y(5) / MW_O
        C(6) = rho * Y(6) / MW_H
        C(7) = sum(C(1:6)) 

        ! Rates (SI Units: m, kmol, s)
        kf(1) = 2.24e11_dp * exp(-8455.0_dp/T)         
        kf(2) = 1.74e10_dp * exp(-4756.0_dp/T)         
        kf(3) = 1.55e10_dp * exp(-2667.0_dp/T)         
        kf(4) = 6.02e9_dp  * exp(-550.0_dp/T)          
        
        kf(5) = 5.50e12_dp * (T**(-1.0_dp)) * exp(-51987.0_dp/T) 
        kf(6) = 7.20e12_dp * (T**(-1.0_dp)) * exp(-59340.0_dp/T) 
        kf(7) = 5.20e15_dp * (T**(-1.5_dp)) * exp(-59386.0_dp/T) 
        kf(8) = 7.10e12_dp * (T**(-1.0_dp))                      

        Rate(1) = kf(1)*C(6)*C(2)
        Rate(2) = kf(2)*C(5)*C(1)
        Rate(3) = kf(3)*C(4)*C(1)
        Rate(4) = kf(4)*C(4)*C(4)
        Rate(5) = kf(5)*C(1)*C(7)
        Rate(6) = kf(6)*C(2)*C(7)
        Rate(7) = kf(7)*C(3)*C(7)
        Rate(8) = kf(8)*C(6)*C(5)*C(7)

        w_dot(1) = -Rate(2) - Rate(3) - Rate(5)
        w_dot(2) = -Rate(1) - Rate(6)
        w_dot(3) = Rate(3) + Rate(4) - Rate(7)
        w_dot(4) = Rate(1) + Rate(2) - Rate(3) - 2.0_dp*Rate(4) + Rate(7) + Rate(8)
        w_dot(5) = Rate(1) - Rate(2) + Rate(4) + 2.0_dp*Rate(6) - Rate(8)
        w_dot(6) = -Rate(1) + Rate(2) + Rate(3) + 2.0_dp*Rate(5) + Rate(7) - Rate(8)
    end subroutine get_production_rates

end module reaction_mod