! ====================================================================
! Module: roe_solver_module
! Purpose: Contains the Roe's Approximate Riemann Solver for the
!          1D Euler equations.
! ====================================================================
module roe_solver_module
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: gamma = 1.4d0

contains

    ! --------------------------------------------------------------------
    ! Subroutine: roe_flux
    ! Purpose: Calculates the numerical flux at an interface using
    !          Roe's approximate Riemann solver.
    ! --------------------------------------------------------------------
    subroutine roe_flux(rho_L, u_L, p_L, rho_R, u_R, p_R, flux)
        implicit none
        ! --- Arguments ---
        real(dp), intent(in)  :: rho_L, u_L, p_L, rho_R, u_R, p_R
        real(dp), intent(out) :: flux(3)

        ! --- Local Variables ---
        ! Conservative and derived states (Left and Right)
        real(dp) :: e_int_L, E_tot_L, H_L, F_L(3)
        real(dp) :: e_int_R, E_tot_R, H_R, F_R(3)
        ! Roe-averaged quantities
        real(dp) :: rho_hat, u_hat, H_hat, a_hat
        ! Eigenvalues (wave speeds)
        real(dp) :: lambda(3)
        ! Eigenvectors
        real(dp) :: e_vec(3,3)
        ! Wave strengths
        real(dp) :: alpha(3)
        ! Jumps in primitive variables
        real(dp) :: drho, du, dp_jump

        ! --- 1. Compute derived quantities for L/R states ---
        ! Left State
        e_int_L = p_L / ((gamma - 1.0_dp) * rho_L)
        E_tot_L = e_int_L + 0.5_dp * u_L**2
        H_L = E_tot_L + p_L / rho_L
        F_L(1) = rho_L * u_L
        F_L(2) = rho_L * u_L**2 + p_L
        F_L(3) = (rho_L * E_tot_L + p_L) * u_L
        
        ! Right State
        e_int_R = p_R / ((gamma - 1.0_dp) * rho_R)
        E_tot_R = e_int_R + 0.5_dp * u_R**2
        H_R = E_tot_R + p_R / rho_R
        F_R(1) = rho_R * u_R
        F_R(2) = rho_R * u_R**2 + p_R
        F_R(3) = (rho_R * E_tot_R + p_R) * u_R
        
        ! --- 2. Compute Roe Averages (Algorithm Step 1) ---
        rho_hat = sqrt(rho_L * rho_R)
        u_hat   = (sqrt(rho_L)*u_L + sqrt(rho_R)*u_R) / (sqrt(rho_L) + sqrt(rho_R))
        H_hat   = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R) / (sqrt(rho_L) + sqrt(rho_R))
        a_hat   = sqrt((gamma - 1.0_dp) * (H_hat - 0.5_dp * u_hat**2))

        ! --- 3. Compute Eigenvalues and Eigenvectors (Algorithm Step 2) ---
        lambda(1) = u_hat - a_hat
        lambda(2) = u_hat
        lambda(3) = u_hat + a_hat

        e_vec(:,1) = [ 1.0_dp, u_hat - a_hat, H_hat - u_hat*a_hat ]
        e_vec(:,2) = [ 1.0_dp, u_hat, 0.5_dp*u_hat**2 ]
        e_vec(:,3) = [ 1.0_dp, u_hat + a_hat, H_hat + u_hat*a_hat ]
        ! Correction for eigenvector 2 (enthalpy term needs to be added)
        e_vec(3,2) = e_vec(3,2) + (a_hat**2)/(gamma - 1.0_dp)

        ! --- 4. Compute Wave Strengths (Algorithm Step 3) ---
        drho    = rho_R - rho_L
        du      = u_R - u_L
        dp_jump = p_R - p_L

        alpha(2) = drho - dp_jump / a_hat**2
        alpha(1) = (dp_jump - rho_hat * a_hat * du) / (2.0_dp * a_hat**2)
        alpha(3) = (dp_jump + rho_hat * a_hat * du) / (2.0_dp * a_hat**2)

        ! --- 5. Compute the Interface Flux (Algorithm Step 4) ---
        ! Form: F_interface = 0.5*(F_L + F_R) - 0.5*Sum( |lambda|*alpha*e_vec )
        flux = 0.5_dp * (F_L + F_R) - 0.5_dp * ( &
               abs(lambda(1))*alpha(1)*e_vec(:,1) + &
               abs(lambda(2))*alpha(2)*e_vec(:,2) + &
               abs(lambda(3))*alpha(3)*e_vec(:,3) )

    end subroutine roe_flux

end module roe_solver_module
