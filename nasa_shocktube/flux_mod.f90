module flux_mod
    use constants_mod, only: dp, gamma
    implicit none
    
    private
    public :: get_flux_vectors

contains

    !> Computes Inviscid Flux Vectors E (x-flux) and F (y-flux)
    !> Input: Q_cons (Conservative Variables: Rho, RhoU, RhoV, RhoE)
    subroutine get_flux_vectors(Q_cons, E_flux, F_flux)
        real(dp), intent(in)  :: Q_cons(4)
        real(dp), intent(out) :: E_flux(4)
        real(dp), intent(out) :: F_flux(4)
        
        real(dp) :: rho, u, v, E_t, P
        real(dp) :: rhou, rhov
        
        ! 1. Decode Conservative Variables
        rho  = Q_cons(1)
        rhou = Q_cons(2)
        rhov = Q_cons(3)
        E_t  = Q_cons(4)
        
        ! Safety check to prevent divide-by-zero
        if (rho < 1.0e-10_dp) rho = 1.0e-10_dp
        
        u = rhou / rho
        v = rhov / rho
        
        ! 2. Compute Pressure (Equation of State)
        ! P = (gamma - 1) * (Internal Energy)
        P = (gamma - 1.0_dp) * (E_t - 0.5_dp * rho * (u**2 + v**2))
        
        ! 3. Compute Flux Vector E (Flow in X direction)
        ! [ rho*u        ]
        ! [ rho*u^2 + P  ]
        ! [ rho*u*v      ]
        ! [ (Et + P)*u   ]
        E_flux(1) = rhou
        E_flux(2) = rhou * u + P
        E_flux(3) = rhou * v
        E_flux(4) = (E_t + P) * u
        
        ! 4. Compute Flux Vector F (Flow in Y direction)
        ! [ rho*v        ]
        ! [ rho*u*v      ]
        ! [ rho*v^2 + P  ]
        ! [ (Et + P)*v   ]
        F_flux(1) = rhov
        F_flux(2) = rhov * u
        F_flux(3) = rhov * v + P
        F_flux(4) = (E_t + P) * v
        
    end subroutine get_flux_vectors

end module flux_mod