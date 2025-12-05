import numpy as np
import matplotlib.pyplot as plt

class Quasi1DNozzle:
    def __init__(self, nx, L, geometry_type='converging-diverging'):
        """
        Initialize the nozzle grid and flow variables.
        """
        self.nx = nx
        self.dx = L / (nx - 1)
        self.x = np.linspace(0, L, nx)
        
        # 1. Grid Generation (Area Profile)
        # Standard test case profile (e.g., Anderson's CFD text)
        self.A = np.zeros(nx)
        if geometry_type == 'converging-diverging':
            # Area = 1 + 2.2(x - 1.5)^2  (Example profile)
            # You should adjust this to match your specific nozzle dimensions
            midpoint = L / 2
            self.A = 1.0 + 2.2 * (self.x - 1.5)**2 
            
        # Geometric source term dA/dx (pre-calculated for efficiency)
        self.dA_dx = np.gradient(self.A, self.dx)
        
        # 2. State Variables (Conservation form: U = [rho*A, rho*u*A, rho*E*A])
        # Using 3 components for Euler; add N_species components for Chemistry later
        self.U = np.zeros((nx, 3)) 
        
        # Primitive variables for easy boundary condition application
        self.rho = np.zeros(nx)
        self.u = np.zeros(nx)
        self.p = np.zeros(nx)
        self.T = np.zeros(nx)
        
        # Physics constants
        self.gamma = 1.4
        self.R = 287.0
        
    def initialize_flow(self, rho_in, T_in, p_exit_ratio=0.1):
        """
        Smart Initialization: Linearly ramps pressure to avoid vacuum creation.
        """
        # 1. Inlet Stagnation Conditions
        self.rho[0] = rho_in
        self.T[0] = T_in
        self.p[0] = rho_in * self.R * T_in
        self.u[0] = 10.0 # Subsonic inlet guess
        
        # 2. Linearly distribute Pressure and Temperature through the nozzle
        # This acts as a 'guide' for the solver
        p_inlet = self.p[0]
        p_exit = p_inlet * p_exit_ratio
        
        # Linear ramp from Inlet Pressure to Exit Pressure
        self.p = np.linspace(p_inlet, p_exit, self.nx)
        
        # Assume T scales roughly with P (Polytropic-ish guess)
        # T / T0 = (P / P0)^((g-1)/g)
        self.T = T_in * (self.p / p_inlet)**((self.gamma - 1) / self.gamma)
        
        # Calculate Density from Ideal Gas Law
        self.rho = self.p / (self.R * self.T)
        
        # 3. Calculate Velocity based on constant Mass Flow guess
        # m_dot = rho * u * A. Let's guess m_dot based on the throat.
        # Assume throat is choked (Mach 1 approx)
        # T_star = T_in * (2 / (gamma + 1))
        # p_star = p_inlet * (2 / (gamma + 1))^(gamma/(gamma-1))
        # rho_star = p_star / (R * T_star)
        # u_star = sqrt(gamma * R * T_star)
        # A_throat = np.min(self.A)
        # m_dot_guess = rho_star * u_star * A_throat
        
        # For simplicity, just use a constant mass flow estimate derived from inlet
        m_dot_guess = 0.6 * self.rho[0] * 340.0 * np.min(self.A) # Rough guess
        
        self.u = m_dot_guess / (self.rho * self.A)
        
        # 4. Populate Conservative Vector U
        self.U[:, 0] = self.rho * self.A
        self.U[:, 1] = self.rho * self.u * self.A
        e = self.p / ((self.gamma - 1) * self.rho)
        E = e + 0.5 * self.u**2
        self.U[:, 2] = self.rho * E * self.A

    def decode_state(self):
        """
        Convert U back to primitives (rho, u, p, T) with SAFETY CHECKS.
        """
        rho_A = self.U[:, 0]
        rho_u_A = self.U[:, 1]
        rho_E_A = self.U[:, 2]
        
        # 1. Recover Density (and floor it to avoid division by zero)
        self.rho = rho_A / self.A
        self.rho = np.maximum(self.rho, 1e-6) # Safety Floor
        
        # 2. Recover Velocity
        self.u = rho_u_A / (self.rho * self.A)
        
        # 3. Recover Pressure
        E = rho_E_A / (self.rho * self.A)
        e = E - 0.5 * self.u**2
        
        # If kinetic energy is too high, internal energy 'e' can become negative numerically.
        # We must floor 'e' or 'p' to prevent sqrt errors.
        self.p = self.rho * (self.gamma - 1) * e
        
        # --- CRITICAL FIX: Pressure Floor ---
        if np.any(self.p <= 0):
            # print("Warning: Negative pressure detected. Correcting...") 
            self.p = np.maximum(self.p, 1e-6)
            
            # OPTIONAL: Consistency correction (update Energy to match new positive Pressure)
            # This keeps the solver stable by adding energy artificially where it crashes.
            e_new = self.p / (self.rho * (self.gamma - 1))
            E_new = e_new + 0.5 * self.u**2
            self.U[:, 2] = self.rho * E_new * self.A

        # 4. Recover Temperature
        self.T = self.p / (self.rho * self.R)

    def compute_flux(self):
        """
        PLACEHOLDER: Implement your Riemann Solver here (Roe, HLLC, or MacCormack).
        Reference Guidelines Section 2: 'Completion of the Cold CFD Solver'.
        """
        flux = np.zeros_like(self.U)
        # TODO: Calculate Flux F = [rho*u*A, (rho*u^2 + p)*A, (rho*E + p)*u*A]
        # TODO: Add artificial viscosity or use an upwind scheme
        return flux

    def compute_source_term(self):
        """
        Calculate Source S = [0, p * dA/dx, 0]
        """
        S = np.zeros_like(self.U)
        S[:, 1] = self.p * self.dA_dx
        return S

    def chemistry_step(self, dt):
        """
        PLACEHOLDER: Implement Operator Splitting here.
        Reference Guidelines Section 3: 'Coupling with a Reactive Flow Routine'.
        """
        # TODO: 
        # 1. Update Temperature/Species based on reaction rates for time dt
        # 2. Update self.U based on new T/Species
        pass

    def time_step(self, dt):
        """
        Main solver loop integration (e.g., RK2 or Euler Explicit)
        """
        # 1. Decode U -> Primitives
        self.decode_state()
        
        # 2. Cold Flow Step
        # F_flux = self.compute_flux()
        # S_term = self.compute_source_term()
        # dU_dt = - (dF/dx) + S
        # self.U += dU_dt * dt
        
        # 3. Reaction Step (Strang Splitting)
        # self.chemistry_step(dt)
        
        # 4. Boundary Conditions
        # Apply inlet/outlet BCs here
        pass

    def compute_flux_roe(self):
        """
        Computes the flux at cell interfaces using the Roe Approximate Riemann Solver.
        Returns:
            flux (numpy array): Flux at cell boundaries. Size (nx-1, 3).
        """
        # 1. Reconstruct Left (L) and Right (R) states at interfaces
        # Simple 1st order: L is node i, R is node i+1
        # Indices: 0 to nx-2 for Left, 1 to nx-1 for Right
        rho_L = self.rho[:-1]
        u_L   = self.u[:-1]
        p_L   = self.p[:-1]
        h_L   = (self.U[:-1, 2] / self.rho[:-1] + self.p[:-1]/self.rho[:-1]) # Total Enthalpy H = E + p/rho
        
        rho_R = self.rho[1:]
        u_R   = self.u[1:]
        p_R   = self.p[1:]
        h_R   = (self.U[1:, 2] / self.rho[1:] + self.p[1:]/self.rho[1:])

        # 2. Roe Averages (The "hat" variables)
        sqrt_rho_L = np.sqrt(rho_L)
        sqrt_rho_R = np.sqrt(rho_R)
        denom = sqrt_rho_L + sqrt_rho_R
        
        u_hat = (sqrt_rho_L * u_L + sqrt_rho_R * u_R) / denom
        h_hat = (sqrt_rho_L * h_L + sqrt_rho_R * h_R) / denom
        
        # Sound speed at the interface (check for negative values during transients)
        a_hat_sq = (self.gamma - 1) * (h_hat - 0.5 * u_hat**2)
        a_hat = np.sqrt(np.maximum(a_hat_sq, 1e-6))

        # 3. Wave Speeds (Eigenvalues)
        lambda1 = u_hat - a_hat
        lambda2 = u_hat
        lambda3 = u_hat + a_hat
        
        # Entropy Fix (Harten's Fix) to prevent expansion shocks
        # If |lambda| < epsilon, replace with (lambda^2 + epsilon^2) / (2*epsilon)
        epsilon = 0.1 * a_hat
        for lmb in [lambda1, lambda2, lambda3]:
            mask = np.abs(lmb) < epsilon
            lmb[mask] = (lmb[mask]**2 + epsilon[mask]**2) / (2 * epsilon[mask])

        # 4. Wave Strengths (Alpha) and Right Eigenvectors (K)
        # Differences across the interface
        drho = rho_R - rho_L
        du   = u_R - u_L
        dp   = p_R - p_L
        
        # Alphas
        alpha2 = (self.gamma - 1) / a_hat**2 * (drho * (h_hat - u_hat**2) + u_hat*du - dp/(self.gamma-1)) # Correction for alpha2 formulation
        # Alternative standard formulation for ideal gas:
        alpha2_alt = drho - dp/a_hat**2
        
        # Let's use the standard diagonalized form for clarity:
        term1 = (dp - rho_L*a_hat*du) / (2*a_hat**2) # This is actually incorrect for Roe average density scaling.
        
        # CORRECT STANDARD ROE ALGORITHM:
        drho = rho_R - rho_L
        dp   = p_R - p_L
        du   = u_R - u_L
        
        alpha1 = (dp - (rho_L * sqrt_rho_R + rho_R * sqrt_rho_L)/denom * a_hat * du) / (2*a_hat**2) # Approximate scaling
        # Actually, using the Roe-averaged density directly for scaling is safer:
        rho_hat = np.sqrt(rho_L * rho_R)
        
        alpha1 = (dp - rho_hat * a_hat * du) / (2 * a_hat**2)
        alpha2 = drho - dp / a_hat**2
        alpha3 = (dp + rho_hat * a_hat * du) / (2 * a_hat**2)

        # Eigenvectors K
        # K1 = [1, u-a, H-ua]^T
        # K2 = [1, u,   0.5u^2]^T  <-- Note: This is for K2 corresponding to entropy wave.
        # K3 = [1, u+a, H+ua]^T
        
        # Vectorized assembly of Dissipation Term: sum(alpha * |lambda| * K)
        # Component 1 (Mass)
        dissip_0 = (np.abs(lambda1) * alpha1 * 1.0 + 
                    np.abs(lambda2) * alpha2 * 1.0 + 
                    np.abs(lambda3) * alpha3 * 1.0)
        
        # Component 2 (Momentum)
        dissip_1 = (np.abs(lambda1) * alpha1 * (u_hat - a_hat) + 
                    np.abs(lambda2) * alpha2 * u_hat + 
                    np.abs(lambda3) * alpha3 * (u_hat + a_hat))
        
        # Component 3 (Energy)
        dissip_2 = (np.abs(lambda1) * alpha1 * (h_hat - u_hat*a_hat) + 
                    np.abs(lambda2) * alpha2 * (0.5 * u_hat**2) + 
                    np.abs(lambda3) * alpha3 * (h_hat + u_hat*a_hat))

        # 5. Average Flux
        # Physical fluxes at L and R nodes
        F_L = np.zeros((self.nx-1, 3))
        F_R = np.zeros((self.nx-1, 3))
        
        # F = [rho*u, rho*u^2 + p, rho*H*u] * A  <-- Wait, A is handled outside or inside?
        # Usually for Quasi-1D, we solve 1D Euler flux then multiply by Area at interface.
        
        # Flux L
        F_L[:,0] = (rho_L * u_L)
        F_L[:,1] = (rho_L * u_L**2 + p_L)
        F_L[:,2] = (rho_L * h_L * u_L)
        
        # Flux R
        F_R[:,0] = (rho_R * u_R)
        F_R[:,1] = (rho_R * u_R**2 + p_R)
        F_R[:,2] = (rho_R * h_R * u_R)
        
        # 6. Final Roe Flux (at interfaces)
        # F_interface = 0.5 * (F_L + F_R) - 0.5 * Dissipation
        
        Flux_bar = np.zeros((self.nx-1, 3))
        Flux_bar[:,0] = 0.5 * (F_L[:,0] + F_R[:,0]) - 0.5 * dissip_0
        Flux_bar[:,1] = 0.5 * (F_L[:,1] + F_R[:,1]) - 0.5 * dissip_1
        Flux_bar[:,2] = 0.5 * (F_L[:,2] + F_R[:,2]) - 0.5 * dissip_2
        
        # Multiply by Area at the interface (midpoint)
        A_interface = 0.5 * (self.A[:-1] + self.A[1:])
        Flux_final = Flux_bar * A_interface[:, np.newaxis]
        
        return Flux_final
    
    def apply_bcs(self):
        """
        Applies boundary conditions for Subsonic Inlet and Supersonic Outlet.
        """
        # --- INLET (Index 0) ---
        # Fixed Reservoir Conditions (Stagnation properties)
        rho_0 = 1.2    # Example stagnation density
        T_0 = 300.0    # Example stagnation temperature
        p_0 = rho_0 * self.R * T_0
        
        # Extrapolate velocity (u) from the first interior point (index 1)
        # (Simple 0th order extrapolation)
        u_inlet = self.u[1] 
        
        # Calculate static T and rho using isentropic relations with fixed T0, P0
        # T0 = T + u^2 / (2*Cp)  -> T = T0 - u^2 / (2*Cp)
        cp = self.gamma * self.R / (self.gamma - 1)
        T_inlet = T_0 - 0.5 * u_inlet**2 / cp
        
        # p / p0 = (T / T0)^(gamma/(gamma-1))
        p_inlet = p_0 * (T_inlet / T_0)**(self.gamma / (self.gamma - 1))
        rho_inlet = p_inlet / (self.R * T_inlet)
        
        # Update Primitive Variables at Inlet
        self.rho[0] = rho_inlet
        self.u[0]   = u_inlet
        self.p[0]   = p_inlet
        self.T[0]   = T_inlet
        
        # Update Conservative Variables U at Inlet
        self.U[0, 0] = self.rho[0] * self.A[0]
        self.U[0, 1] = self.rho[0] * self.u[0] * self.A[0]
        E_inlet = (self.R * self.T[0])/(self.gamma - 1) + 0.5 * self.u[0]**2
        self.U[0, 2] = self.rho[0] * E_inlet * self.A[0]

        # --- OUTLET (Index -1) ---
        # Supersonic outflow: All characteristics go out.
        # Simple Extrapolation (Copy from neighbor)
        self.U[-1, :] = self.U[-2, :]
        self.rho[-1]  = self.rho[-2]
        self.u[-1]    = self.u[-2]
        self.p[-1]    = self.p[-2]
        self.T[-1]    = self.T[-2]


    def time_step(self, CFL=0.5):
        """
        Performs one explicit time step (Forward Euler for simplicity).
        """
        # 1. Calculate time step size (dt) based on CFL condition
        # dt = CFL * dx / (|u| + a)
        a = np.sqrt(self.gamma * self.p / self.rho)
        max_speed = np.max(np.abs(self.u) + a)
        dt = CFL * self.dx / max_speed
        
        # 2. Compute Fluxes at interfaces (Uses the Roe Scheme you added)
        # Flux_final size: (nx-1, 3)
        fluxes = self.compute_flux_roe()
        
        # 3. Compute Source Term (Pressure acting on walls)
        # S = [0, p * dA/dx, 0]
        # We need p at cell centers.
        S = np.zeros_like(self.U)
        S[:, 1] = self.p * self.dA_dx
        
        # 4. Update Interior Cells (Indices 1 to nx-2)
        # dU/dt = - (F_i+1/2 - F_i-1/2) / dx  + S
        # Note on indices: 
        # fluxes[i] is at interface i+1/2. 
        # fluxes[i-1] is at interface i-1/2.
        
        # Flux divergence
        dF = fluxes[1:] - fluxes[:-1] # This matches interior cells 1 to nx-2
        
        # Source term for interior cells
        S_interior = S[1:-1]
        
        # Update
        self.U[1:-1] -= (dt / self.dx) * dF
        self.U[1:-1] += dt * S_interior
        
        # 5. Decode U -> Primitives (rho, u, p) for next step
        self.decode_state()
        
        # 6. Apply Boundary Conditions
        self.apply_bcs()
        
        return dt

if __name__ == "__main__":
    # 1. Setup Simulation
    # 101 points is good for a quick test. L=3.0 meters.
    solver = Quasi1DNozzle(nx=61, L=3.0) 
    
    # Initialize with Reservoir conditions: rho=1.0, T=300 K
    # We initialize the interior somewhat arbitrarily, the BCs will drive the physics.
    solver.initialize_flow(rho_in=1.0, T_in=300.0)
    
    # 2. Simulation Loop
    t = 0.0
    t_end = 0.2  # Enough time for waves to traverse the domain multiple times
    iteration = 0
    max_iter = 5000 # Safety stop
    
    print(f"{'Iter':<10} | {'Time':<10} | {'Max Mach':<10} | {'dt':<10}")
    print("-" * 50)

    try:
        while t < t_end and iteration < max_iter:
            # Run one time step (CFL 0.1 for stability during startup)
            dt = solver.time_step(CFL=0.1)
            t += dt
            iteration += 1
            
            # Monitor progress every 100 steps
            if iteration % 100 == 0:
                # Calculate Mach number for monitoring
                a = np.sqrt(solver.gamma * solver.p / solver.rho)
                mach = solver.u / a
                print(f"{iteration:<10} | {t:<10.4f} | {np.max(mach):<10.4f} | {dt:<10.4e}")

    except KeyboardInterrupt:
        print("\nSimulation stopped by user.")
    except Exception as e:
        print(f"\nError occurred: {e}")
        raise

    print("\nSimulation Complete.")

    # 3. Post-Processing & Validation Plot
    # Calculate final Mach number and Entropy
    a = np.sqrt(solver.gamma * solver.p / solver.rho)
    Mach = solver.u / a
    Entropy = solver.p / (solver.rho ** solver.gamma) # Should be constant in shock-free flow
    
    # Create the visualization (Figure 1 for your project)
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Mach Number (Left Axis)
    color = 'tab:blue'
    ax1.set_xlabel('Position x (m)')
    ax1.set_ylabel('Mach Number', color=color)
    ax1.plot(solver.x, Mach, color=color, linewidth=2, label='CFD Mach')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    
    # Highlight the throat location
    throat_idx = np.argmin(solver.A)
    ax1.axvline(solver.x[throat_idx], color='k', linestyle='--', alpha=0.5, label='Throat')

    # Plot Normalized Pressure (Right Axis) to show expansion
    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('Pressure / P0', color=color)
    P0 = solver.U[0, 2] # Rough estimate of reservoir P based on inlet
    ax2.plot(solver.x, solver.p / solver.p[0], color=color, linestyle='--', linewidth=2, label='Pressure')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title(f'Quasi-1D Nozzle Cold Flow Validation\n(Roe Solver, {iteration} iterations)')
    fig.tight_layout()
    
    # Save the figure for your report
    plt.savefig('nozzle_cold_flow_validation.png', dpi=300)
    print("Plot saved as 'nozzle_cold_flow_validation.png'")
    plt.show()