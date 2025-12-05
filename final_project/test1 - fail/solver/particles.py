import numpy as np

class ParticleSolver:
    """
    Lagrangian particle tracker for molten aluminum.
    Phase 3: Handles Injection, Drag, AND Combustion (Hermsen).
    Reference: NASA-OK-199183
    """

    def __init__(self, grid,
                rho_p=2700.0,   # Aluminum density
                Cd=0.44,        # Approx drag coeff for high Re
                d0=150e-6):     # Initial diameter 150 microns
        self.grid = grid
        self.rho_p = rho_p
        self.Cd = Cd

        # Particle state arrays
        self.x = np.array([])      
        self.r = np.array([])      
        self.u = np.array([])      
        self.v = np.array([])      
        self.d = np.array([])      # Diameter
        self.alive = np.array([])  

        self.default_diameter = d0
        
        # Combustion Constants (Hermsen Correlation)
        # d(D^2)/dt = -K
        # K depends on Pressure and Oxidizer conc. 
        # We use a simplified constant K for stability in this project.
        # K ~ 5.0e-7 m^2/s is typical for Aluminized propellants.
        self.burn_rate_constant = 5.0e-7 
        
        # Heat of Combustion (Al -> Al2O3)
        # Approx 31 MJ/kg released into the gas
        self.H_combustion = 3.1e7 

    def inject_from_inlet(self, n_new, u0=0.0, v0=0.0):
        # Inject slightly downstream (i=1) to catch the flow
        i_inj = 1
        x0 = self.grid.X[i_inj, 0]
        r_wall = self.grid.Y[i_inj, -1]

        x_new = np.full(n_new, x0)
        # Random radial distribution
        r_new = np.random.rand(n_new) * r_wall * 0.95

        # Initial kick (gas velocity is high, particles need a start)
        u_new = np.full(n_new, max(u0, 50.0))
        v_new = np.full(n_new, v0)
        d_new = np.full(n_new, self.default_diameter)
        alive_new = np.ones(n_new, dtype=bool)

        self.x = np.concatenate([self.x, x_new]) if self.x.size else x_new
        self.r = np.concatenate([self.r, r_new]) if self.r.size else r_new
        self.u = np.concatenate([self.u, u_new]) if self.u.size else u_new
        self.v = np.concatenate([self.v, v_new]) if self.v.size else v_new
        self.d = np.concatenate([self.d, d_new]) if self.d.size else d_new
        self.alive = np.concatenate([self.alive, alive_new]) if self.alive.size else alive_new

    def _sample_gas(self, gas):
        # Safe interpolation logic (Same as Phase 2)
        ni_gas, nj_gas = gas.u.shape
        X = gas.grid.X
        Y = gas.grid.Y

        i_indices = np.argmin(np.abs(X[:, 0][:, None] - self.x[None, :]), axis=0)
        
        # Approximate J index
        r_wall_local = np.interp(self.x, X[:,0], Y[:,-1])
        eta = self.r / (r_wall_local + 1e-9)
        j_indices = (eta * (nj_gas - 1)).astype(int)

        # Clamp to avoid IndexErrors
        i_indices = np.clip(i_indices, 0, ni_gas - 1)
        j_indices = np.clip(j_indices, 0, nj_gas - 1)

        u_g = gas.u[i_indices, j_indices]
        v_g = gas.v[i_indices, j_indices]
        rho_g = gas.rho[i_indices, j_indices]
        p_g = gas.p[i_indices, j_indices]

        return u_g, v_g, rho_g, p_g, i_indices, j_indices

    def step(self, dt, gas):
        """
        Move particles AND Calculate Coupling Source Terms.
        Returns dictionary of source grids to add to Gas Solver.
        """
        if self.x.size == 0:
            return None

        alive = self.alive
        
        # 1. Sample Gas (Returns properties for ALL particles)
        u_g, v_g, rho_g, p_g, idx_i, idx_j = self._sample_gas(gas)
        
        # 2. Drag & Movement (Vectorized on ALL, then masked update)
        du = u_g - self.u
        dv = v_g - self.v
        mag_rel = np.sqrt(du**2 + dv**2) + 1e-9
        
        Re_p = (rho_g * mag_rel * self.d) / 1.0e-5
        Cd = np.full_like(Re_p, 0.44)
        mask_lam = Re_p < 1000
        Cd[mask_lam] = (24.0 / (Re_p[mask_lam]+1e-9)) * (1.0 + 0.15 * Re_p[mask_lam]**0.687)
        
        mass_p = (np.pi/6.0) * self.rho_p * self.d**3
        drag_force = 0.5 * rho_g * (np.pi * (self.d/2)**2) * Cd * mag_rel
        
        accel_coeff = drag_force / (mass_p + 1e-12)
        
        # Update only ALIVE particles
        self.u[alive] += accel_coeff[alive] * du[alive] * dt
        self.v[alive] += accel_coeff[alive] * dv[alive] * dt
        self.x[alive] += self.u[alive] * dt
        self.r[alive] += self.v[alive] * dt
        
        # 3. COMBUSTION (Pressure Coupled)
        # Calculate K for ALL particles first
        p_psi = p_g * 0.000145038
        p_psi_safe = np.maximum(p_psi, 1.0) 
        K_pressure = 5.0e-7 * (p_psi_safe / 1000.0)**0.27
        
        # --- FIX: Apply [alive] mask to K_pressure so shapes match ---
        d_change = (K_pressure[alive] * dt) / (2.0 * self.d[alive] + 1e-9)
        
        old_mass = mass_p[alive]
        
        # Update Diameter
        self.d[alive] -= d_change
        self.d[alive] = np.maximum(self.d[alive], 1e-7)
        
        new_mass = (np.pi/6.0) * self.rho_p * self.d[alive]**3
        dm = old_mass - new_mass 
        
        # 4. Aggregate Source Terms
        ni_gas, nj_gas = gas.u.shape
        
        S_rho = np.zeros((ni_gas, nj_gas))
        S_en  = np.zeros((ni_gas, nj_gas))
        S_mx  = np.zeros((ni_gas, nj_gas))
        S_my  = np.zeros((ni_gas, nj_gas))
        
        valid = alive
        
        # Mass Source
        mass_rate = dm / dt
        np.add.at(S_rho, (idx_i[valid], idx_j[valid]), mass_rate)
        
        # Energy Source
        heat_rate = mass_rate * self.H_combustion
        np.add.at(S_en, (idx_i[valid], idx_j[valid]), heat_rate)
        
        # Momentum Source
        force_x = mass_p[valid] * (accel_coeff[valid] * du[valid])
        force_y = mass_p[valid] * (accel_coeff[valid] * dv[valid])
        
        np.add.at(S_mx, (idx_i[valid], idx_j[valid]), -force_x)
        np.add.at(S_my, (idx_i[valid], idx_j[valid]), -force_y)
        
        # Normalize
        S_rho /= (gas.vol + 1e-9)
        S_en  /= (gas.vol + 1e-9)
        S_mx  /= (gas.vol + 1e-9)
        S_my  /= (gas.vol + 1e-9)

        # 5. Boundary cleanup
        L = self.grid.length
        r_max_local = np.interp(self.x, self.grid.X[:, 0], self.grid.Y[:, -1])
        out_domain = (self.x < 0) | (self.x > L) | (self.r > r_max_local) | (self.d < 1e-6)
        self.alive[out_domain] = False
        
        return {'rho': S_rho, 'mom_x': S_mx, 'mom_y': S_my, 'energy': S_en}

    def get_alive_particles(self):
        mask = self.alive
        return (self.x[mask], self.r[mask], self.u[mask], self.v[mask], self.d[mask])