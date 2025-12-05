import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# ==========================================
# 1. CONFIGURATION
# ==========================================
class Config:
    # Grid
    N_POINTS = 200   
    L_NOZZLE = 0.5   
    THROAT_LOC = 0.2 
    
    # Time stepping
    CFL = 0.4        
    T_END = 0.015    # 15 ms to allow full startup
    
    # Ramps
    RAMP_TIME_PARTICLES = 0.005 # Wait 5ms before adding particles
    
    # Gas Properties (CEA)
    GAMMA = 1.125
    MW = 26.17  
    R_UNIV = 8314.46
    R_GAS = R_UNIV / MW  
    
    # Chamber (Stagnation) Conditions
    P_0 = 64.0 * 1e5  
    T_0 = 3729.0      
    RHO_0 = P_0 / (R_GAS * T_0)
    
    # Ambient
    P_AMB = 1.0 * 1e5
    
    # Particle Properties
    RHO_PARTICLE = 2700.0   
    D_INITIAL = 150e-6      
    BURNING_RATE_CONST = 1.5e-7 

# ==========================================
# 2. GEOMETRY
# ==========================================
def generate_nozzle_geometry(n, L, throat_x):
    x = np.linspace(0, L, n)
    
    r_inlet = np.sqrt(0.05 / np.pi)
    r_throat = np.sqrt(0.01 / np.pi)
    r_exit = np.sqrt(0.08 / np.pi)
    
    R = np.zeros_like(x)
    
    for i, xi in enumerate(x):
        if xi < throat_x:
            theta = np.pi * (xi / throat_x)
            R[i] = r_throat + (r_inlet - r_throat) * 0.5 * (1 + np.cos(theta))
        else:
            theta = np.pi * (xi - throat_x) / (L - throat_x)
            R[i] = r_throat + (r_exit - r_throat) * 0.5 * (1 - np.cos(theta))
            
    A = np.pi * R**2
    return x, A

# ==========================================
# 3. GAS SOLVER
# ==========================================
class GasSolver:
    def __init__(self, x, A, cfg):
        self.x = x
        self.A = A
        self.dx = x[1] - x[0]
        self.cfg = cfg
        self.n = len(x)
        
        print("Initializing with Linear Pressure Gradient (Engine Start)...")
        # Linear Pressure Drop: High in chamber, Low at exit
        # This drives the flow naturally
        self.P = np.linspace(cfg.P_0, cfg.P_AMB, self.n)
        
        # Assume Temperature drops slightly (Adiabatic-ish guess)
        self.T = np.linspace(cfg.T_0, cfg.T_0 * 0.6, self.n)
        
        # Density from State Eq
        self.rho = self.P / (cfg.R_GAS * self.T)
        
        # Low initial velocity (Let pressure gradient do the work)
        self.u = np.linspace(10.0, 100.0, self.n)
        
        self.U = np.zeros((3, self.n))
        self.update_conservative()

    def get_energy(self, T, u):
        Cv = self.cfg.R_GAS / (self.cfg.GAMMA - 1)
        return Cv * T + 0.5 * u**2

    def update_conservative(self):
        E = self.get_energy(self.T, self.u)
        self.U[0, :] = self.rho * self.A
        self.U[1, :] = self.rho * self.u * self.A
        self.U[2, :] = self.rho * E * self.A

    def enforce_physicality(self, U_state):
        eps = 1e-16
        min_rho = 1e-4
        min_P = 100.0
        
        # Clamp Density
        U_state[0, :] = np.maximum(U_state[0, :], min_rho * self.A)
        
        rho = U_state[0, :] / self.A
        rho_safe = rho + eps
        
        # Clamp Velocity (Prevent runaway during transient)
        u = U_state[1, :] / (U_state[0, :] + eps)
        max_u = 5000.0
        if np.any(np.abs(u) > max_u):
            u = np.clip(u, -max_u, max_u)
            U_state[1, :] = u * U_state[0, :]

        # Clamp Internal Energy
        E_tot = U_state[2, :] / rho_safe
        ke = 0.5 * u**2
        e_int = E_tot - ke
        
        min_e_int = min_P / ((self.cfg.GAMMA - 1) * rho_safe)
        
        if np.any(e_int < min_e_int):
            # If internal energy is non-physical, heat it up to minimum
            e_int = np.maximum(e_int, min_e_int)
            U_state[2, :] = rho * (e_int + ke) * self.A
            
        return U_state

    def decode_conservative(self):
        self.U = self.enforce_physicality(self.U)
        eps = 1e-16
        
        self.rho = self.U[0, :] / self.A
        self.u = self.U[1, :] / (self.U[0, :] + eps)
        E = self.U[2, :] / (self.U[0, :] + eps)
        
        Cv = self.cfg.R_GAS / (self.cfg.GAMMA - 1)
        self.T = (E - 0.5 * self.u**2) / Cv
        self.P = self.rho * self.cfg.R_GAS * self.T

    def apply_boundary_conditions(self, U_vec):
        U_vec = self.enforce_physicality(U_vec)
        
        # --- INLET (Subsonic Stagnation) ---
        # Fix Stagnation P0, T0
        # U[0] (Density) is extrapolated/calculated from P0, T0 and u_extrap
        
        u_inner = U_vec[1, 1] / U_vec[0, 1]
        u_in = max(0.0, u_inner)
        
        # Isentropic relation for static state from stagnation
        T_in = self.cfg.T_0 - (self.cfg.GAMMA - 1)/(2*self.cfg.GAMMA * self.cfg.R_GAS) * u_in**2
        # Clamp T_in
        T_in = max(T_in, 100.0)
        
        P_in = self.cfg.P_0 * (T_in / self.cfg.T_0)**(self.cfg.GAMMA/(self.cfg.GAMMA-1))
        rho_in = P_in / (self.cfg.R_GAS * T_in)
        
        U_vec[0, 0] = rho_in * self.A[0]
        U_vec[1, 0] = rho_in * u_in * self.A[0]
        E_in = self.get_energy(T_in, u_in)
        U_vec[2, 0] = rho_in * E_in * self.A[0]

        # --- OUTLET (Supersonic Extrapolation) ---
        # We assume the nozzle starts and stays supersonic at exit
        U_vec[:, -1] = 2*U_vec[:, -2] - U_vec[:, -3]
            
        return U_vec

    def step(self, dt, source_terms, ramp_factor):
        U_old = self.U.copy()
        eps = 1e-16
        
        # 1. Flux Function
        def get_flux(U_in):
            rhoA = U_in[0, :]
            rhouA = U_in[1, :]
            rhoEA = U_in[2, :]
            
            rho = rhoA / self.A
            u = rhouA / (rhoA + eps)
            
            # P = (g-1)*(E - 0.5*rho*u^2)
            # PA = (g-1)*(rhoEA - 0.5 * (rhouA)^2 / rhoA)
            PA = (self.cfg.GAMMA - 1) * (rhoEA - 0.5 * rhouA**2 / (rhoA + eps))
            
            F = np.zeros_like(U_in)
            F[0, :] = rhouA
            F[1, :] = rhouA**2 / (rhoA + eps) + PA
            F[2, :] = (rhoEA + PA) * u
            return F, PA / self.A

        # 2. Artificial Viscosity (Jameson-style)
        # Pressure sensor
        P = self.P
        nu_k = np.abs(np.diff(P, append=P[-1]) - np.diff(P, prepend=P[0])) / (4*P + eps)
        nu_k = np.clip(nu_k, 0, 1.0)
        # Adaptive viscosity: High near shocks, low elsewhere
        visc_coeff = 0.5 # Base
        
        smoothing = np.zeros_like(self.U)
        smoothing[:, 1:-1] = (self.U[:, 2:] - 2*self.U[:, 1:-1] + self.U[:, :-2])
        S_visc = visc_coeff * smoothing

        # 3. Predictor
        F, P_val = get_flux(self.U)
        
        # Source term H (Pressure Area force)
        H = np.zeros_like(self.U)
        dA_dx = np.gradient(self.A, self.x)
        H[1, :] = P_val * dA_dx
        
        dFdx = np.zeros_like(self.U)
        dFdx[:, :-1] = (F[:, 1:] - F[:, :-1]) / self.dx
        
        U_bar = U_old - dt * (dFdx - H - source_terms*ramp_factor) + S_visc
        U_bar = self.apply_boundary_conditions(U_bar)

        # 4. Corrector
        F_bar, P_bar = get_flux(U_bar)
        
        H_bar = np.zeros_like(H)
        H_bar[1, :] = P_bar * dA_dx
        
        dFdx_c = np.zeros_like(self.U)
        dFdx_c[:, 1:] = (F_bar[:, 1:] - F_bar[:, :-1]) / self.dx
        
        self.U = 0.5 * (U_old + U_bar - dt * (dFdx_c - H_bar - source_terms*ramp_factor)) + S_visc
        self.U = self.apply_boundary_conditions(self.U)
        
        self.decode_conservative()

# ==========================================
# 4. PARTICLES
# ==========================================
class ParticleCloud:
    def __init__(self, cfg):
        self.cfg = cfg
        self.particles = [] 
        
    def inject(self, dt):
        n_inject = 10 
        x_inj = 2.0 * (0.5/200)
        for _ in range(n_inject):
            p = {
                'x': x_inj,
                'u': 50.0, 
                'T': self.cfg.T_0,
                'D': self.cfg.D_INITIAL,
                'm': (4/3) * np.pi * (self.cfg.D_INITIAL/2)**3 * self.cfg.RHO_PARTICLE
            }
            self.particles.append(p)

    def update(self, dt, gas_solver):
        S_gas = np.zeros_like(gas_solver.U)
        alive_particles = []
        
        for p in self.particles:
            if p['x'] >= gas_solver.x[-2] or p['x'] < 0:
                continue 
            
            idx = int(p['x'] / gas_solver.dx)
            idx = min(idx, gas_solver.n - 1)
            
            u_g = gas_solver.u[idx]
            T_g = gas_solver.T[idx]
            rho_g = gas_solver.rho[idx]

            # 1. Mass
            dm_dt_gas_add = 0.0
            if p['D'] > 50e-6:
                dD_dt = - self.cfg.BURNING_RATE_CONST / p['D'] 
                p['D'] += dD_dt * dt
                if p['D'] < 1e-7: p['D'] = 1e-7
                
                m_new = (4/3) * np.pi * (p['D']/2)**3 * self.cfg.RHO_PARTICLE
                dm_loss = (m_new - p['m']) / dt 
                p['m'] = m_new
                dm_dt_gas_add = -dm_loss 
                S_gas[0, idx] -= dm_dt_gas_add / gas_solver.dx

            # 2. Momentum
            rel_vel = abs(u_g - p['u'])
            Re = rho_g * rel_vel * p['D'] / 1.8e-5
            Cd = 24.0/(Re+1e-9) + 0.44 if Re > 0.1 else 240.0
            tau_d = (4 * self.cfg.RHO_PARTICLE * p['D']) / (3 * Cd * rho_g * rel_vel + 1e-9)
            du_dt = (u_g - p['u']) / tau_d
            p['u'] += du_dt * dt
            p['x'] += p['u'] * dt
            S_gas[1, idx] -= p['m'] * du_dt / gas_solver.dx

            # 3. Energy
            Nu = 2.0 + 0.6 * np.sqrt(Re) * (0.7)**(1/3)
            h = Nu * 0.1 / p['D']
            dT_dt = (h * 6 / (self.cfg.RHO_PARTICLE * p['D'] * 1000)) * (T_g - p['T'])
            p['T'] += dT_dt * dt
            q_dot = -p['m'] * 1000 * dT_dt 
            
            Cp_gas = gas_solver.cfg.GAMMA * gas_solver.cfg.R_GAS / (gas_solver.cfg.GAMMA - 1)
            h_gas_total = Cp_gas * T_g + 0.5 * u_g**2
            energy_from_mass = dm_dt_gas_add * h_gas_total
            
            S_gas[2, idx] -= (q_dot + energy_from_mass) / gas_solver.dx

            alive_particles.append(p)
            
        self.particles = alive_particles
        return S_gas

# ==========================================
# 5. MAIN
# ==========================================
def main():
    cfg = Config()
    x, A = generate_nozzle_geometry(cfg.N_POINTS, cfg.L_NOZZLE, cfg.THROAT_LOC)
    
    gas = GasSolver(x, A, cfg)
    particles = ParticleCloud(cfg)
    
    t = 0.0
    dt = 1e-8 
    
    print(f"Starting Simulation with Linear Init (Engine Start)...")
    
    while t < cfg.T_END:
        # 1. Particles (Delayed)
        source_terms = np.zeros_like(gas.U)
        ramp = 0.0
        
        if t > cfg.RAMP_TIME_PARTICLES:
            if int(t/1e-5) % 20 == 0: 
                particles.inject(dt)
            source_terms = particles.update(dt, gas)
            ramp = min((t - cfg.RAMP_TIME_PARTICLES)/0.001, 1.0)
        
        # 2. Time Step
        a = np.sqrt(cfg.GAMMA * np.abs(gas.P) / (gas.rho + 1e-16))
        max_wave = np.max(np.abs(gas.u) + a)
        
        # Safety check
        if max_wave > 1e5 or np.any(np.isnan(gas.U)):
            print(f"Instability detected at t={t*1000:.3f}ms")
            break
        
        dt_cfl = cfg.CFL * np.min(gas.dx / (max_wave + 1e-9))
        dt = min(dt_cfl, 1e-6)
        if t < 0.001: dt = min(dt, 1e-7) # Slow start
        
        # 3. Step
        gas.step(dt, source_terms, ramp)
        t += dt
        
        if int(t/dt) % 2000 == 0:
            p_ex = gas.P[-1]/1e5
            mach_ex = gas.u[-1]/a[-1]
            print(f"Time: {t*1000:.2f} ms | MachExit: {mach_ex:.2f} | P_Exit: {p_ex:.2f} bar")

    # ==========================================
    # 6. PLOTS
    # ==========================================
    print("Generating Plots...")
    fig = plt.figure(figsize=(16, 12))
    
    # 1D Plots
    mach_1d = gas.u / np.sqrt(cfg.GAMMA * gas.P / gas.rho)
    
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(x, mach_1d, 'b-')
    ax1.axvline(cfg.THROAT_LOC, color='k', linestyle='--', label='Throat')
    ax1.set_ylabel('Mach'); ax1.legend(); ax1.grid()
    ax1.set_title('Mach Number')
    
    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(x, gas.P/1e5, 'r-')
    ax2.set_ylabel('Pressure (bar)'); ax2.grid()
    ax2.set_title('Pressure')
    
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(x, gas.T, 'orange')
    ax3.set_ylabel('Temp (K)'); ax3.grid()
    ax3.set_title('Temperature')
    
    # 2D Contour
    Ny = 50
    Radius = np.sqrt(A / np.pi)
    X_grid = np.zeros((2*Ny, cfg.N_POINTS))
    Y_grid = np.zeros((2*Ny, cfg.N_POINTS))
    Mach_grid = np.zeros((2*Ny, cfg.N_POINTS))
    
    for i in range(cfg.N_POINTS):
        r_wall = Radius[i]
        y_col = np.linspace(-r_wall, r_wall, 2*Ny)
        X_grid[:, i] = x[i]
        Y_grid[:, i] = y_col
        Mach_grid[:, i] = mach_1d[i] 
        
    ax4 = plt.subplot(2, 2, 4)
    levels = np.linspace(0, np.nanmax(mach_1d)*1.05, 30)
    cp = ax4.contourf(X_grid, Y_grid, Mach_grid, levels=levels, cmap='jet')
    fig.colorbar(cp, ax=ax4, label='Mach')
    ax4.plot(x, Radius, 'k-', linewidth=2)
    ax4.plot(x, -Radius, 'k-', linewidth=2)
    
    if len(particles.particles) > 0:
        px = [p['x'] for p in particles.particles[::10]]
        py = []
        for val_x in px:
            idx = int(val_x / gas.dx)
            idx = min(idx, cfg.N_POINTS-1)
            r_local = Radius[idx]
            py.append(np.random.uniform(-r_local*0.9, r_local*0.9))
        ax4.scatter(px, py, s=3, c='white', alpha=0.7)
        
    ax4.set_title('Mach Contour')
    ax4.axis('equal')
    
    plt.tight_layout()
    plt.savefig('srm_cfd_results.png')
    print("Done.")

if __name__ == "__main__":
    main()