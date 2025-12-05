import numpy as np
import matplotlib.pyplot as plt
import time

# ==========================================
# 1. CONFIGURATION
# ==========================================
class Config:
    # Grid
    NX = 200          # Higher resolution for better shock
    NY = 100
    L_X = 0.5         # Length [m]
    L_Y = 0.2         # Height [m]
    WEDGE_ANGLE = 15.0 # Degrees
    
    # Physics
    GAMMA = 1.4
    R_GAS = 287.0
    
    # Simulation
    CFL = 0.1         # Stability factor (Keep < 1.0)
    T_FINAL = 2.0e-4  # End time [s]
    
    # Chemistry (1-Step Westbrook-Dryer Model)
    # This approximates the heat release of the full 7-step model
    # but is much more stable for a student project.
    A_PRE = 2.0e9    # Pre-exponential factor
    TA = 12000.0     # Activation Temp [K]
    Q_HEAT = 4.0e6   # Heat Release [J/kg]

# ==========================================
# 2. GRID GENERATION
# ==========================================
def generate_grid(cfg):
    x = np.linspace(0, cfg.L_X, cfg.NX)
    y_base = np.linspace(0, cfg.L_Y, cfg.NY)
    X, Y = np.meshgrid(x, y_base, indexing='ij')
    
    # Wedge Geometry
    wedge_start = 0.05
    angle_rad = np.deg2rad(cfg.WEDGE_ANGLE)
    
    for i in range(cfg.NX):
        if X[i, 0] > wedge_start:
            y_wall = (X[i, 0] - wedge_start) * np.tan(angle_rad)
            # Compress grid above wedge
            Y[i, :] = y_wall + (Y[i, :] / cfg.L_Y) * (cfg.L_Y - y_wall)
            
    # Calculate Cell Sizes (Simple approx for now)
    dX = np.diff(X, axis=0)
    dY = np.diff(Y, axis=1)
    # Pad to match shape
    dx_grid = np.pad(dX, ((0,1),(0,0)), 'edge')
    dy_grid = np.pad(dY, ((0,0),(0,1)), 'edge')
    
    return X, Y, dx_grid, dy_grid

# ==========================================
# 3. SOLVER CLASS (Manual RK3)
# ==========================================
class ODWESolver:
    def __init__(self, cfg):
        self.cfg = cfg
        # 1. Generate Grid
        self.X, self.Y, _, _ = generate_grid(cfg)
        self.nx, self.ny = self.X.shape
        
        # 2. State Vector U
        self.U = np.zeros((5, self.nx, self.ny))
        
        # 3. ROBUST METRIC CALCULATION (Analytic)
        # Instead of np.gradient, we calculate metrics based on the known wedge geometry.
        # This prevents numerical noise at the sharp corner.
        
        # Grid spacings in Computational Space (assumed 1 for simplicity)
        d_xi = 1.0
        d_eta = 1.0
        
        # x_xi (dx/di): Uniform spacing in X
        x_xi = (cfg.L_X / (cfg.NX - 1)) * np.ones_like(self.X)
        x_eta = np.zeros_like(self.X) # Grid lines are vertical
        
        # y_xi (dy/di): This is the slope of the grid lines
        # y_eta (dy/dj): This is the local cell height
        
        wedge_start = 0.05
        angle_rad = np.deg2rad(cfg.WEDGE_ANGLE)
        
        y_xi = np.zeros_like(self.Y)
        y_eta = np.zeros_like(self.Y)
        
        for i in range(self.nx):
            # Calculate Wall Height and Slope at this i
            x_loc = self.X[i, 0]
            if x_loc > wedge_start:
                h_wall = (x_loc - wedge_start) * np.tan(angle_rad)
                dh_wall_dx = np.tan(angle_rad) # Exact slope
            else:
                h_wall = 0.0
                dh_wall_dx = 0.0
            
            # Y = h_wall + (j/ny) * (L_Y - h_wall)
            # Differentiate w.r.t x (which is xi)
            # dy/dx = dh/dx + (j/ny) * (0 - dh/dx)
            # dy/dx = dh/dx * (1 - Y_normalized)
            
            # Normalized Y coordinate (0 at wall, 1 at top)
            y_norm = (self.Y[i, :] - h_wall) / (cfg.L_Y - h_wall)
            
            y_xi[i, :] = dh_wall_dx * (1.0 - y_norm) * x_xi[i, :] 
            
            # Differentiate w.r.t eta (j)
            # dy/d_eta = (L_Y - h_wall) / (NY - 1)
            y_eta[i, :] = (cfg.L_Y - h_wall) / (cfg.NY - 1)

        # Jacobian and Inverse Metrics
        self.J = 1.0 / (x_xi * y_eta - x_eta * y_xi)
        
        self.xi_x = self.J * y_eta
        self.xi_y = -self.J * x_eta
        self.eta_x = -self.J * y_xi
        self.eta_y = self.J * x_xi
        
    def initialize(self):
        # Mach 5 Inlet
        T_inf = 300.0
        p_inf = 101325.0 
        rho_inf = p_inf / (self.cfg.R_GAS * T_inf)
        a_inf = np.sqrt(self.cfg.GAMMA * self.cfg.R_GAS * T_inf)
        u_inf = 5.0 * a_inf 
        v_inf = 0.0
        Yf_inf = 1.0
        
        E = p_inf/((self.cfg.GAMMA-1)*rho_inf) + 0.5*(u_inf**2 + v_inf**2)
        
        self.U[0,:,:] = rho_inf
        self.U[1,:,:] = rho_inf * u_inf
        self.U[2,:,:] = rho_inf * v_inf
        self.U[3,:,:] = rho_inf * E
        self.U[4,:,:] = rho_inf * Yf_inf

    def get_primitives(self, U_in):
        rho = U_in[0]
        rho = np.maximum(rho, 1e-6) 
        u = U_in[1] / rho
        v = U_in[2] / rho
        E_total = U_in[3] / rho
        E_kin = 0.5 * (u**2 + v**2)
        p = rho * (self.cfg.GAMMA - 1) * (E_total - E_kin)
        p = np.maximum(p, 1e-6) # Safety Floor
        Yf = U_in[4] / rho
        return rho, u, v, p, Yf

    def compute_residuals(self, U_curr):
        rho, u, v, p, Yf = self.get_primitives(U_curr)
        a = np.sqrt(self.cfg.GAMMA * p / rho)
        
        # --- 1. PHYSICAL FLUXES ---
        F = np.zeros_like(U_curr)
        F[0] = U_curr[1]
        F[1] = U_curr[1]*u + p
        F[2] = U_curr[1]*v
        F[3] = (U_curr[3] + p)*u
        F[4] = U_curr[4]*u
        
        G = np.zeros_like(U_curr)
        G[0] = U_curr[2]
        G[1] = U_curr[2]*u
        G[2] = U_curr[2]*v + p
        G[3] = (U_curr[3] + p)*v
        G[4] = U_curr[4]*v
        
        # --- 2. TRANSFORMED FLUXES ---
        F_prime = self.xi_x * F + self.xi_y * G
        G_prime = self.eta_x * F + self.eta_y * G
        
        # --- 3. NUMERICAL FLUX (Rusanov) ---
        U_contra = self.xi_x * u + self.xi_y * v
        V_contra = self.eta_x * u + self.eta_y * v
        
        grad_xi = np.sqrt(self.xi_x**2 + self.xi_y**2)
        grad_eta = np.sqrt(self.eta_x**2 + self.eta_y**2)
        
        lam_xi = np.abs(U_contra) + a * grad_xi
        lam_eta = np.abs(V_contra) + a * grad_eta
        
        Residual = np.zeros_like(U_curr)
        
        # XI Flux
        eig_xi = np.maximum(lam_xi[:-1, :], lam_xi[1:, :])
        Flux_Xi = 0.5 * (F_prime[:, :-1, :] + F_prime[:, 1:, :]) - 0.5 * eig_xi * (U_curr[:, 1:, :] - U_curr[:, :-1, :])
        Residual[:, 1:-1, :] -= (Flux_Xi[:, 1:, :] - Flux_Xi[:, :-1, :])
        
        # ETA Flux
        eig_eta = np.maximum(lam_eta[:, :-1], lam_eta[:, 1:])
        Flux_Eta = 0.5 * (G_prime[:, :, :-1] + G_prime[:, :, 1:]) - 0.5 * eig_eta * (U_curr[:, :, 1:] - U_curr[:, :, :-1])
        Residual[:, :, 1:-1] -= (Flux_Eta[:, :, 1:] - Flux_Eta[:, :, :-1])
        
        Residual *= self.J
        
        # --- 4. SOURCE TERM (Chemistry) ---
        T = p / (rho * self.cfg.R_GAS)
        mask = np.where(T > 800.0, 1.0, 0.0) 
        w_dot = -self.cfg.A_PRE * rho * Yf * np.exp(-self.cfg.TA / T) * mask
        
        Residual[3, :, :] += -w_dot * self.cfg.Q_HEAT
        Residual[4, :, :] += w_dot

        # --- 5. BCs ---
        Residual[:, 0, :] = 0 
        Residual[:, -1, :] = 0 
        Residual[:, :, 0] = 0
        Residual[:, :, -1] = 0
        
        return Residual, np.max(lam_xi + lam_eta)

    def solve(self):
        print(f"Starting Manual TVD-RK3 Simulation (Analytic Metrics)...")
        self.initialize()
        
        t = 0.0
        iter_count = 0
        
        while t < self.cfg.T_FINAL:
            _, max_lam = self.compute_residuals(self.U)
            
            # Robust dt calculation
            if max_lam < 1e-6: max_lam = 1.0
            dt = self.cfg.CFL / max_lam
            
            # Safety clamp for dt
            dt = min(dt, 1e-7) 
            
            if t + dt > self.cfg.T_FINAL: dt = self.cfg.T_FINAL - t
            
            # RK3 Steps
            U0 = self.U.copy()
            R0, _ = self.compute_residuals(U0)
            U1 = U0 + dt * R0
            
            R1, _ = self.compute_residuals(U1)
            U2 = 0.75*U0 + 0.25*U1 + 0.25*dt*R1
            
            R2, _ = self.compute_residuals(U2)
            self.U = (1/3)*U0 + (2/3)*U2 + (2/3)*dt*R2
            
            t += dt
            iter_count += 1
            
            if iter_count % 50 == 0:
                rho, _, _, p, _ = self.get_primitives(self.U)
                T = p / (rho * self.cfg.R_GAS)
                print(f"Iter {iter_count:4d} | t = {t:.4e}s | dt = {dt:.4e} | Max T = {np.max(T):.1f} K")

        print("Simulation Complete.")
        self.plot_result(t)

    def plot_result(self, t):
        rho, u, v, p, Yf = self.get_primitives(self.U)
        T = p / (rho * self.cfg.R_GAS)
        plt.figure(figsize=(10, 6))
        cp = plt.contourf(self.X, self.Y, T, levels=50, cmap='jet')
        plt.colorbar(cp, label='Temperature [K]')
        plt.plot(self.X[:,0], self.Y[:,0], 'k-', linewidth=2)
        plt.title(f'Oblique Detonation Wave (Analytic Metrics, t={t:.2e}s)')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('odwe_final.png', dpi=150)
        plt.show()

# ==========================================
# 4. RUNNER
# ==========================================
if __name__ == "__main__":
    cfg = Config()
    solver = ODWESolver(cfg)
    solver.solve()