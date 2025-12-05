import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. CONFIGURATION & PHYSICS CONSTANTS
# ==========================================
class RocketConfig:
    # Based on NASA CEA Output for AP/Al/Butadiene
    MW = 24.63  
    R_GAS = 8314.46 / MW  # ~337.5 J/(kg K)
    GAMMA = 1.13
    
    # Target Chamber Conditions
    P_CHAMBER_TARGET = 64.0 * 1e5  # 64 bar
    T_CHAMBER = 3729.08            # K
    
    # Simulation Settings
    NX = 100
    NY = 20
    L_X = 1.0
    L_Y = 0.2
    CFL = 0.2  # Conservative CFL
    T_FINAL = 0.002

# ==========================================
# 2. NUMERICAL KERNELS (HLLC & CHEMISTRY)
# ==========================================
def conservative_to_primitive(U):
    """Robust conversion with Safety Floors."""
    MIN_RHO = 1e-4
    MIN_P = 1e-4
    
    rho = U[0]
    if rho < MIN_RHO: rho = MIN_RHO
    
    u = U[1] / rho
    v = U[2] / rho
    E = U[3]
    
    kinetic = 0.5 * (u**2 + v**2)
    p = (RocketConfig.GAMMA - 1.0) * (E - rho * kinetic)
    
    if p < MIN_P: p = MIN_P
        
    return rho, u, v, p

def get_sound_speed(rho, p):
    return np.sqrt(RocketConfig.GAMMA * p / rho)

def hllc_flux_5eq(UL, UR, normal):
    """5-Equation HLLC Riemann Solver (2D Rotated)."""
    nx, ny = normal
    
    # 1. Reconstruct Primitives
    rho_L, u_L_phys, v_L_phys, p_L = conservative_to_primitive(UL)
    rho_R, u_R_phys, v_R_phys, p_R = conservative_to_primitive(UR)
    
    # 2. Rotate to Face Normal
    u_L = u_L_phys * nx + v_L_phys * ny
    v_L = -u_L_phys * ny + v_L_phys * nx
    u_R = u_R_phys * nx + v_R_phys * ny
    v_R = -u_R_phys * ny + v_R_phys * nx
    
    # 3. Wave Speeds
    a_L = get_sound_speed(rho_L, p_L)
    a_R = get_sound_speed(rho_R, p_R)
    
    # Roe-averaged or Davis estimates
    S_L = min(u_L - a_L, u_R - a_R)
    S_R = max(u_L + a_L, u_R + a_R)
    
    # S_star (Contact Wave)
    denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R)
    if abs(denom) < 1e-9: denom = 1e-9 # Prevent div/0
        
    num = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)
    S_star = num / denom
    
    # 4. Flux Assembly
    Y_L = UL[4] / rho_L
    Y_R = UR[4] / rho_R
    
    def get_flux_vector(rho, u, v, p, E, Y_mass):
        return np.array([
            rho * u,
            rho * u**2 + p,
            rho * u * v,
            u * (E + p),
            u * Y_mass
        ])
    
    # Physical fluxes (rotated inputs)
    F_L_n = get_flux_vector(rho_L, u_L, v_L, p_L, UL[3], UL[4])
    F_R_n = get_flux_vector(rho_R, u_R, v_R, p_R, UR[3], UR[4])
    
    # 5. HLLC Logic
    if S_L >= 0:
        F_n = F_L_n
    elif S_L < 0 and S_star >= 0:
        # Star Left
        factor = rho_L * (S_L - u_L) / (S_L - S_star)
        U_star_L = factor * np.array([
            1.0, S_star, v_L,
            UL[3]/rho_L + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))),
            Y_L
        ])
        F_n = F_L_n + S_L * (U_star_L - UL)
    elif S_star < 0 and S_R >= 0:
        # Star Right
        factor = rho_R * (S_R - u_R) / (S_R - S_star)
        U_star_R = factor * np.array([
            1.0, S_star, v_R,
            UR[3]/rho_R + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))),
            Y_R
        ])
        F_n = F_R_n + S_R * (U_star_R - UR)
    else:
        F_n = F_R_n
        
    # 6. Rotate Back
    F_phys = np.zeros(5)
    F_phys[0] = F_n[0]
    F_phys[1] = F_n[1] * nx - F_n[2] * ny
    F_phys[2] = F_n[1] * ny + F_n[2] * nx
    F_phys[3] = F_n[3]
    F_phys[4] = F_n[4]
    
    return F_phys

def reaction_step(U, dt):
    """Stabilized Finite Rate Chemistry"""
    # Parameters
    A_pre = 2.0e7      # Tuned for project scale
    Ta = 15000.0
    Q_heat = 1.5e6     # Moderate heat release
    
    rho = U[:,:,0]
    rhoY = U[:,:,4]
    E = U[:,:,3]
    
    # Recover T safely
    rho_safe = np.maximum(rho, 1e-4)
    u = U[:,:,1] / rho_safe
    v = U[:,:,2] / rho_safe
    p = (RocketConfig.GAMMA - 1.0) * (E - 0.5 * rho * (u**2 + v**2))
    p = np.maximum(p, 1e-4)
    T = p / (rho_safe * RocketConfig.R_GAS)
    
    # Reaction Masks (Where to burn?)
    # 1. Must have fuel (Y > 0)
    # 2. Must be hot enough (T > 600)
    # 3. Must not be too hot already (T < 4000) -> Prevents Runaway
    mask = (rhoY > 1e-6) & (T > 600.0) & (T < 4000.0) & (rho > 0.01)
    
    if not np.any(mask):
        return U
        
    # Calculate Kinetics
    T_loc = T[mask]
    rho_loc = rho[mask]
    Y_loc = rhoY[mask] / rho_loc
    
    # Arrhenius
    k = A_pre * np.exp(-Ta / T_loc)
    
    # Analytical Update
    Y_new = Y_loc * np.exp(-k * dt)
    dY_mass = (Y_loc - Y_new) * rho_loc
    
    # Apply changes
    U[mask, 4] -= dY_mass           # Consume Fuel
    U[mask, 3] += dY_mass * Q_heat  # Add Heat
    
    return U

# ==========================================
# 3. MAIN SIMULATION LOOP
# ==========================================
def main():
    # Setup Grid
    dx = RocketConfig.L_X / RocketConfig.NX
    dy = RocketConfig.L_Y / RocketConfig.NY
    
    # Initialize Domain (Ambient Air - 10 Bar cushion)
    U = np.zeros((RocketConfig.NX, RocketConfig.NY, 5))
    P_amb = 10.0 * 1e5 
    T_amb = 300.0
    Rho_amb = P_amb / (RocketConfig.R_GAS * T_amb)
    E_amb = P_amb / (RocketConfig.GAMMA - 1.0)
    
    for i in range(RocketConfig.NX):
        for j in range(RocketConfig.NY):
            U[i,j,:] = [Rho_amb, 0.0, 0.0, E_amb, 0.0]

    t = 0.0
    it = 0
    print(f"Starting Simulation with Soft Start...")
    
    # --- TIME LOOP ---
    while t < RocketConfig.T_FINAL:
        
        # A. Stability: Vacuum Cleaner
        # Reset cells with unphysical density
        bad_cells = U[:,:,0] < 1e-3
        if np.any(bad_cells):
            U[bad_cells, :] = [Rho_amb, 0.0, 0.0, E_amb, 0.0]
            
        # B. Calculate Time Step
        rho = np.maximum(U[:,:,0], 1e-4)
        u = U[:,:,1] / rho
        v = U[:,:,2] / rho
        p = (RocketConfig.GAMMA - 1.0) * (U[:,:,3] - 0.5 * rho * (u**2 + v**2))
        p = np.maximum(p, 1e-4)
        a = np.sqrt(RocketConfig.GAMMA * p / rho)
        
        max_wave = np.max(np.sqrt(u**2 + v**2) + a)
        max_wave = max(max_wave, 500.0) # Floor for initial step
        
        dt = RocketConfig.CFL * min(dx, dy) / max_wave
        
        # C. Flux Arrays
        Flux_X = np.zeros((RocketConfig.NX+1, RocketConfig.NY, 5))
        Flux_Y = np.zeros((RocketConfig.NX, RocketConfig.NY+1, 5))
        
        # D. Internal Fluxes
        for i in range(1, RocketConfig.NX):
            for j in range(RocketConfig.NY):
                Flux_X[i,j,:] = hllc_flux_5eq(U[i-1,j,:], U[i,j,:], [1,0])
                
        for i in range(RocketConfig.NX):
            for j in range(1, RocketConfig.NY):
                Flux_Y[i,j,:] = hllc_flux_5eq(U[i,j-1,:], U[i,j,:], [0,1])
                
        # E. Boundary Conditions (Soft Start)
        # -----------------------------------
        # Ramp Pressure from 10 bar to 64 bar over 1000 steps
        ramp_factor = min(it / 1000.0, 1.0)
        P_current = P_amb + ramp_factor * (RocketConfig.P_CHAMBER_TARGET - P_amb)
        
        # Injector Properties (Hot for ignition)
        T_inj = 2500.0
        Rho_inj = P_current / (RocketConfig.R_GAS * T_inj)
        E_inj = P_current / (RocketConfig.GAMMA - 1.0)
        
        # Inlet BC
        j_start = int(RocketConfig.NY * 0.4)
        j_end   = int(RocketConfig.NY * 0.6)
        
        for j in range(RocketConfig.NY):
            if j_start <= j <= j_end:
                # Jet Inlet
                U_inlet = np.array([Rho_inj, 0.0, 0.0, E_inj, Rho_inj * 1.0])
                Flux_X[0,j,:] = hllc_flux_5eq(U_inlet, U[0,j,:], [1,0])
            else:
                # Wall
                U_ghost = U[0,j,:].copy(); U_ghost[1] *= -1
                Flux_X[0,j,:] = hllc_flux_5eq(U_ghost, U[0,j,:], [1,0])
                
        # Outlet BC
        for j in range(RocketConfig.NY):
            Flux_X[RocketConfig.NX,j,:] = hllc_flux_5eq(U[-1,j,:], U[-1,j,:], [1,0])
            
        # Top/Bottom Walls
        for i in range(RocketConfig.NX):
            U_bot = U[i,0,:].copy(); U_bot[2] *= -1
            Flux_Y[i,0,:] = hllc_flux_5eq(U_bot, U[i,0,:], [0,1])
            
            U_top = U[i,-1,:].copy(); U_top[2] *= -1
            Flux_Y[i,RocketConfig.NY,:] = hllc_flux_5eq(U[i,-1,:], U_top, [0,1])
            
        # F. Update
        U -= (dt/dx) * (Flux_X[1:,:,:] - Flux_X[:-1,:,:])
        U -= (dt/dy) * (Flux_Y[:,1:,:] - Flux_Y[:,:-1,:])
        
        # G. Reaction
        U = reaction_step(U, dt)
        
        t += dt
        it += 1
        
        if it % 50 == 0:
            rho_d = np.maximum(U[:,:,0], 1e-4)
            T_d = (RocketConfig.GAMMA - 1.0)*(U[:,:,3]-0.5*(U[:,:,1]**2+U[:,:,2]**2)/rho_d)/rho_d/RocketConfig.R_GAS
            u_mag = np.sqrt((U[:,:,1]/rho_d)**2 + (U[:,:,2]/rho_d)**2)
            a_d = np.sqrt(RocketConfig.GAMMA * RocketConfig.R_GAS * T_d)
            mach = u_mag / a_d
            
            print(f"Iter {it:4d} | P_inlet: {P_current/1e5:4.1f} bar | Max Mach: {np.max(mach):5.2f} | Max T: {np.max(T_d):6.1f} K")

    # --- PLOTTING ---
    print("Simulation Complete. Generating Plot...")
    rho_f = np.maximum(U[:,:,0], 1e-4)
    u_f = U[:,:,1]/rho_f
    v_f = U[:,:,2]/rho_f
    speed = np.sqrt(u_f**2 + v_f**2)
    p_f = (RocketConfig.GAMMA - 1.0)*(U[:,:,3] - 0.5*rho_f*speed**2)
    T_f = p_f / (rho_f * RocketConfig.R_GAS)
    
    plt.figure(figsize=(10,6))
    plt.subplot(2,1,1)
    plt.imshow(speed.T, origin='lower', extent=[0,1,0,0.2], cmap='jet', aspect='auto')
    plt.colorbar(label='Velocity (m/s)')
    plt.title('Final Velocity Field')
    
    plt.subplot(2,1,2)
    plt.imshow(T_f.T, origin='lower', extent=[0,1,0,0.2], cmap='inferno', aspect='auto')
    plt.colorbar(label='Temperature (K)')
    plt.title('Final Temperature Field')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()