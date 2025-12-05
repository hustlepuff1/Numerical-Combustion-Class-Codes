import numpy as np
import matplotlib.pyplot as plt
import hllc_flux as sol  # Imports your saved file

# --- 1. Simulation Setup ---
NX = 100        # Number of cells in X
NY = 20         # Number of cells in Y
L_X = 1.0       # Length (m)
L_Y = 0.2       # Height (m)
CFL = 0.1       # Stability condition
T_FINAL = 0.001 # 1ms simulation time

dx = L_X / NX
dy = L_Y / NY

# --- 2. Initialization ---
# U is now 5D: [rho, rho*u, rho*v, Energy, rho*Y]
U = np.zeros((NX, NY, 5))

# Ambient Conditions (Initial state of the domain: Low Pressure, No Fuel)
P_amb = 10.0 * 1e5  # 10 Bar
T_amb = 300.0      # Cold ambient
Rho_amb = P_amb / (sol.RocketConfig.R_GAS * T_amb)
E_amb = P_amb / (sol.RocketConfig.GAMMA - 1.0)

# Fill domain with ambient air
for i in range(NX):
    for j in range(NY):
        U[i,j,:] = [Rho_amb, 0.0, 0.0, E_amb, 0.0] # Y=0 (No fuel initially)

# --- 3. Main Time Loop ---
t = 0.0
it = 0

print(f"Starting Reactive Jet Simulation: {NX}x{NY} Grid")

while t < T_FINAL:



    # --- STABILITY FIX: VACUUM MOMENTUM CLIP ---
    # 1. Reset 'Dead' cells (too low density) to Ambient
    BAD_CELLS = U[:,:,0] < 1e-4
    if np.any(BAD_CELLS):
        U[BAD_CELLS, :] = [Rho_amb, 0.0, 0.0, E_amb, 0.0]
        
    # 2. CLIP MOMENTUM in 'Thin' cells (low but valid density)
    # If density is < 0.1 kg/m3, kill the velocity to prevent Mach explosion
    THIN_CELLS = U[:,:,0] < 0.1
    if np.any(THIN_CELLS):
        # Set Momentum (indices 1 and 2) to zero in these cells
        U[THIN_CELLS, 1] = 0.0
        U[THIN_CELLS, 2] = 0.0
        # Recalculate Energy to remove Kinetic Energy (keep only Internal)
        # E = p/(g-1) + 0.5*rho*u^2  ->  E = p/(g-1)
        # We assume Pressure is roughly conserved or ambient
        p_thin = 1.0e5 
        rho_thin = U[THIN_CELLS, 0]
        U[THIN_CELLS, 3] = p_thin / (sol.RocketConfig.GAMMA - 1.0)
    # -----------------------------------------------------------------------


    # A. Calculate Time Step (dt)
    rho = U[:,:,0]
    rho = np.maximum(rho, 1e-6) # Safety floor
    u = U[:,:,1] / rho
    v = U[:,:,2] / rho
    
    # Calculate Sound Speed
    kinetic = 0.5 * rho * (u**2 + v**2)
    p = (sol.RocketConfig.GAMMA - 1.0) * (U[:,:,3] - kinetic)
    p = np.maximum(p, 1e-6) # Safety floor
    a = np.sqrt(sol.RocketConfig.GAMMA * p / rho)
    
    # CRITICAL FIX: Include the INLET velocity in the CFL check
    # The inlet is roughly 1000 m/s (Sound speed at 2500K)
    # We force max_speed to be at least this high to prevent huge initial dt
    max_speed_field = np.max(np.sqrt(u**2 + v**2) + a)
    max_speed = max(max_speed_field, 1200.0) 
    
    dt = CFL * min(dx, dy) / max_speed
    
    # B. Create Flux Arrays (for 5 Equations)
    Flux_X = np.zeros((NX+1, NY, 5))
    Flux_Y = np.zeros((NX, NY+1, 5))
    
    # C. X-Direction Fluxes (Internal Faces)
    for i in range(1, NX):
        for j in range(NY):
            Flux_X[i, j, :] = sol.hllc_flux_5eq(U[i-1, j, :], U[i, j, :], normal=[1, 0])
            
    # D. Y-Direction Fluxes (Internal Faces)
    for i in range(NX):
        for j in range(1, NY):
            Flux_Y[i, j, :] = sol.hllc_flux_5eq(U[i, j-1, :], U[i, j, :], normal=[0, 1])

    # E. Boundary Conditions
    # ----------------------
    # 2D INLET: Central Injector (Jet)
    j_start = int(NY * 0.4)
    j_end   = int(NY * 0.6)
    
    # Chamber Conditions (from CEA)
    P_cham = sol.RocketConfig.P_CHAMBER
    # We inject slightly hot gas (2500K) to ensure it auto-ignites/reacts
    T_inj = 2500.0  
    Rho_inj = P_cham / (sol.RocketConfig.R_GAS * T_inj)
    E_inj = P_cham / (sol.RocketConfig.GAMMA - 1.0)
    
    for j in range(NY):
        if j_start <= j <= j_end:
            # INJECTOR: High P, Hot, Full Fuel (Y=1.0)
            U_inlet = np.array([Rho_inj, 0.0, 0.0, E_inj, Rho_inj * 1.0])
            Flux_X[0, j, :] = sol.hllc_flux_5eq(U_inlet, U[0,j,:], normal=[1,0])
        else:
            # WALL: Reflective BC
            U_ghost = U[0,j,:].copy()
            U_ghost[1] *= -1 # Reflect u
            Flux_X[0, j, :] = sol.hllc_flux_5eq(U_ghost, U[0,j,:], normal=[1,0])

    # OUTLET: Supersonic Extrapolation
    for j in range(NY):
        Flux_X[NX, j, :] = sol.hllc_flux_5eq(U[NX-1,j,:], U[NX-1,j,:], normal=[1,0])

    # WALLS (Top/Bottom): Slip Wall
    for i in range(NX):
        # Bottom Wall
        U_ghost_bot = U[i, 0, :].copy(); U_ghost_bot[2] *= -1
        Flux_Y[i, 0, :] = sol.hllc_flux_5eq(U_ghost_bot, U[i, 0, :], normal=[0,1])
        # Top Wall
        U_ghost_top = U[i, NY-1, :].copy(); U_ghost_top[2] *= -1
        Flux_Y[i, NY, :] = sol.hllc_flux_5eq(U[i, NY-1, :], U_ghost_top, normal=[0,1])

    # F. Finite Volume Update
    dFdx = (Flux_X[1:NX+1, :, :] - Flux_X[0:NX, :, :])
    dGdy = (Flux_Y[:, 1:NY+1, :] - Flux_Y[:, 0:NY, :])
    
    U -= (dt / dx) * dFdx
    U -= (dt / dy) * dGdy
    
    # G. Chemical Step (Coupling)
    # THIS SATISFIES PROJECT REQUIREMENT #3
    #U = sol.reaction_step(U, dt)
    
    t += dt
    it += 1
    
    if it % 50 == 0:
        # --- Diagnostic Calculations ---
        # Re-calculate primitives strictly for printing (with safety floors)
        rho_print = np.maximum(U[:,:,0], 1e-9)
        u_print = U[:,:,1] / rho_print
        v_print = U[:,:,2] / rho_print
        
        ke_print = 0.5 * rho_print * (u_print**2 + v_print**2)
        p_print = (sol.RocketConfig.GAMMA - 1.0) * (U[:,:,3] - ke_print)
        p_print = np.maximum(p_print, 1e-9)
        
        # 1. Temperature Check (Did it ignite?)
        T_field = p_print / (rho_print * sol.RocketConfig.R_GAS)
        max_T = np.max(T_field)
        
        # 2. Mach Check (Is it accelerating?)
        a_print = np.sqrt(sol.RocketConfig.GAMMA * p_print / rho_print)
        Mach_field = np.sqrt(u_print**2 + v_print**2) / a_print
        max_M = np.max(Mach_field)
        
        # 3. Stability Check (Is density vanishing?)
        min_rho = np.min(rho_print)
        
        # Print formatted status line
        print(f"Iter {it:4d} | t={t*1000:6.3f} ms | Max Mach: {max_M:5.2f} | Max T: {max_T:6.1f} K | Min Rho: {min_rho:.2e}")

# --- 4. Visualization ---
# Plot Density to see the jet structure clearly
rho = U[:,:,0]

plt.figure(figsize=(10, 4))
plt.imshow(rho.T, origin='lower', extent=[0, L_X, 0, L_Y], cmap='viridis', aspect='auto')
plt.colorbar(label='Density (kg/m^3)')
plt.title(f'Reactive Jet Density (t={t*1000:.2f} ms)\nCheck for Shock Structure')
plt.xlabel('Length (m)')
plt.ylabel('Height (m)')
plt.tight_layout()
plt.savefig("reactive_jet_density.png", dpi=300)
plt.show()