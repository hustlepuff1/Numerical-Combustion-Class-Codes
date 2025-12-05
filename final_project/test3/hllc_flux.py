import numpy as np

class RocketConfig:
    """
    Configuration based on NASA CEA Output 'final_porj_nasa_cea.txt'
    """
    # Average Molecular Weight from Chamber to Exit (approx)
    # MW ranges from 24.63 to 25.58 [cite: 10]
    MW = 24.63  
    
    # Specific Gas Constant R = R_univ / MW
    R_GAS = 8314.46 / MW  # ~337.5 J/(kg K)
    
    # Heat Capacity Ratio (Gamma)
    # CEA shows 1.125 (chamber) to 1.155 (exit) [cite: 11]
    # We use a constant average for the "Simple" solver
    GAMMA = 1.13
    
    # Chamber Conditions (for Boundary Conditions)
    P_CHAMBER = 64.0 * 1e5  # 64 bar to Pa [cite: 1, 8]
    T_CHAMBER = 3729.08     # K [cite: 9]
    RHO_CHAMBER = 5.4027    # kg/m^3 [cite: 9]

def get_sound_speed(rho, p):
    """Computes speed of sound a = sqrt(gamma * p / rho)"""
    return np.sqrt(RocketConfig.GAMMA * p / rho)

def conservative_to_primitive(U):
    """
    Converts Conservative [rho, rho*u, rho*v, E] to Primitive [rho, u, v, p]
    Includes FLOORS to prevent negative pressure/density crashes.
    """
    MIN_VAL = 1e-6  # Safety floor
    
    rho = U[0]
    
    # 1. Protect Density
    if rho < MIN_VAL:
        rho = MIN_VAL
        
    u = U[1] / rho
    v = U[2] / rho
    E = U[3]
    
    # Internal Energy calc
    kinetic = 0.5 * (u**2 + v**2)
    p = (RocketConfig.GAMMA - 1.0) * (E - rho * kinetic)
    
    # 2. Protect Pressure
    if p < MIN_VAL:
        p = MIN_VAL
    
    return rho, u, v, p

def hllc_flux(UL, UR, normal):
    """
    Computes HLLC Flux for 2D Euler Equations.
    
    Parameters:
    UL, UR : Left and Right Conservative State Vectors (4,)
    normal : Normal vector [nx, ny] of the face
    """
    # 1. Rotate velocities to face normal frame
    nx, ny = normal
    
    # Get primitives
    rho_L, u_L_phys, v_L_phys, p_L = conservative_to_primitive(UL)
    rho_R, u_R_phys, v_R_phys, p_R = conservative_to_primitive(UR)
    
    # Normal and Tangential velocities
    u_L = u_L_phys * nx + v_L_phys * ny
    v_L = -u_L_phys * ny + v_L_phys * nx
    
    u_R = u_R_phys * nx + v_R_phys * ny
    v_R = -u_R_phys * ny + v_R_phys * nx
    
    # Sound speeds
    a_L = get_sound_speed(rho_L, p_L)
    a_R = get_sound_speed(rho_R, p_R)
    
    # 2. Wave Speed Estimates (Davis approximation for simplicity)
    # You can also use Roe-averaged speeds for higher accuracy
    S_L = min(u_L - a_L, u_R - a_R)
    S_R = max(u_L + a_L, u_R + a_R)
    
    # Contact wave speed S_star
    num = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)
    den = rho_L * (S_L - u_L) - rho_R * (S_R - u_R)
    S_star = num / den
    
    # 3. HLLC Flux Assembly
    # Flux vector F(U) in normal direction
    def get_flux(rho, u, v, p, E):
        return np.array([
            rho * u,
            rho * u**2 + p,
            rho * u * v,
            u * (E + p)
        ])
    
    F_L = get_flux(rho_L, u_L, v_L, p_L, UL[3])
    F_R = get_flux(rho_R, u_R, v_R, p_R, UR[3])
    
    # Logic to select the correct flux
    if S_L >= 0:
        F_n = F_L
    elif S_L < 0 and S_star >= 0:
        # Star Region Left
        # We need to construct the Star State U*_L
        factor_L = rho_L * (S_L - u_L) / (S_L - S_star)
        U_star_L = factor_L * np.array([
            1.0,
            S_star,
            v_L, # Tangential velocity preserved
            UL[3]/rho_L + (S_star - u_L) * (S_star + p_L/(rho_L*(S_L - u_L)))
        ])
        F_n = F_L + S_L * (U_star_L - np.array([rho_L, rho_L*u_L, rho_L*v_L, UL[3]]))
        
    elif S_star < 0 and S_R >= 0:
        # Star Region Right
        factor_R = rho_R * (S_R - u_R) / (S_R - S_star)
        U_star_R = factor_R * np.array([
            1.0,
            S_star,
            v_R, # Tangential velocity preserved
            UR[3]/rho_R + (S_star - u_R) * (S_star + p_R/(rho_R*(S_R - u_R)))
        ])
        F_n = F_R + S_R * (U_star_R - np.array([rho_R, rho_R*u_R, rho_R*v_R, UR[3]]))
        
    else: # S_R < 0
        F_n = F_R

    # 4. Rotate Flux back to physical (x,y) frame
    # F_phys = [mass, mom_x, mom_y, energy]
    F_phys = np.zeros(4)
    F_phys[0] = F_n[0]
    F_phys[1] = F_n[1] * nx - F_n[2] * ny
    F_phys[2] = F_n[1] * ny + F_n[2] * nx
    F_phys[3] = F_n[3]
    
    return F_phys

def reaction_step(U, dt):
    """
    Stabilized Finite-Rate Chemistry Step
    """
    # 1. Constants 
    # Reduced Pre-exponential factor to slow down reaction slightly
    A_pre = 1.0e8        
    Ta = 15000.0         
    Q_heat = 2.0e6       # Keep heat release moderate
    
    rho = U[:,:,0]
    rhoY = U[:,:,4]
    E = U[:,:,3]
    
    # Safety: Avoid division by zero
    rho_safe = np.maximum(rho, 1e-6)
    
    # Recover Temperature
    u = U[:,:,1] / rho_safe
    v = U[:,:,2] / rho_safe
    ke = 0.5 * (u**2 + v**2)
    p = (RocketConfig.GAMMA - 1.0) * (E - rho * ke)
    # Clamp pressure to be positive
    p = np.maximum(p, 1e-6)
    T = p / (rho_safe * RocketConfig.R_GAS)
    
    # --- STABILITY MASKS ---
    # Mask 1: Vacuum Cutoff (Don't burn if density is too low)
    valid_cells = (rho > 0.05) 
    
    # Mask 2: Temperature Cutoff (Don't burn if already too hot)
    # This prevents the 33,000 K explosions
    cool_enough = (T < 4000.0)
    
    # Combine masks
    active_mask = np.logical_and(valid_cells, cool_enough)
    
    # 2. Calculate Rates only where safe
    # We use a temporary rate array initialized to zero
    w_dot = np.zeros_like(rho)
    
    # Extract values for active cells only
    T_active = T[active_mask]
    rho_active = rho[active_mask]
    rhoY_active = rhoY[active_mask]
    
    # Arrhenius Rate
    # Limit T to avoid overflow in exp
    T_active = np.maximum(T_active, 300.0) 
    k_rate = A_pre * np.exp(-Ta / T_active)
    
    # Analytical Integration: Y_new = Y_old * exp(-k * dt)
    Y_old = rhoY_active / rho_active
    Y_new = Y_old * np.exp(-k_rate * dt)
    
    # Calculate consumption
    dY = Y_old - Y_new
    w_dot[active_mask] = dY * rho_active # kg/m^3 consumed per step
    
    # 3. Update Conservative Variables
    # Remove reactant mass
    U[:,:,4] -= w_dot 
    # Add Heat Release to Energy
    U[:,:,3] += w_dot * Q_heat
    
    return U

def hllc_flux_5eq(UL, UR, normal):
    """
    Upgraded HLLC Solver for 5 Equations (Fluid + 1 Species)
    """
    # ... (Keep the exact same logic as before for rho, u, v, p, S_L, S_R, S_star) ...
    # COPY the rotation and wave speed code from your previous hllc_flux here
    # UNTIL you calculate 'S_star'. 
    
    # --- INTERNAL RECALCULATION FOR CONTEXT (You can paste this over old logic) ---
    nx, ny = normal
    rho_L, u_L_phys, v_L_phys, p_L = conservative_to_primitive(UL)
    rho_R, u_R_phys, v_R_phys, p_R = conservative_to_primitive(UR)
    
    u_L = u_L_phys * nx + v_L_phys * ny
    v_L = -u_L_phys * ny + v_L_phys * nx
    u_R = u_R_phys * nx + v_R_phys * ny
    v_R = -u_R_phys * ny + v_R_phys * nx
    
    a_L = get_sound_speed(rho_L, p_L)
    a_R = get_sound_speed(rho_R, p_R)
    
    S_L = min(u_L - a_L, u_R - a_R)
    S_R = max(u_L + a_L, u_R + a_R)
    
    num = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)
    den = rho_L * (S_L - u_L) - rho_R * (S_R - u_R)
    S_star = num / den
    # --------------------------------------------------------------------------

    # Assemble Fluxes for 5 Equations
    # UL[4] is rho*Y
    Y_L = UL[4] / rho_L
    Y_R = UR[4] / rho_R

    # F_L (Physical Flux)
    F_L = np.array([
        rho_L * u_L,
        rho_L * u_L**2 + p_L,
        rho_L * u_L * v_L,
        u_L * (UL[3] + p_L),
        rho_L * u_L * Y_L  # Species Flux
    ])
    
    F_R = np.array([
        rho_R * u_R,
        rho_R * u_R**2 + p_R,
        rho_R * u_R * v_R,
        u_R * (UR[3] + p_R),
        rho_R * u_R * Y_R
    ])

    # HLLC Logic with Species Transport
    if S_L >= 0:
        F_n = F_L
    elif S_L < 0 and S_star >= 0:
        # Left Star State
        factor = rho_L * (S_L - u_L) / (S_L - S_star)
        # Species Mass Fraction is preserved across Contact Wave -> Y*_L = Y_L
        U_star_L = factor * np.array([
            1.0, S_star, v_L,
            UL[3]/rho_L + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))),
            Y_L  # Species density scales with fluid density
        ])
        F_n = F_L + S_L * (U_star_L - UL)
        
    elif S_star < 0 and S_R >= 0:
        # Right Star State
        factor = rho_R * (S_R - u_R) / (S_R - S_star)
        U_star_R = factor * np.array([
            1.0, S_star, v_R,
            UR[3]/rho_R + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))),
            Y_R 
        ])
        F_n = F_R + S_R * (U_star_R - UR)
    else:
        F_n = F_R

    # Rotate back (Handling 5 components)
    F_phys = np.zeros(5)
    F_phys[0] = F_n[0]
    F_phys[1] = F_n[1] * nx - F_n[2] * ny
    F_phys[2] = F_n[1] * ny + F_n[2] * nx
    F_phys[3] = F_n[3]
    F_phys[4] = F_n[4] # Scalar flux doesn't rotate
    
    return F_phys