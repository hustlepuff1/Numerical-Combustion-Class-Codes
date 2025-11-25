import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ==========================================
# 1. PHYSICAL CONSTANTS
# ==========================================
R_univ = 8.314        # J/(mol K)
rho0   = 200.0        # kg/m^3

# Pre-exponential Factors (1/s)
A1 = 4.4e10
A2 = 4.0e10
A3 = 3.0e10

# Activation Energies (J/mol)
E1 = 190000.0
E2 = 180000.0
E3 = 140000.0

# Heat of Reaction (J/kg)
dH1 = 2.7e5
dH2 = -8.0e5
dH3 = -4.0e6

# Specific Heats (J/kg K)
Cv_A = 1200.0
Cv_B = 1200.0
Cv_C = 2000.0
Cv_D = 2000.0

# ==========================================
# 2. THE ODE SYSTEM
# ==========================================
def combustion_ode(t, U):
    Y_A, Y_B, Y_C, Y_D, T = U
    
    # Safety Check
    if T < 200: T = 200
    
    # Reaction Rates (Arrhenius)
    k1 = A1 * np.exp(-E1 / (R_univ * T))
    r1 = k1 * rho0 * Y_A
    
    k2 = A2 * np.exp(-E2 / (R_univ * T))
    r2 = k2 * rho0 * Y_B
    
    k3 = A3 * np.exp(-E3 / (R_univ * T))
    r3 = k3 * rho0 * Y_C
    
    # Mixture Cv
    Cv_mix = Y_A*Cv_A + Y_B*Cv_B + Y_C*Cv_C + Y_D*Cv_D
    
    # Source Terms
    dYA_dt = -r1 / rho0
    dYB_dt = (r1 - r2) / rho0
    dYC_dt = (r2 - r3) / rho0
    dYD_dt = r3 / rho0
    
    heat_release = r1*dH1 + r2*dH2 + r3*dH3
    dT_dt = -heat_release / (rho0 * Cv_mix)
    
    return [dYA_dt, dYB_dt, dYC_dt, dYD_dt, dT_dt]

# ==========================================
# 3. SOLVER & PLOTTING
# ==========================================
def run_and_plot():
    print("Running SciPy Stiff Solver (BDF)...")
    
    # Initial Conditions & Time
    y0 = [1.0, 0.0, 0.0, 0.0, 1500.0]
    t_span = (0.0, 0.3e-4) # 30 microseconds
    
    # Solve using BDF (Implicit)
    sol = solve_ivp(combustion_ode, t_span, y0, method='BDF', rtol=1e-6, atol=1e-10)
    
    if not sol.success:
        print("Solver Failed!")
        return

    # --- FIGURE 1: LINEAR SCALE (Temp + Species) ---
    plt.figure(figsize=(10, 8))
    
    # Subplot 1: Temperature
    plt.subplot(2, 1, 1)
    plt.plot(sol.t, sol.y[4], 'r-', linewidth=2, label='Temperature')
    plt.ylabel('Temperature (K)', fontweight='bold')
    plt.title('Figure 1: Linear Time Scale Results', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.legend()
    
    # Subplot 2: Species
    plt.subplot(2, 1, 2)
    plt.plot(sol.t, sol.y[0], label='Y_A (Fuel)')
    plt.plot(sol.t, sol.y[1], label='Y_B')
    plt.plot(sol.t, sol.y[2], label='Y_C')
    plt.plot(sol.t, sol.y[3], label='Y_D (Prod)')
    plt.xlabel('Time (s)', fontweight='bold')
    plt.ylabel('Mass Fraction', fontweight='bold')
    plt.grid(True, linestyle='--')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('linear_results_open_source.png')
    
    # --- FIGURE 2: LOG SCALE (Species ONLY) ---
    plt.figure(figsize=(10, 6))
    
    plt.plot(sol.t, sol.y[0], label='Y_A (Fuel)')
    plt.plot(sol.t, sol.y[1], label='Y_B')
    plt.plot(sol.t, sol.y[2], label='Y_C')
    plt.plot(sol.t, sol.y[3], label='Y_D (Prod)')
    
    # LOG MAGIC HERE
    plt.xscale('log')
    
    plt.title('Figure 2: Species Mass Fraction (Log Time Scale)', fontsize=14)
    plt.xlabel('Time (s) - Log Scale', fontweight='bold')
    plt.ylabel('Mass Fraction', fontweight='bold')
    plt.grid(True, which="both", linestyle='--', alpha=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('log_species_open_source.png')
    
    print("Graphs saved as 'linear_results_open_source.png' and 'log_species_open_source.png'")
    # plt.show() # Uncomment to see pop-ups

if __name__ == "__main__":
    run_and_plot()