import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Ensure local modules are importable (project root)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from mesh.grid_gen import GridGenerator
from solver.flow import GasFlowSolver
from solver.particles import ParticleSolver
from utils.post_processing import calculate_mach


def main():
    print("--- Phase 3: Fully Coupled Reactive Flow (FVM) ---")
    
    # 1. Grid
    print("Generating Grid...")
    grid = GridGenerator(ni=161, nj=61)
    grid.generate()
    
    # 2. Solvers
    print("Initializing Gas Flow Solver...")
    solver = GasFlowSolver(grid)
    
    print("Initializing Particle Solver...")
    particles = ParticleSolver(grid)
    
    # 3. Conditions
    P_CHAMBER = 6.4e6   # Pa
    T_CHAMBER = 3729.0  # K
    P_EXIT    = 101325  # Pa (1 atm)
    
    solver.initialize_solution(
        p_chamber=P_CHAMBER,
        t_chamber=T_CHAMBER,
        p_exit=P_EXIT
    )
    
    # 4. Time loop
    t = 0.0
    t_end = 0.01  # 10 ms
    iter_count = 0
    
    # Particle injection settings
    inject_interval = 20     # every 20 iterations
    n_inject_each = 50       # particles per injection
    
    print(f"Starting Simulation: t_end={t_end}s")
    
    try:
        while t < t_end:
            dt = solver.get_dt(cfl=0.8)
            
            # --- 1. Particle Phase (Move & Burn) ---
            if iter_count % inject_interval == 0:
                particles.inject_from_inlet(
                    n_new=n_inject_each,
                    u0=50.0, # Give them a kick to match gas inlet
                    v0=0.0
                )
            
            # STEP particles and get the Source Terms (Two-Way Coupling)
            particle_sources = particles.step(dt, solver)
            
            # --- 2. Gas Phase (Flow & React) ---
            # Pass the sources into the gas solver to apply coupling
            solver.step(dt, coupling_terms=particle_sources)
            
            t += dt
            iter_count += 1
            
            if iter_count % 100 == 0:
                mach_exit = calculate_mach(solver)[-1, 0]
                n_p = particles.get_alive_particles()[0].size
                print(
                    f"Iter {iter_count:05d} | t={t:.6f}s | dt={dt:.2e} | "
                    f"Exit Mach={mach_exit:.2f} | Np={n_p}"
                )

    except KeyboardInterrupt:
        print("\nSimulation stopped by user.")
    
    print(f"Simulation completed in {iter_count} iterations, t={t:.6f}s")
    
    # --- PERFORMANCE ANALYSIS ---
    # We need to extract data from the centerline (j=0) and the exit plane (i=-1)
    
    # 1. Mach Number Profile (Centerline)
    mach_field = calculate_mach(solver)
    mach_centerline = mach_field[:, 0]
    
    # Find Throat Location (Minimum physical area)
    # Note: In our grid gen, throat is at index where Y (wall) is min
    throat_idx = np.argmin(solver.grid.Y[:, -1])
    mach_throat = mach_centerline[throat_idx]
    
    # 2. Exit Conditions (Averaged across the exit face)
    # We mass-average the properties at the exit to get "Bulk" values
    rho_exit = solver.rho[-1, :]
    u_exit   = solver.u[-1, :]
    p_exit   = solver.p[-1, :]
    T_exit   = solver.T[-1, :]

    # Compute face areas for the exit plane (axisymmetric rings)
    # area_j = 2*pi * r_j * dy_j  (dy_j approximated from grid spacing)
    r_j = solver.grid.Y[-1, :]
    if solver.grid.nj > 1:
        dy = solver.grid.Y[-1, 1] - solver.grid.Y[-1, 0]
    else:
        dy = 1.0
    area_exit = 2.0 * np.pi * r_j * dy
    
    # Mass Flow Rate (mdot) = Sum(rho * u * Area)
    mdot_local = rho_exit * u_exit * area_exit
    mdot_total = np.sum(mdot_local)
    
    # Bulk Velocity and Pressure
    u_avg = np.sum(u_exit * mdot_local) / (mdot_total + 1e-12)
    p_avg = np.sum(p_exit * area_exit) / (np.sum(area_exit) + 1e-12) # Area-weighted pressure
    T_avg = np.sum(T_exit * mdot_local) / (mdot_total + 1e-12)
    
    # 3. Thrust and Isp Calculation
    # Thrust F = mdot * u_avg + (p_avg - p_ambient) * A_exit
    p_ambient = 101325.0 # Sea level
    A_exit_total = np.sum(area_exit)
    
    thrust = mdot_total * u_avg + (p_avg - p_ambient) * A_exit_total
    g0 = 9.80665
    Isp = thrust / (mdot_total * g0 + 1e-12)
    
    # --- PRINT REPORT ---
    print("\n" + "="*40)
    print("      SOLID ROCKET PERFORMANCE REPORT      ")
    print("="*40)
    print(f"Chamber Pressure : {P_CHAMBER/1e6:.2f} MPa")
    print(f"Throat Mach      : {mach_throat:.3f}  (Target ~1.0)")
    print(f"Exit Mach (Max)  : {np.max(mach_field[-1,:]):.3f}")
    print("-" * 40)
    print(f"Mass Flow Rate   : {mdot_total:.4f} kg/s")
    print(f"Exit Velocity    : {u_avg:.1f} m/s")
    print(f"Exit Temperature : {T_avg:.1f} K")
    print(f"Exit Pressure    : {p_avg/1e5:.2f} bar")
    print("-" * 40)
    print(f"Thrust           : {thrust:.1f} N")
    print(f"Specific Impulse : {Isp:.1f} s")
    print("="*40 + "\n")

    # 5. Visualization
    print("Plotting Results...")
    
    mach_contours = mach_field # Shape is (ni, nj)
    
    x_p, r_p, _, _, _ = particles.get_alive_particles()
    
    # Subsample particles for cleaner plot
    max_plot = 2000
    if x_p.size > max_plot:
        idx = np.random.choice(x_p.size, size=max_plot, replace=False)
        x_p = x_p[idx]
        r_p = r_p[idx]
    
    plt.figure(figsize=(10, 5))
    
    cp = plt.contourf(
        solver.grid.X, solver.grid.Y, mach_contours,
        levels=np.linspace(0, 5.5, 50), cmap='jet'
    )
    cbar = plt.colorbar(cp)
    cbar.set_label('Mach Number')
    
    # Overlay particles
    if x_p.size > 0:
        plt.scatter(x_p, r_p, s=2, c='black', alpha=0.5)
    
    plt.title(f"Phase 3: Reactive Flow (Coupled)\n(Pc={P_CHAMBER/1e6:.1f} MPa)")
    plt.xlabel("Axial Position (m)")
    plt.ylabel("Radial Position (m)")
    plt.axis('equal')
    plt.tight_layout()
    
    plt.savefig('phase3_reactive_result.png')
    print("Result saved to 'phase3_reactive_result.png'")
    # plt.show()

if __name__ == "__main__":
    main()