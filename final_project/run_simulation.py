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
    print("--- Phase 1+2: Cold Gas + Passive Particles ---")
    
    # 1. Grid
    print("Generating Grid...")
    grid = GridGenerator(ni=81, nj=31, length=1.0)
    grid.generate()
    
    # 2. Solvers
    print("Initializing Gas Flow Solver...")
    solver = GasFlowSolver(grid)
    
    print("Initializing Particle Solver...")
    particles = ParticleSolver(grid)
    
    # 3. Conditions
    P_CHAMBER = 6.4e6   # Pa
    T_CHAMBER = 3000.0  # K
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
            dt = solver.get_dt(cfl=0.6)
            solver.step(dt)
            
            # --- Particle phase (gas -> particles only) ---
            if iter_count % inject_interval == 0:
                particles.inject_from_inlet(
                    n_new=n_inject_each,
                    u0=0.0,
                    v0=0.0
                )
            
            particles.advance(dt, solver)
            
            t += dt
            iter_count += 1
            
            if iter_count % 100 == 0:
                mach_exit = calculate_mach(solver)[-1, 0]
                print(
                    f"Iter {iter_count:05d} | t={t:.6f}s | dt={dt:.2e} | "
                    f"Exit Mach={mach_exit:.2f} | Np={particles.get_alive_particles()[0].size}"
                )
    except KeyboardInterrupt:
        print("\nSimulation stopped by user.")
    
    print(f"Simulation completed in {iter_count} iterations, t={t:.6f}s")
    
    # Throat / inlet / exit Mach check
    mach = calculate_mach(solver)
    x_centerline = solver.grid.X[:, 0]
    x_throat_target = 0.5 * solver.grid.length
    i_throat = np.argmin(np.abs(x_centerline - x_throat_target))
    print(f"Mach at inlet centerline  : {mach[0, 0]:.3f}")
    print(f"Mach at throat centerline : {mach[i_throat, 0]:.3f}")
    print(f"Mach at exit centerline   : {mach[-1, 0]:.3f}")
    x_p, r_p, u_p, v_p, _ = particles.get_alive_particles()
    if x_p.size > 0:
        vmag = np.sqrt(u_p**2 + v_p**2)
        print(f"Mean particle speed = {vmag.mean():.1f} m/s, "
            f"max = {vmag.max():.1f} m/s")

    
    # 5. Visualization
    print("Simulation Complete. Plotting Results...")
    
    mach_contours = mach
    x_p, r_p, _, _, _ = particles.get_alive_particles()
    # Subsample for plotting (e.g., at most 1000 points)
    max_plot = 1000
    if x_p.size > max_plot:
        idx = np.random.choice(x_p.size, size=max_plot, replace=False)
        x_p = x_p[idx]
        r_p = r_p[idx]
    
    plt.figure(figsize=(10, 5))
    cp = plt.contourf(
        solver.grid.X, solver.grid.Y, mach_contours,
        levels=20, cmap='jet'
    )
    cbar = plt.colorbar(cp)
    cbar.set_label('Mach Number')
    
    # Overlay particles
    if x_p.size > 0:
        plt.scatter(x_p, r_p, s=5, alpha=0.6)
    
    plt.title(f"Phase 2: Gas Mach & Particle Trajectories\n(Pc={P_CHAMBER/1e6:.1f} MPa)")
    plt.xlabel("Axial Position (m)")
    plt.ylabel("Radial Position (m)")
    plt.axis('equal')
    plt.tight_layout()
    
    plt.savefig('phase2_mach_particles.png')
    print("Result saved to 'phase2_mach_particles.png'")
    plt.show()


if __name__ == "__main__":
    main()
