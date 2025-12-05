# main.py

import numpy as np
import matplotlib.pyplot as plt


from grid import make_1d_grid
from state import prim_to_cons, cons_to_prim
from solver import tvd_rk3_step, compute_dt_cfl
from post import save_snapshot_1d

from grid import make_2d_rect_grid
from state import prim_to_cons, cons_to_prim
from solver import tvd_rk3_step_2d, compute_dt_cfl_2d
from flux import residual_2d
from post import save_snapshot_1d

from grid import make_wedge_grid
from state import prim_to_cons, cons_to_prim
from flux import residual_2d_odwe, max_wave_speed_2d
from solver import tvd_rk3_step_2d



gamma = 1.4

def run_sod():
    # grid
    ni = 201
    x, dx = make_1d_grid(ni, 0.0, 1.0)

    # Sod initial condition
    rhoL, pL, uL = 1.0, 1.0, 0.0
    rhoR, pR, uR = 0.125, 0.1, 0.0

    U = np.zeros((4, ni))
    for i in range(ni):
        if x[i] < 0.5:
            rho = rhoL; p = pL; u = uL
        else:
            rho = rhoR; p = pR; u = uR
        U[:, i] = prim_to_cons(rho, u, p, gamma)

    t      = 0.0
    t_end  = 0.2
    CFL    = 0.5
    it     = 0
    output_every = 50

    while t < t_end:
        dt = compute_dt_cfl(U, dx, gamma, CFL)
        if t + dt > t_end:
            dt = t_end - t

        U = tvd_rk3_step(U, dt, dx, gamma)

        t  += dt
        it += 1

        if it % output_every == 0 or abs(t - t_end) < 1e-10:
            rho, u, p = cons_to_prim(U, gamma)
            print(f"it={it}, t={t:.4f}, dt={dt:.4e}")
            save_snapshot_1d(x, rho, u, p, it)

    # final solution for plotting
    rho, u, p = cons_to_prim(U, gamma)
    return x, rho, u, p

def run_2d_sod():
    print("\n=== Running 2D Sod Tube Test ===")

    # grid
    ni, nj = 201, 5    # thin strip in y direction
    x, y, dx, dy = make_2d_rect_grid(ni, nj, 0.0, 1.0, 0.0, 0.05)

    # IC
    rhoL, pL, uL = 1.0, 1.0, 0.0
    rhoR, pR, uR = 0.125, 0.1, 0.0
    v0 = 0.0

    U = np.zeros((4, ni, nj))
    for i in range(ni):
        for j in range(nj):
            if x[i, j] < 0.5:
                rho, p, u = rhoL, pL, uL
            else:
                rho, p, u = rhoR, pR, uR
            U[:, i, j] = prim_to_cons(rho, u, v0, p, gamma)

    # time integration
    t      = 0.0
    t_end  = 0.2
    CFL    = 0.1
    it     = 0
    print_every = 10

    while t < t_end:
        dt = compute_dt_cfl_2d(U, dx, dy, gamma, CFL)
        if t + dt > t_end:
            dt = t_end - t

        # --- HERE: define the RHS wrapper ---
        def rhs(U_in):
            return residual_2d(U_in, dx, dy, gamma)

        # --- Do RK3 step using this rhs ---
        U = tvd_rk3_step_2d(U, dt, dx, dy, gamma, rhs)
        t  += dt
        it += 1

        # convert to primitive for debugging
        rho, u, v, p = cons_to_prim(U, gamma)
        velmag = np.sqrt(u*u + v*v)

        if it % print_every == 0 or abs(t - t_end) < 1e-12:
            print(f"[2D SOD] it={it:04d}, t={t:.5f}, dt={dt:.3e}")
            print(f"    rho: [{rho.min():.4f}, {rho.max():.4f}]")
            print(f"    p:   [{p.min():.4f}, {p.max():.4f}]")
            print(f"    vel: [{velmag.min():.4f}, {velmag.max():.4f}]")

            # Stability checks
            if np.isnan(rho).any() or np.isnan(p).any():
                print("ERROR: NaN detected! Simulation blowing up.")
                break

    # return centerline as 1D slice
    jmid = nj // 2
    rho, u, v, p = cons_to_prim(U, gamma)
    return x[:, jmid], rho[:, jmid], u[:, jmid], p[:, jmid]


def run_odwe_cold():


    print("\n=== Running cold ODWE (Euler only) ===")

    # grid & wedge
    ni, nj = 121, 41
    x, y, dx, dy, wall_theta = make_wedge_grid(
        ni, nj,
        x_left=0.0, x_right=1.0,
        y_bottom=0.0, y_top=0.2,
        x_ramp=0.3, theta_deg=15.0
    )

    # inflow conditions (supersonic premix; here just air for cold run)
    M_inf = 5.0
    T_inf = 800.0
    p_inf = 1.0e5
    Rgas  = 287.0

    a_inf   = np.sqrt(gamma * Rgas * T_inf)
    u_inf   = M_inf * a_inf
    v_inf   = 0.0
    rho_inf = p_inf / (Rgas * T_inf)

    inflow_state = (rho_inf, u_inf, v_inf, p_inf)

    # initial condition: uniform inflow state everywhere
    U = np.zeros((4, ni, nj))
    for i in range(ni):
        for j in range(nj):
            U[:, i, j] = prim_to_cons(rho_inf, u_inf, v_inf, p_inf, gamma)

    # time-stepping
    t      = 0.0
    t_end  = 1.0e-4 
    CFL    = 0.4
    it     = 0
    print_every = 5

    while t < t_end:
        # CFL dt
        smax = max_wave_speed_2d(U, gamma)
        dt   = CFL * min(dx, dy) / smax
        if t + dt > t_end:
            dt = t_end - t

        # build RHS via ODWE residual
        R = residual_2d_odwe(U, dx, dy, gamma, Rgas, inflow_state, wall_theta)
        # RK3 stage wrapper
        def rhs(U_in):
            return residual_2d_odwe(U_in, dx, dy, gamma, Rgas, inflow_state, wall_theta)

        U = tvd_rk3_step_2d(U, dt, dx, dy, gamma, rhs_func=rhs)  # see note below
        t += dt
        it += 1

        if it % print_every == 0 or abs(t - t_end) < 1e-12:
            rho, u, v, p = cons_to_prim(U, gamma)
            print(f"[ODWE cold] it={it:04d}, t={t:.5e}, "
                  f"rho:[{rho.min():.3e},{rho.max():.3e}], "
                  f"p:[{p.min():.3e},{p.max():.3e}]")

    return x, y, U



if __name__ == "__main__":
    x, y, U = run_odwe_cold()

    # Convert to primitive vars
    rho, u, v, p = cons_to_prim(U, gamma)
    a = np.sqrt(gamma * 287.0 * (p / (rho * 287.0)))  # or directly sqrt(gamma*p/rho)
    Vmag = np.sqrt(u**2 + v**2)
    M = Vmag / np.sqrt(gamma * p / rho)

    # Mach contour
    plt.figure(figsize=(8, 3))
    levels = np.linspace(0, np.max(M), 40)
    cp = plt.contourf(x, y, M, levels=levels)
    plt.colorbar(cp, label="Mach number")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Cold ODWE - Mach contours")
    plt.tight_layout()
    plt.savefig("odwe_cold_Mach.png", dpi=200)

    j_wall = 2  # a few cells above the wall
    plt.figure()
    plt.plot(x[:, j_wall], p[:, j_wall])
    plt.xlabel("x")
    plt.ylabel("Pressure")
    plt.title("Cold ODWE - pressure along near-wall line")
    plt.grid(True)
    plt.savefig("odwe_cold_pressure.png", dpi=200)
