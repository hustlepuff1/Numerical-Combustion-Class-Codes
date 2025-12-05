import numpy as np

# ==============================
# 0. Physical / model parameters
# ==============================

gamma = 1.4
R_gas = 287.0            # J/(kg K) (air-like)
cv = R_gas / (gamma - 1)

# Chemistry model (placeholder values – you will tune these)
Q_release = 3.0e6        # J/kg per unit lambda (heat of reaction)
A_preexp = 1.0e6         # 1/s (pre-exponential factor)
E_act = 8.0e4            # J/mol (activation energy)
R_univ = 8.314           # J/(mol K)
MW_mix = 0.029           # kg/mol (air-ish), for rough scaling


# ==============================
# 1. Geometry: nozzle A(x)
# ==============================

def nozzle_area(x, L=1.0, A_exit=2.5):
    """
    Symmetric C-D nozzle with min area = 1 at x = L/2.
    A(x) = 1 + (A_exit - 1) * xi^2,  xi in [-1, 1].
    """
    xi = (x - 0.5 * L) / (0.5 * L)
    return 1.0 + (A_exit - 1.0) * xi**2

def nozzle_dA_dx(x, L=1.0, A_exit=2.5):
    """
    dA/dx = 4*(A_exit - 1)*xi / L, with xi as in nozzle_area.
    """
    xi = (x - 0.5 * L) / (0.5 * L)
    return 4.0 * (A_exit - 1.0) * xi / L


# ========================================
# 2. Conversions between primitive/cons
# ========================================

def cons_from_prim(rho, u, p, lam):
    """
    Primitive -> Conservative, including chemical energy lam*Q_release in E.
    """
    e_thermal = p / ((gamma - 1.0) * rho)
    e_total = e_thermal + lam * Q_release
    E = e_total + 0.5 * u**2
    U1 = rho
    U2 = rho * u
    U3 = rho * E
    return np.array([U1, U2, U3])

def prim_from_cons(U, lam):
    """
    Conservative -> Primitive, given lambda.
    """
    rho = U[0]
    u = U[1] / rho
    E = U[2] / rho
    e_total = E - 0.5 * u**2
    e_thermal = e_total - lam * Q_release
    p = (gamma - 1.0) * rho * e_thermal
    return rho, u, p

def temperature_from_prim(rho, p):
    """
    Ideal-gas temperature: p = rho * R * T
    """
    return p / (rho * R_gas)


# ==============================
# 3. Flux and wave speed
# ==============================

def compute_flux(U, lam):
    """
    Physical flux F(U) for 1D Euler (no area/source terms here).
    """
    rho, u, p = prim_from_cons(U, lam)
    E = U[2] / rho
    F1 = rho * u
    F2 = rho * u**2 + p
    F3 = u * rho * E + u * p
    return np.array([F1, F2, F3])

def max_wave_speed(U, lam):
    rho, u, p = prim_from_cons(U, lam)
    a = np.sqrt(gamma * p / rho)
    return np.abs(u) + a

def rusanov_flux(UL, UR, lamL, lamR):
    """
    Rusanov (local Lax-Friedrichs) numerical flux.
    """
    FL = compute_flux(UL, lamL)
    FR = compute_flux(UR, lamR)
    sL = max_wave_speed(UL, lamL)
    sR = max_wave_speed(UR, lamR)
    smax = max(sL, sR)
    return 0.5 * (FL + FR) - 0.5 * smax * (UR - UL)


# ==============================
# 4. Flow step (no chemistry)
# ==============================

def flow_step(U, lam, x, A, dA_dx, CFL, dt_max=None):
    """
    One explicit Euler step for flow:
    - Rusanov flux
    - Quasi-1D area source term in momentum
    """
    N = x.size

    # CFL time step
    smax = 0.0
    for i in range(N):
        s = max_wave_speed(U[:, i], lam[i])
        smax = max(smax, s)
    dx = x[1] - x[0]
    dt_cfl = CFL * dx / smax
    dt = min(dt_cfl, dt_max) if dt_max is not None else dt_cfl

    # Ghost cells (simple transmissive BCs for now)
    U_ext = np.zeros((3, N + 2))
    lam_ext = np.zeros(N + 2)
    U_ext[:, 1:-1] = U
    lam_ext[1:-1] = lam
    U_ext[:, 0] = U[:, 0]
    lam_ext[0] = lam[0]
    U_ext[:, -1] = U[:, -1]
    lam_ext[-1] = lam[-1]

    # Interface fluxes
    F_half = np.zeros((3, N + 1))
    for i in range(N + 1):
        UL = U_ext[:, i]
        UR = U_ext[:, i + 1]
        lamL = lam_ext[i]
        lamR = lam_ext[i + 1]
        F_half[:, i] = rusanov_flux(UL, UR, lamL, lamR)

    # Cell updates
    U_new = np.zeros_like(U)
    for i in range(N):
        # flux divergence term, note A is at cell centers
        dF = (A[i+1] * F_half[:, i+1] - A[i] * F_half[:, i]) / (dx * A[i])
        # area source term
        S = np.zeros(3)
        rho, u, p = prim_from_cons(U[:, i], lam[i])
        S[1] = p * dA_dx[i] / A[i]
        U_new[:, i] = U[:, i] - dt * dF + dt * S

    return U_new, dt


# ==============================
# 5. Chemistry step
# ==============================

def reaction_rate(T, lam):
    """
    Simple Arrhenius progress-variable model:
    dlam/dt = k(T) * (1 - lam)
    """
    # VERY rough Arrhenius, with MW_mix for scaling.
    # You can tweak this; the main point is T dependence.
    k = A_preexp * np.exp(-E_act / (R_univ * T))  # treat E_act/R_univ/T only
    return k * (1.0 - lam)

def reaction_step(U, lam, dt_reac, n_sub=1):
    """
    Operator-split reaction update:
    - constant rho, u during chemistry step
    - update lambda and energy (via Q_release)
    """
    N = lam.size
    lam_new = lam.copy()
    U_new = U.copy()
    dt_local = dt_reac / n_sub

    for _ in range(n_sub):
        for i in range(N):
            rho, u, p = prim_from_cons(U_new[:, i], lam_new[i])
            T = temperature_from_prim(rho, p)

            dlam_dt = reaction_rate(T, lam_new[i])
            dlam = dt_local * dlam_dt
            lam_star = lam_new[i] + dlam
            lam_star = max(0.0, min(1.0, lam_star))

            # energy update: add heat Q * dλ
            dE_chem = Q_release * (lam_star - lam_new[i])  # per kg
            U_new[2, i] += rho * dE_chem

            lam_new[i] = lam_star

    return U_new, lam_new


# ==============================
# 6. Driver
# ==============================

def run_nozzle_sim():
    # Grid
    N = 201
    L = 1.0
    x = np.linspace(0.0, L, N)
    dx = x[1] - x[0]

    # Geometry (A and dA/dx on cell centers)
    A_center = nozzle_area(x, L=L)
    dA_center = nozzle_dA_dx(x, L=L)

    # For flux divergence we need A at i and i+1.
    # Simple padding at ends:
    A = np.zeros(N + 1)
    A[:-1] = A_center
    A[-1] = A_center[-1]

    # Initial condition (you will refine this)
    p0 = 5.0e5      # Pa
    T0 = 1000.0     # K
    rho0 = p0 / (R_gas * T0)
    u0 = 50.0       # m/s
    lam0 = 0.0      # start unreacted (or 1.0 if your story is different)

    U = np.zeros((3, N))
    lam = np.full(N, lam0)
    for i in range(N):
        U[:, i] = cons_from_prim(rho0, u0, p0, lam[i])

    CFL = 0.5
    t = 0.0
    t_end = 0.01
    it = 0
    output_every = 50

    while t < t_end:
        U, dt_flow = flow_step(U, lam, x, A, dA_center, CFL)
        U, lam = reaction_step(U, lam, dt_flow, n_sub=2)

        t += dt_flow
        it += 1

        if it % output_every == 0:
            print(f"[it={it}] t={t:.6e}, dt={dt_flow:.3e}")

    # Convert to primitive for postprocessing
    rho = np.zeros(N)
    u = np.zeros(N)
    p = np.zeros(N)
    T = np.zeros(N)
    Mach = np.zeros(N)
    for i in range(N):
        rho[i], u[i], p[i] = prim_from_cons(U[:, i], lam[i])
        T[i] = temperature_from_prim(rho[i], p[i])
        a = np.sqrt(gamma * p[i] / rho[i])
        Mach[i] = u[i] / a

    return x, A_center, rho, u, p, T, Mach, lam


if __name__ == "__main__":
    x, A, rho, u, p, T, Mach, lam = run_nozzle_sim()
    # TODO: Add matplotlib plots, e.g.:
    # import matplotlib.pyplot as plt
    # plt.figure(); plt.plot(x, Mach); plt.xlabel("x"); plt.ylabel("Mach"); plt.show()
