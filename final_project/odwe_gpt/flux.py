# flux.py

import numpy as np
from state import cons_to_prim
from bc import slip_wall_reflect
from state import prim_to_cons, cons_to_prim

def compute_flux(U, gamma):
    """
    Physical Euler flux F(U) for 1D.
    U: shape (4, N)
    Returns F with same shape.
    """
    rho, u, p = cons_to_prim(U, gamma)
    E = U[2] / rho

    F = np.zeros_like(U)
    F[0] = rho * u
    F[1] = rho * u**2 + p
    F[2] = u * rho * E + u * p
    # F[3] stays zero for now
    return F

def max_wave_speed(U, gamma):
    rho, u, p = cons_to_prim(U, gamma)
    a = np.sqrt(gamma * p / rho)
    return np.max(np.abs(u) + a)

def rusanov_flux(UL, UR, gamma):
    """
    Local Lax-Friedrichs (Rusanov) numerical flux.
    UL, UR: (4,) arrays at left/right of interface.
    """
    FL = compute_flux(UL[:, None], gamma)[:, 0]
    FR = compute_flux(UR[:, None], gamma)[:, 0]

    # wave speeds
    rhoL, uL, pL = cons_to_prim(UL[:, None], gamma)
    rhoR, uR, pR = cons_to_prim(UR[:, None], gamma)
    aL = np.sqrt(gamma * pL / rhoL)
    aR = np.sqrt(gamma * pR / rhoR)
    smax = np.max([np.abs(uL) + aL, np.abs(uR) + aR])

    return 0.5 * (FL + FR) - 0.5 * smax * (UR - UL)

def residual_1d(U, dx, gamma):
    """
    Compute spatial residual dU/dt = -dF/dx for 1D Euler.
    Simple transmissive BC: zero-gradient at both ends.
    U: (4, N)
    Returns R: (4, N)
    """
    nvar, N = U.shape

    # ghost cells
    U_ext = np.zeros((nvar, N + 2))
    U_ext[:, 1:-1] = U
    U_ext[:, 0]    = U[:, 0]
    U_ext[:, -1]   = U[:, -1]

    # interface fluxes
    F_half = np.zeros((nvar, N + 1))
    for i in range(N + 1):
        UL = U_ext[:, i    ]
        UR = U_ext[:, i + 1]
        F_half[:, i] = rusanov_flux(UL, UR, gamma)

    # residual
    R = np.zeros_like(U)
    for i in range(N):
        R[:, i] = - (F_half[:, i+1] - F_half[:, i]) / dx

    return R


def compute_flux_x(U, gamma):
    """
    F(U) in x-direction for 2D Euler.
    U: (4, ni, nj)
    Returns F: (4, ni, nj)
    """
    rho, u, v, p = cons_to_prim(U, gamma)
    E = U[3] / rho

    F = np.zeros_like(U)
    F[0] = rho * u
    F[1] = rho * u**2 + p
    F[2] = rho * u * v
    F[3] = u * (rho * E + p)
    return F

def compute_flux_y(U, gamma):
    """
    G(U) in y-direction for 2D Euler.
    U: (4, ni, nj)
    Returns G: (4, ni, nj)
    """
    rho, u, v, p = cons_to_prim(U, gamma)
    E = U[3] / rho

    G = np.zeros_like(U)
    G[0] = rho * v
    G[1] = rho * u * v
    G[2] = rho * v**2 + p
    G[3] = v * (rho * E + p)
    return G

def max_wave_speed_2d(U, gamma):
    rho, u, v, p = cons_to_prim(U, gamma)
    a = np.sqrt(gamma * p / rho)
    return np.max(np.sqrt(u**2 + v**2) + a)

def rusanov_flux_normal(UL, UR, n, gamma):
    """
    Rusanov flux in direction given by unit normal n = (nx, ny).
    UL, UR: (4,) left/right state at interface.
    """
    nx, ny = n
    # left/right primitive
    rhoL, uL, vL, pL = cons_to_prim(UL[:, None, None], gamma)
    rhoR, uR, vR, pR = cons_to_prim(UR[:, None, None], gamma)

    uL_n = uL*nx + vL*ny
    uR_n = uR*nx + vR*ny
    aL = np.sqrt(gamma * pL / rhoL)
    aR = np.sqrt(gamma * pR / rhoR)
    smax = float(max(abs(uL_n) + aL, abs(uR_n) + aR))

    # fluxes in x,y separately
    FL = compute_flux_x(UL[:, None, None], gamma)[:, 0, 0]
    GL = compute_flux_y(UL[:, None, None], gamma)[:, 0, 0]
    FR = compute_flux_x(UR[:, None, None], gamma)[:, 0, 0]
    GR = compute_flux_y(UR[:, None, None], gamma)[:, 0, 0]

    FnL = nx*FL + ny*GL
    FnR = nx*FR + ny*GR

    return 0.5*(FnL + FnR) - 0.5*smax*(UR - UL)

def residual_2d(U, dx, dy, gamma):
    """
    2D residual for dU/dt = -dF/dx - dG/dy with Rusanov flux.
    Transmissive (zero-gradient) BCs on all sides.
    U: (4, ni, nj)
    Returns R: (4, ni, nj)
    """
    nvar, ni, nj = U.shape

    # 1) Build ghost cells
    U_ext = np.zeros((nvar, ni + 2, nj + 2))
    U_ext[:, 1:-1, 1:-1] = U

    # left/right copy
    U_ext[:, 0,      1:-1] = U[:, 0,     :]
    U_ext[:, -1,     1:-1] = U[:, -1,    :]
    # bottom/top copy
    U_ext[:, 1:-1,   0   ] = U[:, :, 0   ]
    U_ext[:, 1:-1,  -1   ] = U[:, :, -1  ]
    # corners
    U_ext[:, 0,   0   ] = U[:, 0,   0  ]
    U_ext[:, 0,  -1   ] = U[:, 0,  -1  ]
    U_ext[:, -1,  0   ] = U[:, -1,  0  ]
    U_ext[:, -1, -1   ] = U[:, -1, -1  ]

    # 2) Fluxes at interfaces
    # F_half_x: interfaces in x, shape (nvar, ni+1, nj)
    # G_half_y: interfaces in y, shape (nvar, ni, nj+1)
    F_half_x = np.zeros((nvar, ni + 1, nj))
    G_half_y = np.zeros((nvar, ni, nj + 1))

    # x-interfaces: between (i, j) and (i+1, j) in extended grid
    nx = (1.0, 0.0)
    for i in range(ni + 1):
        for j in range(nj):
            UL = U_ext[:, i,     j + 1]
            UR = U_ext[:, i + 1, j + 1]
            F_half_x[:, i, j] = rusanov_flux_normal(UL, UR, nx, gamma)

    # y-interfaces: between (i, j) and (i, j+1) in extended grid
    ny = (0.0, 1.0)
    for i in range(ni):
        for j in range(nj + 1):
            UL = U_ext[:, i + 1, j    ]
            UR = U_ext[:, i + 1, j + 1]
            G_half_y[:, i, j] = rusanov_flux_normal(UL, UR, ny, gamma)

    # 3) Residual in each cell from flux differences
    R = np.zeros_like(U)
    for i in range(ni):
        for j in range(nj):
            F_ip = F_half_x[:, i + 1, j]   # i+1/2
            F_im = F_half_x[:, i,     j]   # i-1/2
            G_jp = G_half_y[:, i, j + 1]   # j+1/2
            G_jm = G_half_y[:, i, j    ]   # j-1/2

            R[:, i, j] = - (F_ip - F_im) / dx - (G_jp - G_jm) / dy

    return R


# flux.py (add below residual_2d)



def residual_2d_odwe(U, dx, dy, gamma, Rgas, inflow_state, wall_theta):
    """
    Residual for ODWE case with physical BCs:
      - left: fixed supersonic inflow
      - right: extrapolating outflow
      - bottom: slip wall with wedge angle (theta[i])
      - top: extrapolating far field

    inflow_state: (rho_inf, u_inf, v_inf, p_inf)
    wall_theta  : (ni,) wedge angle at bottom wall
    """
    nvar, ni, nj = U.shape
    rho_inf, u_inf, v_inf, p_inf = inflow_state

    # 1) build ghost layer
    U_ext = np.zeros((nvar, ni + 2, nj + 2))
    U_ext[:, 1:-1, 1:-1] = U

    # --- Left boundary (supersonic inflow) ---
    U_inflow = prim_to_cons(rho_inf, u_inf, v_inf, p_inf, gamma)
    for j in range(nj):
        U_ext[:, 0, j+1] = U_inflow

    # --- Right boundary (extrapolating outflow) ---
    U_ext[:, -1, 1:-1] = U[:, -1, :]

    # --- Top boundary (copy / far field) ---
    U_ext[:, 1:-1, -1] = U[:, :, -1]

    # --- Bottom boundary (slip wall with wedge angle) ---
    for i in range(ni):
        # interior cell just above bottom wall
        U_int = U[:, i, 0]
        theta = wall_theta[i]
        U_ghost = slip_wall_reflect(U_int, theta, gamma)
        U_ext[:, i+1, 0] = U_ghost

    # corners (simple copies)
    U_ext[:, 0,   0]  = U_ext[:, 0,   1]
    U_ext[:, 0,  -1]  = U_ext[:, 0,  -2]
    U_ext[:, -1, 0]   = U_ext[:, -1, 1]
    U_ext[:, -1, -1]  = U_ext[:, -1, -2]

    # 2) fluxes at interfaces (same pattern as residual_2d)
    F_half_x = np.zeros((nvar, ni + 1, nj))
    G_half_y = np.zeros((nvar, ni, nj + 1))

    nx = (1.0, 0.0)
    for i in range(ni + 1):
        for j in range(nj):
            UL = U_ext[:, i,     j + 1]
            UR = U_ext[:, i + 1, j + 1]
            F_half_x[:, i, j] = rusanov_flux_normal(UL, UR, nx, gamma)

    ny = (0.0, 1.0)
    for i in range(ni):
        for j in range(nj + 1):
            UL = U_ext[:, i + 1, j    ]
            UR = U_ext[:, i + 1, j + 1]
            G_half_y[:, i, j] = rusanov_flux_normal(UL, UR, ny, gamma)

    # 3) residual in each cell
    R = np.zeros_like(U)
    for i in range(ni):
        for j in range(nj):
            F_ip = F_half_x[:, i + 1, j]
            F_im = F_half_x[:, i,     j]
            G_jp = G_half_y[:, i, j + 1]
            G_jm = G_half_y[:, i, j    ]
            R[:, i, j] = - (F_ip - F_im) / dx - (G_jp - G_jm) / dy

    return R
