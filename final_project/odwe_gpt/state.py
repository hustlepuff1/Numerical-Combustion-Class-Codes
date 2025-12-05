# state.py

import numpy as np

def prim_to_cons(rho, u, v, p, gamma):
    """
    Primitive -> conservative for 2D Euler.
    rho, u, v, p can be scalars or arrays of same shape.
    Returns U with shape (4, ...) : [rho, rho*u, rho*v, rho*E]
    """
    E = p / ((gamma - 1.0) * rho) + 0.5 * (u**2 + v**2)
    U = np.zeros((4,) + np.shape(rho))
    U[0] = rho
    U[1] = rho * u
    U[2] = rho * v
    U[3] = rho * E
    return U

def cons_to_prim(U, gamma):
    """
    Conservative -> primitive (2D Euler).
    U shape: (4, N) or (4, ni, nj)
    Returns rho, u, v, p (same shape as rho).
    """
    rho = U[0]
    u   = U[1] / rho
    v   = U[2] / rho
    E   = U[3] / rho
    q2  = 0.5 * (u**2 + v**2)
    p   = (gamma - 1.0) * rho * (E - q2)
    return rho, u, v, p
