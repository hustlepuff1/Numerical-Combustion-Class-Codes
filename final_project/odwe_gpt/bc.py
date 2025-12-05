# bc.py

import numpy as np
from state import cons_to_prim, prim_to_cons

def slip_wall_reflect(U_cell, theta, gamma):
    """
    Given conservative state U_cell at a boundary-adjacent cell,
    reflect the velocity at a wall with angle theta (radians).

    theta = 0   -> horizontal wall (along x-axis)
    theta > 0   -> ramp up (compression corner)

    Returns conservative state for the ghost cell.
    """
    rho, u, v, p = cons_to_prim(U_cell[:, None, None], gamma)

    # Convert to scalars
    rho = float(rho); u = float(u); v = float(v); p = float(p)

    # Wall tangent & normal
    # tangent t = (cosθ, sinθ)
    # normal  n = (-sinθ, cosθ) pointing into the fluid
    t = np.array([np.cos(theta), np.sin(theta)])
    n = np.array([-np.sin(theta), np.cos(theta)])

    v_vec = np.array([u, v])

    # Decompose into normal and tangential components
    v_n = np.dot(v_vec, n) * n
    v_t = v_vec - v_n

    # Reflect normal component: v_n -> -v_n
    v_ref = v_t - v_n
    u_ref, v_ref_y = v_ref[0], v_ref[1]

    # Build conservative ghost state (same rho, p)
    U_ghost = prim_to_cons(rho, u_ref, v_ref_y, p, gamma)
    return U_ghost
