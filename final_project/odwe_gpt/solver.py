# solver.py

import numpy as np
from flux import residual_1d, max_wave_speed
from flux import residual_2d, max_wave_speed_2d

def tvd_rk3_step(U, dt, dx, gamma):
    """
    Shu-Osher TVD RK3 for 1D Euler.
    U: (4, N)
    """
    # stage 1
    R1 = residual_1d(U, dx, gamma)
    U1 = U + dt * R1

    # stage 2
    R2 = residual_1d(U1, dx, gamma)
    U2 = 0.75 * U + 0.25 * (U1 + dt * R2)

    # stage 3
    R3 = residual_1d(U2, dx, gamma)
    U3 = (1.0/3.0) * U + (2.0/3.0) * (U2 + dt * R3)

    return U3

def compute_dt_cfl(U, dx, gamma, CFL):
    """
    Compute time step from CFL condition.
    """
    smax = max_wave_speed(U, gamma)
    return CFL * dx / smax


def tvd_rk3_step_2d(U, dt, dx, dy, gamma, rhs_func):
    """
    Generic TVD RK3 with user-supplied RHS(U).
    rhs_func(U) must return same shape as U.
    """
    R1 = rhs_func(U)
    U1 = U + dt * R1

    R2 = rhs_func(U1)
    U2 = 0.75 * U + 0.25 * (U1 + dt * R2)

    R3 = rhs_func(U2)
    U3 = (1.0/3.0)*U + (2.0/3.0)*(U2 + dt * R3)
    return U3

def compute_dt_cfl_2d(U, dx, dy, gamma, CFL):
    smax = max_wave_speed_2d(U, gamma)
    # simple choice: use min(dx,dy) for CFL
    h = min(dx, dy)
    return CFL * h / smax