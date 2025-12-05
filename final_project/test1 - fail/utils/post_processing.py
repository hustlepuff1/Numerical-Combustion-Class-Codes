import numpy as np

def calculate_mach(solver):
    """
    Calculates the Mach Number field from the solver state.
    Used for validation against NASA Paper Figure 3.
    """
    epsilon = 1e-9
    
    a = np.sqrt(solver.gamma * np.maximum(solver.p, 1e-5) / (solver.rho + epsilon))
    vel_mag = np.sqrt(solver.u**2 + solver.v**2)
    
    return vel_mag / (a + epsilon)
