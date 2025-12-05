import numpy as np

def axisymmetric_sources(rho, u, v, p, energy, Y,
                         r_min=1e-3, clip_val=1e5, scale=0.3,
                         use_energy=False):
    """
    Axisymmetric geometric source terms for (x,r) coordinates, no swirl.

    Governing equations in conservative form (simplified):

        Continuity:
            ρ_t + ... + ρ v / r = 0
            -> S_rho = - ρ v / r

        Axial momentum:
            (ρu)_t + ... + ρ u v / r = 0
            -> S_mx = - ρ u v / r

        Radial momentum:
            (ρv)_t + ... + (ρ v^2 - p) / r = 0
            -> S_my = - (ρ v^2 - p) / r

        Energy (optional):
            E_t + ... + (E + p) v / r = 0
            -> S_en = - (E + p) v / r

    We return S_rho, S_mx, S_my, S_en with some clipping and an overall
    scaling factor `scale` for stability.
    """
    # Safe radius to avoid division by zero near the axis
    r_safe = np.maximum(np.abs(Y), r_min)

    # Continuity and axial momentum depend on v
    S_rho = - (rho * v) / r_safe
    S_mx  = - (rho * u * v) / r_safe

    # Radial momentum has a pressure term that can generate v even if v=0
    S_my  = - (rho * v * v - p) / r_safe

    if use_energy:
        S_en = - (energy + p) * v / r_safe
    else:
        S_en = np.zeros_like(rho)

    # Symmetry on centerline: no geometric source at r = 0
    S_rho[:, 0] = 0.0
    S_mx[:,  0] = 0.0
    S_my[:,  0] = 0.0
    S_en[:,  0] = 0.0

    if clip_val is not None:
        S_rho = np.clip(S_rho, -clip_val, clip_val)
        S_mx  = np.clip(S_mx,  -clip_val, clip_val)
        S_my  = np.clip(S_my,  -clip_val, clip_val)
        S_en  = np.clip(S_en,  -clip_val, clip_val)

    # Global scaling to keep things tame (tune if needed)
    S_rho *= scale
    S_mx  *= scale
    S_my  *= scale
    S_en  *= scale

    return S_rho, S_mx, S_my, S_en
