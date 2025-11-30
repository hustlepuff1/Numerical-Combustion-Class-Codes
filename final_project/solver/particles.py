import numpy as np

class ParticleSolver:
    """
    Lagrangian particle tracker for molten aluminum / Al2O3 droplets.

    Phase 2, Step 1: 'Cold' particles with drag only (no combustion, no
    mass/heat coupling yet). Particles are advanced in time inside a
    gas field provided by GasFlowSolver.

    Positions are stored in (x, r) coordinates consistent with grid.X, grid.Y.
    """

    def __init__(self, grid,
                 rho_p=1500.0,   # particle density (kg/m^3), ~aluminum
                 Cd=2.0,         # drag coefficient (order 1)
                 d0=80e-6):     # default diameter 150 micron
        self.grid = grid
        self.rho_p = rho_p
        self.Cd = Cd

        # Particle state arrays (1D: each index = one particle)
        self.x = np.array([])      # axial position
        self.r = np.array([])      # radial position
        self.u = np.array([])      # axial velocity
        self.v = np.array([])      # radial velocity
        self.d = np.array([])      # diameter
        self.alive = np.array([])  # bool mask for active particles

        self.default_diameter = d0

    # ------------------------------------------------------------------
    # Particle injection
    # ------------------------------------------------------------------
    def inject_from_inlet(self, n_new, u0=0.0, v0=0.0):
        """
        Inject n_new particles slightly downstream of the inlet boundary
        (in the first interior cell, i=1), so they feel the accelerating
        gas flow instead of the stagnant reservoir BC at i=0.
        """
        # Use first interior cell i = 1
        i_inj = 1
        x0 = self.grid.X[i_inj, 0]
        r_wall = self.grid.Y[i_inj, -1]

        # New particle positions
        x_new = np.full(n_new, x0)           # exactly at that interior x
        r_new = np.random.rand(n_new) * r_wall

        u_new = np.full(n_new, u0 if u0 != 0.0 else 50.0)
        v_new = np.full(n_new, v0)
        d_new = np.full(n_new, self.default_diameter)
        alive_new = np.ones(n_new, dtype=bool)

        # Append to existing arrays
        self.x = np.concatenate([self.x, x_new]) if self.x.size else x_new
        self.r = np.concatenate([self.r, r_new]) if self.r.size else r_new
        self.u = np.concatenate([self.u, u_new]) if self.u.size else u_new
        self.v = np.concatenate([self.v, v_new]) if self.v.size else v_new
        self.d = np.concatenate([self.d, d_new]) if self.d.size else d_new
        self.alive = np.concatenate([self.alive, alive_new]) if self.alive.size else alive_new

    # ------------------------------------------------------------------
    # Helper: interpolate gas field at particle positions
    # ------------------------------------------------------------------
    def _sample_gas(self, gas):
        """
        Sample gas velocity and density at each particle position using
        nearest-cell lookup for simplicity.

        gas: instance of GasFlowSolver
        Returns (u_g, v_g, rho_g, p_g)
        """
        X = gas.grid.X
        Y = gas.grid.Y

        # Find nearest grid indices in i (x) and j (r)
        # (vectorized nearest-neighbor)
        i_indices = np.argmin(np.abs(X[:, 0][:, None] - self.x[None, :]), axis=0)
        j_indices = np.zeros_like(i_indices)

        # For each i, pick nearest j in radius
        for idx, i in enumerate(i_indices):
            r_col = Y[i, :]
            j_indices[idx] = np.argmin(np.abs(r_col - self.r[idx]))

        u_g = gas.u[i_indices, j_indices]
        v_g = gas.v[i_indices, j_indices]
        rho_g = gas.rho[i_indices, j_indices]
        p_g = gas.p[i_indices, j_indices]

        return u_g, v_g, rho_g, p_g

    # ------------------------------------------------------------------
    # Time integration
    # ------------------------------------------------------------------
    def advance(self, dt, gas, g=0.0):
        """
        Advance particles by dt using a simple explicit scheme with drag.

        gas: GasFlowSolver (provides u, v, rho, p)
        g  : gravitational acceleration (m/s^2) in -r direction, default 0.
        """
        if self.x.size == 0:
            return  # no particles yet

        alive = self.alive

        # Sample gas field at particle locations
        u_g, v_g, rho_g, _ = self._sample_gas(gas)

        # Relative velocity
        du = u_g - self.u
        dv = v_g - self.v
        mag_rel = np.sqrt(du**2 + dv**2) + 1e-12

        # Drag acceleration (spherical particle)
        # F_d = 0.5 * Cd * rho_g * A * |rel| * rel_dir
        # a_p = F_d / m_p = (3 Cd rho_g |rel|)/(4 rho_p d) * rel
        coeff = 3.0 * self.Cd * rho_g * mag_rel / (4.0 * self.rho_p * (self.d + 1e-12))
        a_u = coeff * du
        a_v = coeff * dv - g  # gravity along -r if g>0

        # Explicit Euler update
        self.u[alive] += a_u[alive] * dt
        self.v[alive] += a_v[alive] * dt

        self.x[alive] += self.u[alive] * dt
        self.r[alive] += self.v[alive] * dt

        # Simple particle "death" conditions: outside domain
        L = self.grid.length
        r_max_local = np.interp(self.x, self.grid.X[:, 0], self.grid.Y[:, -1])

        out_left  = self.x < 0.0
        out_right = self.x > L
        out_rad   = (self.r < 0.0) | (self.r > r_max_local + 1e-6)

        self.alive[out_left | out_right | out_rad] = False

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------
    def get_alive_particles(self):
        """Return views of active particles only."""
        mask = self.alive
        return (self.x[mask], self.r[mask],
                self.u[mask], self.v[mask],
                self.d[mask])
