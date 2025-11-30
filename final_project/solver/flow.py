import numpy as np

class GasFlowSolver:
    """
    Eulerian Solver for Compressible Gas Flow.
    Solves the inviscid Navier-Stokes equations (Euler) using a
    finite-volume Rusanov scheme on a 2D nozzle grid.
    """

    def __init__(self, grid, gamma=1.4, R=287.0):
        self.grid = grid
        self.gamma = gamma
        self.R = R
        self.ni = grid.ni
        self.nj = grid.nj

        # State Variables (Conserved)
        self.rho   = np.zeros((self.ni, self.nj))
        self.mom_x = np.zeros((self.ni, self.nj))
        self.mom_y = np.zeros((self.ni, self.nj))
        self.energy = np.zeros((self.ni, self.nj))

        # Primitive Variables
        self.u = np.zeros((self.ni, self.nj))
        self.v = np.zeros((self.ni, self.nj))
        self.p = np.zeros((self.ni, self.nj))
        self.T = np.zeros((self.ni, self.nj))

        # Boundary condition parameters (reservoir & back pressure)
        self.p_inlet = None
        self.T_inlet = None
        self.p_outlet = None

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------
    def initialize_solution(self, p_chamber, t_chamber, p_exit):
        """
        Initialize the flow field as a uniform stagnant reservoir at
        chamber conditions. The outlet boundary will enforce p_exit,
        which drives the nozzle to choke and go supersonic.
        """
        self.p_inlet  = p_chamber
        self.T_inlet  = t_chamber
        self.p_outlet = p_exit

        # Uniform stagnant reservoir in the interior
        self.p[:, :] = p_chamber
        self.T[:, :] = t_chamber
        self.rho[:, :] = self.p / (self.R * self.T)
        self.u[:, :] = 0.0
        self.v[:, :] = 0.0

        # Build conserved variables and sync primitives
        self.update_conserved_variables()
        self.update_primitive_variables()

    # ------------------------------------------------------------------
    # Primitive / Conserved conversions
    # ------------------------------------------------------------------
    def update_conserved_variables(self):
        self.mom_x = self.rho * self.u
        self.mom_y = self.rho * self.v
        kinetic = 0.5 * self.rho * (self.u**2 + self.v**2)
        self.energy = self.p / (self.gamma - 1.0) + kinetic

    def update_primitive_variables(self):
        eps = 1e-9
        self.u = self.mom_x / (self.rho + eps)
        self.v = self.mom_y / (self.rho + eps)
        kinetic = 0.5 * self.rho * (self.u**2 + self.v**2)
        self.p = (self.gamma - 1.0) * (self.energy - kinetic)
        self.T = self.p / (self.rho * self.R + eps)

    # ------------------------------------------------------------------
    # Time step (CFL)
    # ------------------------------------------------------------------
    def get_dt(self, cfl=0.5):
        """
        Compute stable timestep using a CFL estimate.
        """
        a = np.sqrt(self.gamma * np.maximum(self.p, 1e-5) / (self.rho + 1e-9))
        vel = np.sqrt(self.u**2 + self.v**2)
        wave_speed = vel + a

        dx = self.grid.length / max((self.ni - 1), 1)
        max_wave = np.max(wave_speed) + 1e-9

        return cfl * dx / max_wave

    # ------------------------------------------------------------------
    # Main update (Rusanov FVM)
    # ------------------------------------------------------------------
    def step(self, dt):
        """
        Perform one time step of the Euler equations.
        """

        # ---------- I-direction (axial) flux ----------
        F_rho   = self.rho * self.u
        F_mom_x = self.rho * self.u**2 + self.p
        F_mom_y = self.rho * self.u * self.v
        F_eng   = (self.energy + self.p) * self.u

        a = np.sqrt(self.gamma * np.maximum(self.p, 1e-5) / (self.rho + 1e-9))
        lambda_max_i = np.abs(self.u) + a

        def rusanov_flux_i(U, F, lam):
            U_L = U[:-1, :]
            U_R = U[1:, :]
            F_L = F[:-1, :]
            F_R = F[1:, :]
            lam_face = 0.5 * (lam[:-1, :] + lam[1:, :])
            return 0.5 * (F_L + F_R) - 0.5 * lam_face * (U_R - U_L)

        flux_rho_i   = rusanov_flux_i(self.rho,   F_rho,   lambda_max_i)
        flux_mom_x_i = rusanov_flux_i(self.mom_x, F_mom_x, lambda_max_i)
        flux_mom_y_i = rusanov_flux_i(self.mom_y, F_mom_y, lambda_max_i)
        flux_eng_i   = rusanov_flux_i(self.energy, F_eng,  lambda_max_i)

        # ---------- J-direction (radial) flux ----------
        G_rho   = self.rho * self.v
        G_mom_x = self.rho * self.v * self.u
        G_mom_y = self.rho * self.v**2 + self.p
        G_eng   = (self.energy + self.p) * self.v

        lambda_max_j = np.abs(self.v) + a

        def rusanov_flux_j(U, G, lam):
            U_L = U[:, :-1]
            U_R = U[:, 1:]
            G_L = G[:, :-1]
            G_R = G[:, 1:]
            lam_face = 0.5 * (lam[:, :-1] + lam[:, 1:])
            return 0.5 * (G_L + G_R) - 0.5 * lam_face * (U_R - U_L)

        flux_rho_j   = rusanov_flux_j(self.rho,   G_rho,   lambda_max_j)
        flux_mom_x_j = rusanov_flux_j(self.mom_x, G_mom_x, lambda_max_j)
        flux_mom_y_j = rusanov_flux_j(self.mom_y, G_mom_y, lambda_max_j)
        flux_eng_j   = rusanov_flux_j(self.energy, G_eng,  lambda_max_j)

        # ---------- Finite volume update ----------
        dx = self.grid.length / max((self.ni - 1), 1)
        dy_avg = np.mean(self.grid.Y[:, -1]) / max((self.nj - 1), 1)

        dF_rho = flux_rho_i[1:,   1:-1] - flux_rho_i[:-1, 1:-1]
        dG_rho = flux_rho_j[1:-1, 1:  ] - flux_rho_j[1:-1, :-1]

        dF_mx = flux_mom_x_i[1:,   1:-1] - flux_mom_x_i[:-1, 1:-1]
        dG_mx = flux_mom_x_j[1:-1, 1:  ] - flux_mom_x_j[1:-1, :-1]

        dF_my = flux_mom_y_i[1:,   1:-1] - flux_mom_y_i[:-1, 1:-1]
        dG_my = flux_mom_y_j[1:-1, 1:  ] - flux_mom_y_j[1:-1, :-1]

        dF_en = flux_eng_i[1:,   1:-1] - flux_eng_i[:-1, 1:-1]
        dG_en = flux_eng_j[1:-1, 1:  ] - flux_eng_j[1:-1, :-1]

        self.rho[1:-1, 1:-1]   -= dt * (dF_rho/dx + dG_rho/dy_avg)
        self.mom_x[1:-1, 1:-1] -= dt * (dF_mx/dx  + dG_mx/dy_avg)
        self.mom_y[1:-1, 1:-1] -= dt * (dF_my/dx  + dG_my/dy_avg)
        self.energy[1:-1, 1:-1]-= dt * (dF_en/dx  + dG_en/dy_avg)

        # Apply BCs and update primitives
        self.apply_boundary_conditions()
        self.update_primitive_variables()

    # ------------------------------------------------------------------
    # Boundary conditions
    # ------------------------------------------------------------------
    def apply_boundary_conditions(self):
        """
        Inlet: fixed reservoir at (p_inlet, T_inlet, u=0, v=0)
        Outlet: fixed static pressure p_outlet (subsonic outlet BC)
        Centerline: symmetry
        Wall: slip wall
        """
        eps = 1e-9
        # ---------- Inlet (i = 0): Reservoir ----------
        if self.p_inlet is not None and self.T_inlet is not None:
            rho_in = self.p_inlet / (self.R * self.T_inlet)
            u_in   = 0.0
            v_in   = 0.0
            E_in   = self.p_inlet / (self.gamma - 1.0) + 0.5 * rho_in * (u_in**2 + v_in**2)

            self.rho[0, :]   = rho_in
            self.mom_x[0, :] = rho_in * u_in
            self.mom_y[0, :] = rho_in * v_in
            self.energy[0, :] = E_in

        # ---------- Outlet (i = ni-1): Fixed static pressure ----------
        if self.p_outlet is not None:
            # Extrapolate density and velocity from interior
            rho_R   = self.rho[-2, :]
            momx_R  = self.mom_x[-2, :]
            momy_R  = self.mom_y[-2, :]
            u_R = momx_R / (rho_R + eps)
            v_R = momy_R / (rho_R + eps)

            p_R = self.p_outlet
            E_R = p_R / (self.gamma - 1.0) + 0.5 * rho_R * (u_R**2 + v_R**2)

            self.rho[-1, :]   = rho_R
            self.mom_x[-1, :] = rho_R * u_R
            self.mom_y[-1, :] = rho_R * v_R
            self.energy[-1, :] = E_R

        # ---------- Centerline (j = 0): Symmetry ----------
        self.rho[:, 0]    = self.rho[:, 1]
        self.mom_x[:, 0]  = self.mom_x[:, 1]
        self.mom_y[:, 0]  = -self.mom_y[:, 1]  # v -> -v so average is 0
        self.energy[:, 0] = self.energy[:, 1]

        # ---------- Wall (j = nj-1): Slip wall ----------
        self.rho[:, -1]    = self.rho[:, -2]
        self.mom_x[:, -1]  = self.mom_x[:, -2]
        self.mom_y[:, -1]  = 0.0
        self.energy[:, -1] = self.energy[:, -2]
