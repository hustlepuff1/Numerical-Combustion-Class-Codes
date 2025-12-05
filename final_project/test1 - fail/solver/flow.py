import numpy as np
from solver.source_terms import axisymmetric_sources
from utils.gas_properties import GasProperties

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

        # --- Thermodynamics: NASA CEA-based properties ---
        self.thermo = GasProperties()   # uses your CEA table internally

        # Cell-wise gamma and R (start with something reasonable)
        self.gamma_field = np.full((self.ni, self.nj), 1.25)
        self.R_field     = np.full((self.ni, self.nj), 320.0)

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

        # Approximate grid spacings for CFL
        self.dx = self.grid.length / max((self.ni - 1), 1)
        self.dy = np.mean(self.grid.Y[:, -1]) / max((self.nj - 1), 1)
        
        # --- FIX: Initialize Volume Array for Particle Coupling ---
        self.vol = np.zeros((self.ni, self.nj))
        self._compute_approx_volumes()

        # Boundary condition parameters (reservoir & back pressure)
        self.p_inlet = None
        self.T_inlet = None
        self.p_outlet = None

    def _compute_approx_volumes(self):
        """
        Calculate approximate cell volumes (Annular Rings) for source term scaling.
        V_ij approx = 2 * pi * r * dx * dy
        """
        dx = self.grid.length / max((self.ni - 1), 1)
        
        for i in range(self.ni):
            # Local radial spacing (dy) changes with x because nozzle contours
            y_wall = self.grid.Y[i, -1]
            dy = y_wall / max((self.nj - 1), 1)
            
            # Radius at each node
            r = np.abs(self.grid.Y[i, :])
            # Safety clamp for r to avoid zero volume at centerline
            r_safe = np.maximum(r, 1e-6)
            
            # Volume = Circumference * Area
            self.vol[i, :] = 2.0 * np.pi * r_safe * dx * dy

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------
    def initialize_solution(self, p_chamber, t_chamber, p_exit):
        self.p_inlet  = p_chamber
        self.T_inlet  = t_chamber
        self.p_outlet = p_exit

        # Uniform stagnant reservoir in the interior
        self.p[:, :] = p_chamber
        self.T[:, :] = t_chamber
        self.rho[:, :] = self.p / (self.R * self.T)
        self.u[:, :] = 0.0
        self.v[:, :] = 0.0

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
        """
        Update primitive variables (u, v, p, T) from conserved ones,
        using *cell-wise* gamma(T) and R(T) based on NASA CEA data.
        """
        eps = 1e-9

        # Density must stay positive
        self.rho = np.clip(self.rho, 1e-4, 1e4)

        # Velocities
        self.u = self.mom_x / (self.rho + eps)
        self.v = self.mom_y / (self.rho + eps)

        # Kinetic energy
        kinetic = 0.5 * self.rho * (self.u**2 + self.v**2)

        # Pressure from E = p/(γ-1) + ½ρu²
        self.p = (self.gamma_field - 1.0) * (self.energy - kinetic)
        # keep pressure positive and bounded
        self.p = np.clip(self.p, 1.0e3, 1.0e9)

        # Temperature from p = ρ R T (cell-wise R)
        self.T = self.p / (self.rho * self.R_field + eps)
        self.T = np.clip(self.T, 200.0, 5000.0)  # don't let T go crazy

        # --- Update gamma and R from NASA CEA table ---
        gamma_new, R_new = self.thermo.get_properties(self.T)

        # Mild under-relaxation for stability
        relax = 0.5
        self.gamma_field = (1.0 - relax) * self.gamma_field + relax * gamma_new
        self.R_field     = (1.0 - relax) * self.R_field     + relax * R_new


    # ------------------------------------------------------------------
    # Time step (CFL)
    # ------------------------------------------------------------------
    def get_dt(self, cfl=0.5):
        """
        Compute time step from CFL condition using local γ(T), R(T).
        """
        eps = 1e-9
        # Primitives should be up to date
        self.update_primitive_variables()

        a  = np.sqrt(self.gamma_field * self.p / (self.rho + eps))
        ws = np.sqrt(self.u**2 + self.v**2) + a
        ws = np.nan_to_num(ws, nan=1.0, posinf=1.0e6, neginf=1.0e6)

        max_ws = np.max(ws) + 1e-9
        dx_min = min(self.dx, self.dy)

        return cfl * dx_min / max_ws

    # ------------------------------------------------------------------
    # Main update (Rusanov FVM)
    # ------------------------------------------------------------------
    def step(self, dt, coupling_terms=None):
        """
        Perform one time step of the Euler equations.
        coupling_terms: Dict of source arrays from ParticleSolver (rho, mom_x, etc.)
        """

        # ---------- I-direction (axial) flux ----------
        F_rho   = self.rho * self.u
        F_mom_x = self.rho * self.u**2 + self.p
        F_mom_y = self.rho * self.u * self.v
        F_eng   = (self.energy + self.p) * self.u

        eps = 1e-9
        # sound speed using local γ(T), not constant γ
        a = np.sqrt(self.gamma_field * np.maximum(self.p, 1e-5) / (self.rho + eps))
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

        # --- Axisymmetric Source Terms ---
        S_rho, S_mx, S_my, S_en = axisymmetric_sources(
            self.rho,
            self.u,
            self.v,
            self.p,
            self.energy,
            self.grid.Y,
            r_min=1e-3,
            clip_val=1e5,
            scale=1.0,
            use_energy=True
        )

        # --- Add Particle Coupling Terms ---
        if coupling_terms is not None:
            # Particle source terms are already Density Rates (kg/m^3/s)
            # if they were divided by gas.vol in particles.py
            if 'rho' in coupling_terms: S_rho += coupling_terms['rho']
            if 'mom_x' in coupling_terms: S_mx += coupling_terms['mom_x']
            if 'mom_y' in coupling_terms: S_my += coupling_terms['mom_y']
            if 'energy' in coupling_terms: S_en += coupling_terms['energy']

        # Apply Total Source Terms
        self.rho[1:-1, 1:-1]   += dt * S_rho[1:-1, 1:-1]
        self.mom_x[1:-1, 1:-1] += dt * S_mx[1:-1, 1:-1]
        self.mom_y[1:-1, 1:-1] += dt * S_my[1:-1, 1:-1]
        self.energy[1:-1, 1:-1] += dt * S_en[1:-1, 1:-1]

        # Apply BCs and update primitives
        self.apply_boundary_conditions()
        self.update_primitive_variables()

    # ------------------------------------------------------------------
    # Boundary conditions
    # ------------------------------------------------------------------
    def apply_boundary_conditions(self):
        eps = 1e-9
        # Inlet
        if self.p_inlet is not None and self.T_inlet is not None:
            rho_in = self.p_inlet / (self.R * self.T_inlet)
            u_in = 50.0 # Diode
            v_in = 0.0
            E_in = self.p_inlet / (self.gamma - 1.0) + 0.5 * rho_in * (u_in**2 + v_in**2)

            self.rho[0, :]   = rho_in
            self.mom_x[0, :] = rho_in * u_in
            self.mom_y[0, :] = rho_in * v_in
            self.energy[0, :] = E_in

        # Outlet
        if self.p_outlet is not None:
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

        # Centerline
        self.rho[:, 0]    = self.rho[:, 1]
        self.mom_x[:, 0]  = self.mom_x[:, 1]
        self.mom_y[:, 0]  = -self.mom_y[:, 1]
        self.energy[:, 0] = self.energy[:, 1]

        # Wall
        self.rho[:, -1]    = self.rho[:, -2]
        self.mom_x[:, -1]  = self.mom_x[:, -2]
        self.mom_y[:, -1]  = 0.0
        self.energy[:, -1] = self.energy[:, -2]