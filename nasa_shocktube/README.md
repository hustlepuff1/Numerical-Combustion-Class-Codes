# Development and Validation of a Compressible Euler Solver for Supersonic Hydrogen–Air Combustion Using a 7-Species Reduced Reaction Mechanism

![Fortran](https://img.shields.io/badge/Fortran-Modern-purple)
![Python](https://img.shields.io/badge/Python-Visualization-blue)
![Validation](https://img.shields.io/badge/Validation-NASA_NPARC-green)

This project implements a parallelized 2D compressible flow solver for supersonic hydrogen–air combustion. It was developed as a final project for the graduate-level course **"Numerical Combustion"**.

The solver is validated against the **NASA NPARC "Channel Combustion" Study #1**, demonstrating the ability to capture shock waves, ignition delay, and species evolution in high-speed flows.

---

## 📖 Physics & Mathematical Formulation

### 1. Governing Equations (Vector Form)
The solver integrates the **Reactive Compressible Euler Equations** using a Finite Volume Method (FVM). The system is solved in strong conservative vector form:

$$
\frac{\partial \mathbf{Q}}{\partial t} + \frac{\partial \mathbf{F}(\mathbf{Q})}{\partial x} + \frac{\partial \mathbf{G}(\mathbf{Q})}{\partial y} = \mathbf{S}(\mathbf{Q})
$$

Where the state vector $\mathbf{Q}$, inviscid fluxes $\mathbf{F}, \mathbf{G}$, and source vector $\mathbf{S}$ are defined as:

$$
\mathbf{Q} = \begin{bmatrix} \rho \\ \rho u \\ \rho v \\ \rho E \\ \rho Y_1 \\ \dots \\ \rho Y_{N_s} \end{bmatrix}
$$
$$
\mathbf{F} = \begin{bmatrix} \rho u \\ \rho u^2 + p \\ \rho u v \\ u(\rho E + p) \\ \rho u Y_1 \\ \dots \\ \rho u Y_{N_s} \end{bmatrix}
$$
$$
\mathbf{G} = \begin{bmatrix} \rho v \\ \rho v u \\ \rho v^2 + p \\ v(\rho E + p) \\ \rho v Y_1 \\ \dots \\ \rho v Y_{N_s} \end{bmatrix}
$$
$$
\mathbf{S} = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \\ \dot{\omega}_1 \\ \dots \\ \dot{\omega}_{N_s} \end{bmatrix}
$$

* **Conservation:** Mass, X-Momentum, Y-Momentum, Total Energy.
* **Species:** Transport equations for $N_s$ species with chemical source terms $\dot{\omega}_k$.

### 2. Numerical Flux (Rusanov Scheme)
To capture strong shocks stably, the **Rusanov (Local Lax-Friedrichs)** scheme is used as an approximate Riemann solver at cell interfaces:

$$
\hat{\mathbf{F}}_{i+1/2} = \frac{1}{2} \left( \mathbf{F}_L + \mathbf{F}_R \right) - \frac{1}{2} \lambda_{\max} \left( \mathbf{Q}_R - \mathbf{Q}_L \right)
$$

Where the maximum local wave speed is:

$$
\lambda_{\max} = \max(|u_L| + a_L, |u_R| + a_R)
$$

### 3. Thermodynamics & Equation of State
The flow is modeled as a **Thermally Perfect Gas Mixture**:
* **NASA Polynomials:** Specific heats $C_p(T)$ and enthalpies $h(T)$ are calculated dynamically for each species using **NASA 7-Coefficient Polynomials** (from GRI-Mech 3.0).
* **Temperature Recovery:** Temperature $T$ is recovered from the conserved internal energy $e$ using a **Newton-Raphson iterative solver**.
* **Note:** A constant $\gamma = 1.4$ is used *only* for the acoustic wave speed calculation ($a = \sqrt{\gamma p / \rho}$) in the numerical flux dissipation term.

### 4. Chemical Kinetics
* **Mechanism:** A **7-species, 8-reaction reduced mechanism** for $H_2$-Air combustion.
* **Sources:** Rate coefficients are derived from the benchmark configurations of Evans & Schexnayder (1980) and Mani et al. (1991).
* **Implementation:** Arrhenius rates are integrated using **Operator Splitting** with sub-cycling (20 sub-steps per fluid step) to handle stiffness.

---

## 📂 Project Structure & Input

* **Grid:** Structured Cartesian grid ($N_i \times N_j$) covering domain $[0, L] \times [0, H]$.
* **State Vector:** 10 Variables $\to$ $[\rho, \rho u, \rho v, \rho E, \rho Y_{H_2}, \rho Y_{O_2}, \rho Y_{H_2O}, \rho Y_{OH}, \rho Y_{O}, \rho Y_{H}]$.
    * $N_2$ *is calculated as the inert balance species.*

### Input File (`input.txt`)
The solver reads simulation parameters from `input.txt` in the following strict order:

| Line | Parameter | Description |
| :--- | :--- | :--- |
| 1 | `ni` | Number of grid points in X (e.g., 201) |
| 2 | `nj` | Number of grid points in Y (e.g., 15) |
| 3 | `L` | Domain Length (meters) |
| 4 | `H` | Domain Height (meters) |
| 5 | `P_init` | Initial Pressure (Pa) |
| 6 | `T_init` | Initial Temperature (K) |
| 7 | `M_init` | Initial Mach Number |
| 8 | `Y_H2` | Initial Mass Fraction of $H_2$ |
| 9 | `Y_O2` | Initial Mass Fraction of $O_2$ |
| 10 | `t_final` | Final Simulation Time (s) |
| 11 | `max_steps` | Maximum Iterations |
| 12 | `chem_switch` | 0 = Frozen Flow, 1 = Finite-Rate Chemistry |
| 13 | `CFL` | CFL Number for Time Step Stability |

---

## 🛠 Core Routines

| File | Description |
| :--- | :--- |
| **`main.f90`** | Driver loop. Calls initialization, time stepping, and I/O routines. |
| **`solver_mod.f90`** | Contains the core PDE solver: `compute_residuals` (FVM), `update_solution` (Euler), `apply_boundary_conditions`. |
| **`flux_mod.f90`** | Computes inviscid fluxes using the Rusanov Approximate Riemann Solver. |
| **`thermo_mod.f90`** | Handles thermodynamics: NASA polynomials and Newton-Raphson $T$ recovery. |
| **`reaction_mod.f90`** | Computes chemical source terms using Arrhenius kinetics and handles sub-cycling. |
| **`read_input`** | Parses `input.txt` and initializes the grid/history file. |

---

## 🧪 Validation Benchmark

**Case:** NASA NPARC "Channel Combustion" Study #1.
* **Scenario:** Supersonic ($M=1.82$), Premixed ($H_2/O_2$) flow in a 2D channel.
* **Conditions:** $T_{in} = 1388$ K, $P_{in} = 1$ atm.
* **Targets:**
    * Ignition location ($\approx 0.02$ m).
    * Vertical ignition front topology.
    * Steady-state species mass fractions (Centerline).

---

## 🚀 Usage

### 1. Build
Prerequisites: `CMake`, Fortran Compiler (`ifx` or `gfortran`).
```bash
sh run.sh
````

*This script compiles the project, runs the solver with OpenMP (8 threads), and organizes output files.*

### 2\. Visualization

Python scripts provided for post-processing:

```bash
# Plot latest contours
python3 plot_results.py

# Plot contours at specific frame (e.g., flow_009.dat)
python3 plot_results.py 9
```

*Generates contours (`.png`) and unsteady animations (`.mp4`).*

-----

## 🔗 References

1.  NPARC Alliance Validation Archive, "Channel Combustion: Study \#1."
2.  Mani, M., Bush, R.H., Vogel, P.G., "Implicit Equilibrium and Finite-Rate Chemistry Models For High Speed Flow Applications," *AIAA Paper 91-3299-CP*, Jan. 1991.
3.  Jachimowski, C. J., "An Analytical Study of the Hydrogen-Air Reaction Mechanism With Application to Scramjet Combustion," *NASA TP-2791*, 1988.
4.  Evans, J. S., and Schexnayder, C. J., "Influence of Chemical Kinetics and Unmixedness on Burning in Supersonic Hydrogen Flames," *AIAA Journal*, Vol.18, No. 2, 1980.
5.  GRI-Mech 3.0 Thermodynamics Database.