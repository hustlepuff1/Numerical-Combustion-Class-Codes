import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# ---------------------------------------------------
# Global config
# ---------------------------------------------------
RESULTS_DIR = "results"

# ---------------------------------------------------
# 1. Read simple Fortran input file
# ---------------------------------------------------
def read_input(fname="odwe_input.dat"):
    if not os.path.exists(fname):
        raise FileNotFoundError(fname)

    vals = []
    with open(fname, "r") as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            vals.extend(line.split())

    it = iter(vals)
    ni   = int(next(it)); nj = int(next(it))
    x_left = float(next(it)); x_right = float(next(it))
    y_bottom = float(next(it)); y_top = float(next(it))
    x_ramp = float(next(it)); theta_deg = float(next(it))
    gamma = float(next(it)); Rgas = float(next(it))
    rho_inf = float(next(it)); u_inf = float(next(it)); v_inf = float(next(it))
    p_inf = float(next(it))
    t_end = float(next(it)); CFL = float(next(it))
    chemistry_on = int(next(it))

    return {
        "ni": ni, "nj": nj,
        "x_left": x_left, "x_right": x_right,
        "y_bottom": y_bottom, "y_top": y_top,
        "x_ramp": x_ramp, "theta_deg": theta_deg,
        "gamma": gamma, "Rgas": Rgas,
        "rho_inf": rho_inf, "u_inf": u_inf, "v_inf": v_inf, "p_inf": p_inf,
        "t_end": t_end, "CFL": CFL,
        "chemistry_on": chemistry_on,
    }

# ---------------------------------------------------
# 2. Snapshot utilities
# ---------------------------------------------------
def find_snapshot_files():
    snaps = []
    for root in [".", "build"]:
        if not os.path.isdir(root):
            continue
        for fname in os.listdir(root):
            if fname.startswith("snap_") and fname.endswith(".txt"):
                core = fname[len("snap_"):-4]
                try:
                    num = int(core)
                except ValueError:
                    continue
                snaps.append((num, os.path.join(root, fname)))
    if not snaps:
        raise FileNotFoundError("No snapshot files found.")
    snaps.sort(key=lambda x: x[0])
    return snaps

def latest_snapshot():
    snaps = find_snapshot_files()
    return snaps[-1]

def load_snapshot_ascii(fname_txt, ni, nj):
    data = np.loadtxt(fname_txt, usecols=(2, 3, 4, 5))
    if data.shape[0] != ni * nj:
        print(f"Warning: {fname_txt} has {data.shape[0]} rows, expected {ni*nj}. Skipping.")
        return None
    U = data.reshape((nj, ni, 4))
    return np.transpose(U, (2, 1, 0))

def load_species_snapshot(fname_txt, ni, nj, nspec=None):
    if not os.path.exists(fname_txt):
        return None
    # Quick check of file length
    with open(fname_txt, 'r') as f:
        lines = f.readlines()
        if len(lines) != ni * nj:
             return None

    with open(fname_txt, "r") as f:
        first = f.readline().strip()
    ncols = len(first.split())
    if nspec is None:
        nspec = ncols - 2
    cols = tuple(range(2, 2 + nspec))
    
    try:
        data = np.loadtxt(fname_txt, usecols=cols)
    except:
        return None

    if data.shape[0] != ni * nj:
        return None
    Y = data.reshape((nj, ni, nspec))
    return np.transpose(Y, (2, 1, 0))

def cons_to_prim(U, gamma):
    rho = U[0]
    u   = U[1] / rho
    v   = U[2] / rho
    E   = U[3] / rho
    p   = (gamma - 1.0) * rho * (E - 0.5 * (u*u + v*v))
    return rho, u, v, p

# ---------------------------------------------------
# 3. Build wedge-mapped grid
# ---------------------------------------------------
def build_wedge_grid(params):
    ni = params["ni"]; nj = params["nj"]
    x_left = params["x_left"]; x_right = params["x_right"]
    y_bottom = params["y_bottom"]; y_top = params["y_top"]
    x_ramp = params["x_ramp"]; theta = np.deg2rad(params["theta_deg"])

    x = np.linspace(x_left, x_right, ni)
    y_comp = np.linspace(y_bottom, y_top, nj)
    Xc, Yc = np.meshgrid(x, y_comp, indexing="ij")

    Y_phys = np.empty_like(Yc)
    mask_flat = Xc < x_ramp
    Y_phys[mask_flat] = Yc[mask_flat]

    mask_ramp = ~mask_flat
    x_r = Xc[mask_ramp]
    y_comp_r = Yc[mask_ramp]
    y_ramp = (x_r - x_ramp) * np.tan(theta)
    Y_phys[mask_ramp] = y_ramp + (y_top - y_ramp) * (y_comp_r / y_top)

    return x, Xc, Y_phys

# --- IMPROVED SAFE_LEVELS ---
def safe_levels(field, n=40):
    """
    Generate contour levels safely, even for constant fields.
    """
    # Filter out NaNs for min/max calculation
    valid_data = field[np.isfinite(field)]
    if len(valid_data) == 0:
        return np.linspace(0, 1, n)

    vmin = float(np.min(valid_data))
    vmax = float(np.max(valid_data))

    # If the field is constant (or nearly constant)
    if vmax <= vmin + 1e-14:
        # Create a small artificial range centered on the value
        offset = abs(vmin) * 0.1 if abs(vmin) > 1e-14 else 0.1
        return np.linspace(vmin - offset, vmin + offset, n)
    
    return np.linspace(vmin, vmax, n)

# ---------------------------------------------------
# 4. PLOTTING FUNCTIONS
# ---------------------------------------------------

def make_grid_plot(params):
    print("Generating Grid Mesh Plot...")
    ni = params["ni"]; nj = params["nj"]
    _, X, Y = build_wedge_grid(params)

    plt.figure(figsize=(10, 5))
    plt.plot(X, Y, 'k-', linewidth=0.3, alpha=0.5)
    plt.plot(X.T, Y.T, 'k-', linewidth=0.3, alpha=0.5)
    plt.title(f"Computational Mesh ({ni} x {nj})")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.axis("equal")
    plt.tight_layout()
    outpath = os.path.join(RESULTS_DIR, "grid_mesh.png")
    plt.savefig(outpath, dpi=300)
    print(f"Saved: {outpath}")
    plt.close()

def make_static_plots(params):
    gamma = params["gamma"]; Rgas = params["Rgas"]
    ni = params["ni"]; nj = params["nj"]

    snap_num, snap_path = latest_snapshot()
    print(f"Processing static plots for snapshot: {snap_path}")

    U = load_snapshot_ascii(snap_path, ni, nj)
    if U is None: return

    rho, u, v, p = cons_to_prim(U, gamma)
    a = np.sqrt(gamma * p / rho)
    Vmag = np.sqrt(u*u + v*v)
    M = Vmag / a
    T = p / (rho * Rgas)

    _, X, Y = build_wedge_grid(params)

    Yspec = None
    spec_path = os.path.join(
        os.path.dirname(snap_path),
        os.path.basename(snap_path).replace("snap_", "snapY_"),
    )
    if os.path.exists(spec_path):
        Yspec = load_species_snapshot(spec_path, ni, nj)
        print(f"Loaded species data from {spec_path}")

    # Plot Mach
    plt.figure(figsize=(8, 3))
    levelsM = safe_levels(M, n=40)
    cp = plt.contourf(X, Y, M, levels=levelsM, cmap="jet")
    plt.colorbar(cp, label="Mach")
    plt.plot(X[:, 0], Y[:, 0], "k-", linewidth=2)
    plt.title("Mach Number Contour")
    plt.xlabel("x"); plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "mach_snapshot.png"), dpi=200)
    plt.close()

    # Plot Temperature
    plt.figure(figsize=(8, 3))
    levelsT = safe_levels(T, n=40)
    ct = plt.contourf(X, Y, T, levels=levelsT, cmap="inferno")
    plt.colorbar(ct, label="Temperature [K]")
    plt.plot(X[:, 0], Y[:, 0], "k-", linewidth=2)
    plt.title("Temperature Contour")
    plt.xlabel("x"); plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "T_snapshot.png"), dpi=200)
    plt.close()

    if Yspec is not None:
        Y_H2O = Yspec[4]
        plt.figure(figsize=(8, 3))
        levelsY = safe_levels(Y_H2O, n=40)
        cY = plt.contourf(X, Y, Y_H2O, levels=levelsY, cmap="gist_heat_r") 
        plt.colorbar(cY, label="Mass Fraction H2O")
        plt.plot(X[:, 0], Y[:, 0], "k-", linewidth=2)
        plt.title("H2O Mass Fraction")
        plt.xlabel("x"); plt.ylabel("y")
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, "Y_H2O_contour.png"), dpi=200)
        plt.close()

def save_animation(anim, basename):
    mp4_path = os.path.join(RESULTS_DIR, basename + ".mp4")
    gif_path = os.path.join(RESULTS_DIR, basename + ".gif")
    try:
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save(mp4_path, writer=writer)
        print(f"Saved: {mp4_path}")
    except Exception as e:
        print(f"FFmpeg failed, saving GIF: {gif_path}")
        try:
            anim.save(gif_path, writer="pillow", fps=10)
            print(f"Saved: {gif_path}")
        except Exception as e2:
            print(f"GIF save failed too: {e2}")

def make_field_animation(params, field_name, basename, max_frames=None):
    gamma = params["gamma"]; Rgas = params["Rgas"]
    ni = params["ni"]; nj = params["nj"]
    _, X, Y = build_wedge_grid(params)

    snaps = find_snapshot_files()
    if max_frames: snaps = snaps[:max_frames]
    
    print(f"Animating {field_name}...")
    fields = []
    
    for num, path in snaps:
        U = load_snapshot_ascii(path, ni, nj)
        if U is None: continue
        
        rho, u, v, p = cons_to_prim(U, gamma)
        
        if field_name == "Mach":
            a = np.sqrt(gamma * p / rho)
            F = np.sqrt(u*u+v*v)/a
        elif field_name == "Temperature":
            F = p / (rho * Rgas)
        elif field_name == "Y_H2O":
            spath = path.replace("snap_", "snapY_")
            Yspec = load_species_snapshot(spath, ni, nj)
            if Yspec is None: continue
            F = Yspec[4] # H2O
        elif field_name == "Y_H2":
            spath = path.replace("snap_", "snapY_")
            Yspec = load_species_snapshot(spath, ni, nj)
            if Yspec is None: continue
            F = Yspec[0] # H2
        
        fields.append(F)

    if not fields: 
        print(f"No valid data for {field_name}")
        return

    # Use the improved safe_levels to avoid crashes on flat data
    levels = safe_levels(np.array(fields), n=40)
    
    fig, ax = plt.subplots(figsize=(8, 3))

    def init():
        ax.clear()
        return []

    def update(i):
        ax.clear()
        cf = ax.contourf(X, Y, fields[i], levels=levels, cmap="jet" if "Y_H2O" not in field_name else "gist_heat_r")
        ax.plot(X[:, 0], Y[:, 0], "k-", linewidth=2)
        ax.set_title(f"{field_name} - Frame {i}")
        return []

    anim = FuncAnimation(fig, update, frames=len(fields), init_func=init)
    save_animation(anim, basename)
    plt.close(fig)

if __name__ == "__main__":
    os.makedirs(RESULTS_DIR, exist_ok=True)
    params = read_input("odwe_input.dat")

    make_grid_plot(params)
    make_static_plots(params)
    make_field_animation(params, "Mach", "mach_anim")
    make_field_animation(params, "Temperature", "T_anim")
    make_field_animation(params, "Y_H2", "Y_H2_anim")
    make_field_animation(params, "Y_H2O", "Y_H2O_anim")