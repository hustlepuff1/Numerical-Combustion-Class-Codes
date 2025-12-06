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
    """
    Parser matching the Fortran read sequence, but ignoring '#'-comments.

      1: ni nj
      2: x_left x_right y_bottom y_top
      3: x_ramp theta_deg
      4: gamma Rgas
      5: rho_inf u_inf v_inf
      6: p_inf
      7: t_end CFL
      8: chemistry_on  (1=on, 0=off)
    """
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
    """Return sorted list of all (index, path) for snap_XXXXX.txt."""
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
    return snaps[-1]   # (num, path)


def load_snapshot_ascii(fname_txt, ni, nj):
    """Load U (4,ni,nj) from snap_XXXXX.txt."""
    data = np.loadtxt(fname_txt, usecols=(2, 3, 4, 5))
    if data.shape[0] != ni * nj:
        raise ValueError(f"Unexpected number of rows in snapshot {fname_txt}.")
    U = data.reshape((nj, ni, 4))
    return np.transpose(U, (2, 1, 0))


def load_species_snapshot(fname_txt, ni, nj, nspec=None):
    """
    Load species mass fractions from snapY_XXXXX.txt.

    Fortran format: i, j, Y1, Y2, ..., YN.
    If nspec is None, infer N from number of columns.
    """
    if not os.path.exists(fname_txt):
        return None

    # First read one line to know how many columns total
    with open(fname_txt, "r") as f:
        first = f.readline().strip()
    ncols = len(first.split())
    # 2 columns are i, j
    if nspec is None:
        nspec = ncols - 2

    cols = tuple(range(2, 2 + nspec))
    data = np.loadtxt(fname_txt, usecols=cols)
    if data.shape[0] != ni * nj:
        raise ValueError(f"Unexpected number of rows in species snapshot {fname_txt}.")
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
# 3. Build wedge-mapped grid (like Fortran geometry)
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


# ---------------------------------------------------
# Helper for field levels (avoid constant-field crash)
# ---------------------------------------------------
def safe_levels(field, n=40):
    vmin = float(np.nanmin(field))
    vmax = float(np.nanmax(field))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("Field has NaN/inf values.")
    if vmax <= vmin:
        eps = 1e-12 if abs(vmin) < 1e-6 else 1e-6 * abs(vmin)
        return np.linspace(vmin, vmin + eps, n)
    return np.linspace(vmin, vmax, n)


# ---------------------------------------------------
# 4. Static plots for the latest snapshot
# ---------------------------------------------------
def make_static_plots(params):
    gamma = params["gamma"]; Rgas = params["Rgas"]
    ni = params["ni"]; nj = params["nj"]

    snap_num, snap_path = latest_snapshot()
    print(f"Using latest snapshot: {snap_path}")

    U = load_snapshot_ascii(snap_path, ni, nj)
    rho, u, v, p = cons_to_prim(U, gamma)
    a = np.sqrt(gamma * p / rho)
    Vmag = np.sqrt(u*u + v*v)
    M = Vmag / a
    T = p / (rho * Rgas)

    x, X, Y = build_wedge_grid(params)

    # Load species (if file exists)
    Yspec = None
    spec_path = os.path.join(
        os.path.dirname(snap_path),
        os.path.basename(snap_path).replace("snap_", "snapY_"),
    )
    if os.path.exists(spec_path):
        Yspec = load_species_snapshot(spec_path, ni, nj)
        print(f"Loaded species from {spec_path}")
    else:
        print("Species snapshot not found; skipping Y plots.")

    # ----- Mach contour -----
    plt.figure(figsize=(8, 3))
    levelsM = safe_levels(M, n=40)
    cp = plt.contourf(X, Y, M, levels=levelsM, cmap="jet")
    plt.colorbar(cp, label="Mach")
    plt.plot(x, Y[:, 0], "k-", linewidth=2)
    plt.title("Fortran ODWE - Mach Contour (wedge-mapped grid)")
    plt.xlabel("x"); plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "mach_snapshot.png"), dpi=200)
    plt.close()

    # ----- Temperature contour -----
    plt.figure(figsize=(8, 3))
    levelsT = safe_levels(T, n=40)
    ct = plt.contourf(X, Y, T, levels=levelsT, cmap="jet")
    plt.colorbar(ct, label="Temperature [K]")
    plt.plot(x, Y[:, 0], "k-", linewidth=2)
    plt.title("Fortran ODWE - Temperature Contour (wedge-mapped grid)")
    plt.xlabel("x"); plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "T_snapshot.png"), dpi=200)
    plt.close()

    # ----- Wall pressure -----
    j_wall = 0
    p_wall = p[:, j_wall]
    plt.figure(figsize=(6, 4))
    plt.plot(x, p_wall)
    plt.xlabel("x"); plt.ylabel("Pressure [Pa]")
    plt.title("Fortran ODWE - Wall Pressure")
    plt.grid(True); plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "wall_pressure.png"), dpi=200)
    plt.close()

    # ----- Centerline Mach -----
    j_mid = nj // 2
    M_center = M[:, j_mid]
    plt.figure(figsize=(6, 4))
    plt.plot(x, M_center)
    plt.xlabel("x"); plt.ylabel("Mach")
    plt.title("Fortran ODWE - Centerline Mach (y mid, wedge-mapped)")
    plt.grid(True); plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, "centerline_Mach.png"), dpi=200)
    plt.close()

    # ----- Species plots -----
    if Yspec is not None:
        # index: 0 = H2, 4 = H2O  (matches your Fortran state_mod)
        Y_H2  = Yspec[0]
        Y_H2O = Yspec[4]

        # H2 contour
        plt.figure(figsize=(8, 3))
        levelsY = safe_levels(Y_H2, n=40)
        cY = plt.contourf(X, Y, Y_H2, levels=levelsY, cmap="jet")
        plt.colorbar(cY, label="Y_H2")
        plt.plot(x, Y[:, 0], "k-", linewidth=2)
        plt.title("Fortran ODWE - H2 Mass Fraction")
        plt.xlabel("x"); plt.ylabel("y")
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, "Y_H2_contour.png"), dpi=200)
        plt.close()

        # H2O contour
        plt.figure(figsize=(8, 3))
        levelsY2 = safe_levels(Y_H2O, n=40)
        cY2 = plt.contourf(X, Y, Y_H2O, levels=levelsY2, cmap="jet")
        plt.colorbar(cY2, label="Y_H2O")
        plt.plot(x, Y[:, 0], "k-", linewidth=2)
        plt.title("Fortran ODWE - H2O Mass Fraction")
        plt.xlabel("x"); plt.ylabel("y")
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, "Y_H2O_contour.png"), dpi=200)
        plt.close()

        # H2 along the wall
        Y_H2_wall = Y_H2[:, j_wall]
        plt.figure(figsize=(6, 4))
        plt.plot(x, Y_H2_wall)
        plt.xlabel("x"); plt.ylabel("Y_H2 at wall")
        plt.title("Fortran ODWE - H2 Mass Fraction along Wall")
        plt.grid(True); plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, "Y_H2_wall.png"), dpi=200)
        plt.close()


# ---------------------------------------------------
# 5. Generic animation helper
# ---------------------------------------------------
def save_animation(anim, basename):
    """
    Try saving as MP4 with FFmpeg; if that fails, fall back to GIF.
    """
    mp4_path = os.path.join(RESULTS_DIR, basename + ".mp4")
    gif_path = os.path.join(RESULTS_DIR, basename + ".gif")

    try:
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save(mp4_path, writer=writer)
        print(f"Saved animation: {mp4_path}")
    except Exception as e:
        print(f"MP4 save failed ({e}), falling back to GIF...")
        anim.save(gif_path, writer="pillow", fps=10)
        print(f"Saved animation: {gif_path}")


def make_field_animation(params, field_name, basename, max_frames=None):
    """
    Generic field animation over time using contourf + wedge mapping.

    field_name: "Mach", "Temperature", "Y_H2", "Y_H2O"
    """
    gamma = params["gamma"]; Rgas = params["Rgas"]
    ni = params["ni"]; nj = params["nj"]
    x, X, Y = build_wedge_grid(params)

    snaps = find_snapshot_files()
    if max_frames is not None:
        snaps = snaps[:max_frames]
    nframes = len(snaps)
    print(f"Animating {field_name} over {nframes} frames")

    # Precompute 2D fields for all frames + global min/max
    fields = []
    fmin, fmax = None, None

    for num, path in snaps:
        U = load_snapshot_ascii(path, ni, nj)
        rho, u, v, p = cons_to_prim(U, gamma)
        a = np.sqrt(gamma * p / rho)
        M = np.sqrt(u*u + v*v) / a
        T = p / (rho * Rgas)

        # Initialise as default None
        F = None

        if field_name == "Mach":
            F = M
        elif field_name == "Temperature":
            F = T
        elif field_name in ("Y_H2", "Y_H2O"):
            spec_path = os.path.join(
                os.path.dirname(path),
                os.path.basename(path).replace("snap_", "snapY_"),
            )
            Yspec = load_species_snapshot(spec_path, ni, nj)
            if Yspec is None:
                # If no species snapshot, skip this frame
                continue
            if field_name == "Y_H2":
                F = Yspec[0]   # H2
            else:  # Y_H2O
                F = Yspec[4]   # H2O

        if F is None:
            continue

        fields.append(F)
        fmin_frame = float(np.nanmin(F))
        fmax_frame = float(np.nanmax(F))
        if fmin is None or fmin_frame < fmin:
            fmin = fmin_frame
        if fmax is None or fmax_frame > fmax:
            fmax = fmax_frame

    if not fields:
        print(f"No data to animate for field {field_name}.")
        return

    levels = safe_levels(np.array(fields), n=40)

    fig, ax = plt.subplots(figsize=(8, 3))

    def init():
        ax.clear()
        F0 = fields[0]
        cf = ax.contourf(X, Y, F0, levels=levels, cmap="jet")
        plt.colorbar(cf, ax=ax, label=field_name)
        ax.plot(x, Y[:, 0], "k-", linewidth=2)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"{field_name} (frame 1/{len(fields)})")
        return []

    def update(i):
        ax.clear()
        F = fields[i]
        cf = ax.contourf(X, Y, F, levels=levels, cmap="jet")
        if i == 0:
            plt.colorbar(cf, ax=ax, label=field_name)
        ax.plot(x, Y[:, 0], "k-", linewidth=2)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"{field_name} (frame {i+1}/{len(fields)})")
        return []

    anim = FuncAnimation(fig, update, frames=len(fields), init_func=init,
                         blit=False)
    save_animation(anim, basename)
    plt.close(fig)


# ---------------------------------------------------
# Main
# ---------------------------------------------------
if __name__ == "__main__":
    # Make sure results directory exists
    os.makedirs(RESULTS_DIR, exist_ok=True)

    params = read_input("odwe_input.dat")

    # Static plots (flow + chemistry)
    make_static_plots(params)

    # Animations
    make_field_animation(params, field_name="Mach",        basename="mach_anim",        max_frames=None)
    make_field_animation(params, field_name="Temperature", basename="T_anim",           max_frames=None)
    make_field_animation(params, field_name="Y_H2",        basename="Y_H2_anim",        max_frames=None)
    make_field_animation(params, field_name="Y_H2O",       basename="Y_H2O_anim",       max_frames=None)
