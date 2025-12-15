import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import matplotlib.animation as animation
from pathlib import Path

# --- Configuration ---
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / 'results_fortran'
OUTPUT_DIR = BASE_DIR / 'results_python'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Variable Mapping ---
VAR_MAP = {
    'Rho': 2, 'U': 3, 'V': 4, 'P': 5, 'T': 6,
    'Y_H2': 7, 'Y_O2': 8, 'Y_H2O': 9, 
    'Y_OH': 10, 'Y_O': 11, 'Y_H': 12, 'Mach': 13
}

def load_dat(filepath):
    try:
        raw = np.loadtxt(filepath, skiprows=1)
        x_col = raw[:, 0]
        # Auto-detect grid dimensions
        jumps = np.where(np.diff(x_col) < 0)[0]
        Nx = jumps[0] + 1 if len(jumps) > 0 else 201
        Ny = len(raw) // Nx
        return raw.reshape((Ny, Nx, raw.shape[1]))
    except Exception as e:
        return None

# --- 0. Grid Plotter (Added) ---
def plot_grid(data, title, filename):
    X = data[:, :, 0]
    Y = data[:, :, 1]
    
    fig, ax = plt.subplots(figsize=(10, 3))
    
    # Plot vertical lines (columns)
    ax.plot(X.T, Y.T, 'k-', linewidth=0.5, alpha=0.5)
    # Plot horizontal lines (rows)
    ax.plot(X, Y, 'k-', linewidth=0.5, alpha=0.5)
    
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect('equal')
    
    plt.savefig(OUTPUT_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated: {filename}")

# --- 1. Contour Plotter ---
def plot_contour(data, var_name, title, filename, cmap='jet'):
    if var_name not in VAR_MAP: return
    idx = VAR_MAP[var_name]
    
    Z = data[:, :, idx]
    X = data[:, :, 0]
    Y = data[:, :, 1]
    
    fig, ax = plt.subplots(figsize=(10, 3))
    cf = ax.contourf(X, Y, Z, levels=50, cmap=cmap)
    cbar = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(var_name)
    
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect('equal')
    
    plt.savefig(OUTPUT_DIR / filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Generated: {filename}")

# --- 2. Centerline Plotter ---
def plot_centerline(data, var_name, title, filename):
    if var_name not in VAR_MAP: return
    idx = VAR_MAP[var_name]
    ny = data.shape[0]
    mid_j = ny // 2
    
    x = data[mid_j, :, 0]
    y_val = data[mid_j, :, idx]
    
    plt.figure(figsize=(8, 4))
    plt.plot(x, y_val, 'b-', linewidth=2, label='Simulation')
    plt.title(title)
    plt.xlabel("Distance (m)")
    plt.ylabel(var_name)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.savefig(OUTPUT_DIR / filename, dpi=200)
    plt.close()
    print(f"Generated: {filename}")

# --- 3. NASA Validation Plot (Fig 6 Style) ---
def plot_centerline_validation(data, filename):
    ny = data.shape[0]
    mid_j = ny // 2
    x = data[mid_j, :, 0]
    
    y_h2  = data[mid_j, :, VAR_MAP['Y_H2']]
    y_h2o = data[mid_j, :, VAR_MAP['Y_H2O']]
    y_o2  = data[mid_j, :, VAR_MAP['Y_O2']]
    
    fig, ax1 = plt.subplots(figsize=(9, 5))
    
    # Left Axis (H2, H2O)
    line1, = ax1.plot(x, y_h2, 'g--', linewidth=2, label='H2 (Left)')
    line2, = ax1.plot(x, y_h2o, 'r-', linewidth=2, label='H2O (Left)')
    ax1.set_xlabel("Distance (m)")
    ax1.set_ylabel("Mass Fraction (H2, H2O)", color='k')
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Right Axis (O2)
    ax2 = ax1.twinx()
    line3, = ax2.plot(x, y_o2, 'b-.', linewidth=2, label='O2 (Right)')
    ax2.set_ylabel("Mass Fraction (O2)", color='b')
    ax2.tick_params(axis='y', labelcolor='b')
    
    lines = [line1, line2, line3]
    labs = [l.get_label() for l in lines]
    ax1.legend(lines, labs, loc='center right')
    
    plt.title("Centerline Species Distribution (NASA Fig 6 Comparison)")
    plt.savefig(OUTPUT_DIR / filename, dpi=200)
    plt.close()
    print(f"Generated: {filename}")

# --- 4. Convergence History ---
def plot_convergence():
    hist_file = DATA_DIR / 'history.dat'
    if not hist_file.exists(): 
        print("Warning: history.dat not found.")
        return

    try:
        try: hist = np.loadtxt(hist_file, skiprows=1, delimiter=',')
        except: hist = np.loadtxt(hist_file, skiprows=1)

        if hist.ndim < 2: return
        
        plt.figure(figsize=(8, 5))
        plt.semilogy(hist[:, 0], hist[:, 2], 'r-', linewidth=1.5)
        plt.title("Convergence History (L2 Norm)")
        plt.xlabel("Iterations")
        plt.ylabel("L2 Norm (Rho)")
        plt.grid(True, which="both", ls="-", alpha=0.4)
        plt.savefig(OUTPUT_DIR / 'fig2_convergence.png', dpi=200)
        plt.close()
        print("Generated: fig2_convergence.png")
    except Exception as e:
        print(f"Failed to plot convergence: {e}")

# --- 5. Movie Maker (No Colorbar, Fixed Limits) ---
def get_global_limits(file_list, var_name):
    idx = VAR_MAP[var_name]
    g_min, g_max = 1e20, -1e20
    # Scan every 5th frame for speed
    for f in file_list[::5]: 
        d = load_dat(f)
        if d is not None:
            val = d[:, :, idx]
            g_min = min(g_min, val.min())
            g_max = max(g_max, val.max())
    if g_max <= g_min: g_max = g_min + 1e-9
    return g_min, g_max

def make_movie_fixed(file_list, var_name, filename, cmap='jet'):
    if not file_list: return
    print(f"Rendering {filename}...")
    
    vmin, vmax = get_global_limits(file_list, var_name)
    data0 = load_dat(file_list[0])
    if data0 is None: return

    fig, ax = plt.subplots(figsize=(10, 3))
    idx = VAR_MAP[var_name]
    
    # Plot without colorbar
    X = data0[:, :, 0]
    Y = data0[:, :, 1]
    Z = data0[:, :, idx]
    ax.contourf(X, Y, Z, levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_aspect('equal')
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    
    def update(frame_file):
        ax.clear()
        d = load_dat(frame_file)
        if d is not None:
            ax.contourf(X, Y, d[:, :, idx], levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_title(f"{var_name} - {os.path.basename(frame_file)}")
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            ax.set_aspect('equal')

    ani = animation.FuncAnimation(fig, update, frames=file_list, interval=100)
    ani.save(OUTPUT_DIR / filename, writer='ffmpeg', dpi=150)
    plt.close()

# ====== MAIN EXECUTION ======
files = sorted(glob.glob(str(DATA_DIR / 'flow_*.dat')))
if not files:
    print("No data found.")
    exit()

# --- A. STATIC PLOTS (Selectable Frame) ---
target_file = files[-1] # Default: Last Frame
if len(sys.argv) > 1:
    arg = sys.argv[1]
    # Check if integer index
    try:
        idx = int(arg)
        potential = DATA_DIR / f"flow_{idx:03d}.dat"
        if potential.exists():
            target_file = str(potential)
        else:
            print(f"Warning: Frame {idx} not found. Using last frame.")
    except ValueError:
        # Check if direct filename
        potential = DATA_DIR / arg
        if potential.exists():
            target_file = str(potential)
        else:
            print(f"Warning: File {arg} not found. Using last frame.")

print(f"--- Plotting Static Validation (Frame: {os.path.basename(str(target_file))}) ---")
data = load_dat(target_file)

if data is not None:
    # 0. Grid Plot
    plot_grid(data, 'Computational Grid (201 x 15)', 'fig1_grid.png')

    # 1. Contours (Figs 9-17)
    plot_contour(data, 'Mach', 'Mach Number (Fig 9)', 'fig9_mach.png')
    plot_contour(data, 'P', 'Static Pressure (Fig 10)', 'fig10_pressure.png')
    plot_contour(data, 'T', 'Static Temperature (Fig 11)', 'fig11_temp.png')
    plot_contour(data, 'Y_H2', 'H2 Mass Fraction (Fig 12)', 'fig12_h2.png', cmap='inferno')
    plot_contour(data, 'Y_O2', 'O2 Mass Fraction (Fig 13)', 'fig13_o2.png')
    plot_contour(data, 'Y_H2O', 'H2O Mass Fraction (Fig 14)', 'fig14_h2o.png', cmap='viridis')
    plot_contour(data, 'Y_OH', 'OH Mass Fraction (Fig 15)', 'fig15_oh.png', cmap='magma')
    plot_contour(data, 'Y_O', 'O Mass Fraction (Fig 16)', 'fig16_o.png', cmap='plasma')
    plot_contour(data, 'Y_H', 'H Mass Fraction (Fig 17)', 'fig17_h.png', cmap='plasma')

    # 2. Centerlines (Figs 4, 5, 7, 8)
    plot_centerline(data, 'T', 'Centerline Temperature (Fig 7)', 'fig7_center_temp.png')
    plot_centerline(data, 'P', 'Centerline Pressure (Fig 8)', 'fig8_center_p.png')
    plot_centerline(data, 'Y_H2O', 'Centerline H2O (Fig 4)', 'fig4_center_h2o.png')
    plot_centerline(data, 'Y_O2', 'Centerline O2 (Fig 5)', 'fig5_center_o2.png')
    
    # 3. Special Validation (Fig 6 Comparison)
    plot_centerline_validation(data, 'fig6_validation.png')

# 4. Convergence History
plot_convergence()

# 5. Movies (Always Full History)
print("--- Generating Unsteady Animations (Full History) ---")
movie_files = files 

# Physics
make_movie_fixed(movie_files, 'T',    'movie_T.mp4',    cmap='jet')
make_movie_fixed(movie_files, 'P',    'movie_P.mp4',    cmap='jet')
make_movie_fixed(movie_files, 'Mach', 'movie_Mach.mp4', cmap='jet')
make_movie_fixed(movie_files, 'Rho',  'movie_Rho.mp4',  cmap='jet')

# Species
make_movie_fixed(movie_files, 'Y_H2',  'movie_species_H2.mp4',  cmap='inferno')
make_movie_fixed(movie_files, 'Y_O2',  'movie_species_O2.mp4',  cmap='Blues')
make_movie_fixed(movie_files, 'Y_H2O', 'movie_species_H2O.mp4', cmap='viridis')
make_movie_fixed(movie_files, 'Y_OH',  'movie_species_OH.mp4',  cmap='magma')
make_movie_fixed(movie_files, 'Y_O',   'movie_species_O.mp4',   cmap='plasma')
make_movie_fixed(movie_files, 'Y_H',   'movie_species_H.mp4',   cmap='plasma')