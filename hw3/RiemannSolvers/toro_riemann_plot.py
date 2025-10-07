import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def parse_toro_output(filepath, num_cells, scale_factor=1.0):
    """
    Custom function to parse the specific output format of the Toro FORTRAN code.
    It skips the time header and reads the second column of data.
    """
    try:
        # Use pandas to read the file, skipping the header row
        # The data has an index column and a value column, separated by spaces
        data = pd.read_csv(filepath, header=None, sep='\s+', skiprows=1, names=['index', 'value'])
        # Return the 'value' column, scaled if necessary
        return data['value'].values * scale_factor
    except FileNotFoundError:
        print(f"Error: The data file '{filepath}' was not found.")
        # Return an array of zeros if the file doesn't exist, so the script doesn't crash
        return np.zeros(num_cells)
    except Exception as e:
        print(f"An error occurred while reading '{filepath}': {e}")
        return np.zeros(num_cells)

# --- 1. Define Simulation Parameters ---
# These must match the values used in init.dat
domain_length = 10.0
num_cells = 400
gamma = 1.4
R = 287.0

# --- 2. Create the Spatial Grid ---
# The Toro code uses a cell-centered grid, so we calculate the x-coordinates
# for the center of each of the 400 cells.
dx = domain_length / num_cells
x_coords = np.linspace(dx/2, domain_length - dx/2, num_cells)

# --- 3. Load the Data from Output Files ---
# Update filenames to match the new output.
# Note: Pressure is scaled by 1.0E+6 in the output, so we scale it back.
density_loc = Path(__file__).parent / 'density.dat'
velocity_loc = Path(__file__).parent / 'velocity.dat'
pressure_loc = Path(__file__).parent / 'pressure.dat'


density = parse_toro_output(density_loc, num_cells)
velocity = parse_toro_output(velocity_loc, num_cells)
pressure = parse_toro_output(pressure_loc, num_cells, scale_factor=1.0E+6)

# --- 4. Calculate Derived Quantities ---
# Calculate the other variables needed for the plots
mass_flow = density * velocity
speed_of_sound = np.sqrt(gamma * pressure / density)
# Avoid division by zero for Mach number
mach_number = np.divide(velocity, speed_of_sound, out=np.zeros_like(velocity), where=speed_of_sound!=0)
entropy = (R / (gamma - 1.0)) * np.log(pressure / (density**gamma))

# --- 5. Create the 2x3 Subplot Grid ---
fig, axes = plt.subplots(2, 3, figsize=(22, 12))
fig.suptitle('Open-Source Toro Code Solution at t = 3.9 ms', fontsize=20, y=0.97)

# --- Plot 1: Pressure ---
axes[0, 0].plot(x_coords, pressure, '.-', color='C0')
axes[0, 0].set_title('Pressure', fontsize=14)
axes[0, 0].set_ylabel('Pressure (Pa)', fontsize=12)
axes[0, 0].grid(True)

# --- Plot 2: Entropy ---
axes[0, 1].plot(x_coords, entropy, '.-', color='C1')
axes[0, 1].set_title('Entropy', fontsize=14)
axes[0, 1].set_ylabel('Entropy (s)', fontsize=12)
axes[0, 1].grid(True)

# --- Plot 3: Velocity ---
axes[0, 2].plot(x_coords, velocity, '.-', color='C2')
axes[0, 2].set_title('Velocity', fontsize=14)
axes[0, 2].set_ylabel('Velocity (m/s)', fontsize=12)
axes[0, 2].grid(True)

# --- Plot 4: Mach Number ---
axes[1, 0].plot(x_coords, mach_number, '.-', color='C3')
axes[1, 0].set_title('Mach Number', fontsize=14)
axes[1, 0].set_ylabel('Mach Number', fontsize=12)
axes[1, 0].set_xlabel('Position (x)', fontsize=12)
axes[1, 0].grid(True)

# --- Plot 5: Density ---
axes[1, 1].plot(x_coords, density, '.-', color='C4')
axes[1, 1].set_title('Density', fontsize=14)
axes[1, 1].set_ylabel('Density (kg/m^3)', fontsize=12)
axes[1, 1].set_xlabel('Position (x)', fontsize=12)
axes[1, 1].grid(True)
    
# --- Plot 6: Mass Flow ---
axes[1, 2].plot(x_coords, mass_flow, '.-', color='C5')
axes[1, 2].set_title('Mass Flow', fontsize=14)
axes[1, 2].set_ylabel('Mass Flow (kg/m^2/s)', fontsize=12)
axes[1, 2].set_xlabel('Position (x)', fontsize=12)
axes[1, 2].grid(True)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.show()

