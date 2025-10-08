# ====================================================================
# Program: FVM Roe Solver Failure Visualizer
# Purpose: Loads the data from the time-marching FVM solver that
#          failed and plots the results to show where the instability
#          (NaN values) occurred.
# ====================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Load the Data File ---
# We explicitly define the column names to ensure robust parsing.
# The FVM code output has columns: x, pressure, density, velocity
column_names = ['x', 'pressure', 'density', 'velocity']
file_to_plot = 'results_roe.dat'

try:
    # Load the data, skipping the header row and using our defined names
    df = pd.read_csv(file_to_plot, sep='\s+', header=None, names=column_names, skiprows=1)
    print(f"Successfully loaded '{file_to_plot}'")

    # --- 2. Calculate Derived Quantities ---
    gamma = 1.4
    R = 287.0
    
    df['mass_flow'] = df['density'] * df['velocity']
    df['speed_of_sound'] = np.sqrt(gamma * df['pressure'] / df['density'])
    df['mach_number'] = np.divide(df['velocity'], df['speed_of_sound'], out=np.zeros_like(df['velocity']), where=df['speed_of_sound']!=0)
    df['entropy'] = (R / (gamma - 1.0)) * np.log(df['pressure'] / (df['density']**gamma))

    # --- 3. Create the 2x3 Subplot Grid ---
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle('FVM Roe Solver Results Showing Instability in Expansion Fan (t = 3.9 ms)', fontsize=20, y=0.97)

    # --- Plot 1: Pressure ---
    axes[0, 0].plot(df['x'], df['pressure'], '.-', color='C0')
    axes[0, 0].set_title('Pressure', fontsize=14)
    axes[0, 0].set_ylabel('Pressure (Pa)', fontsize=12)
    axes[0, 0].grid(True)

    # --- Plot 2: Entropy ---
    axes[0, 1].plot(df['x'], df['entropy'], '.-', color='C1')
    axes[0, 1].set_title('Entropy', fontsize=14)
    axes[0, 1].set_ylabel('Entropy (s)', fontsize=12)
    axes[0, 1].grid(True)

    # --- Plot 3: Velocity ---
    axes[0, 2].plot(df['x'], df['velocity'], '.-', color='C2')
    axes[0, 2].set_title('Velocity', fontsize=14)
    axes[0, 2].set_ylabel('Velocity (m/s)', fontsize=12)
    axes[0, 2].grid(True)

    # --- Plot 4: Mach Number ---
    axes[1, 0].plot(df['x'], df['mach_number'], '.-', color='C3')
    axes[1, 0].set_title('Mach Number', fontsize=14)
    axes[1, 0].set_ylabel('Mach Number', fontsize=12)
    axes[1, 0].set_xlabel('Position (x)', fontsize=12)
    axes[1, 0].grid(True)

    # --- Plot 5: Density ---
    axes[1, 1].plot(df['x'], df['density'], '.-', color='C4')
    axes[1, 1].set_title('Density', fontsize=14)
    axes[1, 1].set_ylabel('Density (kg/m^3)', fontsize=12)
    axes[1, 1].set_xlabel('Position (x)', fontsize=12)
    axes[1, 1].grid(True)
    
    # --- Plot 6: Mass Flow ---
    axes[1, 2].plot(df['x'], df['mass_flow'], '.-', color='C5')
    axes[1, 2].set_title('Mass Flow', fontsize=14)
    axes[1, 2].set_ylabel('Mass Flow (kg/m^2/s)', fontsize=12)
    axes[1, 2].set_xlabel('Position (x)', fontsize=12)
    axes[1, 2].grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.show()

except FileNotFoundError:
    print(f"Error: Data file not found at '{file_to_plot}'")
except Exception as e:
    print(f"An error occurred: {e}")
