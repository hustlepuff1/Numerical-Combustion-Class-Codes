import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_simulation():
    output_file = 'results.csv'
    
    print(f"[Plot] Reading {output_file}...")
    
    try:
        # skipinitialspace=True handles headers written like "Time, Temp"
        df = pd.read_csv(output_file, skipinitialspace=True)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # 1. Clean Column Names
    df.columns = df.columns.str.strip()

    # 2. Check for 'Time'
    if 'Time' not in df.columns:
        for col in df.columns:
            if col.lower().strip() == 'time':
                df.rename(columns={col: 'Time'}, inplace=True)
                break
        else:
            print(f"Error: 'Time' column not found. Found: {df.columns.tolist()}")
            sys.exit(1)

    # 3. Check for Divergence (NaNs)
    if df.isnull().values.any():
        print("\n[WARNING] NaNs detected in the data! The simulation likely diverged.")
        nan_rows = df[df.isna().any(axis=1)]
        first_fail_time = nan_rows['Time'].iloc[0] if not nan_rows.empty else 'Unknown'
        print(f"Simulation failed at Time = {first_fail_time}")
        df = df.dropna()

    if len(df) < 2:
        print("\n[ERROR] Not enough valid data points to plot.")
        return

    # 4. Plotting - UPDATED to 2 Subplots
    # nrows=2, ncols=1, sharex=True links the time axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # --- Plot 1: Temperature ---
    if 'Temp' in df.columns:
        color = 'tab:red'
        ax1.plot(df['Time'], df['Temp'], color=color, linewidth=2, label='Temperature')
        ax1.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.legend(loc='best')
        ax1.set_title('Combustion Simulation Results', fontsize=14)
    else:
        ax1.text(0.5, 0.5, "Temperature Data Missing", ha='center')

    # --- Plot 2: Species ---
    species_cols = [col for col in df.columns if col.startswith('Y_')]
    
    if species_cols:
        # Loop through species columns (e.g., Y_A, Y_B, Y_fuel)
        for i, species in enumerate(species_cols):
            # Cycle through different line styles and colors automatically
            ax2.plot(df['Time'], df[species], linewidth=2, label=species)
            
        ax2.set_ylabel('Mass Fraction', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend(loc='best')
    else:
        ax2.text(0.5, 0.5, "Species Data Missing", ha='center')

    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    save_name = 'results_plot.png'
    plt.savefig(save_name, dpi=300)
    print(f"[Plot] Graph saved to {save_name}")

if __name__ == "__main__":
    plot_simulation()