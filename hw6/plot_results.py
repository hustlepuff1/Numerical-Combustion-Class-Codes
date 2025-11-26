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

    # 3. Check for Divergence
    if df.isnull().values.any():
        print("\n[WARNING] NaNs detected in the data!")
        df = df.dropna()

    if len(df) < 2:
        print("\n[ERROR] Not enough valid data points to plot.")
        return

    # ==========================================
    # FIGURE 1: LINEAR SCALE
    # ==========================================
    print("[Plot] Generating Linear Plots...")
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Temperature
    if 'Temp' in df.columns:
        color = 'tab:red'
        ax1.plot(df['Time'], df['Temp'], color=color, linewidth=2, label='Temperature')
        ax1.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.legend(loc='best')
        ax1.set_title('Combustion Results (Linear Scale)', fontsize=14)

    # Species
    species_cols = [col for col in df.columns if col.startswith('Y_')]
    if species_cols:
        for species in species_cols:
            ax2.plot(df['Time'], df[species], linewidth=2, label=species)
        ax2.set_ylabel('Mass Fraction', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend(loc='center right')

    plt.tight_layout()
    plt.savefig('results_plot.png', dpi=300)
    plt.close(fig1)

    # ==========================================
    # FIGURE 2: LOG SCALE (As requested)
    # ==========================================
    if species_cols:
        print("[Plot] Generating Log-Time Scale Plot...")
        fig2, ax3 = plt.subplots(figsize=(10, 6))
        
        for species in species_cols:
            ax3.plot(df['Time'], df[species], linewidth=2, label=species)
            
        # LOG MAGIC HERE
        ax3.set_xscale('log') 
        
        ax3.set_ylabel('Mass Fraction', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Time (s) - Log Scale', fontsize=12, fontweight='bold')
        ax3.set_title('Species Mass Fraction (Log Time Scale)', fontsize=14)
        ax3.grid(True, which="both", linestyle='--', alpha=0.5)
        ax3.legend(loc='best')
        
        plt.tight_layout()
        plt.savefig('results_plot_log.png', dpi=300)
        plt.close(fig2)
        print(f"[Plot] Saved 'results_plot.png' and 'results_plot_log.png'")

if __name__ == "__main__":
    plot_simulation()