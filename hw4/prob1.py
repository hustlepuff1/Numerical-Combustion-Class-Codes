import pandas as pd
import matplotlib.pyplot as plt

def parse_cea_output(filename):
    """
    Parses a NASA CEA output file to extract Phi and Adiabatic Flame Temp (T).
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data_lines = []
    # Find the start of the data table
    for i, line in enumerate(lines):
        if "PHI" in line and "T(K)" in line:
            data_lines = lines[i+2:] # Skip header lines
            break
            
    phi_vals = []
    temp_vals = []
    
    for line in data_lines:
        if line.strip() == "": # Stop at the end of the table
            break
        parts = line.split()
        if len(parts) > 2:
            phi_vals.append(float(parts[0]))
            temp_vals.append(float(parts[2]))
            
    return pd.Series(temp_vals, index=phi_vals)

# --- Load Data ---
files_q1 = {
    300: 'prob1_300k_output.txt',
    400: 'prob1_400k_output.txt',
    500: 'prob1_500k_output.txt'
}

data_q1 = {T: parse_cea_output(f) for T, f in files_q1.items()}

# --- Plot Data ---
plt.figure(figsize=(12, 8))
colors = ['blue', 'green', 'red']
linestyles = ['-', '--', ':']

# Plot P=1 atm
for i, T_init in enumerate([300, 400, 500]):
    data_p1 = data_q1[T_init].iloc[:14] # First 14 points are for P=1
    plt.plot(data_p1.index, data_p1.values, 
             label=f'P=1 atm, T_init={T_init}K', 
             color=colors[i], linestyle='-')

# Plot P=3 atm
for i, T_init in enumerate([300, 400, 500]):
    data_p3 = data_q1[T_init].iloc[14:28] # Next 14 are for P=3
    plt.plot(data_p3.index, data_p3.values, 
             label=f'P=3 atm, T_init={T_init}K', 
             color=colors[i], linestyle='--')

# Plot P=5 atm
for i, T_init in enumerate([300, 400, 500]):
    data_p5 = data_q1[T_init].iloc[28:] # Last 14 are for P=5
    plt.plot(data_p5.index, data_p5.values, 
             label=f'P=5 atm, T_init={T_init}K', 
             color=colors[i], linestyle=':')

plt.title('Adiabatic Flame Temperature vs. Equivalence Ratio for H2/O2', fontsize=16)
plt.xlabel('Equivalence Ratio ($\phi$)', fontsize=12)
plt.ylabel('Adiabatic Flame Temperature (K)', fontsize=12)
plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()