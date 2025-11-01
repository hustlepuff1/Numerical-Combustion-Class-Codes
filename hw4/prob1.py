# ---------------------------------------------------
# Plot Adiabatic Flame Temperature (H2/O2, HP case)
# ---------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd

# Assume 'prob1_data' DataFrame already exists or load from CSV
# Columns: T0_K, phi, P_atm, Tad_K
# If you already saved it:
prob1_data = pd.read_csv("prob1_H2O2_Tad_table.csv")

# Plot: one chart per initial temperature
for T0 in [300, 400, 500]:
    dfT = prob1_data[prob1_data["T0_K"] == T0].sort_values(["phi", "P_atm"])
    plt.figure(figsize=(7, 5))
    for P in [1, 3, 5]:
        sub = dfT[dfT["P_atm"] == P]
        plt.plot(sub["phi"], sub["Tad_K"], marker='o', label=f"{P} atm")

    plt.xlabel("Equivalence Ratio (ϕ)")
    plt.ylabel("Adiabatic Flame Temperature (K)")
    plt.title(f"H₂/O₂ Adiabatic Flame Temperature (HP) — T₀ = {T0} K")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"Tad_vs_phi_T0_{T0}K.png", dpi=200)
    plt.show()
