import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- parameters ---
phi_values = np.arange(0.7, 1.21, 0.1)
P = 5 * ct.one_atm
T_in = 400  # K
gas = ct.Solution("gri30.yaml")

flame_speeds = []
prev_flame = None  # for continuation (speeds up successive runs)

for phi in phi_values:
    gas.set_equivalence_ratio(phi, "C3H8", "O2:1.0, N2:3.76")
    gas.TP = T_in, P

    flame = ct.FreeFlame(gas, width=0.03)  # narrower domain
    flame.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
    flame.transport_model = "mixture-averaged"  # updated keyword

    if prev_flame is not None:
        flame.set_initial_guess(data=prev_flame.to_array())  # reuse previous solution

    flame.solve(loglevel=1, auto=True)
    prev_flame = flame

    flame_speeds.append(flame.velocity[0])
    print(f"phi={phi:.2f}, S_L={flame.velocity[0]:.3f} m/s")

# --- save & plot ---
df = pd.DataFrame({"phi": phi_values, "S_L_m_per_s": flame_speeds})
df.to_csv("prob4_flame_speed.csv", index=False)

plt.figure(figsize=(7,5))
plt.plot(df["phi"], df["S_L_m_per_s"], "o-", lw=2)
plt.xlabel("Equivalence Ratio (ϕ)")
plt.ylabel("Laminar Flame Speed [m/s]")
plt.title("Propane–Air Laminar Flame Speed (P=5 atm, T₀=400 K)")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("prob4_propane_flame_speed.png", dpi=200)
plt.show()
