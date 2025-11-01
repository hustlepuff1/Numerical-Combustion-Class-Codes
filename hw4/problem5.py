# ------------------------------------------------------------
# Problem 5: DME/Air ignition delay (CV vs CP) using Cantera
# Requires: dme.yaml (converted from LLNL CHEMKIN)
# ------------------------------------------------------------
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MECH_FILE = "dme.yaml"     # produced by ck2yaml
P0 = 1.0 * ct.one_atm # change to 5 if NaN at low T
phi = 1.0
T_min, T_max, dT = 800.0, 1400.0, 50.0
# Allow more time at low T; you can raise this if 800–900 K doesn’t ignite
max_time_map = [(0, 0.5), (900, 0.2), (1100, 0.05), (1300, 0.02)]  # (T_threshold, max_time s)
rtol, atol = 1e-9, 1e-15

gas = ct.Solution(MECH_FILE)

# Detect DME species name
fuel_candidates = ["ch3och3", "CH3OCH3", "DME"]
fuel = next((s for s in fuel_candidates if s in gas.species_names), None)
if fuel is None:
    raise ValueError(f"Could not find a DME species among {fuel_candidates} in {MECH_FILE}")

air = "O2:1.0, N2:3.76"

def set_phi_TP(g, phi_val, T0, P0):
    g.set_equivalence_ratio(phi_val, fuel=fuel, oxidizer=air)
    g.TP = T0, P0

def pick_max_time(T0):
    mt = max_time_map[0][1]
    for thr, val in max_time_map:
        if T0 >= thr: mt = val
    return mt

def ignition_delay(times, temps):
    dTdt = np.gradient(temps, times, edge_order=2)
    k = np.argmax(dTdt)
    return times[k], k

def run_one(T0, is_cp):
    g = ct.Solution(MECH_FILE)
    set_phi_TP(g, phi, T0, P0)
    r = ct.IdealGasConstPressureReactor(g) if is_cp else ct.IdealGasReactor(g)
    sim = ct.ReactorNet([r])
    sim.rtol = rtol; sim.atol = atol

    t_final = pick_max_time(T0)
    t = 0.0
    times, temps = [], []
    while t < t_final:
        t = sim.step()
        times.append(t); temps.append(r.T)
        # early stop: large temp jump → take a few extra steps and stop
        if len(times) > 5 and temps[-1] - temps[0] > 1000:
            for _ in range(5):
                t = sim.step(); times.append(t); temps.append(r.T)
            break

    times = np.array(times); temps = np.array(temps)
    if len(times) < 5 or (temps.max() - temps[0]) < 10:
        return np.nan
    tau, _ = ignition_delay(times, temps)
    return tau

T_list = np.arange(T_min, T_max + 1e-9, dT)
rows = []
print(f"Mechanism: {MECH_FILE}, fuel='{fuel}'")
for T0 in T_list:
    tau_cv = run_one(T0, is_cp=False)
    tau_cp = run_one(T0, is_cp=True)
    rows.append((T0, 1.0/T0, tau_cv, tau_cp))
    print(f"T0={T0:.1f} K  ->  tau_CV={tau_cv:.6g} s,  tau_CP={tau_cp:.6g} s")

df = pd.DataFrame(rows, columns=["T0_K", "invT_1perK", "tau_CV_s", "tau_CP_s"])
df.to_csv("prob5_DME_ignition_delays.csv", index=False)
print("Saved: prob5_DME_ignition_delays.csv")

plt.figure(figsize=(7,5))
plt.semilogy(df["invT_1perK"], df["tau_CV_s"], "o-", label="Constant-Volume")
plt.semilogy(df["invT_1perK"], df["tau_CP_s"], "s-", label="Constant-Pressure")
plt.xlabel("Inverse Temperature 1/T₀ [1/K]")
plt.ylabel("Ignition Delay τ [s]")
plt.title(f"DME/Air Ignition Delay (ϕ={phi:.1f}, P₀={P0/ct.one_atm:.1f} atm)")
plt.grid(True, which="both", linestyle="--", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig("prob5_DME_tau_vs_invT.png", dpi=200)
plt.show()
