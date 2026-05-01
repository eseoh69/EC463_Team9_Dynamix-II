import matplotlib
matplotlib.use("Agg")  # non-interactive backend — no display needed
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("output", exist_ok=True)

summary_path = "output/energy_summary.txt"
summary = open(summary_path, "w")
summary.write("=== ENERGY SUMMARY REPORT ===\n\n")

def load_energy_csv(path):
    if not os.path.exists(path):
        return None
    steps, K, U, E = [], [], [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("step"):
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            try:
                steps.append(int(parts[0]))
                K.append(float(parts[1]))
                U.append(float(parts[2]))
                E.append(float(parts[3]))
            except:
                continue
    if not steps:
        return None
    return np.array(steps), np.array(K), np.array(U), np.array(E)


files = {
    "full":                   "./output/full_energies.csv",
    "cell":                   "./output/cell_energies.csv",
    "cell_half_mt":           "./output/cell_half_mt_energies.csv",
    "nbl":                    "./output/nbl_energies.csv",
    "full_pthread":           "./output/full_pthread_energies.csv",
    "cell_half_pthread":      "./output/cell_half_mt_energies.csv",
    "cell_manhattan_pthread": "./output/manhattan_pthread_energies.csv",
    "nbl_pthread":            "./output/nbl_pthread_energies.csv",
    "manhat":                 "./output/cell_manhat_energies.csv",
    "final":                  "./output/final_energies.csv",
}

# ============================================================
# Per-method figure: K on left axis, U on right axis, ΔE% inset
# ============================================================
n_methods = len(files)
fig, axes = plt.subplots(n_methods, 1, figsize=(12, 3 * n_methods))

final_data = None

for ax, (name, path) in zip(axes, files.items()):
    data = load_energy_csv(path)
    if data is None:
        ax.text(0.5, 0.5, f"No data for {name}", ha='center', va='center',
                transform=ax.transAxes)
        summary.write(f"{name.upper()}: MISSING FILE ({path})\n\n")
        continue

    steps, K, U, E = data
    if name == "final":
        final_data = data

    # Compute changes from initial value
    dK = K - K[0]
    dU = U - U[0]
    dE = E - E[0]

    E0  = E[0]
    pct = (dE[-1] / abs(E0)) * 100 if E0 != 0 else 0.0

    summary.write(f"--- {name.upper()} ENERGY CHECK ---\n")
    summary.write(f"  Initial E0 : {E0:.6e}\n")
    summary.write(f"  Final Ef   : {E[-1]:.6e}\n")
    summary.write(f"  ΔE         : {dE[-1]:.6e}\n")
    summary.write(f"  ΔE/|E0|    : {pct:.6f} %\n\n")

    # ΔK and ΔU on same axis — should mirror each other if energy is conserved
    ax.plot(steps, dK, color="tab:blue",   alpha=0.85, label="ΔK (kinetic)",   linewidth=1.2)
    ax.plot(steps, dU, color="tab:orange", alpha=0.85, label="ΔU (potential)", linewidth=1.2)
    ax.plot(steps, dE, color="tab:green",  alpha=0.7,  label="ΔE (total)",     linewidth=1.0, linestyle="--")
    ax.axhline(0, color="black", linewidth=0.6, linestyle=":")

    ax.set_ylabel("Energy change", fontsize=9)
    ax.set_title(f"{name.upper()}  —  ΔK, ΔU, ΔE from step 0  |  ΔE/|E₀| = {pct:.6f}%", fontsize=9)
    ax.set_xlabel("Step")
    ax.grid(True, alpha=0.25)
    ax.legend(loc='upper right', fontsize=8)

axes[-1].set_xlabel("Step")
plt.tight_layout()
plt.savefig("output/energy_all_comparison.png", dpi=100)
plt.close(fig)
print("Comparison saved: output/energy_all_comparison.png")

# ============================================================
# Summary panel: ΔE/|E0| drift across all methods (bar chart)
# ============================================================
drift_names = []
drift_vals  = []

for name, path in files.items():
    data = load_energy_csv(path)
    if data is None:
        continue
    _, _, _, E = data
    if E[0] != 0:
        pct = abs((E[-1] - E[0]) / abs(E[0])) * 100
        drift_names.append(name)
        drift_vals.append(pct)

if drift_names:
    fig2, ax3 = plt.subplots(figsize=(12, 5))
    bars = ax3.bar(drift_names, drift_vals, color='steelblue', edgecolor='black', alpha=0.8)
    ax3.set_ylabel("|ΔE/E₀| (%)")
    ax3.set_title("Energy Drift (%) per Method — lower is better")
    ax3.set_xticklabels(drift_names, rotation=30, ha='right', fontsize=9)
    ax3.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, drift_vals):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                 f"{val:.4f}%", ha='center', va='bottom', fontsize=8)
    plt.tight_layout()
    plt.savefig("output/energy_drift_comparison.png", dpi=100)
    plt.close(fig2)
    print("Drift comparison saved: output/energy_drift_comparison.png")

# ============================================================
# Log-scale plot for FINAL
# ============================================================
if final_data is not None:
    steps, K, U, E = final_data
    eps = 1e-12

    dK = K - K[0]
    dU = U - U[0]
    dE = E - E[0]

    fig3, (axL, axR) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: ΔK and ΔU — should mirror each other
    axL.plot(steps, dK, color='tab:blue',   label='ΔK', linewidth=1.2)
    axL.plot(steps, dU, color='tab:orange', label='ΔU', linewidth=1.2)
    axL.plot(steps, dE, color='tab:green',  label='ΔE (total)', linewidth=1.0, linestyle='--')
    axL.axhline(0, color='black', linewidth=0.6, linestyle=':')
    axL.set_title("FINAL: ΔK, ΔU, ΔE from step 0\n(ΔK ≈ −ΔU means energy conserved)")
    axL.set_xlabel("Step")
    axL.set_ylabel("Energy change")
    axL.legend(fontsize=9)
    axL.grid(True, alpha=0.25)

    # Right: |ΔE| on log scale — conservation quality
    dE_abs = np.abs(dE) + eps
    axR.semilogy(steps, dE_abs, color='tab:green', linewidth=1.2)
    axR.set_title("FINAL: |ΔE| drift (log scale)\nlower = better conservation")
    axR.set_xlabel("Step")
    axR.set_ylabel("|E(t) − E(0)|")
    axR.grid(True, alpha=0.25)

    plt.tight_layout()
    plt.savefig("output/energy_final_logplot.png", dpi=100)
    plt.close(fig3)
    print("Final log-plot saved: output/energy_final_logplot.png")

summary.close()
print(f"Energy summary written to: {summary_path}")
