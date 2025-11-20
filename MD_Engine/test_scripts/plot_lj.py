import matplotlib.pyplot as plt
import numpy as np
import os

# =============================
# Ensure output folder exists
# =============================
os.makedirs("output", exist_ok=True)

# =============================
# Load energy CSV
# =============================
def load_energy_csv(path):
    if not os.path.exists(path):
        print(f" Missing: {path}")
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

    if len(steps) == 0:
        print(f"⚠ No data in {path}")
        return None

    return np.array(steps), np.array(K), np.array(U), np.array(E)


# =============================
# File definitions (NOW 4)
# =============================
files = {
    "full":  "output/full_energies.csv",
    "cell":  "output/cell_energies.csv",
    "nbl":   "output/nbl_energies.csv",
    "final": "output/final_energies.csv"
}

# =============================
# Create main 4-panel figure
# =============================
fig, axes = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

final_data = None

for ax, (name, path) in zip(axes, files.items()):
    data = load_energy_csv(path)
    if data is None:
        ax.text(0.5, 0.5, f"No data for {name}", ha='center', va='center')
        continue

    steps, K, U, E = data

    # Save FINAL dataset for log plot later
    if name == "final":
        final_data = data

    # =============================
    # Energy conservation check
    # =============================
    E0 = E[0]
    Ef = E[-1]
    dE = Ef - E0
    pct = (dE / E0) * 100

    print(f"\n=== {name.upper()} ENERGY CHECK ===")
    print(f"Initial E0: {E0:.6e}")
    print(f"Final Ef : {Ef:.6e}")
    print(f"ΔE       : {dE:.6e}")
    print(f"ΔE/E0    : {pct:.4f} %")

    # =============================
    # Plot energies
    # =============================
    ax.plot(steps, K, label="Kinetic", alpha=0.8)
    ax.plot(steps, U, label="Potential", alpha=0.8)
    ax.plot(steps, E, label="Total", linewidth=2.0)

    ax.text(
        0.02, 0.95,
        f"Drift = {pct:.3f}%",
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.set_title(f"{name.upper()} Energy vs Step")
    ax.set_ylabel("Energy")
    ax.grid(True)
    ax.legend()

axes[-1].set_xlabel("Step")

plt.tight_layout()
plt.savefig("output/energy_four_panel.png", dpi=300)
print("\n Saved → output/energy_four_panel.png\n")

# ==========================================================
# EXTRA PLOT: Log plot of FINAL energies
# ==========================================================
if final_data is not None:
    steps, K, U, E = final_data

    eps = 1e-12
    Klog = np.log10(np.abs(K) + eps)
    Ulog = np.log10(np.abs(U) + eps)
    Elog = np.log10(np.abs(E) + eps)

    # =============================
    # Print log-scale info to terminal
    # =============================
    print("\n===== FINAL RUN LOG-SCALE ENERGY INFO =====")
    print(f"log10|K| range: {np.min(Klog):.4f} to {np.max(Klog):.4f}")
    print(f"log10|U| range: {np.min(Ulog):.4f} to {np.max(Ulog):.4f}")
    print(f"log10|E| range: {np.min(Elog):.4f} to {np.max(Elog):.4f}")
    print(f"Max K spike magnitude (orders): {np.max(Klog) - np.min(Klog):.4f}")
    print(f"Max U spike magnitude (orders): {np.max(Ulog) - np.min(Ulog):.4f}")
    print(f"Total energy flatness (orders): {np.max(Elog) - np.min(Elog):.4f}")

    # =============================
    # Plot log10 energies
    # =============================
    plt.figure(figsize=(12, 6))
    plt.plot(steps, Klog, label="log10|K|", alpha=0.8)
    plt.plot(steps, Ulog, label="log10|U|", alpha=0.8)
    plt.plot(steps, Elog, label="log10|E|", linewidth=2.0)

    plt.title("FINAL Simulation: Log-Scale Energies (log10|E|)")
    plt.xlabel("Step")
    plt.ylabel("log10(|Energy|)")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig("output/energy_final_logplot.png", dpi=300)
    print("\n Saved → output/energy_final_logplot.png\n")

else:
    print("\n⚠ FINAL energy file missing — log plot skipped.\n")
