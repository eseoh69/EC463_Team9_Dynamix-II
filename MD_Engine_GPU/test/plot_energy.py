import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("output", exist_ok=True)

# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

def load_energy_csv(path):
    if not os.path.exists(path):
        return None
    steps, K, U, E = [], [], [], []
    with open(path) as f:
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
            except ValueError:
                continue
    if not steps:
        return None
    return np.array(steps), np.array(K), np.array(U), np.array(E)

# ---------------------------------------------------------------------------
# Energy drift summary
# ---------------------------------------------------------------------------

def energy_drift(label, E, f):
    E0, Ef = E[0], E[-1]
    dE = Ef - E0
    pct = (dE / E0) * 100 if E0 != 0 else float("nan")
    f.write(f"--- {label} ---\n")
    f.write(f"  Initial E0 : {E0:.6e}\n")
    f.write(f"  Final   Ef : {Ef:.6e}\n")
    f.write(f"  ΔE         : {dE:.6e}\n")
    f.write(f"  ΔE/E0      : {pct:.4f}%\n\n")

# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------

datasets = [
    ("Full O(N²)",    "output/full_energies.csv"),
    ("Cell-list",     "output/cell_energies.csv"),
    ("Neighbor-list", "output/nbl_energies.csv"),
    ("Manhattan-cell-list", "output/manhattan_cell_energies.csv"),
    ("Final run",     "output/final_energies.csv"),
]

loaded = [(label, load_energy_csv(path)) for label, path in datasets]
loaded = [(label, data) for label, data in loaded if data is not None]

if not loaded:
    print("No energy CSV files found in output/. Run the simulation first.")
    exit(1)

# ---------------------------------------------------------------------------
# Drift summary
# ---------------------------------------------------------------------------

with open("output/energy_summary.txt", "w") as f:
    f.write("=== GPU MD ENGINE — ENERGY DRIFT SUMMARY ===\n\n")
    for label, (steps, K, U, E) in loaded:
        energy_drift(label, E, f)

print("Saved: output/energy_summary.txt")

# ---------------------------------------------------------------------------
# Per-algorithm plots (K, U, E on one figure)
# ---------------------------------------------------------------------------

for label, (steps, K, U, E) in loaded:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(steps, K, label="Kinetic (K)")
    ax.plot(steps, U, label="Potential (U)")
    ax.plot(steps, E, label="Total (E)", linewidth=2)
    ax.set_title(f"{label} — Energy vs Timestep")
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    fname = "output/energy_" + label.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("²", "2") + ".png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"Saved: {fname}")
