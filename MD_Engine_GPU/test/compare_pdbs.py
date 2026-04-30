import numpy as np
import os
from itertools import combinations

os.makedirs("output", exist_ok=True)

def load_pdb(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                coords.append((x, y, z))
    return np.array(coords)

datasets = [
    ("Full O(N2)",    "output/full_positions.pdb"),
    ("Cell-list",     "output/cell_positions.pdb"),
    ("Neighbor-list", "output/nbl_positions.pdb"),
    ("Manhattan-cell-list", "output/manhattan_positions.pdb"),
    ("Final run",     "output/final_positions.pdb"),
]

loaded = [(label, load_pdb(path)) for label, path in datasets if os.path.exists(path)]

if not loaded:
    print("No PDB files found in output/. Run the simulation first.")
    exit(1)

with open("output/pdb_comparison.txt", "w") as f:
    f.write("=== PDB POSITION COMPARISON ===\n\n")
    for (la, ca), (lb, cb) in combinations(loaded, 2):
        if len(ca) != len(cb):
            f.write(f"{la} vs {lb}: particle count mismatch ({len(ca)} vs {len(cb)})\n\n")
            continue
        diff = np.sqrt(np.sum((ca - cb)**2, axis=1))
        ref  = np.sqrt(np.sum(ca**2, axis=1))
        pct  = (diff / np.where(ref > 0, ref, 1.0)) * 100.0
        f.write(f"{la} vs {lb}:\n")
        f.write(f"  mean deviation : {pct.mean():.4f}%\n")
        f.write(f"  max  deviation : {pct.max():.4f}%\n")
        f.write(f"  std  deviation : {pct.std():.4f}%\n\n")

print("Saved: output/pdb_comparison.txt")
