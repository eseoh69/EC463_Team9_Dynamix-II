import numpy as np
import re

# ============================================================
# Load coordinates from PDB
# ============================================================
def load_pdb(path):
    coords = []
    atom_pattern = re.compile(r'(ATOM|HETATM)')
    try:
        with open(path, 'r') as f:
            for line in f:
                if not atom_pattern.search(line):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    parts = line.split()
                    nums = [p for p in parts if re.match(r'^-?\d+(\.\d*)?$', p)]
                    if len(nums) >= 3:
                        coords.append([float(nums[-3]), float(nums[-2]), float(nums[-1])])
    except FileNotFoundError:
        return np.zeros((0, 3))
    return np.array(coords)


# ============================================================
# RMSD metrics
# ============================================================
def compute_metrics(A, B):
    n = min(len(A), len(B))
    if n == 0:
        return None
    diff = A[:n] - B[:n]
    dists = np.linalg.norm(diff, axis=1)
    return {
        "N":              n,
        "RMSD":           float(np.sqrt((dists ** 2).mean())),
        "Mean abs error": float(dists.mean()),
        "Max deviation":  float(dists.max()),
    }


# ============================================================
# All output files
# ============================================================
outputs = {
    "full":             "output/full_positions.pdb",
    "cell":             "output/cell_positions.pdb",
    "nbl":              "output/nbl_positions.pdb",
    "full_pthread":     "output/full_pthread_positions.pdb",
    "cell_half_mt":     "output/cell_half_mt_positions.pdb",
    "nbl_pthread":      "output/nbl_pthread_positions.pdb",
    "manhattan_pthread":"output/manhattan_pthread_positions.pdb",
    "manhat":           "output/manhat_positions.pdb",
    "final":            "output/final_positions.pdb",
}

print("=== Loading PDB files ===")
data = {}
for name, path in outputs.items():
    arr = load_pdb(path)
    status = f"{len(arr)} atoms" if len(arr) > 0 else "MISSING"
    print(f"  {name:22s}: {status}  ({path})")
    data[name] = arr
print()

# ============================================================
# Section 1: All methods vs FULL (reference)
# ============================================================
print("=" * 60)
print("1) All methods vs FULL O(N^2) (reference)")
print("=" * 60)
ref = "full"
for name, arr in data.items():
    if name == ref:
        continue
    m = compute_metrics(data[ref], arr)
    if m is None:
        print(f"  {name:22s}: MISSING or empty")
        continue
    print(f"\n  {name:22s} vs {ref}:")
    for k, v in m.items():
        print(f"    {k:16s}: {v:.6f}" if isinstance(v, float) else f"    {k:16s}: {v}")

# ============================================================
# Section 2: Serial vs Pthread pairs (correctness check)
# ============================================================
print()
print("=" * 60)
print("2) Serial vs Pthread equivalents (correctness check)")
print("=" * 60)

pairs = [
    ("full",        "full_pthread",      "FULL serial vs FULL pthread"),
    ("cell",        "cell_half_mt",      "CELL serial vs CELL half-shell pthread"),
    ("nbl",         "nbl_pthread",       "NBL serial vs NBL pthread"),
    ("manhattan_pthread", "manhat",      "Manhattan pthread vs Manhattan serial"),
]

for A, B, label in pairs:
    m = compute_metrics(data[A], data[B])
    if m is None:
        print(f"\n  {label}: MISSING data")
        continue
    print(f"\n  {label}:")
    for k, v in m.items():
        print(f"    {k:16s}: {v:.6f}" if isinstance(v, float) else f"    {k:16s}: {v}")

# ============================================================
# Section 3: Final vs all methods
# ============================================================
print()
print("=" * 60)
print("3) FINAL simulation vs all autotuning outputs")
print("=" * 60)

for name, arr in data.items():
    if name == "final":
        continue
    m = compute_metrics(data["final"], arr)
    if m is None:
        print(f"  final vs {name:22s}: MISSING or empty")
        continue
    print(f"\n  final vs {name:22s}:")
    for k, v in m.items():
        print(f"    {k:16s}: {v:.6f}" if isinstance(v, float) else f"    {k:16s}: {v}")

print("\nDone.\n")
