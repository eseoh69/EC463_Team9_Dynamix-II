#!/usr/bin/env python3
import os
import math
import argparse
from typing import List, Tuple, Optional

Vec3 = Tuple[float, float, float]

RUN_PARAMS = {
    "N": 1024,
    "input": "input/random_particles-1024.pdb",
    "bounds": "x[0.540,299.724] y[0.561,299.627] z[0.044,299.864]",
    "L": 299.820,
    "sigma": 1.000,
    "epsilon": 1.000,
    "rc": 30.000,
    "dt": 0.001000,
    "steps": 4000,
    "gpu_timing_total_s": 1.519203,
    "gpu_ms_per_step": 0.37980075,
}

def _try_parse_xyz_columns(line: str) -> Optional[Vec3]:
    if len(line) < 54:
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return (x, y, z)
    except ValueError:
        return None

def _try_parse_xyz_floats(line: str) -> Optional[Vec3]:
    vals = []
    s = line.replace(",", " ")
    for tok in s.split():
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    if len(vals) < 3:
        return None
    return (vals[0], vals[1], vals[2])

def read_pdb_coords(path: str) -> List[Vec3]:
    coords: List[Vec3] = []
    with open(path, "r") as f:
        for line in f:
            rec = line[0:6].strip()  # robust against leading spaces
            if rec not in ("ATOM", "HETATM"):
                continue

            xyz = _try_parse_xyz_columns(line)
            if xyz is None:
                xyz = _try_parse_xyz_floats(line)
            if xyz is None:
                continue

            coords.append(xyz)

    if not coords:
        with open(path, "r") as f:
            head = "".join([next(f, "") for _ in range(8)])
        raise RuntimeError(
            f"No ATOM/HETATM coords parsed from {path}\n"
            f"First lines:\n{head}"
        )
    return coords

def _min_image(d: float, L: float) -> float:
    return d - round(d / L) * L

def compare_coords(a: List[Vec3], b: List[Vec3], pbc_L: Optional[float] = None) -> dict:
    if len(a) != len(b):
        raise RuntimeError(f"Atom count mismatch: {len(a)} vs {len(b)}")

    n = len(a)
    sum_sq = 0.0
    percent_error_sum = 0.0
    max_dist = 0.0

    for (ax, ay, az), (bx, by, bz) in zip(a, b):
        dx = ax - bx
        dy = ay - by
        dz = az - bz

        if pbc_L is not None:
            dx = _min_image(dx, pbc_L)
            dy = _min_image(dy, pbc_L)
            dz = _min_image(dz, pbc_L)

        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        sum_sq += dist * dist
        if dist > max_dist:
            max_dist = dist

        pos1 = math.sqrt(ax*ax + ay*ay + az*az)
        pos2 = math.sqrt(bx*bx + by*by + bz*bz)
        if pos1 > 1e-12:
            percent_error_sum += abs(pos1 - pos2) / pos1

    rmsd = math.sqrt(sum_sq / n)
    mean_percent_error = percent_error_sum / n
    return {
        "n": n,
        "rmsd": rmsd,
        "mean_percent_error": mean_percent_error,
        "max_dist": max_dist,
    }

def print_table(rows: List[dict]) -> None:
    headers = ["file", "atoms", "RMSD", "Mean%Err", "MaxDist"]
    col0 = max(len(headers[0]), *(len(r["file"]) for r in rows))
    fmt = f"{{:<{col0}}}  {{:>5}}  {{:>12}}  {{:>12}}  {{:>12}}"

    print(fmt.format(*headers))
    print("-" * (col0 + 5 + 12 + 12 + 12 + 10))

    for r in rows:
        print(fmt.format(
            r["file"],
            r["n"],
            f"{r['rmsd']:.6g}",
            f"{(r['mean_percent_error']*100):.6g}",
            f"{r['max_dist']:.6g}",
        ))

def main():
    ap = argparse.ArgumentParser(
        description="Compare reference PDB to all output PDBs in output_dir. "
                    "Use --pbc L to compare with periodic minimum-image distances."
    )
    ap.add_argument("--ref", default="output/full_positions.pdb",
                    help="Reference PDB path (default: output/full_positions.pdb)")
    ap.add_argument("--output_dir", default="output",
                    help="Directory containing PDBs (default: output)")
    ap.add_argument("--pbc", type=float, default=None,
                    help="Periodic box length L for minimum-image comparison (optional)")
    ap.add_argument("--all_positions", action="store_true",
                    help="If set, compares ALL '*positions*.pdb' files found in output_dir (recommended).")
    args = ap.parse_args()

    print("=== Run Params (for context) ===")
    for k in ["N", "input", "bounds", "L", "sigma", "epsilon", "rc", "dt", "steps", "gpu_timing_total_s", "gpu_ms_per_step"]:
        if k in RUN_PARAMS:
            print(f"{k}: {RUN_PARAMS[k]}")
    if args.pbc is not None:
        print(f"PBC: ON (L={args.pbc})")
    else:
        print("PBC: OFF")
    print("===============================\n")

    if not os.path.exists(args.ref):
        raise SystemExit(f"Reference not found: {args.ref}")

    ref = read_pdb_coords(args.ref)

    # Default list (includes GPU output now)
    default_candidates = [
        "full_positions.pdb",
        "cell_positions.pdb",
        "cell_half_mt_positions.pdb",
        "nbl_positions.pdb",
        "nbl_pthread_positions.pdb",
        "manhattan_pthread_positions.pdb",
        "final_positions.pdb",
        "gpu_quick_positions.pdb",
    ]

    if args.all_positions:
        candidates = sorted([
            f for f in os.listdir(args.output_dir)
            if f.lower().endswith(".pdb") and "positions" in f.lower()
        ])
    else:
        candidates = default_candidates

    rows = []
    for fname in candidates:
        path = os.path.join(args.output_dir, fname)
        if not os.path.exists(path):
            continue
        try:
            coords = read_pdb_coords(path)
            stats = compare_coords(ref, coords, pbc_L=args.pbc)
            stats["file"] = fname
            rows.append(stats)
        except Exception as e:
            print(f"[warn] Failed comparing {fname}: {e}")

    if not rows:
        raise SystemExit("No comparable PDBs found. Did you generate output/*.pdb?")

    # Sort with reference first if it's inside output_dir, otherwise alphabetical
    ref_base = os.path.basename(args.ref)
    rows.sort(key=lambda r: (r["file"] != ref_base, r["file"]))
    print_table(rows)

if __name__ == "__main__":
    main()