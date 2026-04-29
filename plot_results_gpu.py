#!/usr/bin/env python3
"""
plot_results_gpu.py – Generate plots from run_sweep_gpu.sh CSV.

4 GPU algorithms timed per particle count:
  Full O(N^2), Cell-List, Neighbor-List, Manhattan

Usage:
    python3 plot_results_gpu.py ./results/timings.csv ./results/plots
"""

import sys
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 3:
    print("Usage: python3 plot_results_gpu.py <timings.csv> <plots_dir>")
    sys.exit(1)

csv_path  = sys.argv[1]
plots_dir = sys.argv[2]
os.makedirs(plots_dir, exist_ok=True)

rows = []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        rows.append(row)

if not rows:
    print("ERROR: CSV is empty.")
    sys.exit(1)

print(f"Loaded {len(rows)} rows from {csv_path}")

def fval(row, key):
    try:
        return float(row.get(key, ""))
    except (ValueError, TypeError):
        return None

particle_counts = sorted(int(r["N"]) for r in rows)
print(f"Particle counts: {particle_counts}")

METHODS = [
    ("time_full", "su_full", "Full O(N²) GPU",    "steelblue"),
    ("time_cell", "su_cell", "Cell-List GPU",      "darkorange"),
    ("time_nbl",  "su_nbl",  "Neighbor-List GPU",  "seagreen"),
    ("time_man",  "su_man",  "Manhattan GPU",       "crimson"),
]

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# ── PLOT 1: Wall-clock time vs N (log-log, all 4 methods) ────────────────────
fig, ax = plt.subplots(figsize=(10, 6))

for time_col, _, label, color in METHODS:
    Ns, ts = [], []
    for r in rows:
        v = fval(r, time_col)
        if v is not None:
            Ns.append(int(r["N"])); ts.append(v)
    if Ns:
        ax.plot(Ns, ts, marker="o", label=label, color=color, linewidth=2, markersize=6)

# O(N) / O(N^2) reference lines anchored at first full-method point
full_pts = [(int(r["N"]), fval(r, "time_full")) for r in rows if fval(r, "time_full") is not None]
if full_pts:
    n0, t0 = full_pts[0]
    ref_ns = np.array(particle_counts, dtype=float)
    ax.plot(ref_ns, t0 * (ref_ns / n0),    "k--", linewidth=0.8, label="O(N) ref")
    ax.plot(ref_ns, t0 * (ref_ns / n0)**2, "k:",  linewidth=0.8, label="O(N²) ref")

ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(particle_counts)
ax.set_xticklabels([str(n) for n in particle_counts])
ax.set_xlabel("Number of particles (N)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title("GPU Methods — Wall-clock Time vs N  (log-log)")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
p = os.path.join(plots_dir, "gpu_walltime_vs_N.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ── PLOT 2: Speedup vs N (cell, nbl, manhattan over full O(N^2)) ─────────────
fig, ax = plt.subplots(figsize=(10, 6))

SPEEDUP_METHODS = [
    ("su_cell", "Cell-List GPU",     "darkorange"),
    ("su_nbl",  "Neighbor-List GPU", "seagreen"),
    ("su_man",  "Manhattan GPU",     "crimson"),
]

for su_col, label, color in SPEEDUP_METHODS:
    Ns, sus = [], []
    for r in rows:
        v = fval(r, su_col)
        if v is not None:
            Ns.append(int(r["N"])); sus.append(v)
    if Ns:
        ax.plot(Ns, sus, marker="o", label=label, color=color, linewidth=2, markersize=6)

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="Baseline (Full O(N²) GPU)")
ax.set_xscale("log", base=2)
ax.set_xticks(particle_counts)
ax.set_xticklabels([str(n) for n in particle_counts])
ax.set_xlabel("Number of particles (N)")
ax.set_ylabel("Speedup over Full O(N²) GPU")
ax.set_title("GPU Methods — Speedup vs Particle Count")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
p = os.path.join(plots_dir, "gpu_speedup_vs_N.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ── PLOT 3: Fastest method table (rows = N) ───────────────────────────────────
fig, ax = plt.subplots(figsize=(8, max(3, len(particle_counts) * 0.7)))
ax.axis("off")

row_labels = [f"N={n}" for n in particle_counts]
table_data = []
col_labels = ["Fastest Method", "Time (s)"]

N_to_row = {int(r["N"]): r for r in rows}
for N in particle_counts:
    r = N_to_row.get(N)
    if r:
        method = r.get("fastest_method", "N/A").strip('"')
        ft     = r.get("fastest_time",   "N/A")
        table_data.append([method, ft])
    else:
        table_data.append(["—", "—"])

tbl = ax.table(
    cellText=table_data,
    rowLabels=row_labels,
    colLabels=col_labels,
    loc="center",
    cellLoc="center"
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.6, 1.8)
ax.set_title("Fastest GPU Method per Particle Count", fontsize=11, pad=12)
plt.tight_layout()
p = os.path.join(plots_dir, "gpu_fastest_method_table.png")
plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
print(f"Saved: {p}")

print("\nAll GPU plots generated successfully.")
