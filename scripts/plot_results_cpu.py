#!/usr/bin/env python3
"""
plot_results.py  –  Generate speedup/timing plots from run_sweep.sh CSV.

8 algorithms:
  Serial:   Full O(N^2), Cell-List (half-shell), Cell-List Manhattan, Neighbor-List
  Parallel: Full Pthread, Cell Half-Shell Pthread, Cell Manhattan Pthread, NBL Pthread

Usage:
    python3 plot_results.py ./results/speedups.csv ./results/plots
"""

import sys
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

if len(sys.argv) < 3:
    print("Usage: python3 plot_results.py <speedups.csv> <plots_dir>")
    sys.exit(1)

csv_path  = sys.argv[1]
plots_dir = sys.argv[2]
os.makedirs(plots_dir, exist_ok=True)

# ── load CSV ──────────────────────────────────────────────────────────────────
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

particle_counts = sorted(set(int(r["N"])       for r in rows))
thread_counts   = sorted(set(int(r["threads"])  for r in rows))
print(f"Particle counts : {particle_counts}")
print(f"Thread counts   : {thread_counts}")

# ── method tables ─────────────────────────────────────────────────────────────
# (key, display label, speedup_col, time_col)
SERIAL_METHODS = [
    ("cell",         "Cell-List (serial)",           "su_cell",         "time_cell"),
    ("cell_manhat",  "Cell-List Manhattan (serial)",  "su_cell_manhat",  "time_cell_manhat"),
    ("nbl",          "Neighbor-List (serial)",        "su_nbl",          "time_nbl"),
]
PARALLEL_METHODS = [
    ("full_pt",   "Full O(N²) Pthread",           "su_full_pt",   "time_full_pt"),
    ("cell_mt",   "Cell Half-Shell Pthread",       "su_cell_mt",   "time_cell_mt"),
    ("manh_pt",   "Cell Manhattan Pthread",        "su_manh_pt",   "time_manh_pt"),
    ("nbl_pt",    "NBL Pthread",                   "su_nbl_pt",    "time_nbl_pt"),
]
ALL_METHODS = SERIAL_METHODS + PARALLEL_METHODS

ALL_TIME_COLS = [
    ("time_full",        "Full O(N²) serial"),
    ("time_cell",        "Cell-List (serial)"),
    ("time_cell_manhat", "Cell-List Manhattan (serial)"),
    ("time_nbl",         "Neighbor-List (serial)"),
    ("time_full_pt",     "Full O(N²) Pthread"),
    ("time_cell_mt",     "Cell Half-Shell Pthread"),
    ("time_manh_pt",     "Cell Manhattan Pthread"),
    ("time_nbl_pt",      "NBL Pthread"),
]

colors        = plt.rcParams["axes.prop_cycle"].by_key()["color"]
thread_colors = cm.viridis(np.linspace(0.15, 0.85, len(thread_counts)))

max_t      = max(thread_counts)
t_max_rows = [r for r in rows if int(r["threads"]) == max_t]
t_max_by_N = {int(r["N"]): r for r in t_max_rows}

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 1 – Speedup vs N  (all 7 methods vs full serial, threads=max)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))

for idx, (key, label, su_col, _) in enumerate(ALL_METHODS):
    Ns, sus = [], []
    for N in particle_counts:
        row = t_max_by_N.get(N)
        if row is None:
            continue
        v = fval(row, su_col)
        if v is not None:
            Ns.append(N)
            sus.append(v)
    if Ns:
        ax.plot(Ns, sus, marker="o", label=label, color=colors[idx % len(colors)])

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="Baseline (Full O(N²) serial)")
ax.set_xscale("log", base=2)
ax.set_xticks(particle_counts)
ax.set_xticklabels([str(n) for n in particle_counts])
ax.set_xlabel("Number of particles (N)")
ax.set_ylabel("Speedup over Full O(N²) serial")
ax.set_title(f"Speedup vs Particle Count  [threads = {max_t}]")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
p = os.path.join(plots_dir, f"speedup_vs_N_t{max_t}.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 2 – Speedup vs N per parallel method, one line per thread count
# ═════════════════════════════════════════════════════════════════════════════
ncols = 2
nrows = (len(PARALLEL_METHODS) + 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(13, 5 * nrows), sharex=True)
axes = axes.flatten()

for ax_i, (key, label, su_col, _) in enumerate(PARALLEL_METHODS):
    ax = axes[ax_i]
    for t_i, T in enumerate(thread_counts):
        by_N = {int(r["N"]): r for r in rows if int(r["threads"]) == T}
        Ns, sus = [], []
        for N in particle_counts:
            row = by_N.get(N)
            if row is None:
                continue
            v = fval(row, su_col)
            if v is not None:
                Ns.append(N); sus.append(v)
        if Ns:
            ax.plot(Ns, sus, marker="o", label=f"{T} threads",
                    color=thread_colors[t_i])

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xscale("log", base=2)
    ax.set_xticks(particle_counts)
    ax.set_xticklabels([str(n) for n in particle_counts], rotation=30, ha="right")
    ax.set_title(label)
    ax.set_ylabel("Speedup")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3)

for ax_i in range(len(PARALLEL_METHODS), len(axes)):
    axes[ax_i].set_visible(False)

fig.suptitle("Parallel Method Speedup vs N  (by thread count)", fontsize=13)
fig.text(0.5, 0.01, "Number of particles (N)", ha="center")
plt.tight_layout()
p = os.path.join(plots_dir, "speedup_vs_N_by_threads.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 3 – Speedup vs Thread count, one subplot per N
# ═════════════════════════════════════════════════════════════════════════════
ncols = 3
nrows = (len(particle_counts) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows))
axes = axes.flatten()

for ax_i, N in enumerate(particle_counts):
    ax = axes[ax_i]
    by_T = {int(r["threads"]): r for r in rows if int(r["N"]) == N}

    for idx, (key, label, su_col, _) in enumerate(PARALLEL_METHODS):
        Ts, sus = [], []
        for T in thread_counts:
            row = by_T.get(T)
            if row is None:
                continue
            v = fval(row, su_col)
            if v is not None:
                Ts.append(T); sus.append(v)
        if Ts:
            ax.plot(Ts, sus, marker="o", label=label,
                    color=colors[idx % len(colors)])

    ax.plot(thread_counts, thread_counts, "k--", linewidth=0.8, label="Ideal linear")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xscale("log", base=2)
    ax.set_xticks(thread_counts)
    ax.set_xticklabels([str(t) for t in thread_counts])
    ax.set_title(f"N = {N}")
    ax.set_xlabel("Threads")
    ax.set_ylabel("Speedup")
    ax.legend(fontsize=6)
    ax.grid(True, which="both", alpha=0.3)

for ax_i in range(len(particle_counts), len(axes)):
    axes[ax_i].set_visible(False)

fig.suptitle("Parallel Speedup vs Thread Count  (per particle count)", fontsize=13)
plt.tight_layout()
p = os.path.join(plots_dir, "speedup_vs_threads_by_N.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 4 – Serial method speedup vs N (threads=1 run, clean comparison)
# ═════════════════════════════════════════════════════════════════════════════
t1_rows = [r for r in rows if int(r["threads"]) == 1]
t1_by_N = {int(r["N"]): r for r in t1_rows}

fig, ax = plt.subplots(figsize=(9, 6))
for idx, (key, label, su_col, _) in enumerate(SERIAL_METHODS):
    Ns, sus = [], []
    for N in particle_counts:
        row = t1_by_N.get(N)
        if row is None:
            continue
        v = fval(row, su_col)
        if v is not None:
            Ns.append(N); sus.append(v)
    if Ns:
        ax.plot(Ns, sus, marker="o", label=label, color=colors[idx % len(colors)])

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="Baseline")
ax.set_xscale("log", base=2)
ax.set_xticks(particle_counts)
ax.set_xticklabels([str(n) for n in particle_counts])
ax.set_xlabel("Number of particles (N)")
ax.set_ylabel("Speedup over Full O(N²) serial")
ax.set_title("Serial Algorithm Speedup vs Particle Count")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
p = os.path.join(plots_dir, "speedup_serial_vs_N.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 5 – Raw wall-clock times vs N, log-log  (threads=max)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
for idx, (col, label) in enumerate(ALL_TIME_COLS):
    Ns, times = [], []
    for N in particle_counts:
        row = t_max_by_N.get(N)
        if row is None:
            continue
        v = fval(row, col)
        if v is not None:
            Ns.append(N); times.append(v)
    if Ns:
        ax.plot(Ns, times, marker="o", label=label, color=colors[idx % len(colors)])

ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(particle_counts)
ax.set_xticklabels([str(n) for n in particle_counts])
ax.set_xlabel("Number of particles (N)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(f"Wall-clock Time vs N  [threads = {max_t}, log-log]")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
p = os.path.join(plots_dir, f"walltime_vs_N_t{max_t}.png")
plt.savefig(p, dpi=150); plt.close()
print(f"Saved: {p}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 6 – Fastest method chosen per (N, threads) as a heatmap/table
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 4))
ax.axis("off")

col_labels = [str(t) + "T" for t in thread_counts]
row_labels  = [f"N={n}" for n in particle_counts]
table_data  = []

for N in particle_counts:
    row_data = []
    for T in thread_counts:
        match = [r for r in rows if int(r["N"]) == N and int(r["threads"]) == T]
        if match:
            method = match[0].get("fastest_method", "N/A").strip('"')
            # shorten long names for table readability
            method = method.replace("CELL-LIST", "Cell").replace("NEIGHBOR-LIST", "NBL") \
                           .replace("MANHATTAN", "Manh").replace("PTHREAD", "PT") \
                           .replace("(Serial)", "ser").replace("(Pthread)", "PT") \
                           .replace("HALF-SHELL", "HS").replace("O(N^2)", "O(N²)")
            row_data.append(method)
        else:
            row_data.append("—")
    table_data.append(row_data)

tbl = ax.table(
    cellText=table_data,
    rowLabels=row_labels,
    colLabels=col_labels,
    loc="center",
    cellLoc="center"
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
tbl.scale(1.4, 1.6)
ax.set_title("Fastest Method Selected by Autotuner  (rows=N, cols=threads)",
             fontsize=11, pad=12)
plt.tight_layout()
p = os.path.join(plots_dir, "fastest_method_table.png")
plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
print(f"Saved: {p}")

print("\nAll plots generated successfully.")
