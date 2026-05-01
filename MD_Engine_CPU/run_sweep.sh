#!/usr/bin/env bash
# =============================================================================
# run_sweep.sh
#
# Sweeps the MD autotuning framework over:
#   - Particle counts : 1024 2048 4096 8192 16384 32768
#   - Thread counts   : 1 2 4 8 16
#
# 8 algorithms timed per run:
#   Serial:   Full O(N^2), Cell-List, Cell-List Manhattan, Neighbor-List
#   Parallel: Full Pthread, Cell Half-Shell Pthread, Cell Manhattan Pthread, NBL Pthread
#
# Usage:
#   chmod +x run_sweep.sh
#   ./run_sweep.sh [--binary PATH] [--input-dir PATH] [--output-dir PATH]
#
# Defaults:
#   --binary     ./md_run
#   --input-dir  ./input
#   --output-dir ./results
# =============================================================================

set -euo pipefail

# ── defaults ─────────────────────────────────────────────────────────────────
BINARY="./md_run"
INPUT_DIR="./input"
OUTPUT_DIR="./results"

PARTICLE_COUNTS=(1024 2048 4096 8192 16384 32768)
THREAD_COUNTS=(1 2 4 8 16)

# ── argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --binary)     BINARY="$2";     shift 2 ;;
        --input-dir)  INPUT_DIR="$2";  shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── sanity checks ─────────────────────────────────────────────────────────────
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: binary not found or not executable: $BINARY"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/plots"
mkdir -p output   # md_run writes energy CSVs and PDBs here

CSV="$OUTPUT_DIR/speedups.csv"

# ── CSV header ────────────────────────────────────────────────────────────────
# Serial:   full, cell, cell_manhat, nbl
# Parallel: full_pt, cell_mt, manh_pt, nbl_pt
cat > "$CSV" << 'EOF'
N,threads,time_full,time_cell,time_cell_manhat,time_nbl,time_full_pt,time_cell_mt,time_manh_pt,time_nbl_pt,su_cell,su_cell_manhat,su_nbl,su_full_pt,su_cell_mt,su_manh_pt,su_nbl_pt,fastest_method,fastest_time
EOF

# ── helpers ───────────────────────────────────────────────────────────────────
parse_time() {
    # $1 = log file, $2 = grep regex
    # returns first float on the matching line (the wall-clock time)
    grep -m1 -E "$2" "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | head -1
}

parse_speedup() {
    # returns second float on matching line (the Yx speedup value)
    grep -m1 -E "$2" "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | sed -n '2p'
}

parse_fastest_method() {
    grep -m1 "FASTEST METHOD:" "$1" \
        | sed 's/.*FASTEST METHOD: //' \
        | sed 's/ (.*//'
}

parse_fastest_time() {
    grep -m1 "FASTEST METHOD:" "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | head -1
}

# ── main sweep ────────────────────────────────────────────────────────────────
TOTAL=$(( ${#PARTICLE_COUNTS[@]} * ${#THREAD_COUNTS[@]} ))
RUN=0

for N in "${PARTICLE_COUNTS[@]}"; do
    PDB="$INPUT_DIR/random_particles-${N}.pdb"
    if [[ ! -f "$PDB" ]]; then
        echo "WARNING: input file not found, skipping: $PDB"
        continue
    fi

    for T in "${THREAD_COUNTS[@]}"; do
        RUN=$(( RUN + 1 ))
        LOG="$OUTPUT_DIR/logs/N${N}_T${T}.log"

        echo "──────────────────────────────────────────────────────"
        echo "Run ${RUN}/${TOTAL}  |  N=${N}  threads=${T}"
        echo "  binary : $BINARY"
        echo "  input  : $PDB"
        echo "  log    : $LOG"
        echo "──────────────────────────────────────────────────────"

        "$BINARY" "$PDB" --threads "$T" 2>&1 | tee "$LOG"

        # ── parse the 8 method lines ──────────────────────────────
        # Patterns match your exact printf strings:
        # "FULL O(N^2):                      %.6f seconds (1.00x speedup)"
        # "CELL-LIST:                        %.6f seconds (%.2fx speedup)"
        # "CELL-LIST MANHATTAN:              %.6f seconds (%.2fx speedup)"
        # "NEIGHBOR-LIST:                    %.6f seconds (%.2fx speedup)"
        # "FULL O(N^2) PTHREADS              %.6f seconds (%.2fx speedup)"
        # "CELL-LIST HALF-SHELL PTHREADS:    %.6f seconds (%.2fx speedup)"
        # "CELL-LIST MANHATTAN PTHREADS:     %.6f seconds (%.2fx speedup)"
        # "NBL PTHREADS:                     %.6f seconds (%.2fx speedup)"

        T_FULL=$(        parse_time "$LOG" "FULL O\(N\^2\):" )
        T_CELL=$(        parse_time "$LOG" "CELL-LIST:" )
        T_CELL_MANHAT=$( parse_time "$LOG" "CELL-LIST MANHATTAN:" )
        T_NBL=$(         parse_time "$LOG" "NEIGHBOR-LIST:" )
        T_FULL_PT=$(     parse_time "$LOG" "FULL O\(N\^2\) PTHREADS" )
        T_CELL_MT=$(     parse_time "$LOG" "CELL-LIST HALF-SHELL PTHREADS" )
        T_MANH_PT=$(     parse_time "$LOG" "CELL-LIST MANHATTAN PTHREADS" )
        T_NBL_PT=$(      parse_time "$LOG" "NBL PTHREADS" )

        SU_CELL=$(        parse_speedup "$LOG" "CELL-LIST:" )
        SU_CELL_MANHAT=$( parse_speedup "$LOG" "CELL-LIST MANHATTAN:" )
        SU_NBL=$(         parse_speedup "$LOG" "NEIGHBOR-LIST:" )
        SU_FULL_PT=$(     parse_speedup "$LOG" "FULL O\(N\^2\) PTHREADS" )
        SU_CELL_MT=$(     parse_speedup "$LOG" "CELL-LIST HALF-SHELL PTHREADS" )
        SU_MANH_PT=$(     parse_speedup "$LOG" "CELL-LIST MANHATTAN PTHREADS" )
        SU_NBL_PT=$(      parse_speedup "$LOG" "NBL PTHREADS" )

        FASTEST_METHOD=$(parse_fastest_method "$LOG")
        FASTEST_TIME=$(  parse_fastest_time   "$LOG")

        # Fallback to N/A if any parse failed
        T_FULL=${T_FULL:-N/A};             T_CELL=${T_CELL:-N/A}
        T_CELL_MANHAT=${T_CELL_MANHAT:-N/A}; T_NBL=${T_NBL:-N/A}
        T_FULL_PT=${T_FULL_PT:-N/A};       T_CELL_MT=${T_CELL_MT:-N/A}
        T_MANH_PT=${T_MANH_PT:-N/A};       T_NBL_PT=${T_NBL_PT:-N/A}
        SU_CELL=${SU_CELL:-N/A};           SU_CELL_MANHAT=${SU_CELL_MANHAT:-N/A}
        SU_NBL=${SU_NBL:-N/A};             SU_FULL_PT=${SU_FULL_PT:-N/A}
        SU_CELL_MT=${SU_CELL_MT:-N/A};     SU_MANH_PT=${SU_MANH_PT:-N/A}
        SU_NBL_PT=${SU_NBL_PT:-N/A}
        FASTEST_METHOD=${FASTEST_METHOD:-N/A}; FASTEST_TIME=${FASTEST_TIME:-N/A}

        echo "$N,$T,$T_FULL,$T_CELL,$T_CELL_MANHAT,$T_NBL,$T_FULL_PT,$T_CELL_MT,$T_MANH_PT,$T_NBL_PT,$SU_CELL,$SU_CELL_MANHAT,$SU_NBL,$SU_FULL_PT,$SU_CELL_MT,$SU_MANH_PT,$SU_NBL_PT,\"$FASTEST_METHOD\",$FASTEST_TIME" \
            >> "$CSV"

        echo "  → full=${T_FULL}s | cell_su=${SU_CELL}x | nbl_su=${SU_NBL}x | nbl_pt_su=${SU_NBL_PT}x | fastest=${FASTEST_METHOD}"
    done
done

echo ""
echo "============================================================"
echo "Sweep complete. Results: $CSV"
echo "Generating plots..."
echo "============================================================"

./venv/bin/python3 plot_results.py "$CSV" "$OUTPUT_DIR/plots"

echo ""
echo "Done. Plots saved to: $OUTPUT_DIR/plots/"
