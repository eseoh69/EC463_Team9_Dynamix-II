#!/usr/bin/env bash
# =============================================================================
# run_sweep_gpu.sh
#
# Sweeps the GPU MD engine over particle counts:
#   1024  2048  4096  8192  16384  32768
#
# 4 GPU algorithms timed per run (AUTOTUNE_N_TIMESTEPS each):
#   Full O(N^2), Cell-List, Neighbor-List, Manhattan
#
# Usage:
#   chmod +x run_sweep_gpu.sh
#   ./run_sweep_gpu.sh [--binary PATH] [--input-dir PATH] [--output-dir PATH]
#
# Defaults:
#   --binary     ./md_run
#   --input-dir  ./input
#   --output-dir ./results
# =============================================================================

set -uo pipefail

BINARY="./md_run"
INPUT_DIR="./input"
OUTPUT_DIR="./results"

PARTICLE_COUNTS=(1024 2048 4096 8192 16384 32768)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --binary)     BINARY="$2";     shift 2 ;;
        --input-dir)  INPUT_DIR="$2";  shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: binary not found or not executable: $BINARY"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/plots"
mkdir -p output   # md_run writes energy CSVs and PDBs here

CSV="$OUTPUT_DIR/timings.csv"

cat > "$CSV" << 'EOF'
N,time_full,time_cell,time_nbl,time_man,su_cell,su_nbl,su_man,fastest_method,fastest_time
EOF

# ── helpers ───────────────────────────────────────────────────────────────────
parse_time() {
    # $1 = log file, $2 = fixed-string pattern
    grep -m1 -F "$2" "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | head -1
}

parse_speedup() {
    # Second float on the matching line (the Yx value)
    grep -m1 -F "$2" "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | sed -n '2p'
}

parse_fastest_method() {
    grep -m1 "^FASTEST:" "$1" \
        | sed 's/^FASTEST: //'
}

parse_fastest_time() {
    grep -m1 -F "Full run done." "$1" \
        | grep -oE '[0-9]+\.[0-9]+' \
        | head -1
}

# ── main sweep ────────────────────────────────────────────────────────────────
TOTAL=${#PARTICLE_COUNTS[@]}
RUN=0

for N in "${PARTICLE_COUNTS[@]}"; do
    PDB="$INPUT_DIR/random_particles-${N}.pdb"
    if [[ ! -f "$PDB" ]]; then
        echo "WARNING: input file not found, skipping: $PDB"
        continue
    fi

    RUN=$(( RUN + 1 ))
    LOG="$OUTPUT_DIR/logs/N${N}.log"

    echo "──────────────────────────────────────────────────────"
    echo "Run ${RUN}/${TOTAL}  |  N=${N}"
    echo "  binary : $BINARY"
    echo "  input  : $PDB"
    echo "  log    : $LOG"
    echo "──────────────────────────────────────────────────────"

    if "$BINARY" "$PDB" 2>&1 | tee "$LOG"; then
        T_FULL=$(parse_time "$LOG" "Full O(N^2) done.")
        T_CELL=$(parse_time "$LOG" "Cell-list done.")
        T_NBL=$( parse_time "$LOG" "Neighbor-list done.")
        T_MAN=$(  parse_time "$LOG" "Manhattan done.")

        SU_CELL=$(parse_speedup "$LOG" "Cell-list:     ")
        SU_NBL=$( parse_speedup "$LOG" "Neighbor-list: ")
        SU_MAN=$(  parse_speedup "$LOG" "Manhattan:     ")

        FASTEST_METHOD=$(parse_fastest_method "$LOG")
        FASTEST_TIME=$(  parse_fastest_time   "$LOG")

        # Preserve per-run output PDBs
        for f in full_positions cell_positions nbl_positions manhattan_positions final_positions; do
            [[ -f "output/${f}.pdb" ]] && \
                cp "output/${f}.pdb" "$OUTPUT_DIR/logs/N${N}_${f}.pdb"
        done
    else
        echo "  WARNING: binary exited with error for N=${N}, recording N/A"
        T_FULL="N/A"; T_CELL="N/A"; T_NBL="N/A"; T_MAN="N/A"
        SU_CELL="N/A"; SU_NBL="N/A"; SU_MAN="N/A"
        FASTEST_METHOD="N/A"; FASTEST_TIME="N/A"
    fi

    T_FULL=${T_FULL:-N/A};  T_CELL=${T_CELL:-N/A}
    T_NBL=${T_NBL:-N/A};    T_MAN=${T_MAN:-N/A}
    SU_CELL=${SU_CELL:-N/A}; SU_NBL=${SU_NBL:-N/A}; SU_MAN=${SU_MAN:-N/A}
    FASTEST_METHOD=${FASTEST_METHOD:-N/A}; FASTEST_TIME=${FASTEST_TIME:-N/A}

    echo "$N,$T_FULL,$T_CELL,$T_NBL,$T_MAN,$SU_CELL,$SU_NBL,$SU_MAN,\"$FASTEST_METHOD\",$FASTEST_TIME" \
        >> "$CSV"

    echo "  → full=${T_FULL}s | su_cell=${SU_CELL}x | su_nbl=${SU_NBL}x | su_man=${SU_MAN}x | fastest=${FASTEST_METHOD}"
done

echo ""
echo "============================================================"
echo "Sweep complete. Results: $CSV"
echo "Generating plots..."
echo "============================================================"

python3 test/plot_results_gpu.py "$CSV" "$OUTPUT_DIR/plots"

echo ""
echo "Done. Plots saved to: $OUTPUT_DIR/plots/"
