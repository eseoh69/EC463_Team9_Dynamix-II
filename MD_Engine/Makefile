# ================================
#   Makefile for LJ MD Engine
#   With Pthreads Support (CPU only)
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude -pthread
LDFLAGS = -lm -lpthread

CXX     = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra

PYTHON  = python3     # for plot_lj.py

# Number of threads for pthread neighbor list (default: 4)
# Override with: make NTHREADS=8
NTHREADS ?= 4
CFLAGS += -DNBL_NUM_THREADS=$(NTHREADS)

SRC = \
    src/main.c \
    src/md.c \
    src/cell_list.c \
    src/neighbor_list.c \
    src/neighbor_list_pthread.c \
    src/md_io.c \
    src/pdb_importer.c \
    src/cell_mt.c \
    src/cell_manhattan_pthread.c \
    src/manhat.c

OBJ = $(SRC:.c=.o)

TARGET      = md_run
COMPARE_BIN = compare_pdb
COMPARE_SRC = compare.cpp

# -----------------------------------------
# Default build
# -----------------------------------------
all: output $(TARGET)
	@echo ""
	@echo "Build complete with $(NTHREADS) threads for pthread neighbor list"
	@echo "To change: make clean && make NTHREADS=8"
	@echo ""

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

# -----------------------------------------
# Build comparison binary
# -----------------------------------------
$(COMPARE_BIN): $(COMPARE_SRC)
	$(CXX) $(CXXFLAGS) $(COMPARE_SRC) -o $(COMPARE_BIN)

# -----------------------------------------
# Create output directory automatically
# -----------------------------------------
output:
	mkdir -p output

# -----------------------------------------
# Build object files
# -----------------------------------------
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

# -----------------------------------------
# Run simulation with default PDB
# -----------------------------------------
run: $(TARGET) output
	./md_run input/random_particles-1024.pdb

# -----------------------------------------
# Run comparison tests (writes to text file)
# -----------------------------------------
compare: $(COMPARE_BIN)
	@echo "Writing comparison report to output/compare_report.txt"
	$(eval INPUT_PDB := $(shell cat output/input_used.txt 2>/dev/null || echo "input/random_particles-1024.pdb"))
	@echo "Using input: $(INPUT_PDB)"
	@rm -f output/compare_report.txt

	@echo "========================================" >> output/compare_report.txt
	@echo "PDB COMPARISON REPORT" >> output/compare_report.txt
	@echo "Energy Conservation & Numerical Stability Test" >> output/compare_report.txt
	@echo "========================================" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "1) Compare all methods against each other (after autotuning steps)" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs CELL-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs FULL (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/full_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs CELL-LIST (half shell) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs CELL-LIST (Manhattan) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs CELL-LIST MANHATTAN (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/cell_manhat_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt


	@echo "FULL (Pthread) vs CELL-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_pthread_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL (pthread) vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_pthread_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FULL (pthread) vs CELL-LIST (half shell) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_pthread_positions.pdb output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt
	
	@echo "FULL (pthread) vs CELL-LIST (Manhattan) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_pthread_positions.pdb output/cell_manhat_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt	

	@echo "FULL (pthread) vs NEIGHBOR-LIST (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_pthread_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Serial) vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Serial) vs CELL-LIST (half-shell) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Serial) vs CELL-LIST (Manhattan) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Serial) vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "NEIGHBOR-LIST (Serial) vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/nbl_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "NEIGHBOR-LIST (Serial) vs CELL-LIST (half shell) (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/nbl_positions.pdb output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "NEIGHBOR-LIST (Serial) vs CELL-LIST (Manhattan) (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/nbl_positions.pdb output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (half shell) (pthread) vs CELL-LIST (Manhattan) (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_half_mt_positions.pdb output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (half shell) (pthread) vs NEIGHBOR-LIST (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_half_mt_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Manhattan) (pthread) vs NEIGHBOR-LIST (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/manhattan_pthread_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST (Manhattan) (pthread) vs CELL LIST MANHATTAN (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/manhattan_pthread_positions.pdb output/cell_manhat_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST MANHATTAN (Serial) vs FULL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_manhat_positions.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "CELL-LIST MANHATTAN (Serial) vs NEIGHBOR-LIST (pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_manhat_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt


	@echo "2) Final positions vs all autotuning outputs" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs FULL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs CELL-LIST:" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs FULL (PTHREAD):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/full_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt


	@echo "FINAL vs CELL LIST (HALFSHELL) (PTHREAD):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs CELL LIST MANHATTAN (PTHREAD):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "FINAL vs MANHATTAN (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/cell_manhat_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	


	@echo "3) Input vs all outputs (drift from initial)" >> output/compare_report.txt
	@echo "----------------------------------------" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs FULL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs CELL-LIST:" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs NEIGHBOR-LIST (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs FULL (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/full_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs CELL LIST (HALFSHELL) (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/cell_half_mt_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs CELL LIST MANHATTAN (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/manhattan_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs NEIGHBOR-LIST (Pthread):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/nbl_pthread_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs MANHATTAN (Serial):" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/cell_manhat_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "INPUT vs FINAL:" >> output/compare_report.txt
	@./$(COMPARE_BIN) $(INPUT_PDB) output/final_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt




	@echo "========================================" >> output/compare_report.txt
	@echo "Report complete." >> output/compare_report.txt
	@echo "Comparison report written to output/compare_report.txt"
	@cat output/compare_report.txt

# -----------------------------------------
# Plot energy graphs (Python script)
# -----------------------------------------
plot:
	$(PYTHON) test_scripts/plot_lj.py

# -----------------------------------------
# Basic clean
# -----------------------------------------
clean:
	rm -f $(OBJ) $(TARGET) $(COMPARE_BIN)
	rm -rf output
	rm -f energy_four_panel.png
	rm -f energy_final_logplot.png




# ── configurable variables (override on command line if needed) ───────────────
BINARY      ?= ./md_run
INPUT_DIR   ?= ./input
OUTPUT_DIR  ?= ./results
PYTHON      ?= ./venv/bin/python3
SWEEP_SCRIPT := run_sweep.sh
PLOT_SCRIPT  := plot_results.py
 
# ── sweep: run all simulations then plot ─────────────────────────────────────
.PHONY: sweep
sweep: $(BINARY) $(SWEEP_SCRIPT) $(PLOT_SCRIPT)
	@echo "Starting sweep: binary=$(BINARY) input=$(INPUT_DIR) output=$(OUTPUT_DIR)"
	@chmod +x $(SWEEP_SCRIPT)
	./$(SWEEP_SCRIPT) --binary $(BINARY) --input-dir $(INPUT_DIR) --output-dir $(OUTPUT_DIR)
 
# ── insights: replot from existing CSV without rerunning simulations ──────────
.PHONY: insights
insights: $(PLOT_SCRIPT) $(OUTPUT_DIR)/speedups.csv
	@echo "Regenerating plots from $(OUTPUT_DIR)/speedups.csv ..."
	$(PYTHON) $(PLOT_SCRIPT) $(OUTPUT_DIR)/speedups.csv $(OUTPUT_DIR)/plots
	@echo "Plots saved to $(OUTPUT_DIR)/plots/"
 
# ── sweep-clean: wipe results directory ──────────────────────────────────────
.PHONY: sweep-clean
sweep-clean:
	@echo "Removing $(OUTPUT_DIR)/ ..."
	rm -rf $(OUTPUT_DIR)
 
# ── guard: remind user if CSV is missing when running 'make insights' ─────────
$(OUTPUT_DIR)/speedups.csv:
	@echo "ERROR: $(OUTPUT_DIR)/speedups.csv not found. Run 'make sweep' first."
	@exit 1




# -----------------------------------------
# Full clean — also remove MD output files
# -----------------------------------------
distclean: clean
	rm -f output/full_positions.pdb
	rm -f output/cell_positions.pdb
	rm -f output/nbl_positions.pdb
	rm -f output/nbl_pthread_positions.pdb
	rm -f output/final_positions.pdb

	rm -f output/full_energies.csv
	rm -f output/cell_energies.csv
	rm -f output/nbl_energies.csv
	rm -f output/nbl_pthread_energies.csv
	rm -f output/final_energies.csv

	rm -f energy_final_logplot.png
	rm -f energy_four_panel.png

	rm -rf output/*

.PHONY: all run clean distclean output compare plot