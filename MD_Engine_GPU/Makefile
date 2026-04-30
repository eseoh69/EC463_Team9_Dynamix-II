# ================================
#   Makefile for MD Engine GPU
# ================================
NVCC      = nvcc
NVCCFLAGS = -O3 -arch=sm_86 -Iinclude -rdc=true
TARGET    = md_run
DLINK_OBJ = device_link.o

CU_SRC  = src/main.cu src/cuda_cell_list.cu src/cuda_nbl.cu src/cuda_full.cu \
           src/gpu_manhattan_cell_cap.cu
C_SRC   = src/md_io.c src/cell_list.c

CU_OBJ  = $(CU_SRC:src/%.cu=src/%.o)
C_OBJ   = $(C_SRC:src/%.c=src/%.o)
OBJ     = $(CU_OBJ) $(C_OBJ)

# -----------------------------------------
# Default build (no parallel — avoids rdc race)
# -----------------------------------------
.NOTPARALLEL:

all: output $(TARGET)

$(TARGET): $(OBJ) $(DLINK_OBJ)
	$(NVCC) -arch=sm_86 $(OBJ) $(DLINK_OBJ) -o $(TARGET) -lm

# -----------------------------------------
# Device link step (CUDA objects only, not C)
# -----------------------------------------
$(DLINK_OBJ): $(CU_OBJ)
	$(NVCC) -arch=sm_86 -dlink $(CU_OBJ) -o $(DLINK_OBJ)

# -----------------------------------------
# Compile CUDA sources
# -----------------------------------------
src/%.o: src/%.cu include/*.h include/*.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# -----------------------------------------
# Compile C sources through nvcc
# -----------------------------------------
src/%.o: src/%.c include/*.h include/*.cuh
	$(NVCC) $(NVCCFLAGS) -x c -c $< -o $@

# -----------------------------------------
# Create output directory automatically
# -----------------------------------------
output:
	mkdir -p output

# -----------------------------------------
# Run simulation
# -----------------------------------------
BINARY     ?= ./md_run
INPUT_FILE ?= input/random_particles-1024.pdb
INPUT_DIR  ?= ./input
OUTPUT_DIR ?= ./results
PYTHON     ?= python3

run: $(TARGET) output
	$(BINARY) $(INPUT_FILE)

# -----------------------------------------
# Plot energy from output/ CSVs
# -----------------------------------------
PLOT_ENERGY_SCRIPT := test/plot_energy.py

.PHONY: plot
plot:
	@echo "Plotting energy from output/ ..."
	$(PYTHON) $(PLOT_ENERGY_SCRIPT)
	@echo "Energy plots saved to output/"

# -----------------------------------------
# Compare final PDB positions in output/
# -----------------------------------------
COMPARE_SCRIPT := test/compare_pdbs.py

.PHONY: compare
compare:
	@echo "Comparing PDB positions in output/ ..."
	$(PYTHON) $(COMPARE_SCRIPT)
	@echo "Comparison saved to output/pdb_comparison.txt"

# -----------------------------------------
# Clean
# -----------------------------------------
clean:
	rm -f $(OBJ) $(DLINK_OBJ) $(TARGET)
	rm -rf output

# -----------------------------------------
# Sweep & plotting
# -----------------------------------------
SWEEP_SCRIPT := run_sweep_gpu.sh
PLOT_SCRIPT  := test/plot_results_gpu.py

.PHONY: sweep
sweep: $(TARGET) $(SWEEP_SCRIPT) $(PLOT_SCRIPT)
	@echo "Starting GPU sweep: binary=$(BINARY) input=$(INPUT_DIR) output=$(OUTPUT_DIR)"
	@chmod +x $(SWEEP_SCRIPT)
	./$(SWEEP_SCRIPT) \
	    --binary     $(BINARY)     \
	    --input-dir  $(INPUT_DIR)  \
	    --output-dir $(OUTPUT_DIR)

# -----------------------------------------
# Insights: replot from existing CSV
# -----------------------------------------
.PHONY: insights
insights: $(PLOT_SCRIPT) $(OUTPUT_DIR)/timings.csv
	@echo "Regenerating plots from $(OUTPUT_DIR)/timings.csv ..."
	$(PYTHON) $(PLOT_SCRIPT) $(OUTPUT_DIR)/timings.csv $(OUTPUT_DIR)/plots
	@echo "Plots saved to $(OUTPUT_DIR)/plots/"

# -----------------------------------------
# Sweep clean
# -----------------------------------------
.PHONY: sweep-clean
sweep-clean:
	@echo "Removing $(OUTPUT_DIR)/ ..."
	rm -rf $(OUTPUT_DIR)

# guard: error if CSV missing when running 'make insights'
$(OUTPUT_DIR)/timings.csv:
	@echo "ERROR: $(OUTPUT_DIR)/timings.csv not found. Run 'make sweep' first."
	@exit 1

.PHONY: all run clean output plot compare