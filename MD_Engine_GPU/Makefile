# ================================
#   Makefile for MD Engine GPU
# ================================

NVCC      = nvcc
NVCCFLAGS = -O3 -arch=sm_86 -Iinclude -rdc=true

TARGET    = md_run
DLINK_OBJ = device_link.o

SRC = src/main.cu src/cuda_cell_list.cu src/cuda_nbl.cu src/cuda_full.cu \
      src/gpu_manhattan_cell_cap.cu src/cell_list.c src/md_io.c
OBJ = $(patsubst src/%.cu,src/%.o,$(patsubst src/%.c,src/%.o,$(SRC)))

# -----------------------------------------
# Default build
# -----------------------------------------
MAKEFLAGS += -j4

all: output $(TARGET)

$(TARGET): $(OBJ) $(DLINK_OBJ)
	$(NVCC) -arch=sm_86 $(OBJ) $(DLINK_OBJ) -o $(TARGET) -lm

# -----------------------------------------
# Device link step
# -----------------------------------------
$(DLINK_OBJ): $(OBJ)
	$(NVCC) -arch=sm_86 -dlink $(OBJ) -o $(DLINK_OBJ)

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
run: $(TARGET) output

# -----------------------------------------
# Clean
# -----------------------------------------
clean:
	rm -f $(OBJ) $(DLINK_OBJ) $(TARGET)
	rm -rf output

# -----------------------------------------
# Sweep & plotting
# -----------------------------------------
# Override on command line: make sweep INPUT_DIR=./myinput
BINARY      ?= ./md_run
INPUT_DIR   ?= ./input
OUTPUT_DIR  ?= ./results
PYTHON      ?= python3
SWEEP_SCRIPT := run_sweep_gpu.sh
PLOT_SCRIPT  := plot_results_gpu.py

# sweep: build, run all particle counts, then plot
.PHONY: sweep
sweep: $(TARGET) $(SWEEP_SCRIPT) $(PLOT_SCRIPT)
	@echo "Starting GPU sweep: binary=$(BINARY) input=$(INPUT_DIR) output=$(OUTPUT_DIR)"
	@chmod +x $(SWEEP_SCRIPT)
	./$(SWEEP_SCRIPT) \
	    --binary     $(BINARY)     \
	    --input-dir  $(INPUT_DIR)  \
	    --output-dir $(OUTPUT_DIR)

# insights: replot from existing CSV without rerunning
.PHONY: insights
insights: $(PLOT_SCRIPT) $(OUTPUT_DIR)/timings.csv
	@echo "Regenerating plots from $(OUTPUT_DIR)/timings.csv ..."
	$(PYTHON) $(PLOT_SCRIPT) $(OUTPUT_DIR)/timings.csv $(OUTPUT_DIR)/plots
	@echo "Plots saved to $(OUTPUT_DIR)/plots/"

# sweep-clean: wipe results directory
.PHONY: sweep-clean
sweep-clean:
	@echo "Removing $(OUTPUT_DIR)/ ..."
	rm -rf $(OUTPUT_DIR)

# guard: error if CSV missing when running 'make insights'
$(OUTPUT_DIR)/timings.csv:
	@echo "ERROR: $(OUTPUT_DIR)/timings.csv not found. Run 'make sweep' first."
	@exit 1

.PHONY: all run clean output
