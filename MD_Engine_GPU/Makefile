# ================================
#   Makefile for MD Engine GPU
# ================================

NVCC      = nvcc
NVCCFLAGS = -O3 -arch=sm_86 -Iinclude -rdc=true

TARGET    = md_run
DLINK_OBJ = device_link.o

SRC = src/main.cu src/cuda_cell_list.cu src/cuda_nbl.cu src/cuda_full.cu src/md_io.c
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

.PHONY: all run clean output
