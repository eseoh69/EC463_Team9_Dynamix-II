# ================================
#   Makefile.gpu — GPU-only build
#   (CUDA Manhattan runner)
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude

NVCC    = nvcc
NVFLAGS = -O3 -std=c++17 -Iinclude

# Use NVCC to link so cudart gets pulled in cleanly
TARGET  = md_run_gpu

# --- Minimal CPU sources needed by main_gpu.c ---
SRC_C = \
    src/main_gpu_manhattan.cu \
    src/md_io.c \
    src/pdb_importer.c \
    src/cell_list.c

# --- CUDA source ---
SRC_CU = src/gpu_manhattan_cell_cap.cu

OBJ_C  = $(SRC_C:.c=.o)
OBJ_CU = $(SRC_CU:.cu=.o)

all: output $(TARGET)

output:
	mkdir -p output

# Compile C
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

# Compile CUDA
src/%.o: src/%.cu include/*.h
	$(NVCC) $(NVFLAGS) -c $< -o $@

# Link (NVCC)
$(TARGET): $(OBJ_C) $(OBJ_CU)
	$(NVCC) -o $(TARGET) $(OBJ_C) $(OBJ_CU) -lm

clean:
	rm -f $(OBJ_C) $(OBJ_CU) $(TARGET)

.PHONY: all clean output