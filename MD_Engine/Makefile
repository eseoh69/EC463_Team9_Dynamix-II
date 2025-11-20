# ================================
#   Makefile for LJ MD Engine
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude
LDFLAGS = -lm

PYTHON  = /opt/homebrew/bin/python3    # <--- ADDED FOR plot_lj.py

SRC = \
    src/main.c \
    src/md.c \
    src/cell_list.c \
    src/neighbor_list.c \
    src/md_io.c \
    src/pdb_importer.c

OBJ = $(SRC:.c=.o)

TARGET = md_run

# -----------------------------------------
# Default build
# -----------------------------------------
all: output $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

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
run: $(TARGET)
	./md_run input.pdb

# -----------------------------------------
# NEW: Generate energy plots (4-panel)
# -----------------------------------------
plot:
	$(PYTHON) test_scripts/plot_lj.py

# -----------------------------------------
# Basic clean (build artifacts only)
# -----------------------------------------
clean:
	rm -f $(OBJ) $(TARGET)
	rm -rf output
	rm -f energy_four_panel.png
	rm -f energy_final_logplot.png

# -----------------------------------------
# Full clean: remove MD outputs
# -----------------------------------------
distclean: clean
	# Remove MD output PDBs ONLY
	rm -f output/full_positions.pdb
	rm -f output/cell_positions.pdb
	rm -f output/nbl_positions.pdb
	rm -f output/final_positions.pdb

	# Remove MD output CSVs ONLY
	rm -f output/full_energies.csv
	rm -f output/cell_energies.csv
	rm -f output/nbl_energies.csv
	rm -f output/final_energies.csv

	rm -f energy_final_logplot.png
	rm -f energy_four_panel.png

	# Remove everything inside output/
	rm -rf output/*

.PHONY: all run clean distclean output plot
