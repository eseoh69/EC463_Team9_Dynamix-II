# ================================
#   Makefile for LJ MD Engine
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude
LDFLAGS = -lm

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
# Basic clean (build artifacts only)
# -----------------------------------------
clean:
	rm -f $(OBJ) $(TARGET)
	rm -rf output

# -----------------------------------------
# Full clean: remove old stray outputs too
# -----------------------------------------
distclean: clean
	# Remove MD output PDBs ONLY
	rm -f full_positions.pdb cell_positions.pdb nbl_positions.pdb
	rm -f final_positions.pdb

	# Remove MD output CSVs ONLY
	rm -f full_energies.csv cell_energies.csv nbl_energies.csv

	# Remove anything inside output/ folder
	rm -rf output/*

	# DO NOT delete input.pdb or any user-created PDBs
	# DO NOT delete *.pdb in project root

.PHONY: all run clean distclean output
