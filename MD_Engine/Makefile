# ================================
#   Makefile for LJ MD Engine
# ================================

CC      = gcc
CFLAGS  = -O2 -std=c11 -Wall -Wextra -Iinclude
LDFLAGS = -lm

CXX     = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra

PYTHON  = /opt/homebrew/bin/python3     # for plot_lj.py

SRC = \
    src/main.c \
    src/md.c \
    src/cell_list.c \
    src/neighbor_list.c \
    src/md_io.c \
    src/pdb_importer.c \
	src/cell_mt.c

OBJ = $(SRC:.c=.o)

TARGET      = md_run
COMPARE_BIN = compare_pdb
COMPARE_SRC = compare.cpp

# -----------------------------------------
# Default build
# -----------------------------------------
all: output $(TARGET)

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
	./md_run input.pdb

# -----------------------------------------
# Run comparison tests (writes to text file)
# -----------------------------------------
compare: $(COMPARE_BIN)
	@echo "Writing comparison report to output/compare_report.txt"
	@rm -f output/compare_report.txt

	@echo "Running PDB comparison tests" >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "1) full_positions, cell_positions, nbl_positions vs each other" >> output/compare_report.txt
	@echo "Comparing output/full_positions.pdb vs output/cell_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing output/full_positions.pdb vs output/nbl_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/full_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing output/cell_positions.pdb vs output/nbl_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/cell_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "2) final_positions vs all other outputs" >> output/compare_report.txt
	@echo "Comparing output/final_positions.pdb vs output/full_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing output/final_positions.pdb vs output/cell_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing output/final_positions.pdb vs output/nbl_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) output/final_positions.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "3) input.pdb vs all outputs" >> output/compare_report.txt
	@echo "Comparing input.pdb vs output/full_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/full_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing input.pdb vs output/cell_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/cell_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing input.pdb vs output/nbl_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/nbl_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

	@echo "Comparing input.pdb vs output/final_positions.pdb" >> output/compare_report.txt
	@./$(COMPARE_BIN) input.pdb output/final_positions.pdb >> output/compare_report.txt
	@echo "" >> output/compare_report.txt

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

# -----------------------------------------
# Full clean — also remove MD output files
# -----------------------------------------
distclean: clean
	rm -f output/full_positions.pdb
	rm -f output/cell_positions.pdb
	rm -f output/nbl_positions.pdb
	rm -f output/final_positions.pdb

	rm -f output/full_energies.csv
	rm -f output/cell_energies.csv
	rm -f output/nbl_energies.csv
	rm -f output/final_energies.csv

	rm -f energy_final_logplot.png
	rm -f energy_four_panel.png

	rm -rf output/*

.PHONY: all run clean distclean output compare plot
