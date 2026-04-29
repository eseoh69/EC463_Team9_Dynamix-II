# Dynamix-II — GPU MD Engine Usage Guide

## 1. Connect to SCC

```bash
ssh username@scc1.bu.edu
```

## 2. Clone the Repository

```bash
git clone https://github.com/eseoh69/EC463_Team9_Dynamix-II.git
cd EC463_Team9_Dynamix-II
```

## 3. Checkout the GPU Branch

```bash
git checkout md_engine_gpu
cd md_engine_gpu
```

## 4. Request a GPU Node

```bash
qrsh -l gpus=1 -l gpu_type=A40 -P ece601
```

> After the session starts, navigate back to the `MD_Engine_GPU` directory.

## 5. Set Up Python Environment

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install matplotlib MDAnalysis
```

## 6. Load CUDA

```bash
module load cuda/12.8
```

## 7. Build

```bash
make clean
make
```

The build runs with `-j4` parallel compilation. The `output/` directory is created automatically if it doesn't exist.

## 8. Run

Runs the simulation on a single input file. You can override the input file on the command line:

```bash
make run
```

```bash
make run INPUT_FILE=input/random_particles-4096.pdb
```

The default input is `input/random_particles-1024.pdb`. Output files are written to `output/`:
- `full_energies.csv`, `cell_energies.csv`, `nbl_energies.csv`, `manhattan_cell_energies.csv`, `final_energies.csv`
- `full_positions.pdb`, `cell_positions.pdb`, `nbl_positions.pdb`, `final_positions.pdb`

## 9. Visualize Energy

Builds and runs the simulation (if not already done), then plots kinetic, potential, and total energy over time for each algorithm. Plots and a drift summary are saved to `output/`.

```bash
make plot
```

You can override the input file the same way as `make run`:

```bash
make plot INPUT_FILE=input/random_particles-4096.pdb
```

Output saved to `output/`:
- `energy_<method>.png` — per-algorithm energy plots
- `energy_summary.txt` — energy drift summary across all methods

## 10. Compare Final Positions

Builds and runs the simulation (if not already done), then computes per-particle positional deviation between all algorithm pairs.

```bash
make compare
```

Output saved to `output/`:
- `pdb_comparison.txt` — mean/max/std deviation between each pair of methods

## 11. Run GPU Sweep Benchmark

Runs the GPU MD engine across all input particle counts in the `input/` directory, saves timing results to `results/`, and generates plots.

```bash
make sweep
```

You can override the defaults on the command line:

```bash
make sweep BINARY=./md_run INPUT_DIR=./input OUTPUT_DIR=./results
```

| Variable | Default | Description |
|---|---|---|
| `BINARY` | `./md_run` | Path to the compiled binary |
| `INPUT_DIR` | `./input` | Directory containing `.pdb` input files |
| `OUTPUT_DIR` | `./results` | Directory where timing CSV and plots are saved |
| `PYTHON` | `python3` | Python interpreter used for plotting |

## 12. Regenerate Plots from Existing Results

If you've already run a sweep and only want to regenerate the plots without re-running the simulation:

```bash
make insights
```

This reads from `results/timings.csv` and writes updated plots to `results/plots/`. Requires a prior `make sweep` to have produced the CSV.

## 13. Clean Sweep Results

To delete the `results/` directory:

```bash
make sweep-clean
```

## 14. Clean Build Artifacts

To remove compiled objects and the binary (does not affect `results/`):

```bash
make clean
```