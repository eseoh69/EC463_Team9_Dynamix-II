# Team_09 Dynamix II

<p align="center">
<img src="./images/team_9_4-9.png" width="50%">
</p>
<p align="center">
Team Photo
</p>

## Blurb

Our product is a molecular dynamics simulator that uses runtime autotuning to empirically determine the fastest MD algorithm for a user’s platform. Our product simulates the movement of particles over time like other MD simulators (LAMMPS, GROMACS [1]) by computing 
The Lennard-Jones force between particles. This is the only force we compute at this time because it is the most computationally expensive and we focused on its optimization. For a more realistic simulation, other forces would have to be calculated as well.

What makes our product unique is that we use autotuning on the backend to choose the fastest simulation algorithm for a user’s platform. There are multiple MD algorithms that can be used for simulation and they each have tradeoffs in terms of speed, scalability and memory. The specs of a user’s hardware determine which of these tradeoffs are worth making to achieve maximum performance. For example, one algorithm may be the fastest on one CPU core but get overtaken by a better scaling algorithm running on multiple cores. We follow an autotuning approach [2] where we evaluate MD implementations on a user’s input over a short period of time to identify the fastest algorithm.

## Team links
- [Team Google Drive](https://drive.google.com/drive/folders/1uDh5uxmoZppGSzFuEcW462L6dL4mUN0_?usp=drive_link)

## Course links
- [ECE Senior Design Piazza Site](https://piazza.com/bu/fall2025/ec463/home)
- [Blackboard](http://learn.bu.edu/)

## Project structure
Below is a high-level file tree for the repository:

- `README.md`
- `cpu/`
  - `Makefile` - CPU build and run targets
  - `md_run` - CPU autotuning executable output after build
  - `src/` - CPU source files and algorithm implementations
  - `include/` - CPU headers
  - `run_sweep.sh` - CPU sweep script
- `gpu/`
  - `Makefile` - GPU build and run targets
  - `md_run` - GPU autotuning executable output after build
  - `src/` - GPU source files and CUDA kernels
  - `include/` - GPU headers
  - `run_sweep_gpu.sh` - GPU sweep script
- `scripts/` scripts for output validation and performance evaluation
  - `plot_results_cpu.py`, `plot_results_gpu.py`, plot speedup charts for cpu and gpu
   `plot_energy.py`, plot kinetic, potential, total energy charts over course of simulation
  - `compare_pdbs.py` compare outputs
  - `compare_pdbs_all.py` compare outputs, more features
  - ...
- `images/` - 
- `input/` - shared sample input files
- `results/` - our results from our poster

## Setup and build
1. Clone the repository:

```bash
git clone https://github.com/<your-org>/EC463_Team9_Dynamix-II.git
cd EC463_Team9_Dynamix-II
```

2. Install dependencies

The MDAnalyis Python package is a dependency for CPU and GPU autotuners.
Install your package to your preferred location (system wide install, venv, conda env)
```
pip install MDAnalysis==2.10.0 # note: MDAnalysis 2.10.0 needs Python 3.11.0 or later
```
Using a text editor, set the ```PYTHON``` variable in ```cpu/Makefile``` and ```gpu/Makefile``` to the python path where MDAnalysis is installed.
By default, it is set to ```python3```. If using a venv, your path would be ```<abs path to venv>/bin/python3```.


3. Build the CPU version:

```bash
cd cpu
make
```

Optional: change the number of pthread threads for the CPU build:

```bash
make clean
make NTHREADS=8
```

4. Build the GPU version (if you have CUDA support):

```bash
cd ../gpu
make
```

If your GPU architecture differs from the default `sm_86`, update `NVCCFLAGS` in `gpu/Makefile` before building.

## Usage examples
### CPU run
Run the CPU binary from the `cpu/` directory:

```bash
cd cpu
./md_run ../input/random_particles-1024.pdb
```

From the repository root, you can also run:

```bash
cpu/md_run ./input/random_particles-1024.pdb
```

outputs will be written to whatever directory the command is run from.

### GPU run
Run the GPU binary from the `gpu/` directory:

```bash
cd gpu
./md_run input/random_particles-1024.pdb
```

From the repository root, you can also run:

```bash
gpu/md_run ../input/random_particles-1024.pdb
```

### Sweep and benchmark scripts
run_sweep, will run all algorithms at various sizes and record results to ```results/``` in your current directory. 
These results can be plotted with: ```scripts/plot_results_cpy.py```, ```scripts/plot_results_gpu.py```.

CPU sweep:

```bash
cd cpu
./run_sweep.sh
```

GPU sweep:

```bash
cd gpu
./run_sweep_gpu.sh
```

### Plotting and analysis
The repository includes plotting scripts in `scripts/`.

For CPU performance graphs:

```bash
python3 scripts/plot_results_cpu.py
```

For GPU performance graphs:

```bash
python3 scripts/plot_results_gpu.py
```

For conservation of energy plots:

```bash
python3 scripts/plot_energy.py
```

Output validation. By default checks CPU outputs only, can be changed in the file.
```bash
python3 scripts/compare_all.py
```

Notes:
- Output files are written to your current directory's `output/` and `results/` directories.
