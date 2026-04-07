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
cd MD_Engine_GPU
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

## 8. Run

```bash
./md_run input/random_particles-1024.pdb
```

Output files are written to `output/`:
- `full_energies.csv`, `cell_energies.csv`, `nbl_energies.csv`, `final_energies.csv`
- `full_positions.pdb`, `cell_positions.pdb`, `nbl_positions.pdb`, `final_positions.pdb`

## 9. Visualize Energy

```bash
python3 plot_energy.py
```

## 10. Compare Final Positions

```bash
python3 compare_pdbs.py
```
