# Usage

## 1. Setup the venv
Setup Python path in Makefile

We use MD_Engine uses the MD_Analysis python package to parse input. We use a venv to contain dependencies This only has to be done once per install.
Python 3.12.4

```
cd MD_Engine
python3 -m venv venv
source venv/bin/activate
pip install MDAnalysis==2.10.0
deactivate
```

## 2. Example run
```
make
./md_run input/2FBD.pdb
./md_run input/random_particles-1024.pdb
```