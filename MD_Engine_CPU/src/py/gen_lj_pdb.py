#!/usr/bin/env python3
"""
Generate a BCC lattice PDB in LJ reduced units at target density.
BCC has 2 atoms per unit cell: N = 2*n^3.
Default: n=8 -> N=1024, rho*=0.8 -> L=(N/rho*)^(1/3) ~ 10.86
All particles start near the LJ minimum distance (stable initialization).
"""
import sys, math

n     = int(sys.argv[1])   if len(sys.argv) > 1 else 8      # cells per side
rho   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8
out   = sys.argv[3]        if len(sys.argv) > 3 else "input/lj_1024.pdb"

N = 2 * n * n * n
L = (N / rho) ** (1.0 / 3.0)
a = L / n   # unit cell size

print(f"BCC lattice: n={n}, N={N}, rho*={rho:.2f}, L={L:.4f}, a={a:.4f} LJ units -> {out}")

# BCC basis: corner (0,0,0) and body-center (0.5,0.5,0.5)*a
basis = [(0.0, 0.0, 0.0), (0.5 * a, 0.5 * a, 0.5 * a)]

atoms = []
for ix in range(n):
    for iy in range(n):
        for iz in range(n):
            for (bx, by, bz) in basis:
                x = ix * a + bx
                y = iy * a + by
                z = iz * a + bz
                atoms.append((x, y, z))

assert len(atoms) == N, f"Expected {N} atoms, got {len(atoms)}"

with open(out, "w") as f:
    # CRYST1 record tells MDAnalysis (and main.c) the true periodic box size
    f.write(f"CRYST1{L:9.3f}{L:9.3f}{L:9.3f}  90.00  90.00  90.00 P 1           1\n")
    for i, (x, y, z) in enumerate(atoms, start=1):
        f.write(f"ATOM  {i:5d}  X   XXX A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           X  \n")
    f.write("END\n")

print(f"Done. Nearest-neighbor distance = {a * math.sqrt(3)/2:.4f} sigma (LJ min = {2**(1/6):.4f})")
