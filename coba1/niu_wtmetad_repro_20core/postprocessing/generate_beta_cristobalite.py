#!/usr/bin/env python3
"""Generate idealized beta-cristobalite SiO2 supercell for Niu-style tests.
Default: 4x4x4 conventional cells = 512 Si + 1024 O = 1536 atoms.
Si atoms form a diamond lattice. O atoms are placed at unique Si-Si bond midpoints.
This is an idealized structural target for XRD and WTMETAD testing, not a relaxed DFT structure.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np

SI_MASS = 28.0855
O_MASS = 15.9994

def min_image(delta: np.ndarray, box: float) -> np.ndarray:
    return delta - box * np.round(delta / box)

def write_lammps_data(path: Path, atoms, box: float):
    # atoms: list of (type, charge, x,y,z, comment)
    with path.open('w') as f:
        f.write('# idealized beta-cristobalite SiO2; atom_style charge\n\n')
        f.write(f'{len(atoms)} atoms\n')
        f.write('2 atom types\n\n')
        f.write(f'0.0 {box:.10f} xlo xhi\n0.0 {box:.10f} ylo yhi\n0.0 {box:.10f} zlo zhi\n\n')
        f.write('Masses\n\n')
        f.write(f'1 {SI_MASS:.6f} # Si\n')
        f.write(f'2 {O_MASS:.6f} # O\n\n')
        f.write('Atoms # charge\n\n')
        for i,(typ,q,x,y,z,c) in enumerate(atoms, start=1):
            f.write(f'{i:8d} {typ:2d} {q: .8f} {x:15.8f} {y:15.8f} {z:15.8f} # {c}\n')

def write_xyz(path: Path, atoms, box: float):
    sym = {1:'Si',2:'O'}
    with path.open('w') as f:
        f.write(f'{len(atoms)}\n')
        f.write(f'Lattice="{box} 0 0 0 {box} 0 0 0 {box}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for typ,q,x,y,z,c in atoms:
            f.write(f'{sym[typ]} {x:.8f} {y:.8f} {z:.8f}\n')

def generate(a=7.15, nrep=4):
    # conventional diamond basis: 8 Si / cell
    fcc = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]], dtype=float)
    basis = np.vstack([fcc, fcc + 0.25]) % 1.0
    box = a*nrep
    si_pos = []
    for i in range(nrep):
      for j in range(nrep):
       for k in range(nrep):
        shift = np.array([i,j,k], float)
        for b in basis:
            si_pos.append((shift+b)*a)
    si_pos = np.array(si_pos)
    # Find unique nearest-neighbor Si-Si bonds in periodic cell.
    target = a*math.sqrt(3)/4
    bonds = []
    nsi = len(si_pos)
    for i in range(nsi):
        for j in range(i+1, nsi):
            dvec = min_image(si_pos[j]-si_pos[i], box)
            d = np.linalg.norm(dvec)
            if abs(d-target) < 0.05:
                mid = (si_pos[i] + 0.5*dvec) % box
                bonds.append(mid)
    # In diamond, unique bonds should be 2*Nsi = 1024 for Nsi=512.
    atoms = []
    for p in si_pos:
        atoms.append((1,0.0,float(p[0]),float(p[1]),float(p[2]),'Si'))
    for p in bonds:
        atoms.append((2,0.0,float(p[0]),float(p[1]),float(p[2]),'O'))
    return atoms, box, len(si_pos), len(bonds)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--a', type=float, default=7.15, help='beta-cristobalite conventional-cell parameter in Angstrom')
    ap.add_argument('--nrep', type=int, default=4, help='replication per direction; 4 gives 1536 atoms')
    ap.add_argument('--outdir', default='../structures')
    args = ap.parse_args()
    atoms, box, nsi, no = generate(args.a, args.nrep)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    stem = f'beta_cristobalite_a{args.a:g}_N{len(atoms)}'
    write_lammps_data(outdir / f'{stem}.data', atoms, box)
    write_xyz(outdir / f'{stem}.xyz', atoms, box)
    print(f'Wrote {stem}: total={len(atoms)}, Si={nsi}, O={no}, box={box:.4f} A')

if __name__ == '__main__':
    main()
