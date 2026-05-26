#!/usr/bin/env python3
"""ASE version of the beta-cristobalite generator.

Requires:
    pip install ase

The dependency-free generator `make_beta_cristobalite.py` is recommended if ASE
is unavailable. This script uses ASE to write XYZ/PDB and then calls the same
LAMMPS-data writer logic for atom_style charge.
"""
from __future__ import annotations

import argparse
from pathlib import Path

try:
    from ase import Atoms
    from ase.io import write
except ImportError as exc:
    raise SystemExit("ASE is not installed. Install with `pip install ase`, or use scripts/make_beta_cristobalite.py instead.") from exc

# Reuse the dependency-free implementation for robust geometry and LAMMPS data output.
import make_beta_cristobalite as base


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=7.15)
    ap.add_argument("--reps", type=int, nargs=3, default=(2, 2, 2))
    ap.add_argument("--jitter", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out-prefix", default="structures/beta_cristobalite_192_ase")
    args = ap.parse_args()

    symbols, positions, box = base.generate(args.a, tuple(args.reps), args.jitter, args.seed)
    atoms = Atoms(symbols=symbols, positions=positions, cell=box, pbc=True)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    write(prefix.with_suffix(".xyz"), atoms)
    write(prefix.with_suffix(".pdb"), atoms)
    base.write_lammps_data(prefix.with_suffix(".data"), symbols, positions, box)
    print(f"Wrote ASE structure to {prefix.with_suffix('.xyz')} and LAMMPS data to {prefix.with_suffix('.data')}")


if __name__ == "__main__":
    main()
