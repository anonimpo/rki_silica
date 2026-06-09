#!/usr/bin/env python3
"""
run_md_chgnet.py

MD runner for FTIR sampling.

Purpose
- Run ASE MD with CHGNet forces.
- Save only a small dipole time series for FTIR, plus optional sparse trajectory.
- Designed to scale from a few hundred atoms to about 1500 atoms without writing huge files.

Important
- CHGNet provides energy, forces, stress, magnetic moments, and site energies through ASE.
- It does not provide IR dipoles or dynamic charges.
- This script computes an approximate dipole M(t) = sum_i q_i r_i(t) using fixed charges.
- Use charges from the LAMMPS data file, an external charge file, or a manual charge map.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np


def parse_key_value_map(text: Optional[str], key_type=str, value_type=float) -> Dict:
    """Parse maps like 'Si:2.4,O:-1.2,H:0.6' or '1:Si,2:O'."""
    if text is None or text.strip() == "":
        return {}
    out = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Map item must contain ':', got {item!r}")
        k, v = item.split(":", 1)
        out[key_type(k.strip())] = value_type(v.strip())
    return out


def load_charge_file(path: str | Path, n_atoms: int) -> np.ndarray:
    """Load charges from a text file.

    Accepted formats:
    - one column: q, already in atom order
    - two or more columns: id q ...; rows are sorted by id
    """
    arr = np.loadtxt(path, comments="#")
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        if arr.size != n_atoms:
            raise ValueError(f"Charge file has {arr.size} values, expected {n_atoms}")
        return arr.copy()
    if arr.shape[0] != n_atoms:
        raise ValueError(f"Charge file has {arr.shape[0]} rows, expected {n_atoms}")
    if arr.shape[1] < 2:
        raise ValueError("Charge file table must have at least two columns: id q")
    ids = arr[:, 0].astype(int)
    q = arr[:, 1].astype(float)
    order = np.argsort(ids)
    ids_sorted = ids[order]
    if not np.array_equal(ids_sorted, np.arange(1, n_atoms + 1)):
        raise ValueError("Charge-file IDs must be exactly 1..N")
    return q[order]


def composition_dict(symbols) -> Dict[str, int]:
    comp: Dict[str, int] = {}
    for s in symbols:
        comp[s] = comp.get(s, 0) + 1
    return dict(sorted(comp.items()))


def get_charges(atoms, args) -> np.ndarray:
    """Choose charges for dipole calculation."""
    n = len(atoms)
    symbols = np.asarray(atoms.get_chemical_symbols())

    if args.charge_file:
        charges = load_charge_file(args.charge_file, n)
        source = f"charge_file:{args.charge_file}"
    elif args.charge_map:
        cmap = parse_key_value_map(args.charge_map, str, float)
        missing = sorted(set(symbols) - set(cmap))
        if missing:
            raise ValueError(f"Missing charge-map entries for symbols: {missing}")
        charges = np.array([cmap[s] for s in symbols], dtype=float)
        source = f"charge_map:{args.charge_map}"
    else:
        charges = np.asarray(atoms.get_initial_charges(), dtype=float).copy()
        source = "input_initial_charges"

    if args.neutralize:
        charges -= charges.sum() / n

    if np.allclose(charges, 0.0, atol=1e-14) and not args.allow_zero_charges:
        raise ValueError(
            "All dipole charges are zero. FTIR intensity will be zero. "
            "Use --charge-map, --charge-file, or --allow-zero-charges for a force-only test."
        )

    net_q = float(charges.sum())
    if abs(net_q) > args.net_charge_tol:
        raise ValueError(
            f"Net dipole charge is {net_q:.8g} e. Dipole in a periodic cell is origin dependent. "
            "Use neutral charges or pass --neutralize if you intentionally want mean-charge correction."
        )

    atoms.set_initial_charges(charges)
    return charges


def compute_dipole_eang(atoms, charges: np.ndarray) -> np.ndarray:
    """Return dipole vector in e Angstrom using continuous ASE positions."""
    pos = atoms.get_positions(wrap=False)
    return charges @ pos


def build_atoms(args):
    try:
        from ase.io import read
        from ase.data import atomic_numbers
    except Exception as exc:
        raise RuntimeError(
            "ASE is required. Install with: pip install ase"
        ) from exc

    read_kwargs = {}
    if args.format == "lammps-data":
        read_kwargs["atom_style"] = args.atom_style
        read_kwargs["sort_by_id"] = True
        if args.z_of_type:
            raw = parse_key_value_map(args.z_of_type, int, str)
            read_kwargs["Z_of_type"] = {k: atomic_numbers[v] for k, v in raw.items()}

    atoms = read(args.input, format=args.format, **read_kwargs)
    atoms.pbc = [args.pbc_x, args.pbc_y, args.pbc_z]
    return atoms


def build_calculator(args):
    try:
        import torch
        from chgnet.model.dynamics import CHGNetCalculator
    except Exception as exc:
        raise RuntimeError(
            "CHGNet and torch are required. Install with: pip install chgnet torch"
        ) from exc

    if hasattr(torch, "set_float32_matmul_precision"):
        torch.set_float32_matmul_precision(args.matmul_precision)
    torch.set_num_threads(args.torch_threads)

    if args.model_path:
        return CHGNetCalculator.from_file(args.model_path, use_device=args.device)
    return CHGNetCalculator(
        use_device=args.device,
        check_cuda_mem=args.check_cuda_mem,
        on_isolated_atoms=args.on_isolated_atoms,
    )


def initialize_velocities(atoms, temperature_K: float, seed: int) -> None:
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation

    rng = np.random.default_rng(seed)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K, rng=rng)
    Stationary(atoms)
    ZeroRotation(atoms)


def maybe_relax(atoms, args, out_dir: Path) -> None:
    if args.relax_steps <= 0:
        return
    from ase.optimize import FIRE
    from ase.io import write

    print(f"Relaxing with FIRE: steps={args.relax_steps}, fmax={args.relax_fmax} eV/A", flush=True)
    opt = FIRE(atoms, logfile=str(out_dir / "relax.log"))
    opt.run(fmax=args.relax_fmax, steps=args.relax_steps)
    write(out_dir / "relaxed.extxyz", atoms)


def run_md(args) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    atoms = build_atoms(args)
    symbols = atoms.get_chemical_symbols()
    charges = get_charges(atoms, args)
    atoms.calc = build_calculator(args)

    maybe_relax(atoms, args, out_dir)
    initialize_velocities(atoms, args.temperature_K, args.seed)

    from ase import units
    from ase.md.langevin import Langevin
    from ase.io import write
    from ase.io.trajectory import Trajectory

    dyn = Langevin(
        atoms,
        timestep=args.timestep_fs * units.fs,
        temperature_K=args.temperature_K,
        friction=(1.0 / args.langevin_tau_fs) / units.fs,
        fixcm=True,
    )

    metadata = {
        "input": str(args.input),
        "n_atoms": len(atoms),
        "composition": composition_dict(symbols),
        "temperature_K": args.temperature_K,
        "timestep_fs": args.timestep_fs,
        "steps": args.steps,
        "sample_interval": args.sample_interval,
        "traj_interval": args.traj_interval,
        "charge_sum_e": float(charges.sum()),
        "charge_min_e": float(charges.min()),
        "charge_max_e": float(charges.max()),
        "pbc": atoms.pbc.tolist(),
        "cell_A": np.asarray(atoms.cell).tolist(),
        "created_unix": time.time(),
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    dipole_path = out_dir / "dipoles.csv"
    thermo_path = out_dir / "thermo.csv"

    dipole_f = open(dipole_path, "w", newline="", encoding="utf-8")
    thermo_f = open(thermo_path, "w", newline="", encoding="utf-8")
    dipole_writer = csv.writer(dipole_f)
    thermo_writer = csv.writer(thermo_f)
    dipole_writer.writerow(["step", "time_fs", "Mx_eA", "My_eA", "Mz_eA"])
    thermo_writer.writerow(["step", "time_fs", "T_K", "Epot_eV", "Ekin_eV", "Etot_eV"])

    traj = None
    if args.traj_interval > 0:
        traj = Trajectory(str(out_dir / "md.traj"), "w", atoms)
        dyn.attach(traj.write, interval=args.traj_interval)

    def record():
        step = int(dyn.nsteps)
        time_fs = step * args.timestep_fs
        M = compute_dipole_eang(atoms, charges)
        dipole_writer.writerow([step, f"{time_fs:.8f}", f"{M[0]:.12e}", f"{M[1]:.12e}", f"{M[2]:.12e}"])

        epot = float(atoms.get_potential_energy())
        ekin = float(atoms.get_kinetic_energy())
        temp = float(atoms.get_temperature())
        thermo_writer.writerow([step, f"{time_fs:.8f}", f"{temp:.8f}", f"{epot:.12e}", f"{ekin:.12e}", f"{epot + ekin:.12e}"])
        dipole_f.flush()
        thermo_f.flush()
        if args.print_interval > 0 and step % args.print_interval == 0:
            print(f"step={step:8d} time_fs={time_fs:12.3f} T={temp:9.2f} Epot={epot:16.6f}", flush=True)

    dyn.attach(record, interval=args.sample_interval)
    record()

    print("Starting MD", flush=True)
    remaining = args.steps
    while remaining > 0:
        n = min(args.chunk_steps, remaining)
        dyn.run(n)
        remaining -= n
        write(out_dir / "restart.extxyz", atoms)

    dipole_f.close()
    thermo_f.close()
    if traj is not None:
        traj.close()
    print(f"Done. Output directory: {out_dir}", flush=True)


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run CHGNet ASE MD and save dipoles for FTIR.")
    p.add_argument("--input", required=True, help="Input structure. Example: precursor_pH7p0_N192.data")
    p.add_argument("--format", default="lammps-data", help="ASE format. Use lammps-data for LAMMPS .data files.")
    p.add_argument("--atom-style", default="charge", help="LAMMPS atom_style for lammps-data.")
    p.add_argument("--z-of-type", default=None, help="Optional LAMMPS type map. Example: 1:Si,2:O,3:H")
    p.add_argument("--out-dir", default="run_chgnet_md", help="Output directory.")

    p.add_argument("--steps", type=int, default=50000, help="MD steps after optional relaxation.")
    p.add_argument("--timestep-fs", type=float, default=0.5, help="Time step in fs. Use 0.25 to 0.5 fs when H is present.")
    p.add_argument("--temperature-K", type=float, default=300.0, help="Target temperature in K.")
    p.add_argument("--langevin-tau-fs", type=float, default=100.0, help="Langevin damping time in fs. Friction = 1/tau.")
    p.add_argument("--sample-interval", type=int, default=1, help="Save dipole every N MD steps.")
    p.add_argument("--traj-interval", type=int, default=100, help="Save ASE trajectory every N steps. Set 0 to disable.")
    p.add_argument("--print-interval", type=int, default=100, help="Print status every N steps. Set 0 to disable.")
    p.add_argument("--chunk-steps", type=int, default=1000, help="Run MD in chunks and write restart.extxyz after each chunk.")

    p.add_argument("--relax-steps", type=int, default=200, help="FIRE relaxation steps before MD. Set 0 to skip.")
    p.add_argument("--relax-fmax", type=float, default=0.05, help="Relaxation force threshold in eV/A.")

    p.add_argument("--device", default=None, choices=[None, "cpu", "cuda", "mps"], help="Torch device for CHGNet. Default lets CHGNet choose.")
    p.add_argument("--model-path", default=None, help="Optional custom CHGNet model path.")
    p.add_argument("--check-cuda-mem", action="store_true", help="Let CHGNet choose CUDA device with available memory.")
    p.add_argument("--on-isolated-atoms", default="warn", choices=["ignore", "warn", "error"], help="CHGNet isolated atom behavior.")
    p.add_argument("--torch-threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="Torch CPU threads.")
    p.add_argument("--matmul-precision", default="high", choices=["highest", "high", "medium"], help="Torch float32 matmul precision.")

    p.add_argument("--charge-map", default=None, help="Fixed charges by element. Example: Si:2.4,O:-1.2,H:0.6")
    p.add_argument("--charge-file", default=None, help="Text file with charges: either q or id q.")
    p.add_argument("--neutralize", action="store_true", help="Subtract mean charge so total charge is exactly zero.")
    p.add_argument("--net-charge-tol", type=float, default=1e-6, help="Allowed net charge in e.")
    p.add_argument("--allow-zero-charges", action="store_true", help="Allow zero charges for force-only debugging. FTIR will be zero.")

    p.add_argument("--pbc-x", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--pbc-y", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--pbc-z", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--seed", type=int, default=12345)
    return p


if __name__ == "__main__":
    run_md(make_parser().parse_args())
