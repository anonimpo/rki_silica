#!/usr/bin/env python3
"""Generate pH-conditioned hydrated beta-cristobalite / silica precursors.

Purpose
-------
The supplied ReaxFF file contains Si/O/H only. Therefore this generator does not
insert explicit HCl/NaOH/Na+/Cl- species. Instead, pH is represented by a
controlled amount of initial water/hydroxylation capacity, which is then allowed
to react during ReaxFF annealing. This is a mechanistic proxy for the sol-gel
precursor state: acid-side and strong-base-side samples are assigned higher
hydration/defect load, while pH 7.0 is assigned the lowest defect load.

Atom type order in the generated LAMMPS data files is fixed as:
  1 Si, 2 O, 3 H
and must be paired in LAMMPS with:
  pair_coeff * * potentials/ffield.reax.SiOH Si O H
"""
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path
from typing import Iterable, List, Tuple

Vec = Tuple[float, float, float]

SI_MASS = 28.0855
O_MASS = 15.9994
H_MASS = 1.0080

# Conventional diamond-cubic Si sublattice used as beta-cristobalite template.
SI_BASIS_FRAC = [
    (0.0, 0.0, 0.0),
    (0.25, 0.25, 0.25),
    (0.0, 0.5, 0.5),
    (0.25, 0.75, 0.75),
    (0.5, 0.0, 0.5),
    (0.75, 0.25, 0.75),
    (0.5, 0.5, 0.0),
    (0.75, 0.75, 0.25),
]

# Relative number of H2O molecules per Si atom. Symmetric around pH 7.0 to
# reflect the thesis trend: pH 7.0 gives the most ordered post-calcination
# product, whereas both acid and stronger base sides introduce more defects.
PH_WATER_PER_SI = {
    "pH6p0": 0.1875,
    "pH6p5": 0.15625,
    "pH7p0": 0.1250,
    "pH7p5": 0.15625,
    "pH8p0": 0.1875,
}
PH_NUMERIC = {
    "pH6p0": 6.0,
    "pH6p5": 6.5,
    "pH7p0": 7.0,
    "pH7p5": 7.5,
    "pH8p0": 8.0,
}


def wrap(x: float, L: float) -> float:
    y = x % L
    if abs(y - L) < 1.0e-10:
        return 0.0
    return y


def min_image_delta(a: Vec, b: Vec, box: Vec) -> Vec:
    d = [b[k] - a[k] for k in range(3)]
    for k, L in enumerate(box):
        d[k] -= round(d[k] / L) * L
    return (d[0], d[1], d[2])


def dist_pbc(a: Vec, b: Vec, box: Vec) -> float:
    d = min_image_delta(a, b, box)
    return math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def norm(v: Vec) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(v: Vec, s: float) -> Vec:
    return (v[0] * s, v[1] * s, v[2] * s)


def unit_random(rng: random.Random) -> Vec:
    z = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    r = math.sqrt(max(0.0, 1.0 - z * z))
    return (r * math.cos(phi), r * math.sin(phi), z)


def perpendicular_unit(v: Vec, rng: random.Random) -> Vec:
    # Generate a vector approximately perpendicular to v.
    for _ in range(100):
        u = unit_random(rng)
        dot = v[0] * u[0] + v[1] * u[1] + v[2] * u[2]
        w = (u[0] - dot * v[0], u[1] - dot * v[1], u[2] - dot * v[2])
        n = norm(w)
        if n > 1.0e-8:
            return (w[0] / n, w[1] / n, w[2] / n)
    return (1.0, 0.0, 0.0)


def generate_beta_cristobalite(a: float, reps: Tuple[int, int, int]) -> Tuple[List[str], List[Vec], Vec]:
    nx, ny, nz = reps
    box = (a * nx, a * ny, a * nz)
    si_positions: List[Vec] = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for fx, fy, fz in SI_BASIS_FRAC:
                    si_positions.append(((ix + fx) * a, (iy + fy) * a, (iz + fz) * a))

    nn = a * math.sqrt(3.0) / 4.0
    cutoff = nn * 1.08
    o_positions: List[Vec] = []
    for i in range(len(si_positions)):
        for j in range(i + 1, len(si_positions)):
            d = min_image_delta(si_positions[i], si_positions[j], box)
            if norm(d) < cutoff:
                m = add(si_positions[i], scale(d, 0.5))
                o_positions.append((wrap(m[0], box[0]), wrap(m[1], box[1]), wrap(m[2], box[2])))

    expected_o = 2 * len(si_positions)
    if len(o_positions) != expected_o:
        raise RuntimeError(f"Expected {expected_o} bridging O, found {len(o_positions)}")

    symbols = ["Si"] * len(si_positions) + ["O"] * len(o_positions)
    positions = si_positions + o_positions
    return symbols, positions, box


def scale_box(symbols: List[str], positions: List[Vec], box: Vec, factor: float) -> Tuple[List[str], List[Vec], Vec]:
    new_box = (box[0] * factor, box[1] * factor, box[2] * factor)
    scaled = [(x * factor, y * factor, z * factor) for x, y, z in positions]
    return list(symbols), scaled, new_box


def jitter_positions(positions: List[Vec], box: Vec, rng: random.Random, amp: float) -> List[Vec]:
    if amp <= 0.0:
        return positions
    out = []
    for x, y, z in positions:
        out.append((
            wrap(x + rng.uniform(-amp, amp), box[0]),
            wrap(y + rng.uniform(-amp, amp), box[1]),
            wrap(z + rng.uniform(-amp, amp), box[2]),
        ))
    return out


def add_waters(symbols: List[str], positions: List[Vec], box: Vec, n_water: int, rng: random.Random,
               min_ow_dist: float = 1.35) -> Tuple[List[str], List[Vec]]:
    """Append water molecules as O, H, H with simple gas-phase geometry.

    The positions are only initial guesses for ReaxFF relaxation. The molecule
    topology is not fixed; ReaxFF decides bond order dynamically.
    """
    symbols_out = list(symbols)
    positions_out = list(positions)
    oh = 0.96
    angle = math.radians(104.5)
    half = 0.5 * angle

    for w in range(n_water):
        placed = False
        min_try_dist = min_ow_dist
        for attempt in range(8000):
            o = (rng.random() * box[0], rng.random() * box[1], rng.random() * box[2])
            if all(dist_pbc(o, p, box) >= min_try_dist for p in positions_out):
                axis = unit_random(rng)
                perp = perpendicular_unit(axis, rng)
                h1_dir = add(scale(axis, math.cos(half)), scale(perp, math.sin(half)))
                h2_dir = add(scale(axis, math.cos(half)), scale(perp, -math.sin(half)))
                h1 = (wrap(o[0] + oh * h1_dir[0], box[0]), wrap(o[1] + oh * h1_dir[1], box[1]), wrap(o[2] + oh * h1_dir[2], box[2]))
                h2 = (wrap(o[0] + oh * h2_dir[0], box[0]), wrap(o[1] + oh * h2_dir[1], box[1]), wrap(o[2] + oh * h2_dir[2], box[2]))
                symbols_out.extend(["O", "H", "H"])
                positions_out.extend([o, h1, h2])
                placed = True
                break
            if attempt in (2500, 5000, 7000):
                min_try_dist *= 0.90
        if not placed:
            raise RuntimeError(f"Could not place water molecule {w+1}/{n_water}. Increase --scale-box.")
    return symbols_out, positions_out


def reorder_by_type(symbols: List[str], positions: List[Vec]) -> Tuple[List[str], List[Vec]]:
    order = {"Si": 0, "O": 1, "H": 2}
    items = sorted(zip(symbols, positions), key=lambda x: order[x[0]])
    return [s for s, _ in items], [p for _, p in items]


def write_lammps_data(path: Path, symbols: List[str], positions: List[Vec], box: Vec) -> None:
    type_map = {"Si": 1, "O": 2, "H": 3}
    masses = {1: (SI_MASS, "Si"), 2: (O_MASS, "O"), 3: (H_MASS, "H")}
    with path.open("w", encoding="utf-8") as f:
        f.write("# pH-conditioned hydrated silica precursor for ReaxFF (atom_style charge)\n\n")
        f.write(f"{len(symbols)} atoms\n")
        f.write("3 atom types\n\n")
        f.write(f"0.0 {box[0]:.10f} xlo xhi\n")
        f.write(f"0.0 {box[1]:.10f} ylo yhi\n")
        f.write(f"0.0 {box[2]:.10f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in (1, 2, 3):
            mass, name = masses[t]
            f.write(f"{t} {mass:.6f} # {name}\n")
        f.write("\nAtoms # charge\n\n")
        for i, (sym, (x, y, z)) in enumerate(zip(symbols, positions), start=1):
            f.write(f"{i:8d} {type_map[sym]:2d} {0.0: .8f} {x:16.8f} {y:16.8f} {z:16.8f} # {sym}\n")


def write_xyz(path: Path, symbols: List[str], positions: List[Vec], box: Vec) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(symbols)}\n")
        f.write(f'Lattice="{box[0]:.10f} 0 0 0 {box[1]:.10f} 0 0 0 {box[2]:.10f}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for sym, (x, y, z) in zip(symbols, positions):
            f.write(f"{sym:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def write_plumed_localq6(path: Path, ph_label: str, base_atoms: int, n_si: int, temp: float) -> None:
    """Write stock-PLUMED Local-Q6 WTMetaD input.

    This is not the exact Niu XRD/Debye CV. It is a stock PLUMED fallback for
    exploratory enhanced sampling when the custom XRD CV is unavailable.
    """
    with path.open("w", encoding="utf-8") as f:
        f.write("# Stock PLUMED fallback: Local-Q6 order parameter on the Si sublattice.\n")
        f.write("# Requires PLUMED built with the crystallization module.\n")
        f.write("# This is NOT the exact Niu et al. XRD/Debye collective variable.\n")
        f.write("UNITS LENGTH=A TIME=fs ENERGY=kcal/mol\n\n")
        f.write(f"# {ph_label}: Si atoms are contiguous because data files are written Si/O/H.\n")
        f.write(f"q6: Q6 SPECIES=1-{n_si} SWITCH={{RATIONAL D_0=2.7 R_0=0.6}}\n")
        f.write("lq6: LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=2.7 R_0=0.6} MEAN\n\n")
        f.write(f"metad: METAD ARG=lq6.mean SIGMA=0.02 HEIGHT=0.50 PACE=2000 BIASFACTOR=50 TEMP={temp:.1f} FILE=runs/{ph_label}_N{base_atoms}/05_wtmetad_plumed/HILLS_{ph_label}_N{base_atoms} GRID_MIN=-0.30 GRID_MAX=1.00 GRID_BIN=650 CALC_RCT\n")
        f.write(f"\nPRINT ARG=lq6.mean,metad.bias,metad.rct,metad.rbias STRIDE=100 FILE=runs/{ph_label}_N{base_atoms}/05_wtmetad_plumed/COLVAR\n")
        f.write("FLUSH STRIDE=100\n")


def write_manifest(path: Path, rows: List[dict]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("ph_label,ph_value,base_atoms,n_si,n_o_base,n_h2o,n_o_total,n_h_total,total_atoms,box_A\n")
        for r in rows:
            f.write(
                f"{r['ph_label']},{r['ph_value']:.1f},{r['base_atoms']},{r['n_si']},{r['n_o_base']},"
                f"{r['n_h2o']},{r['n_o_total']},{r['n_h_total']},{r['total_atoms']},{r['box_A']:.6f}\n"
            )


def reps_from_base_atoms(base_atoms: int) -> Tuple[int, int, int]:
    # One conventional cell = 24 atoms. We support 192=2^3*24 and 1536=4^3*24 by default.
    cells = base_atoms // 24
    root = round(cells ** (1.0 / 3.0))
    if 24 * root ** 3 != base_atoms:
        raise ValueError("base_atoms must be 24*n^3, e.g., 192 or 1536")
    return (root, root, root)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-atoms", type=int, default=1536, help="SiO2 template atoms before water addition: 192 or 1536 recommended")
    ap.add_argument("--a", type=float, default=7.15, help="beta-cristobalite conventional-cell lattice parameter in Angstrom")
    ap.add_argument("--scale-box", type=float, default=1.18, help="scale box/coordinates to create a porous hydrated precursor")
    ap.add_argument("--jitter", type=float, default=0.08, help="random displacement amplitude in Angstrom")
    ap.add_argument("--seed", type=int, default=20260601)
    ap.add_argument("--structures-dir", default="structures")
    ap.add_argument("--plumed-dir", default="plumed/stock_localq6")
    ap.add_argument("--temp-meta", type=float, default=2300.0)
    args = ap.parse_args()

    reps = reps_from_base_atoms(args.base_atoms)
    rng_master = random.Random(args.seed)
    rows = []
    structures_dir = Path(args.structures_dir)
    plumed_dir = Path(args.plumed_dir)
    structures_dir.mkdir(parents=True, exist_ok=True)
    plumed_dir.mkdir(parents=True, exist_ok=True)

    base_symbols, base_pos, base_box = generate_beta_cristobalite(args.a, reps)
    base_symbols, base_pos, box = scale_box(base_symbols, base_pos, base_box, args.scale_box)
    n_si = base_symbols.count("Si")
    n_o_base = base_symbols.count("O")

    for ph_label, water_per_si in PH_WATER_PER_SI.items():
        rng = random.Random(rng_master.randint(1, 10**9))
        n_water = int(round(water_per_si * n_si))
        symbols = list(base_symbols)
        pos = jitter_positions(list(base_pos), box, rng, args.jitter)
        symbols, pos = add_waters(symbols, pos, box, n_water, rng)
        symbols, pos = reorder_by_type(symbols, pos)

        stem = f"precursor_{ph_label}_N{args.base_atoms}"
        write_lammps_data(structures_dir / f"{stem}.data", symbols, pos, box)
        write_xyz(structures_dir / f"{stem}.xyz", symbols, pos, box)
        write_plumed_localq6(plumed_dir / f"plumed_localq6_{ph_label}_N{args.base_atoms}.dat", ph_label, args.base_atoms, n_si, args.temp_meta)

        rows.append({
            "ph_label": ph_label,
            "ph_value": PH_NUMERIC[ph_label],
            "base_atoms": args.base_atoms,
            "n_si": n_si,
            "n_o_base": n_o_base,
            "n_h2o": n_water,
            "n_o_total": symbols.count("O"),
            "n_h_total": symbols.count("H"),
            "total_atoms": len(symbols),
            "box_A": box[0],
        })
        print(f"Wrote {stem}: total={len(symbols)} Si={n_si} O={symbols.count('O')} H={symbols.count('H')} H2O={n_water}")

    write_manifest(structures_dir / f"precursor_manifest_N{args.base_atoms}.csv", rows)
    print(f"Manifest: {structures_dir / f'precursor_manifest_N{args.base_atoms}.csv'}")


if __name__ == "__main__":
    main()
