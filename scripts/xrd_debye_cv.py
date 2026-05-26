#!/usr/bin/env python3
"""Post-process a LAMMPS dump or XYZ trajectory and compute Debye XRD-like CVs.

Computes the quantity used by Niu et al. in spirit:
    I(Q) = (1/N) sum_i sum_j f_i(Q) f_j(Q) sinc(Q r_ij) W(r_ij)
with a Lorch window W(r)=sinc(pi r/Rc) and minimum-image PBC.

This is useful for validating trajectories and for the two-executable workflow
(ReaxFF in LAMMPS; PLUMED/analysis separately). It does NOT apply a bias to MD.
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

Vec = Tuple[float, float, float]

# Cromer-Mann neutral-atom coefficients. f = sum a_i exp[-b_i (Q/4pi)^2] + c.
CM = {
    "O":  ([3.0485, 2.2868, 1.5463, 0.8670], [13.2771, 5.7011, 0.3239, 32.9089], 0.2508),
    "Si": ([6.2915, 3.0353, 1.9891, 1.5410], [2.4386, 32.3337, 0.6785, 81.6937], 1.1407),
}
Z = {"O": 8.0, "Si": 14.0}


def sinc(x: float) -> float:
    return 1.0 if abs(x) < 1.0e-12 else math.sin(x) / x


def ff(symbol: str, q: float, mode: str) -> float:
    if mode == "constant":
        return Z.get(symbol, 1.0)
    if mode == "none":
        return 1.0
    if symbol not in CM:
        return Z.get(symbol, 1.0)
    a, b, c = CM[symbol]
    s2 = (q / (4.0 * math.pi)) ** 2
    return sum(ai * math.exp(-bi * s2) for ai, bi in zip(a, b)) + c


def minimum_image(dx: float, L: float) -> float:
    return dx - round(dx / L) * L


def dist_pbc(a: Vec, b: Vec, box: Vec) -> float:
    dx = minimum_image(b[0] - a[0], box[0])
    dy = minimum_image(b[1] - a[1], box[1])
    dz = minimum_image(b[2] - a[2], box[2])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def debye(symbols: List[str], pos: List[Vec], box: Vec, q: float, rc: Optional[float], formfactor: str) -> float:
    n = len(pos)
    if rc is None:
        rc = 0.5 * min(box)
    f = [ff(s, q, formfactor) for s in symbols]
    total = 0.0
    for i in range(n):
        # self term
        total += f[i] * f[i]
        for j in range(i + 1, n):
            r = dist_pbc(pos[i], pos[j], box)
            if r > rc:
                continue
            term = f[i] * f[j] * sinc(q * r) * sinc(math.pi * r / rc)
            total += 2.0 * term
    return total / n


def parse_xyz(path: Path) -> Iterator[Tuple[int, List[str], List[Vec], Vec]]:
    with path.open("r", encoding="utf-8") as f:
        frame = 0
        while True:
            line = f.readline()
            if not line:
                return
            line = line.strip()
            if not line:
                continue
            n = int(line)
            comment = f.readline().strip()
            m = re.search(r'Lattice="([^"]+)"', comment)
            if not m:
                raise ValueError("XYZ comment lacks Lattice=\"...\" box information")
            vals = [float(x) for x in m.group(1).split()]
            box = (vals[0], vals[4], vals[8])
            symbols: List[str] = []
            pos: List[Vec] = []
            for _ in range(n):
                parts = f.readline().split()
                symbols.append(parts[0])
                pos.append((float(parts[1]), float(parts[2]), float(parts[3])))
            yield frame, symbols, pos, box
            frame += 1


def parse_lammpstrj(path: Path, type_map: Dict[int, str]) -> Iterator[Tuple[int, List[str], List[Vec], Vec]]:
    with path.open("r", encoding="utf-8") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(f.readline().strip())
            assert f.readline().startswith("ITEM: NUMBER")
            n = int(f.readline().strip())
            box_header = f.readline().strip()
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Unsupported dump: missing BOX BOUNDS")
            bounds = []
            for _ in range(3):
                lo, hi, *_ = f.readline().split()
                bounds.append((float(lo), float(hi)))
            box = tuple(hi - lo for lo, hi in bounds)  # type: ignore[assignment]
            atom_header = f.readline().strip().split()[2:]
            col = {name: i for i, name in enumerate(atom_header)}
            needed = ["id", "type"]
            for name in needed:
                if name not in col:
                    raise ValueError(f"Dump lacks {name} column")
            # Accept x/y/z or xu/yu/zu or xs/ys/zs.
            if all(k in col for k in ("x", "y", "z")):
                coord_kind = "x"
            elif all(k in col for k in ("xu", "yu", "zu")):
                coord_kind = "xu"
            elif all(k in col for k in ("xs", "ys", "zs")):
                coord_kind = "xs"
            else:
                raise ValueError("Dump needs x y z, xu yu zu, or xs ys zs columns")
            rows = []
            for _ in range(n):
                parts = f.readline().split()
                rows.append(parts)
            rows.sort(key=lambda p: int(p[col["id"]]))
            symbols: List[str] = []
            pos: List[Vec] = []
            for parts in rows:
                t = int(parts[col["type"]])
                symbols.append(type_map[t])
                if coord_kind == "x":
                    x, y, z = (float(parts[col["x"]]), float(parts[col["y"]]), float(parts[col["z"]]))
                elif coord_kind == "xu":
                    x, y, z = (float(parts[col["xu"]]), float(parts[col["yu"]]), float(parts[col["zu"]]))
                else:
                    x = float(parts[col["xs"]]) * box[0]
                    y = float(parts[col["ys"]]) * box[1]
                    z = float(parts[col["zs"]]) * box[2]
                pos.append((x % box[0], y % box[1], z % box[2]))
            yield timestep, symbols, pos, box  # type: ignore[arg-type]


def q_from_hkl(a: float, hkl: Tuple[int, int, int]) -> float:
    h, k, l = hkl
    return 2.0 * math.pi * math.sqrt(h * h + k * k + l * l) / a


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", help=".xyz or LAMMPS .lammpstrj/.dump")
    ap.add_argument("--q", type=float, nargs="*", help="Q values in inverse Angstrom")
    ap.add_argument("--a", type=float, default=7.15, help="lattice a for default q111 and q022, Angstrom")
    ap.add_argument("--rc", type=float, default=None, help="cutoff Rc in Angstrom; default 0.5*min(box)")
    ap.add_argument("--formfactor", choices=["cromer", "constant", "none"], default="cromer")
    ap.add_argument("--type-map", default="1:Si,2:O", help="LAMMPS dump type map, e.g. 1:Si,2:O")
    ap.add_argument("--stride", type=int, default=1)
    args = ap.parse_args()

    qs = args.q if args.q else [q_from_hkl(args.a, (1, 1, 1)), q_from_hkl(args.a, (0, 2, 2))]
    p = Path(args.trajectory)
    if p.suffix.lower() == ".xyz":
        frames = parse_xyz(p)
    else:
        type_map = {int(k): v for k, v in (item.split(":") for item in args.type_map.split(","))}
        frames = parse_lammpstrj(p, type_map)

    print("# step " + " ".join(f"I_Q{q:.5f}" for q in qs))
    for iframe, (step, symbols, pos, box) in enumerate(frames):
        if iframe % args.stride != 0:
            continue
        vals = [debye(symbols, pos, box, q, args.rc, args.formfactor) for q in qs]
        print(str(step) + " " + " ".join(f"{v:.10g}" for v in vals))


if __name__ == "__main__":
    main()
