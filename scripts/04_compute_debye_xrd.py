#!/usr/bin/env python3
"""Compute powder-like XRD intensity from LAMMPS data/dump using the Debye equation.

This post-processing script is intentionally independent from PLUMED. It provides
an observable comparable with the experimental XRD trend: amorphous halo vs.
cristobalite/tridymite-like peak sharpening.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

TYPE_TO_ELEM = {1: "Si", 2: "O", 3: "H"}
LAMBDA_CUKA = 1.5406  # Angstrom

# Cromer-Mann X-ray atomic form factor coefficients. q units: Angstrom^-1,
# s = q / (4*pi). H contributes weakly but is included for hydrated systems.
CM = {
    "H": {
        "a": [0.489918, 0.262003, 0.196767, 0.049879],
        "b": [20.6593, 7.74039, 49.5519, 2.20159],
        "c": 0.001305,
    },
    "O": {
        "a": [3.0485, 2.2868, 1.5463, 0.8670],
        "b": [13.2771, 5.7011, 0.3239, 32.9089],
        "c": 0.2508,
    },
    "Si": {
        "a": [6.2915, 3.0353, 1.9891, 1.5410],
        "b": [2.4386, 32.3337, 0.6785, 81.6937],
        "c": 1.1407,
    },
}


def f_xray(elem: str, q: np.ndarray) -> np.ndarray:
    s2 = (q / (4.0 * np.pi)) ** 2
    pars = CM[elem]
    out = np.full_like(q, pars["c"], dtype=float)
    for a, b in zip(pars["a"], pars["b"]):
        out += a * np.exp(-b * s2)
    return out


def read_data(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    box = np.zeros((3, 2), dtype=float)
    for line in lines:
        parts = line.split()
        if len(parts) >= 4 and parts[-2:] == ["xlo", "xhi"]:
            box[0] = [float(parts[0]), float(parts[1])]
        elif len(parts) >= 4 and parts[-2:] == ["ylo", "yhi"]:
            box[1] = [float(parts[0]), float(parts[1])]
        elif len(parts) >= 4 and parts[-2:] == ["zlo", "zhi"]:
            box[2] = [float(parts[0]), float(parts[1])]
    atom_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atom_start = i + 2
            break
    if atom_start is None:
        raise ValueError("Atoms section not found")
    rows = []
    for line in lines[atom_start:]:
        if not line.strip() or line.strip().startswith("#"):
            continue
        if line.split()[0].isalpha():
            break
        parts = line.split("#", 1)[0].split()
        if len(parts) < 6:
            break
        # atom_style charge: id type q x y z
        rows.append((int(parts[1]), float(parts[3]), float(parts[4]), float(parts[5])))
    types = np.array([r[0] for r in rows], dtype=int)
    coords = np.array([[r[1], r[2], r[3]] for r in rows], dtype=float)
    lengths = box[:, 1] - box[:, 0]
    coords -= box[:, 0]
    coords %= lengths
    return types, coords, lengths


def read_last_dump(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    idxs = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    if not idxs:
        raise ValueError("No ITEM: TIMESTEP blocks found")
    i = idxs[-1]
    n = int(lines[i + 3].strip())
    box_lines = lines[i + 5:i + 8]
    bounds = np.array([[float(x) for x in line.split()[:2]] for line in box_lines])
    lengths = bounds[:, 1] - bounds[:, 0]
    header = lines[i + 8].split()[2:]
    data = lines[i + 9:i + 9 + n]
    col = {name: k for k, name in enumerate(header)}
    types = np.empty(n, dtype=int)
    coords = np.empty((n, 3), dtype=float)
    for r, line in enumerate(data):
        p = line.split()
        types[r] = int(p[col["type"]])
        if {"x", "y", "z"}.issubset(col):
            coords[r] = [float(p[col["x"]]), float(p[col["y"]]), float(p[col["z"]])]
            coords[r] -= bounds[:, 0]
        elif {"xs", "ys", "zs"}.issubset(col):
            coords[r] = [float(p[col["xs"]]) * lengths[0], float(p[col["ys"]]) * lengths[1], float(p[col["zs"]]) * lengths[2]]
        else:
            raise ValueError("Dump must contain x y z or xs ys zs columns")
    coords %= lengths
    return types, coords, lengths


def pair_histograms(types: np.ndarray, coords: np.ndarray, box: np.ndarray, dr: float, rmax: float | None) -> Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]:
    if rmax is None:
        rmax = 0.5 * float(np.min(box))
    edges = np.arange(0.0, rmax + dr, dr)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist: Dict[Tuple[int, int], np.ndarray] = {}
    n = len(types)
    for i in range(n - 1):
        d = coords[i + 1:] - coords[i]
        d -= np.round(d / box) * box
        r = np.linalg.norm(d, axis=1)
        mask = r < rmax
        if not np.any(mask):
            continue
        js_types = types[i + 1:][mask]
        rsel = r[mask]
        for tj in np.unique(js_types):
            pair = tuple(sorted((int(types[i]), int(tj))))
            if pair not in hist:
                hist[pair] = np.zeros(len(centers), dtype=float)
            h, _ = np.histogram(rsel[js_types == tj], bins=edges)
            hist[pair] += h
    return {k: (centers, v) for k, v in hist.items()}


def compute_xrd(types: np.ndarray, coords: np.ndarray, box: np.ndarray, qmin: float, qmax: float, nq: int, dr: float, rmax: float | None) -> Tuple[np.ndarray, np.ndarray]:
    q = np.linspace(qmin, qmax, nq)
    hist = pair_histograms(types, coords, box, dr=dr, rmax=rmax)
    counts = {t: int(np.count_nonzero(types == t)) for t in np.unique(types)}

    intensity = np.zeros_like(q)
    # Self-scattering term.
    for t, count in counts.items():
        elem = TYPE_TO_ELEM.get(int(t), "O")
        intensity += count * f_xray(elem, q) ** 2

    # Pair terms: 2 * f_i f_j sin(qr)/(qr).
    for (ta, tb), (r, h) in hist.items():
        if np.sum(h) == 0:
            continue
        fa = f_xray(TYPE_TO_ELEM.get(ta, "O"), q)
        fb = f_xray(TYPE_TO_ELEM.get(tb, "O"), q)
        qr = np.outer(q, r)
        sinc = np.ones_like(qr)
        nz = qr != 0.0
        sinc[nz] = np.sin(qr[nz]) / qr[nz]
        pair_sum = sinc @ h
        intensity += 2.0 * fa * fb * pair_sum

    # Normalize for easier comparison between pH systems with different H2O contents.
    denom = max(1.0, float(len(types)))
    intensity /= denom
    intensity -= np.nanmin(intensity)
    mx = np.nanmax(intensity)
    if mx > 0:
        intensity /= mx
    return q, intensity


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="LAMMPS data file or dump trajectory")
    ap.add_argument("--format", choices=["data", "dump"], default="data")
    ap.add_argument("--out", required=True, help="CSV output")
    ap.add_argument("--plot", default=None, help="Optional PNG plot")
    ap.add_argument("--qmin", type=float, default=0.5)
    ap.add_argument("--qmax", type=float, default=4.5)
    ap.add_argument("--nq", type=int, default=801)
    ap.add_argument("--dr", type=float, default=0.02)
    ap.add_argument("--rmax", type=float, default=None)
    args = ap.parse_args()

    path = Path(args.input)
    if args.format == "data":
        types, coords, box = read_data(path)
    else:
        types, coords, box = read_last_dump(path)
    q, intensity = compute_xrd(types, coords, box, args.qmin, args.qmax, args.nq, args.dr, args.rmax)
    theta2 = 2.0 * np.degrees(np.arcsin(np.clip(q * LAMBDA_CUKA / (4.0 * np.pi), -1.0, 1.0)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    arr = np.column_stack([q, theta2, intensity])
    np.savetxt(out, arr, delimiter=",", header="q_Ainv,two_theta_deg,I_norm", comments="")
    print(f"Wrote {out}")

    if args.plot:
        import matplotlib.pyplot as plt
        p = Path(args.plot)
        p.parent.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(7, 4))
        plt.plot(theta2, intensity, linewidth=1.2)
        plt.xlabel(r"2$\theta$ (degree, Cu K$\alpha$)")
        plt.ylabel("Normalized Debye intensity")
        plt.title(path.name)
        plt.xlim(float(np.nanmin(theta2)), float(np.nanmax(theta2)))
        plt.tight_layout()
        plt.savefig(p, dpi=300)
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
