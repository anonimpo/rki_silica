#!/usr/bin/env python3
"""Distance-based structural metrics for hydrated silica ReaxFF outputs.

Metrics are approximate but useful for comparing pH trends:
- Si tetrahedral fraction from Si-O coordination.
- Bridging oxygen and non-bridging oxygen fractions.
- Silanol proxy: O bonded to one Si and at least one H.
- Water proxy: O bonded to zero Si and at least two H.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

from importlib.machinery import SourceFileLoader
_xrd = SourceFileLoader("debye", str(Path(__file__).with_name("04_compute_debye_xrd.py"))).load_module()


def pbc_dist_matrix(a: np.ndarray, b: np.ndarray, box: np.ndarray) -> np.ndarray:
    d = a[:, None, :] - b[None, :, :]
    d -= np.round(d / box) * box
    return np.linalg.norm(d, axis=2)


def analyze(types: np.ndarray, coords: np.ndarray, box: np.ndarray, cutoff_sio: float, cutoff_oh: float) -> dict:
    si = coords[types == 1]
    oxy = coords[types == 2]
    hyd = coords[types == 3]

    if len(si) == 0 or len(oxy) == 0:
        raise ValueError("Si and O atoms are required")

    d_sio = pbc_dist_matrix(si, oxy, box)
    si_coord = np.sum(d_sio < cutoff_sio, axis=1)
    o_si_coord = np.sum(d_sio < cutoff_sio, axis=0)

    if len(hyd) > 0:
        d_oh = pbc_dist_matrix(oxy, hyd, box)
        o_h_coord = np.sum(d_oh < cutoff_oh, axis=1)
    else:
        o_h_coord = np.zeros(len(oxy), dtype=int)

    n_si = len(si)
    n_o = len(oxy)
    n_h = len(hyd)
    out = {
        "n_si": n_si,
        "n_o": n_o,
        "n_h": n_h,
        "si_o_coord_mean": float(np.mean(si_coord)),
        "si_o_coord_std": float(np.std(si_coord)),
        "frac_si_coord4": float(np.mean(si_coord == 4)),
        "frac_si_undercoord": float(np.mean(si_coord < 4)),
        "frac_si_overcoord": float(np.mean(si_coord > 4)),
        "frac_o_bridging_2si_0h": float(np.mean((o_si_coord == 2) & (o_h_coord == 0))),
        "frac_o_nonbridging_1si": float(np.mean(o_si_coord == 1)),
        "frac_o_silanol_proxy_1si_ge1h": float(np.mean((o_si_coord == 1) & (o_h_coord >= 1))),
        "frac_o_water_proxy_0si_ge2h": float(np.mean((o_si_coord == 0) & (o_h_coord >= 2))),
        "siloxane_bridge_count": int(np.sum((o_si_coord == 2) & (o_h_coord == 0))),
        "silanol_proxy_count": int(np.sum((o_si_coord == 1) & (o_h_coord >= 1))),
        "water_proxy_count": int(np.sum((o_si_coord == 0) & (o_h_coord >= 2))),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--format", choices=["data", "dump"], default="data")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ph", default="NA")
    ap.add_argument("--cutoff-sio", type=float, default=2.0)
    ap.add_argument("--cutoff-oh", type=float, default=1.25)
    args = ap.parse_args()

    path = Path(args.input)
    if args.format == "data":
        types, coords, box = _xrd.read_data(path)
    else:
        types, coords, box = _xrd.read_last_dump(path)
    metrics = analyze(types, coords, box, args.cutoff_sio, args.cutoff_oh)
    metrics = {"ph_label": args.ph, "source": str(path), **metrics}

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = list(metrics.keys())
    with out.open("w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        f.write(",".join(str(metrics[k]) for k in keys) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
