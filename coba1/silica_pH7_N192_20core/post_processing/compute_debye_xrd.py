#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np

TYPE_TO_ELEM = {1: "Si", 2: "O", 3: "H"}
LAMBDA_CUKA = 1.5406
CM = {
    "H": {"a": [0.489918, 0.262003, 0.196767, 0.049879], "b": [20.6593, 7.74039, 49.5519, 2.20159], "c": 0.001305},
    "O": {"a": [3.0485, 2.2868, 1.5463, 0.8670], "b": [13.2771, 5.7011, 0.3239, 32.9089], "c": 0.2508},
    "Si": {"a": [6.2915, 3.0353, 1.9891, 1.5410], "b": [2.4386, 32.3337, 0.6785, 81.6937], "c": 1.1407},
}

def f_xray(elem: str, q: np.ndarray) -> np.ndarray:
    s2 = (q / (4.0 * np.pi)) ** 2
    pars = CM[elem]
    out = np.full_like(q, pars["c"], dtype=float)
    for a, b in zip(pars["a"], pars["b"]):
        out += a * np.exp(-b * s2)
    return out

def read_data(path: Path):
    lines = path.read_text(errors="ignore").splitlines()
    box = np.zeros((3, 2))
    for line in lines:
        p = line.split()
        if len(p) >= 4 and p[-2:] == ["xlo", "xhi"]: box[0] = [float(p[0]), float(p[1])]
        if len(p) >= 4 and p[-2:] == ["ylo", "yhi"]: box[1] = [float(p[0]), float(p[1])]
        if len(p) >= 4 and p[-2:] == ["zlo", "zhi"]: box[2] = [float(p[0]), float(p[1])]
    start = next((i + 2 for i, l in enumerate(lines) if l.strip().startswith("Atoms")), None)
    if start is None: raise ValueError(f"Atoms section not found in {path}")
    rows = []
    for line in lines[start:]:
        s = line.strip()
        if not s or s.startswith("#"): continue
        p = line.split("#", 1)[0].split()
        if len(p) < 6: break
        try:
            rows.append((int(p[1]), float(p[3]), float(p[4]), float(p[5])))
        except ValueError:
            break
    types = np.array([r[0] for r in rows], int)
    coords = np.array([[r[1], r[2], r[3]] for r in rows], float)
    lengths = box[:, 1] - box[:, 0]
    coords = (coords - box[:, 0]) % lengths
    return types, coords, lengths

def read_dump_frames(path: Path):
    lines = path.read_text(errors="ignore").splitlines()
    idxs = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    for i in idxs:
        step = int(lines[i + 1].strip())
        n = int(lines[i + 3].strip())
        bounds = np.array([[float(x) for x in lines[i + 5 + j].split()[:2]] for j in range(3)])
        lengths = bounds[:, 1] - bounds[:, 0]
        header = lines[i + 8].split()[2:]
        col = {name: k for k, name in enumerate(header)}
        types = np.empty(n, int); coords = np.empty((n, 3), float)
        for r, line in enumerate(lines[i + 9:i + 9 + n]):
            p = line.split(); types[r] = int(p[col["type"]])
            if {"x", "y", "z"}.issubset(col):
                coords[r] = [float(p[col["x"]]), float(p[col["y"]]), float(p[col["z"]])] - bounds[:, 0]
            else:
                coords[r] = [float(p[col["xs"]]) * lengths[0], float(p[col["ys"]]) * lengths[1], float(p[col["zs"]]) * lengths[2]]
        yield step, types, coords % lengths, lengths

def read_last_dump(path: Path):
    last = None
    for last in read_dump_frames(path):
        pass
    if last is None: raise ValueError(f"No dump frames found in {path}")
    return last[1], last[2], last[3]

def pair_histograms(types, coords, box, dr=0.02, rmax=None):
    if rmax is None: rmax = 0.5 * float(np.min(box))
    edges = np.arange(0.0, rmax + dr, dr)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist = {}
    n = len(types)
    for i in range(n - 1):
        d = coords[i + 1:] - coords[i]
        d -= np.round(d / box) * box
        r = np.linalg.norm(d, axis=1)
        mask = r < rmax
        for tj in np.unique(types[i + 1:][mask]):
            pair = tuple(sorted((int(types[i]), int(tj))))
            hist.setdefault(pair, np.zeros(len(centers)))
            h, _ = np.histogram(r[mask][types[i + 1:][mask] == tj], bins=edges)
            hist[pair] += h
    return {k: (centers, v) for k, v in hist.items()}

def compute_xrd(types, coords, box, qmin=0.5, qmax=4.5, nq=801, dr=0.02, rmax=None):
    q = np.linspace(qmin, qmax, nq)
    hist = pair_histograms(types, coords, box, dr, rmax)
    intensity = np.zeros_like(q)
    for t in np.unique(types):
        elem = TYPE_TO_ELEM.get(int(t), "O")
        intensity += np.count_nonzero(types == t) * f_xray(elem, q) ** 2
    for (ta, tb), (r, h) in hist.items():
        fa = f_xray(TYPE_TO_ELEM.get(ta, "O"), q); fb = f_xray(TYPE_TO_ELEM.get(tb, "O"), q)
        qr = np.outer(q, r)
        sinc = np.ones_like(qr)
        nz = qr != 0.0
        sinc[nz] = np.sin(qr[nz]) / qr[nz]
        intensity += 2.0 * fa * fb * (sinc @ h)
    intensity /= max(1, len(types))
    intensity -= intensity.min()
    if intensity.max() > 0: intensity /= intensity.max()
    theta2 = 2.0 * np.degrees(np.arcsin(np.clip(q * LAMBDA_CUKA / (4.0 * np.pi), -1, 1)))
    return q, theta2, intensity

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--format", choices=["data", "dump"], default="data")
    ap.add_argument("--out", required=True)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()
    if args.format == "data": types, coords, box = read_data(Path(args.input))
    else: types, coords, box = read_last_dump(Path(args.input))
    q, t2, I = compute_xrd(types, coords, box)
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, np.column_stack([q, t2, I]), delimiter=",", header="q_Ainv,two_theta_deg,I_norm", comments="")
    if args.plot:
        import matplotlib.pyplot as plt
        p = Path(args.plot); p.parent.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(7, 4))
        plt.plot(t2, I, linewidth=1.2)
        plt.xlabel(r"2$\theta$ (degree, Cu K$\alpha$)")
        plt.ylabel("Normalized Debye intensity")
        plt.title(Path(args.input).name)
        plt.xlim(10, 65)
        plt.tight_layout(); plt.savefig(p, dpi=300)
    print(f"Wrote {out}")
if __name__ == "__main__": main()
