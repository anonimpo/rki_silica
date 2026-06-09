#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
TRAJ = ROOT / "runs/pH7p0_N192/04_highT_unbiased_2300K/traj_highT_unbiased.lammpstrj"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

def read_frames(path):
    lines = path.read_text(errors="ignore").splitlines()
    idx = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    for i in idx:
        step = int(lines[i + 1]); n = int(lines[i + 3])
        bounds = np.array([[float(x) for x in lines[i + 5 + j].split()[:2]] for j in range(3)])
        L = bounds[:, 1] - bounds[:, 0]
        header = lines[i + 8].split()[2:]; col = {v: k for k, v in enumerate(header)}
        types = np.empty(n, int); coords = np.empty((n, 3), float)
        for r, line in enumerate(lines[i + 9:i + 9 + n]):
            p = line.split(); types[r] = int(p[col["type"]])
            coords[r] = [float(p[col["x"]]), float(p[col["y"]]), float(p[col["z"]])] - bounds[:, 0]
        yield step, types, coords % L, L

def pbc_dist(a, b, L):
    d = a[:, None, :] - b[None, :, :]
    d -= np.round(d / L) * L
    return np.linalg.norm(d, axis=2)

def order_proxy(types, coords, L):
    si = coords[types == 1]
    o = coords[types == 2]
    d_so = pbc_dist(si, o, L)
    coord4 = np.sum(d_so < 2.0, axis=1) == 4
    d_ss = pbc_dist(si, si, L)
    np.fill_diagonal(d_ss, np.inf)
    nearest4 = np.sort(d_ss, axis=1)[:, :4]
    # beta-cristobalite-like Si-Si first-neighbor distance is roughly around 3.1 A near this scale.
    dev = np.mean(np.abs(nearest4 - 3.1), axis=1)
    dev_norm = np.clip(dev / 1.2, 0, 1)
    # Map to Niu-style color range only for visualization.
    sbar = -1.75 + 0.78 * dev_norm
    solid_like = coord4 & (dev < 1.1)
    return si, sbar, solid_like

if not TRAJ.exists():
    raise SystemExit(f"Missing trajectory: {TRAJ}. Run stage 04 first.")
frames = list(read_frames(TRAJ))
sel = np.linspace(0, len(frames) - 1, 6, dtype=int)
letters = list("ABCDEF")
fig = plt.figure(figsize=(11, 7))
for k, idx in enumerate(sel):
    step, types, coords, L = frames[idx]
    si, sbar, solid_like = order_proxy(types, coords, L)
    ax = fig.add_subplot(2, 3, k + 1, projection="3d")
    pts = si[solid_like]
    vals = sbar[solid_like]
    if len(pts) == 0:
        pts = si; vals = sbar
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=vals, s=18, vmin=-1.75, vmax=-0.97, alpha=0.85)
    # Draw simple simulation box.
    for a in [0, L[0]]:
        for b in [0, L[1]]:
            ax.plot([a, a], [b, b], [0, L[2]], linewidth=0.5)
    for a in [0, L[0]]:
        for c in [0, L[2]]:
            ax.plot([a, a], [0, L[1]], [c, c], linewidth=0.5)
    for b in [0, L[1]]:
        for c in [0, L[2]]:
            ax.plot([0, L[0]], [b, b], [c, c], linewidth=0.5)
    ax.text2D(0.02, 0.92, letters[k], transform=ax.transAxes, fontsize=18, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=16, azim=-65)
    ax.set_title(f"step {step}", fontsize=8)
fig.suptitle("Figure-7-like snapshots: crystal-like Si atoms, pH 7 N192", y=0.98)
plt.tight_layout()
out = OUT / "figure7_like_pH7_N192.png"
plt.savefig(out, dpi=300)
print(f"Wrote {out}")
