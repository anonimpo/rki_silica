#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

pairs = [
    ("before_calcination_pH7", ROOT / "kode_lammps/structures/precursor_pH7p0_N192.data"),
    ("after_calcination_1173K_pH7", ROOT / "runs/pH7p0_N192/03_calcination_1173K/calcined_1173K.data"),
]

csvs = []
for label, path in pairs:
    if not path.exists():
        print(f"Missing {path}; run LAMMPS stage first.")
        continue
    csv = OUT / f"xrd_{label}.csv"
    subprocess.check_call([sys.executable, str(ROOT / "post_processing/compute_debye_xrd.py"), "--input", str(path), "--format", "data", "--out", str(csv)])
    csvs.append((label, csv))

if csvs:
    plt.figure(figsize=(7.2, 4.2))
    for label, csv in csvs:
        df = pd.read_csv(csv)
        plt.plot(df["two_theta_deg"], df["I_norm"], linewidth=1.2, label=label.replace("_", " "))
    plt.xlabel(r"2$\theta$ (degree, Cu K$\alpha$)")
    plt.ylabel("Normalized Debye intensity")
    plt.title("pH 7, N192: XRD before vs after calcination")
    plt.xlim(10, 65)
    plt.legend(frameon=False)
    plt.tight_layout()
    fig = OUT / "xrd_before_after_pH7_N192.png"
    plt.savefig(fig, dpi=300)
    print(f"Wrote {fig}")
