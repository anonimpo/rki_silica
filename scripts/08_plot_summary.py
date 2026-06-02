#!/usr/bin/env python3
"""Aggregate structural metrics and XRD curves across pH cases."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PH_ORDER = ["pH6p0", "pH6p5", "pH7p0", "pH7p5", "pH8p0"]
PH_VALUE = {"pH6p0": 6.0, "pH6p5": 6.5, "pH7p0": 7.0, "pH7p5": 7.5, "pH8p0": 8.0}


def read_metrics(metrics_dir: Path) -> pd.DataFrame:
    rows = []
    for f in metrics_dir.glob("metrics_*.csv"):
        rows.append(pd.read_csv(f))
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    df["pH"] = df["ph_label"].map(PH_VALUE)
    return df.sort_values("pH")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    analysis_dir = Path(args.analysis_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_metrics(analysis_dir / "metrics")
    if not df.empty:
        df.to_csv(analysis_dir / "metrics_all.csv", index=False)

        plt.figure(figsize=(6.5, 4.2))
        plt.plot(df["pH"], df["frac_si_coord4"], marker="o")
        plt.xlabel("pH proxy")
        plt.ylabel("Fraction of tetrahedral Si (SiO4)")
        plt.title("MD metric: keteraturan tetrahedral SiO4")
        plt.ylim(0, 1.05)
        plt.tight_layout()
        plt.savefig(outdir / "fig_md_tetrahedral_si_fraction_vs_pH.png", dpi=300)

        plt.figure(figsize=(6.5, 4.2))
        plt.plot(df["pH"], df["frac_o_bridging_2si_0h"], marker="o", label="Bridging O / siloxane")
        plt.plot(df["pH"], df["frac_o_silanol_proxy_1si_ge1h"], marker="s", label="Silanol proxy")
        plt.plot(df["pH"], df["frac_o_water_proxy_0si_ge2h"], marker="^", label="Water proxy")
        plt.xlabel("pH proxy")
        plt.ylabel("Fraction of O atoms")
        plt.title("MD metric: siloxane, silanol, dan air terikat")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "fig_md_oxygen_speciation_vs_pH.png", dpi=300)

    # XRD overlay.
    xrd_files = sorted((analysis_dir / "xrd").glob("xrd_*.csv"))
    if xrd_files:
        plt.figure(figsize=(7.0, 4.5))
        offset = 0.0
        for f in xrd_files:
            label = next((p for p in PH_ORDER if p in f.name), f.stem)
            x = pd.read_csv(f)
            plt.plot(x["two_theta_deg"], x["I_norm"] + offset, linewidth=1.0, label=label)
            offset += 1.05
        plt.xlabel(r"2$\theta$ (degree, Cu K$\alpha$)")
        plt.ylabel("Normalized intensity + offset")
        plt.title("MD Debye-XRD overlay across pH proxies")
        plt.xlim(10, 60)
        plt.legend(ncol=2, fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / "fig_md_xrd_overlay_vs_pH.png", dpi=300)

    print(f"Wrote summaries to {outdir}")


if __name__ == "__main__":
    main()
