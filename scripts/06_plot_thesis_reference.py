#!/usr/bin/env python3
"""Plot thesis reference data used as experimental anchors for MD interpretation."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

PH = np.array([6.0, 6.5, 7.0, 7.5, 8.0])
IK_IT = np.array([0.639, 0.834, 0.909, 0.656, 0.707])
FWHM_K = np.array([0.3958, 0.3684, 0.3349, 0.4136, 0.4131])
SIZE_K = np.array([34, 38, 46, 41, 36])
FWHM_T = np.array([0.3128, 0.3345, 0.2803, 0.3618, 0.3669])
SIZE_T = np.array([48, 44, 59, 58, 52])


def main() -> None:
    outdir = Path("outputs/thesis_reference_plots")
    outdir.mkdir(parents=True, exist_ok=True)

    # Save table.
    table = np.column_stack([PH, IK_IT, FWHM_K, SIZE_K, FWHM_T, SIZE_T])
    np.savetxt(outdir / "thesis_xrd_reference_table.csv", table, delimiter=",",
               header="pH,Ik_over_It,FWHM_cristobalite,size_cristobalite_nm,FWHM_tridymite,size_tridymite_nm", comments="")

    plt.figure(figsize=(6.5, 4.2))
    plt.plot(PH, SIZE_K, marker="o", label="Cristobalite")
    plt.plot(PH, SIZE_T, marker="s", label="Tridymite")
    plt.xlabel("pH sintesis sol-gel")
    plt.ylabel("Ukuran kristalit Scherrer (nm)")
    plt.title("Acuan eksperimen: ukuran kristalit setelah kalsinasi 900 °C")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_thesis_crystallite_size_vs_pH.png", dpi=300)

    plt.figure(figsize=(6.5, 4.2))
    plt.plot(PH, IK_IT, marker="o")
    plt.xlabel("pH sintesis sol-gel")
    plt.ylabel("Rasio intensitas Ik/It")
    plt.title("Acuan eksperimen: dominasi cristobalite relatif terhadap tridymite")
    plt.tight_layout()
    plt.savefig(outdir / "fig_thesis_Ikratio_vs_pH.png", dpi=300)

    plt.figure(figsize=(6.5, 4.2))
    plt.plot(PH, FWHM_K, marker="o", label="FWHM cristobalite")
    plt.plot(PH, FWHM_T, marker="s", label="FWHM tridymite")
    plt.xlabel("pH sintesis sol-gel")
    plt.ylabel("FWHM (degree 2θ)")
    plt.title("Acuan eksperimen: pelebaran puncak XRD setelah kalsinasi")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_thesis_fwhm_vs_pH.png", dpi=300)

    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
