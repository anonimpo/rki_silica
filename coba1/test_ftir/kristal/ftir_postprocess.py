#!/usr/bin/env python3
"""
ftir_postprocess.py

Post-process dipole time series from MD into a relative FTIR spectrum.

Input
- dipoles.csv from run_md_chgnet.py with columns:
  step,time_fs,Mx_eA,My_eA,Mz_eA

Method
- Uses the dipole-derivative power spectrum.
- This gives a relative IR spectrum from M(t) = sum_i q_i r_i(t).
- For absolute intensities, use a validated dipole model and prefactors for your material.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Tuple

import numpy as np

C_CM_S = 2.99792458e10
H_J_S = 6.62607015e-34
KB_J_K = 1.380649e-23


def load_dipoles_csv(path: str | Path) -> Tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError("Dipole file is empty")
    time_fs = np.asarray(data["time_fs"], dtype=float)
    M = np.column_stack([
        np.asarray(data["Mx_eA"], dtype=float),
        np.asarray(data["My_eA"], dtype=float),
        np.asarray(data["Mz_eA"], dtype=float),
    ])
    if M.ndim != 2 or M.shape[1] != 3:
        raise ValueError("Dipole array must have shape N x 3")
    if len(time_fs) < 8:
        raise ValueError("Need at least 8 frames for a spectrum")
    return time_fs, M


def detrend_columns(x: np.ndarray) -> np.ndarray:
    """Remove linear trend from each component without requiring scipy."""
    n = x.shape[0]
    t = np.arange(n, dtype=float)
    out = np.empty_like(x, dtype=float)
    for j in range(x.shape[1]):
        coeff = np.polyfit(t, x[:, j], deg=1)
        out[:, j] = x[:, j] - np.polyval(coeff, t)
    return out


def gaussian_smooth(y: np.ndarray, sigma_points: float) -> np.ndarray:
    if sigma_points <= 0:
        return y
    try:
        from scipy.ndimage import gaussian_filter1d
        return gaussian_filter1d(y, sigma=sigma_points, mode="nearest")
    except Exception:
        radius = max(1, int(4 * sigma_points + 0.5))
        grid = np.arange(-radius, radius + 1, dtype=float)
        kernel = np.exp(-0.5 * (grid / sigma_points) ** 2)
        kernel /= kernel.sum()
        return np.convolve(y, kernel, mode="same")


def welch_spectrum(signal: np.ndarray, dt_s: float, segment_points: int, overlap: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return frequency in Hz and summed 3D power spectrum."""
    fs_hz = 1.0 / dt_s
    n = signal.shape[0]
    nperseg = int(max(8, min(segment_points, n)))
    noverlap = int(np.clip(overlap, 0.0, 0.95) * nperseg)

    try:
        from scipy.signal import welch
        p_sum = None
        f_hz = None
        for j in range(signal.shape[1]):
            f_hz, pxx = welch(
                signal[:, j],
                fs=fs_hz,
                window="hann",
                nperseg=nperseg,
                noverlap=noverlap,
                detrend="constant",
                scaling="density",
                return_onesided=True,
            )
            p_sum = pxx if p_sum is None else p_sum + pxx
        return f_hz, p_sum
    except Exception:
        step = max(1, nperseg - noverlap)
        window = np.hanning(nperseg)
        norm = fs_hz * np.sum(window**2)
        acc = None
        count = 0
        for start in range(0, n - nperseg + 1, step):
            block = signal[start:start + nperseg]
            block = block - block.mean(axis=0, keepdims=True)
            block = block * window[:, None]
            fft = np.fft.rfft(block, axis=0)
            pxx = (np.abs(fft) ** 2).sum(axis=1) / norm
            if nperseg % 2 == 0:
                pxx[1:-1] *= 2.0
            else:
                pxx[1:] *= 2.0
            acc = pxx if acc is None else acc + pxx
            count += 1
        if count == 0:
            raise ValueError("Not enough frames for selected segment length")
        f_hz = np.fft.rfftfreq(nperseg, d=dt_s)
        return f_hz, acc / count


def compute_ftir(
    time_fs: np.ndarray,
    M_eA: np.ndarray,
    temperature_K: float,
    drop_first_ps: float,
    segment_ps: float,
    overlap: float,
    smooth_cm1: float,
    min_cm1: float,
    max_cm1: float,
    quantum_correction: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    keep = time_fs >= (time_fs[0] + 1000.0 * drop_first_ps)
    time_fs = time_fs[keep]
    M_eA = M_eA[keep]
    if len(time_fs) < 16:
        raise ValueError("Too few frames after drop-first-ps")

    dt_fs_values = np.diff(time_fs)
    dt_fs = float(np.median(dt_fs_values))
    if not np.allclose(dt_fs_values, dt_fs, rtol=1e-4, atol=1e-8):
        raise ValueError("Time spacing is not uniform enough for FFT/Welch")
    dt_s = dt_fs * 1e-15

    M = detrend_columns(M_eA)
    dM_dt = np.gradient(M, dt_s, axis=0)
    dM_dt -= dM_dt.mean(axis=0, keepdims=True)

    segment_points = int(round((segment_ps * 1000.0) / dt_fs))
    segment_points = max(16, segment_points)
    f_hz, power = welch_spectrum(dM_dt, dt_s, segment_points, overlap)
    wn_cm1 = f_hz / C_CM_S

    intensity = np.asarray(power, dtype=float)
    if quantum_correction:
        x = (H_J_S * C_CM_S * wn_cm1) / (KB_J_K * temperature_K)
        qc = np.ones_like(x)
        nonzero = x > 1e-12
        qc[nonzero] = x[nonzero] / (1.0 - np.exp(-x[nonzero]))
        intensity *= qc

    if len(wn_cm1) > 1 and smooth_cm1 > 0:
        dwn = float(np.median(np.diff(wn_cm1)))
        sigma_points = (smooth_cm1 / 2.354820045) / max(dwn, 1e-30)
        intensity = gaussian_smooth(intensity, sigma_points)

    mask = (wn_cm1 >= min_cm1) & (wn_cm1 <= max_cm1)
    wn_cm1 = wn_cm1[mask]
    intensity = intensity[mask]
    if intensity.size == 0:
        raise ValueError("No spectrum points inside requested cm-1 range")
    max_i = float(np.max(intensity))
    if max_i > 0:
        intensity = intensity / max_i
    return wn_cm1, intensity


def write_spectrum_csv(path: str | Path, wn_cm1: np.ndarray, intensity: np.ndarray) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["wavenumber_cm-1", "relative_intensity"])
        for x, y in zip(wn_cm1, intensity):
            w.writerow([f"{x:.8f}", f"{y:.12e}"])


def plot_spectrum(path: str | Path, wn_cm1: np.ndarray, intensity: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.5, 4.5), dpi=160)
    ax = fig.add_subplot(111)
    ax.plot(wn_cm1, intensity, linewidth=1.3)
    ax.set_xlabel("Wavenumber / cm$^{-1}$")
    ax.set_ylabel("Relative intensity")
    ax.set_xlim(float(wn_cm1.min()), float(wn_cm1.max()))
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_test_dipoles(path: str | Path, dt_fs: float = 0.5, total_ps: float = 5.0) -> None:
    """Create synthetic dipoles with peaks near 1000 and 3400 cm-1."""
    n = int(round(total_ps * 1000.0 / dt_fs)) + 1
    t_fs = np.arange(n) * dt_fs
    t_s = t_fs * 1e-15
    f1 = 1000.0 * C_CM_S
    f2 = 3400.0 * C_CM_S
    Mx = np.sin(2 * np.pi * f1 * t_s)
    My = 0.4 * np.sin(2 * np.pi * f2 * t_s + 0.3)
    Mz = 0.2 * np.sin(2 * np.pi * f1 * t_s + 1.1)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["step", "time_fs", "Mx_eA", "My_eA", "Mz_eA"])
        for i in range(n):
            w.writerow([i, f"{t_fs[i]:.8f}", f"{Mx[i]:.12e}", f"{My[i]:.12e}", f"{Mz[i]:.12e}"])


def main(args) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.self_test:
        dipole_path = out_dir / "synthetic_dipoles.csv"
        make_test_dipoles(dipole_path, dt_fs=0.5, total_ps=8.0)
    else:
        dipole_path = Path(args.dipoles)

    time_fs, M = load_dipoles_csv(dipole_path)
    wn, intensity = compute_ftir(
        time_fs=time_fs,
        M_eA=M,
        temperature_K=args.temperature_K,
        drop_first_ps=args.drop_first_ps,
        segment_ps=args.segment_ps,
        overlap=args.overlap,
        smooth_cm1=args.smooth_cm1,
        min_cm1=args.min_cm1,
        max_cm1=args.max_cm1,
        quantum_correction=args.quantum_correction,
    )

    spectrum_path = out_dir / args.spectrum_name
    plot_path = out_dir / args.plot_name
    write_spectrum_csv(spectrum_path, wn, intensity)
    plot_spectrum(plot_path, wn, intensity)

    top = np.argsort(intensity)[-8:][::-1]
    print(f"Wrote: {spectrum_path}")
    print(f"Wrote: {plot_path}")
    print("Top peak candidates, cm-1:")
    for idx in top:
        print(f"  {wn[idx]:10.3f}  {intensity[idx]:.4f}")


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert dipole time series to relative FTIR spectrum.")
    p.add_argument("--dipoles", default="run_chgnet_md/dipoles.csv", help="Dipole CSV from MD.")
    p.add_argument("--out-dir", default="ftir_out", help="Output directory.")
    p.add_argument("--spectrum-name", default="ftir_spectrum.csv")
    p.add_argument("--plot-name", default="ftir_spectrum.png")
    p.add_argument("--temperature-K", type=float, default=300.0)
    p.add_argument("--drop-first-ps", type=float, default=0.0, help="Discard initial time before FTIR analysis.")
    p.add_argument("--segment-ps", type=float, default=4.0, help="Welch segment length in ps. Longer gives finer resolution, shorter gives smoother spectra.")
    p.add_argument("--overlap", type=float, default=0.5, help="Welch segment overlap fraction.")
    p.add_argument("--smooth-cm1", type=float, default=20.0, help="Gaussian FWHM smoothing in cm-1. Set 0 to disable.")
    p.add_argument("--min-cm1", type=float, default=0.0)
    p.add_argument("--max-cm1", type=float, default=4500.0)
    p.add_argument("--quantum-correction", action="store_true", help="Apply harmonic quantum correction factor.")
    p.add_argument("--self-test", action="store_true", help="Generate synthetic dipoles with peaks near 1000 and 3400 cm-1.")
    return p


if __name__ == "__main__":
    main(make_parser().parse_args())
