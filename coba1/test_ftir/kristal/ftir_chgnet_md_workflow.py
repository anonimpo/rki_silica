#!/usr/bin/env python3
"""
ftir_chgnet_md_workflow.py

Python workflow hasil konversi dari notebook ftir_chgnet_md_workflow.ipynb.

File ini tidak otomatis menjalankan MD 50.000 step saat dibuka. Gunakan argumen CLI:

    python ftir_chgnet_md_workflow.py --inspect
    python ftir_chgnet_md_workflow.py --run-md
    python ftir_chgnet_md_workflow.py --run-ftir
    python ftir_chgnet_md_workflow.py --self-test
    python ftir_chgnet_md_workflow.py --run-md --run-ftir

Instalasi dependency jika belum tersedia:

    pip install ase chgnet torch scipy matplotlib

Catatan:
- CHGNet digunakan untuk energi dan gaya pada MD.
- FTIR dihitung dari sinyal dipol M(t) = sum(q_i r_i).
- Karena CHGNet tidak memberi dipol IR langsung, muatan dipol berasal dari charge_map, charge_file, atau muatan input.
"""

from __future__ import annotations

import argparse
from argparse import Namespace
from pathlib import Path

from run_md_chgnet import build_atoms, composition_dict, run_md
from ftir_postprocess import main as run_ftir_postprocess


BASE_DIR = Path.cwd()
INPUT_FILE = Path("cristobalite_1536.data")
MD_OUT_DIR = BASE_DIR / "run_pH7_chgnet"
FTIR_OUT_DIR = BASE_DIR / "ftir_pH7"


def make_md_args() -> Namespace:
    return Namespace(
        input=str(INPUT_FILE),
        format="lammps-data",
        atom_style="charge",
        z_of_type="1:Si,2:O,3:H",
        out_dir=str(MD_OUT_DIR),
        steps=50_000,
        timestep_fs=0.5,
        temperature_K=300.0,
        langevin_tau_fs=100.0,
        sample_interval=1,
        traj_interval=100,
        print_interval=100,
        chunk_steps=1_000,
        relax_steps=200,
        relax_fmax=0.05,
        device="cuda",          # ganti ke "cpu" jika tidak ada GPU
        model_path=None,
        check_cuda_mem=False,
        on_isolated_atoms="warn",
        torch_threads=8,
        matmul_precision="high",
        charge_map="Si:2.4,O:-1.2,H:0.6",
        charge_file=None,
        neutralize=False,
        net_charge_tol=1e-6,
        allow_zero_charges=False,
        pbc_x=True,
        pbc_y=True,
        pbc_z=True,
        seed=12345,
    )


def make_ftir_args() -> Namespace:
    return Namespace(
        dipoles=str(MD_OUT_DIR / "dipoles.csv"),
        out_dir=str(FTIR_OUT_DIR),
        spectrum_name="ftir_spectrum.csv",
        plot_name="ftir_spectrum.png",
        temperature_K=300.0,
        drop_first_ps=2.0,
        segment_ps=8.0,
        overlap=0.5,
        smooth_cm1=20.0,
        min_cm1=100.0,
        max_cm1=4500.0,
        quantum_correction=True,
        self_test=False,
    )


def make_self_test_args() -> Namespace:
    return Namespace(
        dipoles="",
        out_dir=str(BASE_DIR / "ftir_test_notebook"),
        spectrum_name="ftir_spectrum.csv",
        plot_name="ftir_spectrum.png",
        temperature_K=300.0,
        drop_first_ps=0.0,
        segment_ps=4.0,
        overlap=0.5,
        smooth_cm1=20.0,
        min_cm1=100.0,
        max_cm1=4500.0,
        quantum_correction=True,
        self_test=True,
    )


def inspect_structure(md_args: Namespace) -> None:
    print("BASE_DIR:", BASE_DIR)
    print("INPUT_FILE:", Path(md_args.input))
    print("Input exists:", Path(md_args.input).exists())

    try:
        atoms_preview = build_atoms(md_args)
        input_charges = atoms_preview.get_initial_charges()

        print("n_atoms:", len(atoms_preview))
        print("composition:", composition_dict(atoms_preview.get_chemical_symbols()))
        print("cell_A:\n", atoms_preview.cell)
        print("pbc:", atoms_preview.pbc)
        print("input_charge_sum:", float(input_charges.sum()))
        print("input_charge_min_max:", float(input_charges.min()), float(input_charges.max()))
    except Exception as exc:
        print("Inspeksi struktur belum bisa dijalankan:", exc)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Workflow MD CHGNet dan post-processing FTIR.")
    parser.add_argument("--inspect", action="store_true", help="Cek input struktur sebelum MD.")
    parser.add_argument("--run-md", action="store_true", help="Jalankan simulasi MD utama.")
    parser.add_argument("--run-ftir", action="store_true", help="Jalankan post-processing FTIR dari dipoles.csv.")
    parser.add_argument("--self-test", action="store_true", help="Jalankan uji post-processing FTIR dengan dipol sintetis.")
    parser.add_argument("--input", default=str(INPUT_FILE), help="Path file input struktur, misalnya /content/calcined_1173K.data.")
    parser.add_argument("--device", default="cuda", choices=["cpu", "cuda", "mps"], help="Device untuk CHGNet.")
    parser.add_argument("--steps", type=int, default=50_000, help="Jumlah step MD.")
    return parser


def main() -> None:
    cli_args = make_parser().parse_args()

    md_args = make_md_args()
    md_args.input = cli_args.input
    md_args.device = cli_args.device
    md_args.steps = cli_args.steps

    if cli_args.inspect:
        inspect_structure(md_args)

    if cli_args.run_md:
        run_md(md_args)

    if cli_args.run_ftir:
        run_ftir_postprocess(make_ftir_args())

    if cli_args.self_test:
        run_ftir_postprocess(make_self_test_args())

    if not any([cli_args.inspect, cli_args.run_md, cli_args.run_ftir, cli_args.self_test]):
        print("Tidak ada aksi dijalankan. Gunakan --inspect, --run-md, --run-ftir, atau --self-test.")


if __name__ == "__main__":
    main()
