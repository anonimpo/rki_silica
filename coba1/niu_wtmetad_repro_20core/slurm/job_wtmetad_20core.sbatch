#!/bin/bash
#SBATCH --job-name=niu_xrd_wtmetad
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# module load lammps plumed   # edit for your cluster
export NCORE=20
export LMP=${LMP:-lmp}
export DT=${DT:-0.25}
export DUMP_STRIDE=${DUMP_STRIDE:-160000}
export RUN_STEPS=${RUN_STEPS:-960000}

./runs/00_check_environment.sh
./runs/02_smoke_test_plumed_xrd_20core.sh
./runs/03_make_liquid_from_beta_20core.sh
./runs/05_run_wtmetad_xrd111_20core.sh
./runs/07_postprocess.sh
