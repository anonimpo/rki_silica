#!/bin/bash
#SBATCH --job-name=ph_silica_reaxff
#SBATCH --array=0-4
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

# Edit these for your cluster.
# module load lammps
# module load plumed

set -euo pipefail
PH_LIST=(pH6p0 pH6p5 pH7p0 pH7p5 pH8p0)
PH=${PH_LIST[$SLURM_ARRAY_TASK_ID]}
N=${N:-1536}
TPROD=${TPROD:-2300}
export LMP=${LMP:-lmp}
export N_EQ=${N_EQ:-200000}
export N_HEAT=${N_HEAT:-400000}
export N_HOLD=${N_HOLD:-800000}
export N_PROD=${N_PROD:-2000000}
export RUN_META=${RUN_META:-0}

./scripts/02_run_one_pipeline.sh "$PH" "$N"
