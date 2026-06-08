#!/bin/bash
#SBATCH --job-name=pH6p0_1536_stage_04
#SBATCH --array=0-4
#SBATCH --ntasks=16
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

# Activate conda environment (adjust conda installation path if different)
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate rki
    echo "Activated conda environment: rki"
else
    echo "Warning: conda.sh not found. Attempting to use current Python environment."
fi

lmp -var phlabel pH6p0 -var baseatoms 1536 -var Tprod 2300 -var nheat 80000 -var nprod 400000 -in lammps/in.04_highT_unbiased_reaxff
