#!/bin/bash
#SBATCH --job-name=chgnet_gpu
#SBATCH --partition=gpu          # Change to your cluster's GPU partition name if different
#SBATCH --gres=gpu:1             # Request 1 GPU
#SBATCH --ntasks=1               # Single task (python script)
#SBATCH --cpus-per-task=4        # Number of CPU cores for data loading/processing
#SBATCH --time=04:00:00          # 4 hours limit (simulation should take <1 hour on GPU)
#SBATCH --output=logs/chgnet_gpu_%j.out
#SBATCH --error=logs/chgnet_gpu_%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

set -euo pipefail

echo "Starting CHGNet GPU simulation job..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CUDA Devices: $CUDA_VISIBLE_DEVICES"

# Load CUDA module if required by your HPC cluster
# module load cuda

# Activate conda environment (adjust conda installation path if different)
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate rki
    echo "Activated conda environment: rki"
else
    echo "Warning: conda.sh not found. Attempting to use current Python environment."
fi

# Run the GPU-enabled molecular dynamics simulation
python3 scripts/custom/read_data_gpu.py

echo "Job finished successfully."
echo "End Date: $(date)"
