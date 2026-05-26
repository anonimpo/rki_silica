#!/usr/bin/env bash
#SBATCH --job-name=sio2_crystallization
#SBATCH --output=logs/sio2_crystallization_%j.out
#SBATCH --error=logs/sio2_crystallization_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --partition=standard

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# ==============================================================================
# Environment Configuration
# ==============================================================================
# Load necessary modules for your HPC environment. Examples:
# module load intel/2021.4.0
# module load openmpi/4.1.2
# module load python/3.9.7

# Define the LAMMPS binary.
# Set this to your compiled LAMMPS executable. If running Path A (Metadynamics),
# this binary MUST be built with both REAXFF and PLUMED packages.
LMP="lmp"

# MPI command prefix (e.g., srun, mpirun, mpiexec)
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"

# Select the crystallization pathway:
# "A" = Well-Tempered Metadynamics (requires LAMMPS compiled with PLUMED)
# "B" = Production ReaxFF + Offline Python XRD post-processing
PATHWAY="B"

# ==============================================================================
# Workflow Execution
# ==============================================================================

echo "======================================================================"
echo "Starting SiO2 Crystallization Workflow"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Running on ${SLURM_NNODES:-1} node(s) with ${SLURM_NTASKS:-16} MPI task(s)"
echo "Selected Pathway: Pathway ${PATHWAY}"
echo "======================================================================"

# 1. Ensure output and log directories exist
mkdir -p logs outputs structures

# 2. Step 0: Generate initial crystal structure
echo "[Stage 0/4] Generating beta-cristobalite initial structure..."
python3 scripts/make_beta_cristobalite.py --a 7.15 --reps 2 2 2 --out-prefix structures/beta_cristobalite_192

# 3. Step 1: Energy Minimization
echo "[Stage 1/4] Running ReaxFF minimization..."
${MPI_RUN} ${LMP} -in lammps/in.00_minimize_reaxff | tee logs/00_minimize.log

# 4. Step 2: High-temperature Melt (4000 K)
echo "[Stage 2/4] Running ReaxFF melt simulation..."
${MPI_RUN} ${LMP} -in lammps/in.01_melt_reaxff | tee logs/01_melt.log

# 5. Step 3: Liquid Equilibration (2400 K)
echo "[Stage 3/4] Running ReaxFF equilibration..."
${MPI_RUN} ${LMP} -in lammps/in.02_equilibrate_liquid_reaxff | tee logs/02_equilibrate.log

# 6. Step 4: Production / Crystallization Stage
if [ "${PATHWAY}" = "A" ]; then
    echo "[Stage 4/4] Running WELL-TEMPERED METADYNAMICS (Pathway A)..."
    echo "Using PLUMED input: plumed/plumed_saxs_metad_approx.dat"
    # Ensure PLUMED kernel can be detected if shared library is used
    # export PLUMED_KERNEL=/path/to/libplumedKernel.so
    
    ${MPI_RUN} ${LMP} -in lammps/in.03_reaxff_plumed_wtmetad | tee logs/03_wtmetad.log
    echo "Pathway A simulation finished."
    
elif [ "${PATHWAY}" = "B" ]; then
    echo "[Stage 4/4] Running REAXFF PRODUCTION (Pathway B)..."
    ${MPI_RUN} ${LMP} -in lammps/in.04_reaxff_production_no_plumed | tee logs/04_production_no_plumed.log
    
    echo "[Stage 4/4] Post-processing trajectory for Debye XRD CVs..."
    python3 scripts/xrd_debye_cv.py outputs/production_no_plumed.lammpstrj --stride 1 > outputs/xrd_cv_from_production.dat
    echo "Pathway B simulation and offline post-processing finished."
    echo "XRD CV output written to outputs/xrd_cv_from_production.dat"
else
    echo "ERROR: Invalid PATHWAY selection: ${PATHWAY}. Choose 'A' or 'B'." >&2
    exit 1
fi

echo "======================================================================"
echo "SiO2 Crystallization Workflow Completed Successfully"
echo "======================================================================"
