#!/usr/bin/env bash
set -euo pipefail

# Run from the project root: ./kode_lammps/run_pH7_N192_20core.sh
# Override executable if needed, e.g. LMP=lmp_mpi ./kode_lammps/run_pH7_N192_20core.sh

LMP=${LMP:-lmp}
NP=${NP:-20}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

run_lmp() {
  echo "[RUN] mpirun -np ${NP} ${LMP} $*"
  mpirun -np "${NP}" "${LMP}" "$@"
}

PH=pH7p0
BASE=192

# Short trial settings. Increase these for a real production run.
run_lmp -var phlabel ${PH} -var baseatoms ${BASE} \
  -in kode_lammps/in.01_minimize_precursor

run_lmp -var phlabel ${PH} -var baseatoms ${BASE} -var nsteps 40000 \
  -in kode_lammps/in.02_equilibrate_300K

run_lmp -var phlabel ${PH} -var baseatoms ${BASE} -var nheat 80000 -var nhold 160000 \
  -in kode_lammps/in.03_calcination_1173K

run_lmp -var phlabel ${PH} -var baseatoms ${BASE} -var Tprod 2300 -var nheat 80000 -var nprod 240000 \
  -in kode_lammps/in.04_highT_unbiased_2300K

# Optional PLUMED stage; comment this block if PLUMED is not compiled into LAMMPS.
run_lmp -var phlabel ${PH} -var baseatoms ${BASE} -var Tprod 2300 \
  -var nheat 80000 -var nmeta 240000 \
  -var plumedfile kode_lammps/plumed/plumed_pH7p0_N192_localq6.dat \
  -in kode_lammps/in.05_wtmetad_plumed_FIXED
