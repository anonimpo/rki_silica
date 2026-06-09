#!/usr/bin/env bash
set -euo pipefail
LMP=${LAMMPS_BIN:-lmp}
mkdir -p logs outputs structures
"$LMP" -in lammps/in.04_reaxff_production_no_plumed | tee logs/04_production_no_plumed.log
python3 scripts/xrd_debye_cv.py outputs/production_no_plumed.lammpstrj --stride 1 > outputs/xrd_cv_from_production.dat
echo "Wrote outputs/xrd_cv_from_production.dat"
