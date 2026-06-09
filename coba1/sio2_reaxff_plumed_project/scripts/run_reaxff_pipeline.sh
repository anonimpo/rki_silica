#!/usr/bin/env bash
set -euo pipefail
LMP=${LAMMPS_BIN:-lmp}
mkdir -p logs outputs structures

"$LMP" -in lammps/in.00_minimize_reaxff | tee logs/00_minimize.log
"$LMP" -in lammps/in.01_melt_reaxff     | tee logs/01_melt.log
"$LMP" -in lammps/in.02_equilibrate_liquid_reaxff | tee logs/02_equilibrate.log

echo "ReaxFF preparation done. Next options:"
echo "  1) LAMMPS+PLUMED executable: $LMP -in lammps/in.03_reaxff_plumed_wtmetad"
echo "  2) ReaxFF-only production:    $LMP -in lammps/in.04_reaxff_production_no_plumed"
