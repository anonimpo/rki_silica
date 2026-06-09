#!/usr/bin/env bash
set -euo pipefail
LMP=${LAMMPS_BIN:-lmp}
mkdir -p logs outputs structures
"$LMP" -in lammps/in.03_reaxff_plumed_wtmetad | tee logs/03_wtmetad.log
