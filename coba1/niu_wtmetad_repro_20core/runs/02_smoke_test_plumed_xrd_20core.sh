#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LMP=${LMP:-lmp}
NCORE=${NCORE:-20}
DT=${DT:-0.25}
DATA=${DATA:-structures/beta_cristobalite_a7.15_N1536.data}
POTENTIAL=${POTENTIAL:-potentials/ffield.reax.SiOH}
mkdir -p logs outputs structures
mpirun -np "$NCORE" "$LMP" -in input_lammps/in.00_plumed_xrd_smoke_test \
  -var DT "$DT" -var DATA "$DATA" -var POTENTIAL "$POTENTIAL" | tee logs/00_smoke_test.log
