#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LMP=${LMP:-lmp}
NCORE=${NCORE:-20}
DT=${DT:-0.25}
DATA=${DATA:-structures/liquid_from_beta_2300K_N1536.data}
POTENTIAL=${POTENTIAL:-potentials/ffield.reax.SiOH}
DUMP_STRIDE=${DUMP_STRIDE:-160000}
RUN_STEPS=${RUN_STEPS:-800000}     # 200 ps if DT=0.25 fs
mkdir -p logs outputs structures
mpirun -np "$NCORE" "$LMP" -in input_lammps/in.02_unbiased_control_2300K \
  -var DT "$DT" -var DATA "$DATA" -var POTENTIAL "$POTENTIAL" \
  -var DUMP_STRIDE "$DUMP_STRIDE" -var RUN_STEPS "$RUN_STEPS" | tee logs/02_unbiased_control.log
