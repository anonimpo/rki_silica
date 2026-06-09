#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LMP=${LMP:-lmp}
NCORE=${NCORE:-20}
DT=${DT:-0.25}
DATA=${DATA:-structures/beta_cristobalite_a7.15_N1536.data}
POTENTIAL=${POTENTIAL:-potentials/ffield.reax.SiOH}
# ReaxFF-safe default: 0.25 fs. DUMP_STRIDE=160000 means 40 ps.
DUMP_STRIDE=${DUMP_STRIDE:-160000}
MELT_STEPS=${MELT_STEPS:-200000}   # 50 ps if DT=0.25 fs
COOL_STEPS=${COOL_STEPS:-200000}   # 50 ps
EQ_STEPS=${EQ_STEPS:-200000}       # 50 ps
mkdir -p logs outputs structures
mpirun -np "$NCORE" "$LMP" -in input_lammps/in.01_melt_beta_make_liquid_2300K \
  -var DT "$DT" -var DATA "$DATA" -var POTENTIAL "$POTENTIAL" \
  -var DUMP_STRIDE "$DUMP_STRIDE" -var MELT_STEPS "$MELT_STEPS" -var COOL_STEPS "$COOL_STEPS" -var EQ_STEPS "$EQ_STEPS" | tee logs/01_melt_make_liquid.log
