#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   LMP=lmp ./scripts/02_run_one_pipeline.sh pH7p0 1536
# Optional environment variables:
#   TPROD=2300 N_EQ=20000 N_HEAT=80000 N_HOLD=120000 N_PROD=400000 N_META=400000 RUN_META=0 DRYRUN=0

PH=${1:-pH7p0}
N=${2:-1536}
LMP=${LMP:-lmp}
TPROD=${TPROD:-2300}
RUN_META=${RUN_META:-0}
DRYRUN=${DRYRUN:-0}

run_cmd() {
  echo "+ $*"
  if [[ "$DRYRUN" != "1" ]]; then
    "$@"
  fi
}

run_cmd "$LMP" -var phlabel "$PH" -var baseatoms "$N" \
  -in lammps/in.01_minimize_precursor

run_cmd "$LMP" -var phlabel "$PH" -var baseatoms "$N" -var nsteps "${N_EQ:-20000}" \
  -in lammps/in.02_equilibrate_300K

run_cmd "$LMP" -var phlabel "$PH" -var baseatoms "$N" -var nheat "${N_HEAT:-80000}" -var nhold "${N_HOLD:-120000}" \
  -in lammps/in.03_calcination_1173K

run_cmd "$LMP" -var phlabel "$PH" -var baseatoms "$N" -var Tprod "$TPROD" -var nheat "${N_HEAT:-80000}" -var nprod "${N_PROD:-400000}" \
  -in lammps/in.04_highT_unbiased_reaxff

if [[ "$RUN_META" == "1" ]]; then
  PLUMED_FILE="plumed/stock_localq6/plumed_localq6_${PH}_N${N}.dat"
  run_cmd "$LMP" -var phlabel "$PH" -var baseatoms "$N" -var Tprod "$TPROD" \
    -var nheat "${N_HEAT:-80000}" -var nmeta "${N_META:-400000}" -var plumedfile "$PLUMED_FILE" \
    -in lammps/in.05_wtmetad_plumed
fi
