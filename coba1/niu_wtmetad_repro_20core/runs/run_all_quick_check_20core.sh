#!/usr/bin/env bash
set -euo pipefail
# Quick functional check. This is short and may NOT crystallize; use longer RUN_STEPS for science.
export NCORE=${NCORE:-20}
export DT=${DT:-0.25}
export MELT_STEPS=${MELT_STEPS:-20000}
export COOL_STEPS=${COOL_STEPS:-20000}
export EQ_STEPS=${EQ_STEPS:-20000}
export RUN_STEPS=${RUN_STEPS:-80000}
export DUMP_STRIDE=${DUMP_STRIDE:-20000}
./runs/00_check_environment.sh
./runs/02_smoke_test_plumed_xrd_20core.sh
./runs/03_make_liquid_from_beta_20core.sh
./runs/04_run_unbiased_control_20core.sh
./runs/05_run_wtmetad_xrd111_20core.sh
./runs/07_postprocess.sh
