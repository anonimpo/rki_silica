#!/usr/bin/env bash
set -euo pipefail

# Run all pH cases sequentially. Use this for local testing; for HPC prefer slurm/job_array.sbatch.
# Example quick smoke-test:
#   LMP=lmp DRYRUN=1 ./scripts/03_run_all_ph.sh 192
# Example full baseline:
#   LMP=lmp N_EQ=200000 N_HEAT=400000 N_HOLD=400000 N_PROD=2000000 ./scripts/03_run_all_ph.sh 1536

N=${1:-1536}
for PH in pH6p0 pH6p5 pH7p0 pH7p5 pH8p0; do
  echo "=== Running $PH N=$N ==="
  ./scripts/02_run_one_pipeline.sh "$PH" "$N"
done
