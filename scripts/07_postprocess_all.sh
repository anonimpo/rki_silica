#!/usr/bin/env bash
set -euo pipefail

# Post-process the calcination stage for all pH cases.
# Usage:
#   ./scripts/07_postprocess_all.sh 1536 03_calcination_1173K calcined_1173K.data
# or for WTMetaD final:
#   ./scripts/07_postprocess_all.sh 1536 05_wtmetad_plumed wtmetad_final.data

N=${1:-1536}
STAGE=${2:-03_calcination_1173K}
DATA=${3:-calcined_1173K.data}
OUTROOT="analysis/${STAGE}_N${N}"
mkdir -p "$OUTROOT/xrd" "$OUTROOT/metrics" "outputs/${STAGE}_N${N}"

for PH in pH6p0 pH6p5 pH7p0 pH7p5 pH8p0; do
  IN="runs/${PH}_N${N}/${STAGE}/${DATA}"
  if [[ ! -f "$IN" ]]; then
    echo "Skipping $PH: missing $IN" >&2
    continue
  fi
  python3 scripts/04_compute_debye_xrd.py --input "$IN" --format data \
    --out "$OUTROOT/xrd/xrd_${PH}_N${N}.csv" \
    --plot "outputs/${STAGE}_N${N}/xrd_${PH}_N${N}.png"
  python3 scripts/05_structural_metrics.py --input "$IN" --format data --ph "$PH" \
    --out "$OUTROOT/metrics/metrics_${PH}_N${N}.csv"
done

python3 scripts/08_plot_summary.py --analysis-dir "$OUTROOT" --outdir "outputs/${STAGE}_N${N}"
