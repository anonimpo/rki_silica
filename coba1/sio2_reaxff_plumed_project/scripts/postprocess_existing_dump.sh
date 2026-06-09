#!/usr/bin/env bash
set -euo pipefail
if [ $# -lt 1 ]; then
  echo "Usage: $0 path/to/dump.lammpstrj [output.dat]" >&2
  exit 1
fi
DUMP=$1
OUT=${2:-outputs/xrd_cv_from_dump.dat}
python3 scripts/xrd_debye_cv.py "$DUMP" > "$OUT"
echo "Wrote $OUT"
