#!/usr/bin/env bash
set -euo pipefail
LMP=${1:-${LAMMPS_BIN:-lmp}}

echo "Using LAMMPS command: $LMP"
if ! command -v "$LMP" >/dev/null 2>&1; then
  echo "ERROR: LAMMPS executable not found. Set LAMMPS_BIN=/path/to/lmp or pass it as argument." >&2
  exit 1
fi

"$LMP" -help > logs/lammps_help.txt 2>&1 || true

echo "Checking styles/packages from lammps_help.txt"
for token in reaxff qeq/reaxff plumed; do
  if grep -qi "$token" logs/lammps_help.txt; then
    echo "OK: found $token"
  else
    echo "MISSING or not listed: $token"
  fi
done

echo "\nIf 'plumed' is missing, you can still use lammps/in.04_reaxff_production_no_plumed and Python post-processing."
echo "For metadynamics bias, you need a single LAMMPS executable with both reaxff and plumed."
