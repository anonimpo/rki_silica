#!/usr/bin/env bash
set -euo pipefail
LMP=${LMP:-lmp}

echo "== LAMMPS executable =="
if command -v "$LMP" >/dev/null 2>&1; then
  "$LMP" -h | head -80 || true
else
  echo "LAMMPS executable '$LMP' not found. Set LMP=/path/to/lmp."
fi

echo "\n== PLUMED executable =="
if command -v plumed >/dev/null 2>&1; then
  plumed info --version || true
else
  echo "plumed command not found. Needed only for Stage 05."
fi

echo "\n== Python packages =="
python3 - <<'PY'
for pkg in ['numpy','pandas','matplotlib']:
    try:
        __import__(pkg)
        print(f'{pkg}: OK')
    except Exception as exc:
        print(f'{pkg}: missing ({exc})')
PY
