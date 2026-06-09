#!/usr/bin/env bash
set -euo pipefail
LMP=${LMP:-lmp}
NCORE=${NCORE:-20}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "== LAMMPS binary =="
command -v "$LMP" || { echo "LAMMPS binary '$LMP' not found. Set LMP=/path/to/lmp"; exit 1; }
"$LMP" -help | head -40 || true

echo "== Check PLUMED in LAMMPS styles =="
if "$LMP" -help 2>/dev/null | grep -qi "plumed"; then
  echo "LAMMPS help mentions PLUMED. Good sign."
else
  echo "WARNING: LAMMPS help did not show PLUMED. If in.00 fails with 'Unknown fix style plumed', rebuild LAMMPS with PLUMED package."
fi

echo "== PLUMED binary =="
if command -v plumed >/dev/null 2>&1; then
  plumed --version || true
else
  echo "WARNING: 'plumed' CLI not found in PATH. LAMMPS may still be linked, but postprocessing with plumed driver may fail."
fi

echo "== Project files =="
[[ -f potentials/ffield.reax.SiOH ]] && echo "ReaxFF potential found: potentials/ffield.reax.SiOH"
[[ -f structures/beta_cristobalite_a7.15_N1536.data ]] && echo "beta target found"
[[ -f custom_cv/plumed_xrd111_si_matrix_WTMETAD_v210.dat ]] && echo "PLUMED XRD-like CV found"
