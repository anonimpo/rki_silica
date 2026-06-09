#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
python3 postprocessing/generate_beta_cristobalite.py --a 7.15 --nrep 4 --outdir structures
python3 postprocessing/compute_xrd_debye.py structures/beta_cristobalite_a7.15_N1536.data --out outputs/xrd_beta_cristobalite_target_SiOnly.dat --species Si --rcut 12 --npts 1201
