#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p outputs

# 1) XRD before: use the liquid made from beta if available; otherwise use supplied amorphous pH7 sanity structure.
BEFORE="structures/liquid_from_beta_2300K_N1536.data"
if [[ ! -f "$BEFORE" ]]; then BEFORE="structures/amorphous_pH7_N1536.data"; fi
python3 postprocessing/compute_xrd_debye.py "$BEFORE" --out outputs/xrd_before_unbiased_or_amorphous_SiOnly.dat --species Si --rcut 12 --npts 1201

# 2) XRD after WTMETAD: prefer final WTMETAD data; fallback to ideal beta target for preview.
AFTER="structures/wtmetad_xrd111_after_2300K_N1536.data"
if [[ ! -f "$AFTER" ]]; then AFTER="structures/beta_cristobalite_a7.15_N1536.data"; fi
python3 postprocessing/compute_xrd_debye.py "$AFTER" --out outputs/xrd_after_wtmetad_or_beta_SiOnly.dat --species Si --rcut 12 --npts 1201
python3 postprocessing/plot_xrd_compare.py --before outputs/xrd_before_unbiased_or_amorphous_SiOnly.dat --after outputs/xrd_after_wtmetad_or_beta_SiOnly.dat --out outputs/xrd_before_after_wtmetad_check.png --title "WTMETAD crystallization check: before vs after"

# 3) Fig. 7-like six snapshots if trajectory exists.
if [[ -f outputs/dump_03_wtmetad_xrd111_2300K.lammpstrj ]]; then
  python3 postprocessing/make_fig7_like.py outputs/dump_03_wtmetad_xrd111_2300K.lammpstrj --out outputs/fig7_like_wtmetad_xrd111.png --dt-fs 0.25 --dump-stride 160000
else
  echo "No WTMETAD dump found yet; run in.03 first."
fi
