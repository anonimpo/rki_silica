#!/usr/bin/env bash
set -euo pipefail
python post_processing/plot_xrd_before_after.py
python post_processing/make_figure7_like_snapshots.py
