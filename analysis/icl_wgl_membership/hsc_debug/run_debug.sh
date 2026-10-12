#!/bin/bash
# Reproduce the ICL stack debug up to the mass gate.
# Steps 01-04 always run. Step 05 writes the menu stack only after M_200c
# and global-sky cutouts are supplied. Steps 06-08 match CAMIRA members,
# probe one HSC cutout, and average BCG-masked frames on a kpc grid when
# cutouts/frames.tsv names those frames. Missing inputs record a stop.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${ASTRAFEK_DATA:-/gpfs/kjhan/Astrafek/data}
PYTHON=${ASTRAFEK_PYTHON:-/home/kjhan/miniconda3/bin/python3}
export ASTRAFEK_DATA="$ROOT"
export ASTRAFEK_PYTHON="$PYTHON"
bash "$HERE/01_fetch_catalog.sh"
"$PYTHON" "$HERE/02_select_clusters.py"
"$PYTHON" "$HERE/03_mass_gate.py"
"$PYTHON" "$HERE/04_cutout_status.py"
"$PYTHON" "$HERE/05_stack.py"
"$PYTHON" "$HERE/06_fetch_members.py"
"$PYTHON" "$HERE/07_fetch_cutouts.py"
"$PYTHON" "$HERE/08_mask_stack.py"
echo "debug steps written under $ROOT/logs"
