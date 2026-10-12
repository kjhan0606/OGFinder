#!/bin/bash
# Reproduce the ICL stack debug up to the mass gate.
# Steps 01-04 always run. Step 05 writes a stack only after M_200c and
# global-sky cutouts are supplied; otherwise it records a stop.
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
echo "debug steps written under $ROOT/logs"
