#!/bin/bash
# Product source extraction is ogfmeas/sextract.py.
# This script does not compile sep_src or ds9/library/ds9_sextract.c
# and does not write bin/ds9_sextract. sep_src is LGPL-3.0+ and stays
# in the tree; it is not linked into a product binary.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [ "${1:-}" = "clean" ]; then
    rm -f "$SCRIPT_DIR"/sep_src/*.o "$SCRIPT_DIR"/bin/ds9_sextract
    echo "Removed leftover SEP objects. No product binary is built."
    exit 0
fi

echo "The product source extractor is ogfmeas/sextract.py." >&2
echo "build_sextract.sh does not link sep_src (LGPL-3.0+) into a product binary." >&2
exit 2
