#!/bin/bash
# Download DES DR2 r, i, z band coadd tiles for Tanoglidis+2021 reproduction.
# Tile DES0125-0124 covers the densest LSBG patch (RA~21.2-21.8, Dec~-1.55--1.2).
#
# Usage: cd OGFinder && bash scripts/fetch_des_riz.sh
#
# Total download: ~500 MB (3 x ~165 MB fpacked FITS)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATADIR="$SCRIPT_DIR/../data"

BASE="https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/DES0125-0124"

mkdir -p "$DATADIR"

for band in r i z; do
    OUTFILE="$DATADIR/des_tanoglidis2021_${band}.fits.fz"
    URL="${BASE}/DES0125-0124_r4920p01_${band}.fits.fz"

    if [ -f "$OUTFILE" ]; then
        SIZE=$(du -h "$OUTFILE" | cut -f1)
        echo "Already exists: $OUTFILE ($SIZE)" >&2
        continue
    fi

    echo "Downloading DES DR2 ${band}-band tile..." >&2
    echo "  URL: $URL" >&2
    wget --no-verbose -O "$OUTFILE" "$URL" || {
        echo "ERROR: Failed to download ${band}-band tile" >&2
        rm -f "$OUTFILE"
        exit 1
    }
    SIZE=$(du -h "$OUTFILE" | cut -f1)
    echo "  Saved: $OUTFILE ($SIZE)" >&2
done

echo "" >&2
echo "All r/i/z tiles downloaded." >&2
echo "Next: python3 scripts/make_riz_coadd.py" >&2
