#!/usr/bin/env bash
set -euo pipefail

# Capture the real OGFinder runtime after the NEO menu has marked a FITS
# sequence. This intentionally refuses to create a substitute image.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DS9="$ROOT/bin/ds9"
DETECTOR="${CODES_ROOT:-/home/kjhan/BACKUP/CODES}/../ds9/OGFinder/ds9/library/ds9_neo_detect.py"
OUT="${1:-$ROOT/../3.5ST/appendixC_assets/ogfinder_neo_runtime.png}"
REGIONS="${NEO_REGIONS:-/home/kjhan/BACKUP/CODES/output/dad_mpcid_14941/detection4/neo_candidates.reg}"
XPA="$ROOT/bin/xpaset"

if [[ -z "${DISPLAY:-}" ]]; then
  echo "DISPLAY is not set. Start a real X display and rerun this script." >&2
  exit 2
fi
if [[ ! -x "$DS9" ]]; then
  echo "OGFinder binary not found: $DS9" >&2
  exit 2
fi
if [[ ! -f "$REGIONS" ]]; then
  echo "NEO region file not found: $REGIONS" >&2
  exit 2
fi
if ! "$ROOT/bin/xpaaccess" ds9 >/dev/null 2>&1; then
  echo "No running OGFinder DS9 XPA endpoint. Start the custom binary first." >&2
  exit 2
fi

"$XPA" -p ds9 "regions file $REGIONS" >/dev/null
"$XPA" -p ds9 "regions select all" >/dev/null
"$XPA" -p ds9 "regions color green" >/dev/null

if command -v import >/dev/null 2>&1; then
  import -window root "$OUT"
elif command -v xwd >/dev/null 2>&1 && command -v convert >/dev/null 2>&1; then
  xwd -root | convert xwd:- "$OUT"
else
  echo "No X screenshot utility found (import or xwd+convert)." >&2
  exit 2
fi

echo "Wrote genuine OGFinder screen capture: $OUT"
