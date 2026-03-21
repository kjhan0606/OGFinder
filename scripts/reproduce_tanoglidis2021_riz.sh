#!/bin/bash
# Reproduce Tanoglidis et al. (2021, ApJS 252 18) LSBG detection
# using r+i+z coadd for detection + g-band forced photometry.
#
# Tanoglidis+2021 used r+i+z coadd for detection, then measured
# photometry on individual bands. This script follows the same approach.
#
# Key insight: LSBGs are NOT low surface brightness in riz (~0.6 mag brighter).
# Therefore, we use simple SEP detection (no aggressive masking/background)
# on the riz coadd, then apply LSBG selection criteria on g-band photometry.
#
# Usage: cd OGFinder && bash scripts/reproduce_tanoglidis2021_riz.sh
#
# Prerequisites:
#   bash scripts/fetch_des_riz.sh          # download r/i/z tiles
#   python3 scripts/make_riz_coadd.py      # create coadd

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

DATADIR="../data"
OUTDIR="results/tanoglidis2021_riz"
LSBG_PY="ds9/library/ds9_lsbg.py"
SEP_DETECT="scripts/sep_detect_riz.py"

# Detection image: r+i+z coadd
RIZ_COADD="$DATADIR/des_tanoglidis2021_riz_coadd.fits"
# Forced photometry image: g-band cutout (same pixel grid)
G_FITS="$DATADIR/des_tanoglidis2021_g_cutout.fits"

# -- DES parameters --
PIXEL_SCALE=0.263     # arcsec/pixel
MAG_ZP=30.0           # DES DR2 MAGZERO

# -- Detection parameters --
DETECT_THRESH=1.5
DETECT_MINAREA=50
DEBLEND_NTHRESH=32
DEBLEND_MINCONT=0.005
FILTER_FWHM=5.0       # Gaussian kernel FWHM in pixels
BKG_SIZE=256

# -- Size/shape pre-filter (applied on riz detection) --
R_EFF_MIN=2.5
R_EFF_MAX=20.0
ELLIPTICITY_MAX=0.7

# -- Tanoglidis+2021 g-band selection criteria --
MU_EFF_MIN=24.2
MU_EFF_MAX=28.5

# ---- Validate ----

if [ ! -f "$RIZ_COADD" ]; then
    echo "ERROR: r+i+z coadd not found: $RIZ_COADD" >&2
    echo "Run:" >&2
    echo "  bash scripts/fetch_des_riz.sh" >&2
    echo "  python3 scripts/make_riz_coadd.py" >&2
    exit 1
fi

if [ ! -f "$G_FITS" ]; then
    echo "ERROR: g-band cutout not found: $G_FITS" >&2
    echo "Run: python3 scripts/fetch_tanoglidis2021_data.py" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
echo "=== Tanoglidis+2021 LSBG Detection (r+i+z Coadd) ===" >&2
echo "Detection image: $RIZ_COADD (r+i+z inverse-variance coadd)" >&2
echo "Forced phot:     $G_FITS (g-band)" >&2
echo "Pixel scale: ${PIXEL_SCALE}\"/pixel, ZP=${MAG_ZP}" >&2
echo "Output: $OUTDIR/" >&2
echo "" >&2

# ---- Step 1: Simple SEP detection on r+i+z coadd ----
# No aggressive masking or iterative background subtraction.
# LSBGs appear as moderate-brightness extended sources in riz.
# Only apply size/shape filter; surface brightness cuts are for g-band.

echo ">>> Step 1: SEP detection on r+i+z coadd..." >&2
echo "  detect-thresh=${DETECT_THRESH}sigma, minarea=${DETECT_MINAREA}pix" >&2
echo "  r_eff: ${R_EFF_MIN}-${R_EFF_MAX}\", e<${ELLIPTICITY_MAX}" >&2
echo "" >&2

python3 "$SEP_DETECT" "$RIZ_COADD" \
    --pixel-scale $PIXEL_SCALE \
    --mag-zeropoint $MAG_ZP \
    --detect-thresh $DETECT_THRESH \
    --detect-minarea $DETECT_MINAREA \
    --deblend-nthresh $DEBLEND_NTHRESH \
    --deblend-mincont $DEBLEND_MINCONT \
    --filter-fwhm $FILTER_FWHM \
    --bkg-size $BKG_SIZE \
    --r-eff-min $R_EFF_MIN \
    --r-eff-max $R_EFF_MAX \
    --ellipticity-max $ELLIPTICITY_MAX \
    --output "$OUTDIR/lsbg_riz_detected.tsv"

echo "" >&2

if [ ! -f "$OUTDIR/lsbg_riz_detected.tsv" ]; then
    echo "ERROR: Detection failed, no output catalog" >&2
    exit 1
fi

N_RIZ=$(tail -n +2 "$OUTDIR/lsbg_riz_detected.tsv" | wc -l)
echo "  Detected on r+i+z coadd: $N_RIZ candidates" >&2

# ---- Step 2: g-band forced photometry ----

echo "" >&2
echo ">>> Step 2: g-band forced photometry at detected positions..." >&2

python3 "$LSBG_PY" "$G_FITS" --mode forced \
    --pixel-scale $PIXEL_SCALE \
    --mag-zeropoint $MAG_ZP \
    --catalog "$OUTDIR/lsbg_riz_detected.tsv" \
    --band-name g \
    > "$OUTDIR/lsbg_forced_g.tsv"

echo "  Forced photometry saved: $OUTDIR/lsbg_forced_g.tsv" >&2

# ---- Step 3: Merge detection catalog + forced g-band photometry ----

echo "" >&2
echo ">>> Step 3: Merging catalogs and applying g-band selection..." >&2

python3 -c "
import sys

# Load detection catalog
det_lines = open('$OUTDIR/lsbg_riz_detected.tsv').read().strip().split('\n')
det_header = det_lines[0].split('\t')
det_rows = [line.split('\t') for line in det_lines[1:]]

# Load forced photometry
forced_lines = open('$OUTDIR/lsbg_forced_g.tsv').read().strip().split('\n')
forced_header = forced_lines[0].split('\t')
forced_rows = [line.split('\t') for line in forced_lines[1:]]

# Build merged header with MU_EFF_G column
merge_header = det_header + [c for c in forced_header if c != 'NUMBER'] + ['MU_EFF_G']
print('\t'.join(merge_header))

import math

n_pass = 0
for i, det_row in enumerate(det_rows):
    if i >= len(forced_rows):
        break

    # Parse forced g-band magnitude
    forced_dict = dict(zip(forced_header, forced_rows[i]))
    mag_g = float(forced_dict.get('MAG_g', '99'))

    # Two r_eff measures:
    # - R_EFF_ARCSEC (moment-based): stable, used for size filtering
    # - FLUX_RADIUS_ARCSEC (half-light): accurate, used for mu_eff_g
    fr_idx = det_header.index('FLUX_RADIUS_ARCSEC') if 'FLUX_RADIUS_ARCSEC' in det_header else -1
    r_eff_idx = det_header.index('R_EFF_ARCSEC') if 'R_EFF_ARCSEC' in det_header else -1

    # Moment r_eff for size filter
    r_eff_moment = float(det_row[r_eff_idx]) if r_eff_idx >= 0 else 5.0
    # Flux radius for mu_eff_g calculation
    r_eff_mu = float(det_row[fr_idx]) if fr_idx >= 0 else r_eff_moment
    # Clamp flux_radius to reasonable range for mu_eff calculation
    r_eff_mu = max(r_eff_mu, r_eff_moment * 0.5)
    r_eff_mu = min(r_eff_mu, 20.0)

    # Compute g-band mu_eff using flux_radius
    if mag_g < 90 and r_eff_mu > 0:
        mu_eff_g = mag_g + 2.5 * math.log10(2 * math.pi * r_eff_mu**2)
    else:
        mu_eff_g = 99.0

    # Apply Tanoglidis+2021 g-band selection
    # Size filter uses moment-based r_eff (stable for blended sources)
    if mu_eff_g < $MU_EFF_MIN or mu_eff_g > $MU_EFF_MAX:
        continue
    if r_eff_moment < $R_EFF_MIN or r_eff_moment > $R_EFF_MAX:
        continue

    # Merge columns
    forced_vals = [forced_dict[c] for c in forced_header if c != 'NUMBER']
    merged = det_row + forced_vals + [f'{mu_eff_g:.3f}']
    print('\t'.join(merged))
    n_pass += 1

print(f'  g-band selection: {n_pass}/{len(det_rows)} passed '
      f'(mu_eff_g {$MU_EFF_MIN}-{$MU_EFF_MAX}, r_eff {$R_EFF_MIN}-{$R_EFF_MAX}\")',
      file=sys.stderr)
" > "$OUTDIR/lsbg_detected.tsv"

N_FINAL=$(tail -n +2 "$OUTDIR/lsbg_detected.tsv" | wc -l)
echo "  Final LSBG candidates: $N_FINAL" >&2

# ---- Summary ----

echo "" >&2
echo "=== Pipeline Complete ===" >&2
echo "  r+i+z detection: $N_RIZ candidates" >&2
echo "  After g-band selection: $N_FINAL LSBGs" >&2

echo "" >&2
echo "Output files:" >&2
ls -lh "$OUTDIR/"*.tsv 2>/dev/null | while read line; do
    echo "  $line" >&2
done
echo "" >&2
echo "Next: python3 scripts/compare_tanoglidis2021.py --riz" >&2
