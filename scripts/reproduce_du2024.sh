#!/bin/bash
# Reproduce Du et al. (2024, RAA 24 055015) LSBG detection
# using OGFinder LSBG pipeline on DESI Legacy Survey DR9 data.
#
# Usage: cd OGFinder && bash scripts/reproduce_du2024.sh
#
# Reference: "Low Surface Brightness Galaxies from BASS+MzLS
#             with Machine Learning"
# Survey: DESI Legacy Survey DR9 North (BASS+MzLS)
# Brick: 2828p530 (RA=282.8, Dec=+53.0)

set -euo pipefail

# ---- Step 0: Configuration ----

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

DATADIR="../data"
OUTDIR="results/du2024"
LSBG_PY="ds9/library/ds9_lsbg.py"

# Primary data (r-band coadd)
FITS="$DATADIR/desi_dr9_2828p530_r.fits"

# -- DESI Legacy Survey DR9 parameters --
PIXEL_SCALE=0.262     # arcsec/pixel (ground-based, BASS+MzLS)
MAG_ZP=22.5           # nanomaggy zeropoint
DETECT_THRESH=1.0     # low threshold for LSB detection
DETECT_MINAREA=100    # larger area (PSF ~1.4" FWHM -> ~5 px)
FILTER_KERNEL=gauss21x21  # large kernel matched to LSBG size (sigma=5px=1.3")
DEBLEND_NTHRESH=32    # deblending sub-thresholds
DEBLEND_MINCONT=0.01  # slightly higher to avoid over-deblending

# -- Du+2024 selection criteria --
MU_EFF_MIN=24.0       # standard LSBG cutoff (slightly wider than paper's 24.2)
MU_EFF_MAX=29.0       # faint end (wider than paper's 28.8 for completeness)
R_EFF_MIN=2.5         # arcsec (Du+2024 criterion)
R_EFF_MAX=60.0        # arcsec (upper limit)
ELLIPTICITY_MAX=0.7   # Du+2024 criterion

# -- Masking parameters (mask bright sources, keep LSB) --
MASK_THRESH=1.5       # detection threshold for masking
MASK_MAG_THRESH=20.0  # only mask bright sources (relaxed from 22)
EXPAND_FACTOR=1.0     # dilation factor (relaxed from 1.5)
MAX_DILATE=10         # max dilation radius (relaxed from 30)
BRIGHT_MAG=16.0       # bright star masking limit (relaxed from 18)

# -- Background parameters --
BKG_METHOD=sep_large  # SEP mesh-based background
BKG_MESH=256          # large mesh for gradients
BKG_NITER=3           # iterative refinement
BKG_TOL=0.01          # convergence tolerance

# -- Multi-scale detection --
MULTISCALE_FACTORS="1,2,4"

# Parallel workers (0 = auto)
NWORKERS=0

# ---- Validate data ----

if [ ! -f "$FITS" ]; then
    echo "ERROR: DESI image not found: $FITS" >&2
    echo "Run: python3 scripts/fetch_du2024_data.py" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
echo "=== Du et al. (2024) LSBG Detection Reproduction ===" >&2
echo "Data: $FITS" >&2
echo "Survey: DESI Legacy Survey DR9 (BASS+MzLS)" >&2
echo "Brick: 2828p530 (RA=282.8, Dec=+53.0)" >&2
echo "Output: $OUTDIR/" >&2
echo "" >&2

# ---- Option A: Full pipeline in one command ----

echo ">>> Running full LSBG pipeline..." >&2
echo "  pixel-scale=${PIXEL_SCALE}\"/pix, ZP=${MAG_ZP}" >&2
echo "  detect-thresh=${DETECT_THRESH}sigma, minarea=${DETECT_MINAREA}pix" >&2
echo "  mu_eff: ${MU_EFF_MIN}-${MU_EFF_MAX} mag/arcsec^2" >&2
echo "  r_eff: ${R_EFF_MIN}-${R_EFF_MAX} arcsec" >&2
echo "  ellipticity < ${ELLIPTICITY_MAX}" >&2
echo "" >&2

python3 "$LSBG_PY" "$FITS" --mode run \
    --pixel-scale $PIXEL_SCALE \
    --mag-zeropoint $MAG_ZP \
    --mask-detect-thresh $MASK_THRESH \
    --mask-mag-threshold $MASK_MAG_THRESH \
    --mask-expand-factor $EXPAND_FACTOR \
    --max-dilate-radius $MAX_DILATE \
    --bright-star-mag-limit $BRIGHT_MAG \
    --lsb-protect --lsb-mu-threshold 24.0 \
    --bkg-method $BKG_METHOD \
    --bkg-mesh-size $BKG_MESH \
    --bkg-n-iterations $BKG_NITER \
    --bkg-convergence-tol $BKG_TOL \
    --detect-thresh $DETECT_THRESH \
    --detect-minarea $DETECT_MINAREA \
    --detect-filter-kernel $FILTER_KERNEL \
    --deblend-nthresh $DEBLEND_NTHRESH \
    --deblend-mincont $DEBLEND_MINCONT \
    --multiscale --multiscale-factors "$MULTISCALE_FACTORS" \
    --mu-eff-min $MU_EFF_MIN --mu-eff-max $MU_EFF_MAX \
    --r-eff-min $R_EFF_MIN --r-eff-max $R_EFF_MAX \
    --ellipticity-max $ELLIPTICITY_MAX \
    --sersic-fit \
    --n-workers $NWORKERS \
    --catalog-output "$OUTDIR/lsbg_detected.tsv" \
    --mask-output "$OUTDIR/lsbg_mask.fits" \
    --cleaned-output "$OUTDIR/lsbg_cleaned.fits" \
    --segmap-output "$OUTDIR/lsbg_segmap.fits"

echo "" >&2

# ---- Summary ----

echo "=== Pipeline Complete ===" >&2

if [ -f "$OUTDIR/lsbg_detected.tsv" ]; then
    N_DET=$(tail -n +2 "$OUTDIR/lsbg_detected.tsv" | wc -l)
    echo "  Detected LSBGs: $N_DET" >&2
    echo "" >&2

    # Quick statistics
    echo "  Grade distribution:" >&2
    for grade in A B C D; do
        N_GRADE=$(tail -n +2 "$OUTDIR/lsbg_detected.tsv" | \
                  awk -F'\t' '{print $27}' | grep -c "^${grade}$" || true)
        echo "    Grade $grade: $N_GRADE" >&2
    done
    echo "" >&2
fi

echo "Output files:" >&2
ls -lh "$OUTDIR/"*.fits "$OUTDIR/"*.tsv 2>/dev/null | while read line; do
    echo "  $line" >&2
done
echo "" >&2
echo "Next: python3 scripts/compare_du2024.py" >&2
