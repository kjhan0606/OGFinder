#!/bin/bash
# Reproduce Tanoglidis et al. (2021, ApJS 252 18) LSBG detection
# using OGFinder LSBG pipeline on DES DR2 g-band coadd.
#
# Usage: cd OGFinder && bash scripts/reproduce_tanoglidis2021.sh
#
# Reference: "Shadows in the Dark: Low-Surface-Brightness Galaxies
#             Discovered in the Dark Energy Survey"
# Original: DES Y3 SourceExtractor → SVM → 23,790 LSBGs
# Our test: SEP detection on DES DR2 g-band coadd → same selection cuts

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

DATADIR="../data"
OUTDIR="results/tanoglidis2021"
LSBG_PY="ds9/library/ds9_lsbg.py"

# Primary data: DES DR2 g-band coadd (cutout from tile DES0125-0124)
# Covers all 19 reference LSBGs (RA~21.2-21.8, Dec~-1.55--1.2)
FITS="$DATADIR/des_tanoglidis2021_g_cutout.fits"

# -- DES parameters --
PIXEL_SCALE=0.263     # arcsec/pixel (DES standard)
MAG_ZP=30.0           # DES DR2 MAGZERO

# -- Detection parameters --
# DES Y3 Gold standard SExtractor uses ~1.5sigma
DETECT_THRESH=1.5      # sigma (match DES standard)
DETECT_MINAREA=50      # min connected pixels
FILTER_KERNEL=gauss7x7 # sigma=2px, FWHM~4.7px (~1.2" at 0.263"/px), match DES SExtractor
DEBLEND_NTHRESH=32
DEBLEND_MINCONT=0.005

# -- Tanoglidis+2021 selection criteria --
MU_EFF_MIN=24.2        # g-band mu_eff lower limit (mag/arcsec^2)
MU_EFF_MAX=28.5        # Tanoglidis+2021 effective limit
R_EFF_MIN=2.5          # arcsec
R_EFF_MAX=20.0         # arcsec (Tanoglidis upper limit)
ELLIPTICITY_MAX=0.7

# -- Masking (conservative) --
MASK_THRESH=2.0
MASK_MAG_THRESH=19.0
EXPAND_FACTOR=1.0
MAX_DILATE=10
BRIGHT_MAG=15.0

# -- Background --
BKG_METHOD=sep_large
BKG_MESH=256
BKG_NITER=3
BKG_TOL=0.01

# -- Multi-scale --
MULTISCALE_FACTORS="1,2,4"

# Parallel workers
NWORKERS=0

# No FITS extension needed (cutout is standard PrimaryHDU)

# ---- Validate ----

if [ ! -f "$FITS" ]; then
    echo "ERROR: DES g-band cutout not found: $FITS" >&2
    echo "Run: python3 scripts/fetch_tanoglidis2021_data.py" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
echo "=== Tanoglidis+2021 LSBG Detection Reproduction ===" >&2
echo "Data: $FITS" >&2
echo "Survey: DES DR2 g-band coadd (tile DES0125-0124, EXPTIME=1170s)" >&2
echo "Pixel scale: ${PIXEL_SCALE}\"/pixel, ZP=${MAG_ZP}" >&2
echo "Output: $OUTDIR/" >&2
echo "" >&2

# ---- Run LSBG pipeline ----

echo ">>> Running LSBG pipeline..." >&2
echo "  detect-thresh=${DETECT_THRESH}sigma, minarea=${DETECT_MINAREA}pix" >&2
echo "  filter-kernel=${FILTER_KERNEL}" >&2
echo "  mu_eff: ${MU_EFF_MIN}-${MU_EFF_MAX} mag/arcsec^2" >&2
echo "  r_eff: ${R_EFF_MIN}-${R_EFF_MAX} arcsec" >&2
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

    echo "  Grade distribution:" >&2
    for grade in A B C D; do
        N_GRADE=$(tail -n +2 "$OUTDIR/lsbg_detected.tsv" | \
                  grep -c "$grade" || true)
        echo "    Grade $grade: $N_GRADE" >&2
    done
    echo "" >&2
fi

echo "Output files:" >&2
ls -lh "$OUTDIR/"*.fits "$OUTDIR/"*.tsv 2>/dev/null | while read line; do
    echo "  $line" >&2
done
echo "" >&2
echo "Next: python3 scripts/compare_tanoglidis2021.py" >&2
