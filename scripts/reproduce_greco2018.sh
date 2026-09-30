#!/bin/bash
# Reproduce Greco et al. (2018, ApJ 857 104) LSBG detection
# using OGFinder LSBG pipeline on HSC-SSP data.
#
# Usage: cd OGFinder && bash scripts/reproduce_greco2018.sh
#
# Reference: "Illuminating Low Surface Brightness Galaxies with the
#             Hyper Suprime-Cam Survey"
# Survey: HSC-SSP Wide (S16A → PDR3)
# Method: SExtractor-based detection — fair comparison with OGFinder/SEP

set -euo pipefail

# ---- Step 0: Configuration ----

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

DATADIR="../data"
OUTDIR="results/greco2018"
LSBG_PY="ds9/library/ds9_lsbg.py"

# Primary data (i-band coadd from HSC-SSP)
FITS="$DATADIR/hsc_greco2018_i.fits"

# -- HSC-SSP parameters --
PIXEL_SCALE=0.168     # arcsec/pixel (HSC)
MAG_ZP=27.0           # HSC-SSP coadd zeropoint (auto-detect if header has it)

# -- Greco+2018 SExtractor parameters (reproduced) --
# Original: DETECT_THRESH=0.7, DETECT_MINAREA=100
# Kernel: gauss_6.0_31x31.conv (FWHM=6px → σ≈2.55px)
# Our closest: gauss9x9 (σ=2.5px, FWHM≈5.9px)
DETECT_THRESH=0.7      # very low threshold for LSB detection
DETECT_MINAREA=100     # minimum area (larger than typical PSF)
FILTER_KERNEL=gauss9x9 # matches Greco FWHM=6px kernel (σ=2.5px)
DEBLEND_NTHRESH=32     # deblending sub-thresholds
DEBLEND_MINCONT=0.001  # Greco used 0.001

# -- Greco+2018 selection criteria --
MU_EFF_MIN=24.3        # g-band mu_eff lower limit (mag/arcsec^2)
MU_EFF_MAX=28.8        # g-band mu_eff upper limit
R_EFF_MIN=2.5          # arcsec (Greco criterion)
R_EFF_MAX=14.0         # arcsec (Greco upper limit)
ELLIPTICITY_MAX=0.7    # Greco criterion

# -- Masking parameters (conservative, protect diffuse sources) --
MASK_THRESH=2.0        # detection threshold for masking
MASK_MAG_THRESH=19.0   # only mask bright sources
EXPAND_FACTOR=1.0      # minimal dilation (don't mask LSBGs)
MAX_DILATE=10          # max dilation radius
BRIGHT_MAG=15.0        # bright star masking limit

# -- Background parameters --
# Greco used BACK_SIZE=128, BACK_FILTERSIZE=5
BKG_METHOD=sep_large   # SEP mesh-based background
BKG_MESH=128           # matches Greco BACK_SIZE
BKG_NITER=3            # iterative refinement
BKG_TOL=0.01           # convergence tolerance

# -- Multi-scale detection (Greco+2018 inspired) --
MULTISCALE_FACTORS="1,2,4"

# Parallel workers (0 = auto)
NWORKERS=0

# ---- Validate data ----

if [ ! -f "$FITS" ]; then
    echo "ERROR: HSC-SSP image not found: $FITS" >&2
    echo "Run: python3 scripts/fetch_greco2018_data.py" >&2
    echo "" >&2
    echo "HSC-SSP requires free registration:" >&2
    echo "  https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
echo "=== Greco et al. (2018) LSBG Detection Reproduction ===" >&2
echo "Data: $FITS" >&2
echo "Survey: HSC-SSP Wide (PDR3)" >&2
echo "Pixel scale: ${PIXEL_SCALE}\"/pixel" >&2
echo "Output: $OUTDIR/" >&2
echo "" >&2

# ---- Run full LSBG pipeline ----

echo ">>> Running full LSBG pipeline..." >&2
echo "  pixel-scale=${PIXEL_SCALE}\"/pix, ZP=${MAG_ZP}" >&2
echo "  detect-thresh=${DETECT_THRESH}sigma, minarea=${DETECT_MINAREA}pix" >&2
echo "  filter-kernel=${FILTER_KERNEL}" >&2
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
echo "Next: python3 scripts/compare_greco2018.py" >&2
