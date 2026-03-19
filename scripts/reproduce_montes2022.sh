#!/bin/bash
# Reproduce Montes & Trujillo (2022, ApJ 940 L51) ICL analysis of SMACS J0723
# using OGFinder CLI pipeline.
#
# Usage: cd OGFinder && bash scripts/reproduce_montes2022.sh
#
# Reference: "A New Era of Intracluster Light Studies with JWST"
# Cluster: SMACS J0723.3-7327 (z=0.39), JWST ERO Program 2736

set -euo pipefail

# ---- Step 0: Configuration ----

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

DATADIR="../data"
OUTDIR="results/montes2022"
ICL_PY="ds9/library/ds9_icl.py"

# Primary data (F444W cutout, 115 MB)
FITS="$DATADIR/jwst_smacs0723_nircam_f444w_cutout.fits"
# Secondary data (F090W, 1.7 GB, optional for color)
F090W="$DATADIR/jwst_smacs0723_f090w.fits"

# -- Paper parameters (Montes & Trujillo 2022) --
PIXEL_SCALE=0.063    # arcsec/pixel (NIRCam LW, Module B)
MAG_ZP=auto          # auto-compute from FITS header (PIXAR_SR)
EXPAND_FACTOR=1.5    # moderate dilation (preserve ICL signal)
MAX_DILATE=20        # capped dilation to prevent over-masking
DETECT_THRESH=5.0    # detect only bright compact sources
BRIGHT_MAG=18.0      # bright star masking limit
NSTEPS=51            # paper: 51 logarithmic radial bins
RMIN=3.0             # ~0.19" (inner cutoff)
RMAX=1143            # ~72" ≈ 380 kpc at z=0.39
BKG_ORDER=2          # 2D polynomial plane (paper: polynomial fit)

# BCG center (brightest galaxy, found via smoothed peak detection)
CENTER="1004,1067"

# Parallel workers (0 = auto)
NWORKERS=0

# ---- Validate data ----

if [ ! -f "$FITS" ]; then
    echo "ERROR: F444W cutout not found: $FITS" >&2
    echo "Download from MAST (JWST ERO Program 2736) or check data/ directory" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
echo "=== Montes & Trujillo 2022 ICL Reproduction ===" >&2
echo "Data: $FITS" >&2
echo "Output: $OUTDIR/" >&2
echo "" >&2

# ---- Step 1: Source Masking ----

echo ">>> Step 1: Source masking (SEP-based, expand=${EXPAND_FACTOR}x)..." >&2

python3 "$ICL_PY" "$FITS" --mode mask \
    --expand-factor $EXPAND_FACTOR \
    --max-dilate-radius $MAX_DILATE \
    --detect-thresh $DETECT_THRESH \
    --bright-star-mag-limit $BRIGHT_MAG \
    --interp-method linear \
    --mask-output "$OUTDIR/icl_mask.fits" \
    --masked-output "$OUTDIR/icl_masked.fits"

echo "  Mask: $OUTDIR/icl_mask.fits" >&2
echo "  Masked image: $OUTDIR/icl_masked.fits" >&2
echo "" >&2

# ---- Step 2: Background Model (2D Polynomial) ----

echo ">>> Step 2: Background model (polynomial order=${BKG_ORDER}, iterative)..." >&2

python3 "$ICL_PY" "$FITS" --mode background \
    --mask "$OUTDIR/icl_mask.fits" \
    --bkg-method polynomial \
    --bkg-order $BKG_ORDER \
    --bkg-sigma-clip 3.0 \
    --iterative \
    --bkg-n-iterations 3 \
    --bkg-convergence-tol 0.01 \
    --bkg-output "$OUTDIR/icl_background.fits" \
    --bgsub-output "$OUTDIR/icl_bgsub.fits"

echo "  Background: $OUTDIR/icl_background.fits" >&2
echo "  Bg-subtracted: $OUTDIR/icl_bgsub.fits" >&2
echo "" >&2

# ---- Step 3: Circular SB Profile (main result) ----

BGSUB="$OUTDIR/icl_bgsub.fits"

echo ">>> Step 3: Circular SB profile (${NSTEPS} log bins, r=${RMIN}-${RMAX} pix)..." >&2

python3 "$ICL_PY" "$BGSUB" --mode profile \
    --center "$CENTER" \
    --rmin $RMIN \
    --rmax $RMAX \
    --nsteps $NSTEPS \
    --spacing log \
    --mag-zeropoint $MAG_ZP \
    --pixel-scale $PIXEL_SCALE \
    --profile-output "$OUTDIR/sb_profile.tsv" \
    > "$OUTDIR/sb_profile_stdout.tsv"

echo "  Profile: $OUTDIR/sb_profile.tsv" >&2
echo "" >&2

# ---- Step 4: Sector Profiles (East / West) ----

echo ">>> Step 4: Sector profiles (East PA=0, West PA=180, width=20 deg)..." >&2

python3 "$ICL_PY" "$BGSUB" --mode profile \
    --center "$CENTER" \
    --rmin $RMIN --rmax $RMAX --nsteps $NSTEPS --spacing log \
    --sector-width 20.0 --sector-pa 0.0 \
    --mag-zeropoint $MAG_ZP --pixel-scale $PIXEL_SCALE \
    --profile-output "$OUTDIR/sb_profile_east.tsv" \
    > /dev/null

python3 "$ICL_PY" "$BGSUB" --mode profile \
    --center "$CENTER" \
    --rmin $RMIN --rmax $RMAX --nsteps $NSTEPS --spacing log \
    --sector-width 20.0 --sector-pa 180.0 \
    --mag-zeropoint $MAG_ZP --pixel-scale $PIXEL_SCALE \
    --profile-output "$OUTDIR/sb_profile_west.tsv" \
    > /dev/null

echo "  East: $OUTDIR/sb_profile_east.tsv" >&2
echo "  West: $OUTDIR/sb_profile_west.tsv" >&2
echo "" >&2

# ---- Step 5: ICL Fraction (Multi-Threshold Sweep) ----

echo ">>> Step 5: ICL fraction (multi-threshold 25.0-29.0, bootstrap N=500)..." >&2

python3 "$ICL_PY" "$BGSUB" --mode measure-multi \
    --profile-file "$OUTDIR/sb_profile.tsv" \
    --mu-min 25.0 --mu-max 29.0 --mu-step 0.5 \
    --n-bootstrap 500 \
    --pixel-scale $PIXEL_SCALE \
    > "$OUTDIR/icl_fraction_multi.tsv"

echo "  Result: $OUTDIR/icl_fraction_multi.tsv" >&2
echo "" >&2

# ---- Step 6: BCG+ICL Decomposition (Double Sérsic) ----

echo ">>> Step 6: BCG+ICL decomposition (double Sérsic fit)..." >&2

python3 "$ICL_PY" "$BGSUB" --mode decompose \
    --profile-file "$OUTDIR/sb_profile.tsv" \
    --pixel-scale $PIXEL_SCALE \
    --mag-zeropoint $MAG_ZP \
    > "$OUTDIR/decompose_result.tsv"

echo "  Result: $OUTDIR/decompose_result.tsv" >&2
echo "" >&2

# ---- Step 7: Multi-band Color Profile (F090W + F444W, optional) ----

if [ -f "$F090W" ]; then
    echo ">>> Step 7: Color profile (F090W - F444W)..." >&2

    python3 "$ICL_PY" "$FITS" --mode color \
        --mask "$OUTDIR/icl_mask.fits" \
        --center "$CENTER" \
        --bands "F090W:$F090W,F444W:$FITS" \
        --rmin $RMIN --rmax $RMAX --nsteps $NSTEPS --spacing log \
        --mag-zeropoint $MAG_ZP --pixel-scale $PIXEL_SCALE \
        > "$OUTDIR/color_profile.tsv"

    echo "  Result: $OUTDIR/color_profile.tsv" >&2
else
    echo ">>> Step 7: SKIPPED (F090W not found: $F090W)" >&2
fi
echo "" >&2

# ---- Summary ----

echo "=== Pipeline Complete ===" >&2
echo "Output files:" >&2
ls -lh "$OUTDIR/"*.fits "$OUTDIR/"*.tsv 2>/dev/null | while read line; do
    echo "  $line" >&2
done
echo "" >&2
echo "Next: python3 scripts/compare_montes2022.py" >&2
