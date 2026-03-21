#!/bin/bash
# Train LSBG SVM classifier using Tanoglidis+2021 data
# Usage: cd OGFinder && bash scripts/train_lsbg_svm.sh
#
# Prerequisites:
#   1. Run the LSBG pipeline on Tanoglidis+2021 data first:
#      bash scripts/reproduce_tanoglidis2021.sh
#   2. Ensure results/tanoglidis2021/lsbg_detected.tsv exists
#   3. Ensure ../data/tanoglidis2021_lsbg_patch.csv exists

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/.."

DETECTIONS="results/tanoglidis2021/lsbg_detected.tsv"
REFERENCE="../data/tanoglidis2021_lsbg_patch.csv"
LABELED="lsbg/data/training/labeled.tsv"
CHECKPOINT="lsbg/data/checkpoints/svm_lsbg_best.joblib"

# Check prerequisites
if [ ! -f "$DETECTIONS" ]; then
    echo "ERROR: Detection catalog not found: $DETECTIONS"
    echo "Run: bash scripts/reproduce_tanoglidis2021.sh first"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference catalog not found: $REFERENCE"
    exit 1
fi

echo "=== Step 1: Generate labeled training data ==="
python3 -m lsbg.training.make_labels \
    --detections "$DETECTIONS" \
    --reference "$REFERENCE" \
    --ra-col ra_se --dec-col dec_se \
    --match-radius 3.0 \
    --output "$LABELED"

echo ""
echo "=== Step 2: Train SVM classifier ==="
python3 -m lsbg.training.train_svm \
    --data "$LABELED" \
    --output "$CHECKPOINT"

echo ""
echo "=== Done ==="
echo "Model saved: $CHECKPOINT"
echo ""
echo "To use with LSBG pipeline:"
echo "  python3 ds9/library/ds9_lsbg.py FITS --mode run --svm-classify"
