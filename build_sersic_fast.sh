#!/bin/bash
# Build libsersic_fast.so — fast 2D Sérsic model evaluation with analytical Jacobian.
#
# Usage:
#   ./build_sersic_fast.sh          # build
#   ./build_sersic_fast.sh clean    # remove artifact
#
# Requires: gcc (or CC env var), libm. OpenMP optional but recommended.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC="$SCRIPT_DIR/sersic_fit/sersic_fast.c"
OUTDIR="$SCRIPT_DIR/lib"
OUTPUT="$OUTDIR/libsersic_fast.so"

if [ "$1" = "clean" ]; then
    rm -f "$OUTPUT"
    echo "Cleaned $OUTPUT"
    exit 0
fi

if [ ! -f "$SRC" ]; then
    echo "Error: source not found: $SRC" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

CC="${CC:-gcc}"

# Detect OpenMP support
OMP_FLAGS=""
if $CC -fopenmp -E -x c /dev/null >/dev/null 2>&1; then
    OMP_FLAGS="-fopenmp"
    echo "OpenMP: enabled"
else
    echo "OpenMP: not available, building single-threaded"
fi

echo "Compiling $SRC -> $OUTPUT"
$CC -O3 -march=native -fno-math-errno $OMP_FLAGS -shared -fPIC \
    -o "$OUTPUT" "$SRC" -lm

echo "Built: $OUTPUT"
ls -lh "$OUTPUT"
