#!/bin/bash
# Build ds9_sextract binary from sep_src/ sources and ds9/library/ds9_sextract.c
#
# Prerequisites:
#   - cfitsio (brew install cfitsio on macOS, yum/dnf install cfitsio-devel on Linux)
#   - C compiler (gcc or clang)
#
# Usage:
#   ./build_sextract.sh          # auto-detect platform
#   ./build_sextract.sh clean    # remove build artifacts

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SEP_DIR="$SCRIPT_DIR/sep_src"
SRC="$SCRIPT_DIR/ds9/library/ds9_sextract.c"
BINDIR="$SCRIPT_DIR/bin"
OUTPUT="$BINDIR/ds9_sextract"

if [ "$1" = "clean" ]; then
    rm -f "$SEP_DIR"/*.o "$OUTPUT"
    echo "Cleaned."
    exit 0
fi

# Check prerequisites
if [ ! -f "$SRC" ]; then
    echo "ERROR: $SRC not found" >&2
    exit 1
fi
if [ ! -d "$SEP_DIR" ]; then
    echo "ERROR: $SEP_DIR not found" >&2
    exit 1
fi

# Detect compiler
CC="${CC:-cc}"

# Detect cfitsio
CFITSIO_CFLAGS=""
CFITSIO_LIBS="-lcfitsio"

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists cfitsio 2>/dev/null; then
    CFITSIO_CFLAGS="$(pkg-config --cflags cfitsio)"
    CFITSIO_LIBS="$(pkg-config --libs cfitsio)"
else
    # Search common paths: conda, Homebrew, system
    FOUND=0
    for prefix in "$CONDA_PREFIX" "$HOME/miniconda3" "$HOME/anaconda3" \
                  /opt/homebrew /usr/local /usr; do
        if [ -n "$prefix" ] && [ -f "$prefix/include/fitsio.h" ]; then
            CFITSIO_CFLAGS="-I$prefix/include"
            CFITSIO_LIBS="-L$prefix/lib -lcfitsio"
            FOUND=1
            break
        fi
    done
    if [ "$FOUND" = "0" ]; then
        echo "ERROR: cfitsio not found. Install via:" >&2
        echo "  conda install cfitsio    (or)" >&2
        echo "  brew install cfitsio     (macOS)" >&2
        echo "  dnf install cfitsio-devel (RHEL/Fedora)" >&2
        exit 1
    fi
fi
echo "cfitsio: $CFITSIO_CFLAGS $CFITSIO_LIBS"

# Compile sep source files
SEP_OBJS=""
echo "Compiling SEP library..."
for src in "$SEP_DIR"/analyse.c "$SEP_DIR"/aperture.c "$SEP_DIR"/background.c \
           "$SEP_DIR"/convolve.c "$SEP_DIR"/deblend.c "$SEP_DIR"/extract.c \
           "$SEP_DIR"/lutz.c "$SEP_DIR"/util.c; do
    obj="${src%.c}.o"
    echo "  $(basename $src)"
    $CC -O2 -I"$SEP_DIR" -c "$src" -o "$obj"
    SEP_OBJS="$SEP_OBJS $obj"
done

# Compile and link ds9_sextract
mkdir -p "$BINDIR"
echo "Compiling ds9_sextract..."
$CC -O2 -I"$SEP_DIR" $CFITSIO_CFLAGS -o "$OUTPUT" "$SRC" $SEP_OBJS $CFITSIO_LIBS -lm
echo "Built: $OUTPUT"
