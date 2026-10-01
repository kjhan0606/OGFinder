#!/bin/bash
# Behaviour-capture test for the ::ogf::cat decoupling (see scripts/verify_cat_behavior.tcl).
#   scripts/verify_cat_behavior.sh            compare with scripts/golden/cat_behavior.golden  (exit 0 = identical)
#   scripts/verify_cat_behavior.sh --update   rewrite the golden file (only do this deliberately!)
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
DISP="${DISPLAY_OVERRIDE:-:77}"; GOLD="$HERE/golden/cat_behavior.golden"
H=$(mktemp -d /tmp/ogf_cat_home.XXXXXX); OUT="$H/cat_out.txt"
cd "$ROOT"; unset OGFINDER_PYTHON   # the golden was captured with the default interpreter name (python3)
DISPLAY=$DISP HOME=$H OGF_CAT_OUT=$OUT timeout -s KILL 300 bin/ds9 "${OGF_TEST_FITS:-/workspace/fits}/m51.fits" -geometry 1300x950 -source scripts/verify_cat_behavior.tcl > $H/stdout.txt 2>&1
rc=$?
sed -e "s#$H#<HOME>#g" "$OUT" > "$H/cat_norm.txt"
if [ "$1" = "--update" ]; then cp "$H/cat_norm.txt" "$GOLD"; echo "golden updated ($(wc -l < $GOLD) lines)"; exit 0; fi
if ! grep -q SUMMARY-DONE "$OUT"; then echo "FAIL: harness did not finish (rc=$rc); see $H"; tail -5 "$OUT"; exit 1; fi
if diff -u "$GOLD" "$H/cat_norm.txt" > "$H/diff.txt"; then echo "cat behaviour identical to golden ($(grep -c '^===' $GOLD) features, $(wc -l < $GOLD) lines)"; rm -rf "$H"; exit 0
else echo "FAIL: behaviour differs from golden ($(grep -c '^[-+][^-+]' $H/diff.txt) changed lines); diff in $H/diff.txt"; head -40 "$H/diff.txt"; exit 1; fi
