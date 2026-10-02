#!/bin/bash
# Run scripts/verify_newplugins.tcl in a real ds9 (scratch HOME), then replay the exported session script headless and require that every
# step's stdout and the catalog after every step are identical to the GUI run (the session template compares CRC32s itself).
#   scripts/verify_newplugins.sh [workdir]        env: OGF_NP_ONLY=isophote,...  DISPLAY_OVERRIDE=:77
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_np.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home"
cd "$ROOT" || exit 2
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_NP_OUT="$W/np.txt" OGF_NP_DIR="$W/gui" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 900 \
  bin/ds9 "${OGF_TEST_FITS:-/workspace/fits}/m51.fits" -geometry 1300x950 -source scripts/verify_newplugins.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/np.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR' "$W/np.txt" | head -8; grep '^SUMMARY' "$W/np.txt"
ok=1; grep -q '^SUMMARY failures=0' "$W/np.txt" && ! grep -q '^FAIL' "$W/np.txt" || ok=0
if [ -f "$W/gui/ogfinder_session.py" ]; then
  "$PY" "$W/gui/ogfinder_session.py" --mode replay --outdir "$W/replay" "${OGF_TEST_FITS:-/workspace/fits}/m51.fits" > "$W/replay.txt" 2>&1; rc=$?
  nst=$(grep -c 'step [0-9]*/[0-9]*' "$W/replay.txt")
  nid=$(grep -c 'catalog after step identical' "$W/replay.txt"); nsd=$(grep -c 'stdout identical' "$W/replay.txt")
  nw=$(grep -c 'WARNING: .*differs' "$W/replay.txt")
  echo "replay: rc=$rc steps=$nst stdout-identical=$nsd catalog-identical=$nid differ-warnings=$nw"
  [ $rc = 0 ] && [ "$nw" = 0 ] && [ "$nid" -ge 1 ] || ok=0
else echo "no exported session script"; ok=0; fi
[ $ok = 1 ]
