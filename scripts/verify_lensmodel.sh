#!/bin/bash
# Real-ds9 test of plugins/lensmodel on a synthetic lens, then headless replay of the exported session (catalog and stdout identical).
#   scripts/verify_lensmodel.sh [workdir]        env: DISPLAY_OVERRIDE=:77
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_lm.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home" "$W/sim"
cd "$ROOT" || exit 2
"$PY" plugins/lensmodel/lensmodel.py --task simulate --work "$W/sim" --seed 3 > "$W/sim.txt" 2>&1 || { cat "$W/sim.txt"; exit 2; }
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_LM_OUT="$W/lm.txt" OGF_LM_DIR="$W/gui" OGF_LM_SIM="$W/sim" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 600 \
  bin/ds9 "$W/sim/sim_lens.fits" -geometry 1300x950 -source scripts/verify_lensmodel.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/lm.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR' "$W/lm.txt" | head -8; grep '^SUMMARY' "$W/lm.txt"
ok=1; grep -q '^SUMMARY failures=0' "$W/lm.txt" && ! grep -q '^FAIL' "$W/lm.txt" || ok=0
if [ -f "$W/gui/ogfinder_session.py" ]; then
  "$PY" "$W/gui/ogfinder_session.py" --mode replay --outdir "$W/replay" "$W/sim/sim_lens.fits" > "$W/replay.txt" 2>&1; rc=$?
  nst=$(grep -c 'step [0-9]*/[0-9]*' "$W/replay.txt"); nid=$(grep -c 'catalog after step identical' "$W/replay.txt"); nw=$(grep -c 'WARNING: .*differs' "$W/replay.txt")
  echo "replay: rc=$rc steps=$nst catalog-identical=$nid differ-warnings=$nw"
  [ $rc = 0 ] && [ "$nw" = 0 ] && [ "$nid" -ge 1 ] || ok=0
else echo "no exported session script"; ok=0; fi
[ $ok = 1 ]
