#!/bin/bash
# Real-ds9 test of plugins/trails ("Remove Trails" button) on a synthetic field with two injected trails (one crossing a bright galaxy), then headless
# replay of the exported session.   scripts/verify_trails.sh [workdir]        env: DISPLAY_OVERRIDE=:77
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_tr.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home" "$W/sim"
cd "$ROOT" || exit 2
"$PY" - "$W/sim" > "$W/sim.txt" 2>&1 <<'PYEOF' || { cat "$W/sim.txt"; exit 2; }
import sys, json
sys.path.insert(0, '.')
import numpy as np
from astropy.io import fits
from ogfkit import trailsynth as S
img, truth = S.test_field(700, seed=5, sky=100.0, sigma=3.0, n_gal=10, n_star=50)
g = max((t for t in truth if t[3] == 'galaxy'), key=lambda t: t[2])
nx = ny = 700
# trail 1 through the brightest galaxy: line through (gx, gy) at 35 degrees
import math
th = 35.0; rho = -(g[0] - (nx - 1) / 2.0) * math.sin(math.radians(th)) + (g[1] - (ny - 1) / 2.0) * math.cos(math.radians(th))
img = img + S.trail_image(img.shape, th, rho, 5, 4.0) + S.trail_image(img.shape, 110.0, 160.0, 3, 3.0, s0=-250, s1=200)
fits.writeto(sys.argv[1] + '/trails_test.fits', img.astype(np.float32), overwrite=True)
json.dump(dict(trail1=dict(theta=th, rho=rho, fwhm=5), trail2=dict(theta=110.0, rho=160.0, fwhm=3, s0=-250, s1=200), galaxy=list(map(float, g[:3]))), open(sys.argv[1] + '/truth.json', 'w'))
PYEOF
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_TR_OUT="$W/tr.txt" OGF_TR_DIR="$W/gui" OGF_TR_SIM="$W/sim" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 900 \
  bin/ds9 "$W/sim/trails_test.fits" -geometry 1300x950 -source scripts/verify_trails.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/tr.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR' "$W/tr.txt" | head -8; grep '^SUMMARY' "$W/tr.txt"
ok=1; grep -q '^SUMMARY failures=0' "$W/tr.txt" && ! grep -q '^FAIL' "$W/tr.txt" || ok=0
if [ -f "$W/gui/ogfinder_session.py" ]; then
  "$PY" "$W/gui/ogfinder_session.py" --mode replay --outdir "$W/replay" "$W/sim/trails_test.fits" > "$W/replay.txt" 2>&1; rc=$?
  nst=$(grep -c 'step [0-9]*/[0-9]*' "$W/replay.txt"); nw=$(grep -c 'WARNING: .*differs' "$W/replay.txt")
  echo "replay: rc=$rc steps=$nst differ-warnings=$nw"
  [ $rc = 0 ] && [ "$nst" -ge 1 ] || ok=0
else echo "no exported session script"; ok=0; fi
[ $ok = 1 ]
