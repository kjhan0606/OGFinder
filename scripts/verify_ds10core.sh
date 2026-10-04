#!/bin/bash
# Real-ds9 test of plugins/ds10core (stand-alone shell of the shared ds10 core): chip + menu entries, header calibration, region statistics / mask, pixel table,
# multi-band forced photometry on synthetic bands, ds10-script run, session records; optional screenshots (docs/shots/standalone_*.png).
#   scripts/verify_ds10core.sh [workdir]        env: DISPLAY_OVERRIDE=:77  OGF_SHOTS=docs/shots (screenshots on when set)  OGF_DS_GEOM=1300x950 (ds9 window)
# SKIP (exit 77) when the ds10core package of ds10-web cannot be found.
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_ds10.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home" "$W/in"
cd "$ROOT" || exit 2
"$PY" plugins/ds10core/ds10.py --where >/dev/null 2>&1 || { echo "ds10core not found (set DS10_CORE)"; exit 77; }
"$PY" - "$W/in" > "$W/sim.txt" 2>&1 <<'PYEOF' || { cat "$W/sim.txt"; exit 2; }
import sys, math
import numpy as np
from astropy.io import fits
d = sys.argv[1]
rng = np.random.default_rng(11); n = 400
yy, xx = np.mgrid[0:n, 0:n]
src = [(50 + 40 * i, 60 + 37 * ((i * 3) % 8), 10.0 + 4.0 * i) for i in range(9)]
spec = {"F105W": (1.0, 1.0551e-19, 10551.0), "F125W": (1.4, 2.2486e-20, 12486.0), "F160W": (2.0, 1.9276e-20, 15369.0)}
for lab, (scale, flam, plam) in spec.items():
    a = rng.normal(0.0, 0.02, (n, n))
    for cx, cy, f in src:
        a += scale * f / (2 * math.pi * 1.8 ** 2) * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.8 ** 2)) * 0.05
    h = fits.PrimaryHDU(a.astype(np.float32))
    h.header.update({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRPIX1": n / 2 + .5, "CRPIX2": n / 2 + .5, "CRVAL1": 53.1, "CRVAL2": -27.8, "CD1_1": -0.06 / 3600, "CD2_2": 0.06 / 3600,
                     "CD1_2": 0.0, "CD2_1": 0.0, "FILTER": lab, "BUNIT": "ELECTRONS/S", "PHOTFLAM": flam, "PHOTPLAM": plam, "PHOTZPT": -21.1, "EXPTIME": 1000.0, "INSTRUME": "WFC3"})
    h.writeto(f"{d}/syn_{lab.lower()}.fits", overwrite=True)
open(f"{d}/r.reg", "w").write("# Region file format: DS9 version 4.1\nimage\ncircle(%g,%g,12.3)\nbox(200.3,200.2,40,24,0)\nannulus(%g,%g,20.3,32.3)\n" % (src[0][0], src[0][1], src[0][0], src[0][1]))
PYEOF
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_DS_OUT="$W/ds.txt" OGF_DS_DIR="$W/in" OGF_DS_SHOTS="${OGF_SHOTS:-}" OGF_DS_DISPLAY="${DISPLAY_OVERRIDE:-:77}" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 600 \
  bin/ds9 "$W/in/syn_f160w.fits" -geometry "${OGF_DS_GEOM:-1300x950}" -source scripts/verify_ds10core.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/ds.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR' "$W/ds.txt" | head -8; grep '^SUMMARY' "$W/ds.txt"
grep -q '^SUMMARY failures=0' "$W/ds.txt" && ! grep -q '^FAIL' "$W/ds.txt"
