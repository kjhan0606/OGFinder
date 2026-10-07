#!/bin/bash
# Real-ds9 test of the shared-core star aperture photometry steps of plugins/ds10core (aper-phot, aper-series): stars picked on the image
# (real X clicks with xdotool when available, otherwise point regions created in the frame), the run with the GUI parameters, the catalogue
# table, the aperture + sky-annulus overlay, and a light curve over a synthetic frame sequence with an injected variable star.
#   scripts/verify_aper_phot.sh [workdir]     env: DISPLAY_OVERRIDE=:77  OGF_SHOTS=dir (screenshots aper_desktop_*.png when set)  OGF_DS_GEOM=1300x950
# SKIP (exit 77) when the ds10core package of Astrafex Web cannot be found.
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_aper.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home" "$W/in"
cd "$ROOT" || exit 2
"$PY" plugins/ds10core/ds10.py --where >/dev/null 2>&1 || { echo "ds10core not found (set DS10_CORE)"; exit 77; }
"$PY" - "$W/in" > "$W/sim.txt" 2>&1 <<'PYEOF' || { cat "$W/sim.txt"; exit 2; }
import sys, math
import numpy as np
from astropy.io import fits
d = sys.argv[1]
rng = np.random.default_rng(7); ny, nx = 360, 480
yy, xx = np.mgrid[0:ny, 0:nx]
pts = []
while len(pts) < 18:
    p = rng.uniform([40, 40], [nx - 40, ny - 40])
    if all(math.hypot(*(p - q)) > 45 for q in pts):
        pts.append(p)
xy = np.array(pts); flux = 10 ** rng.uniform(3.8, 5.0, len(xy)); flux[:4] = [6e4, 8e4, 7e4, 5e4]
with open(f"{d}/stars.txt", "w") as fh:
    for (x, y), f in zip(xy, flux):
        fh.write(f"{x + 1:.3f} {y + 1:.3f} {f:.1f}\n")                       # 1-based image pixels
for k in range(6):
    f = flux.copy(); f[0] *= 1 + 0.08 * math.sin(2 * math.pi * k / 6)
    s = (3.0 + 0.1 * k) / 2.3548; dx = 0.6 * k
    img = np.full((ny, nx), 150.0)
    for (x, y), fl in zip(xy, f):
        img += fl / (2 * math.pi * s * s) * np.exp(-((xx - x - dx) ** 2 + (yy - y) ** 2) / (2 * s * s))
    img = rng.poisson(img * 2.0) / 2.0 + rng.normal(0, 3, img.shape)
    h = fits.PrimaryHDU(img.astype(np.float32))
    h.header.update({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRPIX1": nx / 2 + .5, "CRPIX2": ny / 2 + .5, "CRVAL1": 210.0, "CRVAL2": 20.0, "CD1_1": -0.5 / 3600,
                     "CD2_2": 0.5 / 3600, "CD1_2": 0.0, "CD2_1": 0.0, "GAIN": 2.0, "RDNOISE": 6.0, "EXPTIME": 60.0, "FILTER": "V", "MAGZERO": 25.0,
                     "MJD-OBS": 60600.2 + 0.01 * k, "AIRMASS": 1.1 + 0.02 * k, "OBJECT": "SYN-VAR"})
    h.writeto(f"{d}/seq_{k}.fits", overwrite=True)
PYEOF
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_DS_OUT="$W/ds.txt" OGF_DS_DIR="$W/in" OGF_DS_SHOTS="${OGF_SHOTS:-}" OGF_DS_DISPLAY="${DISPLAY_OVERRIDE:-:77}" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 600 \
  bin/ds9 "$W/in/seq_0.fits" -geometry "${OGF_DS_GEOM:-1300x950}" -source scripts/verify_aper_phot.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/ds.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR|^INFO' "$W/ds.txt" | head -12; grep '^SUMMARY' "$W/ds.txt"
grep -q '^SUMMARY failures=0' "$W/ds.txt" && ! grep -q '^FAIL' "$W/ds.txt"
