#!/bin/bash
# Real-ds9 test of the spectra "Fit Kinematics" step and the Kinematics Viewer on synthetic IFU cube / long-slit data with known truth.
#   scripts/verify_spectra_kin.sh [workdir]        env: DISPLAY_OVERRIDE=:77
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
W="${1:-$(mktemp -d /tmp/ogf_sk.XXXXXX)}"; rm -rf "$W"; mkdir -p "$W/home" "$W/spec"
cd "$ROOT" || exit 2
"$PY" - "$W/spec" <<'PYEOF' || exit 2
import sys, numpy as np
sys.path.insert(0, '.')
from astropy.io import fits
from ogfkit import spectrasynth as SS
d = sys.argv[1]
rng = np.random.default_rng(1)
fits.PrimaryHDU(rng.normal(100, 3, (120, 120)).astype(np.float32)).writeto(d + '/../img.fits')
w, s, t = SS.disc_slit(seed=3, noise=0.05, vc=220.0, rt=3.5, inc=60.0, pa=40.0, sigma0=45.0, z=0.12)
h = fits.PrimaryHDU(s.astype(np.float32)); h.header['CRVAL1'] = w[0]; h.header['CDELT1'] = w[1] - w[0]; h.header['CRPIX1'] = 1; h.writeto(d + '/spec_1.fits')
w, c, t = SS.disc_cube(seed=4, noise=0.05, pa=130.0, inc=55.0, vc=200.0, rt=4.0, sigma0=40.0, z=0.12)
h = fits.PrimaryHDU(c.astype(np.float32)); h.header['CRVAL3'] = w[0]; h.header['CDELT3'] = w[1] - w[0]; h.header['CRPIX3'] = 1; h.writeto(d + '/spec_2.fits')
open(d + '/../cat.tsv', 'w').write('NUMBER\tX_IMAGE\tY_IMAGE\tZ\n1\t30\t30\t0.12\n2\t21\t21\t0.12\n3\t90\t90\t0.12\n')
PYEOF
rm -rf ~/ds9.auto ~/ds9.auto.dir
HOME="$W/home" OGF_SK_OUT="$W/sk.txt" OGF_SK_CAT="$W/cat.tsv" OGF_SK_SPEC="$W/spec" DISPLAY="${DISPLAY_OVERRIDE:-:77}" timeout -s KILL 600 \
  bin/ds9 "$W/img.fits" -geometry 1300x950 -source scripts/verify_spectra_kin.tcl > "$W/raw.txt" 2>&1
grep -c '^PASS' "$W/sk.txt" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^BGERROR' "$W/sk.txt" | head -8; grep '^SUMMARY' "$W/sk.txt"
grep -q '^SUMMARY failures=0' "$W/sk.txt" && ! grep -q '^FAIL' "$W/sk.txt"
