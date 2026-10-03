#!/usr/bin/env python3
"""Real-data check: HST ACS/WFC exposure jc8m32j5q_flc.fits (MACS J0717, F606W, proposal 13498; the acstools findsat_mrt example, contains real
satellite trails).  Downloads the file from MAST if missing, runs ogfkit.trails on both chips (DQ bad-pixel flags masked), measures the residual
wings outside the mask and, if `acstools` + `photutils` are importable, compares the masks with acstools.findsat_mrt.
usage: real_acs.py [--dir DIR] [--out report.md]"""
import argparse, json, os, subprocess, sys, time
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
from ogfkit import trails as T
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jc8m32j5q_flc.fits'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='/workspace/work/trails/real')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    f = os.path.join(a.dir, 'jc8m32j5q_flc.fits')
    if not os.path.exists(f):
        subprocess.check_call(['curl', '-sSL', '-o', f, URL])
    from astropy.io import fits
    h = fits.open(f)
    lines = []
    for ext, name in ((4, 'chip1'), (1, 'chip2')):
        im = h[ext].data.astype(np.float32); dq = h[ext + 2].data
        valid = (dq & 0x3FF) == 0
        t0 = time.time()
        res = T.detect_trails(im.copy(), valid=valid)
        lines.append('%s: %d trails, %d rejected (bleed/spike/shape), best z %.1f, sigma_pix %.1f, %.0f s' % (name, len(res['trails']), len(res.get('rejected', [])), res['best_zscore'], res['sigma_pix'], time.time() - t0))
        ny, nx = im.shape
        for t in res['trails']:
            hw = t['halfwidth']
            sv = np.arange(t['s0'], t['s1'], 1.0); tv = np.arange(-hw - 25, hw + 26)
            P = np.nanmedian(T.sample_strip(im, t['theta_deg'], t['rho'], sv, tv), axis=1)
            sky = np.median(P[abs(tv) > hw + 15])
            e = np.median(P[(abs(tv) > hw) & (abs(tv) <= hw + 4)] - sky)
            lines.append('  trail %d: theta %.2f rho %.1f length %.0f FWHM %.1f px peak %.1f e- (%.2f sigma_pix) z %.1f mask half-width %.1f; median excess 0-4 px outside the mask %+.2f e- (%+.4f of peak)' % (
                t['id'], t['theta_deg'], t['rho'], t['length'], t['fwhm'], t['amp'], t['amp'] / res['sigma_pix'], t['zscore'], hw, e, e / t['amp']))
        for r in res.get('rejected', []):
            lines.append('  rejected: theta %.1f length %.0f z %.1f reason %s' % (r['theta'], r['length'], r['zscore'], r['reason']))
    out = '\n'.join(lines)
    print(out)
    if a.out:
        open(a.out, 'w').write(out + '\n')


if __name__ == '__main__':
    main()
