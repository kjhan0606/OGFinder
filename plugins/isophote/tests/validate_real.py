#!/usr/bin/env python3
"""Real-data validation of the isophote plugin (not a pytest: prints numbers).  HUDF F160W crop + M51.
For each of the brightest extended sources of the crop (SEP detection) the isophote result is compared with independent SEP numbers:
  * total flux inside the last isophote vs SEP elliptical-aperture flux of the same ellipse (same centre, SMA, eps, PA of that isophote)
  * eps / PA at the half-light radius vs the second-moment ellipse of SEP (A_IMAGE, B_IMAGE, THETA_IMAGE)
usage: validate_real.py [--n 6]"""
import argparse, os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
import isophote as I
from ogfkit import imageio

ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=6); a = ap.parse_args()
import sep
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')
data, hdr = imageio.load_image(FITS + '/hudf_f160w.fits')
x0, y0, W = 1500, 1500, 1000
cut = np.ascontiguousarray(data[y0:y0 + W, x0:x0 + W])
bkg = sep.Background(cut)
sub = cut - bkg.back()
objs = sep.extract(sub, 2.0, err=bkg.globalrms, minarea=30)
flux, ferr, fl = sep.sum_ellipse(sub, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], 2.5, err=bkg.globalrms)
order = np.argsort(-flux)
zp = 25.0
print('HUDF F160W crop %dx%d at (%d,%d): %d sources, rms %.4g' % (W, W, x0, y0, len(objs), bkg.globalrms))
print('%4s %7s %7s | %6s %6s %6s %6s | %6s %6s | %8s %8s %7s' % ('id', 'x', 'y', 'e_iso', 'e_sep', 'pa_iso', 'pa_sep', 'R50', 'rconv', 'F_iso', 'F_sep', 'dF/F'))
dFs, des, dpas = [], [], []
nok = 0
for i in order:
    if nok >= a.n:
        break
    o = objs[i]
    if o['a'] < 6 or o['b'] / o['a'] < 0.3:
        continue
    src = dict(num=int(i), x=float(o['x']), y=float(o['y']), a=float(o['a']), b=float(o['b']), theta=float(np.rad2deg(o['theta'])), kron=3.5, mag=-flux[i])
    res, model, org = I.fit_source(sub, src, maxsma_scale=2.0, step=0.12, zp=zp, mask_neighbours=2.5,
                                   others=[dict(num=int(j), x=float(objs['x'][j]), y=float(objs['y'][j]), a=float(objs['a'][j]), b=float(objs['b'][j]), theta=float(np.rad2deg(objs['theta'][j])), kron=3, mag=0) for j in range(len(objs))])
    if res['status'] != 'ok':
        print('%4d  fit %s' % (i, res['status']))
        continue
    s = res['summary']; t = res['table']
    good = [k for k, c in enumerate(t['stop_code']) if c == 0]
    k = good[-1]
    sma, eps, pa = t['sma'][k], t['ellipticity'][k], t['pa'][k]
    xc, yc = t['x0'][k], t['y0'][k]
    f_sep, _, _ = sep.sum_ellipse(sub, [xc], [yc], [sma], [sma * (1 - eps)], [np.deg2rad(((pa + 90) % 180) - 90)], 1.0)
    f_iso = t['growth_flux'][k]
    dF = f_iso / f_sep[0] - 1
    de, dpa = s['eps_hl'] - (1 - o['b'] / o['a']), ((s['pa_hl'] - np.rad2deg(o['theta']) + 90) % 180) - 90
    print('%4d %7.1f %7.1f | %6.3f %6.3f %6.1f %6.1f | %6.2f %6.1f | %8.1f %8.1f %7.3f' % (i, o['x'], o['y'], s['eps_hl'], 1 - o['b'] / o['a'],
          s['pa_hl'], np.rad2deg(o['theta']) % 180, s['r50'], s['r_conv'], f_iso, f_sep[0], dF))
    dFs.append(dF); des.append(de); dpas.append(dpa); nok += 1
print('galaxies %d: growth-curve flux vs SEP elliptical aperture: median dF/F %.4f, max |dF/F| %.4f' % (len(dFs), np.median(dFs), np.max(np.abs(dFs))))
des = [d for d in des if np.isfinite(d)]; dpas = [d for d in dpas if np.isfinite(d)]
print('eps(iso, r50) - eps(SEP moments): median %.3f rms %.3f ;  PA diff: median %.1f deg rms %.1f deg' % (np.median(des), np.std(des), np.median(dpas), np.std(dpas)))
