"""Prints the validation numbers quoted in docs/daophot.md (synthetic crowded field with known truth)."""
import os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import synth, psfmodel as pm, daophot as dp

img, t = synth.star_field((600, 600), n=1200, mag_range=(15.5, 23.0), slope=0.3, noise=0.5, seed=3)
bkg, rms = pm.background(img)
out = {}
for order in (0, 1, 2):
    mdl, info = pm.build_psf_model(img, bkg=bkg, rms=rms, fwhm_prior=3.3, degree=order, snr_min=20)
    res, resid = dp.allstar(img, mdl, bkg=bkg, rms=rms, fwhm=3.3, n_iter=3, n_workers=8)
    S = res['stars']
    d, j = cKDTree(np.c_[[s['x'] for s in S], [s['y'] for s in S]]).query(np.c_[t['x'], t['y']])
    print('== PSF order %d (%d PSF stars): %d stars fitted, %d good, residual rms/noise %.3f' % (order, info['n_used'], res['n'], res['n_good'], np.std(resid[60:540, 60:540]) / 0.5))
    print(' mag        N  recov  compl  dF/F(med)  dF/F(rms)  dpos(med px)  chi(med)')
    for lo, hi in [(15.5, 17), (17, 18), (18, 19), (19, 20), (20, 21), (21, 22)]:
        m = (t['mag'] >= lo) & (t['mag'] < hi)
        mm = m & (d < 1.0)
        if mm.sum() == 0:
            print(' %.1f-%.1f %4d %5d' % (lo, hi, m.sum(), 0)); continue
        fr = np.array([S[k]['flux'] for k in j[mm]]) / t['flux'][mm] - 1
        print(' %.1f-%.1f %4d %5d  %.3f   %+.4f     %.4f      %.3f        %.2f' % (lo, hi, m.sum(), mm.sum(), mm.sum() / m.sum(), np.median(fr), 1.4826 * np.median(np.abs(fr - np.median(fr))), np.median(d[mm]), np.median([S[k]['chi'] for k in j[mm]])))
    print(' unmatched fitted stars:', len(S) - len(set(j[d < 1.5])), 'of', len(S))
