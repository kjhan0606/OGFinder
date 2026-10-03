#!/usr/bin/env python3
"""Compare ogfkit.trails masks with acstools.findsat_mrt on jc8m32j5q_flc.fits chip 1 (needs `pip install acstools photutils` in a separate venv and the venv
python for ogfkit; step 1 writes the MRT mask, step 2 compares).
    ACS_PY compare_acstools.py mrt  FLC.fits MASK.npy     (run with the acstools python)
    OGF_PY compare_acstools.py cmp  FLC.fits MASK.npy     (run with the ogfkit python)"""
import sys, os
import numpy as np
mode, flc, mpath = sys.argv[1:4]
if mode == 'mrt':
    from acstools.findsat_mrt import WfcWrapper
    w = WfcWrapper(flc, extension=4, binsize=2, processes=2, preprocess=True, execute=True, plot=False, output_root=os.path.splitext(mpath)[0])
    np.save(mpath, w.mask)
else:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
    from astropy.io import fits
    from ogfkit import trails as T
    h = fits.open(flc)
    im = h[4].data.astype(np.float32); dq = h[6].data
    res = T.detect_trails(im, valid=(dq & 0x3FF) == 0)
    m = np.kron(np.load(mpath).astype(np.uint8), np.ones((2, 2), np.uint8)) > 0       # MRT mask is binned by 2
    for t in res['trails']:
        o = T.trail_mask(im.shape, [t]); mt = m & T.trail_mask(im.shape, [t], extra_margin=40)
        i = (o & mt).sum()
        print('trail %d (theta %.1f, FWHM %.1f): ours %d px, MRT %d px nearby, IoU %.2f, ours inside MRT %.2f, MRT covered by ours %.2f' % (
            t['id'], t['theta_deg'], t['fwhm'], o.sum(), mt.sum(), i / (o | mt).sum(), i / o.sum(), i / max(mt.sum(), 1)))
