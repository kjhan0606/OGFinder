"""Real-data part of psf_validate.py: held-out star residuals of the PSF models and fitted wings on HUDF F160W (crop) and M51."""
import json
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
from ogfkit import psfmodel as PM, psfext as PX  # noqa: E402


def heldout_residual(mdl, sub, test_xy, r=8.0):
    """For each held-out star: centroid (windowed), linear flux fit, sum |data - flux * model| inside r px / flux.  -> array"""
    import sep
    d = np.ascontiguousarray(sub.astype(np.float64))
    out = []
    ny, nx = d.shape
    h = int(r) + 3
    for x, y in test_xy:
        xw, yw, _ = sep.winpos(d, [x], [y], [1.5], subpix=5)
        x, y = float(xw[0]), float(yw[0])
        ix, iy = int(round(x)), int(round(y))
        if ix - h < 0 or iy - h < 0 or ix + h >= nx or iy + h >= ny:
            continue
        st = d[iy - h:iy + h + 1, ix - h:ix + h + 1]
        m = mdl.stamp(x, y, x - ix, y - iy, 2 * h + 1)
        yy, xx = np.mgrid[-h:h + 1, -h:h + 1]
        sel = np.hypot(xx, yy) <= r
        f = (st[sel] * m[sel]).sum() / (m[sel] ** 2).sum()
        out.append(np.abs(st[sel] - f * m[sel]).sum() / max(f, 1e-30))
    return np.array(out)


def one_field(name, data, fwhm, snr_min, rmax_w, size, sat=None, seed=0):
    bkg, rms = PM.background(data)
    sub = np.nan_to_num(data) - bkg
    mdl0, info = PM.build_psf_model(data, bkg=bkg, rms=rms, fwhm_prior=fwhm, snr_min=snr_min, degree=2, size=size)
    st = [s for s in info['stars'] if s['used']]
    xy = np.array([[s['x'], s['y']] for s in st])
    res = dict(name=name, n_candidates=info['n_candidates'], n_used_all=len(st), models={})
    if len(xy) < 12:
        res['error'] = 'too few stars'
        return res
    rng = np.random.default_rng(seed)
    scores = {}
    for rep in range(4):                                              # 4 random half splits
        perm = rng.permutation(len(xy))
        tr, te = xy[perm[:len(xy) // 2]], xy[perm[len(xy) // 2:]]
        for label, deg, rank in (('deg0', 0, None), ('deg1', 1, None), ('deg2', 2, None), ('deg3', 3, None), ('deg3_rank1', 3, 1), ('deg3_rank2', 3, 2)):
            m, _ = PM.build_psf_model(data, bkg=bkg, rms=rms, xy=tr, fwhm_prior=fwhm, degree=deg, size=size)
            if rank:
                m = PX.reduce_rank(m, rank)
            scores.setdefault(label, []).extend(heldout_residual(m, sub, te).tolist())
    res['models'] = {k: dict(median=float(np.median(v)), mean=float(np.mean(v)), n=len(v)) for k, v in scores.items()}
    # wings of the brightest isolated stars
    flux = np.array([s['flux'] for s in st])
    bright = xy[np.argsort(-flux)[:20]]
    r, prof, ns = PX.measure_wings(sub, bright, rcore=6.0 if fwhm > 2.5 else 3.0, rmax=rmax_w, saturation=sat, min_dist=rmax_w * 0.7)
    w = PX.fit_wings(r, prof, rmin=0.4 * rmax_w, rmax=rmax_w)
    res['wings'] = dict(n_stars=ns, fit=w, flux_beyond_stamp_per_core_flux=float(PX.wing_flux(w, size / 2.0)) if w and w['valid'] else None,
                        profile=[(float(a), float(b)) for a, b in zip(r[::3], prof[::3]) if np.isfinite(b)])
    return res


def real(a):
    out = []
    with fits.open('/workspace/fits/hudf_f160w.fits') as h:
        d = np.array(h[0].data[800:2800, 800:2800], float)
    out.append(one_field('HUDF F160W crop 2000x2000', d, 3.0, 40.0, 30.0, 31))
    with fits.open('/workspace/fits/m51.fits') as h:
        d = np.array(h[0].data, float)
    out.append(one_field('M51 600x600', d, 2.0, 15.0, 12.0, 15))
    json.dump(out, open(a.out, 'w'), indent=1)
    for r in out:
        print(r['name'], r.get('n_used_all'), {k: round(v['median'], 4) for k, v in r.get('models', {}).items()}, r.get('wings', {}).get('fit'), r.get('error'))
