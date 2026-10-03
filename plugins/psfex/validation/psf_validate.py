#!/usr/bin/env python3
"""Validation of the PSF extensions (ogfkit/psfext.py) and of the PSF builder on injected stars and on real HUDF / M51 stars.

    psf_validate.py synth OUT.json [--seeds 4]     spatially varying Moffat PSF (FWHM 3.0 -> 4.0 px, e 0.02 -> 0.14, PA rotating), ~150 stars + 15 % star-like galaxies + 10 % blends, noise rms 3:
                                                   accuracy of the model against the truth PSF at 100 positions for polynomial degree 0-3 and reduced rank, contaminant rejection,
                                                   halo (second Moffat, 2 % of the flux) recovery with measure_wings / fit_wings, error budget on galaxy fits
    psf_validate.py real OUT.json                  HUDF F160W (crop) and M51: held-out star residual for the models, fitted wings of bright stars
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('OMP_NUM_THREADS', '2')
from ogfkit import models as M, psfmodel as PM, psfext as PX  # noqa: E402

N = 800
SIZE = 31


def truth_params(x, y):
    u, v = x / N, y / N
    return dict(fwhm=3.0 + 1.0 * u, e=0.02 + 0.12 * v, pa=20.0 + 50.0 * u)


def truth_psf(x, y, size=SIZE, dx=0.0, dy=0.0, halo=0.0):
    p = truth_params(x, y)
    core = M.moffat_psf(size, p['fwhm'], beta=2.5, e=p['e'], pa_deg=p['pa'], dx=dx, dy=dy)
    if halo > 0:
        h = M.moffat_psf(size, 12.0, beta=1.8, e=0.0, dx=dx, dy=dy, nsub=1)
        return (1 - halo) * core + halo * h, core
    return core


def inject(img, x, y, flux, size, halo=0.0):
    ix, iy = int(round(x)), int(round(y))
    h = size // 2
    p = truth_psf(x, y, size, x - ix, y - iy, halo)
    p = p[0] if isinstance(p, tuple) else p
    y0, x0 = iy - h, ix - h
    ya, xa = max(0, y0), max(0, x0)
    yb, xb = min(img.shape[0], y0 + size), min(img.shape[1], x0 + size)
    img[ya:yb, xa:xb] += flux * p[ya - y0:yb - y0, xa - x0:xb - x0]


def make_field(seed, halo=0.0, nstars=150, big=SIZE, logf=(3.3, 5.3)):
    rng = np.random.default_rng(seed)
    img = rng.normal(0, 3.0, (N, N)) + 100.0
    stars = []
    for _ in range(nstars):
        x, y = rng.uniform(25, N - 25, 2)
        f = 10 ** rng.uniform(*logf)
        inject(img, x, y, f, big, halo)
        stars.append((x, y, f))
    kinds = []
    for _ in range(int(0.15 * nstars)):                                # star-like galaxies
        x, y = rng.uniform(25, N - 25, 2)
        f = 10 ** rng.uniform(3.5, 5.0)
        g = M.render_sersic((31, 31), 15 + (x - int(x)), 15 + (y - int(y)), M.sersic_Ie_from_flux(f, rng.uniform(1.8, 3.0), 1.0, rng.uniform(0.5, 0.9)), rng.uniform(1.8, 3.0), 1.0, rng.uniform(0.5, 0.9), rng.uniform(0, 180))
        ix, iy = int(x), int(y)
        img[iy - 15:iy + 16, ix - 15:ix + 16] += g
        kinds.append(('gal', x, y))
    for (x, y, f) in stars[:int(0.1 * nstars)]:                          # blends: companion 5-8 px away
        a = rng.uniform(0, 2 * math.pi)
        r = rng.uniform(5, 8)
        inject(img, x + r * math.cos(a), y + r * math.sin(a), f * rng.uniform(0.3, 1.0), big, halo)
        kinds.append(('blend', x, y))
    return img, stars, kinds


def eval_model(mdl, pts, size=SIZE):
    err_peak, err_l1, dfw, de = [], [], [], []
    for x, y in pts:
        a = mdl.stamp(x, y, 0.0, 0.0, size)
        b = truth_psf(x, y, size)
        err_peak.append(np.abs(a - b).max() / b.max())
        err_l1.append(np.abs(a - b).sum())
        sa, sb = PM.psf_shape(a), PM.psf_shape(b)
        dfw.append(sa['fwhm'] - sb['fwhm'])
        de.append(sa['e'] - sb['e'])
    return dict(max_over_peak_median=float(np.median(err_peak)), max_over_peak_p90=float(np.percentile(err_peak, 90)), l1_median=float(np.median(err_l1)),
                fwhm_err_rms=float(np.sqrt(np.mean(np.square(dfw)))), e_err_rms=float(np.sqrt(np.mean(np.square(de)))))


def synth(a):
    rng = np.random.default_rng(0)
    pts = list(zip(rng.uniform(40, N - 40, 100), rng.uniform(40, N - 40, 100)))
    out = dict(models={}, selection={}, wings={}, budget={})
    acc = {}
    sel = dict(stars_total=0, stars_used=0, gal_total=0, gal_used=0, blend_total=0, blend_used=0)
    for seed in range(a.seeds):
        img, stars, kinds = make_field(seed)
        cfgs = [('deg0', dict(degree=0)), ('deg1', dict(degree=1)), ('deg2', dict(degree=2)), ('deg3', dict(degree=3)), ('deg3_rank1', dict(degree=3, rank=1)), ('deg3_rank2', dict(degree=3, rank=2)),
                ('deg4_rank2', dict(degree=4, rank=2))]
        for name, c in cfgs:
            mdl, info = PM.build_psf_model(img, fwhm_prior=3.5, snr_min=25, degree=c['degree'], size=SIZE)
            if 'rank' in c:
                mdl = PX.reduce_rank(mdl, c['rank'])
            acc.setdefault(name, []).append(eval_model(mdl, pts))
            if name == 'deg2':
                used = [(s['x'], s['y']) for s in info['stars'] if s['used']]
                allst = [(s['x'], s['y']) for s in info['stars']]
                for lab, items in (('stars', [(x, y) for x, y, f in stars]), ('gal', [(x, y) for k, x, y in kinds if k == 'gal']), ('blend', [(x, y) for k, x, y in kinds if k == 'blend'])):
                    for (x, y) in items:
                        near = [p for p in allst if math.hypot(p[0] - x, p[1] - y) < 2.0]
                        if near:
                            sel[lab + '_total'] += 1
                            if any(math.hypot(p[0] - x, p[1] - y) < 2.0 for p in used):
                                sel[lab + '_used'] += 1
    for name, L in acc.items():
        out['models'][name] = {k: float(np.mean([d[k] for d in L])) for k in L[0]}
    out['selection'] = sel
    # halo recovery
    img, stars, kinds = make_field(11, halo=0.08, nstars=60, big=121, logf=(5.0, 6.3))
    bright = [(x, y) for x, y, f in stars if f > 2e5]
    r, prof, ns = PX.measure_wings(img, bright, bkg=100.0, rcore=6.0, rmax=45.0)
    w = PX.fit_wings(r, prof, rmin=8.0, rmax=40.0)
    # truth halo (second Moffat) fraction beyond r > 8 px (per unit core flux)
    hal = M.moffat_psf(241, 12.0, beta=1.8, nsub=1)
    yy, xx = np.mgrid[-120:121, -120:121]
    rr = np.hypot(xx, yy)
    pt = truth_psf(400, 400, 241, halo=0.08)[0]
    fc = pt[rr < 6].sum()
    true_beyond = float(pt[(rr > 8) & (rr <= 120)].sum() / fc)
    out['wings'] = dict(n_stars=ns, fit=w, true_halo_flux_beyond_8px_per_core_flux=float(true_beyond), fitted_wing_flux_beyond_8px_to_121_per_core_flux=float(PX.wing_flux(w, 8.0, 121.0)) if w else None)
    mdl, _ = PM.build_psf_model(img, fwhm_prior=3.5, snr_min=25, degree=1, size=31)
    ext, inf = PX.extended_psf(mdl, w, size=121)
    ytrue = truth_psf(N / 2, N / 2, 121, halo=0.08)[0]
    ytrue = ytrue / ytrue.sum()
    def curve(p):
        rr = np.hypot(*(np.mgrid[:p.shape[0], :p.shape[1]] - p.shape[0] // 2))
        return [float(p[rr <= R].sum()) for R in (6, 12, 20, 30, 40, 55)]
    core = mdl.stamp(N / 2, N / 2, 0, 0)
    big = np.zeros((121, 121)); big[45:76, 45:76] = core
    out['wings']['enclosed_flux_R=6,12,20,30,40,55'] = dict(truth=curve(ytrue), core_only=curve(big / big.sum()), extended=curve(ext))
    out['wings']['extended_info'] = inf
    # error budget with the true halo PSF for big galaxies
    ptrue = ytrue[40:81, 40:81] / ytrue[40:81, 40:81].sum()
    pext = ext[40:81, 40:81] / ext[40:81, 40:81].sum()
    pcore = core / core.sum()
    out['budget'] = dict(core31_vs_true=PX.error_budget(ptrue, pcore, profiles=((1.0, 8.0), (4.0, 8.0))), extended_vs_true=PX.error_budget(ptrue, pext, profiles=((1.0, 8.0), (4.0, 8.0))))
    base = truth_psf(400, 400, 31)
    cases = {'fwhm+5%': M.moffat_psf(31, truth_params(400, 400)['fwhm'] * 1.05, beta=2.5, e=truth_params(400, 400)['e'], pa_deg=truth_params(400, 400)['pa']),
             'fwhm-5%': M.moffat_psf(31, truth_params(400, 400)['fwhm'] * 0.95, beta=2.5, e=truth_params(400, 400)['e'], pa_deg=truth_params(400, 400)['pa']),
             'e+0.05': M.moffat_psf(31, truth_params(400, 400)['fwhm'], beta=2.5, e=truth_params(400, 400)['e'] + 0.05, pa_deg=truth_params(400, 400)['pa']),
             'shift 0.3px': M.moffat_psf(31, truth_params(400, 400)['fwhm'], beta=2.5, e=truth_params(400, 400)['e'], pa_deg=truth_params(400, 400)['pa'], dx=0.3),
             'stamp 15px': base[8:23, 8:23] / base[8:23, 8:23].sum()}
    out['budget']['mismatch'] = {k: PX.error_budget(base, v, profiles=((1.0, 4.0), (4.0, 4.0))) for k, v in cases.items()}
    json.dump(out, open(a.out, 'w'), indent=1)
    print(json.dumps({k: {m: round(v['max_over_peak_median'], 4) for m, v in out['models'].items()} for k in ['max|d|/peak median']}))
    print('selection', sel)
    print('wings', ns, w and {k: round(v, 3) for k, v in w.items() if not isinstance(v, list)}, 'true beyond8 %.4f' % true_beyond)
    print('enclosed', out['wings']['enclosed_flux_R=6,12,20,30,40,55'])


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['synth', 'real'])
    ap.add_argument('out')
    ap.add_argument('--seeds', type=int, default=4)
    a = ap.parse_args()
    if a.mode == 'synth':
        synth(a)
    else:
        from psf_real import real
        real(a)


if __name__ == '__main__':
    main()
