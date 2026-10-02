#!/usr/bin/env python3
"""Synthetic-truth validation of ogfkit/lensmodel.py  (python3 run_validation.py OUT.json [n_lens])

1. noise-free recovery of SIE+shear quads (should be exact)
2. Monte-Carlo with Gaussian image-position noise 3 mas and 10 mas: recovered theta_E, q, PA, shear, source position, predicted-vs-true images
3. misspecification: truth has an extra PIEMD perturber near an image; fit with and without the perturber in the model
4. time delays and counter-images against the true model for the fitted parameters
"""
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import lensmodel as L, lenssynth as S  # noqa: E402


def stats(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return dict(n=int(v.size), median=float(np.median(v)), mean=float(np.mean(v)), std=float(np.std(v)), p16=float(np.percentile(v, 16)), p84=float(np.percentile(v, 84)),
                max_abs=float(np.max(np.abs(v)))) if v.size else {}


def summarize(res):
    q = [r for r in res if r['quad'] and r['fit'].get('success')]
    d = dict(n_quads=len(q), n_doubles=len([r for r in res if not r['quad']]))
    dth = [(r['fit']['params']['theta_E'] - r['truth']['theta_E']) / r['truth']['theta_E'] for r in q]
    dq = [r['fit']['params']['q'] - r['truth']['q'] for r in q]
    dphi = [S.angle_diff(r['fit']['params']['phi'], r['truth']['phi']) for r in q if r['truth']['q'] < 0.9]
    dg = [r['fit']['params']['gamma'] - r['truth']['gamma'] for r in q]
    dsrc = [math.hypot(r['fit']['source'][0] - r['source'][0], r['fit']['source'][1] - r['source'][1]) for r in q]
    rms = [r['fit']['rms_arcsec'] for r in q]
    chi = [r['fit']['chi2'] for r in q]
    d.update(dtheta_E_rel=stats(dth), dq=stats(dq), dphi_deg=stats(dphi), dgamma=stats(dg), dsource_arcsec=stats(dsrc), rms_arcsec=stats(rms), chi2=stats(chi))
    dbl = [r for r in res if not r['quad'] and r['fit'].get('success')]
    if dbl:
        d['doubles_dtheta_E_rel'] = stats([(r['fit']['params']['theta_E'] - r['truth']['theta_E']) / r['truth']['theta_E'] for r in dbl])
    return d


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    R = {}
    t0 = time.time()
    # 1. noise free
    res = S.validation_sample(n_lens=min(n, 20), noise_arcsec=0.0, quad_fraction=1.0, seed=11)
    # noise 0 -> sigma 0 would divide by zero: validation_sample uses sigma=noise; handle by tiny sigma
    R['noise_free'] = summarize(res)
    for sig in (0.003, 0.010):
        res = S.validation_sample(n_lens=n, noise_arcsec=sig, quad_fraction=0.75, seed=int(sig * 1e4))
        R['noise_%gmas' % (sig * 1000)] = summarize(res)
        R['noise_%gmas' % (sig * 1000)]['expected_rms_for_quad'] = 'sigma * sqrt(dof/ndata) approx %.4f arcsec' % (sig * math.sqrt(1 / 8.0))
        if sig == 0.003:
            # 4. counter-images / delays for the fit vs truth model
            dd = []
            cnt = []
            for r in res:
                if not r['quad'] or not r['fit'].get('success'):
                    continue
                tm = L.build_sie_shear(r['truth'])
                fm = L.build_sie_shear(r['fit']['params'])
                tb = r['source']
                fb = r['fit']['source']
                ti = tm.solve_images(tb[0], tb[1], half=3 * r['truth']['theta_E'], n=301)
                fi = fm.solve_images(fb[0], fb[1], half=3 * r['truth']['theta_E'], n=301)
                cnt.append(len(fi))
                if len(ti) == len(fi) == 4:
                    dt_t = np.array(L.time_delays(tm, ti, tb[0], tb[1], 0.5, 2.0))
                    # match by position
                    fi_m = []
                    for a in ti:
                        k = int(np.argmin([math.hypot(a['x'] - b['x'], a['y'] - b['y']) for b in fi]))
                        fi_m.append(fi[k])
                    dt_f = np.array(L.time_delays(fm, fi_m, fb[0], fb[1], 0.5, 2.0))
                    span = dt_t.max() - dt_t.min()
                    dd.append(float(np.max(np.abs(dt_f - dt_t)) / max(span, 1e-9)))
            R['predicted_images_count_quads'] = dict(counts={str(k): cnt.count(k) for k in sorted(set(cnt))})
            R['time_delay_max_rel_error_over_span'] = stats(dd)
    # 3. misspecification
    rng = np.random.default_rng(77)
    mis = dict(with_perturber_known=[], without_perturber=[])
    for i in range(min(n, 20)):
        lens = S.random_lens(rng, quad=True)
        p = lens['params']
        # perturber truth placed near a random image
        im = lens['images'][int(rng.integers(0, 4))]
        px, py = im['x'] + rng.uniform(-0.8, 0.8), im['y'] + rng.uniform(-0.8, 0.8)
        pert = L.PIEMD(b0=0.12, a=0.05, s=10.0, x0=px, y0=py)
        m = L.build_sie_shear(p, [pert])
        ims = m.solve_images(*lens['source'], half=3 * p['theta_E'], n=301)
        ims = [a for a in ims if abs(a['mu']) > 0.5]
        if len(ims) < 4:
            continue
        obs = np.array([[a['x'], a['y']] for a in ims[:4]]) + rng.normal(0, 0.003, (4, 2))
        for key, pert_list in (('with_perturber_known', [pert]), ('without_perturber', [])):
            f = L.fit_sie_shear(obs, 0.003, (0.0, 0.0), perturbers=pert_list, seed=i)
            if f.get('success'):
                mis[key].append(dict(dtheta_E_rel=(f['params']['theta_E'] - p['theta_E']) / p['theta_E'], rms=f['rms_arcsec'], chi2=f['chi2'], dq=f['params']['q'] - p['q']))
    R['perturber_test'] = {k: dict(n=len(v), rms_arcsec=stats([a['rms'] for a in v]), chi2=stats([a['chi2'] for a in v]), dtheta_E_rel=stats([a['dtheta_E_rel'] for a in v]),
                                   dq=stats([a['dq'] for a in v])) for k, v in mis.items()}
    # 5. source-plane reconstruction with the true model and with a 3 % / 10 % wrong Einstein radius
    rec = {}
    for err in (0.0, 0.03, 0.10):
        rr = [S.source_recon_test(seed=sd, theta_E_error=err) for sd in range(1, 9)]
        rec['theta_E_error_%g' % err] = dict(corr=stats([r['corr'] for r in rr]), flux_ratio=stats([r['flux_ratio'] for r in rr]),
                                             resid_over_noise=stats([r['resid_over_noise'] for r in rr]), centroid_err_arcsec=stats([r['centroid_err'] for r in rr]))
    R['source_reconstruction'] = rec
    # 6. optional independent cross-check against lenstronomy (only when installed)
    try:
        import lenstronomy
        rg = np.random.default_rng(0)
        mx = 0.0
        for _ in range(50):
            p = dict(theta_E=rg.uniform(0.8, 2.5), q=rg.uniform(0.5, 0.98), phi=rg.uniform(0, 180), x0=rg.uniform(-.1, .1), y0=rg.uniform(-.1, .1))
            x, y = rg.uniform(-4, 4, 30), rg.uniform(-4, 4, 30)
            ax, ay = L.SIE(**p).alpha(x, y)
            bx, by = L.lenstronomy_alpha('SIE', p, x, y)
            mx = max(mx, float(np.abs(ax - bx).max()), float(np.abs(ay - by).max()))
        R['lenstronomy_cross_check'] = dict(version=lenstronomy.__version__, sie_max_abs_alpha_diff_arcsec=mx, n_lenses=50, points_per_lens=30)
    except ImportError:
        R['lenstronomy_cross_check'] = dict(skipped='lenstronomy not installed')
    R['elapsed_s'] = time.time() - t0
    json.dump(R, open(out, 'w'), indent=1)
    print(json.dumps(R, indent=1)[:3000])


if __name__ == '__main__':
    main()
