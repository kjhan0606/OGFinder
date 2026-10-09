#!/usr/bin/env python3
"""Synthetic validation of the merger / tidal-feature indicators (morphometry.features) through morphometry.extended.measure_object.

    python features_validate.py OUT.json [--n 40] [--workers 8]

Hosts: Sersic galaxies (n 2-4, R_e 9-14 px, q 0.55-1, random PA, total flux 1e5, Gaussian noise sigma 0.5 per pixel).  Injected: shells (2 concentric arcs at 1.3 and 1.9 R_e, 70-120 deg
wide, radial FWHM 0.25 R_e), tidal tails (curved Gaussian ridge leaving the host at ~1 R_e, 3 R_e long, FWHM 0.4 R_e), each with a fraction f of the host flux; controls (smooth hosts, same noise);
mergers: host + companion of flux ratio 1 / 0.5 / 0.25 / 0.1 at 1.0 / 1.6 / 2.5 R_e separation.  Reported: detection rate vs f with the false-positive rate on the controls, asymmetry and Lotz class of the
mergers, pair classification of a mock catalog."""
def _import_ogfmeas():
    """In-tree measurements. Finds ogfmeas from this file so a script does not need PYTHONPATH."""
    import pathlib
    import sys
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas
    return ogfmeas.measurement_library()

import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path[:0] = [ROOT]
from ogfkit import models as M  # noqa: E402
from morphometry import extended as X, features as FT  # noqa: E402

SHAPE = (321, 321)
XC = YC = 160.0
FTOT = 1e5
SIGMA = 0.5
SIGS = (2.0, 8.0, 24.0)                  # noise multipliers: sigma = 1, 4, 12 per pixel


def host(rng):
    n, re_ = rng.uniform(2, 4), rng.uniform(9, 14)
    q, pa = rng.uniform(0.55, 1.0), rng.uniform(0, 180)
    return dict(n=n, re=re_, q=q, pa=pa, x=XC + rng.uniform(-1, 1), y=YC + rng.uniform(-1, 1))


def render_host(h, flux=FTOT):
    return M.render_sersic(SHAPE, h['x'], h['y'], M.sersic_Ie_from_flux(flux, h['re'], h['n'], h['q']), h['re'], h['n'], h['q'], h['pa'])


def ell_r(h, dx=0.0, dy=0.0):
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]].astype(float)
    t = math.radians(h['pa'])
    u = (xx - h['x'] - dx) * math.cos(t) + (yy - h['y'] - dy) * math.sin(t)
    v = -(xx - h['x'] - dx) * math.sin(t) + (yy - h['y'] - dy) * math.cos(t)
    return np.sqrt(u ** 2 + (v / h['q']) ** 2), np.arctan2(v / h['q'], u)


def shells(h, rng, frac, flux=FTOT):
    img = np.zeros(SHAPE)
    for k, r0 in enumerate((1.3, 1.9)):
        r, phi = ell_r(h, rng.uniform(-0.15, 0.15) * h['re'], rng.uniform(-0.15, 0.15) * h['re'])
        mid = rng.uniform(0, 2 * math.pi)
        half = math.radians(rng.uniform(35, 60))
        dphi = np.abs((phi - mid + math.pi) % (2 * math.pi) - math.pi)
        arc = np.exp(-0.5 * ((r - r0 * h['re']) / (0.25 * h['re'] / 2.355)) ** 2) * np.exp(-0.5 * (dphi / (half / 2)) ** 2)
        img += arc / arc.sum() * frac * flux * (0.6 if k else 0.4) * 1.0
    return img


def tail(h, rng, frac, flux=FTOT):
    img = np.zeros(SHAPE)
    ang = rng.uniform(0, 2 * math.pi)
    curv = rng.uniform(-0.5, 0.5)
    s = np.linspace(0.8, 3.8, 400) * h['re']
    px = h['x'] + s * np.cos(ang + curv * (s / h['re'] - 0.8) * 0.4)
    py = h['y'] + s * np.sin(ang + curv * (s / h['re'] - 0.8) * 0.4)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    w = 0.4 * h['re'] / 2.355
    for a, b in zip(px, py):
        if 2 <= a < SHAPE[1] - 2 and 2 <= b < SHAPE[0] - 2:
            x0, y0 = int(a), int(b)
            sl = (slice(max(0, y0 - 14), y0 + 15), slice(max(0, x0 - 14), x0 + 15))
            img[sl] += np.exp(-0.5 * ((xx[sl] - a) ** 2 + (yy[sl] - b) ** 2) / w ** 2)
    return img / img.sum() * frac * flux


def measure(img, h, seed):
    return X.measure_object(img, h['x'], h['y'], h['q'], math.radians(h['pa']), SIGMA, rmax=130, features=True, seed=seed)


def trial(task):
    kind, seed, f, sig = task
    rng = np.random.default_rng(10007 * seed + (hash(kind) % 997 if False else {'control': 1, 'shell': 2, 'tail': 3, 'merger': 4}[kind]))
    h = host(rng)
    img = render_host(h)
    meta = {}
    if kind == 'shell':
        img = img + shells(h, rng, f)
    elif kind == 'tail':
        img = img + tail(h, rng, f)
    elif kind == 'merger':
        ratio, sep = f
        t = rng.uniform(0, 2 * math.pi)
        h2 = dict(n=1.5, re=h['re'] * math.sqrt(ratio) * 0.8 + 2, q=0.9, pa=0.0, x=h['x'] + sep * h['re'] * math.cos(t), y=h['y'] + sep * h['re'] * math.sin(t))
        img = img + render_host(h2, FTOT * ratio)
    noise_rng = np.random.default_rng(seed + 99)
    img = img + noise_rng.normal(0, SIGMA * sig, SHAPE)
    global SIGMA_USED
    r = X.measure_object(img, h['x'], h['y'], h['q'], math.radians(h['pa']), SIGMA * sig, rmax=130, features=True, seed=seed)
    ft = r['feat']
    return dict(kind=kind, seed=seed, f=f if not isinstance(f, tuple) else list(f), sig=sig, asym=ft['asym'], flags=ft['flags'], lotz=ft['lotz'], n_shell=ft['n_shell'], n_tail=ft['n_tail'], npeak=ft['npeak'],
                shell_frac=ft['shell_frac'], tail_frac=ft['tail_frac'], tidal_frac=ft['tidal_frac'], rp=r['rp'][0.2], re=h['re'])


# ---------------------------------------------------------------------------------------------------------------- injection into a real image
REAL = {}


def _real_setup(img_path, cat_path):
    from astropy.io import fits
    from ogfkit import tsvio
    sep = _import_ogfmeas()
    d = np.ascontiguousarray(np.array(fits.getdata(img_path), float))
    d[~np.isfinite(d)] = 0.0
    bk = sep.Background(d, bw=64, bh=64)
    REAL['img'] = d - bk.back()
    REAL['rms'] = float(np.median(bk.rms()))
    cols, rows = tsvio.read_catalog(cat_path)
    REAL['cat'] = [(float(r['X_IMAGE']) - 1, float(r['Y_IMAGE']) - 1, float(r['A_IMAGE']), float(r['B_IMAGE']), float(r['THETA_IMAGE']), float(r['FLUX_AUTO'])) for r in rows]
    REAL['mask'] = np.zeros(d.shape, bool)
    yy, xx = np.mgrid[:d.shape[0], :d.shape[1]]
    # every catalog object masks a 3 A ellipse (a host is injected on empty sky only, the catalog objects around it are the "neighbours")
    REAL['cat_arr'] = np.array(REAL['cat'])


def real_trial(task):
    kind, seed, f, nsig, img_path, cat_path = task
    if not REAL:
        _real_setup(img_path, cat_path)
    rng = np.random.default_rng(31 * seed + {'control': 1, 'shell': 2, 'tail': 3}[kind])
    d, rms, cat = REAL['img'], REAL['rms'], REAL['cat_arr']
    ny, nx = d.shape
    for _ in range(200):                                       # random position with no catalog object within 60 px (the host is injected on empty sky)
        x, y = rng.uniform(150, nx - 150), rng.uniform(150, ny - 150)
        if np.hypot(cat[:, 0] - x, cat[:, 1] - y).min() > 60:
            break
    h = host(rng)
    h.update(x=x, y=y)
    re_ = h['re']
    Ie = 25.0 * rms                                            # surface brightness at R_e: 25 sigma per pixel
    flux = M.sersic_flux_from_Ie(Ie, re_, h['n'], h['q']) if hasattr(M, 'sersic_flux_from_Ie') else None
    if flux is None:
        flux = float(np.sum(M.render_sersic((401, 401), 200, 200, Ie, re_, h['n'], h['q'], 0.0)))
    global SHAPE, XC, YC
    old = SHAPE
    SHAPE = d.shape
    try:
        m = M.render_sersic(SHAPE, x, y, Ie, re_, h['n'], h['q'], h['pa'])
        if kind == 'shell':
            m = m + shells(h, rng, f, flux)
        elif kind == 'tail':
            m = m + tail(h, rng, f, flux)
    finally:
        SHAPE = old
    x0, y0 = int(x) - 160, int(y) - 160
    cut = (d + m)[y0:y0 + 321, x0:x0 + 321]
    mask = np.zeros(cut.shape, bool)
    yy, xx = np.mgrid[:321, :321]
    for cx, cy, A, B, th, fl in cat:                           # neighbours inside the cutout: 3 A ellipses (as the driver's --mask-neighbours)
        if x0 - 3 * A <= cx <= x0 + 321 + 3 * A and y0 - 3 * A <= cy <= y0 + 321 + 3 * A:
            t = math.radians(th)
            u = (xx - (cx - x0)) * math.cos(t) + (yy - (cy - y0)) * math.sin(t)
            v = -(xx - (cx - x0)) * math.sin(t) + (yy - (cy - y0)) * math.cos(t)
            mask |= (u / (3 * max(A, 1))) ** 2 + (v / (3 * max(B, 1))) ** 2 <= 1.0
    r = X.measure_object(cut, x - x0, y - y0, h['q'], math.radians(h['pa']), rms, rmax=130, mask=mask, features=True, seed=seed, feat_opts=dict(nsig=nsig))
    ft = r.get('feat')
    if ft is None:
        return dict(kind=kind, seed=seed, f=f, nsig=nsig, ok=False)
    return dict(kind=kind, seed=seed, f=f, nsig=nsig, ok=True, n_shell=ft['n_shell'], n_tail=ft['n_tail'], shell_frac=ft['shell_frac'], tail_frac=ft['tail_frac'], asym=ft['asym'], flags=ft['flags'])


def real_main(a):
    tasks = []
    for nsig in (3.0, 4.0, 5.0):
        tasks += [('control', s, 0, nsig, a.real, a.catalog) for s in range(a.n)]
        for f in (0.01, 0.02, 0.05, 0.10):
            tasks += [('shell', s, f, nsig, a.real, a.catalog) for s in range(a.n)] + [('tail', s, f, nsig, a.real, a.catalog) for s in range(a.n)]
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(real_trial, tasks, chunksize=4))
    S = {}
    for nsig in (3.0, 4.0, 5.0):
        rs = [r for r in res if r['nsig'] == nsig and r['ok']]
        ctl = [r for r in rs if r['kind'] == 'control']
        S['nsig_%g' % nsig] = dict(
            n_control=len(ctl), fp_shell=float(np.mean([r['n_shell'] > 0 for r in ctl])), fp_tail=float(np.mean([r['n_tail'] > 0 for r in ctl])),
            shell={str(f): float(np.mean([r['n_shell'] > 0 for r in rs if r['kind'] == 'shell' and r['f'] == f])) for f in (0.01, 0.02, 0.05, 0.10)},
            tail={str(f): float(np.mean([r['n_tail'] > 0 for r in rs if r['kind'] == 'tail' and r['f'] == f])) for f in (0.01, 0.02, 0.05, 0.10)})
    json.dump(dict(summary=S, results=res), open(a.out, 'w'), default=lambda o: None)
    print(json.dumps(S, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--n', type=int, default=40)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--real', default='', help='inject into this real image (background subtracted by sep) instead of the synthetic noise; needs --catalog')
    ap.add_argument('--catalog', default='')
    a = ap.parse_args()
    if a.real:
        return real_main(a)
    tasks = []
    for sig in SIGS:
        tasks += [('control', s, 0, sig) for s in range(a.n)]
        for f in (0.005, 0.01, 0.02, 0.05, 0.10):
            tasks += [('shell', s, f, sig) for s in range(a.n)] + [('tail', s, f, sig) for s in range(a.n)]
    for ratio in (1.0, 0.5, 0.25, 0.1):
        for sep in (1.0, 1.6, 2.5):
            tasks += [('merger', s, (ratio, sep), 2.0) for s in range(max(a.n // 2, 10))]
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(trial, tasks, chunksize=4))
    S = {}
    for sig in SIGS:
        ctl = [r for r in res if r['kind'] == 'control' and r['sig'] == sig]
        S['sigma_%g' % (SIGMA * sig)] = dict(
            control=dict(n=len(ctl), fp_shell=float(np.mean([r['n_shell'] > 0 for r in ctl])), fp_tail=float(np.mean([r['n_tail'] > 0 for r in ctl])), asym_median=float(np.median([r['asym'] for r in ctl])),
                         asym_p95=float(np.percentile([r['asym'] for r in ctl], 95)), merger_flag=float(np.mean([bool(r['flags'] & 3) for r in ctl])), double=float(np.mean([bool(r['flags'] & 4) for r in ctl]))),
            shell={str(f): dict(rate=float(np.mean([r['n_shell'] > 0 for r in res if r['kind'] == 'shell' and r['sig'] == sig and r['f'] == f])),
                                frac_recovered=float(np.nanmedian([r['shell_frac'] / f for r in res if r['kind'] == 'shell' and r['sig'] == sig and r['f'] == f]))) for f in (0.005, 0.01, 0.02, 0.05, 0.10)},
            tail={str(f): dict(rate=float(np.mean([r['n_tail'] > 0 for r in res if r['kind'] == 'tail' and r['sig'] == sig and r['f'] == f])),
                               frac_recovered=float(np.nanmedian([r['tail_frac'] / f for r in res if r['kind'] == 'tail' and r['sig'] == sig and r['f'] == f]))) for f in (0.005, 0.01, 0.02, 0.05, 0.10)})
    mg = {}
    for ratio in (1.0, 0.5, 0.25, 0.1):
        for sep in (1.0, 1.6, 2.5):
            rs = [r for r in res if r['kind'] == 'merger' and r['f'] == [ratio, sep]]
            mg['q%g_sep%g' % (ratio, sep)] = dict(n=len(rs), asym_median=float(np.nanmedian([r['asym'] for r in rs])), merger_flag=float(np.mean([bool(r['flags'] & 3) for r in rs])),
                                                  gm20=float(np.mean([bool(r['flags'] & 1) for r in rs])), asym_flag=float(np.mean([bool(r['flags'] & 2) for r in rs])), double=float(np.mean([bool(r['flags'] & 4) for r in rs])))
    S['mergers'] = mg
    json.dump(dict(summary=S, results=res), open(a.out, 'w'), default=lambda o: None)
    print(json.dumps(S, indent=1)[:6000])


if __name__ == '__main__':
    main()
