#!/usr/bin/env python3
"""Completeness and limiting magnitude by artificial-source injection.

Artificial stars (PSF stamps) or galaxies (Sersic, PSF-convolved) are injected at random positions in magnitude (and size) bins, a
DETECTION CALLBACK is run on the image, and the recovery fraction vs magnitude (and size), the photometric bias and scatter of the recovered
sources and the false-positive rate are measured.  The 50 % and 90 % completeness limits come from a fitted logistic curve (and from direct
interpolation of the data).  The result is a JSON-serialisable dict; the CLI step writes it as JSON + a curve TSV, adds catalog columns and
updates the catalog metadata (`catalog_meta.json`, see ogfkit/meta.py).

Detection callback API (any detection step can be tested this way)::

    detect(image: float32 2-D ndarray, mask: bool ndarray | None) -> dict-of-arrays | list-of-dicts
        required keys x, y (0-based pixel centres); optional mag (calibrated), a, b.   Must be picklable for n_workers > 1.

Built in: sextract_detector (the application's own ds9_sextract binary = "the existing extraction"), sep_detector, and "module:function"
(a factory called with the keyword options of the CLI; see load_detector).

CLI:  completeness.py IMAGE [--detector sextract|sep|pkg.mod:fn] [--kind star|galaxy] [--mag-min 22 --mag-max 29 --n-bins 14 --per-bin 60]
        [--per-image 25] [--mag-zeropoint 25] [--psf FILE | --psf-fwhm 3] [--match-radius 3] [--re-bins 2,4,8,16] [--catalog CAT]
        [--mask MASK] [--crop-size N] [--plot-out F] [--json-out F] [--curve-out F] [--meta-out catalog_meta.json] [--n-workers N] [--seed 1] [detector options]
        With --catalog the TSV `NUMBER COMPL_FRAC COMPL_LIM50 COMPL_LIM90` is printed (completeness at each object's MAG_AUTO).
"""
import argparse
import importlib
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.abspath(os.path.join(_here, '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

from ogfkit import tsvio, imageio, models, meta  # noqa: E402

COLUMNS = ['COMPL_FRAC', 'COMPL_LIM50', 'COMPL_LIM90']


# ------------------------------------------------------------------------------------------------ detectors
def _norm(det):
    """Normalise a detector's return value to dict of float arrays x, y, mag (+ a, b)."""
    if isinstance(det, dict):
        d = {k: np.asarray(v, dtype=float) for k, v in det.items()}
    else:
        det = list(det or [])
        keys = set().union(*[set(r) for r in det]) if det else {'x', 'y'}
        d = {k: np.array([float(r.get(k, np.nan)) for r in det], dtype=float) for k in keys | {'x', 'y'}}
    n = len(d.get('x', []))
    for k in ('x', 'y'):
        d.setdefault(k, np.zeros(0))
    for k in ('mag', 'a', 'b'):
        if k not in d or len(d[k]) != n:
            d[k] = np.full(n, np.nan)
    return d


class SextractDetector:
    """Runs bin/ds9_sextract on the (temporary) image and parses X_IMAGE / Y_IMAGE / MAG_AUTO / A_IMAGE / B_IMAGE."""

    def __init__(self, zp=25.0, thresh=1.5, minarea=5, binary=None, extra=()):
        self.zp, self.thresh, self.minarea, self.extra = zp, thresh, minarea, list(extra)
        self.binary = binary or os.path.join(_root, 'bin', 'ds9_sextract')

    def __call__(self, img, mask=None):
        with tempfile.TemporaryDirectory(prefix='ogf_cmp_') as td:
            fn = os.path.join(td, 'im.fits')
            a = np.asarray(img, dtype=np.float32)
            if mask is not None:
                a = a.copy()
                a[mask] = np.nan
            imageio.save_fits(fn, a)
            p = subprocess.run([self.binary, fn, '--detect-thresh', str(self.thresh), '--detect-minarea', str(self.minarea),
                                '--mag-zeropoint', str(self.zp)] + self.extra, capture_output=True, text=True)
            if p.returncode != 0:
                raise RuntimeError('ds9_sextract failed: ' + p.stderr[-200:])
        lines = [l for l in p.stdout.split('\n') if l.strip() and not l.startswith('#')]
        if len(lines) < 2:
            return dict(x=[], y=[], mag=[], a=[], b=[])
        cols = lines[0].split('\t')
        ix = {c: cols.index(c) for c in ('X_IMAGE', 'Y_IMAGE', 'MAG_AUTO', 'A_IMAGE', 'B_IMAGE')}
        rows = [l.split('\t') for l in lines[1:]]
        out = {}
        for key, c in (('x', 'X_IMAGE'), ('y', 'Y_IMAGE'), ('mag', 'MAG_AUTO'), ('a', 'A_IMAGE'), ('b', 'B_IMAGE')):
            out[key] = np.array([tsvio.fnum(r[ix[c]]) for r in rows])
        out['x'] -= 1.0
        out['y'] -= 1.0                      # to 0-based
        out['mag'] = np.where(out['mag'] > 90, np.nan, out['mag'])
        return out


class SepDetector:
    """SEP (Source Extractor as a library): background + extract + elliptical Kron-like aperture photometry."""

    def __init__(self, zp=25.0, thresh=1.5, minarea=5, back_size=64, deblend_cont=0.005, filter=True, local_rms=False):
        self.zp, self.thresh, self.minarea, self.bs, self.dc = zp, thresh, minarea, back_size, deblend_cont
        self.filter = filter
        self.local_rms = local_rms          # threshold relative to the local rms map instead of the global rms (varying depth)

    def __call__(self, img, mask=None):
        import sep
        a = np.ascontiguousarray(img, dtype=np.float32)
        m = None if mask is None else np.ascontiguousarray(mask.astype(np.uint8))
        bad = ~np.isfinite(a)
        if bad.any():
            a = a.copy()
            a[bad] = 0
            m = bad.astype(np.uint8) if m is None else (m | bad.astype(np.uint8))
        bkg = sep.Background(a, mask=m, bw=self.bs, bh=self.bs)
        sub = a - bkg.back()
        o = sep.extract(sub, self.thresh, err=(bkg.rms() if self.local_rms else bkg.globalrms), minarea=self.minarea, mask=m, deblend_cont=self.dc,
                        **({} if self.filter else {'filter_kernel': None}))
        if len(o):
            o = o[np.isfinite(o['a']) & np.isfinite(o['b']) & np.isfinite(o['theta']) & np.isfinite(o['x']) & np.isfinite(o['y'])]
        if len(o) == 0:
            return dict(x=[], y=[], mag=[], a=[], b=[])
        oa = np.maximum(o['a'], 0.8)
        ob = np.minimum(np.maximum(o['b'], 0.8), oa)             # single-pixel detections have degenerate ellipses
        th = np.clip(o['theta'], -np.pi / 2, np.pi / 2)
        kr, _ = sep.kron_radius(sub, o['x'], o['y'], oa, ob, th, 6.0, mask=m)
        kr = np.where(kr > 0, kr, 1.0)
        fl, _, _ = sep.sum_ellipse(sub, o['x'], o['y'], oa, ob, th, 2.5 * kr, mask=m, subpix=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            mag = np.where(fl > 0, self.zp - 2.5 * np.log10(fl), np.nan)
        return dict(x=o['x'], y=o['y'], mag=mag, a=o['a'], b=o['b'])


def load_detector(spec, **kw):
    """'sextract' | 'sep' | 'package.module:factory' (factory(**kw) -> callback; unknown kw are passed through)."""
    if spec in ('sextract', '', None):
        return SextractDetector(**{k: v for k, v in kw.items() if k in ('zp', 'thresh', 'minarea', 'binary', 'extra')})
    if spec == 'sep':
        return SepDetector(**{k: v for k, v in kw.items() if k in ('zp', 'thresh', 'minarea', 'back_size', 'deblend_cont', 'filter', 'local_rms')})
    mod, _, fn = spec.partition(':')
    if not fn:
        raise ValueError('detector spec must be sextract, sep or module:function')
    return getattr(importlib.import_module(mod), fn)(**kw)


# ------------------------------------------------------------------------------------------------ injection
def default_psf(fwhm, size=None):
    return models.gaussian_psf(size or models.stamp_size_for(fwhm), fwhm)


def _draw_objects(rng, n, shape, mags, kind, exclude, border, min_sep, re_edges, n_sersic, q_range, existing_xy=None, existing_r=None):
    """Random positions (not in `exclude`, >= min_sep from each other and from existing sources) + parameters."""
    ny, nx = shape
    pts = []
    tries = 0
    while len(pts) < n and tries < 400 * n:
        tries += 1
        x = rng.uniform(border, nx - 1 - border)
        y = rng.uniform(border, ny - 1 - border)
        if exclude is not None and exclude[int(round(y)), int(round(x))]:
            continue
        if pts and min((x - p[0]) ** 2 + (y - p[1]) ** 2 for p in pts) < min_sep ** 2:
            continue
        pts.append((x, y))
    objs = []
    for (x, y), m in zip(pts, mags):
        o = dict(x=float(x), y=float(y), mag=float(m))
        if kind == 'galaxy':
            k = int(rng.integers(len(re_edges) - 1))
            o['re'] = float(np.exp(rng.uniform(np.log(re_edges[k]), np.log(re_edges[k + 1]))))
            o['size_bin'] = k
            o['n'] = float(rng.choice(n_sersic))
            o['q'] = float(rng.uniform(*q_range))
            o['pa'] = float(rng.uniform(0, 180))
        objs.append(o)
    return objs


def render_object(o, psf, zp, gain=0.0, rng=None):
    """Stamp (2-D array) and its lower-left corner (y0, x0) for one injected object (flux from the magnitude)."""
    flux = 10.0 ** (-0.4 * (o['mag'] - zp))
    h = psf.shape[0] // 2
    if 're' in o:
        half = int(min(np.ceil(7.0 * o['re'] * (1.0 + 0.5 * (o['n'] > 2)) + h), 160))
        half = max(half, h + 4)
        size = 2 * half + 1
        ix, iy = int(round(o['x'])), int(round(o['y']))
        Ie = models.sersic_Ie_from_flux(flux, o['re'], o['n'], o['q'])
        gal = models.render_sersic((size, size), o['x'] - (ix - half), o['y'] - (iy - half), Ie, o['re'], o['n'], o['q'], o['pa'])
        stamp = models.convolve_same(gal, psf.stamp(o['x'], o['y'], o['x'] - ix, o['y'] - iy) if hasattr(psf, 'stamp') else psf)
        y0, x0 = iy - half, ix - half
    else:
        ix, iy = int(round(o['x'])), int(round(o['y']))
        stamp = flux * (psf.stamp(o['x'], o['y'], o['x'] - ix, o['y'] - iy) if hasattr(psf, 'stamp')       # spatially varying PSFModel (ogfkit.psfmodel)
                        else models.shifted_psf(psf, o['x'] - ix, o['y'] - iy))
        y0, x0 = iy - h, ix - h
    if gain and rng is not None:
        stamp = stamp + rng.normal(0.0, 1.0, stamp.shape) * np.sqrt(np.clip(stamp, 0, None) / gain)
    return stamp, y0, x0


def inject(data, objs, psf, zp, gain=0.0, rng=None):
    out = np.array(data, dtype=np.float32, copy=True)
    ny, nx = out.shape
    for o in objs:
        st, y0, x0 = render_object(o, psf, zp, gain, rng)
        ya, xa = max(y0, 0), max(x0, 0)
        yb, xb = min(y0 + st.shape[0], ny), min(x0 + st.shape[1], nx)
        if ya < yb and xa < xb:
            out[ya:yb, xa:xb] += st[ya - y0:yb - y0, xa - x0:xb - x0].astype(np.float32)
    return out


def match_nearest(ox, oy, dx, dy, radius):
    """Greedy one-to-one nearest-neighbour match of objects (ox, oy) to detections (dx, dy) within radius.
    Returns idx array (-1 = no match)."""
    idx = -np.ones(len(ox), int)
    if len(dx) == 0:
        return idx
    D = np.hypot(ox[:, None] - dx[None, :], oy[:, None] - dy[None, :])
    pairs = np.argwhere(D <= radius)
    order = np.argsort(D[pairs[:, 0], pairs[:, 1]])
    used_d = set()
    for k in order:
        i, j = pairs[k]
        if idx[i] < 0 and j not in used_d:
            idx[i] = j
            used_d.add(j)
    return idx


def _realisation(args):
    (data, mask, detect, objs, psf, zp, gain, seed, base_xy, match_radius) = args
    rng = np.random.default_rng(seed)
    img = inject(data, objs, psf, zp, gain, rng)
    det = _norm(detect(img, mask))
    ox = np.array([o['x'] for o in objs]); oy = np.array([o['y'] for o in objs])
    rad = np.array([max(match_radius, 0.5 * o.get('re', 0.0)) for o in objs])
    # match with a per-object radius: use the largest and filter afterwards
    idx = match_nearest(ox, oy, det['x'], det['y'], float(rad.max()))
    res = []
    used = set()
    for i, o in enumerate(objs):
        j = idx[i]
        ok = j >= 0 and math.hypot(det['x'][j] - o['x'], det['y'][j] - o['y']) <= rad[i]
        r = dict(o)
        r['recovered'] = bool(ok)
        r['mag_out'] = float(det['mag'][j]) if ok else float('nan')
        r['dpos'] = float(math.hypot(det['x'][j] - o['x'], det['y'][j] - o['y'])) if ok else float('nan')
        if ok:
            used.add(int(j))
        res.append(r)
    # detections that are neither a baseline source nor an injected object (spurious / fragments)
    n_unmatched = 0
    for j in range(len(det['x'])):
        if j in used:
            continue
        if base_xy is not None and len(base_xy[0]) and np.min(np.hypot(base_xy[0] - det['x'][j], base_xy[1] - det['y'][j])) <= 2 * match_radius:
            continue
        if np.min(np.hypot(ox - det['x'][j], oy - det['y'][j])) <= 3 * match_radius + 0.5 * max([o.get('re', 0.0) for o in objs] + [0]):
            continue
        n_unmatched += 1
    return res, n_unmatched


# ------------------------------------------------------------------------------------------------ statistics
def wilson(k, n, z=1.0):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def logistic(m, m50, w, fmax=1.0):
    return fmax / (1.0 + np.exp((np.asarray(m, dtype=float) - m50) / max(w, 1e-3)))


def fit_curve(mag, k, n):
    """Weighted binomial least-squares fit of f(m) = fmax / (1 + exp((m - m50)/w)); returns dict(m50, w, fmax, m90, ok)."""
    from scipy.optimize import least_squares
    mag = np.asarray(mag, float); k = np.asarray(k, float); n = np.asarray(n, float)
    use = n > 0
    if use.sum() < 3 or k[use].sum() == 0:
        return dict(ok=False, mh=float('nan'), m50=float('nan'), w=float('nan'), fmax=float('nan'), m90=float('nan'))
    f = k / np.where(n > 0, n, 1)
    m0 = float(np.interp(0.5, f[use][::-1], mag[use][::-1])) if f[use].min() < 0.5 < f[use].max() else float(np.median(mag[use]))

    def resid(p):
        mod = logistic(mag[use], p[0], p[1], p[2])
        sig = np.sqrt(np.clip(mod * (1 - mod), 0.01, None) / n[use])
        return (f[use] - mod) / sig

    sol = least_squares(resid, [m0, 0.4, min(1.0, max(f[use].max(), 0.6))], bounds=([mag.min() - 5, 0.02, 0.3], [mag.max() + 5, 5.0, 1.0]))
    mh, w, fmax = [float(v) for v in sol.x]
    # f = fmax / (1 + exp((m - mh)/w));  f = L  ->  m = mh + w ln(fmax/L - 1)
    lim = lambda L: float(mh + w * math.log(fmax / L - 1.0)) if fmax > L else float('nan')
    return dict(ok=bool(sol.success), mh=mh, w=w, fmax=fmax, m50=lim(0.5), m90=lim(0.9), chi2=float(np.sum(sol.fun ** 2)),
                ndof=int(use.sum() - 3))


def direct_limit(mag, f, level):
    """Faintest-going crossing of `level` by linear interpolation of the binned data (bright -> faint); NaN if never crossed."""
    mag = np.asarray(mag, float); f = np.asarray(f, float)
    ok = np.isfinite(f)
    mag, f = mag[ok], f[ok]
    for i in range(len(mag) - 1):
        if f[i] >= level > f[i + 1]:
            return float(mag[i] + (f[i] - level) / (f[i] - f[i + 1]) * (mag[i + 1] - mag[i]))
    return float('nan')


def summarize(records, edges, n_fp_neg, area_pix, n_unmatched, n_images, pixel_scale=None):
    mags = np.array([r['mag'] for r in records])
    rec = np.array([r['recovered'] for r in records])
    dm = np.array([r['mag_out'] - r['mag'] for r in records])
    centers = 0.5 * (edges[:-1] + edges[1:])
    bins = []
    for i in range(len(centers)):
        sel = (mags >= edges[i]) & (mags < edges[i + 1] if i < len(centers) - 1 else mags <= edges[i + 1])
        n = int(sel.sum()); k = int((sel & rec).sum())
        lo, hi = wilson(k, n)
        d = dm[sel & rec & np.isfinite(dm)]
        bins.append(dict(mag=float(centers[i]), mag_lo=float(edges[i]), mag_hi=float(edges[i + 1]), n_inj=n, n_rec=k,
                         frac=(k / n if n else float('nan')), frac_lo=lo, frac_hi=hi,
                         bias=(float(np.median(d)) if len(d) else float('nan')),
                         scatter=(float(1.4826 * np.median(np.abs(d - np.median(d)))) if len(d) > 2 else float('nan')),
                         n_phot=int(len(d))))
    k = np.array([b['n_rec'] for b in bins]); n = np.array([b['n_inj'] for b in bins])
    c = np.array([b['mag'] for b in bins]); f = np.array([b['frac'] for b in bins])
    fit = fit_curve(c, k, n)
    out = dict(bins=bins, fit=fit, direct=dict(m50=direct_limit(c, f, 0.5), m90=direct_limit(c, f, 0.9)),
               lim50=fit['m50'] if np.isfinite(fit['m50']) else direct_limit(c, f, 0.5),
               lim90=fit['m90'] if np.isfinite(fit['m90']) else direct_limit(c, f, 0.9),
               n_injected=int(n.sum()), n_recovered=int(k.sum()), n_images=int(n_images),
               false_positive=dict(negative_image_detections=int(n_fp_neg), area_pixels=int(area_pix),
                                   per_1e6_pix=(1e6 * n_fp_neg / area_pix if area_pix else float('nan')),
                                   per_arcmin2=(n_fp_neg / (area_pix * pixel_scale ** 2 / 3600.0) if (area_pix and pixel_scale) else None),
                                   unmatched_per_image=(n_unmatched / max(n_images, 1))))
    return out


# ------------------------------------------------------------------------------------------------ driver
def run_completeness(data, detect, mask=None, kind='star', mag_min=22.0, mag_max=29.0, n_bins=14, per_bin=60, per_image=25,
                     zp=25.0, psf=None, psf_fwhm=3.0, match_radius=3.0, re_edges=(2.0, 4.0, 8.0, 16.0), n_sersic=(1.0, 4.0),
                     q_range=(0.4, 1.0), seed=1, gain=0.0, pixel_scale=None, avoid_detected=True, min_sep=None, n_workers=0,
                     progress=None, return_records=False):
    """Full completeness measurement.  `detect` is the detection callback (see module docstring).  Returns a JSON-serialisable dict:
    config, bins (per magnitude bin: n_inj, n_rec, frac + Wilson interval, bias, scatter), by_size (galaxies), fit, lim50, lim90, false_positive."""
    data = np.asarray(data, dtype=np.float32)
    ny, nx = data.shape
    psf = default_psf(psf_fwhm) if psf is None else psf
    rng = np.random.default_rng(seed)
    edges = np.linspace(mag_min, mag_max, n_bins + 1)
    base = _norm(detect(data, mask))
    base_xy = (base['x'], base['y'])
    # exclusion map: masked pixels + discs around baseline detections (so that injections land on empty sky)
    exclude = np.zeros(data.shape, bool) if mask is None else mask.copy()
    exclude |= ~np.isfinite(data)
    if avoid_detected and len(base['x']):
        yy, xx = np.mgrid[:ny, :nx]
        for x, y, a in zip(base['x'], base['y'], np.nan_to_num(base['a'], nan=2.0)):
            r = max(2.0 * match_radius, 3.0 * a) + (0.0 if kind == 'star' else 0.5 * max(re_edges))
            xa, xb, ya, yb = int(max(0, x - r)), int(min(nx, x + r + 1)), int(max(0, y - r)), int(min(ny, y + r + 1))
            if xa < xb and ya < yb:
                exclude[ya:yb, xa:xb] |= ((xx[ya:yb, xa:xb] - x) ** 2 + (yy[ya:yb, xa:xb] - y) ** 2) <= r * r
    h = psf.shape[0] // 2
    border = h + (0 if kind == 'star' else int(min(3 * max(re_edges), 60)))
    border = max(border, 6)
    n_total = n_bins * per_bin
    mags = np.concatenate([rng.uniform(edges[i], edges[i + 1], per_bin) for i in range(n_bins)])
    rng.shuffle(mags)
    if min_sep is None:
        min_sep = max(8.0 * match_radius, 2.0 * h, 6.0 * max(re_edges) if kind == 'galaxy' else 0.0)
    jobs, k = [], 0
    nimg = int(math.ceil(n_total / per_image))
    for i in range(nimg):
        m = mags[i * per_image:(i + 1) * per_image]
        objs = _draw_objects(rng, len(m), data.shape, m, kind, exclude, border, min_sep, list(re_edges), list(n_sersic), q_range)
        jobs.append((data, mask, detect, objs, psf, zp, gain, int(rng.integers(1 << 30)), base_xy, match_radius))
    if n_workers and n_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(n_workers) as ex:
            outs = list(ex.map(_realisation, jobs))
    else:
        outs = []
        for i, j in enumerate(jobs):
            outs.append(_realisation(j))
            if progress:
                progress(i + 1, len(jobs))
    records = [r for o in outs for r in o[0]]
    n_unmatched = sum(o[1] for o in outs)
    # false positives on the negative image: mirror around the median background
    good = np.isfinite(data) & (~mask if mask is not None else True)
    med = float(np.median(data[good])) if good.any() else 0.0
    neg = _norm(detect((2.0 * med - data).astype(np.float32), mask))
    area = int(good.sum())
    res = summarize(records, edges, len(neg['x']), area, n_unmatched, len(jobs), pixel_scale)
    if return_records:
        res['records'] = records          # every injected object (x, y, mag, recovered, mag_out, ...); not JSON-cleaned
    res['config'] = dict(kind=kind, mag_min=mag_min, mag_max=mag_max, n_bins=n_bins, per_bin=per_bin, per_image=per_image, zp=zp,
                         match_radius=match_radius, seed=seed, avoid_detected=bool(avoid_detected), n_baseline_detections=int(len(base['x'])),
                         min_sep=float(min_sep), psf_size=int(psf.shape[0]), detector=type(detect).__name__,
                         re_edges=list(re_edges) if kind == 'galaxy' else None, n_sersic=list(n_sersic) if kind == 'galaxy' else None)
    if kind == 'galaxy':
        by = []
        sb = np.array([r['size_bin'] for r in records])
        for s in range(len(re_edges) - 1):
            rr = [r for r, b in zip(records, sb) if b == s]
            if not rr:
                continue
            sub = summarize(rr, edges, 0, area, 0, len(jobs), pixel_scale)
            by.append(dict(re_lo=float(re_edges[s]), re_hi=float(re_edges[s + 1]), lim50=sub['lim50'], lim90=sub['lim90'],
                           frac=[b['frac'] for b in sub['bins']], n_inj=[b['n_inj'] for b in sub['bins']]))
        res['by_size'] = by
    return res


def json_clean(o):
    """NaN / inf -> None recursively (strict JSON for browsers and other readers)."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: json_clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return json_clean(o.item())
    return o


def write_curve(path, res):
    cols = ['mag', 'n_inj', 'n_rec', 'frac', 'frac_lo', 'frac_hi', 'bias', 'scatter', 'n_phot']
    tsvio.write_table(path, cols, res['bins'])


def plot_curve(res, path):
    """Completeness curve (data with Wilson 68 % intervals, logistic fit, 50/90 % limits), bias and scatter; PNG (Agg)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b = res['bins']
    m = np.array([x['mag'] for x in b]); f = np.array([x['frac'] for x in b])
    lo = np.array([x['frac_lo'] for x in b]); hi = np.array([x['frac_hi'] for x in b])
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].errorbar(m, f, [f - lo, hi - f], fmt='ko', ms=4, label='recovered fraction')
    mm = np.linspace(m.min() - 0.5, m.max() + 0.5, 200)
    ax[0].plot(mm, curve_at(res, mm), 'r-', label='logistic fit')
    for lvl, key in ((0.5, 'lim50'), (0.9, 'lim90')):
        if np.isfinite(res[key]):
            ax[0].axvline(res[key], ls='--', c='gray'); ax[0].axhline(lvl, ls=':', c='gray')
            ax[0].text(res[key], lvl + 0.02, ' %d%%: %.2f' % (lvl * 100, res[key]), fontsize=8)
    ax[0].set_xlabel('input magnitude'); ax[0].set_ylabel('completeness'); ax[0].set_ylim(-0.02, 1.05); ax[0].legend(fontsize=8)
    ax[0].set_title('%s, %d injected' % (res['config']['kind'], res['n_injected']), fontsize=9)
    ax[1].errorbar(m, [x['bias'] for x in b], [x['scatter'] for x in b], fmt='bs-', ms=4)
    ax[1].axhline(0, c='gray'); ax[1].set_xlabel('input magnitude'); ax[1].set_ylabel('measured - input [mag] (median, robust scatter)')
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)
    return path


def curve_at(res, mags):
    """Completeness (0..1) at the magnitudes `mags`: the fitted logistic if the fit worked, else linear interpolation of the binned fraction."""
    mags = np.asarray(mags, float)
    fit = res['fit']
    if fit.get('ok') and np.isfinite(fit.get('mh', np.nan)):
        return np.clip(logistic(mags, fit['mh'], fit['w'], fit['fmax']), 0, 1)
    c = np.array([b['mag'] for b in res['bins']]); f = np.array([b['frac'] for b in res['bins']])
    ok = np.isfinite(f)
    return np.interp(mags, c[ok], f[ok], left=f[ok][0], right=f[ok][-1]) if ok.any() else np.full(len(mags), np.nan)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--detector', default='sextract')
    ap.add_argument('--kind', default='star', choices=['star', 'galaxy'])
    ap.add_argument('--mag-min', type=float, default=22.0)
    ap.add_argument('--mag-max', type=float, default=29.0)
    ap.add_argument('--n-bins', type=int, default=14)
    ap.add_argument('--per-bin', type=int, default=60)
    ap.add_argument('--per-image', type=int, default=25)
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--psf', default='')
    ap.add_argument('--psf-fwhm', type=float, default=3.0)
    ap.add_argument('--match-radius', type=float, default=3.0)
    ap.add_argument('--re-bins', default='2,4,8,16')
    ap.add_argument('--n-sersic', default='1,4')
    ap.add_argument('--gain', type=float, default=0.0)
    ap.add_argument('--pixel-scale', type=float, default=0.0)
    ap.add_argument('--detect-thresh', type=float, default=1.5)
    ap.add_argument('--detect-minarea', type=int, default=5)
    ap.add_argument('--allow-blends', action='store_true', help='inject on top of real sources too (completeness includes blending)')
    ap.add_argument('--catalog', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--crop', default='', help='x0,y0,x1,y1 (0-based pixels) sub-image to test (large images)')
    ap.add_argument('--crop-size', type=int, default=0, help='use the central N x N pixels only (0 = whole image)')
    ap.add_argument('--plot-out', default='')
    ap.add_argument('--json-out', default='')
    ap.add_argument('--curve-out', default='')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--n-workers', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args(argv)

    data, hdr = imageio.load_image(a.image)
    mask = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    if a.crop:
        x0, y0, x1, y1 = [int(v) for v in a.crop.split(',')]
        data = np.ascontiguousarray(data[y0:y1, x0:x1])
        mask = None if mask is None else np.ascontiguousarray(mask[y0:y1, x0:x1])
    if a.crop_size and min(data.shape) > a.crop_size:
        cy, cx = data.shape[0] // 2, data.shape[1] // 2
        h = a.crop_size // 2
        data = np.ascontiguousarray(data[cy - h:cy - h + a.crop_size, cx - h:cx - h + a.crop_size])
        mask = None if mask is None else np.ascontiguousarray(mask[cy - h:cy - h + a.crop_size, cx - h:cx - h + a.crop_size])
    psf = None
    if a.psf and os.path.isfile(a.psf):
        psf, _ = imageio.load_image(a.psf)
        psf = psf / psf.sum()
    det = load_detector(a.detector, zp=a.mag_zeropoint, thresh=a.detect_thresh, minarea=a.detect_minarea)
    res = run_completeness(data, det, mask=mask, kind=a.kind, mag_min=a.mag_min, mag_max=a.mag_max, n_bins=a.n_bins, per_bin=a.per_bin,
                           per_image=a.per_image, zp=a.mag_zeropoint, psf=psf, psf_fwhm=a.psf_fwhm, match_radius=a.match_radius,
                           re_edges=[float(v) for v in a.re_bins.split(',')], n_sersic=[float(v) for v in a.n_sersic.split(',')],
                           seed=a.seed, gain=a.gain, pixel_scale=(a.pixel_scale or None), avoid_detected=not a.allow_blends,
                           n_workers=a.n_workers)
    res['image'] = os.path.basename(a.image)
    if a.json_out:
        with open(a.json_out, 'w') as f:
            json.dump(json_clean(res), f, indent=1, allow_nan=False)
    if a.curve_out:
        write_curve(a.curve_out, res)
    if a.plot_out:
        plot_curve(res, a.plot_out)
    sys.stderr.write('completeness: %s, %d injected, %d recovered; 50%% limit %.2f, 90%% limit %.2f (fit), %s detections on the negative image\n'
                     % (a.kind, res['n_injected'], res['n_recovered'], res['lim50'], res['lim90'], res['false_positive']['negative_image_detections']))
    if not a.catalog:
        print('completeness (%s): 50%% limit %.2f, 90%% limit %.2f' % (a.kind, res['lim50'], res['lim90']))
        print('mag\tn_inj\tn_rec\tfrac\tbias\tscatter')
        for b in res['bins']:
            print('%.2f\t%d\t%d\t%.3f\t%s\t%s' % (b['mag'], b['n_inj'], b['n_rec'], b['frac'], tsvio.fmt(b['bias'], 3), tsvio.fmt(b['scatter'], 3)))
    if a.catalog:
        cols, rows = tsvio.read_catalog(a.catalog)
        mags = np.array([tsvio.fnum(r.get('MAG_AUTO')) for r in rows])
        comp = curve_at(res, np.where(np.isfinite(mags), mags, 99.0))
        out = []
        for r, c, m in zip(rows, comp, mags):
            out.append((r['NUMBER'], dict(COMPL_FRAC=(float(c) if np.isfinite(m) else None), COMPL_LIM50=res['lim50'], COMPL_LIM90=res['lim90'])))
        tsvio.write_columns(COLUMNS, out)
        if a.meta_out:
            meta.update(a.meta_out, 'completeness', json_clean(dict(kind=a.kind, lim50=res['lim50'], lim90=res['lim90'], fit_m50=res['fit'].get('m50'),
                        fit_width=res['fit'].get('w'), fmax=res['fit'].get('fmax'), n_injected=res['n_injected'], detector=a.detector,
                        mag_zeropoint=a.mag_zeropoint, curve_mag=[b['mag'] for b in res['bins']], curve_frac=[b['frac'] for b in res['bins']])),
                        nrows=len(rows))
    return 0


if __name__ == '__main__':
    sys.exit(main())
