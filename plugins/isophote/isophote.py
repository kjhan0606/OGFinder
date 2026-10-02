#!/usr/bin/env python3
"""Isophote ellipse fitting (IRAF `ellipse` style) for the galaxies of the catalog, engine: photutils.isophote.

Per galaxy: radial profile of intensity / surface brightness, ellipticity, position angle, centre, harmonic
deviations A3/B3/A4/B4 (boxy / disky), growth curve (flux inside each isophote), a 2-D model image (cmodel) and the
residual image.  The shared mask (bit-flag mask manager, bool file) and optional masking of the catalog neighbours make the
fit mask-aware.

CLI (catalog step; prints `NUMBER` + new columns as TSV on stdout, all logging on stderr):
  isophote.py IMAGE --catalog CAT.tsv [--numbers 1,5-9 | --max-objects N] [--mask MASK.fits] [--mag-zeropoint 25]
              [--pixel-scale 1] [--maxsma-scale 3] [--step 0.1] [--linear] [--fix-center] [--no-harmonics]
              [--mask-neighbours 2.5] [--n-workers 0] [--model-out M.fits] [--resid-out R.fits]
              [--table-out T.tsv] [--json-out J.json] [--plot-out P.png]
  isophote.py IMAGE --xy X,Y [--sma0 10 --eps0 0.2 --pa0 30] ...        (one object at a 1-based position, no catalog)

Python API (pure functions, JSON-serialisable result): fit_galaxy(), summarize(), growth_radii(), run_catalog().
"""
import argparse
import json
import math
import os
import sys
import warnings

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.abspath(os.path.join(_here, '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

from ogfkit import tsvio, imageio  # noqa: E402

COLUMNS = ['ISO_NISO', 'ISO_RMAX', 'ISO_X0', 'ISO_Y0', 'ISO_R50', 'ISO_R80', 'ISO_MAG_TOT', 'ISO_EPS_HL', 'ISO_PA_HL',
           'ISO_EPS_OUT', 'ISO_PA_OUT', 'ISO_B4_HL', 'ISO_A4_HL', 'ISO_SHAPE', 'ISO_RCONV', 'ISO_FLAG']
PROFILE_COLS = ['NUMBER', 'sma', 'intens', 'intens_err', 'mu', 'mu_err', 'growth_flux', 'ellipticity', 'ellipticity_err', 'pa', 'pa_err',
                'x0', 'y0', 'a3', 'b3', 'a4', 'b4', 'a3_err', 'b3_err', 'a4_err', 'b4_err', 'tflux_e', 'npix_e', 'rms',
                'n_data', 'n_flag', 'stop_code']
DISKY_LIMIT = 0.01      # |B4| above this is called disky (+) / boxy (-), in the units of photutils' B4 (relative amplitude)


# ----------------------------------------------------------------------------------------------- engine
def _isophote_modules():
    from photutils.isophote import EllipseGeometry, Ellipse, build_ellipse_model
    return EllipseGeometry, Ellipse, build_ellipse_model


def _run_ellipse(arr, x0, y0, sma0, eps0, pa0_deg, **kw):
    EllipseGeometry, Ellipse, _ = _isophote_modules()
    g = EllipseGeometry(x0, y0, sma0, eps0, np.deg2rad(pa0_deg))
    try:
        ell = Ellipse(arr, geometry=g)
    except TypeError:                       # photutils < 3.0: positional geometry
        ell = Ellipse(arr, g)
    return ell.fit_image(**kw)


def fit_galaxy(data, x0, y0, mask=None, sma0=None, eps0=0.2, pa0=0.0, maxsma=None, minsma=0.0, step=0.1, linear=False,
               fix_center=False, fix_eps=False, fix_pa=False, conver=0.05, maxit=50, sclip=3.0, n_clip=0, integrmode='bilinear',
               zp=None, pixel_scale=1.0, with_model=True, high_harmonics=True):
    """Fit isophotes to one galaxy.  data: 2-D array; (x0, y0): 0-based pixel centre; mask: bool array (True = ignore).

    Returns (result, model).  result is a JSON-serialisable dict:
      status 'ok' | 'failed' | 'few' (< 6 isophotes), 'attempt' (which start parameters worked), 'table' (dict of lists,
      per isophote: sma, intens, intens_err, mu, mu_err, ellipticity(+err), pa(deg, CCW from +x, +err), x0, y0, a3, b3, a4, b4
      (+errors), growth_flux (mask-aware growth curve: flux inside each isophote), tflux_e (photutils' flux sum of the unmasked pixels), npix_e, rms, n_data, n_flag, stop_code), 'summary'.
    model: float32 array of the 2-D model (None if with_model is False or the fit failed)."""
    arr = np.ma.masked_array(np.asarray(data, dtype=float), mask=(np.zeros(data.shape, bool) if mask is None else mask.copy()))
    arr[~np.isfinite(np.asarray(data, dtype=float))] = np.ma.masked
    if maxsma is None:
        maxsma = 0.45 * min(arr.shape)
    if sma0 is None:
        sma0 = max(3.0, 0.03 * maxsma)
    kw = dict(minsma=minsma, maxsma=maxsma, step=step, linear=bool(linear), conver=conver, maxit=maxit, sclip=sclip,
              n_clip=n_clip, integrmode=integrmode, fix_center=fix_center, fix_eps=fix_eps, fix_pa=fix_pa)
    attempts = [(sma0, eps0, pa0), (sma0 * 2, 0.3, pa0 + 45.0), (sma0 * 0.6, 0.1, pa0 + 90.0)]
    isolist, used = None, None
    for k, (s0, e0, p0) in enumerate(attempts):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                iso = _run_ellipse(arr, x0, y0, s0, e0, p0 % 180.0, **kw)
        except Exception as exc:                                   # photutils raises on a hopeless start
            sys.stderr.write('isophote: attempt %d failed: %s\n' % (k, str(exc)[:80]))
            continue
        if iso is not None and len(iso) > 5:
            isolist, used = iso, k
            break
        if iso is not None and (isolist is None or len(iso) > len(isolist)):
            isolist, used = iso, k
    res = dict(status='failed', attempt=None, x0_in=float(x0), y0_in=float(y0), table={}, summary={})
    if isolist is None or len(isolist) < 2:
        return res, None
    res['status'] = 'ok' if len(isolist) > 5 else 'few'
    res['attempt'] = used
    t = isolist.to_table(columns='all')
    tab = {}
    for c in ('sma', 'intens', 'intens_err', 'ellipticity', 'ellipticity_err', 'pa', 'pa_err', 'x0', 'y0', 'a3', 'b3', 'a4', 'b4',
              'a3_err', 'b3_err', 'a4_err', 'b4_err', 'tflux_e', 'rms', 'n_data', 'n_flag', 'stop_code'):
        if c in t.colnames:
            tab[c] = [None if (v is None or (isinstance(v, float) and not math.isfinite(v))) else float(v) for v in
                      [(x.value if hasattr(x, 'value') else x) for x in t[c]]]
    tab['npix_e'] = [_f(getattr(i, 'npix_e', None)) for i in isolist]
    # surface brightness: mu = zp - 2.5 log10(I / pixel_area)   (I in counts per pixel; pixel_scale in arcsec)
    mu, mue = [], []
    for I, dI in zip(tab['intens'], tab['intens_err']):
        if zp is None or I is None or I <= 0:
            mu.append(None); mue.append(None)
        else:
            mu.append(float(zp - 2.5 * math.log10(I / pixel_scale ** 2)))
            mue.append(None if dI is None else float(1.0857 * dI / I))
    tab['mu'], tab['mu_err'] = mu, mue
    tab['growth_flux'] = growth_curve(tab['sma'], tab['intens'], tab['ellipticity'])
    res['table'] = tab
    res['summary'] = summarize(tab, zp=zp)
    model = None
    if with_model:
        try:
            _, _, bem = _isophote_modules()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                model = np.asarray(bem(arr.shape, isolist, high_harmonics=high_harmonics), dtype=np.float32)
        except Exception as exc:
            sys.stderr.write('isophote: model failed: %s\n' % str(exc)[:100])
    return res, model


def _f(v):
    return None if v is None or not math.isfinite(float(v)) else float(v)


def growth_curve(sma, intens, eps):
    """Mask-aware growth curve: flux inside the isophote of semi-major axis sma[k], obtained by integrating the fitted mean
    intensity over the elliptical annuli (area = pi a^2 (1 - eps)).  Pixels that were masked do not bias it (the isophote
    intensity is the mean of the unmasked pixels), unlike photutils' tflux_e which sums only the unmasked pixels.  List of floats."""
    s = [float(v) for v in sma]
    I = [float('nan') if v is None else float(v) for v in intens]
    e = [0.0 if v is None else float(v) for v in eps]
    out, F, A_prev, I_prev = [], 0.0, 0.0, None
    for k in range(len(s)):
        A = math.pi * s[k] ** 2 * (1.0 - e[k])
        if s[k] <= 0:
            out.append(0.0)
            A_prev, I_prev = 0.0, I[k]
            continue
        Ik = I[k] if math.isfinite(I[k]) else (I_prev if I_prev is not None else 0.0)
        Im = Ik if I_prev is None or not math.isfinite(I_prev) else 0.5 * (Ik + I_prev)
        F += Im * max(A - A_prev, 0.0)
        out.append(F)
        A_prev, I_prev = max(A, A_prev), Ik
    return out


def growth_radii(sma, flux, fractions=(0.5, 0.8), total=None):
    """Semi-major axes enclosing the given fractions of `total` (default: the flux inside the last isophote) by linear
    interpolation of the growth curve; NaN when the curve never reaches the fraction."""
    s = np.array([float(a) for a in sma])
    f = np.array([float('nan') if v is None else float(v) for v in flux])
    ok = np.isfinite(f) & (s > 0)
    s, f = s[ok], f[ok]
    if s.size < 3:
        return [float('nan')] * len(fractions)
    # growth curves from noisy outer isophotes are not monotonic: use the running maximum
    f = np.maximum.accumulate(f)
    tot = f[-1] if total is None else total
    out = []
    for fr in fractions:
        out.append(float(np.interp(fr * tot, f, s)) if f[-1] >= fr * tot and tot > 0 else float('nan'))
    return out


def summarize(tab, zp=None, disky_limit=DISKY_LIMIT):
    """Scalars from the profile table: half-light radius, total magnitude (flux inside the last good isophote), shape at
    the half-light radius and outside, boxy/disky classification."""
    sma = np.array(tab['sma'], dtype=float)
    ell = np.array([np.nan if v is None else v for v in tab['ellipticity']], dtype=float)
    pa = np.array([np.nan if v is None else v for v in tab['pa']], dtype=float)
    tfl = [np.nan if v is None else v for v in (tab.get('growth_flux') or tab['tflux_e'])]
    b4 = np.array([np.nan if v is None else v for v in tab['b4']], dtype=float)
    a4 = np.array([np.nan if v is None else v for v in tab['a4']], dtype=float)
    sc = np.array([-1 if v is None else v for v in tab['stop_code']], dtype=float)
    good = (sc >= 0) & (sc <= 2) & (sma > 0)          # photutils stop codes: 0 converged, 1 too many flagged, 2 too few, 3 gradient...
    r50, r80 = growth_radii(sma, tfl, (0.5, 0.8))
    ftot = float(np.nanmax(np.maximum.accumulate(np.nan_to_num(np.array(tfl, dtype=float), nan=0.0)))) if len(tfl) else float('nan')
    s = dict(n_iso=int(len(sma)), n_good=int(good.sum()), rmax=float(sma.max()), flux_total=ftot, r50=r50, r80=r80,
             mag_total=(float(zp - 2.5 * math.log10(ftot)) if (zp is not None and ftot > 0) else None))

    def around(r, lo, hi, arr, widen=(0.5, 2.0)):
        """mean of arr over the converged isophotes in [lo, hi] * r; widened once when there are none; NaN otherwise"""
        for a_, b_ in ((lo, hi), widen):
            sel = good & (sma >= a_ * r) & (sma <= b_ * r) & np.isfinite(arr)
            if sel.any():
                return float(np.mean(arr[sel]))
        return float('nan')

    ref = r50 if np.isfinite(r50) else float(np.median(sma))
    s['r_conv'] = float(sma[good].max()) if good.any() else float('nan')
    s['eps_hl'] = around(ref, 0.8, 1.25, ell)
    sel_pa = good & (sma >= 0.8 * ref) & (sma <= 1.25 * ref) & np.isfinite(pa)
    if not sel_pa.any():
        sel_pa = good & (sma >= 0.5 * ref) & (sma <= 2.0 * ref) & np.isfinite(pa)
    s['pa_hl'] = _circ_mean_pa(pa[sel_pa])
    out = good & (sma >= 0.75 * s['r_conv'])
    s['eps_out'] = float(np.nanmean(ell[out])) if out.any() else float('nan')
    s['pa_out'] = _circ_mean_pa(pa[out & np.isfinite(pa)])
    s['b4_hl'] = around(ref, 0.5, 1.5, b4, widen=(0.3, 3.0))
    s['a4_hl'] = around(ref, 0.5, 1.5, a4, widen=(0.3, 3.0))
    s['shape'] = 'disky' if s['b4_hl'] > disky_limit else ('boxy' if s['b4_hl'] < -disky_limit else 'elliptical') \
        if np.isfinite(s['b4_hl']) else ''
    return s


def _circ_mean_pa(pa_deg):
    """Mean of position angles (period 180 deg), in [0,180)."""
    if len(pa_deg) == 0:
        return float('nan')
    z = np.exp(2j * np.deg2rad(np.asarray(pa_deg, dtype=float)))
    return float((np.rad2deg(np.angle(z.mean())) / 2.0) % 180.0)


# ----------------------------------------------------------------------------------------------- catalog driver
def _neighbour_mask(shape, rows, skip_number, factor, box=None):
    """Mask ellipses of the other catalog objects (factor * A_IMAGE, B_IMAGE, THETA_IMAGE)."""
    ny, nx = shape
    m = np.zeros(shape, bool)
    for r in rows:
        if r['num'] == skip_number or not np.isfinite(r['a']) or not np.isfinite(r['b']):
            continue
        a, b = max(factor * r['a'], 2.0), max(factor * r['b'], 2.0)
        x0, y0 = r['x'], r['y']
        if box is not None and (x0 < box[0] - a or x0 > box[1] + a or y0 < box[2] - a or y0 > box[3] + a):
            continue
        xa, xb = int(max(0, x0 - a - 1)), int(min(nx, x0 + a + 2))
        ya, yb = int(max(0, y0 - a - 1)), int(min(ny, y0 + a + 2))
        if xa >= xb or ya >= yb:
            continue
        yy, xx = np.mgrid[ya:yb, xa:xb]
        t = np.deg2rad(r['theta'])
        dx, dy = xx - x0, yy - y0
        u = dx * np.cos(t) + dy * np.sin(t)
        v = -dx * np.sin(t) + dy * np.cos(t)
        m[ya:yb, xa:xb] |= (u / a) ** 2 + (v / b) ** 2 <= 1.0
    return m


def read_sources(catalog):
    cols, rows = tsvio.read_catalog(catalog)
    out = []
    for r in rows:
        try:
            num = int(float(r['NUMBER']))
        except (KeyError, ValueError):
            continue
        out.append(dict(num=num, x=tsvio.fnum(r.get('X_IMAGE')) - 1.0, y=tsvio.fnum(r.get('Y_IMAGE')) - 1.0,
                        a=tsvio.fnum(r.get('A_IMAGE'), 5.0), b=tsvio.fnum(r.get('B_IMAGE'), 5.0),
                        theta=tsvio.fnum(r.get('THETA_IMAGE'), 0.0), kron=tsvio.fnum(r.get('KRON_RADIUS'), 3.0),
                        mag=tsvio.fnum(r.get('MAG_AUTO'))))
    return out


def fit_source(data, src, mask=None, others=None, maxsma_scale=3.0, mask_neighbours=0.0, **kw):
    """Fit one catalog source on a cutout.  Returns (result, model_cutout, (y0, x0) of the cutout)."""
    ny, nx = data.shape
    a = src['a'] if np.isfinite(src['a']) else 5.0
    kron = src['kron'] if np.isfinite(src['kron']) else 3.0
    maxsma = kw.pop('maxsma', None) or max(10.0, maxsma_scale * kron * a)
    half = int(math.ceil(maxsma * 1.2)) + 2
    xa, xb = max(0, int(src['x']) - half), min(nx, int(src['x']) + half + 1)
    ya, yb = max(0, int(src['y']) - half), min(ny, int(src['y']) + half + 1)
    cut = data[ya:yb, xa:xb]
    m = np.zeros(cut.shape, bool)
    if mask is not None:
        m |= mask[ya:yb, xa:xb]
    if mask_neighbours and others:
        m |= _neighbour_mask(data.shape, others, src['num'], mask_neighbours, box=(xa, xb, ya, yb))[ya:yb, xa:xb]
    # never mask the galaxy centre itself (the shared mask contains the target as a detected source!)
    cx, cy = src['x'] - xa, src['y'] - ya
    yy, xx = np.mgrid[:cut.shape[0], :cut.shape[1]]
    own = ((xx - cx) ** 2 + (yy - cy) ** 2) <= (max(2.0, 0.5 * a)) ** 2
    m &= ~own
    pa0 = src['theta'] if np.isfinite(src['theta']) else 0.0
    eps0 = 1.0 - src['b'] / src['a'] if (np.isfinite(src['a']) and src['a'] > 0 and np.isfinite(src['b'])) else 0.2
    eps0 = float(min(max(eps0, 0.05), 0.8))
    res, model = fit_galaxy(cut, cx, cy, mask=m, sma0=kw.pop('sma0', None) or max(3.0, min(0.5 * a, 0.1 * maxsma)),
                            eps0=kw.pop('eps0', eps0), pa0=kw.pop('pa0', pa0), maxsma=maxsma, **kw)
    res['cutout'] = [int(ya), int(xa), int(yb), int(xb)]
    res['masked_fraction'] = float(m.mean())
    # profile coordinates back to the full image (0-based)
    for k, off in (('x0', xa), ('y0', ya)):
        if k in res['table']:
            res['table'][k] = [None if v is None else v + off for v in res['table'][k]]
    return res, model, (ya, xa)


def _work(args):
    data, src, mask, others, kw = args
    res, model, org = fit_source(data, src, mask=mask, others=others, **kw)
    res['number'] = src['num']
    return res, model, org


def run_catalog(data, sources, numbers=None, max_objects=10, mask=None, n_workers=0, **kw):
    """Fit several catalog sources; returns (results, model_image, residual_image).  The model is the sum of the per-galaxy
    models (overlapping galaxies add up: documented limitation)."""
    sel = [s for s in sources if (not numbers or s['num'] in set(numbers))]
    if not numbers:
        sel = sorted([s for s in sel if np.isfinite(s['mag'])], key=lambda s: s['mag'])[:max_objects] or sel[:max_objects]
    jobs = [(data, s, mask, sources, dict(kw)) for s in sel]
    if n_workers and n_workers > 1 and len(jobs) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(n_workers, len(jobs))) as ex:
            out = list(ex.map(_work, jobs))
    else:
        out = [_work(j) for j in jobs]
    model = np.zeros(data.shape, np.float32)
    for res, mod, (ya, xa) in out:
        if mod is not None:
            model[ya:ya + mod.shape[0], xa:xa + mod.shape[1]] += mod
    resid = (np.asarray(data, dtype=np.float32) - model)
    return [o[0] for o in out], model, resid


# ----------------------------------------------------------------------------------------------- outputs
def catalog_row(res):
    s = res.get('summary') or {}
    d = {c: None for c in COLUMNS}
    d['ISO_FLAG'] = {'ok': 0, 'few': 1, 'failed': 2}.get(res['status'], 2)
    if not s:
        return d
    t = res['table']
    d.update(ISO_NISO=s['n_iso'], ISO_RMAX=s['rmax'], ISO_X0=(t['x0'][-1] + 1.0 if t.get('x0') else None),
             ISO_Y0=(t['y0'][-1] + 1.0 if t.get('y0') else None), ISO_R50=s['r50'], ISO_R80=s['r80'], ISO_MAG_TOT=s['mag_total'],
             ISO_EPS_HL=s['eps_hl'], ISO_PA_HL=s['pa_hl'], ISO_EPS_OUT=s['eps_out'], ISO_PA_OUT=s['pa_out'],
             ISO_B4_HL=s['b4_hl'], ISO_A4_HL=s['a4_hl'], ISO_SHAPE=s['shape'], ISO_RCONV=s['r_conv'])
    return d


def profile_records(res):
    t = res['table']
    n = len(t.get('sma', []))
    recs = []
    for i in range(n):
        r = {k: t[k][i] for k in t if len(t[k]) == n}
        r['NUMBER'] = res['number']
        r['x0'] = None if r.get('x0') is None else r['x0'] + 1.0       # 1-based like the catalog
        r['y0'] = None if r.get('y0') is None else r['y0'] + 1.0
        recs.append(r)
    return recs


def plot_profiles(results, path, title=''):
    """Four-panel PNG (surface brightness, ellipticity, PA, B4 and the growth curve) for up to 6 galaxies."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(11, 6.2))
    ax = ax.ravel()
    for res in results[:6]:
        t = res['table']
        if not t or 'sma' not in t:
            continue
        s = np.array(t['sma'], dtype=float)
        lab = '#%s' % res['number']

        def arr(k):
            return np.array([np.nan if v is None else v for v in t[k]], dtype=float)
        ok = s > 0
        y = arr('mu') if np.isfinite(arr('mu')).any() else -2.5 * np.log10(np.clip(arr('intens'), 1e-12, None))
        ax[0].plot(s[ok] ** 0.25, y[ok], '.-', ms=3, label=lab)
        ax[1].errorbar(s[ok], arr('ellipticity')[ok], arr('ellipticity_err')[ok], fmt='.-', ms=3, lw=0.8)
        ax[2].errorbar(s[ok], arr('pa')[ok], arr('pa_err')[ok], fmt='.-', ms=3, lw=0.8)
        ax[3].errorbar(s[ok], arr('b4')[ok], arr('b4_err')[ok], fmt='.-', ms=3, lw=0.8)
        ax[4].plot(s[ok], np.maximum.accumulate(np.nan_to_num(arr('growth_flux')))[ok], '.-', ms=3)
        ax[5].plot(s[ok], arr('x0')[ok] - np.nanmedian(arr('x0')), '.-', ms=3, lw=0.8)
        ax[5].plot(s[ok], arr('y0')[ok] - np.nanmedian(arr('y0')), '.--', ms=3, lw=0.8)
    for a, (xl, yl) in zip(ax, [('sma^0.25', 'mu [mag/arcsec2]'), ('sma [pix]', 'ellipticity'), ('sma [pix]', 'PA [deg]'),
                                ('sma [pix]', 'B4 (+disky / -boxy)'), ('sma [pix]', 'flux inside isophote'),
                                ('sma [pix]', 'centre offset [pix]')]):
        a.set_xlabel(xl); a.set_ylabel(yl); a.grid(alpha=0.3)
    ax[0].invert_yaxis()
    ax[0].legend(fontsize=7)
    fig.suptitle(title or 'isophote profiles')
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--catalog')
    ap.add_argument('--xy', help='X,Y (1-based) of one galaxy instead of a catalog')
    ap.add_argument('--numbers', default='')
    ap.add_argument('--max-objects', type=int, default=10)
    ap.add_argument('--mask', default='')
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--pixel-scale', type=float, default=1.0)
    ap.add_argument('--maxsma-scale', type=float, default=3.0)
    ap.add_argument('--maxsma', type=float, default=0.0)
    ap.add_argument('--sma0', type=float, default=0.0)
    ap.add_argument('--eps0', type=float, default=0.0)
    ap.add_argument('--pa0', type=float, default=None)
    ap.add_argument('--step', type=float, default=0.1)
    ap.add_argument('--linear', action='store_true')
    ap.add_argument('--fix-center', action='store_true')
    ap.add_argument('--no-harmonics', action='store_true', help='model without the A3/B3/A4/B4 terms')
    ap.add_argument('--mask-neighbours', type=float, default=0.0, help='mask other catalog objects out to this multiple of A_IMAGE/B_IMAGE (0 = off)')
    ap.add_argument('--n-workers', type=int, default=0)
    ap.add_argument('--model-out', default='')
    ap.add_argument('--resid-out', default='')
    ap.add_argument('--table-out', default='')
    ap.add_argument('--json-out', default='')
    ap.add_argument('--plot-out', default='')
    a = ap.parse_args(argv)

    data, hdr = imageio.load_image(a.image)
    mask = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    if mask is not None:
        sys.stderr.write('isophote: mask %s (%.2f%% masked)\n' % (a.mask, 100.0 * mask.mean()))
    kw = dict(maxsma_scale=a.maxsma_scale, mask_neighbours=a.mask_neighbours, step=a.step, linear=a.linear,
              fix_center=a.fix_center, zp=a.mag_zeropoint, pixel_scale=a.pixel_scale, high_harmonics=not a.no_harmonics)
    if a.maxsma:
        kw['maxsma'] = a.maxsma
    if a.sma0:
        kw['sma0'] = a.sma0
    if a.eps0:
        kw['eps0'] = a.eps0
    if a.pa0 is not None:
        kw['pa0'] = a.pa0
    if a.xy:
        x, y = [float(v) for v in a.xy.split(',')]
        sources = [dict(num=1, x=x - 1.0, y=y - 1.0, a=10.0, b=8.0, theta=0.0, kron=3.0, mag=0.0)]
        numbers = []
    else:
        if not a.catalog:
            ap.error('--catalog or --xy is required')
        sources = read_sources(a.catalog)
        numbers = tsvio.parse_numbers(a.numbers)
    results, model, resid = run_catalog(data, sources, numbers=numbers, max_objects=a.max_objects, mask=mask,
                                        n_workers=a.n_workers, **kw)
    extra = dict(HISTORY='OGFinder isophote fit', OGF_NGAL=len(results))
    if a.model_out:
        imageio.save_fits(a.model_out, model, hdr)
    if a.resid_out:
        imageio.save_fits(a.resid_out, resid, hdr)
    recs = [r for res in results for r in profile_records(res)]
    if a.table_out:
        tsvio.write_table(a.table_out, PROFILE_COLS, recs)
    if a.json_out:
        with open(a.json_out, 'w') as f:
            json.dump({'results': results}, f)
    if a.plot_out and results:
        plot_profiles(results, a.plot_out, os.path.basename(a.image))
    tsvio.write_columns(COLUMNS, [(res['number'], catalog_row(res)) for res in results])
    ok = sum(1 for r in results if r['status'] == 'ok')
    sys.stderr.write('isophote: %d galaxies, %d fitted ok, %d isophotes\n' % (len(results), ok, len(recs)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
