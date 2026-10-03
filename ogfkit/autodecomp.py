"""Automatic structural decomposition: isophote profile -> initial guesses -> multi-component fit -> model selection.

    prof = clean_profile(isophote_table, fwhm)                    # sma, I, dI, eps, pa of an isophote fit (plugins/isophote)
    g = guess_from_profile(prof, fwhm, zp)                        # 1-D single-Sersic and Sersic+exponential fits, structure features, 2-D starting components
    out = decompose(data, x, y, psf=psf, rms=rms, mask=mask, ...) # profile (photutils isophote) + guesses + `multifit.fit_multistart` per candidate + BIC selection

Candidates (in order of complexity): `sersic`; `bulge+disk` (free-n Sersic bulge + exponential disc, common centre); `nucleus+bulge+disk` (adds a PSF component,
tried unless `nucleus='off'`; it must beat bulge+disk by the BIC margin and carry 0.5-50 % of the light).  A more complex model must win by `bic_margin` in the 2-D BIC and pass the sanity
tests (B/T between `bt_min` = 0.05 and `bt_max` = 0.97, bulge smaller than disc, components not stuck on a bound); otherwise the simpler model is kept (`note` says why).
All radii are along the major axis (the multifit convention), `pa` degrees counter-clockwise from +x, fluxes in counts.
Pure numpy / scipy; the isophote engine (photutils) is only needed by `isophote_profile`."""
import math
import os
import sys
import warnings

import numpy as np
from scipy.optimize import least_squares
from scipy.special import gamma, gammaincinv

from . import multifit as mf

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
TYPES = {'sersic': 1, 'bulge+disk': 2, 'nucleus+bulge+disk': 3}
FLAGS = dict(NOPROFILE=1, GUESS_FALLBACK=2, DEGENERATE=4, NOFIT=8, BAR=16, BOUND=32, EDGEON=64)


# ----------------------------------------------------------------------------------------------------------------- profiles
def bn(n):
    return float(gammaincinv(2.0 * n, 0.5))


def sersic_profile(r, ie, re, n):
    b = bn(n)
    return ie * np.exp(-b * ((np.maximum(r, 1e-6) / re) ** (1.0 / n) - 1.0))


def sersic_total(ie, re, n, q):
    """Total flux of a Sersic with intensity `ie` at the (major-axis) radius `re`, axis ratio q."""
    b = bn(n)
    return 2.0 * math.pi * n * math.exp(b) * b ** (-2.0 * n) * gamma(2.0 * n) * re * re * q * ie


def sersic_ie_from_flux(flux, re, n, q):
    return flux / sersic_total(1.0, re, n, q)


def exp_profile(r, i0, h):
    return i0 * np.exp(-np.asarray(r) / h)


def exp_total(i0, h, q):
    return 2.0 * math.pi * h * h * q * i0


def circ_median_pa(pa_deg):
    """Median of axial angles (period 180 deg)."""
    a = np.deg2rad(2.0 * np.asarray(pa_deg, float))
    return float((0.5 * np.rad2deg(np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a))))) % 180.0)


def pa_diff(a, b):
    d = (a - b + 90.0) % 180.0 - 90.0
    return abs(d)


def clean_profile(table, fwhm=3.0):
    """Arrays from the `table` dict of `isophote.fit_galaxy` (lists with None)."""
    def col(k):
        return np.array([np.nan if v is None else v for v in table.get(k, [])], float)
    sma, I, dI, eps, pa = col('sma'), col('intens'), col('intens_err'), col('ellipticity'), col('pa')
    gf = col('growth_flux') if 'growth_flux' in table else np.full(len(sma), np.nan)
    ok = np.isfinite(sma) & np.isfinite(I) & (sma > 0) & np.isfinite(eps)
    dI = np.where(np.isfinite(dI), dI, np.nanmedian(np.abs(dI)) if np.any(np.isfinite(dI)) else 0.0)
    o = np.argsort(sma[ok])
    return dict(sma=sma[ok][o], I=I[ok][o], dI=dI[ok][o], eps=eps[ok][o], pa=pa[ok][o], growth=gf[ok][o], fwhm=float(fwhm))


def half_light_radius(prof):
    """Semi-major axis enclosing half of the total flux of the growth curve (total = plateau = last value)."""
    g = prof['growth']
    if not np.any(np.isfinite(g)):
        return float('nan'), float('nan')
    ok = np.isfinite(g)
    s, g = prof['sma'][ok], g[ok]
    tot = np.max(g)
    if not tot > 0:
        return float('nan'), float(tot)
    k = int(np.argmax(g >= 0.5 * tot))
    if k == 0:
        return float(s[0]), float(tot)
    return float(np.interp(0.5 * tot, g[k - 1:k + 1], s[k - 1:k + 1])), float(tot)


# --------------------------------------------------------------------------------------------------------------- 1-D fits
def _sigma(prof, rel=0.03):
    return np.sqrt(prof['dI'] ** 2 + (rel * np.abs(prof['I'])) ** 2) + 1e-30


def _bic(chi2, k, n):
    return float(chi2 + k * math.log(max(n, 2)))


def fit_single_1d(prof, r50, rmin):
    m = prof['sma'] >= rmin
    r, I, sg = prof['sma'][m], prof['I'][m], _sigma(prof)[m]
    if len(r) < 6:
        return None
    Ie0 = float(np.interp(r50, prof['sma'], prof['I']))
    best = None
    for n0 in (1.0, 2.5, 4.0):
        def res(p):
            ie, re, n, s = np.exp(p[0]), np.exp(p[1]), p[2], p[3]
            return (I - sersic_profile(r, ie, re, n) - s) / sg
        try:
            f = least_squares(res, [math.log(max(Ie0, 1e-12)), math.log(max(r50, 1.0)), n0, 0.0], bounds=([-60, math.log(0.3), 0.4, -np.inf], [60, math.log(20 * max(r50, 2) + 50), 8.0, np.inf]), max_nfev=300)
        except Exception:
            continue
        if best is None or f.cost < best.cost:
            best = f
    if best is None:
        return None
    p = best.x
    chi2 = 2 * best.cost
    return dict(ie=float(np.exp(p[0])), re=float(np.exp(p[1])), n=float(p[2]), sky=float(p[3]), chi2=float(chi2), npts=int(len(r)), bic=_bic(chi2, 4, len(r)))


def fit_double_1d(prof, r50, rmin, single=None):
    m = prof['sma'] >= rmin
    r, I, sg = prof['sma'][m], prof['I'][m], _sigma(prof)[m]
    if len(r) < 9:
        return None
    # disc start: log-linear slope of the outer third of the points
    k0 = int(len(r) * 0.55)
    sl = np.polyfit(r[k0:], np.log(np.maximum(I[k0:], 1e-12 * np.max(I))), 1)
    h0 = float(np.clip(-1.0 / sl[0] if sl[0] < -1e-6 else r50 / 1.678, 0.5 * r50 / 1.678, 4 * r50))
    i00 = float(np.exp(sl[1])) if np.isfinite(sl[1]) else float(np.max(I)) * 0.3
    best = None
    starts = []
    for fb, reb in ((0.15, 0.2), (0.3, 0.35), (0.5, 0.5)):
        for nb in (1.5, 3.5):
            starts.append((fb, reb, nb))
    Imax = float(np.max(I))
    for fb, reb, nb in starts:
        reb_px = max(reb * r50, 0.6)
        ieb = max(float(np.interp(reb_px, prof['sma'], prof['I'])) * fb * 2, 1e-12)
        p0 = [math.log(ieb), math.log(reb_px), nb, math.log(max(min(i00, Imax), 1e-12)), math.log(h0), 0.0]
        lo = [-60, math.log(0.3), 0.5, -60, math.log(max(0.8, 0.1 * r50)), -np.inf]
        hi = [60, math.log(max(r50, 3.0) * 1.2), 7.0, 60, math.log(20 * max(r50, 2) + 50), np.inf]

        def res(p):
            ieb_, reb_, nb_, i0_, h_, s_ = np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3]), np.exp(p[4]), p[5]
            return (I - sersic_profile(r, ieb_, reb_, nb_) - exp_profile(r, i0_, h_) - s_) / sg
        try:
            f = least_squares(res, np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6), bounds=(lo, hi), max_nfev=400)
        except Exception:
            continue
        # a valid decomposition: bulge more concentrated than the disc
        ieb_, reb_, nb_, i0_, h_ = np.exp(f.x[0]), np.exp(f.x[1]), f.x[2], np.exp(f.x[3]), np.exp(f.x[4])
        valid = reb_ < 1.678 * h_ * 1.05
        key = f.cost + (0 if valid else 1e6)
        if best is None or key < best[0]:
            best = (key, f)
    if best is None:
        return None
    f = best[1]
    p = f.x
    chi2 = 2 * f.cost
    return dict(ie_b=float(np.exp(p[0])), re_b=float(np.exp(p[1])), n_b=float(p[2]), i0_d=float(np.exp(p[3])), h_d=float(np.exp(p[4])), sky=float(p[5]),
                chi2=float(chi2), npts=int(len(r)), bic=_bic(chi2, 6, len(r)))


# ------------------------------------------------------------------------------------------------------------- features
def structure_features(prof, fwhm, r50=None):
    """Shape information of the isophote profile: outer / inner ellipticity and PA, concentration, bar-like signature."""
    s, eps, pa = prof['sma'], prof['eps'], prof['pa']
    rmax = float(s.max())
    if r50 is None or not np.isfinite(r50):
        r50 = half_light_radius(prof)[0]
    f = dict(rmax=rmax, r50=r50)
    outer = s >= max(1.2 * r50, 0.5 * rmax)
    if outer.sum() < 2:
        outer = s >= 0.6 * rmax
    f['eps_out'] = float(np.nanmedian(eps[outer])) if outer.any() else float('nan')
    f['pa_out'] = circ_median_pa(pa[outer]) if outer.any() else float('nan')
    inner = (s >= 0.7 * fwhm) & (s <= max(0.5 * r50, 1.5 * fwhm))
    f['eps_in'] = float(np.nanmedian(eps[inner])) if inner.any() else float('nan')
    f['pa_in'] = circ_median_pa(pa[inner]) if inner.any() else float('nan')
    # concentration: I(0.3 r50) / I(r50)
    try:
        f['conc'] = float(np.interp(0.3 * r50, s, prof['I']) / max(np.interp(r50, s, prof['I']), 1e-30))
    except Exception:
        f['conc'] = float('nan')
    # bar signature: ellipticity maximum inside the disc followed by a drop, with a PA twist across it
    f['bar'] = False
    f['eps_max'] = float(np.nanmax(eps)) if len(eps) else float('nan')
    sel = (s > 1.5 * fwhm) & (s < 0.7 * rmax)
    if sel.sum() >= 5:
        k = int(np.nanargmax(np.where(sel, eps, -1)))
        after = (s > s[k]) & (s <= min(1.6 * s[k], rmax))
        if after.any() and eps[k] > 0.35 and eps[k] - np.nanmin(eps[after]) > 0.12:
            tw = pa_diff(pa[k], pa[after][int(np.nanargmin(eps[after]))])
            f['bar_twist'] = float(tw)
            f['bar'] = bool(tw > 20.0)
    f['edgeon'] = bool(np.isfinite(f['eps_out']) and f['eps_out'] > 0.75)
    return f


def nucleus_excess(prof, model_I, fwhm):
    """Ratio of the observed to the modelled intensity in the central isophotes (sma < 1.2 FWHM): a point source makes it > 1 where an extended profile
    (smoothed by the PSF) would make it <= 1.  `model_I` = model profile at prof['sma']."""
    m = prof['sma'] <= 1.2 * fwhm
    if m.sum() < 1:
        return float('nan')
    return float(np.median(prof['I'][m] / np.maximum(model_I[m], 1e-30)))


def guess_from_profile(prof, fwhm, zp=25.0, nucleus='auto'):
    """1-D fits and 2-D starting components (major-axis radii).  Returns dict: ok, r50, flux_tot, single, double, features, nucleus_ratio, comps{name: [...]},
    ties{name: [...]}, bt_guess, bic1d {single, double}, guess_type."""
    r50, ftot = half_light_radius(prof)
    out = dict(ok=False, r50=r50, flux_tot=ftot, fwhm=fwhm)
    if not (np.isfinite(r50) and r50 > 0 and ftot > 0 and len(prof['sma']) >= 8):
        return out
    rmin = 0.8 * fwhm
    feat = structure_features(prof, fwhm, r50)
    out['features'] = feat
    s1 = fit_single_1d(prof, r50, rmin)
    s2 = fit_double_1d(prof, r50, rmin, s1)
    out['single'], out['double'] = s1, s2
    q_out = float(np.clip(1.0 - feat['eps_out'], 0.1, 1.0)) if np.isfinite(feat['eps_out']) else 0.8
    pa_out = feat['pa_out'] if np.isfinite(feat['pa_out']) else 0.0
    q_in = float(np.clip(1.0 - feat['eps_in'], 0.25, 1.0)) if np.isfinite(feat['eps_in']) else min(1.0, q_out + 0.2)
    pa_in = feat['pa_in'] if np.isfinite(feat['pa_in']) else pa_out
    comps, ties = {}, {}
    # single Sersic
    if s1 is not None:
        n1, re1 = float(np.clip(s1['n'], 0.5, 7.5)), s1['re']
        f1 = ftot if ftot > 0 else sersic_total(s1['ie'], re1, n1, q_out)
        comps['sersic'] = [dict(kind='sersic', flux=f1, re=re1, n=n1, q=q_out, pa=pa_out)]
    else:
        comps['sersic'] = [dict(kind='sersic', flux=ftot, re=r50, n=2.0, q=q_out, pa=pa_out)]
    # bulge + disc
    bt = None
    if s2 is not None:
        fb = sersic_total(s2['ie_b'], s2['re_b'], s2['n_b'], q_in)
        fd = exp_total(s2['i0_d'], s2['h_d'], q_out)
        bt = fb / (fb + fd) if fb + fd > 0 else None
    degenerate = s2 is None or bt is None or not (0.03 < bt < 0.97)
    if degenerate:                                                      # heuristic start
        reb, nb, fb_frac, hd = max(0.25 * r50, 0.8 * fwhm), 2.5, 0.3, r50 / 1.678 * 1.1
        bt_g = fb_frac
    else:
        reb, nb, bt_g, hd = s2['re_b'], s2['n_b'], bt, s2['h_d']
    out['bt_guess'] = float(bt_g)
    out['guess_fallback'] = bool(degenerate)
    out['bt_1d'] = None if bt is None else float(bt)
    comps['bulge+disk'] = [dict(kind='sersic', flux=bt_g * ftot, re=float(max(reb, 0.5)), n=float(np.clip(nb, 0.6, 6.5)), q=q_in, pa=pa_in),
                           dict(kind='exp', flux=(1 - bt_g) * ftot, re=float(1.678 * hd), q=q_out, pa=pa_out)]
    ties['bulge+disk'] = [('1.x', '0.x'), ('1.y', '0.y')]
    # nucleus: the candidate is always offered (the central excess of a compact bulge and of a point source are degenerate in the profile)
    nr = float('nan')
    if s2 is not None and not degenerate:
        mI = sersic_profile(prof['sma'], s2['ie_b'], s2['re_b'], s2['n_b']) + exp_profile(prof['sma'], s2['i0_d'], s2['h_d']) + s2['sky']
        nr = nucleus_excess(prof, mI, fwhm)
    elif s1 is not None:
        nr = nucleus_excess(prof, sersic_profile(prof['sma'], s1['ie'], s1['re'], s1['n']) + s1['sky'], fwhm)
    out['nucleus_ratio'] = nr
    want_nuc = nucleus != 'off'
    out['nucleus_candidate'] = bool(want_nuc)
    if want_nuc:
        fn = 0.03 * ftot
        bdc = comps['bulge+disk']
        comps['nucleus+bulge+disk'] = [dict(kind='psf', flux=fn)] + [dict(c, flux=c['flux'] * (1 - 0.03)) for c in bdc]
        ties['nucleus+bulge+disk'] = [('0.x', '1.x'), ('0.y', '1.y'), ('2.x', '1.x'), ('2.y', '1.y')]
    b1, b2 = (s1['bic'] if s1 else float('nan')), (s2['bic'] if s2 else float('nan'))
    out['bic1d'] = dict(single=b1, double=b2)
    out['guess_type'] = 'bulge+disk' if (s2 is not None and s1 is not None and b2 < b1 - 10 and not degenerate) else 'sersic'
    out['comps'], out['ties'] = comps, ties
    out['ok'] = True
    return out


# ------------------------------------------------------------------------------------------------------------- isophote engine
_ISO = None


def _iso_module():
    global _ISO
    if _ISO is None:
        import importlib.util
        p = os.path.join(ROOT, 'plugins', 'isophote', 'isophote.py')
        spec = importlib.util.spec_from_file_location('ogf_isophote_engine', p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _ISO = m
    return _ISO


def isophote_profile(data, x, y, mask=None, sma0=None, eps0=0.2, pa0=0.0, maxsma=None, zp=None, fwhm=3.0, step=0.1):
    """Isophote fit of the object (0-based centre x, y in `data`).  Returns the cleaned profile dict or None."""
    iso = _iso_module()
    maxsma = maxsma or 0.45 * min(data.shape)
    starts = [sma0 or max(2.0, 1.2 * fwhm)] + [f * maxsma for f in (0.25, 0.12)]       # photutils needs a start radius where the first isophote converges
    for s0 in starts:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res, _ = iso.fit_galaxy(data, x, y, mask=mask, sma0=s0, eps0=eps0, pa0=pa0, maxsma=maxsma, step=step, zp=zp, pixel_scale=1.0, with_model=False, high_harmonics=False)
        except Exception:
            continue
        if res.get('status') in ('ok', 'few') and res.get('table'):
            prof = clean_profile(res['table'], fwhm)
            if len(prof['sma']) >= 6:
                return prof
    return None


# ------------------------------------------------------------------------------------------------------------- 2-D decomposition
def _bound_comp(c, half, fwhm):
    c = dict(c)
    b = dict(c.get('bounds', {}))
    if c['kind'] == 'sersic':
        b.setdefault('n', (0.5, 7.5))
        b.setdefault('re', (0.3, 1.2 * half))
        b.setdefault('q', (0.2, 1.0))
    elif c['kind'] == 'exp':
        b.setdefault('re', (0.8, 3.0 * half))
        b.setdefault('q', (0.08, 1.0))
    c['bounds'] = b
    return c


def sanity(name, res, ntarget, half, bt_min=0.05, bt_max=0.97):
    """-> (ok, reason) for a decomposition result (first `ntarget` components are the target's)."""
    cs = res['components'][:ntarget]
    if name == 'sersic':
        return True, ''
    bulge = cs[-2]
    disc = cs[-1]
    fb, fd = bulge['flux'], disc['flux']
    if not (fb > 0 and fd > 0):
        return False, 'negative component flux'
    bt = fb / (fb + fd)
    if bt < bt_min or bt > bt_max:
        return False, 'B/T %.3f outside [%.2f, %.2f]' % (bt, bt_min, bt_max)
    if bulge['re'] > disc['re'] * 1.05:
        return False, 'bulge larger than disc'
    if bulge['re'] < 0.4:
        return False, 'bulge shrunk to a point'
    nb = bulge.get('n', 2.0)
    if nb > 7.3 or nb < 0.52:
        return False, 'bulge Sersic index on its bound'
    if name == 'nucleus+bulge+disk':
        if cs[0]['flux'] < 0.005 * (fb + fd) or cs[0]['flux'] > 0.5 * (fb + fd):
            return False, 'nucleus flux implausible'
    return True, ''


def decompose(data, x, y, psf=None, rms=1.0, mask=None, sky='const', zp=25.0, flux0=None, re0=None, q0=0.8, pa0=0.0, fwhm=3.0, extra=(), tie_extra=(),
              models=('sersic', 'bulge+disk', 'nucleus+bulge+disk'), bic_margin=10.0, nucleus='auto', restarts=1, max_nfev=150, profile=None, use_profile=True,
              maxsma=None, bt_min=0.05, bt_max=0.97, gain=None, **fit_kw):
    """Run the whole chain for one object.  `x`, `y` are 0-based centre coordinates in `data`; `flux0`, `re0`, `q0`, `pa0` catalog start values (used when the
    isophote profile fails).  `extra` / `tie_extra`: neighbour components fitted simultaneously (appended after the target's).
    Returns dict: name (selected), res (multifit result of the selection), ntarget, all {name: result | None}, notes, flag, guess, features, bic {name: bic}."""
    ny, nx = data.shape
    half = 0.5 * min(ny, nx)
    flag = 0
    notes = []
    if maxsma is None:
        maxsma = 0.45 * min(ny, nx)
    prof = profile
    if prof is None and use_profile:
        prof = isophote_profile(data, x, y, mask=mask, sma0=max(2.0, 1.2 * fwhm), eps0=min(max(1 - q0, 0.05), 0.7), pa0=pa0, maxsma=maxsma, zp=zp, fwhm=fwhm)
    g = guess_from_profile(prof, fwhm, zp, nucleus) if prof is not None else dict(ok=False)
    if not g.get('ok'):
        flag |= FLAGS['NOPROFILE']
        notes.append('no usable isophote profile: preset starting values')
        f0 = flux0 if (flux0 and flux0 > 0) else max(float(np.nansum(np.clip(data - np.nanmedian(data), 0, None))), 1.0)
        r0 = re0 if (re0 and re0 > 0) else 3.0
        comps = {'sersic': [dict(kind='sersic', flux=f0, re=r0, n=2.0, q=q0, pa=pa0)],
                 'bulge+disk': [dict(kind='sersic', flux=0.3 * f0, re=0.4 * r0, n=2.5, q=min(1.0, q0 + 0.2), pa=pa0), dict(kind='exp', flux=0.7 * f0, re=1.2 * r0, q=q0, pa=pa0)]}
        ties = {'bulge+disk': [('1.x', '0.x'), ('1.y', '0.y')]}
        g = dict(ok=False, comps=comps, ties=ties)
        if nucleus != 'off':
            comps['nucleus+bulge+disk'] = [dict(kind='psf', flux=0.03 * f0)] + [dict(c, flux=c['flux'] * 0.97) for c in comps['bulge+disk']]
            ties['nucleus+bulge+disk'] = [('0.x', '1.x'), ('0.y', '1.y'), ('2.x', '1.x'), ('2.y', '1.y')]
        feat = {}
    else:
        comps, ties, feat = g['comps'], g['ties'], g['features']
        if g.get('guess_fallback'):
            flag |= FLAGS['GUESS_FALLBACK']
        if feat.get('bar'):
            flag |= FLAGS['BAR']
        if feat.get('edgeon'):
            flag |= FLAGS['EDGEON']
    allr, bic, ntg = {}, {}, {}
    kw = dict(psf=psf, rms=rms, mask=mask, sky=sky, zp=zp, gain=gain, max_nfev=max_nfev, **fit_kw)
    for name in models:
        if name not in comps or (name == 'nucleus+bulge+disk' and nucleus == 'off'):
            continue
        start = []
        for c in comps[name]:
            c = dict(c)
            c['x'], c['y'] = float(x), float(y)
            c['bounds'] = dict(c.get('bounds', {}))
            if c['kind'] != 'psf':
                c = _bound_comp(c, half, fwhm)
            c['bounds'].update(x=(x - 3.0, x + 3.0), y=(y - 3.0, y + 3.0)) if c['kind'] != 'psf' else c['bounds'].update(x=(x - 3.0, x + 3.0), y=(y - 3.0, y + 3.0))
            start.append(c)
        ntg[name] = len(start)
        try:
            r = mf.fit_multistart(data, start + [dict(e) for e in extra], restarts=restarts, tie=list(ties.get(name, ())) + list(tie_extra), **kw)
        except Exception as e:
            allr[name] = None
            notes.append('%s fit failed: %s' % (name, str(e)[:60]))
            continue
        allr[name] = r
        bic[name] = r['bic']
    order = [m for m in models if allr.get(m) is not None]
    if not order:
        return dict(name=None, res=None, ntarget=0, all=allr, notes=notes, flag=flag | FLAGS['NOFIT'], guess=g, features=feat, bic=bic, profile=prof)
    # Real galaxies leave structured residuals (spiral arms, clumps, dust): the pixel errors are then too small and the plain BIC accepts any extra component.
    # Errors are inflated so that the best model has chi2_red = 1 (never deflated): BIC' = chi2 / max(1, min chi2_red) + k ln N.
    scale = max(1.0, min(allr[m]['chi2_red'] for m in order))
    bic_raw = dict(bic)
    for m in order:
        r = allr[m]
        penalty = r['bic'] - r['chi2']
        bic[m] = r['chi2'] / scale + penalty
        r['bic_raw'], r['bic'] = r['bic'], bic[m]
    best = order[0]
    for m in order[1:]:
        if allr[m]['bic'] < allr[best]['bic'] - bic_margin:
            ok, why = sanity(m, allr[m], ntg[m], half, bt_min, bt_max)
            if ok:
                best = m
            else:
                flag |= FLAGS['DEGENERATE']
                notes.append('%s rejected: %s (dBIC %.1f)' % (m, why, allr[best]['bic'] - allr[m]['bic']))
    return dict(name=best, res=allr[best], ntarget=ntg[best], all=allr, notes=notes, flag=flag, guess=g, features=feat, bic=bic, bic_raw=bic_raw, scale=scale, ntg=ntg, profile=prof)
