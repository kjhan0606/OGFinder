"""Bar / ring / spiral detection and fitting on top of a bulge + disc model.

fit_structure() fits a base model (Sersic bulge + exponential disc, common centre), looks at the residual for a bar (second-moment ellipse of the positive residual inside the disc),
a ring (azimuthal profile of the residual in the disc frame) and the spiral arms (disc with GALFIT's logarithmic coordinate rotation + m = 2 mode, started with both winding senses),
refits the base model PLUS the feature with all parameters free, and accepts a feature only if it lowers the error-rescaled BIC by `bic_margin` and passes sanity cuts
(bar: 3-60 % of the light, q < 0.7, length < 2.2 disc R_e; ring: 2-60 % of the light, width > 0.5 px, radius > 1.5 widths; spiral: |rotation| > 25 deg, m = 2 amplitude > 0.03).
The BIC errors are inflated so that the best model has chi2_red >= 1 (see autodecomp): real galaxies leave structured residuals and a plain BIC accepts any component.
Models: base, +bar, +ring, +spiral (bar and ring together are not tried).  Features that are not 'wanted' (features=...) are skipped.
"""
import math

import numpy as np

from . import multifit as MF

TIE3 = [('1.x', '0.x'), ('1.y', '0.y'), ('2.x', '0.x'), ('2.y', '0.y')]
TIE4 = TIE3 + [('3.x', '0.x'), ('3.y', '0.y')]
TIE2 = [('1.x', '0.x'), ('1.y', '0.y')]
def outer_geometry(data, mask, x, y, r_in=0.45, r_out=0.95, iters=6):
    """Axis ratio and position angle (ccw from +x, degrees) of the outer disc from flux-weighted second moments of the (smoothed, masked) image in an elliptical annulus
    r_in-r_out of the half-size of the cutout, iterated so that the annulus follows the ellipse.  -> (q, pa) or None.  Used when the bar / ring dominate the inner light and a
    free bulge + disc fit would take the bar for the disc."""
    ny, nx = data.shape
    half = 0.5 * min(nx, ny)
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x, yy - y
    good = np.isfinite(data) & ~(np.zeros(data.shape, bool) if mask is None else np.asarray(mask, bool))
    sm = _smooth(np.where(good, data, 0.0), 2.0)
    q, pa = 1.0, 0.0
    for _ in range(iters):
        t = math.radians(pa)
        xr = dx * math.cos(t) + dy * math.sin(t)
        yr = (-dx * math.sin(t) + dy * math.cos(t)) / max(q, 0.2)
        r = np.hypot(xr, yr)
        sel = good & (r > r_in * half) & (r < r_out * half) & (np.hypot(dx, dy) < 0.98 * half)
        w = np.clip(sm[sel] - np.median(sm[sel]) * 0.0, 0, None)
        if sel.sum() < 50 or w.sum() <= 0:
            return None
        sx = (w * dx[sel] ** 2).sum() / w.sum()
        sy = (w * dy[sel] ** 2).sum() / w.sum()
        sxy = (w * dx[sel] * dy[sel]).sum() / w.sum()
        ev, evec = np.linalg.eigh(np.array([[sx, sxy], [sxy, sy]]))
        pa = math.degrees(math.atan2(evec[1, 1], evec[0, 1]))
        pa = ((pa + 90.0) % 180.0) - 90.0
        q = float(np.clip(math.sqrt(max(ev[0], 1e-9) / max(ev[1], 1e-9)), 0.2, 1.0))
    return q, pa


FLAG = dict(BASE_FAIL=1, NOFIT=2, BAR_REJECT=4, RING_REJECT=8, SPIRAL_REJECT=16)


def base_comps(x, y, flux, re, q, pa):
    return [dict(kind='sersic', x=x, y=y, flux=0.25 * flux, re=0.35 * re, n=2.0, q=min(1.0, q + 0.15), pa=pa, bounds=dict(q=(0.3, 1.0), n=(0.5, 6.0))),
            dict(kind='exp', x=x, y=y, flux=0.75 * flux, re=1.2 * re, q=q, pa=pa)]


def _smooth(a, s=1.2):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(a, s)


def residual_bar_start(resid, good, comps, re_d):
    """Second-moment ellipse of the positive, smoothed residual within 1.4 disc R_e of the centre -> (pa_deg (ccw from +x), length, q, flux)."""
    ny, nx = resid.shape
    x0, y0 = comps[0]['x'], comps[0]['y']
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0)
    rs = _smooth(np.where(good, resid, 0.0))
    w = np.where((r < 1.4 * re_d) & (r > 1.5), np.clip(rs, 0, None), 0.0)
    if w.sum() <= 0:
        return None
    sx = (w * (xx - x0) ** 2).sum() / w.sum()
    sy = (w * (yy - y0) ** 2).sum() / w.sum()
    sxy = (w * (xx - x0) * (yy - y0)).sum() / w.sum()
    ev, evec = np.linalg.eigh(np.array([[sx, sxy], [sxy, sy]]))
    l1, l2 = max(ev[1], 1e-6), max(ev[0], 1e-6)
    pa = math.degrees(math.atan2(evec[1, 1], evec[0, 1]))
    pa = ((pa + 90.0) % 180.0) - 90.0
    return pa, 2.4 * math.sqrt(l1), min(max(math.sqrt(l2 / l1), 0.12), 0.6), float(w.sum())


def residual_ring_start(resid, good, comps, re_d):
    """Radius of the strongest positive annulus of the residual in the disc frame -> (r_ring, sigma, flux) or None."""
    ny, nx = resid.shape
    d = comps[1]
    x0, y0, q, pa = d['x'], d['y'], max(d.get('q', 1.0), 0.2), math.radians(d.get('pa', 0.0))
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x0, yy - y0
    xr = dx * math.cos(pa) + dy * math.sin(pa)
    yr = (-dx * math.sin(pa) + dy * math.cos(pa)) / q
    r = np.hypot(xr, yr)
    edges = np.arange(1.0, min(2.5 * re_d, 0.45 * min(nx, ny)), 1.5)
    prof, cnt = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = good & (r >= a) & (r < b)
        prof.append(resid[m].mean() if m.sum() > 5 else 0.0)
        cnt.append(m.sum())
    if not prof:
        return None
    prof = np.array(prof)
    k = int(np.argmax(prof))
    if prof[k] <= 0:
        return None
    rr = 0.5 * (edges[k] + edges[k + 1])
    half = prof[k] / 2
    lo = k
    while lo > 0 and prof[lo - 1] > half:
        lo -= 1
    hi = k
    while hi < len(prof) - 1 and prof[hi + 1] > half:
        hi += 1
    sg = max(0.5 * (edges[hi + 1] - edges[lo]) / 2.355 * 1.2, 1.2)
    flux = prof[k] * 2 * math.pi * rr * q * sg * math.sqrt(2 * math.pi)
    return rr, sg, max(float(flux), 1e-3)


def _bic(res, scale):
    return res['chi2'] / scale + res['nfree'] * math.log(max(res['npix'], 2))


def _flux(c):
    return 10 ** (-0.4 * (c['mag']))             # relative (zero point cancels in fractions)


def sanity(name, comps, re_d):
    """-> (ok, reason)."""
    tot = sum(_flux(c) for c in comps if np.isfinite(c.get('mag', np.nan)))
    if name == 'bar':
        b = comps[2]
        fr = _flux(b) / tot
        if not (0.03 <= fr <= 0.6):
            return False, 'bar flux fraction %.2f' % fr
        if b['q'] >= 0.7:
            return False, 'bar q %.2f' % b['q']
        if b['rout'] > 2.2 * max(comps[1]['re'], 1.0):
            return False, 'bar longer than the disc'
        if b['rout'] < 1.5:
            return False, 'bar shorter than 1.5 px'
    elif name == 'bar+ring':
        ok1, w1 = sanity('bar', comps[:3], re_d)
        ok2, w2 = sanity('ring', [comps[0], comps[1], comps[3]], re_d)
        return (ok1 and ok2), w1 or w2
    elif name == 'ring':
        b = comps[2]
        fr = _flux(b) / tot
        if not (0.02 <= fr <= 0.6):
            return False, 'ring flux fraction %.2f' % fr
        if b['sring'] < 0.5 or b['rring'] < 1.5 * b['sring']:
            return False, 'ring width/radius %.1f/%.1f' % (b['sring'], b['rring'])
    elif name == 'spiral':
        d = comps[1]
        if abs(d.get('rot_theta', 0.0)) < 25.0 or abs(d.get('f2a', 0.0)) < 0.03:
            return False, 'spiral rotation %.0f amplitude %.2f' % (d.get('rot_theta', 0.0), d.get('f2a', 0.0))
    b0 = comps[0]
    if comps[0]['re'] > comps[1]['re'] * 1.0 or b0['mag'] != b0['mag']:
        return False, 'bulge larger than disc'
    return True, ''


def fit_structure(data, x, y, flux, re, q, pa, features=('bar', 'ring', 'spiral'), bic_margin=10.0, max_nfev=300, geometry='given', **kw):
    """-> dict(name in base|bar|ring|spiral, res, all, bic, scale, notes, flag).  kw: psf, rms, mask, sky, zp, gain (as multifit.fit)."""
    notes = []
    flag = 0
    if geometry == 'outer':                                                     # disc axis ratio / PA from the outer isophotes (bar and ring do not bias them)
        og = outer_geometry(data, kw.get('mask'), x, y)
        if og is not None:
            q, pa = og
            notes.append('disc geometry from the outer isophotes: q %.2f, PA %.0f' % og)
    try:
        base = MF.fit_multistart(data, base_comps(x, y, flux, re, q, pa), restarts=1, tie=TIE2, max_nfev=max_nfev, **kw)
    except Exception as e:
        return dict(name='none', res=None, all={}, bic={}, scale=1.0, notes=['base fit failed: %s' % e], flag=FLAG['BASE_FAIL'])
    cs = base['components']
    ok_base = cs[0]['re'] <= cs[1]['re']
    if not ok_base:                                                 # swap labels: the compact one is the bulge
        cs = [cs[1], cs[0]]
    re_d = max(cs[1]['re'], 2.0)
    allr = {'base': base}
    cands = {}
    resid, good = base['residual'], base['good']
    f_tot = flux

    def as_start(c):
        d = {k: v for k, v in c.items() if k not in ('errors', 'free', 'mag', 'at_bound')}
        d['flux'] = d.get('flux', 10 ** (-0.4 * (c['mag'] - kw.get('zp', 25.0)))) if 'flux' in d else d.get('flux')
        return d

    start_base = [as_start(c) for c in cs[:2]]
    if 'bar' in features:
        bs = residual_bar_start(resid, good, cs, re_d)
        if bs is not None:
            pa_b, L, qb, fl = bs
            for dpa in (0.0,):
                bar = MF.bar_component(cs[0]['x'], cs[0]['y'], max(fl * 0.8, 1e-3 * f_tot), min(max(L, 3.0), 2.0 * re_d), pa_b, q=qb)
                cands.setdefault('bar', []).append(start_base + [bar])
            bar2 = MF.bar_component(cs[0]['x'], cs[0]['y'], max(0.1 * f_tot, 1e-3), min(max(0.8 * re_d, 3.0), 2.0 * re_d), pa_b + 90.0 if False else cs[1].get('pa', 0.0) + 40.0, q=0.3)
            cands['bar'].append(start_base + [bar2])
    if 'ring' in features:
        rs_ = residual_ring_start(resid, good, cs, re_d)
        if rs_ is not None:
            rr, sg, fl = rs_
            ring = dict(kind='gring', x=cs[0]['x'], y=cs[0]['y'], flux=min(fl, 0.5 * f_tot), rring=rr, sring=sg, q=cs[1]['q'], pa=cs[1]['pa'])
            cands.setdefault('ring', []).append(start_base + [ring])
    if 'bar' in features and 'ring' in features and 'bar' in cands and 'ring' in cands:
        bst, rst = cands['bar'][0][2], cands['ring'][0][2]
        cands['bar+ring'] = [start_base + [bst, dict(rst)]]
    if 'spiral' in features:
        for th in (-200.0, 200.0):
            d = dict(start_base[1])
            sp = MF.spiral_disc(d['x'], d['y'], d['flux'], d['re'], d['q'], d['pa'], theta=th)
            cands.setdefault('spiral', []).append([start_base[0], sp])
    for name, starts in cands.items():
        best = None
        for st in starts:
            try:
                r = MF.fit(data, st, tie=TIE4 if name == 'bar+ring' else TIE3 if name in ('bar', 'ring') else TIE2, max_nfev=max_nfev, **kw)
            except Exception as e:
                notes.append('%s fit failed: %s' % (name, str(e)[:60]))
                continue
            if best is None or r['chi2'] < best['chi2']:
                best = r
        if best is not None:
            allr[name] = best
    chi_min = min(r['chi2'] / max(r['dof'], 1) for r in allr.values())
    scale = max(1.0, chi_min)
    bic = {k: _bic(r, scale) for k, r in allr.items()}
    best_name = 'base'
    for k in ('bar', 'ring', 'bar+ring', 'spiral'):
        if k not in allr:
            continue
        ok, why = sanity(k, allr[k]['components'], re_d)
        if not ok:
            notes.append('%s rejected: %s' % (k, why))
            flag |= FLAG[k.split('+')[0].upper() + '_REJECT']
            continue
        if bic[k] < bic[best_name] - bic_margin:
            best_name = k
    return dict(name=best_name, res=allr[best_name], all=allr, bic=bic, scale=scale, notes=notes, flag=flag, base_swapped=not ok_base)


def bar_summary(res, zp=25.0):
    """Physical summary of the bar component of a '+bar' result: dict(pa, length, q, c0, frac, bt, bar_to_total)."""
    cs = res['components']
    tot = sum(10 ** (-0.4 * c['mag']) for c in cs)
    f = [10 ** (-0.4 * c['mag']) / tot for c in cs]
    b = cs[2]
    return dict(pa=b['pa'], length=b['rout'], q=b['q'], c0=b['c0'], bar_to_total=f[2], bulge_to_total=f[0], disc_to_total=f[1])
