"""Multi-component galaxy fitting (GALFIT-like): PSF-convolved Sersic / exponential / de Vaucouleurs / PSF components with shared sky.

Pure functions, JSON-serialisable inputs and outputs (server-ready):

    comps = [dict(kind='sersic', x=30, y=30, mag=18.2, re=6, n=2, q=0.7, pa=30), dict(kind='psf', x=31, y=30, mag=21)]
    res = fit(cutout, comps, psf=PSFModel | 2-D array | None, rms=2.0, sky='const', zp=25.0)
    res['components'][i]['params'] / ['errors'] / ['mag'] ...; res['chi2_red']; res['model']; res['residual']

Conventions: pixel coordinates are 0-based within the array passed; `pa` is degrees counter-clockwise from +x (the same as `ogfkit.models.render_sersic`),
`q` = b/a, `re` = effective (half-light) radius along the major axis in pixels, magnitudes `mag = zp - 2.5 log10(flux)` with flux = total counts.
Parameters can be fixed (`fixed=['n']`), bounded (`bounds={'n': (0.5, 8)}`) and tied between components (`tie=[('1.x', '0.x'), ('1.y', '0.y')]`, i.e. common centre).
"""
import math

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

from . import models as M, profiles as PF

KINDS = ('sersic', 'exp', 'dev', 'psf') + PF.EXTRA_KINDS + ('trunc',)
PARAMS = {'sersic': ('x', 'y', 'flux', 're', 'n', 'q', 'pa'), 'exp': ('x', 'y', 'flux', 're', 'q', 'pa'),
          'dev': ('x', 'y', 'flux', 're', 'q', 'pa'), 'psf': ('x', 'y', 'flux'),
          'moffat': ('x', 'y', 'flux', 'fwhm', 'beta', 'q', 'pa'), 'ferrer': ('x', 'y', 'i0', 'rout', 'alpha', 'beta', 'q', 'pa'),
          'king': ('x', 'y', 'i0', 'rc', 'rt', 'alpha', 'q', 'pa'), 'nuker': ('x', 'y', 'ib', 'rb', 'alpha', 'beta', 'gamma', 'q', 'pa'),
          'edgedisk': ('x', 'y', 'i0', 'hs', 'rs', 'pa'), 'brokenexp': ('x', 'y', 'i0', 'h1', 'h2', 'rbreak', 'alpha', 'q', 'pa'),
          'gring': ('x', 'y', 'flux', 'rring', 'sring', 'q', 'pa'), 'trunc': ('rbreak', 'dsoft')}
FIXED_N = {'exp': 1.0, 'dev': 4.0}
FIXED_DEFAULT = {'brokenexp': ('alpha',)}                 # parameters that are held fixed unless the component says otherwise
SB_NORM = {'center': 'i0', 're': 'i0', 'break': 'i0'}     # sersic-family flux normalisations other than the total flux ('sersic1/2/3' of GALFIT): the 'flux' parameter becomes the surface brightness `i0` [counts per pixel]
EXTENDED = ('sersic', 'exp', 'dev') + PF.EXTRA_KINDS


def param_names(c):
    """Names of the fit parameters of component dict c: the kind's own, plus the optional shape modifiers (c0, Fourier, bending, rotation) present in c."""
    k = c['kind']
    names = list(PARAMS[k])
    if k == 'trunc':
        names = [n for n in ('x', 'y') if n in c] + names + [n for n in ('q', 'pa') if n in c]
    elif k in ('sersic', 'exp', 'dev') and c.get('norm', 'total') in SB_NORM:
        names = ['i0' if n == 'flux' else n for n in names]
    return names + PF.modifier_names(c)
FLAGS = dict(BOUND=1, NOCONV=2, CHI2=4, NEIGHBOUR=8, EDGE=16, MASKED=32, SINGULAR=64)


# ------------------------------------------------------------------------------------------------------------------ PSF handling
def psf_stamp(psf, x, y, dx=0.0, dy=0.0, size=None):
    """Unit-sum PSF stamp for a source at (x, y) with sub-pixel offset (dx, dy) from the stamp centre pixel."""
    if psf is None:
        return np.ones((1, 1))
    if hasattr(psf, 'stamp'):
        return np.asarray(psf.stamp(x, y, dx, dy, size), float)
    p = np.asarray(psf, float)
    if dx or dy:
        p = M.shifted_psf(p, dx, dy, order=1 if max(abs(dx), abs(dy)) < 1e-9 else 3)
    s = p.sum()
    return p / s if s > 0 else p


def psf_size(psf):
    if psf is None:
        return 1
    return int(psf.size) if hasattr(psf, 'size') and not hasattr(psf, 'shape') or hasattr(psf, 'stamp') else int(np.asarray(psf).shape[0])


def gaussian_psf_array(fwhm, size=None):
    size = size or max(9, int(math.ceil(7 * fwhm)) | 1)
    return M.gaussian_psf(size, fwhm)


# ------------------------------------------------------------------------------------------------------------------ rendering
def _uses_profiles(c):
    return c['kind'] in PF.EXTRA_KINDS or c.get('norm', 'total') != 'total' or bool(PF.modifier_names(c)) or bool(c.get('trunc_in') or c.get('trunc_out'))


def render_component(c, shape, psf=None, pad=None, comps=None):
    """Convolved image of one component (dict with kind, x, y, flux, [re, n, q, pa, shape modifiers]) on `shape`.  `comps` = the full component list when c refers to truncation components."""
    ny, nx = shape
    kind = c['kind']
    if kind == 'trunc':
        return np.zeros(shape)
    if kind == 'psf':
        x, y = c['x'], c['y']
        ix, iy = int(round(x)), int(round(y))
        st = psf_stamp(psf, x, y, x - ix, y - iy) if psf is not None else np.zeros((1, 1))
        if psf is None:
            st = np.ones((1, 1))
        out = np.zeros(shape)
        h = st.shape[0] // 2
        y0, x0 = iy - h, ix - h
        ya, xa = max(0, y0), max(0, x0)
        yb, xb = min(ny, y0 + st.shape[0]), min(nx, x0 + st.shape[1])
        if yb > ya and xb > xa:
            out[ya:yb, xa:xb] = c['flux'] * st[ya - y0:yb - y0, xa - x0:xb - x0]
        return out
    st = None
    if psf is not None:
        st = psf_stamp(psf, c['x'], c['y'])
        if st.size > 1:
            pad = pad if pad is not None else st.shape[0] // 2
        else:
            st = None
    p = st.shape[0] // 2 if st is not None else 0
    if _uses_profiles(c):
        cc = dict(c)
        cc['q'] = min(max(cc.get('q', 1.0), 0.02), 1.0) if 'q' in cc else cc.get('q')
        if cc['q'] is None:
            cc.pop('q')
        if kind in ('exp', 'dev'):
            cc['n'] = FIXED_N[kind]
        truncs = [(comps[j], True) for j in c.get('trunc_in', ())] + [(comps[j], False) for j in c.get('trunc_out', ())] if comps is not None else []
        big = PF.render(cc, (ny + 2 * p, nx + 2 * p), truncs, origin=(-p, -p))
    else:
        n = FIXED_N.get(kind, c.get('n', 1.0))
        q = min(max(c.get('q', 1.0), 0.02), 1.0)
        re = max(c['re'], 0.05)
        Ie = M.sersic_Ie_from_flux(c['flux'], re, n, q)
        if st is None:
            return M.render_sersic(shape, c['x'], c['y'], Ie, re, n, q, c.get('pa', 0.0))
        big = M.render_sersic((ny + 2 * p, nx + 2 * p), c['x'] + p, c['y'] + p, Ie, re, n, q, c.get('pa', 0.0))
        conv = fftconvolve(big, st, mode='same')
        return conv[p:p + ny, p:p + nx]
    if st is None:
        return big
    conv = fftconvolve(big, st, mode='same')
    return conv[p:p + ny, p:p + nx]


def render_model(comps, shape, psf=None, sky=0.0, sky_grad=(0.0, 0.0)):
    out = np.full(shape, float(sky))
    if sky_grad[0] or sky_grad[1]:
        yy, xx = np.mgrid[:shape[0], :shape[1]]
        out += sky_grad[0] * (xx - shape[1] / 2.0) + sky_grad[1] * (yy - shape[0] / 2.0)
    for c in comps:
        out += render_component(c, shape, psf, comps=comps)
    return out


# ------------------------------------------------------------------------------------------------------------------ parameter handling
DEFAULT_BOUNDS = dict(re=(0.3, 500.0), n=(0.3, 8.0), q=(0.05, 1.0), pa=(-360.0, 360.0))
# bounds / scales of the parameters of the advanced components: (kind, name) overrides name
ADV_BOUNDS = {'rring': (0.5, 2000.0), 'sring': (0.3, 1000.0), 'fwhm': (0.3, 500.0), 'rout': (1.0, 2000.0), 'rc': (0.2, 500.0), 'rt': (1.0, 2000.0), 'rb': (0.2, 500.0), 'hs': (0.2, 500.0), 'rs': (0.3, 2000.0), 'h1': (0.3, 2000.0), 'h2': (0.3, 2000.0),
              'rbreak': (0.5, 2000.0), 'dsoft': (0.05, 2000.0), 'gamma': (0.0, 2.5), 'c0': (-1.5, 3.0), 'rot_in': (-500.0, 500.0), 'rot_out': (0.5, 2000.0), 'rot_theta': (-3600.0, 3600.0),
              'rot_alpha': (-3.0, 3.0), 'rot_ws': (0.1, 1000.0), 'rot_incl': (0.0, 85.0), 'rot_pa': (-720.0, 720.0)}
KIND_BOUNDS = {('moffat', 'beta'): (0.6, 20.0), ('ferrer', 'beta'): (-2.0, 1.95), ('ferrer', 'alpha'): (0.02, 10.0), ('king', 'alpha'): (0.5, 10.0), ('nuker', 'alpha'): (0.3, 10.0),
               ('nuker', 'beta'): (0.05, 10.0), ('brokenexp', 'alpha'): (0.02, 5.0)}
FLUXLIKE = ('flux', 'i0', 'ib')


def _normalise(c, zp):
    c = dict(c)
    kind = c.get('kind', 'sersic')
    if kind not in KINDS:
        raise ValueError('unknown component kind %r' % kind)
    c['kind'] = kind
    if kind == 'trunc':
        c.setdefault('dsoft', 5.0)
        c.setdefault('rbreak', 10.0)
    else:
        if c.get('norm', 'total') in SB_NORM:
            if 'i0' not in c:
                c['i0'] = c.get('flux', 10 ** (-0.4 * (c.get('mag', 20.0) - zp)))
            c.setdefault('flux', c['i0'])
        elif kind in ('ferrer', 'king', 'edgedisk', 'brokenexp'):
            c.setdefault('i0', c.get('flux', 10 ** (-0.4 * (c.get('mag', 20.0) - zp))))
        elif kind == 'nuker':
            c.setdefault('ib', c.get('flux', 10 ** (-0.4 * (c.get('mag', 20.0) - zp))))
        elif 'flux' not in c:
            c['flux'] = 10 ** (-0.4 * (c.get('mag', 20.0) - zp))
        c.setdefault('re', 3.0)
        c.setdefault('q', 0.8)
        c.setdefault('pa', 0.0)
        c.setdefault('n', FIXED_N.get(kind, 1.5))
        if kind == 'moffat':
            c.setdefault('fwhm', 3.0); c.setdefault('beta', 2.5)
        if kind == 'gring':
            c.setdefault('rring', 8.0); c.setdefault('sring', 2.0)
        if kind == 'brokenexp':
            c.setdefault('alpha', 0.5); c.setdefault('h1', 5.0); c.setdefault('h2', 2.0); c.setdefault('rbreak', 10.0)
    c['fixed'] = set(c.get('fixed', ())) | set(FIXED_DEFAULT.get(kind, ())) - set(c.get('free', ()))
    c['bounds'] = dict(c.get('bounds', {}))
    for k in ('trunc_in', 'trunc_out'):
        if k in c:
            c[k] = [int(v) for v in c[k]]
    return c


def _free_list(comps, tie):
    """[(comp index, name)] of the free parameters + map of tied parameters -> (master comp, master name, mode, value) (mode 'eq' | 'offset' | 'ratio')."""
    tie_map = {}
    for t in tie or ():
        a, b = t[0], t[1]
        mode = t[2] if len(t) > 2 and t[2] else 'eq'
        val = t[3] if len(t) > 3 else None
        ia, na = _key(a)
        ib, nb = _key(b)
        if mode in ('offset', 'ratio') and val is None:                       # GALFIT semantics: the offset / ratio is the one of the start values
            va, vb = comps[ia][na], comps[ib][nb]
            val = (va - vb) if mode == 'offset' else (va / vb if vb else 1.0)
        tie_map[(ia, na)] = (ib, nb, mode, val)
    free = []
    for i, c in enumerate(comps):
        for name in param_names(c):
            if name in c['fixed'] or (i, name) in tie_map:
                continue
            free.append((i, name))
    return free, tie_map


def _key(s):
    if isinstance(s, str):
        i, n = s.split('.')
        return int(i), n
    return int(s[0]), s[1]


def _scale(c, name, box):
    if name in ('x', 'y'):
        return 1.0
    if name in FLUXLIKE:
        return max(abs(c[name]), 1e-3)
    if name in ('re', 'fwhm', 'rring', 'sring', 'rout', 'rc', 'rt', 'rb', 'hs', 'rs', 'h1', 'h2', 'rbreak', 'rot_out'):
        return max(c[name], 0.5)
    if name in ('dsoft', 'rot_in'):
        return max(abs(c[name]), 1.0)
    if name in ('n', 'beta', 'alpha', 'gamma', 'rot_alpha'):
        return 1.0
    if name == 'q':
        return 0.3
    if name == 'c0' or name == 'rot_ws' or name.startswith('b') and name[1:].isdigit():
        return 0.3 if name != 'rot_ws' else max(c[name], 1.0)
    if name.startswith('f') and name.endswith('a'):
        return 0.1
    if name in ('rot_theta', 'rot_incl'):
        return 30.0
    return 30.0


def _bounds(c, name, box, pos0):
    if name in c['bounds']:
        return tuple(c['bounds'][name])
    ny, nx = box
    if name == 'x':
        return (max(-0.5, pos0[0] - 0.25 * nx), min(nx - 0.5, pos0[0] + 0.25 * nx))
    if name == 'y':
        return (max(-0.5, pos0[1] - 0.25 * ny), min(ny - 0.5, pos0[1] + 0.25 * ny))
    if name in FLUXLIKE:
        return (0.0, max(abs(c[name]), 1.0) * 1e3)
    if name == 're':
        return (DEFAULT_BOUNDS['re'][0], min(DEFAULT_BOUNDS['re'][1], 0.9 * max(nx, ny)))
    if (c['kind'], name) in KIND_BOUNDS:
        return KIND_BOUNDS[(c['kind'], name)]
    if name in DEFAULT_BOUNDS:
        return DEFAULT_BOUNDS[name]
    if name in ADV_BOUNDS:
        return ADV_BOUNDS[name]
    if name.startswith('f') and name.endswith('a'):
        return (-0.9, 0.9)
    if name.startswith('f') and name.endswith('p'):
        return (-720.0, 720.0)
    if name.startswith('b') and name[1:].isdigit():
        return (-50.0, 50.0)
    return (-np.inf, np.inf)


def _apply(comps, free, tie_map, vec):
    out = [dict(c) for c in comps]
    for v, (i, name) in zip(vec, free):
        out[i][name] = float(v)
    for _ in range(2):                                                   # ties to tied parameters (chains) settle in two passes
        for (i, name), (j, name2, mode, val) in tie_map.items():
            m = out[j][name2]
            out[i][name] = m if mode == 'eq' else (m + val if mode == 'offset' else m * val)
    return out


# ------------------------------------------------------------------------------------------------------------------ the fit
def fit(data, comps, psf=None, rms=1.0, mask=None, sky='const', sky_value=None, gain=None, zp=25.0, tie=None, max_nfev=200, reweight=1, ftol=1e-8, sky_grad=None, ls_kw=None):
    """Fit `comps` to `data` (2-D array, usually a cutout) with a common sky.

    comps    list of component dicts (see module doc); `fixed` (iterable of names), `bounds` ({name: (lo, hi)}) per component.
    psf      PSFModel (spatially varying: evaluated at each component's start position), 2-D array (unit sum), or None (no convolution).
    rms      sigma of the background noise (scalar or array); with `gain` the Poisson noise of the model is added (var = rms^2 + max(model, 0) / gain).
    mask     bool array, True = ignore pixel.
    sky      'const' (fitted), 'plane' (constant + x/y gradient, fitted) or 'fixed' (uses `sky_value`).
    sky_grad (gx, gy) per pixel about the array centre: start value for sky='plane', fixed value otherwise (default 0, 0).
    tie      [('1.x', '0.x'), ...]: parameter of the first entry is set equal to the second (common centre, common PA, ...).
    Returns a dict: components (params, errors, mag, magerr), sky, sky_err, chi2, dof, chi2_red, bic, flags, nfev, converged, model, residual, weights."""
    data = np.asarray(data, float)
    ny, nx = data.shape
    comps0 = [_normalise(c, zp) for c in comps]
    good = np.isfinite(data)
    if mask is not None:
        good &= ~np.asarray(mask, bool)
    rms_a = np.broadcast_to(np.asarray(rms, float), data.shape)
    free, tie_map = _free_list(comps0, tie)
    pos0 = [(c.get('x', (nx - 1) / 2.0), c.get('y', (ny - 1) / 2.0)) for c in comps0]
    # sky start: median of the pixels furthest from the centre
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.hypot(xx - (nx - 1) / 2.0, yy - (ny - 1) / 2.0)
    outer = good & (rr > 0.8 * rr.max())
    sky0 = float(np.median(data[outer])) if outer.sum() > 10 else float(np.median(data[good]))
    if sky == 'fixed':
        sky0 = float(sky_value if sky_value is not None else 0.0)
    nsky = {'const': 1, 'plane': 3, 'fixed': 0}[sky]
    g_fix = (float(sky_grad[0]), float(sky_grad[1])) if sky_grad is not None else (0.0, 0.0)
    lo, hi, x0, xs = [], [], [], []
    for i, name in free:
        c = comps0[i]
        b = _bounds(c, name, (ny, nx), pos0[i])
        lo.append(b[0]); hi.append(b[1])
        v = c[name]
        x0.append(min(max(v, b[0] + 1e-9 * max(1, abs(b[0]))), b[1] - 1e-9 * max(1, abs(b[1]))))
        xs.append(_scale(c, name, (ny, nx)))
    if nsky:
        lo.append(-np.inf); hi.append(np.inf); x0.append(sky0); xs.append(max(float(np.median(rms_a[good])), 1e-6))
    if nsky == 3:
        for _ in range(2):
            lo.append(-np.inf); hi.append(np.inf); x0.append(g_fix[_]); xs.append(max(float(np.median(rms_a[good])), 1e-6) / max(nx, ny))
    x0 = np.array(x0, float); lo = np.array(lo, float); hi = np.array(hi, float); xs = np.array(xs, float)
    nfree = len(free) + nsky
    psfs = [psf for _ in comps0]

    def unpack(v):
        cs = _apply(comps0, free, tie_map, v[:len(free)])
        if nsky == 0:
            return cs, sky0, g_fix
        s = v[len(free)]
        g = (v[len(free) + 1], v[len(free) + 2]) if nsky == 3 else g_fix
        return cs, s, g

    sig = np.where(good, rms_a, np.inf)

    def resid_fn(v, sigma):
        cs, s, g = unpack(v)
        mod = render_model(cs, (ny, nx), psf, s, g)
        r = (data - mod) / sigma
        r[~good] = 0.0
        return r.ravel()

    sigma = sig.copy()
    if gain:
        base = np.clip(data - sky0, 0, None)
        sigma = np.where(good, np.sqrt(rms_a ** 2 + base / gain), np.inf)
    v = x0
    res = None
    for it in range(1 + (reweight if gain else 0)):
        kw_ls = dict(x_scale=xs, method='trf', max_nfev=max_nfev, ftol=ftol, xtol=ftol, gtol=ftol, jac='3-point')      # '2-point' crawled along flat valleys on real data (GALFIT example: 150 evaluations, not converged; 3-point: 20, converged)
        kw_ls.update(ls_kw or {})
        res = least_squares(resid_fn, v, args=(sigma,), bounds=(lo, hi), **kw_ls)
        v = res.x
        if gain and it < reweight:
            cs, s, g = unpack(v)
            mod = render_model(cs, (ny, nx), psf, s, g)
            sigma = np.where(good, np.sqrt(rms_a ** 2 + np.clip(mod - s, 0, None) / gain), np.inf)
    cs, s, g = unpack(v)
    model = render_model(cs, (ny, nx), psf, s, g)
    npix = int(good.sum())
    dof = max(npix - nfree, 1)
    chi2 = float(np.sum(res.fun ** 2))
    chi2_red = chi2 / dof
    # parameter errors from the Jacobian (weighted); scaled by sqrt(chi2_red) as GALFIT does not - report both
    J = res.jac
    flags = 0
    errs = np.full(nfree, np.nan)
    try:
        cov = np.linalg.inv(J.T @ J)
        errs = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        try:
            cov = np.linalg.pinv(J.T @ J)
            errs = np.sqrt(np.clip(np.diag(cov), 0, None))
            flags |= FLAGS['SINGULAR']
        except Exception:
            flags |= FLAGS['SINGULAR']
    out_c = []
    errmap = {f: errs[k] for k, f in enumerate(free)}
    hit = False
    for i, c in enumerate(cs):
        d = dict(kind=c['kind'])
        er = {}
        names = param_names(c)
        for name in names:
            d[name] = float(c[name])
            if (i, name) in errmap:
                er[name] = float(errmap[(i, name)])
        for k_ in ('norm', 'trunc_in', 'trunc_out', 'rot_func'):
            if k_ in c:
                d[k_] = c[k_]
        if c['kind'] in FIXED_N:
            d['n'] = FIXED_N[c['kind']]
        if 'flux' in names and c['kind'] != 'trunc':
            fl = max(c['flux'], 1e-30)
            d['mag'] = float(zp - 2.5 * math.log10(fl)) if c['flux'] > 0 else float('nan')
            er['mag'] = float(1.0857 * er['flux'] / fl) if 'flux' in er and c['flux'] > 0 else float('nan')
        elif c['kind'] == 'trunc' or c['kind'] == 'psf':
            d['mag'] = float('nan')
        else:                                                              # surface-brightness normalised profile: magnitude of the (unconvolved) model inside the fitted box
            try:
                tot = float(render_component(dict(c, flux=c.get('flux', 1.0)), (ny, nx), None, comps=cs).sum())
            except Exception:
                tot = 0.0
            d['mag'] = float(zp - 2.5 * math.log10(tot)) if tot > 0 else float('nan')
            er['mag'] = float('nan')
        d['errors'] = er
        d['free'] = [n_ for n_ in names if (i, n_) in errmap]
        for (ii, name) in free:
            if ii == i:
                k = free.index((ii, name))
                b = (lo[k], hi[k])
                if np.isfinite(c[name]) and (abs(c[name] - b[0]) < 1e-6 * max(1, abs(b[0])) + 1e-9 or abs(c[name] - b[1]) < 1e-6 * max(1, abs(b[1])) + 1e-9) and name not in FLUXLIKE:
                    d.setdefault('at_bound', []).append(name)
                    hit = True
        out_c.append(d)
    if hit:
        flags |= FLAGS['BOUND']
    if not res.success or res.status == 0:
        flags |= FLAGS['NOCONV']
    if chi2_red > 3.0:
        flags |= FLAGS['CHI2']
    if npix < 0.5 * data.size:
        flags |= FLAGS['MASKED']
    k = len(free) + nsky
    bic = chi2 + k * math.log(max(npix, 2))
    sky_err = float(errs[len(free)]) if nsky else 0.0
    return dict(components=out_c, sky=float(s), sky_err=sky_err, sky_grad=[float(g[0]), float(g[1])], chi2=chi2, dof=int(dof), chi2_red=float(chi2_red), bic=float(bic), npix=npix,
                nfree=int(nfree), flags=int(flags), nfev=int(res.nfev), converged=bool(res.success and res.status > 0), model=model, residual=np.where(good, data - model, 0.0), good=good)


# ------------------------------------------------------------------------------------------------------------------ presets
def fit_multistart(data, comps, restarts=2, **kw):
    """fit() from the given start and `restarts` alternative starts (R_e x 0.6 / x 1.6 [/ x 0.35 / x 2.5], Sersic n x 1.3 / x 0.75 [...] of the free extended components); the lowest chi2 wins.
    Trust-region fits of steep profiles can end in a wide, shallow local minimum (seen on GALFIT's own example); one extra start is cheap compared with a wrong answer.
    The result carries `starts_chi2` (chi2 of every start, the first is the given one) and `start_used` (index of the winner)."""
    base = [dict(c) for c in comps]
    res0 = fit(data, [dict(c) for c in base], **kw)
    best, k_best, chis = res0, 0, [res0['chi2']]
    facs = [(0.6, 1.3), (1.6, 0.75), (0.35, 1.0), (2.5, 1.0)]
    for k in range(1, restarts + 1):
        fr, fn = facs[(k - 1) % len(facs)]
        alt = []
        for c in base:
            d = dict(c)
            if d.get('kind', 'sersic') != 'psf':
                if 're' in d and 're' not in d.get('fixed', ()):
                    d['re'] = max(0.5, d['re'] * fr)
                if d.get('kind', 'sersic') == 'sersic' and 'n' in d and 'n' not in d.get('fixed', ()):
                    d['n'] = min(max(d['n'] * fn, 0.4), 7.0)
            alt.append(d)
        try:
            r = fit(data, alt, **kw)
        except Exception:
            continue
        chis.append(r['chi2'])
        if r['chi2'] < best['chi2'] - 1e-6:
            best, k_best = r, k
    best['starts_chi2'] = [float(v) for v in chis]
    best['start_used'] = k_best
    return best


def preset(name, x, y, flux, re, q=0.8, pa=0.0, n=2.0, zp=25.0, psf_frac=0.1):
    """Starting components for a catalog object.  Names: sersic, exp, dev, psf, bulge+disk, psf+sersic, and the structure presets bar (Ferrers), ring (Gaussian ring), spiral (disc with logarithmic spiral arms), bulge+disk+bar, bulge+disk+ring."""
    f = flux
    if name == 'psf':
        return [dict(kind='psf', x=x, y=y, flux=f)]
    if name in ('sersic', 'exp', 'dev'):
        return [dict(kind=name, x=x, y=y, flux=f, re=re, q=q, pa=pa, n=n)]
    if name == 'bulge+disk':
        return [dict(kind='dev', x=x, y=y, flux=0.3 * f, re=0.5 * re, q=min(1.0, q + 0.2), pa=pa, bounds=dict(q=(0.3, 1.0))),
                dict(kind='exp', x=x, y=y, flux=0.7 * f, re=1.2 * re, q=q, pa=pa)]
    if name == 'psf+sersic':
        return [dict(kind='psf', x=x, y=y, flux=psf_frac * f), dict(kind='sersic', x=x, y=y, flux=(1 - psf_frac) * f, re=re, q=q, pa=pa, n=n)]
    if name == 'bar':
        return [bar_component(x, y, f, re, pa)]
    if name == 'ring':
        return [dict(kind='gring', x=x, y=y, flux=f, rring=re, sring=0.25 * re, q=q, pa=pa)]
    if name == 'spiral':
        return [spiral_disc(x, y, f, re, q, pa)]
    if name == 'bulge+disk+bar':
        return [dict(kind='sersic', x=x, y=y, flux=0.2 * f, re=0.3 * re, n=2.0, q=min(1.0, q + 0.2), pa=pa, bounds=dict(q=(0.3, 1.0))),
                dict(kind='exp', x=x, y=y, flux=0.55 * f, re=1.3 * re, q=q, pa=pa),
                bar_component(x, y, 0.25 * f, 0.6 * re, pa + 30.0)]
    if name == 'bulge+disk+ring':
        return [dict(kind='sersic', x=x, y=y, flux=0.2 * f, re=0.3 * re, n=2.0, q=min(1.0, q + 0.2), pa=pa, bounds=dict(q=(0.3, 1.0))),
                dict(kind='exp', x=x, y=y, flux=0.6 * f, re=1.3 * re, q=q, pa=pa),
                dict(kind='gring', x=x, y=y, flux=0.2 * f, rring=1.2 * re, sring=0.25 * re, q=q, pa=pa)]
    raise ValueError('unknown preset %r' % name)


def bar_component(x, y, flux, length, pa, q=0.3, c0=1.0):
    """Ferrers bar (GALFIT 'ferrer': alpha = 2, beta = 0 held fixed, truncation radius `length`, boxy/discy shape c0 free in [0, 3]); `flux` is the integrated flux."""
    return dict(kind='ferrer', x=x, y=y, i0=PF.ferrers_i0_for_flux(flux, length, 2.0, 0.0, q, c0), rout=length, alpha=2.0, beta=0.0, q=q, pa=pa, c0=c0, fixed=['alpha', 'beta'],
                bounds=dict(c0=(0.0, 3.0), q=(0.08, 0.7)))


def spiral_disc(x, y, flux, re, q, pa, theta=-200.0, ampl=0.15):
    """Exponential disc with two logarithmic spiral arms: Fourier mode m = 2 (amplitude free) and GALFIT's coordinate rotation with the logarithmic winding function
    (rot_theta = rotation reached at rot_out, both free; r_in = 0, face-on spiral plane: inclination handled by q / pa)."""
    return dict(kind='exp', x=x, y=y, flux=flux, re=re, q=q, pa=pa, f2a=ampl, f2p=0.0, rot_func='log', rot_in=0.0, rot_out=3.0 * re, rot_theta=theta, rot_ws=0.5 * re,
                rot_incl=0.0, rot_pa=0.0, fixed=['rot_in', 'rot_ws', 'rot_incl', 'rot_pa'], bounds=dict(f2a=(0.0, 0.9)))


TIES = {'bulge+disk': [('1.x', '0.x'), ('1.y', '0.y')], 'psf+sersic': [('0.x', '1.x'), ('0.y', '1.y')],
        'bulge+disk+bar': [('1.x', '0.x'), ('1.y', '0.y'), ('2.x', '0.x'), ('2.y', '0.y')], 'bulge+disk+ring': [('1.x', '0.x'), ('1.y', '0.y'), ('2.x', '0.x'), ('2.y', '0.y')]}


def fit_preset(data, name, x, y, flux, re, q, pa, extra=(), tie_extra=(), **kw):
    """Fit one of the presets plus fixed-start `extra` (neighbour) components (their parameters are free too)."""
    comps = preset(name, x, y, flux, re, q, pa, zp=kw.get('zp', 25.0)) + list(extra)
    tie = list(TIES.get(name, ())) + list(tie_extra)
    return fit(data, comps, tie=tie, **kw), comps


def auto_select(data, x, y, flux, re, q, pa, models=('psf', 'sersic', 'bulge+disk'), bic_margin=10.0, **kw):
    """Fit several presets and keep the one with the lowest BIC (a more complex model must win by `bic_margin`).  Returns (name, result, all results)."""
    allr = {}
    for m in models:
        try:
            allr[m], _ = fit_preset(data, m, x, y, flux, re, q, pa, **kw)
        except Exception as e:  # a failing preset must not stop the others
            allr[m] = None
    order = [m for m in models if allr.get(m) is not None]
    best = order[0]
    for m in order[1:]:
        if allr[m]['bic'] < allr[best]['bic'] - bic_margin:
            best = m
    return best, allr[best], allr
