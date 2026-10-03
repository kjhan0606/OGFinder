"""Advanced GALFIT-3 style model components for ``ogfkit.multifit``: azimuthal shape modifiers, spiral coordinate rotation, truncation and the extra radial profiles.

A component dict (see ``multifit``) may carry, as flat numeric keys (all fit parameters like x, y, re, ...):

  c0                       diskiness (< 0) / boxiness (> 0) of the generalised ellipse  r = (|x|^(c0+2) + |y/q|^(c0+2))^(1/(c0+2))
  f<m>a, f<m>p             Fourier mode m: fractional amplitude and phase [deg] (GALFIT F<m>):  r = r0 (1 + sum_m a_m cos(m (theta + phi_m)))
  b<m>                     bending mode m (GALFIT B<m>):  y' = y + sum_m b_m (x / r_scale)^m   (in the rotated frame, r_scale = R_e, R_s, FWHM-radius ... of the profile)
  rot_func                 'power' | 'log' (string, not fitted);  rot_in, rot_out, rot_theta, rot_alpha (power) or rot_ws (log), rot_incl, rot_pa : coordinate rotation (GALFIT R0 ... R10)
  trunc_in / trunc_out     lists of indices of truncation components (kind 'trunc') that taper the profile inside / outside their break radius

Truncation components (kind 'trunc', parameters x, y, rbreak, dsoft, q, pa and the same modifiers) multiply the profile by the hyperbolic-tangent P(r) (Peng et al. 2010, Appendix B).

The module is deliberately independent of the fitter: ``surface_brightness(c, x, y)`` evaluates one component (unconvolved, counts per pixel) at arbitrary points and
``render(c, shape, truncs)`` integrates it on the pixel grid (sub-pixel sampling where the profile is steep).
"""
import math

import numpy as np

from . import models as M

EXTRA_KINDS = ('moffat', 'ferrer', 'king', 'nuker', 'edgedisk', 'brokenexp', 'gring')
SHAPE_KINDS = ('sersic', 'exp', 'dev', 'moffat', 'ferrer', 'king', 'nuker', 'edgedisk', 'brokenexp', 'gring')
EXP_RE = 1.678346990016661              # b_1 : R_e = EXP_RE R_s for the exponential profile
ROT_KEYS = ('rot_in', 'rot_out', 'rot_theta', 'rot_alpha', 'rot_ws', 'rot_incl', 'rot_pa')
TRUNC_PARAMS = ('x', 'y', 'rbreak', 'dsoft', 'q', 'pa')


# ------------------------------------------------------------------------------------------------------------------ modifier bookkeeping
def modes(c, prefix):
    """Mode numbers present as `f<m>a` (prefix 'f') or `b<m>` (prefix 'b')."""
    out = []
    for k in c:
        if len(k) >= 2 and k[0] == prefix and k[1:].rstrip('ap').isdigit() and (prefix == 'b' and k[1:].isdigit() or prefix == 'f' and k.endswith('a') and k[1:-1].isdigit()):
            out.append(int(k[1:] if prefix == 'b' else k[1:-1]))
    return sorted(set(out))


def modifier_names(c):
    """Names of the optional (shape modifier) parameters present in component c, in a fixed order."""
    names = []
    if 'c0' in c:
        names.append('c0')
    for m in modes(c, 'f'):
        names += ['f%da' % m, 'f%dp' % m]
    for m in modes(c, 'b'):
        names.append('b%d' % m)
    if c.get('rot_func'):
        names += ['rot_in', 'rot_out', 'rot_theta', 'rot_alpha' if c.get('rot_func') == 'power' else 'rot_ws', 'rot_incl', 'rot_pa']
    return names


# ------------------------------------------------------------------------------------------------------------------ coordinate pipeline
def _tanh_rot_angle(r, rin, rout, theta_out):
    """Hyperbolic-tangent part of GALFIT's rotation function (Peng et al. 2010, App. A); angle in degrees."""
    cdef = 20.0                                          # the angle (deg) reached at r_in; the paper prints 0.23, the binary behaves as 20 deg (measured)
    a = min(2.0 * cdef / max(abs(theta_out), 1e-9) - 1.00001, 0.99999)
    b = (2.0 - np.arctanh(a)) * rout / (rout - rin)
    return 0.5 * (np.tanh(b * (r / rout - 1.0) + 2.0) + 1.0)


def rotation_angle(r, c):
    """Spiral rotation angle theta(r) [deg] for circular-centric radius r (GALFIT R-function)."""
    rin, rout, th = c['rot_in'], c['rot_out'], c['rot_theta']
    t = th * _tanh_rot_angle(r, rin, rout, th)
    if c.get('rot_func') == 'log':
        ws = c.get('rot_ws', 1.0)
        t = t * (np.log(r / ws + 1.0) / np.log(rout / ws + 1.0))
    else:
        t = t * (0.5 * (r / rout + 1.0)) ** c.get('rot_alpha', 0.0)
    return t


def shape_radius(dx, dy, c, rscale, pa=None, q=None, dsign=1.0):
    """Generalised radius r (same units as `rscale`) of points (dx, dy) = (x - x0, y - y0) for the shape parameters in c (pa [deg, ccw from +x], q, c0, f<m>, b<m>, rot_*)."""
    pa = c.get('pa', 0.0) if pa is None else pa
    q = c.get('q', 1.0) if q is None else q
    if c.get('rot_func'):
        dx, dy = _apply_rotation(dx, dy, c)
        pa = pa + c.get('rot_pa', 0.0)                       # with a rotation function the model PA is measured from the sky position angle R10 (checked against GALFIT)
    t = math.radians(pa)
    ct, st = math.cos(t), math.sin(t)
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    bm = modes(c, 'b')
    if bm:
        yr = yr + sum((-1.0) ** (m + 1) * c['b%d' % m] * (xr / rscale) ** m for m in bm)       # sign convention measured against GALFIT 3.0.5 (odd modes as written in the paper, even modes opposite)
    p = 2.0 + c.get('c0', 0.0)
    r0 = (np.abs(xr) ** p + np.abs(yr / q) ** p) ** (1.0 / p)
    fm = modes(c, 'f')
    if fm:
        th = np.arctan2(yr / q, -xr)                          # GALFIT measures the phase from the major axis with the x axis mirrored (checked against the binary for m = 1..5)
        g = sum(c['f%da' % m] * np.cos(m * (th + math.radians(c.get('f%dp' % m, 0.0)))) for m in fm)
        r0 = r0 * (1.0 + g)
    return r0


def _apply_rotation(dx, dy, c):
    """Spiral coordinate rotation.  (dx, dy) sky offsets -> coordinates in the plane of the spiral (inclination rot_incl about the axis at sky angle rot_pa; face-on: rot_incl = 0)
    -> rotation by theta(r), r = circular-centric distance in that plane (positive rot_theta winds the image clockwise in GALFIT's frame)."""
    incl = c.get('rot_incl', 0.0)
    if incl:
        beta = math.radians(c.get("rot_pa", 0.0))
        cb, sb = math.cos(beta), math.sin(beta)
        u = dx * cb + dy * sb
        v = (-dx * sb + dy * cb) / max(math.cos(math.radians(incl)), 1e-3)
        dx, dy = u * cb - v * sb, u * sb + v * cb
    r = np.hypot(dx, dy)
    th = np.deg2rad(rotation_angle(r, c))
    ct, st = np.cos(th), np.sin(th)
    return dx * ct + dy * st, -dx * st + dy * ct


# ------------------------------------------------------------------------------------------------------------------ radial profiles
def rscale_of(c):
    k = c['kind']
    if k == 'exp':
        return c['re'] / EXP_RE
    if k in ('sersic', 'dev'):
        return c['re']
    if k == 'moffat':
        return moffat_rd(c['fwhm'], c['beta'])
    return {'ferrer': c.get('rout'), 'king': c.get('rc'), 'nuker': c.get('rb'), 'edgedisk': c.get('rs'), 'brokenexp': c.get('h1'), 'gring': c.get('rring')}[k]


def moffat_rd(fwhm, beta):
    return fwhm / (2.0 * math.sqrt(2.0 ** (1.0 / beta) - 1.0))


def _k1(z):
    from scipy.special import k1
    return k1(z)


def radial(c, r):
    """Surface brightness (counts per pixel) as a function of the generalised radius r for component c (no truncation)."""
    k = c['kind']
    if k in ('sersic', 'exp', 'dev'):
        n = {'exp': 1.0, 'dev': 4.0}.get(k, c.get('n', 1.0))
        re = c['re']
        q = min(max(c.get('q', 1.0), 0.02), 1.0)
        bn = M.sersic_bn(n)
        norm = c.get('norm', 'total')
        if norm == 'total':
            Ie = c['flux'] / (flux_factor(c, 1.0) * sersic_unit_flux(re, n, q))
        elif norm == 'center':                                                   # GALFIT 'sersic1': i0 = surface brightness at r = 0
            Ie = c['i0'] * math.exp(-bn)
        elif norm == 're':                                                       # 'sersic2': i0 = surface brightness at r = R_e
            Ie = c['i0']
        else:                                                                    # 'sersic3': i0 = surface brightness at the break radius of the first truncation
            Ie = c['i0'] / math.exp(-bn * ((c['_rref'] / re) ** (1.0 / n) - 1.0))
        return Ie * np.exp(-bn * ((r / re) ** (1.0 / n) - 1.0))
    if k == 'moffat':
        rd = moffat_rd(c['fwhm'], c['beta'])
        q = min(max(c.get('q', 1.0), 0.02), 1.0)
        I0 = c['flux'] * (c['beta'] - 1.0) * flux_factor(c, 1.0) / (math.pi * rd * rd * q)
        return I0 / (1.0 + (r / rd) ** 2) ** c['beta']
    if k == 'ferrer':
        x = np.clip(r / c['rout'], 0.0, None)
        v = 1.0 - np.where(x < 1.0, x, 1.0) ** (2.0 - c['beta'])
        return np.where(x < 1.0, c['i0'] * np.clip(v, 0, None) ** c['alpha'], 0.0)
    if k == 'king':
        rc, rt, al = c['rc'], c['rt'], c['alpha']
        zt = 1.0 / (1.0 + (rt / rc) ** 2) ** (1.0 / al)
        z = 1.0 / (1.0 + (r / rc) ** 2) ** (1.0 / al)
        v = np.clip(z - zt, 0.0, None) ** al
        return np.where(r < rt, c['i0'] * (1.0 - zt) ** (-al) * v, 0.0)
    if k == 'nuker':
        rb, al, be, ga = c['rb'], c['alpha'], c['beta'], c['gamma']
        x = np.maximum(r / rb, 1e-6)
        return c['ib'] * 2.0 ** ((be - ga) / al) * x ** (-ga) * (1.0 + x ** al) ** ((ga - be) / al)
    if k == 'brokenexp':
        h1, h2, rb, al = c['h1'], c['h2'], c['rbreak'], c.get('alpha', 0.5)
        ex = (1.0 / al) * (1.0 / h1 - 1.0 / h2)
        S = (1.0 + math.exp(-al * rb)) ** (-ex)
        u = np.clip(al * (r - rb), -500, 500)
        return c['i0'] * S * np.exp(-r / h1) * (1.0 + np.exp(u)) ** ex
    if k == 'gring':                                                             # Gaussian ring: I(r) = I0 exp(-(r - r_ring)^2 / (2 s^2)) in the generalised radius, normalised to the total flux
        rr, sg = c['rring'], c['sring']
        q = min(max(c.get('q', 1.0), 0.02), 1.0)
        I0 = c['flux'] / (flux_factor(c, 1.0) * gring_unit_flux(rr, sg, q))
        return I0 * np.exp(-0.5 * ((r - rr) / sg) ** 2)
    raise ValueError('unknown kind %r' % k)


def gring_unit_flux(rring, sring, q):
    """Total flux of the Gaussian ring for I0 = 1: 2 pi q int_0^inf r exp(-(r - r_ring)^2 / (2 s^2)) dr (the generalised-radius shape factor is applied by flux_factor)."""
    from scipy.special import erf
    z = rring / (sring * math.sqrt(2.0))
    return 2.0 * math.pi * q * (sring ** 2 * math.exp(-z * z) + rring * sring * math.sqrt(math.pi / 2.0) * (1.0 + erf(z)))


def ferrers_i0_for_flux(flux, rout, alpha, beta, q, c0=0.0):
    """Central surface brightness of a Ferrers bar carrying `flux` (numerical radial integral; shape modifiers through flux_factor)."""
    x = (np.arange(4000) + 0.5) / 4000.0
    v = np.clip(1.0 - x ** (2.0 - beta), 0.0, None) ** alpha
    area = 2.0 * math.pi * q * rout * rout * float(np.sum(x * v) / 4000.0)
    return flux / (area * flux_factor(dict(c0=c0), 1.0))


def sersic_unit_flux(re, n, q):
    return M.sersic_flux_total(1.0, re, n, q)


_AREA_CACHE = {}


def flux_factor(c, _unused=1.0):
    """1 / R(C0; m): area of the modified shape relative to the ellipse of the same q (the total flux of the profile is F_ellipse / R)."""
    c0 = c.get('c0', 0.0)
    fm = modes(c, 'f')
    incl_f = math.cos(math.radians(c['rot_incl'])) if c.get('rot_func') and c.get('rot_incl') else 1.0       # the total flux is that of the inclined (projected) model
    if not fm and not c0:
        return incl_f
    key = (round(c0, 9), round(incl_f, 9)) + tuple((m, round(c['f%da' % m], 9), round(c.get('f%dp' % m, 0.0), 9)) for m in fm)
    v = _AREA_CACHE.get(key)
    if v is None:
        th = (np.arange(4096) + 0.5) * (2 * np.pi / 4096)
        p = 2.0 + c0
        rho2 = (np.abs(np.cos(th)) ** p + np.abs(np.sin(th)) ** p) ** (-2.0 / p)        # squared Euclidean radius of the unit generalised ellipse at angle th (stretched coordinates)
        thm = np.pi - th                                                                  # the same mirrored phase convention as shape_radius
        g = sum(c['f%da' % m] * np.cos(m * (thm + math.radians(c.get('f%dp' % m, 0.0)))) for m in fm) if fm else 0.0
        v = float(np.mean(rho2 / np.clip(1.0 + g, 1e-3, None) ** 2)) * incl_f
        if len(_AREA_CACHE) > 2048:
            _AREA_CACHE.clear()
        _AREA_CACHE[key] = v
    return v


# ------------------------------------------------------------------------------------------------------------------ truncation
def trunc_p(r, rbreak, dsoft, inner):
    """GALFIT hyperbolic-tangent truncation factor.  Inner: P = 0.5 (tanh((2 - B) r / r_break + B) + 1), B = 2.65 - 4.95 r_break / (r_break - r_soft), r_soft = r_break - dsoft (99 % at r_break, 1 % at r_soft).
    Outer: 1 - P of the same function with the roles swapped (break radius r_break + dsoft, softening radius r_break), i.e. 99 % at r_break and 1 % at r_break + dsoft (measured against GALFIT 3.0.5; the paper prints 4.98, the binary uses 4.95)."""
    if inner:
        rb, rs = rbreak, rbreak - dsoft
    else:
        rb, rs = rbreak + abs(dsoft), rbreak
    B = 2.65 - 4.95 * rb / (rb - rs)
    P = 0.5 * (np.tanh((2.0 - B) * (r / rb) + B) + 1.0)
    return P if inner else 1.0 - P


# ------------------------------------------------------------------------------------------------------------------ evaluation
def surface_brightness(c, x, y, truncs=()):
    """Unconvolved surface brightness (counts per pixel) of component c at points (x, y) [0-based pixel coordinates].  truncs: [(trunc_component, inner_bool), ...]."""
    dx, dy = x - c['x'], y - c['y']
    k = c['kind']
    if c.get('norm') == 'break':
        c = dict(c, _rref=(truncs[0][0]['rbreak'] if truncs else 1.0))
    if k == 'edgedisk':
        t = math.radians(c.get('pa', 0.0))
        xr = dx * math.cos(t) + dy * math.sin(t)
        yr = -dx * math.sin(t) + dy * math.cos(t)
        if modes(c, 'b'):
            yr = yr + sum((-1.0) ** (m + 1) * c['b%d' % m] * (xr / c['rs']) ** m for m in modes(c, 'b'))
        z = np.maximum(np.abs(xr) / c['rs'], 1e-8)
        out = c['i0'] * z * _k1(z) / 1.0 * 1.0 / np.cosh(yr / c['hs']) ** 2
    else:
        r = shape_radius(dx, dy, c, rscale_of(c))
        out = radial(c, r)
    for tc, inner in truncs:
        tx, ty = dx + c['x'] - tc.get('x', c['x']), dy + c['y'] - tc.get('y', c['y'])
        tcc = dict(tc)
        tcc.setdefault('pa', c.get('pa', 0.0))
        tcc.setdefault('q', c.get('q', 1.0))
        rt = shape_radius(tx, ty, tcc, tcc['rbreak'])
        out = out * trunc_p(rt, tc['rbreak'], tc['dsoft'], inner)
    return out


def render(c, shape, truncs=(), nsub=7, origin=(0, 0), edge_nsub=5, edge_rmax=8.0):
    """Pixel-integrated image of component c (centre region sub-sampled; GALFIT integrates the central pixels adaptively too)."""
    ny, nx = shape
    y, x = np.mgrid[:ny, :nx].astype(float)
    y += origin[0]
    x += origin[1]
    img = surface_brightness(c, x, y, truncs)
    steep = c['kind'] in ('sersic', 'dev') and c.get('n', 4.0 if c['kind'] == 'dev' else 1.0) > 2.5 or c['kind'] == 'nuker'
    rad = 3.0
    for rr, ns in ((rad, nsub), (1.5, 40 if steep else 21)):
        sel = (x - c['x']) ** 2 + (y - c['y']) ** 2 < (rr + 1.0) ** 2
        if not sel.any() or (rr == 1.5 and not steep and c['kind'] not in ('sersic', 'dev', 'exp', 'nuker', 'brokenexp') and False):
            continue
        o = (np.arange(ns) + 0.5) / ns - 0.5
        oy_, ox_ = np.meshgrid(o, o, indexing='ij')
        xs, ys = x[sel], y[sel]
        acc = np.zeros(xs.size)
        step = max(1, 200000 // (ns * ns))
        for a0 in range(0, xs.size, step):
            sl = slice(a0, a0 + step)
            acc[sl] = surface_brightness(c, xs[sl, None] + ox_.ravel()[None, :], ys[sl, None] + oy_.ravel()[None, :], truncs).sum(axis=1)
        img[sel] = acc / (ns * ns)
    if truncs and edge_nsub > 1:                                                 # pixels on a truncation edge: integrate over the pixel (GALFIT's tanh edge is smoothed the same way)
        edge = np.zeros(img.shape, bool)
        for tc, inner in truncs:
            tcc = dict(tc)
            tcc.setdefault('pa', c.get('pa', 0.0))
            tcc.setdefault('q', c.get('q', 1.0))
            rt = shape_radius(x - tc.get('x', c['x']), y - tc.get('y', c['y']), tcc, tcc['rbreak'])
            pf = trunc_p(rt, tc['rbreak'], tc['dsoft'], inner)
            edge |= (pf > 1e-3) & (pf < 0.999) & ((x - c['x']) ** 2 + (y - c['y']) ** 2 < edge_rmax ** 2)
        if edge.any():
            ns = edge_nsub
            o = (np.arange(ns) + 0.5) / ns - 0.5
            oy_, ox_ = np.meshgrid(o, o, indexing='ij')
            xs, ys = x[edge], y[edge]
            img[edge] = surface_brightness(c, xs[:, None] + ox_.ravel()[None, :], ys[:, None] + oy_.ravel()[None, :], truncs).mean(axis=1)
    return img
