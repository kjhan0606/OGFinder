"""Spatially varying empirical PSF model (PSFEx-like) - pure functions, JSON-serialisable model.

PSF(x, y) = sum_k phi_k(x, y) * C_k   with a polynomial basis phi_k in the normalised image coordinates and pixel-grid
coefficient images C_k sampled at 1/oversample pixel spacing.  A stamp for a star at sub-pixel offset (dx, dy) is made by
separable cubic-convolution interpolation of the fine grid; every stamp is renormalised to unit sum (flux conservation).

    model = build_psf_model(data, ...)         # select stars, iterate, fit coefficient images by sparse least squares
    model.stamp(x, y, dx, dy)                   # unit-sum PSF at an image position and sub-pixel offset
    diagnostics(model)                          # FWHM / ellipticity / PA maps over the field

Fallbacks (few stars): >= 60 stars -> degree 2, >= 20 -> degree 1, >= 4 -> constant, 1-3 -> Moffat fit to those stars,
0 -> Gaussian of the prior FWHM; model.meta['mode'] says which one was used.
"""
import json
import math
import warnings

import numpy as np
from scipy import ndimage as ndi
from scipy import sparse
from scipy.optimize import least_squares
from scipy.sparse.linalg import lsqr

from . import models as _m


# ----------------------------------------------------------------------------------------------------------------- kernel
def _cubic_w(t):
    """Keys cubic-convolution weights (a = -0.5) for taps -1,0,1,2 at fractional position t in [0,1)."""
    t = np.asarray(t, float)
    t2, t3 = t * t, t * t * t
    return np.stack([-0.5 * t3 + t2 - 0.5 * t, 1.5 * t3 - 2.5 * t2 + 1.0, -1.5 * t3 + 2.0 * t2 + 0.5 * t, 0.5 * t3 - 0.5 * t2], axis=-1)


def _taps(size, nf, s, shift):
    """Fine-grid tap indices (size,4) and weights (size,4) for the stamp rows/cols at star offset `shift` (pixels)."""
    c = (nf - 1) / 2.0
    h = (size - 1) // 2
    u = c + s * (np.arange(-h, h + 1) - shift)
    i0 = np.floor(u).astype(int)
    w = _cubic_w(u - i0)
    idx = i0[:, None] + np.arange(-1, 3)[None, :]
    return idx, w


def poly_terms(degree):
    """Exponent pairs (i, j) in total-degree order: 1, u, v, u^2, uv, v^2, ..."""
    return [(i, d - i) for d in range(degree + 1) for i in range(d, -1, -1)]


class PSFModel:
    def __init__(self, size, oversample, degree, coeffs, xr, yr, meta=None):
        self.size = int(size)
        self.oversample = int(oversample)
        self.degree = int(degree)
        self.coeffs = np.asarray(coeffs, float)          # (K, nf, nf)
        self.xr = (float(xr[0]), float(xr[1]))           # image x range (0-based pixels) mapped to [-1, 1]
        self.yr = (float(yr[0]), float(yr[1]))
        self.meta = dict(meta or {})
        self.nf = self.coeffs.shape[-1]
        self.terms = poly_terms(self.degree)
        if len(self.terms) != self.coeffs.shape[0]:
            raise ValueError('coefficient count does not match the polynomial degree')

    @property
    def shape(self):
        return (self.size, self.size)

    # -- evaluation
    def _uv(self, x, y):
        u = 2.0 * (x - self.xr[0]) / max(self.xr[1] - self.xr[0], 1.0) - 1.0
        v = 2.0 * (y - self.yr[0]) / max(self.yr[1] - self.yr[0], 1.0) - 1.0
        return float(np.clip(u, -1.5, 1.5)), float(np.clip(v, -1.5, 1.5))

    def basis(self, x, y):
        u, v = self._uv(x, y)
        return np.array([u ** i * v ** j for i, j in self.terms])

    def fine(self, x, y):
        return np.tensordot(self.basis(x, y), self.coeffs, axes=(0, 0))

    def stamp(self, x, y, dx=0.0, dy=0.0, size=None):
        """Unit-sum PSF stamp centred on pixel (round(x), round(y)) for a star at offset (dx, dy) from that pixel centre."""
        size = int(size or self.size)
        f = self.fine(x, y)
        iy, wy = _taps(size, self.nf, self.oversample, dy)
        ix, wx = _taps(size, self.nf, self.oversample, dx)
        ok_y = (iy >= 0) & (iy < self.nf)
        ok_x = (ix >= 0) & (ix < self.nf)
        iyc = np.clip(iy, 0, self.nf - 1)
        ixc = np.clip(ix, 0, self.nf - 1)
        g = f[iyc[:, None, :, None], ixc[None, :, None, :]]
        ww = (wy * ok_y)[:, None, :, None] * (wx * ok_x)[None, :, None, :]
        out = (g * ww).sum(axis=(2, 3))
        t = out.sum()
        return out / t if t > 0 else out

    def stamp_at(self, x, y, size=None):
        """Stamp for a star at image position (x, y) (0-based, float): returns (stamp, ix0, iy0) with the stamp origin pixel."""
        ix, iy = int(round(x)), int(round(y))
        s = self.stamp(x, y, x - ix, y - iy, size)
        h = (s.shape[0] - 1) // 2
        return s, ix - h, iy - h

    def fwhm_estimate(self, x=None, y=None):
        x = 0.5 * sum(self.xr) if x is None else x
        y = 0.5 * sum(self.yr) if y is None else y
        return psf_shape(self.stamp(x, y, 0.0, 0.0, None))['fwhm']

    # -- (de)serialisation
    def to_dict(self):
        return dict(kind='ogf_psfmodel', version=1, size=self.size, oversample=self.oversample, degree=self.degree, xr=list(self.xr), yr=list(self.yr),
                    coeffs=self.coeffs.tolist(), meta=_jsonable(self.meta))

    @classmethod
    def from_dict(cls, d):
        return cls(d['size'], d['oversample'], d['degree'], np.asarray(d['coeffs'], float), d['xr'], d['yr'], d.get('meta'))

    def save(self, path):
        """JSON (.json) or FITS cube (anything else; PSFEx-like POLDEG/POLZERO/POLSCAL keywords)."""
        if str(path).endswith('.json'):
            with open(path, 'w') as fh:
                json.dump(self.to_dict(), fh)
            return path
        from astropy.io import fits
        h = fits.PrimaryHDU(self.coeffs.astype(np.float32))
        hd = h.header
        hd['OGFPSF'] = (1, 'OGFinder spatially varying PSF model')
        hd['POLDEG1'] = (self.degree, 'polynomial degree')
        hd['PSF_SAMP'] = (1.0 / self.oversample, 'sampling step in pixels')
        hd['PSFSIZE'] = (self.size, 'stamp size in pixels')
        hd['POLZERO1'] = (0.5 * sum(self.xr), 'x zero point')
        hd['POLSCAL1'] = (0.5 * (self.xr[1] - self.xr[0]), 'x scale')
        hd['POLZERO2'] = (0.5 * sum(self.yr), 'y zero point')
        hd['POLSCAL2'] = (0.5 * (self.yr[1] - self.yr[0]), 'y scale')
        hd['OGFMETA'] = (json.dumps(_jsonable(self.meta))[:60000], 'meta json')
        h.writeto(path, overwrite=True)
        return path


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def load_model(path):
    """Load a model from the JSON or FITS file written by PSFModel.save, or wrap a plain PSF image (constant model)."""
    if str(path).endswith('.json'):
        with open(path) as fh:
            return PSFModel.from_dict(json.load(fh))
    from astropy.io import fits
    with fits.open(path) as hl:
        hd = hl[0].header
        data = np.asarray(hl[0].data, float)
        if 'OGFPSF' in hd:
            deg = int(hd['POLDEG1'])
            s = int(round(1.0 / float(hd['PSF_SAMP'])))
            xr = (hd['POLZERO1'] - hd['POLSCAL1'], hd['POLZERO1'] + hd['POLSCAL1'])
            yr = (hd['POLZERO2'] - hd['POLSCAL2'], hd['POLZERO2'] + hd['POLSCAL2'])
            meta = json.loads(hd['OGFMETA']) if 'OGFMETA' in hd else {}
            return PSFModel(int(hd['PSFSIZE']), s, deg, data, xr, yr, meta)
    return from_image(data)


def from_image(psf):
    """Constant model from a plain PSF image (any shape; made odd/square, unit sum)."""
    p = np.asarray(psf, float)
    n = min(p.shape)
    n -= (n + 1) % 2
    cy, cx = p.shape[0] // 2, p.shape[1] // 2
    h = n // 2
    p = p[cy - h:cy + h + 1, cx - h:cx + h + 1]
    p = np.clip(p, 0, None)
    p = p / p.sum()
    return PSFModel(n, 1, 0, p[None], (0, 1), (0, 1), dict(mode='image'))


def analytic(size, fwhm, kind='moffat', beta=2.5, e=0.0, pa_deg=0.0, xr=(0, 1), yr=(0, 1), meta=None):
    p = _m.moffat_psf(size, fwhm, beta=beta, e=e, pa_deg=pa_deg) if kind == 'moffat' else _m.gaussian_psf(size, fwhm, e=e, pa_deg=pa_deg)
    md = dict(mode='analytic_' + kind, fwhm=fwhm)
    md.update(meta or {})
    return PSFModel(size, 1, 0, (p / p.sum())[None], xr, yr, md)


# ------------------------------------------------------------------------------------------------- analytic profiles (DAOPHOT)
ANALYTIC_KINDS = ('gaussian', 'moffat', 'lorentz', 'penny')


def analytic_profile(kind, x, y, fx, fy, pa_deg, beta=2.5, w=0.3):
    """Un-normalised profile at offsets (x, y) from the centre.  fx, fy = FWHM along the two principal axes, pa_deg = angle of
    the fx axis (CCW from +x).  kinds: gaussian, moffat (beta), lorentz (= moffat beta 1), penny (Gaussian core + Lorentzian wing, weight w)."""
    c, sn = math.cos(math.radians(pa_deg)), math.sin(math.radians(pa_deg))
    u = x * c + y * sn
    v = -x * sn + y * c
    if kind == 'gaussian':
        return np.exp(-0.5 * ((u / (fx / 2.3548)) ** 2 + (v / (fy / 2.3548)) ** 2))
    if kind in ('moffat', 'lorentz'):
        b = 1.0 if kind == 'lorentz' else max(beta, 1.05)
        f = 2.0 * math.sqrt(2.0 ** (1.0 / b) - 1.0)
        return (1.0 + (u / (fx / f)) ** 2 + (v / (fy / f)) ** 2) ** (-b)
    if kind == 'penny':
        g = np.exp(-0.5 * ((u / (fx / 2.3548)) ** 2 + (v / (fy / 2.3548)) ** 2))
        f = 2.0 * math.sqrt(2.0 ** (1.0 / 1.5) - 1.0)
        l = (1.0 + (u / (fx / f)) ** 2 + (v / (fy / f)) ** 2) ** (-1.5)
        return (1.0 - w) * g + w * l
    raise ValueError('unknown analytic PSF kind %r' % kind)


def analytic_fine(kind, size, s, fx, fy, pa_deg, beta=2.5, w=0.3, dx=0.0, dy=0.0, nsub=5):
    """Pixel-integrated analytic profile sampled on the model's fine grid (unit sum over the pixel-spaced samples)."""
    nf = (size - 1) * s + 1
    c = (nf - 1) / 2.0
    g = (np.arange(nf) - c) / s
    sub = (np.arange(nsub) - (nsub - 1) / 2.0) / nsub
    acc = np.zeros((nf, nf))
    for oy in sub:
        for ox in sub:
            X, Y = np.meshgrid(g - dx + ox, g - dy + oy)
            acc += analytic_profile(kind, X, Y, fx, fy, pa_deg, beta, w)
    acc /= nsub * nsub
    return acc / (acc.sum() / s ** 2)


def fit_analytic_shape(stack, kind):
    """Fit (fx, fy, pa, beta|w) of an analytic profile to a centred, unit-sum star stack (integer-grid stamp)."""
    m = stack.shape[0]
    st = stack / stack.sum()
    def mod(p):
        fine = analytic_fine(kind, m, 1, abs(p[0]) + 0.3, abs(p[1]) + 0.3, p[2], beta=p[3] if kind == 'moffat' else 2.5, w=float(np.clip(p[3], 0, 0.9)) if kind == 'penny' else 0.3, nsub=3)
        return fine / fine.sum()
    p0 = [3.0, 3.0, 0.0, 2.5 if kind == 'moffat' else 0.2]
    lo = [0.8, 0.8, -90, 1.1 if kind == 'moffat' else 0.0]
    hi = [30, 30, 90, 12 if kind == 'moffat' else 0.9]
    best = None
    for fw0 in (2.5, 4.0):
        for pa0 in (0.0, 45.0):
            q = list(p0)
            q[0] = q[1] = fw0
            q[2] = pa0
            r = least_squares(lambda p: (mod(p) - st).ravel(), q, bounds=(lo, hi), max_nfev=80)
            if best is None or r.cost < best.cost:
                best = r
    p = best.x
    return dict(kind=kind, fx=float(abs(p[0]) + 0.3), fy=float(abs(p[1]) + 0.3), pa=float(p[2]), beta=float(p[3]) if kind == 'moffat' else None,
                w=float(p[3]) if kind == 'penny' else None, cost=float(best.cost))


def analytic_model(shape_p, size, s, xr, yr, meta=None):
    fine = analytic_fine(shape_p['kind'], size, s, shape_p['fx'], shape_p['fy'], shape_p['pa'], beta=shape_p.get('beta') or 2.5, w=shape_p.get('w') or 0.3)
    md = dict(mode='analytic_' + shape_p['kind'], shape=shape_p)
    md.update(meta or {})
    return PSFModel(size, s, 0, fine[None], xr, yr, md)


def growth(model, x=None, y=None, radii=None):
    """Enclosed flux fraction of the model PSF inside circular apertures (pixel-centre sum on a 4x sub-sampled stamp)."""
    x = 0.5 * sum(model.xr) if x is None else x
    y = 0.5 * sum(model.yr) if y is None else y
    big = model.size
    fine = model.fine(x, y)
    s = model.oversample
    c = (model.nf - 1) / 2.0
    yy, xx = np.indices(fine.shape)
    rr = np.hypot(xx - c, yy - c) / s
    radii = np.asarray(radii if radii is not None else np.arange(1, big // 2 + 1), float)
    tot = fine.sum()
    return radii.tolist(), [float(fine[rr <= r].sum() / tot) for r in radii]


# ------------------------------------------------------------------------------------------------------------- shape stats
def psf_shape(stamp, oversample=1):
    """FWHM (half-maximum area equivalent), ellipticity 1-b/a and PA (deg CCW from +x) from windowed moments of a PSF stamp."""
    p = np.asarray(stamp, float)
    p = np.clip(p, 0, None)
    n = p.shape[0]
    yy, xx = np.indices(p.shape)
    pk = p.max()
    if pk <= 0:
        return dict(fwhm=float('nan'), e=float('nan'), pa=float('nan'))
    # half-max area on a 4x upsampled copy
    z = ndi.zoom(p, 4, order=3) if n < 80 else p
    zs = 4.0 if n < 80 else 1.0
    area = (z >= 0.5 * z.max()).sum() / zs ** 2
    fwhm0 = 2.0 * math.sqrt(area / math.pi)
    # windowed moments with a Gaussian window of sigma ~ 1.2 FWHM/2.355 (iterate the centroid)
    cx = cy = (n - 1) / 2.0
    sig = max(fwhm0, 1.5) / 2.355
    for _ in range(4):
        w = np.exp(-0.5 * (((xx - cx) ** 2 + (yy - cy) ** 2) / (1.5 * sig) ** 2)) * p
        sw = w.sum()
        if sw <= 0:
            break
        cx, cy = (w * xx).sum() / sw, (w * yy).sum() / sw
    w = np.exp(-0.5 * (((xx - cx) ** 2 + (yy - cy) ** 2) / (1.5 * sig) ** 2)) * p
    sw = w.sum()
    mxx = (w * (xx - cx) ** 2).sum() / sw
    myy = (w * (yy - cy) ** 2).sum() / sw
    mxy = (w * (xx - cx) * (yy - cy)).sum() / sw
    tr, det = mxx + myy, mxx * myy - mxy ** 2
    d = math.sqrt(max(tr * tr / 4 - det, 0.0))
    l1, l2 = tr / 2 + d, max(tr / 2 - d, 1e-12)
    pa = 0.5 * math.degrees(math.atan2(2 * mxy, mxx - myy))
    return dict(fwhm=float(fwhm0 / oversample), e=float(1.0 - math.sqrt(l2 / l1)), pa=float(pa % 180.0))


def diagnostics(model, grid=5, shape=None):
    """FWHM / ellipticity / PA maps of the model on a grid x grid lattice over the field (JSON-serialisable)."""
    xs = np.linspace(model.xr[0], model.xr[1], grid)
    ys = np.linspace(model.yr[0], model.yr[1], grid)
    out = dict(x=xs.tolist(), y=ys.tolist(), fwhm=[], e=[], pa=[])
    for y in ys:
        r = [psf_shape(model.stamp(x, y, 0, 0)) for x in xs]
        out['fwhm'].append([q['fwhm'] for q in r])
        out['e'].append([q['e'] for q in r])
        out['pa'].append([q['pa'] for q in r])
    for k in ('fwhm', 'e'):
        a = np.array(out[k])
        out[k + '_min'], out[k + '_max'], out[k + '_mean'] = float(a.min()), float(a.max()), float(a.mean())
    return out


# ------------------------------------------------------------------------------------------------------------ star handling
def background(data, mask=None, bw=64):
    """(bkg, rms) maps from sep.Background (mask: True = bad)."""
    import sep
    d = np.ascontiguousarray(np.nan_to_num(data, nan=0.0).astype(np.float32))
    bad = ~np.isfinite(data)
    if mask is not None:
        bad = bad | mask
    bw = int(min(bw, max(16, min(d.shape) // 4)))
    bk = sep.Background(d, mask=bad, bw=bw, bh=bw, fw=3, fh=3)
    return np.asarray(bk.back(), np.float32), np.asarray(bk.rms(), np.float32)


def find_stars(data, bkg, rms, thresh=5.0, minarea=5, mask=None):
    """Candidate point sources: dict of arrays x, y (0-based), peak, snr, flux, fr50 (half-light radius), a, b, flag."""
    import sep
    d = np.ascontiguousarray((np.nan_to_num(data, nan=0.0) - bkg).astype(np.float32))
    bad = ~np.isfinite(data)
    if mask is not None:
        bad = bad | mask
    o = sep.extract(d, thresh, err=np.ascontiguousarray(rms), minarea=minarea, mask=bad, deblend_cont=0.005)
    if len(o) == 0:
        return dict(x=np.zeros(0), y=np.zeros(0), peak=np.zeros(0), snr=np.zeros(0), flux=np.zeros(0), fr50=np.zeros(0), a=np.zeros(0), b=np.zeros(0), flag=np.zeros(0, int))
    rmax = 3.0 * np.sqrt(o['a'] ** 2 + 1)
    fr, _ = sep.flux_radius(d, o['x'], o['y'], np.maximum(rmax, 4.0), 0.5, subpix=5)
    rr = np.asarray(rms)[np.clip(o['y'].astype(int), 0, d.shape[0] - 1), np.clip(o['x'].astype(int), 0, d.shape[1] - 1)]
    return dict(x=o['x'].astype(float), y=o['y'].astype(float), peak=o['peak'].astype(float), snr=(o['peak'] / rr).astype(float),
                flux=o['flux'].astype(float), fr50=np.asarray(fr, float), a=o['a'].astype(float), b=o['b'].astype(float), flag=o['flag'].astype(int))


def select_psf_stars(cand, shape, half, isolation, snr_range=(25.0, 1e9), saturation=None, max_stars=400, locus_tol=0.18):
    """Pick isolated, unsaturated stars on the stellar locus (half-light radius) from `find_stars` candidates.

    Returns (index array, info dict).  The locus is the median half-light radius of the brightest unsaturated candidates;
    members lie within +-max(locus_tol * median, 3 MAD)."""
    n = len(cand['x'])
    if n == 0:
        return np.zeros(0, int), dict(n_cand=0, n_locus=0, fr50_median=float('nan'))
    x, y = cand['x'], cand['y']
    ok = (cand['snr'] >= snr_range[0]) & (cand['snr'] <= snr_range[1]) & (cand['flag'] == 0) & np.isfinite(cand['fr50']) & (cand['fr50'] > 0.3)
    ok &= (x >= half + 1) & (y >= half + 1) & (x < shape[1] - half - 1) & (y < shape[0] - half - 1)
    if saturation is not None:
        ok &= cand['peak'] < saturation
    # isolation
    if n > 1:
        from scipy.spatial import cKDTree
        tr = cKDTree(np.c_[x, y])
        nn, _ = tr.query(np.c_[x, y], k=2)
        ok &= nn[:, 1] > isolation
    sel = np.where(ok)[0]
    if sel.size == 0:
        return sel, dict(n_cand=int(n), n_locus=0, fr50_median=float('nan'))
    order = sel[np.argsort(-cand['snr'][sel])]
    ref = order[:max(5, min(len(order), 40))]
    med = float(np.median(cand['fr50'][ref]))
    for _ in range(3):                                   # iterate the locus centre on the members
        mad = 1.4826 * float(np.median(np.abs(cand['fr50'][sel] - med)))
        tol = max(locus_tol * med, 3 * mad)
        mem = sel[np.abs(cand['fr50'][sel] - med) <= tol]
        if mem.size >= 3:
            med = float(np.median(cand['fr50'][mem]))
    mem = mem[np.argsort(-cand['snr'][mem])][:max_stars]
    return mem, dict(n_cand=int(n), n_isolated_snr=int(sel.size), n_locus=int(mem.size), fr50_median=med, locus_tol=float(tol))


def cut_stamp(img, x, y, half, fill=np.nan):
    ix, iy = int(round(x)), int(round(y))
    out = np.full((2 * half + 1, 2 * half + 1), fill, float)
    y0, x0 = iy - half, ix - half
    ya, xa = max(0, y0), max(0, x0)
    yb, xb = min(img.shape[0], y0 + 2 * half + 1), min(img.shape[1], x0 + 2 * half + 1)
    if yb > ya and xb > xa:
        out[ya - y0:yb - y0, xa - x0:xb - x0] = img[ya:yb, xa:xb]
    return out, ix, iy


# --------------------------------------------------------------------------------------------------------------- fitting
def _star_fit(stamp, w, model, x, y, ix, iy, f0, max_shift=1.5):
    """Fit (flux, dx, dy) of one star stamp (sky-subtracted) with the current model; returns (flux, dx, dy, chi2_red, ok)."""
    m = w > 0
    if m.sum() < 12:
        return f0, 0.0, 0.0, float('nan'), False
    def res(p):
        f, dx, dy = p
        return ((stamp - f * model.stamp(x, y, dx, dy, stamp.shape[0])) * w)[m]
    try:
        r = least_squares(res, [f0, x - ix, y - iy], bounds=([0, x - ix - max_shift, y - iy - max_shift], [np.inf, x - ix + max_shift, y - iy + max_shift]),
                          x_scale=[max(abs(f0), 1.0), 0.3, 0.3], max_nfev=40)
    except Exception:
        return f0, 0.0, 0.0, float('nan'), False
    chi2 = float((r.fun ** 2).sum() / max(m.sum() - 3, 1))
    return float(r.x[0]), float(r.x[1]), float(r.x[2]), chi2, True


def fit_coefficients(stamps, weights, offsets, fluxes, pos, shape_xy, size, oversample, degree, smooth=0.02, ridge_edge=0.5):
    """Sparse weighted least squares for the K coefficient images from sky-subtracted star stamps.

    stamps/weights (N, m, m); offsets (N, 2) = (dx, dy) of each star w.r.t. its centre pixel; fluxes (N,); pos (N, 2) image positions."""
    N, m, _ = stamps.shape
    s = oversample
    nf = (m - 1) * s + 1                     # odd, so that index (nf-1)/2 is the stamp centre
    terms = poly_terms(degree)
    K = len(terms)
    model_tmp = PSFModel(m, s, degree, np.zeros((K, nf, nf)), shape_xy[0], shape_xy[1])
    rows, cols, vals, rhs = [], [], [], []
    r0 = 0
    for i in range(N):
        wi = weights[i]
        sel = (wi > 0).ravel()
        if sel.sum() < 10 or fluxes[i] <= 0:
            continue
        iy, wy = _taps(m, nf, s, offsets[i][1])
        ix, wx = _taps(m, nf, s, offsets[i][0])
        okx = (ix >= 0) & (ix < nf)
        oky = (iy >= 0) & (iy < nf)
        W = (wy * oky)[:, None, :, None] * (wx * okx)[None, :, None, :]                 # (m, m, 4, 4)
        G = np.clip(iy, 0, nf - 1)[:, None, :, None] * nf + np.clip(ix, 0, nf - 1)[None, :, None, :]
        G = np.broadcast_to(G, W.shape).reshape(m * m, 16)
        W = W.reshape(m * m, 16)
        pix = np.where(sel)[0]
        phi = model_tmp.basis(pos[i][0], pos[i][1])
        for k in range(K):
            a = (fluxes[i] * phi[k] * wi.ravel()[pix])[:, None] * W[pix]
            rows.append(np.repeat(np.arange(r0, r0 + pix.size), 16))
            cols.append((G[pix] + k * nf * nf).ravel())
            vals.append(a.ravel())
        rhs.append((stamps[i].ravel()[pix] * wi.ravel()[pix]))
        r0 += pix.size
    if r0 == 0:
        return None
    # regularisation: Laplacian smoothness (on every coefficient image) + weak pull to zero near the stamp edge
    sig = float(np.median(np.concatenate([w[w > 0] for w in weights if (w > 0).any()])))
    sc = smooth * sig * float(np.median(fluxes[fluxes > 0]))
    yy, xx = np.indices((nf, nf))
    rad = np.hypot(yy - (nf - 1) / 2.0, xx - (nf - 1) / 2.0) / ((nf - 1) / 2.0)
    for k in range(K):
        lk = 1.0 if k == 0 else 3.0
        # second differences in x and y
        for axis in (0, 1):
            ii = np.arange(nf * nf).reshape(nf, nf)
            if axis == 0:
                a, b, c = ii[:-2, :], ii[1:-1, :], ii[2:, :]
            else:
                a, b, c = ii[:, :-2], ii[:, 1:-1], ii[:, 2:]
            n = a.size
            rr = np.arange(r0, r0 + n)
            for cc, vv in ((a, 1.0), (b, -2.0), (c, 1.0)):
                rows.append(rr)
                cols.append(cc.ravel() + k * nf * nf)
                vals.append(np.full(n, vv * sc * lk))
            rhs.append(np.zeros(n))
            r0 += n
        edge = np.where(rad.ravel() > ridge_edge)[0]
        rows.append(np.arange(r0, r0 + edge.size))
        cols.append(edge + k * nf * nf)
        vals.append((sc * 3.0 * (rad.ravel()[edge] - ridge_edge) / (1 - ridge_edge + 1e-9) + sc) * lk * 3)
        rhs.append(np.zeros(edge.size))
        r0 += edge.size
    A = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r0, K * nf * nf))
    b = np.concatenate(rhs)
    sol = lsqr(A, b, atol=1e-9, btol=1e-9, iter_lim=400)[0]
    C = sol.reshape(K, nf, nf)
    return C, nf


def _normalise(C, s):
    tot = C[0].sum() / s ** 2
    return C / tot if tot > 0 else C


def auto_degree(n_stars, requested=None):
    cap = requested if (requested is not None and requested >= 0) else 2
    if n_stars >= 150:
        d = 3
    elif n_stars >= 60:
        d = 2
    elif n_stars >= 20:
        d = 1
    else:
        d = 0
    return min(cap, d)


def fit_moffat_to_stars(stamps, weights):
    """Fallback: Moffat (fwhm, beta, e, pa) fitted to the median-combined stars (used with 1-3 stars)."""
    m = stamps.shape[1]
    stack = np.nanmedian(np.where(weights > 0, stamps / np.maximum(stamps.sum(axis=(1, 2)), 1e-9)[:, None, None], np.nan), axis=0)
    stack = np.nan_to_num(stack, nan=0.0)
    ok = stack > -1
    def res(p):
        mod = _m.moffat_psf(m, max(p[0], 0.8), beta=max(p[1], 1.2), e=float(np.clip(p[2], 0, 0.8)), pa_deg=p[3], nsub=3)
        return (stack / max(stack.sum(), 1e-12) - mod).ravel()
    best = None
    for pa0 in (0.0, 60.0, 120.0):
        r = least_squares(res, [3.0, 3.0, 0.1, pa0], bounds=([0.8, 1.2, 0, -180], [20, 20, 0.8, 360]), max_nfev=60)
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    return dict(fwhm=float(p[0]), beta=float(p[1]), e=float(p[2]), pa=float(p[3]), cost=float(best.cost))


# ----------------------------------------------------------------------------------------------------------------- driver
def build_psf_model(data, bkg=None, rms=None, mask=None, xy=None, fwhm_prior=3.0, size=None, oversample='auto', degree=None, snr_min=25.0,
                    saturation=None, n_iter=2, max_stars=400, neighbour_fit=None, smooth=0.02, reject_sigma=3.5, min_dist=None, kind='empirical', lookup=True):
    """Build a spatially varying PSF model from the field stars of `data`.

    xy        optional (N, 2) 0-based star positions; otherwise stars are detected and picked on the stellar locus.
    neighbour_fit(model, data_sub) -> (xs, ys, fluxes) lets the caller (e.g. the crowded-field code) supply the *other* sources fitted with the
              current model so that their light is subtracted from the star stamps (iterative neighbour subtraction); without it neighbours are masked.
    kind      'empirical' (pixel lookup table only) or an analytic DAOPHOT function ('gaussian', 'moffat', 'lorentz', 'penny'); with an analytic kind
              `lookup=True` adds the spatially varying (order `degree`) empirical residual table on top of the analytic base, as DAOPHOT does.
    Returns (PSFModel, info dict).  info['mode'] in {'empirical', 'moffat_fit', 'gaussian_prior'}; info['stars'] = per-star table."""
    data = np.asarray(data, np.float32)
    if bkg is None or rms is None:
        bkg, rms = background(data, mask)
    sub = np.nan_to_num(data - bkg, nan=0.0)
    ny, nx = data.shape
    fw = float(fwhm_prior)
    size = int(size or max(15, 2 * int(math.ceil(3.2 * fw)) + 1))
    size += (size + 1) % 2
    half = size // 2
    xr, yr = (0.0, float(nx - 1)), (0.0, float(ny - 1))
    cand = find_stars(data, bkg, rms, thresh=5.0, minarea=5, mask=mask)
    info = dict(n_candidates=int(len(cand['x'])), fwhm_prior=fw, size=size)
    if xy is None:
        iso = min_dist or max(2.5 * fw, 0.7 * half)
        idx, sinfo = select_psf_stars(cand, data.shape, half, iso, snr_range=(snr_min, 1e9), saturation=saturation, max_stars=max_stars)
        info.update(sinfo)
        pos = np.c_[cand['x'][idx], cand['y'][idx]]
    else:
        pos = np.asarray(xy, float).reshape(-1, 2)
        idx = None
    nst = len(pos)
    info['n_stars'] = int(nst)
    if nst == 0:
        mdl = analytic(size, fw, 'gaussian', xr=xr, yr=yr, meta=dict(mode='gaussian_prior', n_stars=0))
        info.update(mode='gaussian_prior', stars=[])
        return mdl, info
    # initial stamps
    def make_stamps(P, fluxes=None, offs=None, extra_sub=None):
        S = np.zeros((len(P), size, size))
        Wt = np.zeros((len(P), size, size))
        O = np.zeros((len(P), 2))
        for i, (x, y) in enumerate(P):
            st, ix, iy = cut_stamp(sub - (extra_sub if extra_sub is not None else 0.0), x, y, half)
            rr, _, _ = cut_stamp(rms, x, y, half, fill=np.nan)
            mk = np.zeros_like(st, bool) if mask is None else cut_stamp(mask.astype(float), x, y, half, fill=1.0)[0] > 0.5
            w = np.where(np.isfinite(st) & np.isfinite(rr) & (rr > 0) & ~mk, 1.0 / np.where(rr > 0, rr, 1.0), 0.0)
            # local sky correction: median of the outer ring
            yy, xx = np.indices(st.shape)
            ring = (np.hypot(yy - half, xx - half) > half - 1.5) & (w > 0)
            if ring.sum() > 8:
                st = st - np.median(st[ring])
            S[i], Wt[i] = np.nan_to_num(st), w
            O[i] = (x - ix, y - iy)
        return S, Wt, O
    S, Wt, O = make_stamps(pos)
    # mask neighbours (cand list) in the stamps
    if len(cand['x']):
        from scipy.spatial import cKDTree
        tr = cKDTree(np.c_[cand['x'], cand['y']])
        rad = 0.9 * fw + 1.0
        yy, xx = np.indices((size, size))
        for i, (x, y) in enumerate(pos):
            for j in tr.query_ball_point([x, y], half * 1.5):
                if abs(cand['x'][j] - x) < 0.8 and abs(cand['y'][j] - y) < 0.8:
                    continue
                # fainter and brighter neighbours are masked with a radius that grows with their brightness
                rj = rad * (1.0 + 0.25 * max(0.0, math.log10(max(cand['snr'][j], 1.0) / 5.0)))
                ix, iy = int(round(x)), int(round(y))
                m = np.hypot(xx - (cand['x'][j] - ix + half), yy - (cand['y'][j] - iy + half)) < rj
                Wt[i][m] = 0.0
    F = np.array([max(float((S[i] * (Wt[i] > 0)).sum()), 1.0) for i in range(len(S))])
    # fallback: tiny samples
    if nst < 4:
        mfit = fit_moffat_to_stars(S, Wt) if nst >= 1 else None
        mdl = analytic(size, mfit['fwhm'], 'moffat', beta=mfit['beta'], e=mfit['e'], pa_deg=mfit['pa'], xr=xr, yr=yr, meta=dict(mode='moffat_fit', n_stars=int(nst)))
        info.update(mode='moffat_fit', moffat=mfit, stars=[])
        return mdl, info
    deg = auto_degree(nst, degree)
    if oversample == 'auto':
        s = 2 if nst >= 25 else 1       # s=1 + cubic interpolation sharpens the model by ~3 % in FWHM; s=2 is accurate to ~1 %
    else:
        s = int(oversample)
    info.update(degree=deg, oversample=s)
    mdl = None
    keep = np.ones(nst, bool)
    stars_out = []
    base = None
    if kind != 'empirical':
        stack = np.nanmedian(np.array([ndi.shift(S[i] / max(F[i], 1e-9), (-O[i][1], -O[i][0]), order=3) for i in range(nst) if (Wt[i] > 0).mean() > 0.7]), axis=0)
        base = analytic_model(fit_analytic_shape(np.clip(stack, 0, None), kind), size, s, xr, yr, dict(n_stars=int(nst)))
        info['base'] = base.meta['shape']
    for it in range(max(1, n_iter)):
        k = np.where(keep)[0]
        if base is not None and not lookup:
            mdl = base
            mdl.meta.update(iteration=it)
        else:
            tgt = S[k]
            if base is not None:
                tgt = np.array([S[i] - F[i] * base.stamp(pos[i][0], pos[i][1], O[i][0], O[i][1]) for i in k])
            res = fit_coefficients(tgt, Wt[k], O[k], F[k], pos[k], (xr, yr), size, s, deg, smooth=smooth)
            if res is None:
                break
            C, nf = res
            if base is not None:
                C = C.copy()
                C[0] = C[0] + base.coeffs[0]
            mdl = PSFModel(size, s, deg, _normalise(C, s), xr, yr, dict(mode='empirical' if base is None else 'analytic_' + kind + '+lookup', n_stars=int(keep.sum()), iteration=it))
        # re-fit each star's flux / offset and reject outliers
        chis = np.full(nst, np.nan)
        for i in range(nst):
            if not keep[i] and it > 0:
                pass
            f, dx, dy, c2, ok = _star_fit(S[i], Wt[i], mdl, pos[i][0], pos[i][1], int(round(pos[i][0])), int(round(pos[i][1])), F[i])
            if ok:
                F[i] = f
                O[i] = (dx, dy)
                chis[i] = c2
        good = np.isfinite(chis)
        if good.sum() >= 4:
            med = np.nanmedian(chis[good])
            mad = 1.4826 * np.nanmedian(np.abs(chis[good] - med)) + 1e-9
            keep = good & (chis < med + reject_sigma * max(mad, 0.1 * med))
        if neighbour_fit is not None and it < n_iter - 1:
            xs, ys, fs = neighbour_fit(mdl, sub)
            model_img = np.zeros_like(sub)
            for x, y, f in zip(xs, ys, fs):
                st, x0, y0 = mdl.stamp_at(x, y)
                _paste(model_img, st * f, x0, y0)
            # subtract everything except the target star itself
            S2 = np.zeros_like(S)
            for i, (x, y) in enumerate(pos):
                own = np.zeros_like(sub)
                st, x0, y0 = mdl.stamp_at(x, y)
                tgt = _nearest(xs, ys, x, y)
                fo = fs[tgt] if tgt is not None else 0.0
                _paste(own, st * fo, x0, y0)
                s_i, w_i, _ = make_stamps(pos[i:i + 1], extra_sub=model_img - own)
                S2[i] = s_i[0]
                Wt[i] = w_i[0]         # neighbours are subtracted now, so their pixels are used again
            S = S2
    if mdl is None:
        mdl = analytic(size, fw, 'gaussian', xr=xr, yr=yr, meta=dict(mode='gaussian_prior', n_stars=int(nst)))
        info.update(mode='gaussian_prior', stars=[])
        return mdl, info
    # per-star table (positions after refit) + residual statistics
    for i in range(nst):
        sh = psf_shape(S[i] / max(F[i], 1e-9))
        mo = mdl.stamp(pos[i][0], pos[i][1], O[i][0], O[i][1])
        r = (S[i] - F[i] * mo)
        stars_out.append(dict(x=float(pos[i][0]), y=float(pos[i][1]), flux=float(F[i]), fwhm=sh['fwhm'], e=sh['e'], pa=sh['pa'],
                              used=bool(keep[i]), resid_rms_frac=float(np.sqrt(np.mean(r[Wt[i] > 0] ** 2)) / max(F[i] * mo.max(), 1e-9)) if (Wt[i] > 0).any() else float('nan')))
    mdl.meta.update(n_stars=int(keep.sum()), degree=deg, oversample=s, fwhm_center=mdl.fwhm_estimate())
    info.update(mode='empirical', stars=stars_out, n_used=int(keep.sum()), fwhm_center=mdl.meta['fwhm_center'])
    return mdl, info


def _nearest(xs, ys, x, y, tol=1.5):
    if len(xs) == 0:
        return None
    d = np.hypot(np.asarray(xs) - x, np.asarray(ys) - y)
    j = int(np.argmin(d))
    return j if d[j] < tol else None


def _paste(img, stamp, x0, y0):
    h, w = stamp.shape
    ya, xa = max(0, y0), max(0, x0)
    yb, xb = min(img.shape[0], y0 + h), min(img.shape[1], x0 + w)
    if yb > ya and xb > xa:
        img[ya:yb, xa:xb] += stamp[ya - y0:yb - y0, xa - x0:xb - x0]
