"""Analytic image models: Sersic (optionally generalised-ellipse), Gaussian / Moffat PSFs, rendering with sub-pixel integration."""
import numpy as np


_BN_CACHE = {}


def sersic_bn(n):
    """Exact b_n: gamma(2n, b_n) = Gamma(2n)/2 (scipy), Ciotti & Bertin series if scipy is missing."""
    key = float(n)
    v = _BN_CACHE.get(key)
    if v is not None:
        return v
    try:
        from scipy.special import gammaincinv
        v = float(gammaincinv(2.0 * n, 0.5))
    except Exception:
        v = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n * n)
    if len(_BN_CACHE) > 4096:
        _BN_CACHE.clear()
    _BN_CACHE[key] = v
    return v


def sersic_flux_total(Ie, re, n, q):
    """Total flux of a Sersic profile (counts) for effective intensity Ie (per pixel), axis ratio q."""
    from scipy.special import gamma
    bn = sersic_bn(n)
    return float(2.0 * np.pi * n * re * re * Ie * np.exp(bn) * bn ** (-2.0 * n) * gamma(2.0 * n) * q)


def sersic_Ie_from_flux(flux, re, n, q):
    return flux / sersic_flux_total(1.0, re, n, q)


def _sersic_eval(x, y, xc, yc, Ie, re, n, q, pa_deg, c=0.0):
    t = np.deg2rad(pa_deg)
    dx = x - xc
    dy = y - yc
    xr = dx * np.cos(t) + dy * np.sin(t)
    yr = -dx * np.sin(t) + dy * np.cos(t)
    p = 2.0 + c
    r = (np.abs(xr) ** p + np.abs(yr / q) ** p) ** (1.0 / p)
    bn = sersic_bn(n)
    return Ie * np.exp(-bn * ((r / re) ** (1.0 / n) - 1.0))


def render_sersic(shape, xc, yc, Ie, re, n, q, pa_deg, c=0.0, nsub=7, rsub=None, origin=(0, 0)):
    """Sersic model on a pixel grid (0-based pixel centres; origin=(y0,x0) offsets a cutout).  The centre
    (r < rsub, default 3 px or 0.5 re) is integrated on an nsub x nsub sub-pixel grid.  pa_deg is counter-clockwise from +x."""
    ny, nx = shape
    y, x = np.mgrid[:ny, :nx].astype(float)
    y += origin[0]
    x += origin[1]
    img = _sersic_eval(x, y, xc, yc, Ie, re, n, q, pa_deg, c)
    rs = rsub if rsub is not None else max(3.0, 0.5 * re) if n > 1.5 else 3.0
    for rad, ns in ((rs, nsub), (1.5, 41 if n > 2.5 else 21)):       # steep cusps (n >~ 3) need a much finer grid in the last 1.5 px
        if ns <= 1 or (rad == 1.5 and n <= 1.5):
            continue
        sel = (x - xc) ** 2 + (y - yc) ** 2 < (rad + 1.0) ** 2
        if not sel.any():
            continue
        o = (np.arange(ns) + 0.5) / ns - 0.5
        oy_, ox_ = np.meshgrid(o, o, indexing='ij')
        xs, ys = x[sel], y[sel]
        acc = np.zeros(xs.size)
        step = max(1, 200000 // (ns * ns))                       # vectorised over the sub-pixel grid, in chunks of pixels
        for a0 in range(0, xs.size, step):
            sl = slice(a0, a0 + step)
            acc[sl] = _sersic_eval(xs[sl, None] + ox_.ravel()[None, :], ys[sl, None] + oy_.ravel()[None, :], xc, yc, Ie, re, n, q, pa_deg, c).sum(axis=1)
        img[sel] = acc / (ns * ns)
    return img


def gaussian_psf(size, fwhm, e=0.0, pa_deg=0.0):
    """Normalised (sum 1) elliptical Gaussian PSF stamp, size odd, centred on the middle pixel; fwhm = geometric mean FWHM."""
    s = float(fwhm) / 2.3548200450309493
    q = 1.0 - e
    sx = s / np.sqrt(q)
    sy = s * np.sqrt(q)
    h = size // 2
    y, x = np.mgrid[-h:h + 1, -h:h + 1].astype(float)
    t = np.deg2rad(pa_deg)
    xr = x * np.cos(t) + y * np.sin(t)
    yr = -x * np.sin(t) + y * np.cos(t)
    p = np.exp(-0.5 * ((xr / sx) ** 2 + (yr / sy) ** 2))
    return p / p.sum()


def moffat_psf(size, fwhm, beta=2.5, e=0.0, pa_deg=0.0, nsub=5, dx=0.0, dy=0.0):
    """Normalised Moffat PSF stamp (pixel-integrated on an nsub grid), centre offset by (dx, dy) pixels."""
    h = size // 2
    alpha = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    q = 1.0 - e
    t = np.deg2rad(pa_deg)
    o = (np.arange(nsub) + 0.5) / nsub - 0.5
    y0, x0 = np.mgrid[-h:h + 1, -h:h + 1].astype(float)
    acc = np.zeros_like(x0)
    for oy in o:
        for ox in o:
            x = x0 + ox - dx
            y = y0 + oy - dy
            xr = x * np.cos(t) + y * np.sin(t)
            yr = -x * np.sin(t) + y * np.cos(t)
            acc += (1.0 + (xr * xr / q + yr * yr * q) / alpha ** 2) ** (-beta)
    return acc / acc.sum()


def shifted_psf(psf, dx, dy, order=3):
    """psf stamp shifted by (dx, dy) pixels with spline interpolation, renormalised to the original sum."""
    from scipy.ndimage import shift
    s = psf.sum()
    out = shift(psf, (dy, dx), order=order, mode='constant', cval=0.0)
    t = out.sum()
    return out * (s / t) if t > 0 else out


def stamp_size_for(fwhm, mult=7.0):
    n = int(np.ceil(mult * fwhm)) | 1
    return max(n, 15)


def convolve_same(img, psf):
    """FFT convolution of img with a (normalised) psf stamp, same shape."""
    from scipy.signal import fftconvolve
    return fftconvolve(img, psf, mode='same')
