"""Synthetic test images with known truth (used by the plugin tests and validation scripts; no randomness without a seed)."""
import numpy as np

from . import models as _m


def psf_params(x, y, shape, fwhm0=2.8, fwhm1=3.8, e0=0.04, e1=0.14, pa=30.0, beta=3.0):
    """Linearly varying Moffat PSF parameters: FWHM grows with x, ellipticity with y."""
    fx = x / max(shape[1] - 1, 1)
    fy = y / max(shape[0] - 1, 1)
    return dict(fwhm=fwhm0 + (fwhm1 - fwhm0) * fx, e=e0 + (e1 - e0) * fy, pa_deg=pa, beta=beta)


def star_stamp(x, y, shape, size=31, **kw):
    """True pixel-integrated PSF (unit sum) of a star at (x, y), returned with its origin: (stamp, x0, y0)."""
    p = psf_params(x, y, shape, **kw)
    ix, iy = int(round(x)), int(round(y))
    h = size // 2
    st = _m.moffat_psf(size, p['fwhm'], beta=p['beta'], e=p['e'], pa_deg=p['pa_deg'], nsub=5, dx=x - ix, dy=y - iy)
    return st, ix - h, iy - h


def add_stamp(img, st, x0, y0, flux=1.0):
    h, w = st.shape
    ya, xa = max(0, y0), max(0, x0)
    yb, xb = min(img.shape[0], y0 + h), min(img.shape[1], x0 + w)
    if yb > ya and xb > xa:
        img[ya:yb, xa:xb] += flux * st[ya - y0:yb - y0, xa - x0:xb - x0]


def star_field(shape=(500, 500), n=400, mag_range=(16.5, 24.0), slope=0.45, zp=25.0, sky=50.0, noise=1.0, gain=None, seed=1, min_sep=0.0, **psf_kw):
    """Crowded star field with a power-law luminosity function dN/dm ~ 10^(slope m) and a spatially varying Moffat PSF.

    Returns (image float32, truth dict with x, y (0-based), mag, flux)."""
    rng = np.random.default_rng(seed)
    lo, hi = mag_range
    u = rng.random(n)
    a = slope * np.log(10)
    mag = lo + np.log(1 + u * (np.exp(a * (hi - lo)) - 1)) / a
    x = rng.uniform(8, shape[1] - 9, n)
    y = rng.uniform(8, shape[0] - 9, n)
    if min_sep > 0:
        keep = np.ones(n, bool)
        for i in range(n):
            d = np.hypot(x[:i] - x[i], y[:i] - y[i])
            if np.any(d[keep[:i]] < min_sep):
                keep[i] = False
        x, y, mag = x[keep], y[keep], mag[keep]
    flux = 10 ** (-0.4 * (mag - zp))
    img = np.zeros(shape, float)
    for xi, yi, fi in zip(x, y, flux):
        st, x0, y0 = star_stamp(xi, yi, shape, **psf_kw)
        add_stamp(img, st, x0, y0, fi)
    rng2 = np.random.default_rng(seed + 1000)
    clean = img + sky
    if gain:
        clean = rng2.poisson(np.clip(clean, 0, None) * gain) / gain
    out = clean + rng2.normal(0, noise, shape)
    return out.astype(np.float32), dict(x=x, y=y, mag=mag, flux=flux, zp=zp, sky=sky, noise=noise)


def render_stars(shape, x, y, flux, sky=50.0, noise=1.0, seed=2, **psf_kw):
    """Image with the given stars (true spatially varying Moffat PSF) + sky + Gaussian noise."""
    img = np.zeros(shape, float)
    for xi, yi, fi in zip(x, y, flux):
        if -10 < xi < shape[1] + 10 and -10 < yi < shape[0] + 10:
            st, x0, y0 = star_stamp(xi, yi, shape, **psf_kw)
            add_stamp(img, st, x0, y0, fi)
    rng = np.random.default_rng(seed)
    return (img + sky + rng.normal(0, noise, shape)).astype(np.float32)
