"""Single-component elliptical Sérsic fit.

I = Ie * exp(-bn * ((r/re)**(1/n) - 1)), with bn from the Ciotti & Bertin
approximation. r is the elliptical radius. q = b/a. theta is radians,
counter-clockwise from +x. Coordinates are 0-based.

An optional PSF is normalised to unit sum, clipped at zero, and convolved with
the galaxy. The flat background is added after that convolution. Total flux of
the unconvolved profile is

    Ie * exp(bn) * n * 2π * re**2 * q * Γ(2n) / bn**(2n).

This is one component. It is not a multi-component galaxy fit, not a Fourier
mode, and not a substitute for another profile code.
"""
import math

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve


def bn(n):
    n = np.maximum(np.asarray(n, dtype=np.float64), 0.36)
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n ** 2) + 131.0 / (1148175.0 * n ** 3)


def _clipped(re, n, q):
    return float(max(re, 1e-3)), float(max(n, 0.36)), float(np.clip(q, 0.05, 1.0))


def total_flux(ie, re, n, q):
    """Integral of this elliptical profile over the plane."""
    re, n, q = _clipped(re, n, q)
    b = float(bn(n))
    return float(ie) * math.exp(b) * n * (2.0 * math.pi) * re * re * q * math.gamma(2.0 * n) / (b ** (2.0 * n))


def ie_from_flux(flux, re, n, q):
    """Ie that gives ``flux`` for this profile. Non-positive flux returns a small Ie."""
    re, n, q = _clipped(re, n, q)
    scale = total_flux(1.0, re, n, q)
    if scale <= 0.0 or float(flux) <= 0.0:
        return 1e-3
    return float(flux) / scale


def theta_from_galfit_pa(pa_deg):
    """ogfmeas theta from a feedme position angle (degrees, +y toward -x)."""
    return math.radians(float(pa_deg) + 90.0)


def galfit_pa_from_theta(theta):
    """Feedme position angle in degrees, folded into [-90, 90)."""
    return (math.degrees(float(theta)) % 180.0) - 90.0


def convolve(image, psf):
    """Convolve ``image`` with a 2-D PSF. The kernel is clipped at zero and sums to 1."""
    kernel = np.asarray(psf, dtype=np.float64)
    if kernel.ndim != 2 or kernel.size == 0:
        raise ValueError("PSF must be a 2-D array")
    kernel = np.clip(kernel, 0.0, None)
    total = float(kernel.sum())
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("PSF has no positive flux")
    kernel = kernel / total
    return fftconvolve(np.asarray(image, dtype=np.float64), kernel, mode="same")


def render(shape, x0, y0, re, n, ie, q, theta, background=0.0, psf=None):
    """Evaluate one elliptical Sérsic on an image of ``shape`` (ny, nx)."""
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    galaxy = model(xx, yy, x0, y0, re, n, ie, q, theta, 0.0)
    if psf is not None:
        galaxy = convolve(galaxy, psf)
    return galaxy + float(background)


def model(xx, yy, x0, y0, re, n, ie, q, theta, background):
    q = float(np.clip(q, 0.05, 1.0))
    n = float(max(n, 0.36))
    re = float(max(re, 1e-3))
    dx = xx - x0
    dy = yy - y0
    ct = np.cos(theta)
    st = np.sin(theta)
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    radius = np.sqrt(xr * xr + (yr / q) ** 2)
    return ie * np.exp(-bn(n) * ((radius / re) ** (1.0 / n) - 1.0)) + background


def _clip_start(value, lo, hi):
    if math.isfinite(lo) and value <= lo:
        return lo + 1e-6
    if math.isfinite(hi) and value >= hi:
        return hi - 1e-6
    return value


def fit_sersic(image, x0=None, y0=None, mask=None, psf=None, sigma=None,
               re0=None, n0=None, q0=None, theta0=None, ie0=None, background0=None):
    """Least-squares fit. Returns n, re, ie, q, ellip, theta, x, y, background, chi2.

    ``psf=None`` does not convolve. ``sigma`` weights residuals; None weights
    every good pixel equally. chi2 is the mean square residual after that weight.
    """
    data = np.asarray(image, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError("Sérsic image must be 2-D")
    ny, nx = data.shape
    good = np.isfinite(data)
    if mask is not None:
        good &= ~np.asarray(mask, dtype=bool)
    if not np.any(good):
        raise ValueError("Sérsic fit has no good pixels")
    work = np.where(good, data, 0.0)
    if x0 is None or y0 is None:
        score = np.where(good, work, -np.inf)
        y_peak, x_peak = np.unravel_index(int(np.argmax(score)), score.shape)
        if x0 is None:
            x0 = float(x_peak)
        if y0 is None:
            y0 = float(y_peak)
    x0 = float(x0)
    y0 = float(y0)
    border = np.concatenate([work[0, :], work[-1, :], work[:, 0], work[:, -1]])
    bg_border = float(np.median(border[np.isfinite(border)])) if border.size else 0.0
    peak = float(np.max(work[good])) if good.any() else 0.0
    ie_init = max(peak - bg_border, 1e-3) if ie0 is None else max(float(ie0), 1e-3)
    bg0 = bg_border if background0 is None else float(background0)
    re_init = 4.0 if re0 is None else float(re0)
    n_init = 1.5 if n0 is None else float(n0)
    q_init = 0.85 if q0 is None else float(q0)
    th_init = 0.0 if theta0 is None else float(theta0)
    if theta0 is not None:
        th_init = (th_init + math.pi) % (2.0 * math.pi) - math.pi
    lower = np.array([x0 - 5, y0 - 5, 0.4, 0.4, 0.0, 0.15, -np.pi, -np.inf])
    upper = np.array([x0 + 5, y0 + 5, max(nx, ny), 8.0, np.inf, 1.0, np.pi, np.inf])
    p0 = np.array([x0, y0, re_init, n_init, ie_init, q_init, th_init, bg0], dtype=np.float64)
    for i in range(p0.size):
        p0[i] = _clip_start(float(p0[i]), float(lower[i]), float(upper[i]))
    yy, xx = np.mgrid[0:ny, 0:nx]
    xxg = xx[good]
    yyg = yy[good]
    target = work[good]
    if sigma is None:
        sig = 1.0
    else:
        sig = np.asarray(sigma, dtype=np.float64)
        if sig.shape == data.shape:
            sig = sig[good]
        elif sig.ndim != 0 and sig.shape != target.shape:
            raise ValueError("sigma must be a scalar or an image")
        sig = np.maximum(sig, 1e-8)

    if psf is None:
        def residuals(p):
            return (model(xxg, yyg, *p) - target) / sig

        nfev = 200
    else:
        def residuals(p):
            galaxy = model(xx, yy, p[0], p[1], p[2], p[3], p[4], p[5], p[6], 0.0)
            pred = convolve(galaxy, psf) + p[7]
            return (pred[good] - target) / sig

        nfev = 600

    result = least_squares(residuals, p0, bounds=(lower, upper), method="trf", max_nfev=nfev)
    x, y, re, n, ie, q, theta, background = [float(v) for v in result.x]
    dof = max(int(target.size) - 8, 1)
    chi2 = float(np.sum(result.fun ** 2) / dof)
    return {
        "x": x, "y": y, "re": re, "n": n, "ie": ie, "q": q,
        "ellip": 1.0 - q, "theta": theta, "background": background, "chi2": chi2,
    }
