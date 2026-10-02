"""Extended non-parametric morphology: smoothness (clumpiness), Kron radius / flux, growth curves, Petrosian radii for several eta, Petrosian-segmented Gini / M20.

Complements `morphometry.measure` (C, A, Gini, M20, R_PETRO at eta = 0.2); pure functions, numpy / scipy / sep only.

All radii are along the major axis of elliptical apertures of axis ratio q and position angle theta (radians, counter-clockwise from +x, SExtractor convention); q = 1 gives
circular apertures.  `data` must be background subtracted.
Definitions
  growth curve   F(<r): sum inside the elliptical aperture of semi-major axis r (sub-pixel integrated, `sep.sum_ellipse`), masked pixels excluded and replaced by the mean of the rest.
  Petrosian eta  eta(r) = I(r) / <I>(<r),  I(r) = [F(1.25 r) - F(0.8 r)] / [pi q (1.25^2 - 0.8^2) r^2],  <I>(<r) = F(r) / (pi q r^2)    (Bershady et al. 2000; SDSS-like)
  R_P(eta)       first radius where eta falls below the value (linear interpolation of the fine log-spaced grid); F_P = F(< 2 R_P(0.2));  R50 / R90 etc. contain that fraction of F_P.
  concentration  C = 5 log10(R80 / R20) of F_P.
  Kron           r_k = sum(r I) / sum(I) in the ellipse of radius `kron_rmax` (iterated once with the result), flux in kron_scale * r_k (>= `kron_min` px) as SExtractor's MAG_AUTO.
  smoothness     S = sum(|I - I_s| - |B - B_s|) / sum(I) inside 1.5 R_P (Conselice 2003), I_s = I smoothed with a boxcar of width `smooth_frac` * R_P (odd pixels, >= 3), B = noise patch of the same
                 size (Gaussian noise of the measured rms, or a supplied empty-sky patch); the central box (r < width / 2) is excluded from both sums.
  Gini_P, M20_P  Gini and M20 on the pixels inside 1.5 R_P above the surface brightness at R_P (Lotz et al. 2004 segmentation).
"""
import math

import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import CubicSpline
from scipy.special import gammainc, gammaincinv


def _sep():
    try:
        import sep
    except ImportError:
        import sep_pjw as sep
    return sep


# --------------------------------------------------------------------------------------------------------------------- apertures
def growth_curve(data, x, y, q, theta, radii, mask=None, subpix=5):
    """F(<r) for an array of semi-major axes `radii` (pixels)."""
    sep = _sep()
    r = np.atleast_1d(np.asarray(radii, np.float64))
    d = np.ascontiguousarray(data, np.float64)
    n = r.size
    theta = (float(theta) + math.pi / 2) % math.pi - math.pi / 2       # sep wants -pi/2 .. pi/2
    m = None if mask is None else np.ascontiguousarray(mask, np.uint8)
    # sep needs b >= 0.0 and a >= b; very small apertures are fine
    a = np.maximum(r, 1e-3)
    b = np.maximum(a * q, 1e-3)
    flux, err, flag = sep.sum_ellipse(d, np.full(n, x), np.full(n, y), a, b, np.full(n, theta), r=1.0, mask=m, subpix=subpix)
    return np.asarray(flux, float)


def radial_grid(rmax, rmin=0.8, n=160):
    return np.geomspace(rmin, rmax, n)


class Curve:
    """Growth curve on a fine log grid with monotone-safe interpolation F(r), and its inverse."""

    def __init__(self, data, x, y, q, theta, rmax, mask=None, n=160, rmin=0.8):
        self.r = radial_grid(rmax, rmin, n)
        self.F = growth_curve(data, x, y, q, theta, self.r, mask)
        self.q = q
        self.spl = CubicSpline(np.log(self.r), self.F, extrapolate=False)

    def flux(self, r):
        r = np.asarray(r, float)
        out = np.full(r.shape, np.nan)
        ok = (r >= self.r[0]) & (r <= self.r[-1])
        out[ok] = self.spl(np.log(r[ok]))
        return out

    def radius_of(self, f):
        """Radius containing flux f (first crossing, linear interpolation)."""
        F = self.F
        idx = np.where(F >= f)[0]
        if idx.size == 0 or f <= 0:
            return float('nan')
        i = idx[0]
        if i == 0:
            return float(self.r[0])
        t = (f - F[i - 1]) / max(F[i] - F[i - 1], 1e-30)
        return float(self.r[i - 1] + t * (self.r[i] - self.r[i - 1]))


def petrosian_eta(curve):
    """eta(r) on the usable part of the grid (r such that 1.25 r <= rmax)."""
    r = curve.r
    ok = 1.25 * r <= r[-1]
    rr = r[ok]
    F = curve.F[ok]
    Fo = curve.flux(1.25 * rr)
    Fi = curve.flux(0.8 * rr)
    q = curve.q
    ring = (Fo - Fi) / (math.pi * q * (1.25 ** 2 - 0.8 ** 2) * rr ** 2)
    mean = F / (math.pi * q * rr ** 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        eta = np.where(mean > 0, ring / mean, np.nan)
    return rr, eta


def petrosian_radius(rr, eta, level):
    """First radius (beyond the peak region) where eta drops below `level`; linear interpolation; nan when it never does."""
    bad = ~np.isfinite(eta)
    for i in range(1, len(rr)):
        if bad[i] or bad[i - 1]:
            continue
        if eta[i - 1] > level >= eta[i]:
            t = (eta[i - 1] - level) / max(eta[i - 1] - eta[i], 1e-30)
            return float(rr[i - 1] + t * (rr[i] - rr[i - 1]))
    return float('nan')


# --------------------------------------------------------------------------------------------------------------------- Kron
def kron(data, x, y, q, theta, rmax, kron_scale=2.5, kron_min=3.5, mask=None, rmin_fit=1.0):
    """Kron radius and flux.  rmax: semi-major axis of the first-moment ellipse.  Returns (r_kron, flux, rap) with rap the semi-major axis of the flux aperture."""
    sep = _sep()
    d = np.ascontiguousarray(data, np.float64)
    m = None if mask is None else np.ascontiguousarray(mask, np.uint8)
    ny, nx = d.shape
    yy, xx = np.mgrid[:ny, :nx]
    t = theta
    dx, dy = xx - x, yy - y
    u = dx * math.cos(t) + dy * math.sin(t)
    v = -dx * math.sin(t) + dy * math.cos(t)
    rr = np.sqrt(u ** 2 + (v / max(q, 0.05)) ** 2)           # elliptical radius (major axis units)
    good = np.ones(d.shape, bool) if mask is None else ~np.asarray(mask, bool)
    sel = (rr <= rmax) & good & (rr >= 0)
    num = float(np.sum(rr[sel] * d[sel]))
    den = float(np.sum(d[sel]))
    if den <= 0:
        return float('nan'), float('nan'), float('nan')
    rk = num / den
    # iterate once with the aperture 2 r_k (SExtractor-like), keeps it from being dominated by the outskirts
    sel = (rr <= max(2.0 * rk, rmin_fit)) & good
    den2 = float(np.sum(d[sel]))
    if den2 > 0:
        rk = float(np.sum(rr[sel] * d[sel])) / den2
    rap = max(kron_scale * rk, kron_min)
    flux = float(growth_curve(d, x, y, q, theta, [rap], mask)[0])
    return float(rk), flux, float(rap)


# --------------------------------------------------------------------------------------------------------------------- smoothness etc.
def smoothness(cut, mask_ell, sigma_s, noise_cut=None, centre=None):
    """S of the cutout inside the boolean aperture `mask_ell` (True = inside).  `noise_cut`: background patch (same shape) for the noise correction."""
    cut = np.asarray(cut, float)
    box = max(3, int(round(sigma_s)) | 1)                                   # boxcar of width `sigma_s` px (Conselice 2003: 0.25 R_P)
    sm = ndi.uniform_filter(cut, box, mode='nearest')
    if centre is not None:
        yy, xx = np.indices(cut.shape)
        mask_ell = mask_ell & (np.hypot(xx - centre[0], yy - centre[1]) > 0.5 * box)
    tot = float(cut[mask_ell].sum())
    if tot <= 0:
        return float('nan')
    a = float(np.abs(cut - sm)[mask_ell].sum())
    b = 0.0
    if noise_cut is not None:
        nm = ndi.uniform_filter(noise_cut, box, mode='nearest')
        b = float(np.abs(noise_cut - nm)[mask_ell].sum())
    return (a - b) / tot


def gini_m20(pixels, coords=None, xc=None, yc=None):
    """Gini and M20 of a set of pixel values (positive part), as Lotz et al. 2004."""
    p = np.asarray(pixels, float)
    n = p.size
    if n < 3 or p.sum() <= 0:
        return float('nan'), float('nan')
    s = np.sort(np.abs(p))
    mean = s.mean()
    i = np.arange(1, n + 1)
    g = float(np.sum((2 * i - n - 1) * s) / (mean * n * (n - 1)))
    m20 = float('nan')
    if coords is not None:
        xs, ys = coords
        w = np.clip(p, 0, None)
        if xc is None:
            xc = np.sum(w * xs) / w.sum(); yc = np.sum(w * ys) / w.sum()
        r2 = (xs - xc) ** 2 + (ys - yc) ** 2
        mi = w * r2
        order = np.argsort(-w)
        cum = np.cumsum(w[order])
        k = np.searchsorted(cum, 0.2 * w.sum()) + 1
        mtot = mi.sum()
        if mtot > 0:
            m20 = float(np.log10(max(mi[order][:k].sum(), 1e-30) / mtot))
    return g, m20


def fill_masked(data, mask, x, y, q, theta, dr=1.0):
    """Replace masked pixels by the mean of the unmasked pixels at the same elliptical radius (bins of `dr` px along the major axis, linear interpolation between bin
    centres); pixels in bins without data keep 0.  Used instead of sep's area-ratio correction, which is wrong when a masked blob sits in the outskirts."""
    mask = np.asarray(mask, bool)
    if not mask.any():
        return np.asarray(data, float)
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x, yy - y
    u = dx * math.cos(theta) + dy * math.sin(theta)
    v = -dx * math.sin(theta) + dy * math.cos(theta)
    rr = np.sqrt(u ** 2 + (v / max(q, 0.05)) ** 2)
    k = np.floor(rr / dr).astype(int)
    nb = int(k.max()) + 1
    ok = ~mask & np.isfinite(data)
    s = np.bincount(k[ok], weights=np.asarray(data, float)[ok], minlength=nb)
    c = np.bincount(k[ok], minlength=nb)
    prof = np.where(c > 0, s / np.maximum(c, 1), np.nan)
    cen = (np.arange(nb) + 0.5) * dr
    good = np.isfinite(prof)
    if good.sum() < 2:
        return np.where(mask, 0.0, data)
    fill = np.interp(rr, cen[good], prof[good])
    return np.where(mask, fill, data)


# --------------------------------------------------------------------------------------------------------------------- one object
def measure_object(data, x, y, q, theta, rms, rmax=None, mask=None, etas=(0.1, 0.2, 0.3), main_eta=0.2, kron_scale=2.5, kron_min=3.5, kron_rmax=None,
                   smooth_frac=0.25, noise_patch=None, seed=0, curve_fracs=(0.2, 0.5, 0.8, 0.9), want_curve=False):
    """All extended indices of one object; `x, y` 0-based; returns a dict (nan for what could not be measured) + optionally the growth curve."""
    ny, nx = data.shape
    q = float(min(max(q, 0.15), 1.0))
    edge = min(x, y, nx - 1 - x, ny - 1 - y)
    rmax_img = max(edge / max(q, 0.3) if False else edge, 4.0)
    rmax = float(min(rmax if rmax else 150.0, rmax_img))
    out = dict(rp={e: float('nan') for e in etas}, flux_p=float('nan'), r_frac={f: float('nan') for f in curve_fracs}, conc=float('nan'), kron_r=float('nan'), kron_flux=float('nan'),
               smooth=float('nan'), gini_p=float('nan'), m20_p=float('nan'), flag=0)
    if rmax < 6:
        out['flag'] |= 1
        return out
    if mask is not None and np.asarray(mask).any():
        data = fill_masked(data, mask, x, y, q, theta)
        mask_img = np.asarray(mask, bool)
        mask = None                                  # (the filled image is used everywhere below; the mask still excludes pixels from the smoothness / Gini segmentation)
    else:
        mask_img = None
    cv = Curve(data, x, y, q, theta, rmax, mask)
    rr, eta = petrosian_eta(cv)
    for e in etas:
        out['rp'][e] = petrosian_radius(rr, eta, e)
    rp = out['rp'].get(main_eta, float('nan'))
    if not np.isfinite(rp):
        out['flag'] |= 2
    # Petrosian flux and flux-fraction radii
    if np.isfinite(rp) and 2.0 * rp <= cv.r[-1]:
        fp = float(cv.flux(2.0 * rp))
        out['flux_p'] = fp
        if fp > 0:
            for f in curve_fracs:
                out['r_frac'][f] = cv.radius_of(f * fp)
            r20, r80 = out['r_frac'].get(0.2), out['r_frac'].get(0.8)
            if r20 and r80 and np.isfinite(r20) and np.isfinite(r80) and r20 > 0:
                out['conc'] = 5.0 * math.log10(r80 / r20)
    elif np.isfinite(rp):
        out['flag'] |= 4                       # 2 R_P beyond the cutout / image edge
    # Kron
    kr = kron_rmax or (min(6.0 * (rp if np.isfinite(rp) else 3.0), rmax))
    out['kron_r'], out['kron_flux'], _ = kron(data, x, y, q, theta, min(kr, rmax), kron_scale, kron_min, mask)
    # smoothness and Petrosian segmentation
    if np.isfinite(rp) and 1.5 * rp <= rmax and rp >= 2.0:
        h = int(math.ceil(1.5 * rp)) + 3
        x0, y0 = int(round(x)) - h, int(round(y)) - h
        xa, ya = max(0, x0), max(0, y0)
        xb, yb = min(nx, x0 + 2 * h + 1), min(ny, y0 + 2 * h + 1)
        cut = data[ya:yb, xa:xb].astype(float)
        cm = np.zeros(cut.shape, bool) if mask_img is None else mask_img[ya:yb, xa:xb]
        yy, xx = np.mgrid[ya:yb, xa:xb]
        dx, dy = xx - x, yy - y
        u = dx * math.cos(theta) + dy * math.sin(theta)
        v = -dx * math.sin(theta) + dy * math.cos(theta)
        rell = np.sqrt(u ** 2 + (v / q) ** 2)
        ap = (rell <= 1.5 * rp) & ~cm
        if noise_patch is not None and noise_patch.shape == cut.shape:
            nz = noise_patch
        else:
            nz = np.random.default_rng(seed).normal(0.0, rms, cut.shape)
        out['smooth'] = smoothness(cut, ap, max(smooth_frac * rp, 3.0), nz, centre=(x - xa, y - ya))
        # Lotz segmentation: inside 1.5 R_P and above the SB at R_P
        sb_rp = float(cv.flux(rp * 1.0) - cv.flux(rp * 0.8)) / (math.pi * q * (1 - 0.64) * rp ** 2) if np.isfinite(cv.flux(rp)) else float('nan')
        seg = ap & (cut >= sb_rp) if np.isfinite(sb_rp) and sb_rp > 0 else ap & (cut > 0)
        if seg.sum() >= 10:
            out['gini_p'], out['m20_p'] = gini_m20(cut[seg], (xx[seg], yy[seg]))
    if want_curve:
        out['curve'] = (cv.r.tolist(), cv.F.tolist(), (rr.tolist(), eta.tolist()))
    return out


# --------------------------------------------------------------------------------------------------------------------- analytic references (tests / docs)
def sersic_fraction(r_over_re, n):
    """Fraction of the total flux of a Sersic profile inside r (circular, in units of re)."""
    from ogfkit.models import sersic_bn
    bn = sersic_bn(n)
    return gammainc(2 * n, bn * np.asarray(r_over_re, float) ** (1.0 / n))


def sersic_radius_fraction(f, n):
    """r / re containing the fraction f of the total flux."""
    from ogfkit.models import sersic_bn
    bn = sersic_bn(n)
    return (gammaincinv(2 * n, f) / bn) ** n


def sersic_eta(r_over_re, n):
    """Analytic Petrosian eta(r) of a circular Sersic profile (local SB at r over mean SB inside r, point values)."""
    from ogfkit.models import sersic_bn
    bn = sersic_bn(n)
    r = np.asarray(r_over_re, float)
    I = np.exp(-bn * (r ** (1.0 / n) - 1))
    F = gammainc(2 * n, bn * r ** (1.0 / n))
    # total flux in units where Ie = 1: 2 pi n re^2 e^bn bn^{-2n} Gamma(2n)
    from scipy.special import gamma
    Ftot = 2 * math.pi * n * math.exp(bn) * bn ** (-2 * n) * gamma(2 * n)
    mean = F * Ftot / (math.pi * r ** 2)
    return I / mean


def sersic_eta_ring(r_over_re, n, npts=400):
    """eta as measured with the 0.8 r - 1.25 r annulus (what `petrosian_eta` computes) from the analytic profile."""
    from ogfkit.models import sersic_bn
    from scipy.special import gamma
    bn = sersic_bn(n)
    Ftot = 2 * math.pi * n * math.exp(bn) * bn ** (-2 * n) * gamma(2 * n)
    r = np.asarray(r_over_re, float)
    F = lambda s: gammainc(2 * n, bn * s ** (1.0 / n)) * Ftot
    ring = (F(1.25 * r) - F(0.8 * r)) / (math.pi * (1.25 ** 2 - 0.8 ** 2) * r ** 2)
    return ring / (F(r) / (math.pi * r ** 2))
