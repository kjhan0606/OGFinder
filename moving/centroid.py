"""Trail-aware centroiding: fit a uniformly trailed Gaussian PSF (a line segment convolved with a Gaussian) to the difference-flux image.

A mover with apparent rate v leaves a trail of length v * t_exp (1.5-5 arcsec = 30-100 ACS pixels for the fast movers of the benchmark) in every
exposure.  The detection stage finds the trail as a >3-sigma connected component of the smoothed difference image; for faint trails only the
brightest part passes the threshold, cosmic-ray pixels and neighbours bias the moments, and the segment centre (the position that enters the linker)
moves by up to ~0.5 arcsec.  `fit_segment` fits  F(x, y) = A / L * integral_{-L/2}^{L/2} G(x - x0 - s cos(th), y - y0 - s sin(th); sigma) ds
(closed form: erf along the axis x Gaussian across) to the cutout by robust (soft-L1) least squares, with cosmic-ray / bad pixels excluded.
Returns the centre, length, angle, flux and an uncertainty from the Jacobian.  The pixel-sampled model is the exact line integral, not a sum of
sub-steps, so L -> 0 is a point source."""
import numpy as np
from scipy import optimize
from scipy.special import erf

_S2 = np.sqrt(2.0)


def segment_model(shape, x0, y0, L, th, sigma, amp, yy=None, xx=None):
    if yy is None:
        yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    c, s = np.cos(th), np.sin(th)
    dx = xx - x0; dy = yy - y0
    u = dx * c + dy * s                    # along the trail
    v = -dx * s + dy * c                   # across
    L = max(float(L), 1e-3)                # L -> 0 is a point source (the along-trail term / L tends to the Gaussian profile)
    h = 0.5 * L
    along = 0.5 * (erf((u + h) / (_S2 * sigma)) - erf((u - h) / (_S2 * sigma)))      # in [0, 1], = L * (box convolved with gaussian) / L
    across = np.exp(-0.5 * (v / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)
    # integral of the 2-D Gaussian along the segment = across * along * sqrt(2 pi) sigma / L * ... normalised so the total flux is `amp`
    return amp * along * across / L


def fit_segment(img, x0, y0, L0, th0, sigma, noise, mask=None, hw=None, min_hw=10, max_hw=70, free_sigma=False):
    """Fit one segment.  `img` = difference flux (e-/s per pixel), `noise` = per-pixel sigma (scalar or map), `mask` True = ignore.
    (x0, y0, L0, th0) = starting values in pixel coordinates of `img` (th in radians, +x toward +y).  Returns dict(ok, x, y, L, th, amp,
    sx, sy, chi2_red, n) in the coordinates of `img`."""
    ny, nx = img.shape
    if hw is None:
        hw = int(min(max_hw, max(min_hw, 0.5 * L0 + 4 * sigma + 6)))
    xi, yi = int(round(x0)), int(round(y0))
    xa, xb = max(0, xi - hw), min(nx, xi + hw + 1); ya, yb = max(0, yi - hw), min(ny, yi + hw + 1)
    cut = np.asarray(img[ya:yb, xa:xb], float)
    if cut.size < 25:
        return dict(ok=False)
    w = np.ones(cut.shape, bool) if mask is None else ~np.asarray(mask[ya:yb, xa:xb], bool)
    w &= np.isfinite(cut)
    if w.sum() < 25:
        return dict(ok=False)
    sg = np.broadcast_to(np.asarray(noise, float), img.shape)[ya:yb, xa:xb] if np.ndim(noise) else float(noise)
    yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
    A0 = float(np.clip(cut[w].sum(), 1e-6, None))
    p0 = np.array([x0 - xa, y0 - ya, max(L0, 0.0), th0, A0])
    lo = np.array([0, 0, 0.0, th0 - 1.2, 0.0]); hi = np.array([cut.shape[1] - 1, cut.shape[0] - 1, 2.0 * hw, th0 + 1.2, 50.0 * A0 + 1.0])
    xs = np.array([1.0, 1.0, max(2.0, 0.3 * L0), 0.1, A0])
    if free_sigma:
        p0 = np.r_[p0, sigma]; lo = np.r_[lo, 0.6 * sigma]; hi = np.r_[hi, 2.0 * sigma]; xs = np.r_[xs, 0.2]

    def resid(p):
        sig = p[5] if free_sigma else sigma
        m = segment_model(cut.shape, p[0], p[1], p[2], p[3], sig, p[4], yy, xx)
        return ((m - cut) / sg)[w] if np.ndim(sg) == 0 else ((m - cut)[w] / sg[w])
    try:
        r = optimize.least_squares(resid, np.clip(p0, lo, hi), bounds=(lo, hi), x_scale=xs, loss="soft_l1", f_scale=2.0, max_nfev=60)
    except Exception:
        return dict(ok=False)
    p = r.x
    chi2 = float(np.sum(resid(p) ** 2)) / max(int(w.sum()) - len(p), 1)
    try:
        J = r.jac; cov = np.linalg.inv(J.T @ J) * max(chi2, 1.0)
        sx, sy = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
    except Exception:
        sx = sy = np.nan
    return dict(ok=bool(r.success), x=float(p[0] + xa), y=float(p[1] + ya), L=float(p[2]), th=float(p[3]), amp=float(p[4]), sx=sx, sy=sy,
                chi2_red=chi2, n=int(w.sum()))
