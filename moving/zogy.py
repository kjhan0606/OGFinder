"""ZOGY image subtraction (Zackay, Ofek & Gal-Yam 2016, ApJ 830, 27) in numpy, plus a
simple scaled-subtraction fallback.

Inputs are background-subtracted images N (new) and R (reference) on the same pixel grid,
their PSFs (normalised to unit sum, any odd size, centred), the per-image background noise
sigma_n, sigma_r (pixel rms of sky noise) and flux-scaling factors F_n, F_r.

Returns the proper difference image D (unit-variance white noise for background noise),
the matched-filter score image S (flux-weighted correlation with P_D), the point-source flux
image (in the units of N scaled by F_n/F_r... see below), its sigma, and P_D.

Notation follows the paper:  D^ = (F_r P^_r N^ - F_n P^_n R^)/sqrt(s_n^2 F_r^2 |P^_r|^2 + s_r^2 F_n^2 |P^_n|^2)
F_D = F_n F_r / sqrt(s_n^2 F_r^2 + s_r^2 F_n^2),  P^_D = F_n F_r P^_n P^_r / (F_D sqrt(...)),
S^ = F_D D^ conj(P^_D).  The transient flux estimate is  alpha = S / (F_D^2 sum P_D^2)  in units where the
reference flux scale is F_r (we set F_r = 1, F_n = flux ratio new/reference), with sigma_alpha = 1/(F_D sqrt(sum P_D^2)).

Flux scale: `alpha`/`sigma_alpha` are in the reference flux scale (F_r); `alpha_new`/`sigma_alpha_new` in the scale of N
(= alpha * F_n / F_r).  detect.difference_chip reports the latter (flux of a mover in counts of the target exposure).

Extensions (all optional, `zogy(N, R, Pn, Pr, sn, sr)` is unchanged and returns exactly the old keys):

* **Source-noise term** (`Vn`, `Vr`: variance maps of the Poisson noise of the sources in N and R, in the squared units of N / R).  Paper
  eqs. 26-28: S = S_N - S_R with S_N = k_n (x) N, S_R = k_r (x) R, k^_n = F_n F_r^2 conj(P^_n)|P^_r|^2 / den^2, k^_r = F_n^2 F_r conj(P^_r) |P^_n|^2 / den^2;
  V(S_N) = k_n^2 (x) (sigma_n^2 + Vn), V(S_R) = k_r^2 (x) (sigma_r^2 + Vr)  (pixel-squared real-space kernels, convolved with the variance maps).
* **Astrometric-registration terms** (`astrom_n`, `astrom_r` = (sigma_x, sigma_y) in pixels, scalars or maps): V_ast(S_N) = sigma_xn^2 (dS_N/dx)^2 +
  sigma_yn^2 (dS_N/dy)^2 (same for R), derivatives taken in Fourier space (paper eqs. 29-33).  A bright static star that is mis-registered by a few
  mas leaves a residual in S that this term down-weights.
* The corrected score is  S_corr = S / sqrt(V(S_N) + V(S_R) + V_ast(S_N) + V_ast(S_R));  `alpha`/`sigma_alpha` maps are consistent with it
  (sigma_alpha_map = sqrt(V_S) / (F_D^2 sum P_D^2), so alpha / sigma_alpha_map = S_corr).  With neither extension V_S = F_D^2 sum P_D^2 and S_corr = S/(F_D sqrt(sum P_D^2)).
* **Separate template / target PSFs** were always supported (`Pn`, `Pr`); `difference_chip` can now measure the template PSF from the template itself.
* **Spatially varying PSF**: `zogy_tiled` runs the above per tile (with an overlap margin that is discarded) with the PSF of each tile, see
  `imaging.PSFField`; a constant PSF field reproduces `zogy` up to the tile edge treatment.

Remaining limits: the PSF is piecewise constant (tile seams, no smooth interpolation between tiles); the noise is assumed uncorrelated (the template is
a resampled median whose noise is correlated between pixels; `sr` is measured empirically, the correlation is not modelled); the astrometric term uses one
sigma per image, not a distortion-residual map; flux scale F_n/F_r is a single number per chip; no colour-refraction / PSF-chromatic terms.
"""
import numpy as np
from numpy.fft import rfft2, irfft2, ifftshift


def _pad_psf(psf, shape):
    """Place PSF (odd-sized, centre at the middle) in a zero array of `shape` with its centre at [0,0]."""
    out = np.zeros(shape, dtype=np.float64)
    ph, pw = psf.shape
    cy, cx = ph // 2, pw // 2
    ys = (np.arange(ph) - cy) % shape[0]
    xs = (np.arange(pw) - cx) % shape[1]
    out[np.ix_(ys, xs)] = psf
    return out


def _freq_grids(shape):
    """Angular frequencies (rad / pixel) of the rfft2 grid: (qy, qx) broadcastable to (ny, nx//2+1)."""
    qy = 2.0 * np.pi * np.fft.fftfreq(shape[0])[:, None]
    qx = 2.0 * np.pi * np.fft.rfftfreq(shape[1])[None, :]
    return qy, qx


def _as_map(v, shape):
    v = np.asarray(v, np.float64)
    return np.broadcast_to(v, shape) if v.ndim == 0 else v


def zogy(N, R, Pn, Pr, sn, sr, Fn=1.0, Fr=1.0, eps=1e-12, Vn=None, Vr=None, astrom_n=None, astrom_r=None):
    """Compute ZOGY products. All arrays float; N, R equal shape. NaN-free input expected.

    Optional noise terms (see module docstring): `Vn`, `Vr` = source-noise variance maps of N and R (pixel variance of the Poisson noise of the
    sources, same units as N / R squared); `astrom_n`, `astrom_r` = (sigma_x, sigma_y) registration uncertainty in pixels (scalar or maps).
    With none of them the result is the classic sky-noise ZOGY.  Extra keys when any is given: `S_corr` (S / sqrt(V_S)), `V_S` (variance map of S),
    `sigma_alpha_map`, `sigma_alpha_new_map`.  `alpha` is always S / (F_D^2 sum P_D^2), the matched-filter flux of a point source."""
    shape = N.shape
    N = np.asarray(N, np.float64); R = np.asarray(R, np.float64)
    Nh = rfft2(N); Rh = rfft2(R)
    Pnh = rfft2(_pad_psf(Pn / Pn.sum(), shape)); Prh = rfft2(_pad_psf(Pr / Pr.sum(), shape))
    den = np.sqrt(sn ** 2 * Fr ** 2 * np.abs(Prh) ** 2 + sr ** 2 * Fn ** 2 * np.abs(Pnh) ** 2) + eps
    Dh = (Fr * Prh * Nh - Fn * Pnh * Rh) / den
    FD = Fn * Fr / np.sqrt(sn ** 2 * Fr ** 2 + sr ** 2 * Fn ** 2)
    PDh = Fn * Fr * Pnh * Prh / (FD * den)
    D = irfft2(Dh, s=shape)
    S = irfft2(FD * Dh * np.conj(PDh), s=shape)
    PD = irfft2(PDh, s=shape)
    sumP2 = float(np.sum(PD ** 2))
    # D is the "proper difference" normalised such that its pixel noise is 1/ (white) up to factor; by
    # construction (paper eq. 12 with the 1/den) the noise variance of D is sum |.|^2 -> estimate empirically
    alpha = S / (FD ** 2 * sumP2)
    sig_alpha = 1.0 / (FD * np.sqrt(sumP2))
    # `alpha` is in the flux scale of the REFERENCE (F_r): a source of true flux f in N has alpha = f * F_r / F_n.  `alpha_new` /
    # `sigma_alpha_new` convert to the flux scale of N (counts in the new image), which is what photometry of a mover needs.
    to_new = Fn / Fr
    out = dict(D=D, S=S, alpha=alpha, sigma_alpha=sig_alpha, alpha_new=alpha * to_new, sigma_alpha_new=sig_alpha * to_new,
               PD=np.fft.fftshift(PD), FD=FD, sumP2=sumP2)
    if Vn is None and Vr is None and astrom_n is None and astrom_r is None:
        return out
    # ---- corrected score (paper eqs. 26-33)
    den2 = den ** 2
    kn_h = Fn * Fr ** 2 * np.conj(Pnh) * np.abs(Prh) ** 2 / den2          # S_N = k_n (x) N
    kr_h = Fn ** 2 * Fr * np.conj(Prh) * np.abs(Pnh) ** 2 / den2          # S_R = k_r (x) R,  S = S_N - S_R
    kn2 = irfft2(kn_h, s=shape) ** 2; kr2 = irfft2(kr_h, s=shape) ** 2
    ones = np.ones(shape)
    # sky-noise variance of S is the constant  sn^2 sum(kn^2) + sr^2 sum(kr^2)  (== F_D^2 sum P_D^2 in the white-noise limit)
    VS = sn ** 2 * kn2.sum() * ones + sr ** 2 * kr2.sum() * ones
    if Vn is not None:
        VS = VS + irfft2(rfft2(kn2) * rfft2(np.asarray(Vn, np.float64)), s=shape)
    if Vr is not None:
        VS = VS + irfft2(rfft2(kr2) * rfft2(np.asarray(Vr, np.float64)), s=shape)
    if astrom_n is not None or astrom_r is not None:
        qy, qx = _freq_grids(shape)
        for img_h, kh, ast in ((Nh, kn_h, astrom_n), (Rh, kr_h, astrom_r)):
            if ast is None:
                continue
            sx, sy = ast
            SNh = kh * img_h
            dSdx = irfft2(1j * qx * SNh, s=shape); dSdy = irfft2(1j * qy * SNh, s=shape)
            VS = VS + _as_map(sx, shape) ** 2 * dSdx ** 2 + _as_map(sy, shape) ** 2 * dSdy ** 2
    VS = np.maximum(VS, 1e-30)
    out["V_S"] = VS
    out["S_corr"] = S / np.sqrt(VS)
    out["sigma_alpha_map"] = np.sqrt(VS) / (FD ** 2 * sumP2)
    out["sigma_alpha_new_map"] = out["sigma_alpha_map"] * to_new
    return out


def zogy_tiled(N, R, psf_n, psf_r, sn, sr, Fn=1.0, Fr=1.0, tile=(512, 512), margin=32, **kw):
    """ZOGY with spatially varying PSFs.  `psf_n`, `psf_r`: `imaging.PSFField` objects (or plain PSF arrays = constant).  The images are split into
    tiles of `tile` (ny, nx) pixels; each tile plus `margin` pixels on every side (discarded after the transform: FFT wrap-around) is processed with
    the PSF of the tile centre.  `sn`, `sr`, `Fn`, `Fr` and the keywords of `zogy` (Vn, Vr, astrom_n, astrom_r: scalars or full-size maps) are passed through.
    Returns a dict with the full-size maps of the keys of `zogy` (D, S, alpha, alpha_new, sigma_alpha(_new) as MAPS, S_corr/V_S when
    requested) plus `PD` of the centre tile and `FD`; the tile PSFs used are in `psf_used`."""
    shape = N.shape; ny, nx = shape
    keys = ("D", "S", "alpha", "alpha_new", "sigma_alpha", "sigma_alpha_new", "S_corr", "V_S")
    out = {}
    used = {}
    ty, tx = tile
    for y0 in range(0, ny, ty):
        for x0 in range(0, nx, tx):
            y1 = min(ny, y0 + ty); x1 = min(nx, x0 + tx)
            ya, yb = max(0, y0 - margin), min(ny, y1 + margin); xa, xb = max(0, x0 - margin), min(nx, x1 + margin)
            cy, cx = 0.5 * (y0 + y1), 0.5 * (x0 + x1)
            Pn = psf_n.at(cx, cy) if hasattr(psf_n, "at") else psf_n
            Pr = psf_r.at(cx, cy) if hasattr(psf_r, "at") else psf_r
            sub = {}
            for k_, v in kw.items():
                if k_ in ("Vn", "Vr") and v is not None:
                    sub[k_] = np.asarray(v)[ya:yb, xa:xb]
                elif k_ in ("astrom_n", "astrom_r") and v is not None:
                    sub[k_] = tuple(_slice_if_map(c_, ya, yb, xa, xb) for c_ in v)
                else:
                    sub[k_] = v
            o = zogy(N[ya:yb, xa:xb], R[ya:yb, xa:xb], Pn, Pr, sn, sr, Fn=Fn, Fr=Fr, **sub)
            core = (slice(y0 - ya, y1 - ya), slice(x0 - xa, x1 - xa))
            for k_ in keys:
                if k_ not in o:
                    continue
                v = o[k_]
                if k_ not in out:
                    out[k_] = np.zeros(shape)
                out[k_][y0:y1, x0:x1] = v[core] if np.ndim(v) == 2 else v
            used[(y0, x0)] = (Pn, Pr)
            if "PD" not in out or (y0 <= ny // 2 < y1 and x0 <= nx // 2 < x1):
                out["PD"] = o["PD"]; out["FD"] = o["FD"]; out["sumP2"] = o["sumP2"]
    out["psf_used"] = used
    return out


def _slice_if_map(v, ya, yb, xa, xb):
    v = np.asarray(v)
    return v if v.ndim == 0 else v[ya:yb, xa:xb]


def scaled_subtraction(N, R, Pn, Pr, Fn=1.0, Fr=1.0):
    """Fallback: match resolution by convolving each image with the other's PSF then subtract (scaled)."""
    from scipy.signal import fftconvolve
    Nc = fftconvolve(N, Pr, mode="same"); Rc = fftconvolve(R, Pn, mode="same")
    return Nc / Fn - Rc / Fr


def robust_sigma(a, mask=None):
    a = a[np.isfinite(a)] if mask is None else a[np.isfinite(a) & ~mask]
    if a.size == 0:
        return np.nan
    m = np.median(a)
    return 1.4826 * np.median(np.abs(a - m))
