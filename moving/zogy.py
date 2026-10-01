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

Not implemented (documented limitations): the source-noise term of S_corr (sky-noise-only normalisation; around very bright
stars the score noise is higher than 1, which is compensated empirically by the robust score normalisation and the
static-residual / dipole vetoes in classify), the astrometric-error terms of S_corr, and spatially varying PSFs (one PSF per
chip, measured from field stars in imaging.estimate_psf; the template PSF defaults to the target PSF).
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


def zogy(N, R, Pn, Pr, sn, sr, Fn=1.0, Fr=1.0, eps=1e-12):
    """Compute ZOGY products. All arrays float; N, R equal shape. NaN-free input expected."""
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
    return dict(D=D, S=S, alpha=alpha, sigma_alpha=sig_alpha, alpha_new=alpha * to_new, sigma_alpha_new=sig_alpha * to_new,
                PD=np.fft.fftshift(PD), FD=FD, sumP2=sumP2)


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
