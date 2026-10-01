"""Photometry of moving objects and transients, and absolute magnitude.

* aperture_flux(): forced aperture photometry on a (difference or science) image at a position -> e-/s -> AB mag.
* trail_flux(): flux of a trailed source by integrating the ZOGY difference flux map over the trail segment.
* hg_absolute(): H from V-like magnitudes with the IAU H,G system (Bowell et al. 1989 phase function approximated by the
  two-term form of the HG model; G=0.15 assumed unless fitted); the HST band magnitude is converted to V with a
  user-supplied colour term (default: V - F814W(AB)... NOT calibrated; assumed 0.3 for a C/S-type spectrum, FLAGGED).
* diameter_km(): D = 1329/sqrt(p) * 10^(-H/5) with an ASSUMED albedo (flagged).
"""
import numpy as np


def ab_mag(flux_e_s, zp_ab, err=None):
    f = np.asarray(flux_e_s, float)
    with np.errstate(all="ignore"):
        m = np.where(f > 0, zp_ab - 2.5 * np.log10(f), np.nan)
        if err is None:
            return m
        e = 1.0857 * np.asarray(err, float) / np.where(f > 0, f, np.nan)
    return m, e


def aperture_flux(img, x, y, r_pix=3.0, bkg_in=8.0, bkg_out=14.0, mask=None):
    """Sum in a circle, local background from an annulus (median).  Returns flux, flux_err(background scatter only)."""
    ny, nx = img.shape
    yy, xx = np.mgrid[max(0, int(y - bkg_out) - 1):min(ny, int(y + bkg_out) + 2), max(0, int(x - bkg_out) - 1):min(nx, int(x + bkg_out) + 2)]
    r = np.hypot(xx - x, yy - y)
    sub = img[yy, xx]
    ok = np.ones_like(sub, bool) if mask is None else ~mask[yy, xx]
    ann = (r >= bkg_in) & (r <= bkg_out) & ok
    sky = np.median(sub[ann]) if ann.sum() > 10 else 0.0
    sd = 1.4826 * np.median(np.abs(sub[ann] - sky)) if ann.sum() > 10 else np.nan
    ap = (r <= r_pix) & ok
    return float(np.sum(sub[ap] - sky)), float(sd * np.sqrt(ap.sum())) if np.isfinite(sd) else np.nan


def hg_phase(alpha_deg, G=0.15):
    a = np.radians(np.asarray(alpha_deg, float)) 
    # Bowell et al. (1989): phi_i = exp(-A_i * tan(alpha/2)^B_i); A1=3.33,B1=0.63 ; A2=1.87,B2=1.22
    phi1 = np.exp(-3.33 * np.tan(a / 2) ** 0.63)
    phi2 = np.exp(-1.87 * np.tan(a / 2) ** 1.22)
    return -2.5 * np.log10((1 - G) * phi1 + G * phi2)


def hg_absolute(v_mag, r_au, delta_au, alpha_deg, G=0.15):
    """H = V - 5 log10(r*Delta) - phase term."""
    return np.asarray(v_mag) - 5 * np.log10(np.asarray(r_au) * np.asarray(delta_au)) - hg_phase(alpha_deg, G)


def diameter_km(H, albedo=0.14):
    return 1329.0 / np.sqrt(albedo) * 10 ** (-np.asarray(H) / 5.0)
