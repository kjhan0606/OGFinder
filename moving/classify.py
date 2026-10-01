"""Explainable rule-based classification of difference-image detections.

Per-detection class (single exposure) from morphology only:
  artefact_cr      cosmic-ray like: sharper than the PSF (sharp > 1.5 * PSF sharpness), few pixels, no PSF wings
  artefact_dipole  large negative fraction in the neighbourhood (misregistration residual of a bright star/galaxy)
  artefact_static  on top of a significant template source (static residual)
  artefact_edge    next to bad / saturated pixels or the detector edge / template hole
  trail            elongated (axis ratio >= trail_elong) and longer than min_trail_fwhm * PSF FWHM
  point            otherwise (candidate static transient OR a mover seen as a dot: decided later by linking)

Multi-exposure verdict (classify_linked) uses the tracklet fit (see tracklet.py):
  moving           >= 3 exposures consistent with constant sky velocity, rate above the astrometric noise
                   floor, or a single-exposure trail with consistent direction in the others
  static_transient detection at the same sky position (within tol) in >= 2 exposures of different epochs/visits, flux > 0,
                   not moving, associated or not with a host
  variable/AGN     positive and negative residuals at the same position in different epochs, or nuclear (offset < 0.3 arcsec)
  artefact         remaining single-exposure detections
Every decision stores its reasons in the `why` column.
"""
import numpy as np

DEFAULTS = dict(snr_min=6.0, cr_sharp_factor=1.6, cr_max_npix=6, dipole_negfrac=0.12, trail_elong=2.2,
                min_trail_fwhm=1.5, tpl_snr_max=8.0)


def classify_detection(d, fwhm_pix, p=None):
    p = dict(DEFAULTS, **(p or {}))
    why = []
    if d.get("sign", 1) < 0:
        return "negative", ["negative residual (sign<0)"]
    if d.get("channel") == "trail":
        return "trail", ["elongated component of smoothed difference: L=%.2f\" axis ratio %.1f" % (d.get("trail_len_arcsec", np.nan), d.get("elong", np.nan))]
    if d.get("channel") != "trail" and d.get("tpl_snr", 0.0) > p["tpl_snr_max"]:
        return "artefact_static", ["sits on a significant template source (template S/N %.1f > %.0f): residual of a static star/galaxy" % (d["tpl_snr"], p["tpl_snr_max"])]
    if d.get("near_bad") or d.get("near_sat"):
        return "artefact_edge", ["adjacent to bad/saturated pixels"]
    if d.get("neg_frac", 0) > p["dipole_negfrac"]:
        return "artefact_dipole", ["negative fraction %.2f > %.2f" % (d["neg_frac"], p["dipole_negfrac"])]
    sh, shp = d.get("sharp", np.nan), d.get("sharp_psf", np.nan)
    if np.isfinite(sh) and np.isfinite(shp) and sh > p["cr_sharp_factor"] * shp and d.get("npix", 99) <= p["cr_max_npix"] * 4:
        return "artefact_cr", ["sharpness %.2f > %.1f x PSF %.2f" % (sh, p["cr_sharp_factor"], shp)]
    # pipeline CR flag: only decisive for compact, non-elongated detections (a trailed mover is often partly mis-flagged
    # as CR by the archive pipeline, so elongated detections are kept and left to the multi-exposure linking)
    if d.get("on_cr") and not (np.isfinite(d.get("elong", np.nan)) and d["elong"] >= 2.2 and d.get("a_pix", 0) >= 4.0):
        return "artefact_cr", ["on pipeline CR-flagged pixel, compact (a=%.1f px)" % d.get("a_pix", np.nan)]
    if np.isfinite(d.get("elong", np.nan)) and d["elong"] >= p["trail_elong"] and d.get("a_pix", 0) * 2.355 * 2 > p["min_trail_fwhm"] * fwhm_pix:
        return "trail", ["axis ratio %.1f >= %.1f" % (d["elong"], p["trail_elong"])]
    if d["snr"] < p["snr_min"]:
        return "faint", ["S/N %.1f < %.1f" % (d["snr"], p["snr_min"])]
    return "point", ["S/N %.1f, axis ratio %.1f" % (d["snr"], d.get("elong", np.nan))]
