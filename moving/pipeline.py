"""High-level steps used by the CLI driver (ds9_moving.py) and by the validation scripts.

Every step reads/writes explicit files inside a work directory so that it can be replayed from a recorded session:
  wd/chips.json          list of FITS inputs used (+ sha256)
  wd/align/*.json        per-chip alignment sidecars (see align.py)
  wd/diff/*.fits         difference-image flux maps (alpha, e-/s per pixel)
  wd/detections.tsv      all detections with features and single-exposure class
  wd/tracklets.json      linked tracklets
"""
import os, json, time
import numpy as np
from . import imaging as I, align as A, reference as R, detect as D, classify as C, tracklet as T, util
from .util import log


def detect_in_region(chips, center=None, half_pix=700, psf=None, min_inputs=2, snr_det=5.0, progress=None,
                     save_diff_dir=None, extra=None):
    """Difference every chip that covers `center` (ra, dec; or the whole chips if None) against the median of the
    other exposures and detect/classify.  `chips` must already be aligned.  Returns (dets, infos)."""
    if psf is None:
        psf = I.estimate_psf(chips, snr_min=6.0)
    psf_t, fw, n = psf
    ps0 = float(np.median([c.pixscale for c in chips]))
    fw_floor = 0.085 / ps0                      # ACS/WFC F814W FWHM ~0.09 arcsec; a narrower "empirical" PSF is a CR/hot-pixel artefact
    if psf_t is None or not np.isfinite(fw) or fw < fw_floor:
        sts = [A.detect_stars(c, snr_min=6.0) for c in chips[:4]]
        fws = [s_["fwhm"] for s_ in sts if s_ is not None]
        fw = max(float(np.median(fws)) if fws else fw_floor, fw_floor)
        psf_t = I.gaussian_psf(fw, 25); n = 0
        log("PSF: empirical estimate unusable (%s); Gaussian FWHM %.2f px" % ("none" if psf is None else "FWHM %.2f px < floor" % psf[1], fw))
    dets = []; infos = {}
    exs = sorted({c.path for c in chips}, key=lambda p: min(c.mjd for c in chips if c.path == p))
    for c in chips:
        if center is not None:
            x, y = c.wcs.all_world2pix([center[0]], [center[1]], 0)
            ny, nx = c.shape
            if not (-half_pix < x[0] < nx + half_pix and -half_pix < y[0] < ny + half_pix):
                continue
            x0 = int(max(0, x[0] - half_pix)); x1 = int(min(nx, x[0] + half_pix))
            y0 = int(max(0, y[0] - half_pix)); y1 = int(min(ny, y[0] + half_pix))
            if x1 - x0 < 100 or y1 - y0 < 100:
                continue
            cc = c.crop(x0, x1, y0, y1)
        else:
            cc = c
        others = [o for o in chips if o.path != c.path and o.filter == c.filter]
        tpl, nin = R.build_template(cc, others, min_inputs=min_inputs)
        if tpl is None:
            log("no template for", c.name)
            continue
        if progress:
            progress("difference %s" % c.name)
        out = D.difference_chip(cc, tpl, psf_t=psf_t, fw_t=fw, n_t=n, snr_det=snr_det)
        info = out["info"]; infos[c.name] = info
        ex = exs.index(c.path)
        for d in out["dets"]:
            d["ex"] = ex; d["t"] = float(c.mjd); d["chip"] = c.name; d["texp"] = float(c.texp); d["file"] = os.path.basename(c.path)
            d["cls"], d["why"] = C.classify_detection(d, info["fwhm_target"])
            # pixel position in the full chip (crop offset)
            ox, oy = getattr(cc, "origin", (0, 0))
            d["x_chip"] = d["x"] + ox; d["y_chip"] = d["y"] + oy
            d["tpl_snr_"] = d.get("tpl_snr"); d["zp_ab"] = c.zp_ab; d["pixscale"] = c.pixscale; d["filter"] = c.filter
            dets.append(d)
        if save_diff_dir:
            from astropy.io import fits
            os.makedirs(save_diff_dir, exist_ok=True)
            h = cc.wcs.to_header(); h["BUNIT"] = "e-/s"; h["ORIGFILE"] = os.path.basename(c.path); h["COMMENT"] = "ZOGY difference flux alpha"
            fits.writeto(os.path.join(save_diff_dir, "diff_%s.fits" % c.name.replace("[", "_").replace("]", "").replace(",", "_")),
                         out["alpha"].astype(np.float32), h, overwrite=True)
    return dets, infos


def link_detections(dets, chips_by_ex_offsets, snr_min=8.0, tol_arcsec=1.5, min_exposures=3,
                    classes=("trail", "point", "artefact_cr", "artefact_edge"), max_per_exposure=400, **kw):
    """Select detections and link them across exposures (see tracklet.link_exposures).

    The single-exposure class is only a soft veto here (sharpness test instead): in data with a dense archive cosmic-ray flag a real mover is often
    labelled `artefact_cr`/`artefact_edge` in one or two exposures, whereas random CR hits never line up on a constant-motion
    track.  Instead the brightest `max_per_exposure` positive detections of the allowed classes enter the linker, and the
    class of each member is reported with the tracklet (a tracklet needs >= min_exposures members, so a single mislabelled
    epoch is tolerated)."""
    sel = []
    # cosmic-ray / hot-pixel rejection that does not depend on the archive CR flag: a detection is "sharp" when the peak holds
    # a larger fraction of the 7x7 flux than a PSF-shaped source (trail-channel detections are always kept)
    def _ok(d):
        if d["channel"] == "trail":
            return True
        sh = d.get("sharp", np.nan)
        return (not np.isfinite(sh)) or sh < max(0.30, 1.35 * d.get("sharp_psf", 0.3))
    pool = [(i, d) for i, d in enumerate(dets) if d["sign"] > 0 and d["snr"] >= snr_min and _ok(d)
            and d["cls"] in classes + ("artefact_dipole", "faint")]
    keep = set()
    for ex in {d["ex"] for _, d in pool}:
        grp = sorted([(i, d) for i, d in pool if d["ex"] == ex], key=lambda x: -x[1]["snr"])[:max_per_exposure]
        keep.update(i for i, _ in grp)
    for i, d in enumerate(dets):
        if i in keep:
            sel.append(dict(ex=d["ex"], t=d["t"], ra=d["ra"], dec=d["dec"], sig=max(d.get("sig_pos_arcsec", 0.05), 0.05),
                            id=i, flux=d["flux_e_s"], cls=d["cls"], snr=d["snr"], texp=d.get("texp"), channel=d.get("channel", "point"),
                            trail_pa=d.get("pa_deg") if d["channel"] == "trail" else None,
                            trail_len=d.get("trail_len_arcsec") if d["channel"] == "trail" else None))
    trs = T.link_exposures(sel, tol_arcsec=tol_arcsec, min_exposures=min_exposures, obs_off_au=chips_by_ex_offsets, **kw)
    return trs


def hst_offsets(chips):
    """Geocentric offset (AU, equatorial) of HST at each exposure mid-time, from JPL Horizons (target -48)."""
    from . import horizons
    exs = sorted({c.path for c in chips}, key=lambda p: min(c.mjd for c in chips if c.path == p))
    tm = np.array([np.mean([c.mjd for c in chips if c.path == p]) for p in exs])
    jd = util.utc_mjd_to_tdb_jd(tm)
    v = horizons.vectors("-48", jd, center="500@399")[:, :3]
    return {i: v[i] for i in range(len(exs))}, tm
