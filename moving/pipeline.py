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
                     save_diff_dir=None, extra=None, psf_tile=None, psf_min_stars=5, **diff_kw):
    """Difference every chip that covers `center` (ra, dec; or the whole chips if None) against the median of the
    other exposures and detect/classify.  `chips` must already be aligned.  Returns (dets, infos).
    `diff_kw` goes to `detect.difference_chip` (cr_reject, trail_fit, realbogus, source_noise, astrom_sigma, template_psf, ...).
    `psf_tile=N` (pixels) measures a spatially varying PSF per N x N tile from the stars of all exposures of the same detector
    (`imaging.measure_psf_field`, tiles with < `psf_min_stars` stars use the pooled PSF) and runs the tiled ZOGY (`zogy_tiled`); None = one PSF per chip (default)."""
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
        kw_c = dict(diff_kw)
        if psf_tile:
            same = [o for o in chips if getattr(o, "extver", None) == getattr(c, "extver", None)]
            pf = I.measure_psf_field(same, cc.shape, tile=(int(psf_tile), int(psf_tile)), min_stars=psf_min_stars, snr_min=6.0, origin=getattr(cc, "origin", (0, 0)),
                                     constant=(psf_t, fw), min_fwhm_pix=fw_floor)
            kw_c["psf_field_t"] = pf; kw_c["psf_field_r"] = pf
            log("psf field %s: %s" % (c.name, pf.summary()))
        out = D.difference_chip(cc, tpl, nin=nin, psf_t=psf_t, fw_t=fw, n_t=n, snr_det=snr_det, **kw_c)
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


VETO_DEFAULTS = dict(stationary_arcsec=0.2, stationary_min_other=None, tpl_snr=4.0, tpl_min_a_pix=1.5, edge_pix=10)


def _fl(v):
    try:
        v = float(v)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def detection_veto(dets, chip_shapes=None, stationary_arcsec=0.2, stationary_min_other=None, tpl_snr=4.0, tpl_min_a_pix=1.5, edge_pix=10,
                   snr_min=8.0):
    """Detection-level veto masks applied BEFORE the per-exposure cap of the linker input pool (so that vetoed detections no longer
    take the pool slots of real movers).  Returns (veto: bool array over `dets`, counts: dict rule -> n vetoed).

    rules (a detection is vetoed when any applies; trail-channel detections are exempt from the residual rules because a mover's
    own trail legitimately overlaps nothing in the template):
      stationary  the same sky position (within `stationary_arcsec`) has a positive detection in all the other exposures
                  (`stationary_min_other` of them; default: every other exposure, at most 3) -> a static source / persistent
                  residual cannot be a mover whose displacement between exposures is >= ~0.3 arcsec;
      template    non-trail detection on a template source (S/N of the smoothed template > `tpl_snr`) that is extended
                  (a_pix > `tpl_min_a_pix`): residual halo/core of a static star or bright-galaxy core (the stricter template S/N > 8
                  case is already class `artefact_static`);
      edge        within `edge_pix` of the chip border (needs `chip_shapes`: {chip name: (ny, nx)} and `x_chip`, `y_chip`).
    Rules were tuned on inj1..inj8: they remove 0.6-2 % of the false pool detections and 0.0-0.5 % of the injected ones (see docs)."""
    from scipy.spatial import cKDTree
    n = len(dets); veto = np.zeros(n, bool); cnt = {"stationary": 0, "template": 0, "edge": 0}
    if n == 0:
        return veto, cnt
    pos = np.array([float(d["sign"]) > 0 and float(d["snr"]) >= snr_min for d in dets])
    exs = np.array([d["ex"] for d in dets]); nex = len(set(exs[pos])) if pos.any() else 0
    ra = np.array([d["ra"] for d in dets]); de = np.array([d["dec"] for d in dets])
    xy = np.c_[(ra - ra.mean()) * np.cos(np.radians(de.mean())) * 3600.0, (de - de.mean()) * 3600.0]
    need = stationary_min_other if stationary_min_other is not None else min(max(nex - 1, 2), 3)
    if nex >= 3 and stationary_arcsec:
        idx = np.nonzero(pos)[0]; tree = cKDTree(xy[idx])
        for q, nb in zip(idx, tree.query_ball_point(xy[idx], stationary_arcsec)):
            if len({exs[idx[j]] for j in nb if exs[idx[j]] != exs[q]}) >= need:
                veto[q] = True; cnt["stationary"] += 1
    for i, d in enumerate(dets):
        if not pos[i] or d.get("channel") == "trail":
            continue
        try:
            tp = float(d.get("tpl_snr")); a = float(d.get("a_pix"))
        except (TypeError, ValueError):
            tp = a = np.nan
        if tpl_snr and np.isfinite(tp) and tp > tpl_snr and np.isfinite(a) and a > tpl_min_a_pix:
            if not veto[i]:
                cnt["template"] += 1
            veto[i] = True
        if edge_pix and chip_shapes and d.get("chip") in chip_shapes and _fl(d.get("x_chip")) is not None and _fl(d.get("y_chip")) is not None:
            ny, nx = chip_shapes[d["chip"]]
            if min(_fl(d["x_chip"]), nx - 1 - _fl(d["x_chip"]), _fl(d["y_chip"]), ny - 1 - _fl(d["y_chip"])) < edge_pix:
                if not veto[i]:
                    cnt["edge"] += 1
                veto[i] = True
    return veto, cnt


def link_detections(dets, chips_by_ex_offsets, snr_min=8.0, tol_arcsec=1.5, min_exposures=3,
                    classes=("trail", "point", "artefact_cr", "artefact_edge"), max_per_exposure=400, veto=True, chip_shapes=None,
                    veto_stats=None, use_rb=True, rb_min=0.0, rb_snr_min=6.0, use_lac_cr=True, **kw):
    """Select detections and link them across exposures (see tracklet.link_exposures).

    `veto` (default True; or a dict of `detection_veto` options) removes stationary / template-residual / chip-edge detections from the pool
    before the per-exposure cap; `veto_stats` (dict) receives the counts.  `chip_shapes` {chip name: (ny, nx)} enables the edge rule.
    The single-exposure class is only a soft veto here (sharpness test instead): in data with a dense archive cosmic-ray flag a real mover is often
    labelled `artefact_cr`/`artefact_edge` in one or two exposures, whereas random CR hits never line up on a constant-motion
    track.  Instead the brightest `max_per_exposure` positive detections of the allowed classes enter the linker, and the
    class of each member is reported with the tracklet (a tracklet needs >= min_exposures members, so a single mislabelled
    epoch is tolerated)."""
    sel = []

    def _f(v):
        try:
            v = float(v)
            return v if np.isfinite(v) else None
        except (TypeError, ValueError):
            return None
    # cosmic-ray / hot-pixel rejection that does not depend on the archive CR flag: a detection is "sharp" when the peak holds
    # a larger fraction of the 7x7 flux than a PSF-shaped source (trail-channel detections are always kept)
    def _ok(d):
        if d["channel"] == "trail":
            return True
        sh = d.get("sharp", np.nan)
        return (not np.isfinite(sh)) or sh < max(0.30, 1.35 * d.get("sharp_psf", 0.3))
    vmask = np.zeros(len(dets), bool)
    if veto:
        vmask, vc = detection_veto(dets, chip_shapes, snr_min=snr_min, **(veto if isinstance(veto, dict) else {}))
        if veto_stats is not None:
            veto_stats.update(vc)
    snr_floor = min(snr_min, rb_snr_min) if use_rb and any(_f(d.get("rb")) is not None for d in dets) else snr_min
    pool = [(i, d) for i, d in enumerate(dets) if d["sign"] > 0 and d["snr"] >= snr_floor and _ok(d) and not vmask[i]
            and d["cls"] in classes + ("artefact_dipole", "faint")]
    # real/bogus ordering of the pool (detections carry `rb` when the model was available at detection time, `detect.difference_chip(realbogus=True)`):
    # the per-exposure cap then keeps the most real-looking detections instead of the brightest (cosmic-ray residuals are bright).  `rb_min` drops
    # detections below that probability.  Detections without `rb` (older detections.tsv, no model) are ordered by S/N exactly as before.
    def _rb(d):
        v = _f(d.get("rb")) if use_rb else None
        return v
    have_rb = use_rb and any(_rb(d) is not None for _, d in pool)
    if have_rb:
        # S/N threshold of the pool = rb_snr_min for scored detections (the model was trained from S/N 6), `snr_min` for unscored ones
        pool = [(i, d) for i, d in pool if (_rb(d) is None and d["snr"] >= snr_min) or (_rb(d) is not None and d["snr"] >= rb_snr_min and _rb(d) >= rb_min)]
        order_key = lambda x: -(_rb(x[1]) if _rb(x[1]) is not None else 0.0) * 1e6 - x[1]["snr"]
    else:
        order_key = lambda x: -x[1]["snr"]
    keep = set()
    for ex in {d["ex"] for _, d in pool}:
        grp = sorted([(i, d) for i, d in pool if d["ex"] == ex], key=order_key)[:max_per_exposure]
        keep.update(i for i, _ in grp)
    for i, d in enumerate(dets):
        if i in keep:
            sel.append(dict(ex=d["ex"], t=d["t"], ra=d["ra"], dec=d["dec"], sig=max(d.get("sig_pos_arcsec", 0.05), 0.05),
                            id=i, flux=d["flux_e_s"], cls=d["cls"], snr=d["snr"], texp=d.get("texp"), channel=d.get("channel", "point"),
                            trail_pa=d.get("pa_deg") if d["channel"] == "trail" else None,
                            trail_len=d.get("trail_len_arcsec") if d["channel"] == "trail" else None,
                            on_cr=(d.get("on_cr") in (True, 1, "True", "1")) or (use_lac_cr and d.get("lac3") in (1, "1", True, "True")), sharp=_f(d.get("sharp")), tpl_snr=_f(d.get("tpl_snr")),
                            x_chip=_f(d.get("x_chip")), y_chip=_f(d.get("y_chip")), chip=d.get("chip"), neg_frac=_f(d.get("neg_frac")),
                            a_pix=_f(d.get("a_pix")), elong=_f(d.get("elong")), rb=_f(d.get("rb"))))
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
