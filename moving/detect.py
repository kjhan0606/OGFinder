"""Difference imaging per chip and detection / feature extraction.

difference_chip(): target chip vs template (built from other exposures, see reference.py or an external
deep coadd) -> ZOGY (or scaled subtraction) -> matched-filter score S -> sep detection of positive and
negative peaks -> per-detection features.  Positions are in the (aligned) chip WCS.
"""
import numpy as np
import sep
from scipy import ndimage as ndi
from . import imaging as I, zogy as Z
from .util import log


def _star_flux_ratio(target, tplf, tpl_nan, bkg_t, bkg_r, rms_t):
    """Flux scale (target/template) from sources detected in BOTH images (cosmic rays and movers are absent
    from the template and are excluded by the cross-match); aperture sums, sigma-clipped median of ratios."""
    from scipy.spatial import cKDTree
    a = np.ascontiguousarray((target.data - bkg_t).astype(np.float32))
    b = np.ascontiguousarray((tplf - bkg_r).astype(np.float32))
    unusable = ndi.binary_dilation(target.bad | target.saturated | tpl_nan, iterations=8)
    rms_r = float(np.median(I.background(tplf, tpl_nan)[1]))
    ot = sep.extract(a, 10.0, err=np.ascontiguousarray(rms_t), minarea=10, mask=np.ascontiguousarray(target.bad | target.cr))
    orr = sep.extract(b, 10.0, minarea=10, thresh=None, err=np.full(b.shape, rms_r, np.float32), mask=np.ascontiguousarray(tpl_nan)) if False else \
        sep.extract(b, 10.0, err=np.full(b.shape, rms_r, np.float32), minarea=10, mask=np.ascontiguousarray(tpl_nan))
    if len(ot) < 3 or len(orr) < 3:
        return 1.0, 0
    d, j = cKDTree(np.c_[orr["x"], orr["y"]]).query(np.c_[ot["x"], ot["y"]], distance_upper_bound=1.5)
    m = np.isfinite(d)
    if m.sum() < 3:
        return 1.0, int(m.sum())
    x = ot["x"][m]; y = ot["y"][m]
    fl, _, _ = sep.sum_circle(a, x, y, 6.0)
    fr, _, _ = sep.sum_circle(b, x, y, 6.0)
    xi = np.clip(np.round(x).astype(int), 0, a.shape[1] - 1); yi = np.clip(np.round(y).astype(int), 0, a.shape[0] - 1)
    ok = (fl > 0) & (fr > 0) & np.isfinite(fl) & np.isfinite(fr) & ~unusable[yi, xi]
    if ok.sum() < 3:
        return 1.0, int(ok.sum())
    ratio = fl[ok] / fr[ok]
    med = np.median(ratio)
    for _ in range(3):
        sd = 1.4826 * np.median(np.abs(ratio - med)) + 0.01
        keep = np.abs(ratio - med) < 3 * sd
        if keep.sum() < 3:
            break
        med = np.median(ratio[keep])
    return float(med), int(keep.sum())


def _fallback_psf(chip, bkg, rms, size, min_fwhm_arcsec=0.085):
    """Gaussian with FWHM = 2 * median half-light radius of compact stars (valid for a Gaussian)."""
    from . import align as A
    st = A.detect_stars(chip, snr_min=6.0)
    if st is None:
        return I.gaussian_psf(2.2, size), 2.2
    f = max(st["fwhm"], 1.5, min_fwhm_arcsec / chip.pixscale)
    return I.gaussian_psf(f, size), f


def _template_chip(target, tplf, tpl_nan):
    """A Chip-like view of the template on the target grid (for PSF measurement from the template's own stars)."""
    import copy
    c = copy.copy(target)
    c.data = tplf.astype(np.float32); c.bad = np.asarray(tpl_nan) | target.bad
    c.cr = np.zeros(target.shape, bool); c.name = target.name + "[template]"; c.path = target.path + "#template"
    return c


def _astrom_sigma_pix(chip, nin=None, others=None):
    """Registration uncertainty per axis in pixels of a chip from its alignment sidecar (rms of the fit in mas; 0 for the anchor)."""
    a = getattr(chip, "align", None) or {}
    rms = a.get("rms_mas")
    if rms is None or not np.isfinite(rms):
        return None
    return float(rms) / 1000.0 / chip.pixscale


def difference_chip(target, template, nin=None, method="zogy", snr_det=5.0, psf_size=25,
                    min_pix=3, extra_mask=None, psf_t=None, fw_t=None, n_t=0, psf_r=None, fw_r=None, n_r=0,
                    source_noise=False, astrom_sigma=None, psf_field_t=None, psf_field_r=None, template_psf="target",
                    tile=(512, 512), tile_margin=32, astrom_template_sigma=None,
                    cr_reject=True, trail_fit=True, realbogus=True, trail_fit_max=40):
    """Returns dict(alpha=..., S=..., sigma_alpha=..., score=..., mask=..., dets=[...], info=...).

    Optional noise model / PSF extensions (all off by default: the defaults reproduce the previous behaviour exactly):
      source_noise    include the Poisson noise of the sources in N and R in the score variance (zogy `Vn`, `Vr`; image units are e-/s so
                      V = max(smoothed counts, 0) / t_exp; template: / (t_exp * n_inputs)).
      astrom_sigma    registration uncertainty of the target in pixels (float or (sx, sy)); True = take it from the chip's alignment sidecar
                      (`chip.align['rms_mas']`).  The template term (`astrom_template_sigma`, default: the same value / sqrt(n_inputs)) is added
                      too.  Enables the astrometric-error terms of ZOGY.
      psf_field_t / psf_field_r   `imaging.PSFField` (spatially varying PSF of the target / template); runs `zogy_tiled`.
      template_psf    "target" (default, as before) or "measure": measure the template PSF from the template's own stars (constant per chip,
                      or a PSFField if `psf_field_r` is given).

    Detection-stage options (item 1; `cr_reject=trail_fit=realbogus=False` reproduces the previous detections exactly):
      cr_reject       run L.A.Cosmic (`crrej.chip_lac`) on the target and add the features `lac3`, `lac_n7`, `lac_frac` (cosmic-ray evidence
                      independent of the archive DQ flag) to every detection.  It is a *feature*, not a veto (see `crrej`).
      realbogus       score every detection (S/N >= 6) with the shipped real/bogus model (`realbogus.py`, weights `moving/data/realbogus_model.json`):
                      `rb` in [0, 1]; the linker uses it to choose its input pool.  Detections get shape features (flux profile, PSF-fit chi2, fine
                      structure, ...) whether or not a model is available.
      trail_fit       re-measure the centre / length / angle of the `trail_fit_max` best trail-channel detections (by `rb`, else S/N) with a trailed-PSF
                      line-segment fit (`centroid.fit_segment`) that masks cosmic-ray pixels; the centroid error of long faint trails falls from
                      ~5 px to ~1.6 px median (docs/moving_objects.md).  Adds `x_mom`, `y_mom` (previous centroid) and `trail_fit`=1.
    """
    ny, nx = target.shape
    tpl = template.copy()
    valid = np.isfinite(tpl) & ~target.bad
    if extra_mask is not None:
        valid &= ~extra_mask
    tplf = np.where(np.isfinite(tpl), tpl, 0).astype(np.float32)
    bkg_t, rms_t = I.background(target.data, target.bad)
    bkg_r, rms_r = I.background(tplf, ~np.isfinite(tpl))
    if psf_t is None:
        psf_t, fw_t, n_t = I.estimate_psf(target, size=psf_size)
    info = dict(n_psf_stars_target=n_t, n_psf_stars_template=n_r)
    if psf_t is None or n_t < 6:
        psf_t, fw_t = _fallback_psf(target, bkg_t, rms_t, psf_size); info["psf_target"] = "gaussian(2*r50)"
    else:
        info["psf_target"] = "empirical"
    if psf_r is None and template_psf == "measure" and psf_field_r is None:
        tc = _template_chip(target, tplf, ~np.isfinite(tpl))
        pr_ = None; fr_ = np.nan; nr_ = 0
        try:
            st_ = I.collect_star_cutouts([tc], psf_size, 200, 6.0, confirm=False)[0]
            if len(st_) >= 6:
                pr_, fr_ = I._combine_cutouts(st_, psf_size); nr_ = len(st_)
        except Exception:
            pr_ = None
        if pr_ is not None:
            psf_r, fw_r, n_r = pr_, fr_, nr_
            info["n_psf_stars_template"] = nr_
    if psf_r is None:
        # template = median of cubic-spline-resampled exposures: its PSF is the target PSF broadened by the
        # (small) spline/median smoothing; we reuse the target PSF (documented approximation)
        psf_r = psf_t; fw_r = fw_t; info["psf_template"] = "target-psf"
    else:
        info["psf_template"] = "supplied"
    info["fwhm_target"] = float(fw_t); info["fwhm_template"] = float(fw_r)
    Fr_ratio, nst = _star_flux_ratio(target, tplf, ~np.isfinite(tpl), bkg_t, bkg_r, rms_t)
    info["flux_ratio"] = Fr_ratio; info["flux_ratio_nstars"] = nst
    sn = float(Z.robust_sigma((target.data - bkg_t)[valid & (target.data < np.percentile(target.data[valid], 90))]))
    sr = float(Z.robust_sigma((tplf - bkg_r)[valid & (tplf < np.percentile(tplf[valid], 90))]))
    N = np.where(valid, target.data - bkg_t, 0.0); R = np.where(valid, tplf - bkg_r, 0.0)
    info["sigma_target"] = sn; info["sigma_template"] = sr
    # significance of the TEMPLATE (static sky) smoothed with a PSF-sized kernel: detections sitting on a significant template
    # source are residuals of static stars/galaxies (misregistration, PSF mismatch), not movers or new sources.
    _k = max(1.0, fw_t / 2.355)
    tpl_sig = ndi.gaussian_filter(np.where(valid, tplf - bkg_r, 0.0), _k) / (sr / (2.0 * np.sqrt(np.pi) * _k))
    info["template_smooth_sigma"] = float(_k)
    if method == "zogy":
        zkw = {}
        if source_noise:
            nin_mean = float(np.mean(nin[valid])) if nin is not None and valid.any() else 3.0
            Vn = ndi.gaussian_filter(np.clip(N, 0, None), 1.0) / max(target.texp, 1e-3)
            Vr = ndi.gaussian_filter(np.clip(R, 0, None), 1.0) / (max(target.texp, 1e-3) * max(nin_mean, 1.0))
            zkw.update(Vn=Vn, Vr=Vr)
            info["source_noise"] = True
        if astrom_sigma is not None and astrom_sigma is not False:
            if astrom_sigma is True:
                sg = _astrom_sigma_pix(target)
            else:
                sg = astrom_sigma
            if sg is not None:
                sgx, sgy = (sg, sg) if np.isscalar(sg) else sg
                nin_mean = float(np.mean(nin[valid])) if nin is not None and valid.any() else 3.0
                st_ = astrom_template_sigma if astrom_template_sigma is not None else max(sgx, sgy) / np.sqrt(max(nin_mean, 1.0))
                zkw.update(astrom_n=(sgx, sgy), astrom_r=(st_, st_))
                info["astrom_sigma_pix"] = [float(sgx), float(sgy)]; info["astrom_sigma_template_pix"] = float(st_)
        if psf_field_t is not None or psf_field_r is not None:
            out = Z.zogy_tiled(N, R, psf_field_t if psf_field_t is not None else psf_t, psf_field_r if psf_field_r is not None else psf_r, sn, sr,
                               Fn=Fr_ratio, Fr=1.0, tile=tile, margin=tile_margin, **zkw)
            info["psf_variation"] = (psf_field_t.summary() if psf_field_t is not None else "constant")
        else:
            out = Z.zogy(N, R, psf_t, psf_r, sn, sr, Fn=Fr_ratio, Fr=1.0, **zkw)
        if "S_corr" in out:
            # corrected score: unit-variance map including source noise / astrometric terms; flux maps consistent with it
            S = out["S"]; alpha = out["alpha_new"]
            sig_map = out["sigma_alpha_new_map"] if "sigma_alpha_new_map" in out else out["sigma_alpha_new"]
            sig_alpha = float(np.median(sig_map)) if np.ndim(sig_map) else float(sig_map)
            info["PD_sum2"] = out["sumP2"]
            score_raw = out["S_corr"]
            s_emp = Z.robust_sigma(score_raw, ~valid)
            score = score_raw / s_emp
            info["score_sigma_emp"] = float(s_emp); info["sigma_alpha_theory"] = float(sig_alpha)
            info["sigma_alpha_map_p5_p95"] = [float(np.percentile(sig_map[valid], 5)), float(np.percentile(sig_map[valid], 95))] if np.ndim(sig_map) else None
            S = score_raw
        else:
            S = out["S"]; alpha = out["alpha_new"]; sig_alpha = out["sigma_alpha_new"]       # flux in the scale of the target exposure
            if np.ndim(sig_alpha):
                sig_alpha = float(np.median(sig_alpha))
            info["PD_sum2"] = out["sumP2"]
            # empirical score normalisation (robust) -- guards against non-white noise
            s_emp = Z.robust_sigma(S, ~valid)
            score = S / s_emp
            info["score_sigma_emp"] = float(s_emp)
            info["sigma_alpha_theory"] = float(sig_alpha)
        PD = out["PD"]
    else:
        D = Z.scaled_subtraction(N, R, psf_t, psf_r, Fn=Fr_ratio, Fr=1.0)
        from scipy.signal import fftconvolve
        S = fftconvolve(D, psf_t, mode="same"); s_emp = Z.robust_sigma(S, ~valid)
        score = S / s_emp
        alpha = S / np.sum(psf_t ** 2); sig_alpha = s_emp / np.sum(psf_t ** 2)
        PD = psf_t
    # mask edges and invalid area (grown by PSF half width)
    inval = ndi.binary_dilation(~valid, iterations=3)
    score = np.where(inval, 0.0, score).astype(np.float32)
    alpha = np.where(inval, 0.0, alpha).astype(np.float32)
    dets = []
    for sign in (+1, -1):
        sc = np.ascontiguousarray((sign * score).astype(np.float32))
        o = sep.extract(sc, snr_det, minarea=min_pix, deblend_nthresh=16, deblend_cont=0.05)
        for k in range(len(o)):
            dets.append(_det_features(o[k], sign, score, alpha, target, tplf, valid, psf_t, bkg_t, rms_t, sig_alpha, fw_t))
    # trail channel: connected components of the smoothed difference flux (> 3 sigma), kept when clearly
    # elongated (intra-exposure trails of movers).  Point detections inside a trail are merged into it.
    sm = ndi.gaussian_filter(np.where(inval, 0.0, alpha), 1.5)
    sm_sig = Z.robust_sigma(sm, inval)
    lab, nl = ndi.label(ndi.binary_dilation(sm > 3.0 * sm_sig, iterations=1))
    trail_masks = []
    if nl:
        sizes = ndi.sum(np.ones_like(lab), lab, index=np.arange(1, nl + 1))
        for li in np.where(sizes >= 25)[0] + 1:
            sl = ndi.find_objects((lab == li).astype(int))[0]
            m = lab[sl] == li
            # restrict to the core of the component (>= max(3 sigma, 25 % of its peak)): halos of bright movers and
            # merged neighbours otherwise bias the centroid by up to ~1 arcsec
            core = m & (sm[sl] >= max(3.0 * sm_sig, 0.25 * float(sm[sl][m].max())))
            if core.sum() >= 12:
                m = core
            yy, xx = np.nonzero(m)
            w = np.clip(sm[sl][m], 0, None)
            if w.sum() <= 0 or len(xx) < 25:
                continue
            xm = np.average(xx, weights=w); ym = np.average(yy, weights=w)
            C = np.cov(np.vstack([xx - xm, yy - ym]), aweights=w)
            ev, evec = np.linalg.eigh(C)
            L = 4.0 * np.sqrt(max(ev[1], 1e-9)); Wd = 4.0 * np.sqrt(max(ev[0], 1e-9))      # full extents ~ +-2 sigma
            if L / max(Wd, 1.0) < 3.0 or L < 8:
                continue
            ax_ = evec[:, 1]
            proj = (xx - xm) * ax_[0] + (yy - ym) * ax_[1]
            ends = [(xm + sl[1].start + proj.min() * ax_[0], ym + sl[0].start + proj.min() * ax_[1]),
                    (xm + sl[1].start + proj.max() * ax_[0], ym + sl[0].start + proj.max() * ax_[1])]
            fl = float(np.sum(alpha[sl][m]))
            cx = xm + sl[1].start; cy = ym + sl[0].start
            r1, d1 = target.wcs.all_pix2world([ends[0][0], ends[1][0], cx], [ends[0][1], ends[1][1], cy], 0)
            dd = dict(sign=+1, x=float(cx), y=float(cy), snr=float(fl / (sig_alpha * np.sqrt(len(xx)))), flux_e_s=fl, flux_err=float(sig_alpha * np.sqrt(len(xx))),
                      a_pix=float(L / 2), b_pix=float(Wd / 2), theta_pix=float(np.arctan2(ax_[1], ax_[0])), npix=int(len(xx)),
                      psf_sigma_pix=float(fw_t / 2.355), elong=float(L / max(Wd, 1e-3)), neg_frac=0.0, pos_frac=1.0, sharp=0.0, sharp_psf=0.0,
                      on_cr=False, saturated=False, near_bad=False, channel="trail", trail_len_pix=float(L),
                      trail_len_arcsec=float(L * target.pixscale),
                      trail_ends_ra=[float(r1[0]), float(r1[1])], trail_ends_dec=[float(d1[0]), float(d1[1])],
                      ra=float(r1[2]), dec=float(d1[2]), sig_pos_arcsec=float(max(0.1 * Wd, 0.5) * target.pixscale))
            dx_, dy_ = np.cos(dd["theta_pix"]), np.sin(dd["theta_pix"])
            rr_, ee_ = target.wcs.all_pix2world([cx, cx + dx_], [cy, cy + dy_], 0)
            dd["pa_deg"] = float(np.degrees(np.arctan2((rr_[1] - rr_[0]) * np.cos(np.radians(d1[2])), ee_[1] - ee_[0])) % 180.0)
            dd["a_arcsec"] = float(L / 2 * target.pixscale)
            mk = np.zeros(alpha.shape, bool); mk[sl] = m
            trail_masks.append(ndi.binary_dilation(mk, iterations=2))
            dets.append(dd)
    if trail_masks:
        tm = np.any(trail_masks, axis=0)
        dets = [x for x in dets if x.get("channel") == "trail" or not tm[min(max(int(round(x["y"])), 0), ny - 1), min(max(int(round(x["x"])), 0), nx - 1)]]
    sat_d = ndi.binary_dilation(target.saturated, iterations=6)
    for dd in dets:
        xi_, yi_ = int(round(dd["x"])), int(round(dd["y"]))
        dd["tpl_snr"] = float(tpl_sig[min(max(yi_, 0), ny - 1), min(max(xi_, 0), nx - 1)])
        dd["near_sat"] = bool(sat_d[min(max(yi_, 0), ny - 1), min(max(xi_, 0), nx - 1)])
        dd.setdefault("channel", "point")
    dets = _dedupe(dets)
    lac = None
    if cr_reject or realbogus or trail_fit:
        from . import crrej
        if cr_reject:
            lac = crrej.chip_lac(target)
            info["lac_frac"] = float(lac.mean()); info["lac_backend"] = crrej.backend()
        _shape_features(dets, target, bkg_t, rms_t, lac, fw_t, trail_mask=tm if trail_masks else None, alpha=alpha)
        if realbogus:
            from . import realbogus as RB
            info["realbogus"] = RB.score_dets(dets)
        if trail_fit:
            info["trail_fit"] = _refit_trails(dets, alpha, target, lac, fw_t, sig_alpha, trail_fit_max)
    info["n_pos"] = sum(1 for d in dets if d["sign"] > 0); info["n_neg"] = sum(1 for d in dets if d["sign"] < 0)
    return dict(alpha=alpha, score=score, mask=inval, dets=dets, info=info, sigma_alpha=float(sig_alpha),
                psf=psf_t, bkg=bkg_t, rms=rms_t)


RB_MIN_SNR = 6.0


def _shape_features(dets, target, bkg, rms, lac, fwhm_psf, trail_mask=None, alpha=None, r=10):
    """Shape / cosmic-ray features of every detection from the raw (background-subtracted) target cutout around its position.
    point channel: flux fractions in rings around the centroid (`f_r1`: r<=1 px, `f_r2`: 1<r<=2, `f_r3`: 2<r<=3.5 of the total within 3.5 px;
    a cosmic-ray hit has f_r1 -> 1, a PSF with FWHM ~1.7 px about 0.55), `pk_nb` (brightest neighbour of the peak / peak), `psf_chi2` (reduced chi2
    of a Gaussian of the PSF width fitted to the 9x9 cutout, amplitude + constant free), `psf_amp` (fitted amplitude / peak pixel), `fine`
    (L.A.Cosmic fine-structure ratio at the peak: (median3 - median7(median3)) / rms), `lap` (positive Laplacian at the peak / rms), `n_hi`
    (pixels > 3 rms in 7x7), `asym` (flux asymmetry of the 7x7 cutout), `lac3` / `lac_n7` (L.A.Cosmic pixels in 3x3 / 7x7; 0 without a mask),
    `arch3` (archive CR-flag pixels in 3x3), `n_near` (other S/N>=6 detections within 10 px: crowding / fragmentation).
    trail channel: `lac_frac` / `arch_frac` = fraction of the trail core mask on L.A.Cosmic / archive-flagged pixels, `fill` = npix / (L * 2.355 sigma).
    Positions are (x, y) of the detection in the target chip (already cropped)."""
    from scipy.spatial import cKDTree
    ny, nx = target.shape
    sig_psf = max(fwhm_psf / 2.355, 0.8)
    arch = getattr(target, "cr", None)
    data = target.data
    pos = [d for d in dets if d["sign"] > 0 and d["snr"] >= RB_MIN_SNR]
    if pos:
        xy = np.array([[d["x"], d["y"]] for d in pos])
        tree = cKDTree(xy)
        nnear = np.array([len(v) - 1 for v in tree.query_ball_point(xy, 10.0)])
    for i, d in enumerate(pos):
        d["n_near"] = int(nnear[i])
        if d.get("channel") == "trail":
            continue
        xi, yi = int(round(d["x"])), int(round(d["y"]))
        y0, y1, x0, x1 = max(0, yi - r), min(ny, yi + r + 1), max(0, xi - r), min(nx, xi + r + 1)
        cut = (data[y0:y1, x0:x1] - bkg[y0:y1, x0:x1]).astype(np.float64)
        rm = float(np.median(rms[y0:y1, x0:x1])) + 1e-9
        py, px = yi - y0, xi - x0
        if not (2 <= py < cut.shape[0] - 2 and 2 <= px < cut.shape[1] - 2):
            continue
        cx, cy = d["x"] - x0, d["y"] - y0
        yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
        rr = np.hypot(xx - cx, yy - cy)
        tot = cut[rr <= 3.5].sum()
        if tot > 0:
            d["f_r1"] = float(cut[rr <= 1.0].sum() / tot)
            d["f_r2"] = float(cut[(rr > 1.0) & (rr <= 2.0)].sum() / tot)
            d["f_r3"] = float(cut[(rr > 2.0) & (rr <= 3.5)].sum() / tot)
        pk = cut[py, px]
        nb = cut[py - 1:py + 2, px - 1:px + 2].copy(); nb[1, 1] = -np.inf
        d["pk_nb"] = float(np.max(nb) / pk) if pk > 0 else np.nan
        s7 = cut[max(0, py - 3):py + 4, max(0, px - 3):px + 4]
        d["n_hi"] = int((s7 > 3.0 * rm).sum())
        h = s7.shape[0] // 2
        tt = s7.sum()
        if tt > 0 and s7.shape == (7, 7):
            d["asym"] = float((abs(s7[:, :3].sum() - s7[:, 4:].sum()) + abs(s7[:3, :].sum() - s7[4:, :].sum())) / tt)
        # fine structure and Laplacian (L.A.Cosmic statistics at the peak)
        win = cut[max(0, py - 6):py + 7, max(0, px - 6):px + 7]
        m3 = ndi.median_filter(win, size=3, mode="nearest"); m37 = ndi.median_filter(m3, size=7, mode="nearest")
        wy, wx = min(py, 6), min(px, 6)
        d["fine"] = float((m3[wy, wx] - m37[wy, wx]) / rm)
        d["lap"] = float(max(4 * cut[py, px] - cut[py - 1, px] - cut[py + 1, px] - cut[py, px - 1] - cut[py, px + 1], 0.0) / (rm * np.sqrt(20.0)))
        # Gaussian PSF fit at the detection centroid (linear in amplitude and constant)
        g = np.exp(-0.5 * (rr ** 2) / sig_psf ** 2)
        m = rr <= 4.5
        A = np.stack([g[m], np.ones(m.sum())], 1)
        coef, *_ = np.linalg.lstsq(A, cut[m], rcond=None)
        res = cut[m] - A @ coef
        d["psf_chi2"] = float(np.sum(res ** 2) / rm ** 2 / max(m.sum() - 2, 1))
        d["psf_amp"] = float(coef[0] / pk) if pk > 0 else np.nan
        if lac is not None:
            ly0, lx0 = max(0, yi - 1), max(0, xi - 1)
            d["lac3"] = int(lac[ly0:yi + 2, lx0:xi + 2].any())
            d["lac_n7"] = int(lac[max(0, yi - 3):yi + 4, max(0, xi - 3):xi + 4].sum())
        else:
            d["lac3"] = 0; d["lac_n7"] = 0
        if arch is not None:
            d["arch3"] = int(arch[max(0, yi - 1):yi + 2, max(0, xi - 1):xi + 2].any())
    # trail channel: coverage of the component by cosmic-ray pixels
    for d in pos:
        if d.get("channel") != "trail":
            continue
        L = float(d.get("trail_len_pix", 2 * d.get("a_pix", 1.0))); th = d.get("theta_pix", 0.0); th = float(th) if np.isfinite(th) else 0.0
        hl = 0.5 * L
        R = int(np.ceil(hl)) + 4
        xi, yi = int(round(d["x"])), int(round(d["y"]))
        y0, y1, x0, x1 = max(0, yi - R), min(ny, yi + R + 1), max(0, xi - R), min(nx, xi + R + 1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        u = (xx - d["x"]) * np.cos(th) + (yy - d["y"]) * np.sin(th); v = -(xx - d["x"]) * np.sin(th) + (yy - d["y"]) * np.cos(th)
        sel = (np.abs(u) <= hl) & (np.abs(v) <= 2.5)
        if sel.sum() > 0:
            d["lac_frac"] = float(lac[y0:y1, x0:x1][sel].mean()) if lac is not None else 0.0
            d["arch_frac"] = float(arch[y0:y1, x0:x1][sel].mean()) if arch is not None else 0.0
            d["fill"] = float(d.get("npix", 0) / max(L * 2.355 * sig_psf, 1.0))
            if alpha is not None:
                sub = alpha[y0:y1, x0:x1][sel]
                d["trail_peak_frac"] = float(sub.max() / max(sub.sum(), 1e-9)) if sub.sum() > 0 else np.nan


def _refit_trails(dets, alpha, target, lac, fwhm_psf, sig_alpha, nmax):
    """Trailed-PSF segment fit of the best `nmax` trail detections (by `rb`, else S/N) - see `centroid.fit_segment`.  Updates x, y, ra, dec, the
    trail end points, length, position angle and `sig_pos_arcsec` of the accepted fits; keeps the moment centroid in `x_mom`, `y_mom`."""
    from . import centroid as CE
    tr = [d for d in dets if d.get("channel") == "trail" and d["sign"] > 0 and d["snr"] >= 8.0]
    tr.sort(key=lambda d: -(d.get("rb") if d.get("rb") is not None and np.isfinite(d.get("rb", np.nan)) else 0.0) * 1e6 - d["snr"])
    nz = float(1.4826 * np.median(np.abs(alpha[alpha != 0] - np.median(alpha[alpha != 0])))) if np.any(alpha != 0) else 1.0
    sg = max(fwhm_psf / 2.355, 0.8)
    mask = (target.bad | (lac if lac is not None else False))
    ok = 0; moved = 0
    for d in tr[:nmax]:
        L0 = float(d.get("trail_len_pix", 2 * d.get("a_pix", 3.0))); th0 = float(d.get("theta_pix", 0.0))
        best = None
        for f in (1.0, 2.0):
            r = CE.fit_segment(alpha, d["x"], d["y"], L0 * f, th0, sg, nz, mask=mask, max_hw=120)
            if r.get("ok") and (best is None or r["chi2_red"] < best["chi2_red"]):
                best = r
        if best is None or best["L"] < 4.0 or not np.isfinite(best["sx"]) and False:
            continue
        ok += 1
        d["x_mom"], d["y_mom"] = d["x"], d["y"]
        moved += int(np.hypot(best["x"] - d["x"], best["y"] - d["y"]) > 2.0)
        L = best["L"]; th = best["th"]; cx, cy = best["x"], best["y"]
        ax = np.array([np.cos(th), np.sin(th)])
        ends = [(cx - 0.5 * L * ax[0], cy - 0.5 * L * ax[1]), (cx + 0.5 * L * ax[0], cy + 0.5 * L * ax[1])]
        r1, d1 = target.wcs.all_pix2world([ends[0][0], ends[1][0], cx, cx + ax[0]], [ends[0][1], ends[1][1], cy, cy + ax[1]], 0)
        d["x"], d["y"] = float(cx), float(cy)
        d["ra"], d["dec"] = float(r1[2]), float(d1[2])
        d["trail_ends_ra"] = [float(r1[0]), float(r1[1])]; d["trail_ends_dec"] = [float(d1[0]), float(d1[1])]
        d["trail_len_pix"] = float(L); d["trail_len_arcsec"] = float(L * target.pixscale); d["theta_pix"] = float(th)
        d["pa_deg"] = float(np.degrees(np.arctan2((r1[3] - r1[2]) * np.cos(np.radians(d1[2])), d1[3] - d1[2])) % 180.0)
        d["trail_fit"] = 1; d["trail_chi2"] = float(best["chi2_red"])
        d["sig_pos_arcsec"] = float(max(np.nanmax([best["sx"], best["sy"], 0.0]), 0.3) * target.pixscale)
        d["a_pix"] = float(L / 2); d["a_arcsec"] = float(L / 2 * target.pixscale); d["elong"] = float(max(L, 1.0) / (2.355 * sg))
    return dict(candidates=len(tr), fitted=ok, moved_gt2px=moved)


def _det_features(o, sign, score, alpha, target, tpl, valid, psf, bkg, rms, sig_alpha, fwhm_psf):
    from . import align as A
    x, y = float(o["x"]), float(o["y"])
    xi, yi = int(round(x)), int(round(y))
    r = 12
    y0, y1, x0, x1 = max(0, yi - r), min(target.shape[0], yi + r + 1), max(0, xi - r), min(target.shape[1], xi + r + 1)
    sc = score[y0:y1, x0:x1] * sign
    # peak-weighted centroid on the positive part of the score (above 40% of peak) -> astrometry
    pk = sc.max()
    w = np.clip(sc - 0.4 * pk, 0, None)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    if w.sum() > 0:
        cx = float((w * xx).sum() / w.sum()); cy = float((w * yy).sum() / w.sum())
    else:
        cx, cy = x, y
    # second moments of the image cutout (new image minus background) -> size and elongation
    cut = (target.data[y0:y1, x0:x1] - bkg[y0:y1, x0:x1]) * sign
    cut_pos = np.clip(cut - 0.0, 0, None)
    # threshold the cutout at 2 sigma to avoid noise in the moments
    thr = 2.0 * np.median(rms[y0:y1, x0:x1])
    m = cut_pos > thr
    lab, nl = ndi.label(m)
    ll = lab[min(max(yi - y0, 0), lab.shape[0] - 1), min(max(xi - x0, 0), lab.shape[1] - 1)]
    if ll > 0:
        mm = lab == ll
        wts = cut_pos * mm
        s0 = wts.sum()
        mx = (wts * xx).sum() / s0; my = (wts * yy).sum() / s0
        cxx = (wts * (xx - mx) ** 2).sum() / s0; cyy = (wts * (yy - my) ** 2).sum() / s0; cxy = (wts * (xx - mx) * (yy - my)).sum() / s0
        tr = cxx + cyy; det = cxx * cyy - cxy ** 2
        disc = np.sqrt(max(tr * tr / 4 - det, 0))
        l1, l2 = tr / 2 + disc, max(tr / 2 - disc, 1e-6)
        A_ = np.sqrt(l1); B_ = np.sqrt(l2)
        theta_pix = 0.5 * np.arctan2(2 * cxy, cxx - cyy)       # angle of major axis from +x toward +y
        npix = int(mm.sum())
    else:
        A_ = B_ = np.nan; theta_pix = np.nan; npix = 0
    # PSF second moment for reference
    ph = psf.shape[0] // 2
    py, px = np.mgrid[-ph:ph + 1, -ph:ph + 1]
    psf_sig = np.sqrt(np.sum(psf * (px ** 2 + py ** 2)) / 2.0)
    flux = float(alpha[yi, xi]) if 0 <= yi < alpha.shape[0] and 0 <= xi < alpha.shape[1] else np.nan
    # negative pixel fraction of the difference neighbourhood (dipoles / bad subtraction)
    neigh = score[y0:y1, x0:x1] * sign
    negfrac = float(np.mean(neigh < -3.0))
    posfrac = float(np.mean(neigh > 3.0))
    # sharpness: peak of the new image cutout relative to flux within 3 px  (CR have > PSF value)
    ps = np.max(cut[max(0, yi - y0 - 1):yi - y0 + 2, max(0, xi - x0 - 1):xi - x0 + 2]) if cut.size else np.nan
    fl3 = cut[max(0, yi - y0 - 3):yi - y0 + 4, max(0, xi - x0 - 3):xi - x0 + 4].sum()
    sharp = float(ps / fl3) if fl3 > 0 else np.nan
    sharp_psf = float(psf.max() / psf[ph - 3:ph + 4, ph - 3:ph + 4].sum())
    d = dict(sign=sign, x=cx, y=cy, snr=float(pk), flux_e_s=flux, flux_err=float(sig_alpha),
             a_pix=float(A_), b_pix=float(B_), theta_pix=float(theta_pix), npix=npix, psf_sigma_pix=float(psf_sig),
             elong=float(A_ / B_) if np.isfinite(A_) and B_ > 0 else np.nan, neg_frac=negfrac, pos_frac=posfrac,
             sharp=sharp, sharp_psf=sharp_psf, on_cr=False, saturated=bool(target.saturated[min(yi, target.shape[0] - 1), min(xi, target.shape[1] - 1)]),
             near_bad=bool(target.bad[max(0, yi - 3):yi + 4, max(0, xi - 3):xi + 4].any()))
    if getattr(target, "cr", None) is not None:
        d["on_cr"] = bool(target.cr[max(0, yi - 1):yi + 2, max(0, xi - 1):xi + 2].any())
    ra, de = target.wcs.all_pix2world([cx], [cy], 0)
    d["ra"] = float(ra[0]); d["dec"] = float(de[0])
    # centroid error (pixels): FWHM/(2.355*SNR) -> arcsec; floor 3 mas
    sig_pix = max(psf_sig / max(pk, 1.0), 0.003 / target.pixscale)
    d["sig_pos_arcsec"] = float(sig_pix * target.pixscale)
    # sky position angle of the major axis (N through E) from pixel angle using the local WCS jacobian
    if np.isfinite(theta_pix):
        dx, dy = np.cos(theta_pix), np.sin(theta_pix)
        r2, d2 = target.wcs.all_pix2world([cx, cx + dx], [cy, cy + dy], 0)
        dra = (r2[1] - r2[0]) * np.cos(np.radians(de[0])); dde = d2[1] - d2[0]
        d["pa_deg"] = float(np.degrees(np.arctan2(dra, dde)) % 180.0)
    else:
        d["pa_deg"] = np.nan
    d["a_arcsec"] = float(A_ * target.pixscale) if np.isfinite(A_) else np.nan
    return d


def _dedupe(dets, r=1.5):
    """Merge detections of the same sign closer than r pixels (keep the highest S/N; point channel preferred)."""
    if len(dets) < 2:
        return dets
    from scipy.spatial import cKDTree
    xy = np.array([[d["x"], d["y"]] for d in dets])
    tree = cKDTree(xy)
    order = np.argsort([-d["snr"] for d in dets])
    taken = np.zeros(len(dets), bool); out = []
    for i in order:
        if taken[i]:
            continue
        grp = [j for j in tree.query_ball_point(xy[i], r) if not taken[j] and dets[j]["sign"] == dets[i]["sign"]]
        for j in grp:
            taken[j] = True
        out.append(dets[i])
    return out
