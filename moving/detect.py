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


def difference_chip(target, template, nin=None, method="zogy", snr_det=5.0, psf_size=25,
                    min_pix=3, extra_mask=None, psf_t=None, fw_t=None, n_t=0, psf_r=None, fw_r=None, n_r=0):
    """Returns dict(alpha=..., S=..., sigma_alpha=..., score=..., mask=..., dets=[...], info=...)."""
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
        out = Z.zogy(N, R, psf_t, psf_r, sn, sr, Fn=Fr_ratio, Fr=1.0)
        S = out["S"]; alpha = out["alpha_new"]; sig_alpha = out["sigma_alpha_new"]       # flux in the scale of the target exposure
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
    info["n_pos"] = sum(1 for d in dets if d["sign"] > 0); info["n_neg"] = sum(1 for d in dets if d["sign"] < 0)
    return dict(alpha=alpha, score=score, mask=inval, dets=dets, info=info, sigma_alpha=float(sig_alpha),
                psf=psf_t, bkg=bkg_t, rms=rms_t)


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
