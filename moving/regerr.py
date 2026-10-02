"""Registration error of a target exposure against its template, measured from the data.

ZOGY's astrometric-error term (`zogy(astrom_n, astrom_r)`) needs the sigma of the relative registration.  Two ways were available before:
a number typed by the user, or the rms of the alignment fit of the chip (`chip.align['rms_mas']`, one value per chip against the catalogue/anchor,
and zero for the anchor exposure, although its *template* is a median of the other, misregistered exposures).

`measure_registration(target, template)` measures it for the actual pair: sources (stars and compact galaxies) bright in both images
are re-centred on the template image with the same windowed centroid (`sep.winpos`); the per-axis offsets d = x_target - x_template are the
registration error plus the centroid noise of the two images.  The centroid noise (FWHM / 2.355 / S/N per star, the template with the
S/N of the target scaled by sqrt(n_inputs)) is subtracted in quadrature:

    sigma_reg^2 = max(rms(d)^2 - median(sigma_c^2), 0)         per axis, rms^2 = median(d)^2 + (1.4826 MAD(d))^2  (offset + robust scatter, so a constant shift counts)

Returns None when fewer than `min_stars` sources are usable (the caller then falls back to the sidecar value or to no astrometric term).
It does NOT model spatial variation (one sigma per axis per exposure pair) and it measures the *stars*; the registration of the moving/faint
objects and of galaxies with a different centroid definition may be different.
"""
import numpy as np
from . import imaging as I


def _sources(target, tplf, bkg_t, bkg_r, rms_t, snr_min, size_max=6.0, max_n=600):
    """Sources bright in BOTH the target and the template (so cosmic rays and movers drop out): stars *and* compact galaxies are fine here, only the
    relative centroid matters (the same windowed centroid is applied to both images).  Returns x, y, snr, r50 (target frame, 0-based pixels)."""
    import sep
    sub = np.ascontiguousarray((target.data - bkg_t).astype(np.float32))
    o = sep.extract(sub, 5.0, err=np.ascontiguousarray(rms_t), mask=np.ascontiguousarray(target.bad), minarea=5, deblend_nthresh=32, deblend_cont=0.005)
    if len(o) == 0:
        return None
    snr = o["flux"] / (np.sqrt(o["npix"]) * float(np.median(rms_t)) + 1e-12)       # background-limited S/N (flux is in e-/s: a Poisson term would need t_exp)
    ok = (o["flag"] == 0) & (snr > snr_min) & (o["flux"] > 0)
    o = o[ok]; snr = snr[ok]
    if len(o) == 0:
        return None
    r50, _ = sep.flux_radius(sub, o["x"], o["y"], 12.0 * np.ones(len(o)), 0.5)
    good = np.isfinite(r50) & (r50 > 0.5) & (r50 < size_max)
    ny, nx = target.shape
    for j in np.nonzero(good)[0]:
        yi, xi = int(round(o["y"][j])), int(round(o["x"][j]))
        if not (8 <= xi < nx - 8 and 8 <= yi < ny - 8) or target.saturated[yi - 3:yi + 4, xi - 3:xi + 4].any() or target.cr[yi - 2:yi + 3, xi - 2:xi + 3].mean() > 0.3:
            good[j] = False
    o = o[good]; snr = snr[good]; r50 = r50[good]
    if len(o) == 0:
        return None
    # confirmation in the template: aperture flux of the same size must be a significant fraction of the target's
    subr = np.ascontiguousarray((tplf - bkg_r).astype(np.float32))
    fr, _, _ = sep.sum_circle(subr, o["x"], o["y"], 2.0 * np.maximum(r50, 1.0))
    ft, _, _ = sep.sum_circle(sub, o["x"], o["y"], 2.0 * np.maximum(r50, 1.0))
    conf = (fr > 0.4 * ft) & (fr < 2.5 * ft)
    x, y, snr, r50 = o["x"][conf], o["y"][conf], snr[conf], r50[conf]
    top = np.argsort(-snr)[:max_n]
    return x[top], y[top], snr[top], r50[top]


def measure_registration(target, tplf, tpl_bad, bkg_t, bkg_r, rms_t=None, nin_mean=3.0, min_stars=8, snr_min=10.0, max_off_pix=1.5):
    """target: Chip (cropped); tplf: template image on the target grid (NaN already replaced by 0); tpl_bad: bool mask of template NaN;
    bkg_t / bkg_r: background levels; rms_t: background rms map of the target.  Returns a dict (sigma_x, sigma_y in pixels, raw_x/y, n_stars,
    centroid_noise_pix, mean_dx/dy) or None when fewer than `min_stars` usable sources."""
    import sep
    if rms_t is None:
        _, rms_t = I.background(target.data, target.bad)
    src = None
    for sm in (snr_min, 0.6 * snr_min):                  # too few sources at the strict cut: try a lower one once
        src = _sources(target, tplf, bkg_t, bkg_r, rms_t, sm)
        if src is not None and len(src[0]) >= min_stars:
            break
    if src is None or len(src[0]) < min_stars:
        return None
    x, y, snr, r50 = src
    ny, nx = target.shape
    xi = np.clip(np.round(x).astype(int), 0, nx - 1); yi = np.clip(np.round(y).astype(int), 0, ny - 1)
    ok = np.ones(len(x), bool)
    for j in range(len(x)):
        sl = (slice(max(yi[j] - 5, 0), yi[j] + 6), slice(max(xi[j] - 5, 0), xi[j] + 6))
        if tpl_bad[sl].any() or target.bad[sl].any():
            ok[j] = False
    if ok.sum() < min_stars:
        return None
    x, y, snr, r50 = x[ok], y[ok], snr[ok], r50[ok]
    sw = np.maximum(1.2, 1.2 * r50)                                       # winpos sigma ~ source size
    sub_t = np.ascontiguousarray((target.data - bkg_t).astype(np.float32)); sub_r = np.ascontiguousarray((tplf - bkg_r).astype(np.float32))
    xt, yt, ft = sep.winpos(sub_t, x, y, sw); xr, yr, fr = sep.winpos(sub_r, x, y, sw)
    good = (ft == 0) & (fr == 0) & np.isfinite(xt + yt + xr + yr)
    dx = (xt - xr)[good]; dy = (yt - yr)[good]; sg = (sw / np.maximum(snr, 1.0))[good]       # centroid sigma ~ size / (S/N)
    keep = np.hypot(dx, dy) < max_off_pix                                  # a gross offset is a blend / variable / wrong match, not registration
    dx, dy, sg = dx[keep], dy[keep], sg[keep]
    if len(dx) < min_stars:
        return None
    s_c2 = float(np.median(sg ** 2)) * (1.0 + 1.0 / max(nin_mean, 1.0))   # target + template centroid noise
    out = {}
    for ax, d in (("x", dx), ("y", dy)):
        m0 = float(np.median(d)); mad = 1.4826 * float(np.median(np.abs(d - m0)))
        rr = float(np.sqrt(m0 ** 2 + mad ** 2))            # systematic offset + scatter (a plain MAD about zero would over-estimate a constant shift by 1.48)
        out["raw_" + ax] = rr
        out["sigma_" + ax] = float(np.sqrt(max(rr ** 2 - s_c2, 0.0)))
    out.update(n_stars=int(len(dx)), centroid_noise_pix=float(np.sqrt(s_c2)), mean_dx=float(np.mean(dx)), mean_dy=float(np.mean(dy)))
    return out
