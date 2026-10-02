"""Merger / interaction indicators on top of the extended morphology (morphometry.extended): rotational asymmetry, multiple-peak (double nucleus) test, Lotz et al. 2008 Gini-M20
classification, shell and tidal-tail detection in the residual of an azimuthally-median elliptical model, and catalog pair classification.

All images are background subtracted; (x, y) are 0-based; `q, theta` the aperture ellipse (theta in radians, major axis angle from +x); `rp` the Petrosian radius along the major axis."""
import math

import numpy as np
from scipy import ndimage as ndi

FLAGS = dict(MERGER_GM20=1, MERGER_ASYM=2, DOUBLE=4, SHELL=8, TAIL=16, PAIR_MINOR=32, PAIR_MAJOR=64)


def _ell_coords(shape, x, y, q, theta):
    ny, nx = shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    dx, dy = xx - x, yy - y
    u = dx * math.cos(theta) + dy * math.sin(theta)
    v = -dx * math.sin(theta) + dy * math.cos(theta)
    return u, v / max(q, 0.05), np.sqrt(u ** 2 + (v / max(q, 0.05)) ** 2)


# ------------------------------------------------------------------------------------------------------------------------------------------ Gini-M20 class
def lotz_class(gini, m20):
    """Lotz et al. 2008 (eqs. 4-6): 'merger' above G = -0.14 M20 + 0.33; 'E/S0/Sa' above G = 0.14 M20 + 0.80; else 'Sb-Irr'."""
    if not (np.isfinite(gini) and np.isfinite(m20)):
        return 'unknown'
    if gini > -0.14 * m20 + 0.33:
        return 'merger'
    return 'E/S0/Sa' if gini > 0.14 * m20 + 0.80 else 'Sb-Irr'


# ------------------------------------------------------------------------------------------------------------------------------------------ asymmetry
def _rot180(img, xc, yc):
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    return ndi.map_coordinates(img, [2 * yc - yy, 2 * xc - xx], order=1, mode='constant', cval=np.nan)


def _asym_raw(img, ap, xc, yc, norm_img=None):
    r = _rot180(img, xc, yc)
    rm = _rot180(ap.astype(float), xc, yc) > 0.99
    use = ap & rm & np.isfinite(r)
    s = np.abs((img if norm_img is None else norm_img)[use]).sum()
    if use.sum() < 20 or s <= 0:
        return float('nan')
    return float(np.abs(img[use] - r[use]).sum() / (2.0 * s))


def asymmetry(cut, ap, x, y, noise, max_shift=2.0, step=0.5):
    """Conselice (2003) rotational asymmetry in the aperture `ap`: A = min over centres of sum|I - I_180| / (2 sum|I|) minus the same quantity of a noise image (A_bkg).  -> (A, A_raw, A_bkg)."""
    shifts = np.arange(-max_shift, max_shift + 1e-9, step)

    def minimise(img, norm_img=None):
        best = float('inf')
        for dy in shifts:
            for dx in shifts:
                a = _asym_raw(img, ap, x + dx, y + dy, norm_img)
                if np.isfinite(a) and a < best:
                    best = a
        return best

    best = minimise(cut)
    if not np.isfinite(best):
        return float('nan'), float('nan'), float('nan')
    ab = minimise(noise, cut)                                     # the background term is minimised over the same centres, so that the selection bias cancels
    return best - (ab if np.isfinite(ab) else 0.0), best, ab


# ------------------------------------------------------------------------------------------------------------------------------------------ multiple peaks
def peaks(cut, ap, rms, sigma_s, x, y, min_frac=0.15, nsig=5.0):
    """Local maxima of the smoothed image inside `ap` that are >= nsig sigma_smoothed and >= min_frac of the brightest one, at least 2.5 sigma_s apart.  -> list of (value, px, py), brightest first."""
    sm = ndi.gaussian_filter(np.where(ap, cut, 0.0), sigma_s)
    nz = rms / (2.0 * math.sqrt(math.pi) * sigma_s)
    mx = ndi.maximum_filter(sm, footprint=np.ones((int(2 * max(2, round(1.5 * sigma_s)) + 1),) * 2))
    cand = (sm == mx) & ap & (sm >= nsig * nz)
    ys, xs = np.nonzero(cand)
    if len(ys) == 0:
        return []
    vals = sm[ys, xs]
    order = np.argsort(-vals)
    out = []
    for k in order:
        if vals[k] < min_frac * vals[order[0]]:
            break
        if all(math.hypot(xs[k] - p[1], ys[k] - p[2]) >= 2.5 * sigma_s for p in out):
            out.append((float(vals[k]), float(xs[k]), float(ys[k])))
    return out


# ------------------------------------------------------------------------------------------------------------------------------------------ shells / tails
def symmetric_model(cut, mask, x, y, q, theta, rmax, dr=1.0, smooth=2):
    """Azimuthal-median radial profile in elliptical annuli (robust against localized features), linearly interpolated back onto the pixels."""
    u, v, r = _ell_coords(cut.shape, x, y, q, theta)
    edges = np.arange(0.0, rmax + 2 * dr, dr)
    idx = np.digitize(r.ravel(), edges) - 1
    ok = (~mask).ravel() & np.isfinite(cut.ravel())
    prof = np.full(len(edges) - 1, np.nan)
    rc = 0.5 * (edges[1:] + edges[:-1])
    flat = cut.ravel()
    for k in range(len(prof)):
        s = ok & (idx == k)
        if s.sum() >= 8:
            prof[k] = np.median(flat[s])
    good = np.isfinite(prof)
    if good.sum() < 4:
        return np.full(cut.shape, np.nan), r
    prof = np.interp(rc, rc[good], prof[good])
    if smooth:
        prof = ndi.gaussian_filter1d(prof, smooth, mode='nearest')
    return np.interp(r, rc, prof), r


def detect_features(cut, mask, x, y, q, theta, rp, rms, nsig=3.0, smooth_px=None, r_lo=0.4, r_hi=4.0, min_area=25, shell_dphi=45.0, shell_aspect=1.8, tail_elong=2.5, tail_len=0.8, max_frac=0.25, flux_p=None):
    """Residual features of one galaxy.  The galaxy is modelled by the azimuthal median in elliptical annuli (centre / q / theta given); the residual smoothed by a Gaussian (sigma `smooth_px`, default max(1.5, 0.08 R_P))
    is thresholded at `nsig` sigma and the connected positive components between r_lo and r_hi R_P are classified by their geometry about the galaxy centre:
      shell / arc : angular extent >= shell_dphi deg, tangential/radial size ratio >= shell_aspect, direction tangential (> 55 deg from radial), r in [0.5, r_hi] R_P;
      tail        : elongation (sigma_major/sigma_minor) >= tail_elong, extending beyond 1.1 R_P with length >= tail_len R_P, direction within 50 deg of the radial one at its centroid, starting inside 1.6 R_P;
      other       : the rest (companions, clumps, spiral-arm pieces); components with more than `max_frac` of the galaxy flux are companions, never shells / tails.
    -> dict(n_shell, n_tail, n_other, shell_frac, tail_frac, tail_len, tidal_frac, res_rms_frac, components=[...]) fractions relative to flux_p (or the model flux inside r_hi R_P)."""
    out = dict(n_shell=0, n_tail=0, n_other=0, shell_frac=0.0, tail_frac=0.0, tail_len=0.0, tidal_frac=0.0, res_frac=float('nan'), components=[])
    if not (np.isfinite(rp) and rp >= 2.0):
        return out
    ny, nx = cut.shape
    rmax = min(r_hi * rp + 6, 0.5 * math.hypot(nx, ny))
    model, r = symmetric_model(cut, mask, x, y, q, theta, rmax)
    if not np.isfinite(model).any():
        return out
    res = np.where(mask | ~np.isfinite(cut), 0.0, cut - model)
    s = smooth_px or max(1.5, 0.08 * rp)
    sm = ndi.gaussian_filter(res, s)
    wt = ndi.gaussian_filter((~mask & np.isfinite(cut)).astype(float), s)
    sm = np.where(wt > 0.3, sm / np.maximum(wt, 1e-3), 0.0)
    sig = rms / (2.0 * math.sqrt(math.pi) * s)
    u, v, rr = _ell_coords(cut.shape, x, y, q, theta)
    zone = (rr >= r_lo * rp) & (rr <= r_hi * rp) & ~mask
    # empirical noise of the smoothed residual: rms of its negative half in the zone (a feature is a positive excess; correlated pixel noise, clumps and arms enter the negative half as well)
    neg = sm[zone & (sm < 0)]
    sig_neg = float(math.sqrt(np.mean(neg ** 2))) if neg.size > 50 else sig
    out['sigma_smoothed'] = sig = max(sig, sig_neg)
    cand = (sm > nsig * sig) & zone
    lab, nl = ndi.label(cand, structure=np.ones((3, 3)))
    mflux = float(np.nansum(np.where(rr <= r_hi * rp, model, 0.0))) if flux_p is None else float(flux_p)
    mflux = max(mflux, 1e-30)
    out['res_frac'] = float(np.sum(np.abs(res)[zone]) / mflux)
    phi_all = np.arctan2(v, u)
    near_mask = ndi.distance_transform_edt(~mask) <= (s + 3.0) if mask.any() else np.zeros(mask.shape, bool)
    for k in range(1, nl + 1):
        sel = lab == k
        area = int(sel.sum())
        if area < min_area:
            continue
        w = np.clip(sm[sel], 0, None)
        ys, xs = np.nonzero(sel)
        rk, ph = rr[sel], phi_all[sel]
        # angular extent about the centre: complement of the largest gap
        pa = np.sort(ph)
        gaps = np.diff(np.concatenate([pa, [pa[0] + 2 * math.pi]]))
        dphi = math.degrees(2 * math.pi - gaps.max())
        # second moments in the (u, v/q) plane: radial and tangential spread at the centroid
        uu, vv = u[sel], v[sel]
        wsum = w.sum()
        cu, cv = (w * uu).sum() / wsum, (w * vv).sum() / wsum
        cov = np.cov(np.vstack([uu - cu, vv - cv]), aweights=w + 1e-12)
        ev, evec = np.linalg.eigh(cov)
        smin, smaj = math.sqrt(max(ev[0], 1e-9)), math.sqrt(max(ev[1], 1e-9))
        major = evec[:, 1]
        rad = np.array([cu, cv]) / max(math.hypot(cu, cv), 1e-9)
        ang_rad = math.degrees(math.acos(min(1.0, abs(float(major @ rad)))))          # 0 = major axis radial, 90 = tangential
        rmean = float((w * rk).sum() / wsum)
        flux = float(res[sel].sum())
        length = 4.0 * smaj
        comp = dict(area=area, r_mean=rmean / rp, r_max=float(rk.max()) / rp, dphi=dphi, elong=smaj / smin, angle_to_radial=ang_rad, flux_frac=flux / mflux, length=length / rp, kind='other',
                    x=float(xs.mean()), y=float(ys.mean()))
        tangential = ang_rad > 55.0
        if abs(flux) / mflux > max_frac or near_mask[sel].mean() > 0.2:       # companions / halos of masked neighbours
            pass
        elif dphi >= shell_dphi and tangential and smaj / smin >= shell_aspect and 0.5 <= rmean / rp <= r_hi:
            comp['kind'] = 'shell'
        elif smaj / smin >= tail_elong and comp['r_max'] >= 1.1 and length / rp >= tail_len and ang_rad <= 50.0 and float(rk.min()) <= 1.6 * rp:
            comp['kind'] = 'tail'
        out['components'].append(comp)
        out['n_' + comp['kind']] += 1
        if comp['kind'] == 'shell':
            out['shell_frac'] += comp['flux_frac']
        elif comp['kind'] == 'tail':
            out['tail_frac'] += comp['flux_frac']
            out['tail_len'] = max(out['tail_len'], comp['length'])
        out['tidal_frac'] += comp['flux_frac']
    return out


# ------------------------------------------------------------------------------------------------------------------------------------------ pairs
def find_pairs(x, y, rp, flux, sep_max_rp=1.5, major=0.25, minor=0.1):
    """Catalog pair classification.  x, y, rp, flux: arrays (rp may contain nan -> A_IMAGE-like fallback by the caller).  For object i the companions j with separation < sep_max_rp (R_P,i + R_P,j)
    and flux ratio f_j/f_i in [minor, 1/minor] are found; 'major' if the ratio is in [major, 1/major].
    -> dict of arrays n_comp, sep (px, nearest companion), ratio (companion/target flux), kind (0 none, 1 minor, 2 major), touching (apertures overlap: sep < R_P,i + R_P,j)."""
    x, y, rp, flux = [np.asarray(a, float) for a in (x, y, rp, flux)]
    n = len(x)
    out = dict(n_comp=np.zeros(n, int), sep=np.full(n, np.nan), ratio=np.full(n, np.nan), kind=np.zeros(n, int), touching=np.zeros(n, bool))
    if n < 2:
        return out
    from scipy.spatial import cKDTree
    ok = np.isfinite(x) & np.isfinite(y)
    tree = cKDTree(np.where(ok[:, None], np.column_stack([x, y]), -1e9))
    rmax = np.nanmax(rp[np.isfinite(rp)]) if np.isfinite(rp).any() else 10.0
    for i in np.nonzero(ok)[0]:
        for j in tree.query_ball_point([x[i], y[i]], sep_max_rp * (np.nan_to_num(rp[i], nan=5.0) + rmax)):
            if j == i or not ok[j]:
                continue
            d = math.hypot(x[i] - x[j], y[i] - y[j])
            rs = np.nan_to_num(rp[i], nan=5.0) + np.nan_to_num(rp[j], nan=5.0)
            if d > sep_max_rp * rs or not (flux[i] > 0 and flux[j] > 0):
                continue
            fr = flux[j] / flux[i]
            if not (minor <= fr <= 1.0 / minor):
                continue
            out['n_comp'][i] += 1
            if not np.isfinite(out['sep'][i]) or d < out['sep'][i]:
                out['sep'][i], out['ratio'][i] = d, fr
                out['touching'][i] = d < rs
            out['kind'][i] = max(out['kind'][i], 2 if major <= fr <= 1.0 / major else 1)
    return out
