"""Image depth: local noise-law based limiting magnitude and surface-brightness maps (numpy / scipy / sep).

The noise of an aperture sum in a correlated (drizzled / resampled) image is sigma_N = sigma_pix * alpha * N^beta (ogfkit.noise).  The law
is measured once on blank apertures of the whole image; the pixel rms varies over the image (local rms map: sep mesh, or a user rms / weight
map).  Hence

    sigma_ap(x, y)  = rms(x, y) * f_ap(r)     -> m_lim(x, y)  = zp - 2.5 log10(nsig * sigma_ap) - apcorr
    sigma_box(x, y) = rms(x, y) * f_box(L)    -> mu_lim(x, y) = zp - 2.5 log10(nsig * sigma_box / area_arcsec2)

and the per-tile *empirical* blank-aperture / blank-box scatter is compared with the model (`tile_table`).
"""
def _import_ogfmeas():
    """In-tree measurements. Finds ogfmeas from this file so a script does not need PYTHONPATH."""
    import pathlib
    import sys
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas
    return ogfmeas.measurement_library()

import math

import numpy as np
from scipy import ndimage as ndi

from . import noise as nz


def prepare(data, bad=None, bw=64, thresh=2.0, minarea=8, grow=3.0):
    """Source mask, mesh background / rms and the background-subtracted image."""
    data = np.asarray(data, np.float32)
    mask = nz.source_mask(data, bad, thresh=thresh, minarea=minarea, grow=grow, bw=bw)
    d = np.where(np.isfinite(data), data, 0.0)
    bkg, rms = nz.mesh_background(d, mask, bw, 3)
    return dict(mask=mask, bkg=bkg, rms=rms, sub=d - bkg,
                bad=(np.zeros(data.shape, bool) if bad is None else np.asarray(bad, bool)) | ~np.isfinite(data))


def rms_from_file(rms_map, valid, mesh_rms, kind='rms'):
    """User rms map (or weight map: rms ~ 1/sqrt(w)) scaled to the measured noise: the median ratio mesh_rms / user over `valid` pixels."""
    r = np.asarray(rms_map, np.float64)
    if kind == 'weight':
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(r > 0, 1.0 / np.sqrt(r), np.nan)
    ok = valid & np.isfinite(r) & (r > 0) & np.isfinite(mesh_rms) & (mesh_rms > 0)
    if not ok.any():
        raise ValueError('rms / weight map has no valid pixel')
    r = r * float(np.median(mesh_rms[ok] / r[ok]))
    return np.where(np.isfinite(r) & (r > 0), r, np.nan)


def blank_circle_positions(mask, r, n, rng, box=None, max_tries=40):
    """Random aperture centres (x, y) whose circle of radius r has no masked pixel; box = (x0, y0, x1, y1) restricts the centres (shared implementation: noise.blank_positions)."""
    return nz.blank_positions(mask, r, n, rng, box=box, max_tries=max_tries, kmin=100)


def circle_sums(sub, pos, r):
    sep = _import_ogfmeas()
    d = np.ascontiguousarray(sub, np.float64)
    s, _, _ = sep.sum_circle(d, pos[:, 0], pos[:, 1], r, subpix=5)
    return s


def _sample(m, x, y):
    ny, nx = m.shape
    return m[np.clip(np.round(y).astype(int), 0, ny - 1), np.clip(np.round(x).astype(int), 0, nx - 1)]


def circle_factor(sub, mask, rms, r, n=1500, seed=1, box=None):
    """sigma_N(r) / local pixel rms from blank apertures: clipped std of (aperture sum / rms at the aperture centre).
    Returns (factor, n_used).  Using the local rms makes the factor independent of the noise gradient."""
    rng = np.random.default_rng(seed)
    pos = blank_circle_positions(mask, r, n, rng, box=box)
    if len(pos) < 30:
        return float('nan'), int(len(pos))
    s = circle_sums(sub, pos, r)
    loc = _sample(rms, pos[:, 0], pos[:, 1])
    ok = np.isfinite(loc) & (loc > 0)
    return float(nz.clipped_std(s[ok] / loc[ok], 4.0)), int(ok.sum())


def aperture_law(sub, mask, rms, radii, n_aper=800, seed=1):
    """Noise law sigma_N = rms * alpha * N^beta (N = pi r^2) fitted to circle_factor at several radii; also the raw table."""
    N, S, U = [], [], []
    for r in radii:
        f, n = circle_factor(sub, mask, rms, r, n_aper, seed)
        if np.isfinite(f) and f > 0:
            N.append(math.pi * r * r); S.append(f); U.append(n)
    law = nz.fit_noise_law(np.array(N), np.array(S), 1.0) if len(N) >= 3 else dict(alpha=1.0, beta=0.5, ok=False)
    law.update(radii=list(radii), N=N, sigma=S, n_used=U)
    return law


def box_sums(sub, mask, L, stride=None):
    """Sums of L x L px boxes on a regular grid (stride default L/2) whose masked fraction is below 5 %.  Returns (sums, x_centres, y_centres)."""
    L = int(L)
    ny, nx = sub.shape
    if L >= min(ny, nx) // 2:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    d = np.where(mask, 0.0, sub)
    cs = ndi.uniform_filter(d, L, mode='constant') * L * L
    mf = ndi.uniform_filter(mask.astype(float), L, mode='constant')
    stride = max(1, int(stride or L // 2))
    h = L // 2
    ys = np.arange(h, ny - h, stride)
    xs = np.arange(h, nx - h, stride)
    YY, XX = np.meshgrid(ys, xs, indexing='ij')
    ok = mf[YY, XX] < 0.05
    return cs[YY, XX][ok], XX[ok].astype(float), YY[ok].astype(float)


def box_factor(sub, mask, rms, L):
    """sigma_box(L x L px) / local pixel rms from blank boxes (NaN when fewer than 12 independent-ish boxes)."""
    s, x, y = box_sums(sub, mask, L)
    if len(s) < 12:
        return float('nan'), int(len(s))
    loc = _sample(rms, x, y)
    ok = np.isfinite(loc) & (loc > 0)
    return float(nz.clipped_std(s[ok] / loc[ok], 4.0)), int(ok.sum())


def box_law(sub, mask, rms, sizes):
    """sigma_box = rms * gamma * area^delta from blank boxes of the given sizes (px)."""
    N, S = [], []
    for L in sizes:
        f, n = box_factor(sub, mask, rms, L)
        if np.isfinite(f) and n >= 12:
            N.append(float(L * L)); S.append(f)
    law = nz.fit_noise_law(np.array(N), np.array(S), 1.0) if len(N) >= 3 else dict(alpha=1.0, beta=0.5, ok=False)
    return dict(gamma=law['alpha'], delta=law['beta'], ok=law.get('ok', False), N=N, sigma=S)


def mag_limit_map(rms, factor, zp, nsig=5.0, apcorr=0.0):
    """nsig point-source limiting magnitude: zp - 2.5 log10(nsig * rms * factor) - apcorr, factor = sigma_N(r) / rms from circle_factor."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return zp - 2.5 * np.log10(nsig * rms * factor) - apcorr


def sb_limit_map(rms, factor, L_px, pixscale, zp, nsig=3.0):
    """nsig surface-brightness limit (mag / arcsec^2) of an L x L px box: mu = zp - 2.5 log10(nsig * rms * factor / area_arcsec2)."""
    area = (L_px * pixscale) ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        return zp - 2.5 * np.log10(nsig * rms * factor / area)


def tile_table(sub, mask, rms, f_ap, f_box, tile, r_ap, L_box, n_aper=200, seed=1, min_aper=40):
    """Per-tile empirical blank-aperture scatter (clipped std of circle sums) and blank-box scatter versus the model prediction
    from the local mesh rms.  Returns a list of dicts (tile index, centre, n, sigma_emp, err, sigma_model, ratio, ...)."""
    ny, nx = sub.shape
    rng = np.random.default_rng(seed)
    nty, ntx = max(1, ny // tile), max(1, nx // tile)
    out = []
    sb, bx, by = box_sums(sub, mask, L_box) if L_box else (np.zeros(0), np.zeros(0), np.zeros(0))
    for j in range(nty):
        for i in range(ntx):
            x0, x1 = i * nx // ntx, (i + 1) * nx // ntx
            y0, y1 = j * ny // nty, (j + 1) * ny // nty
            pos = blank_circle_positions(mask, r_ap, n_aper, rng, box=(x0, y0, x1 - 1, y1 - 1))
            rec = dict(tile=j * ntx + i, x0=x0, x1=x1, y0=y0, y1=y1, xc=0.5 * (x0 + x1), yc=0.5 * (y0 + y1))
            rr = rms[y0:y1, x0:x1]
            gm = ~mask[y0:y1, x0:x1] & np.isfinite(rr)
            rec['rms_pix'] = float(np.median(rr[gm])) if gm.any() else float('nan')
            rec['valid_frac'] = float(gm.mean())
            rec['n_aper'] = int(len(pos))
            if len(pos) >= min_aper:
                s = circle_sums(sub, pos, r_ap)
                rec['sigma_emp'] = float(nz.clipped_std(s, 4.0))
                rec['sigma_emp_err'] = rec['sigma_emp'] / math.sqrt(2.0 * len(pos))
            else:
                rec['sigma_emp'] = rec['sigma_emp_err'] = float('nan')
            rec['sigma_model'] = rec['rms_pix'] * f_ap
            rec['ratio'] = rec['sigma_emp'] / rec['sigma_model'] if rec['sigma_model'] > 0 else float('nan')
            if len(sb):
                k = (bx >= x0) & (bx < x1) & (by >= y0) & (by < y1)
                rec['n_box'] = int(k.sum())
                if k.sum() >= 12:
                    rec['box_sigma_emp'] = float(nz.clipped_std(sb[k], 4.0))
                    rec['box_sigma_model'] = rec['rms_pix'] * f_box
                    rec['box_ratio'] = rec['box_sigma_emp'] / rec['box_sigma_model']
            out.append(rec)
    return out


def area_depth(depth, valid, edges=None):
    """Cumulative fraction of the valid area that is deeper than each magnitude: list of (mag, fraction)."""
    v = depth[valid & np.isfinite(depth)]
    if not v.size:
        return []
    if edges is None:
        edges = np.arange(np.floor(v.min() * 4) / 4, v.max() + 0.25, 0.25)
    srt = np.sort(v)
    return [(float(e), float(1.0 - np.searchsorted(srt, e, side='left') / len(srt))) for e in edges]


def region_labels(kind, shape, depth=None, valid=None, grid=(3, 3), classes=4, label_file=None):
    """Integer label image (0 = no region): 'grid' (nx x ny rectangles), 'rms' / 'depth' quantile classes of the depth map, or a label image."""
    ny, nx = shape
    lab = np.zeros(shape, np.int32)
    if kind == 'grid':
        gx, gy = grid
        ix = np.minimum(np.arange(nx) * gx // nx, gx - 1)
        iy = np.minimum(np.arange(ny) * gy // ny, gy - 1)
        lab = (iy[:, None] * gx + ix[None, :] + 1).astype(np.int32)
    elif kind in ('rms', 'depth'):
        ok = valid & np.isfinite(depth)
        q = np.quantile(depth[ok], np.linspace(0, 1, classes + 1))
        q[0] -= 1e-6
        q[-1] += 1e-6
        lab[ok] = np.clip(np.digitize(depth[ok], q[1:-1]) + 1, 1, classes)
    elif kind == 'file':
        lab = np.asarray(label_file).astype(np.int32)
        if lab.shape != shape:
            raise ValueError('region map shape differs from the image')
    else:
        raise ValueError('unknown region kind %s' % kind)
    if valid is not None:
        lab = np.where(valid, lab, 0)
    return lab
