"""Stacking of cutouts and mean radial profiles (no GUI, numpy / scipy only).

Pipeline (all steps optional, see `StackConfig`):
  positions -> aligned cutouts (sub-pixel resampling about the true position, optional rotation by the position angle and
  scaling by a per-source size) -> per-cutout background -> masks (input mask, catalog neighbours, per-cutout detection)
  -> normalisation / weights -> combination (mean, median, weighted mean, sigma-clipped mean) -> bootstrap errors
  (resampling the *sources*) -> radial profile / aperture flux of the stack.  `null_positions` draws blank-sky positions that go
  through exactly the same pipeline (null stack) to measure the residual background / mask bias of a stack.

Pixel coordinates in this module are 0-based (x = column).  Cutouts are (2*half+1)^2, the target at index [half, half].
Masked or outside-the-image pixels are NaN in the cube.  Resampling conserves surface brightness (not flux).
"""
import math

import numpy as np

METHODS = ('mean', 'median', 'wmean', 'clipmean')


class StackConfig(dict):
    DEFAULTS = dict(
        half=30,              # cutout half size (output px; units of `scale` when rescaled)
        align='sub',          # 'sub' = spline resampling on the true position, 'int' = integer shift (no resampling)
        order=3,              # spline order of the resampling
        bkg='annulus',        # 'none' | 'annulus' | 'image' (needs bkg_map)
        bkg_in=0.7, bkg_out=1.0,   # annulus as fractions of `half`
        norm='none',          # 'none' | 'flux' (aperture flux within norm_r) | 'peak' | 'column' (norm_values)
        norm_r=5.0,
        min_valid=0.6,        # minimum fraction of valid pixels in the cutout
        core_r=2.0,           # the central core_r px must be valid
        mask_detect=0.0,      # >0: mask pixels > mask_detect sigma (smoothed) in each cutout, outside protect_r
        detect_grow=2.0, detect_minarea=4, protect_r=3.0, detect_smooth=1.5,
        mask_cat_scale=0.0,   # >0: mask other catalog objects with radius max(mask_cat_min, scale * extent)
        mask_cat_min=3.0,
        method='mean', clip=3.0, weights='none',   # weights: 'none' | 'invvar' | 'column'
        n_boot=200, seed=1,
        q=1.0, pa=0.0,        # ellipse of the radial profile (pa in deg from +x of the cutout)
        ap_r=5.0,             # aperture radius of the stack flux (output px)
        estimator='mean',     # radial profile estimator of the stack: 'mean' | 'median'
    )

    def __init__(self, **kw):
        super().__init__(self.DEFAULTS)
        bad = set(kw) - set(self.DEFAULTS)
        if bad:
            raise ValueError('unknown stacking option(s): %s' % ', '.join(sorted(bad)))
        self.update(kw)

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k)


# ---------------------------------------------------------------------------------------------------------------- cutouts
def _sigclip_median(v, nsig=3.0, iters=3):
    v = v[np.isfinite(v)]
    if v.size < 5:
        return float('nan'), float('nan')
    for _ in range(iters):
        m = np.median(v)
        s = 1.4826 * np.median(np.abs(v - m))
        if s <= 0:
            break
        k = np.abs(v - m) < nsig * s
        if k.all() or k.sum() < 5:
            break
        v = v[k]
    m = float(np.median(v))
    return m, float(1.4826 * np.median(np.abs(v - m)))


def _radial_median_model(img, valid):
    """Azimuthal median about the cutout centre in 1-px bins (robust against the neighbours that are to be found)."""
    S = img.shape[0]
    half = S // 2
    y, x = np.indices(img.shape) - half
    ib = np.hypot(x, y).astype(int)
    model = np.zeros(img.shape)
    for k in range(ib.max() + 1):
        m = (ib == k) & valid
        if m.sum() >= 3:
            model[ib == k] = np.median(img[m])
        elif k:
            model[ib == k] = model[ib == k - 1][0] if (ib == k - 1).any() else 0.0
    return model


def _detect_mask(img, valid, nsig, grow, minarea, protect, smooth):
    """Pixels of other sources: smoothed residual of the cutout about its azimuthal median profile above nsig sigma
    (so that the target and its symmetric halo are not detected), grown; nothing within `protect` px of the centre."""
    from scipy import ndimage as ndi
    res = img - _radial_median_model(img, valid)
    fill = np.where(valid, res, 0.0)
    _, sig = _sigclip_median(res[valid])
    if not np.isfinite(sig) or sig <= 0:
        return np.zeros(img.shape, bool)
    sm = ndi.gaussian_filter(fill, smooth)
    sig_s = sig / (2.0 * math.sqrt(math.pi) * smooth)     # white-noise rms of the smoothed image
    w = ndi.gaussian_filter(valid.astype(float), smooth)
    det = (sm > nsig * sig_s) & (w > 0.5)
    if minarea > 1 and det.any():
        lab, n = ndi.label(det)
        sizes = ndi.sum(det, lab, index=np.arange(1, n + 1))
        keep = np.zeros(n + 1, bool)
        keep[1:] = sizes >= minarea
        det = keep[lab]
    if grow > 0 and det.any():
        det = ndi.binary_dilation(det, structure=ndi.generate_binary_structure(2, 1), iterations=int(math.ceil(grow)))
    half = img.shape[0] // 2
    yy, xx = np.indices(img.shape)
    det &= np.hypot(xx - half, yy - half) > protect
    return det


def make_cutouts(data, bad, xs, ys, cfg, scales=None, angles=None, bkg_map=None, norm_values=None, weight_values=None,
                 neigh=None, self_index=None):
    """-> dict(cube (N,S,S) float32 with NaN, ok (N,) bool, reason list, bkg, rms, norm, weight, valid_frac).

    xs, ys   0-based positions; scales: output pixel size in input pixels (default 1); angles: position angle (deg from +x towards +y,
             e.g. THETA_IMAGE) of each source; the output grid is rotated by it so that the major axis lies along +x of the cutout.
    neigh    (M,3) array of (x, y, radius) catalog objects to mask in every cutout, self_index = index of each target in it.
    """
    from scipy import ndimage as ndi
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    N = len(xs)
    half = int(cfg.half)
    S = 2 * half + 1
    scales = np.ones(N) if scales is None else np.asarray(scales, float)
    angles = np.zeros(N) if angles is None else np.asarray(angles, float)
    ny, nx = data.shape
    cube = np.full((N, S, S), np.nan, np.float32)
    ok = np.zeros(N, bool)
    reason = [''] * N
    bk = np.full(N, np.nan)
    rms = np.full(N, np.nan)
    nrm = np.ones(N)
    wgt = np.ones(N)
    vf = np.zeros(N)
    u = np.arange(S) - half
    U, V = np.meshgrid(u, u)
    R = np.hypot(U, V)
    ring = (R >= cfg.bkg_in * half) & (R <= cfg.bkg_out * half)
    bad_all = np.zeros(data.shape, bool) if bad is None else bad
    use_data = np.where(np.isfinite(data), data, 0.0).astype(np.float64)
    base_bad = bad_all | ~np.isfinite(data)
    for i in range(N):
        sc, th = scales[i], math.radians(angles[i])
        c, s = math.cos(th), math.sin(th)
        if cfg.align == 'int':
            sc, c, s = 1.0, 1.0, 0.0
            xc, yc = float(round(xs[i])), float(round(ys[i]))
            xi, yi = xc + U, yc + V
            ext = half + 1
        else:
            xc, yc = xs[i], ys[i]
            xi = xc + sc * (c * U - s * V)
            yi = yc + sc * (s * U + c * V)
            ext = int(math.ceil(half * sc * 1.5)) + 1
        if not (np.isfinite(xc) and np.isfinite(yc)):
            reason[i] = 'bad position'
            continue
        x0, x1 = int(math.floor(xi.min())) - 3, int(math.ceil(xi.max())) + 4
        y0, y1 = int(math.floor(yi.min())) - 3, int(math.ceil(yi.max())) + 4
        if x1 < 0 or y1 < 0 or x0 >= nx or y0 >= ny:
            reason[i] = 'outside'
            continue
        W, H = x1 - x0, y1 - y0
        win = np.zeros((H, W))
        wb = np.ones((H, W), bool)                      # outside the image = bad
        sx0, sx1, sy0, sy1 = max(x0, 0), min(x1, nx), max(y0, 0), min(y1, ny)
        win[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = use_data[sy0:sy1, sx0:sx1]
        wb[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = base_bad[sy0:sy1, sx0:sx1]
        if bkg_map is not None and cfg.bkg == 'image':
            wbk = np.zeros((H, W))
            wbk[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = bkg_map[sy0:sy1, sx0:sx1]
            win = win - wbk
        if neigh is not None and cfg.mask_cat_scale > 0 and len(neigh):
            sel = np.where((neigh[:, 0] > x0 - neigh[:, 2]) & (neigh[:, 0] < x1 + neigh[:, 2]) &
                           (neigh[:, 1] > y0 - neigh[:, 2]) & (neigh[:, 1] < y1 + neigh[:, 2]))[0]
            gy, gx = np.indices((H, W))
            gx = gx + x0
            gy = gy + y0
            for j in sel:
                if self_index is not None and self_index[i] == j:
                    continue
                if self_index is None and math.hypot(neigh[j, 0] - xc, neigh[j, 1] - yc) < 2.0:
                    continue
                wb |= np.hypot(gx - neigh[j, 0], gy - neigh[j, 1]) <= neigh[j, 2]
        good = ~wb
        if good.sum() < 10:
            reason[i] = 'no valid pixels'
            continue
        fillv = float(np.median(win[good]))
        win = np.where(good, win, fillv)
        co = np.array([yi - y0, xi - x0])
        if cfg.align == 'int':
            cut = ndi.map_coordinates(win, co, order=0, mode='nearest')
            bm = ndi.map_coordinates(wb.astype(float), co, order=0, mode='nearest') > 0.5
        else:
            cut = ndi.map_coordinates(win, co, order=int(cfg.order), mode='nearest', prefilter=cfg.order > 1)
            bm = ndi.map_coordinates(wb.astype(float), co, order=1, mode='nearest') > 1e-3
        valid = ~bm
        if cfg.mask_detect > 0:
            # background first (so that the threshold is relative to the local level)
            m0, _ = _sigclip_median(cut[valid & ring])
            m0 = 0.0 if not np.isfinite(m0) else m0
            dm = _detect_mask(cut - m0, valid, cfg.mask_detect, cfg.detect_grow, cfg.detect_minarea, cfg.protect_r, cfg.detect_smooth)
            valid &= ~dm
        # background
        if cfg.bkg == 'annulus':
            b, sg = _sigclip_median(cut[valid & ring])
            if not np.isfinite(b):
                reason[i] = 'no background pixels'
                continue
        else:
            b = 0.0
            _, sg = _sigclip_median(cut[valid & ring])
        cut = cut - b
        bk[i], rms[i] = b, sg
        vf[i] = valid.mean()
        if vf[i] < cfg.min_valid:
            reason[i] = 'masked fraction %.2f' % (1 - vf[i])
            continue
        if (~valid[R <= cfg.core_r]).any():
            reason[i] = 'core masked'
            continue
        nv = 1.0
        if cfg.norm == 'flux':
            a = (R <= cfg.norm_r) & valid
            f = float(cut[a].mean() * (R <= cfg.norm_r).sum()) if a.any() else float('nan')
            if not (np.isfinite(f) and f > 0):
                reason[i] = 'non-positive norm flux'
                continue
            nv = 1.0 / f
        elif cfg.norm == 'peak':
            pk = float(np.nanmax(np.where(valid & (R <= cfg.core_r + 1), cut, np.nan)))
            if not (np.isfinite(pk) and pk > 0):
                reason[i] = 'non-positive peak'
                continue
            nv = 1.0 / pk
        elif cfg.norm == 'column':
            nv = float(norm_values[i])
            if not np.isfinite(nv) or nv == 0:
                reason[i] = 'bad norm value'
                continue
        nrm[i] = nv
        if cfg.weights == 'invvar':
            wgt[i] = 1.0 / max((sg * abs(nv)) ** 2, 1e-30) if np.isfinite(sg) and sg > 0 else 0.0
        elif cfg.weights == 'column':
            wgt[i] = float(weight_values[i])
        cut = cut * nv
        cube[i] = np.where(valid, cut, np.nan).astype(np.float32)
        ok[i] = True
    return dict(cube=cube, ok=ok, reason=reason, bkg=bk, rms=rms, norm=nrm, weight=wgt, valid_frac=vf)


# --------------------------------------------------------------------------------------------------------------- combine
def combine(cube, method='mean', weights=None, clip=3.0):
    """Pixelwise combination of an (N,S,S) cube with NaN = masked.  Returns (S,S) float64 (NaN where no data)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if method == 'median':
            return np.nanmedian(cube, axis=0).astype(np.float64)
        if method == 'mean':
            return np.nanmean(cube, axis=0, dtype=np.float64)
        if method == 'wmean':
            w = np.ones(len(cube)) if weights is None else np.asarray(weights, float)
            ww = np.where(np.isfinite(cube), w[:, None, None], 0.0)
            num = np.nansum(cube * w[:, None, None], axis=0, dtype=np.float64)
            den = ww.sum(axis=0)
            return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        if method == 'clipmean':
            c = cube.astype(np.float64).copy()
            for _ in range(4):
                m = np.nanmean(c, axis=0)
                s = np.nanstd(c, axis=0)
                c = np.where(np.abs(c - m) > clip * np.maximum(s, 1e-30), np.nan, c)
            return np.nanmean(c, axis=0)
    raise ValueError('unknown method %s' % method)


def radial_bins(half, rlin=6.0, fac=1.25, rmax=None):
    """Bin edges: 1 px steps (first bin 0-0.5... centred on integer radii) up to rlin, then geometric steps."""
    rmax = float(half if rmax is None else rmax)
    e = [0.0, 0.5]
    r = 0.5
    while r < min(rlin, rmax):
        r += 1.0
        e.append(r)
    while e[-1] < rmax:
        e.append(e[-1] * fac if e[-1] * fac - e[-1] > 1.0 else e[-1] + 1.0)
    e = np.array(e)
    e[-1] = min(e[-1], rmax + 1e-6)
    return e


class Profiler:
    """Radial (elliptical) profile of a stack image about its centre, with precomputed bin indices."""

    def __init__(self, S, edges, q=1.0, pa=0.0):
        half = S // 2
        y, x = np.indices((S, S)) - half
        t = math.radians(pa)
        xr = x * math.cos(t) + y * math.sin(t)
        yr = -x * math.sin(t) + y * math.cos(t)
        self.r = np.hypot(xr, yr / max(q, 1e-3)) * 1.0
        # semi-major axis radius: the ellipse area scaling is kept (r is the major-axis radius)
        self.edges = np.asarray(edges, float)
        self.idx = np.digitize(self.r.ravel(), self.edges) - 1
        self.nb = len(self.edges) - 1
        self.q = q
        self.rmid = np.zeros(self.nb)
        for b in range(self.nb):
            m = self.idx == b
            self.rmid[b] = self.r.ravel()[m].mean() if m.any() else 0.5 * (self.edges[b] + self.edges[b + 1])

    def __call__(self, img, estimator='mean'):
        f = img.ravel()
        v = np.isfinite(f)
        ib = np.where(v, self.idx, -1)
        out = np.full(self.nb, np.nan)
        cnt = np.zeros(self.nb)
        if estimator == 'mean':
            ok = ib >= 0
            ok &= ib < self.nb
            s = np.bincount(ib[ok], weights=f[ok], minlength=self.nb)[:self.nb]
            cnt = np.bincount(ib[ok], minlength=self.nb)[:self.nb].astype(float)
            out = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
        else:
            for b in range(self.nb):
                m = ib == b
                cnt[b] = m.sum()
                if cnt[b]:
                    out[b] = np.median(f[m])
        return out, cnt


def aperture(img, r, centre=None):
    """Sum of the pixels in a circle (NaN pixels contribute the mean of the valid ones in the aperture)."""
    S = img.shape[0]
    c = S // 2 if centre is None else centre
    y, x = np.indices(img.shape)
    # exact-ish: sub-sample the boundary pixels 5x5
    d = np.hypot(x - c, y - c)
    inside = d <= r - 0.7
    edge = (d > r - 0.7) & (d < r + 0.7)
    f = img.copy()
    m = np.isfinite(f)
    if not m[(d <= r)].any():
        return float('nan')
    fill = float(np.mean(f[m & (d <= r)]))
    f = np.where(m, f, fill)
    tot = f[inside].sum()
    sub = (np.arange(5) - 2) / 5.0
    for yy, xx in zip(*np.where(edge)):
        dd = np.hypot(xx + sub[None, :] - c, yy + sub[:, None] - c)
        tot += f[yy, xx] * (dd <= r).mean()
    return float(tot)


def stack_cutouts(cube, ok, cfg, weights=None, profiler=None, edges=None, keep_boots=False):
    """Combine the accepted cutouts and bootstrap over the sources.

    Returns dict(stack, err_map, profile, profile_err, npix, aper, aper_err, n, boot_profiles)."""
    idx = np.where(ok)[0]
    n = len(idx)
    if n == 0:
        raise ValueError('no valid cutouts to stack')
    cub = cube[idx]
    w = None if weights is None else np.asarray(weights)[idx]
    S = cube.shape[1]
    if profiler is None:
        profiler = Profiler(S, edges if edges is not None else radial_bins(S // 2), cfg.q, cfg.pa)
    st = combine(cub, cfg.method, w, cfg.clip)
    prof, cnt = profiler(st, cfg.estimator)
    ap = aperture(st, cfg.ap_r)
    nb = int(cfg.n_boot)
    rng = np.random.default_rng(int(cfg.seed))
    bp = np.zeros((nb, len(prof)))
    ba = np.zeros(nb)
    s1 = np.zeros_like(st)
    s2 = np.zeros_like(st)
    cnt_ok = np.zeros_like(st)
    for b in range(nb):
        j = rng.integers(0, n, n)
        sb = combine(cub[j], cfg.method, None if w is None else w[j], cfg.clip)
        bp[b], _ = profiler(sb, cfg.estimator)
        ba[b] = aperture(sb, cfg.ap_r)
        g = np.isfinite(sb)
        s1 += np.where(g, sb, 0)
        s2 += np.where(g, sb * sb, 0)
        cnt_ok += g
    if nb > 1:
        mean = s1 / np.maximum(cnt_ok, 1)
        err_map = np.sqrt(np.maximum(s2 / np.maximum(cnt_ok, 1) - mean ** 2, 0)) * math.sqrt(nb / max(nb - 1, 1))
        perr = np.nanstd(bp, axis=0, ddof=1)
        aerr = float(np.nanstd(ba, ddof=1))
    else:
        err_map = np.full_like(st, np.nan)
        perr = np.full_like(prof, np.nan)
        aerr = float('nan')
    out = dict(stack=st, err_map=err_map, profile=prof, profile_err=perr, npix=cnt, aper=ap, aper_err=aerr, n=n,
               r=profiler.rmid, edges=profiler.edges, index=idx)
    if keep_boots:
        out['boot_profiles'] = bp
        out['boot_aper'] = ba
    return out


def run_stack(data, bad, xs, ys, cfg, null=0, avoid=None, **kw):
    """Full run: cutouts -> stack (+ optional null stack of `null` blank positions with the same pipeline).

    kw are passed to make_cutouts.  avoid = (M,2) catalog positions that blank positions must avoid (radius cfg.core_r+3)."""
    cut = make_cutouts(data, bad, xs, ys, cfg, **kw)
    res = stack_cutouts(cut['cube'], cut['ok'], cfg, cut['weight'])
    res['cutouts'] = cut
    if null:
        nx_, ny_ = null_positions(data, bad, int(null), cfg, avoid=avoid)
        nkw = {k: v for k, v in kw.items() if k in ('bkg_map', 'neigh')}
        nc = make_cutouts(data, bad, nx_, ny_, cfg, **nkw)
        if nc['ok'].sum() >= 3:
            res['null'] = stack_cutouts(nc['cube'], nc['ok'], cfg, nc['weight'])
            res['null_cutouts'] = nc
    return res


def null_positions(data, bad, n, cfg, avoid=None, seed=None, source_mask=None, tries=20):
    """Blank-sky random positions: not within the cutout core of a bad pixel / known source."""
    rng = np.random.default_rng(int(cfg.seed) + 777 if seed is None else seed)
    ny, nx = data.shape
    m = int(math.ceil(cfg.half * 1.0)) + 2
    if ny <= 2 * m or nx <= 2 * m:
        raise ValueError('image too small for the cutout size')
    xs, ys = [], []
    bm = np.zeros(data.shape, bool) if bad is None else bad
    bm = bm | ~np.isfinite(data)
    if source_mask is not None:
        bm = bm | source_mask
    tree = None
    if avoid is not None and len(avoid):
        from scipy.spatial import cKDTree
        tree = cKDTree(np.asarray(avoid, float)[:, :2])
    r0 = int(math.ceil(cfg.core_r + 2))
    for _ in range(n * tries):
        x = rng.uniform(m, nx - 1 - m)
        y = rng.uniform(m, ny - 1 - m)
        xi, yi = int(round(x)), int(round(y))
        if bm[yi - r0:yi + r0 + 1, xi - r0:xi + r0 + 1].any():
            continue
        if tree is not None and tree.query([x, y])[0] < cfg.core_r + 3:
            continue
        xs.append(x)
        ys.append(y)
        if len(xs) >= n:
            break
    return np.array(xs), np.array(ys)


# ------------------------------------------------------------------------------------------------------------- products
def surface_brightness(I, Ierr, zp, pixscale, nsig=3.0):
    """counts/pixel -> mag/arcsec^2.  Returns (mu, mu_err, mu_limit) with NaN mu where I <= nsig*err (mu_limit = limit then)."""
    I = np.asarray(I, float)
    Ierr = np.asarray(Ierr, float)
    k = zp + 5.0 * math.log10(pixscale)
    with np.errstate(divide='ignore', invalid='ignore'):
        mu = k - 2.5 * np.log10(np.where(I > 0, I, np.nan))
        mue = 1.0857 * Ierr / np.where(I > 0, I, np.nan)
        lim = k - 2.5 * np.log10(nsig * Ierr)
    sig = I > nsig * Ierr
    return np.where(sig, mu, np.nan), np.where(sig, mue, np.nan), lim


def summarize(res, cfg, zp=None, pixscale=None):
    s = dict(n_stacked=int(res['n']), method=cfg.method, n_boot=int(cfg.n_boot), aperture_radius=float(cfg.ap_r),
             aperture_sum=float(res['aper']), aperture_err=float(res['aper_err']))
    if np.isfinite(res['aper_err']) and res['aper_err'] > 0:
        s['aperture_snr'] = float(res['aper'] / res['aper_err'])
    if zp is not None and np.isfinite(res['aper']) and res['aper'] > 0 and cfg.norm == 'none':
        s['aperture_mag'] = float(zp - 2.5 * math.log10(res['aper']))
        if np.isfinite(res['aper_err']):
            s['aperture_mag_err'] = float(1.0857 * res['aper_err'] / res['aper'])
    if res.get('null') is not None:
        nl = res['null']
        s['null'] = dict(n=int(nl['n']), aperture_sum=float(nl['aper']), aperture_err=float(nl['aper_err']))
    return s
