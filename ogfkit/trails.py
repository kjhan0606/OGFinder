"""Linear trail (satellite / aircraft) detection, masking, filling and catalogue flagging.

Detection (detect_trails):
  1. background subtraction (sep mesh or a block-median fallback), robust noise, smoothing, clipping of bright pixels
     to +-clip sigma (a galaxy can contribute at most `clip` per pixel to a line sum);
  2. clipped Radon transform of the (binned) standardised image, normalised by sqrt(number of valid pixels on the line) ->
     candidate lines (peaks of the sinogram above `cand`);
  3. every candidate is verified at full resolution: the strip of width w along the line is sampled, the best contiguous
     along-track run (maximum-sum segment) is found for several strip widths, and its score is compared with the same
     statistic on parallel *null* lines of the same image (robust z-score) -> `zscore`.  Accepted when zscore >= `threshold`,
     run length >= `min_length`;
  4. the perpendicular profile (median along the run, so crossing galaxies do not matter) is fitted with a box convolved
     with a Gaussian -> centre, half-width a, blur s, amplitude; the mask half-width is a + `edge_sigma` * s + `margin`.
Coordinates: arrays are [row=y, col=x] 0-based internally; the public trail dicts use 1-based FITS pixel coordinates (x1, y1, x2, y2).
Angle theta: direction of the line from +x towards +y, in [0, 180) degrees.
"""
import math
import numpy as np

TRAIL_BIT = 32          # bit of the shared mask manager (ds9_mask.py)


# ------------------------------------------------------------------ preparation
def robust_sigma(a):
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    if a.size < 10:
        return float('nan')
    m = np.median(a)
    return 1.4826 * float(np.median(np.abs(a - m)))


def background(img, valid, mesh=64):
    """smooth background estimate: sep.Background when available, else a block-median mesh interpolated bilinearly"""
    from scipy.ndimage import zoom
    ny, nx = img.shape
    try:
        import sep
        d = np.ascontiguousarray(np.where(valid, img, 0.0), dtype=np.float32)
        bw = max(16, min(mesh, nx // 3, ny // 3))
        b = sep.Background(d, mask=~valid, bw=bw, bh=bw, fw=3, fh=3)
        return np.asarray(b.back(), dtype=np.float32)
    except Exception:
        pass
    m = max(16, min(mesh, nx // 3, ny // 3))
    gy, gx = max(2, ny // m), max(2, nx // m)
    out = np.zeros((gy, gx), np.float32)
    ys = np.linspace(0, ny, gy + 1).astype(int); xs = np.linspace(0, nx, gx + 1).astype(int)
    for j in range(gy):
        for i in range(gx):
            v = img[ys[j]:ys[j + 1], xs[i]:xs[i + 1]][valid[ys[j]:ys[j + 1], xs[i]:xs[i + 1]]]
            out[j, i] = np.median(v) if v.size else 0.0
    return zoom(out, (ny / gy, nx / gx), order=1)[:ny, :nx].astype(np.float32)


def local_noise(a, valid, block=32, hp_sigma=6.0):
    """noise map of the image `a` from its high-pass (a - gaussian(a, hp_sigma)): block-wise 1.4826 * MAD, bilinearly interpolated.
    Returns (map, global) where global is the 25th percentile of the block values (sky level, robust to galaxies filling the field)."""
    from scipy.ndimage import gaussian_filter, zoom
    ny, nx = a.shape
    w = valid.astype(np.float32)
    lp = gaussian_filter(np.where(valid, a, 0.0), hp_sigma) / np.maximum(gaussian_filter(w, hp_sigma), 1e-3)
    hp = np.where(valid, a - lp, np.nan)
    gy, gx = max(1, ny // block), max(1, nx // block)
    ys = np.linspace(0, ny, gy + 1).astype(int); xs = np.linspace(0, nx, gx + 1).astype(int)
    g = np.full((gy, gx), np.nan, np.float32)
    for j in range(gy):
        for i in range(gx):
            v = hp[ys[j]:ys[j + 1], xs[i]:xs[i + 1]]
            v = v[np.isfinite(v)]
            if v.size > 0.3 * block * block:
                g[j, i] = 1.4826 * np.median(np.abs(v - np.median(v)))
    good = np.isfinite(g) & (g > 0)
    if not good.any():
        sg = max(float(np.nanstd(hp)), 1e-30)
        return np.full(a.shape, sg, np.float32), sg
    glob = float(np.percentile(g[good], 25))
    g = np.where(good, g, glob)
    if gy * gx == 1:
        return np.full(a.shape, float(g[0, 0]), np.float32), glob
    m = zoom(g, (ny / gy, nx / gx), order=1)
    pad = np.pad(m, ((0, max(0, ny - m.shape[0])), (0, max(0, nx - m.shape[1]))), mode='edge')[:ny, :nx]
    return pad.astype(np.float32), glob


def standardise(img, valid, smooth=1.0, clip=3.0, mesh=32):
    """-> z (smoothed, background subtracted, divided by the local noise of the smoothed image, clipped to +-clip; 0 where invalid), sigma_pix, bkg.
    The local noise (high-pass MAD per 32-pixel block, floored at 0.7 x the sky value) down-weights galaxies, whose scatter is large."""
    from scipy.ndimage import gaussian_filter
    img = np.asarray(img, np.float32)
    bkg = background(img, valid, mesh)
    r = np.where(valid, img - bkg, 0.0).astype(np.float32)
    nmap0, sig_pix = local_noise(r, valid)
    if smooth and smooth > 0:
        w = gaussian_filter(valid.astype(np.float32), smooth)
        s = gaussian_filter(r, smooth) / np.maximum(w, 1e-3)
        s[w < 0.3] = 0.0
    else:
        s = r
    nmap, glob = local_noise(s, valid)
    nmap = np.maximum(nmap, 0.7 * glob)
    z = np.clip(s / np.maximum(nmap, 1e-30), -clip, clip)
    z[~valid] = 0.0
    return z.astype(np.float32), float(sig_pix), bkg


def bin_image(a, w, b):
    if b <= 1:
        return a, w
    ny, nx = a.shape
    ny2, nx2 = ny // b * b, nx // b * b
    A = a[:ny2, :nx2].reshape(ny2 // b, b, nx2 // b, b)
    W = w[:ny2, :nx2].reshape(ny2 // b, b, nx2 // b, b)
    ws = W.sum((1, 3))
    return (A * W).sum((1, 3)) / np.maximum(ws, 1e-6), (ws / (b * b)).astype(np.float32)


# ------------------------------------------------------------------ geometry
# the line (theta, rho): points  C + (rho + t) n + s d  with d = (cos th, sin th), n = (-sin th, cos th), C = image centre
def chord(theta_deg, rho, nx, ny):
    """along-track interval [s0, s1] of the line inside the image (None if it misses it)"""
    th = math.radians(theta_deg)
    c, sn = math.cos(th), math.sin(th)
    hx, hy = (nx - 1) / 2.0, (ny - 1) / 2.0
    lo, hi = -1e18, 1e18
    for p0, d, h in ((-rho * sn, c, hx), (rho * c, sn, hy)):
        if abs(d) < 1e-12:
            if abs(p0) > h:
                return None
            continue
        a, b = (-h - p0) / d, (h - p0) / d
        if a > b:
            a, b = b, a
        lo, hi = max(lo, a), min(hi, b)
    return (lo, hi) if hi > lo else None


def to_xy(theta_deg, rho, s, t, nx, ny):
    th = math.radians(theta_deg)
    return ((nx - 1) / 2.0 - (rho + t) * math.sin(th) + s * math.cos(th),
            (ny - 1) / 2.0 + (rho + t) * math.cos(th) + s * math.sin(th))


def sample_strip(img, theta, rho, s, t, order=1, cval=np.nan):
    """samples of img at along-track s (1-D) x perpendicular t (1-D): array (len(t), len(s))"""
    from scipy.ndimage import map_coordinates
    ny, nx = img.shape
    S, T = np.meshgrid(s, t)
    X, Y = to_xy(theta, rho, S, T, nx, ny)
    return map_coordinates(img, [Y, X], order=order, mode='constant', cval=cval)


# ------------------------------------------------------------------ stage 1: Radon candidates
def radon(z, w, step=0.5):
    """clipped, validity-weighted Radon transform with linear splatting -> angles, rho axis, sums, counts (n_angles x n_rho)"""
    ny, nx = z.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float32)
    xx -= (nx - 1) / 2.0; yy -= (ny - 1) / 2.0
    half = int(math.ceil(0.5 * math.hypot(nx, ny))) + 2
    nb = 2 * half + 1
    angles = np.arange(0.0, 180.0, step)
    zs = np.zeros((len(angles), nb), np.float64); ws = np.zeros_like(zs)
    sel = w.ravel() > 0
    zf, wf = z.ravel()[sel].astype(np.float64) * w.ravel()[sel], w.ravel()[sel].astype(np.float64)
    xf, yf = xx.ravel()[sel], yy.ravel()[sel]
    for i, a in enumerate(angles):
        th = math.radians(a)
        r = -xf * math.sin(th) + yf * math.cos(th) + half
        i0 = np.floor(r).astype(np.int64); f = r - i0
        for k, wt in ((0, 1 - f), (1, f)):
            zs[i] += np.bincount(i0 + k, weights=zf * wt, minlength=nb + 1)[:nb]
            ws[i] += np.bincount(i0 + k, weights=wf * wt, minlength=nb + 1)[:nb]
    return angles, np.arange(-half, half + 1), zs, ws


def candidates(z, w, binning, nmax=40, cand=4.0, min_count=30, step=0.5):
    """peaks of the sinogram S/N = sum / sqrt(count), robustly normalised per angle -> list of dict(theta, rho (original pixels), snr)"""
    from scipy.ndimage import maximum_filter, uniform_filter1d
    angles, rho, zs, ws = radon(z, w, step)
    zs = uniform_filter1d(zs, 3, axis=1); ws = uniform_filter1d(ws, 3, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        sn = zs / np.sqrt(np.maximum(ws, 1e-9))
    bad = ws < min_count
    sn[bad] = np.nan
    with np.errstate(all='ignore'):
        med = np.nanmedian(sn, axis=1, keepdims=True)
        mad = 1.4826 * np.nanmedian(np.abs(sn - med), axis=1, keepdims=True)
    sn = (sn - med) / np.maximum(mad, 1e-3)
    sn[bad] = 0.0
    sn = np.nan_to_num(sn)
    mx = maximum_filter(sn, size=(9, 15), mode='nearest')
    pk = np.argwhere((sn == mx) & (sn > cand))
    out = [dict(theta=float(angles[i]), rho=float(rho[j]) * binning, snr=float(sn[i, j])) for i, j in pk]
    out.sort(key=lambda d: -d['snr'])
    return out[:nmax], sn


# ------------------------------------------------------------------ stage 2: verification by an along-track scan statistic
WIDTHS = (3, 6, 12, 24)
KAPPAS = (0.4, 0.8, 1.6)


def best_runs(U, min_len, kappa):
    """U (rows, n): per row the maximum-sum segment of (U - kappa) with length >= min_len -> score (sum U / sqrt(n)), i0, i1 arrays"""
    nrow, n = U.shape
    L = int(min_len)
    if n < L + 1:
        return np.zeros(nrow), np.zeros(nrow, int), np.zeros(nrow, int)
    V = U - kappa
    C = np.cumsum(V, axis=1)
    Cp = np.concatenate([np.zeros((nrow, 1)), C[:, :-1]], axis=1)
    mn = np.minimum.accumulate(Cp, axis=1)
    gain = C[:, L - 1:] - mn[:, :n - L + 1]
    j = np.argmax(gain, axis=1) + L - 1
    rows = np.arange(nrow)
    i0 = np.array([int(np.argmin(Cp[r, :j[r] - L + 2])) for r in rows])
    cs = np.concatenate([np.zeros((nrow, 1)), np.cumsum(U, axis=1)], axis=1)
    su = cs[rows, j + 1] - cs[rows, i0]
    return su / np.sqrt(j - i0 + 1), i0, j


class Scanner:
    """verification statistic on the original-resolution standardised image z (weights w)"""

    def __init__(self, z, w, min_len, widths=WIDTHS, kappas=KAPPAS, dtheta=0.075, ntheta=9, drho=6):
        self.z, self.w = z, w
        self.min_len, self.widths, self.kappas = int(min_len), widths, kappas
        self.dth = (np.arange(ntheta) - (ntheta - 1) / 2.0) * dtheta
        self.drho = drho
        self.sig = {wd: 1.0 for wd in widths}
        self.allvalid = bool((w > 0.5).all())

    def strips(self, theta, rho):
        """-> (s, {width: (t0 offsets, U (n_off, n_s))}) or None.  U is the *contrast* of the strip of width w centred on t0: its mean minus the mean of
        the two flanking strips of the same width (centred at t0 +- (w + 2)), so smooth galaxy light and gradients cancel and only a ridge remains."""
        ny, nx = self.z.shape
        ch = chord(theta, rho, nx, ny)
        if ch is None or ch[1] - ch[0] < self.min_len + 4:
            return None
        s = np.arange(math.ceil(ch[0]) + 1, math.floor(ch[1]) - 1, 1.0)
        if len(s) < self.min_len + 2:
            return None
        D = self.drho
        wmax = max(self.widths)
        ext = D + 1.5 * wmax + 4
        t = np.arange(-ext, ext + 1, 1.0)
        Z = sample_strip(self.z, theta, rho, s, t)
        ok = np.isfinite(Z)
        if not self.allvalid:
            ok &= sample_strip(self.w, theta, rho, s, t) > 0.5
        Z = np.where(ok, Z, 0.0)
        cz = np.concatenate([np.zeros((1, len(s))), np.cumsum(Z, axis=0)])
        cw = np.concatenate([np.zeros((1, len(s))), np.cumsum(ok.astype(float), axis=0)])
        t0s = np.arange(-D, D + 1.0, 1.0)

        def mean_at(c0, wd):
            lo = np.round(c0 - wd / 2.0 - t[0]).astype(int); hi = np.round(c0 + wd / 2.0 - t[0]).astype(int) + 1
            sz = cz[hi] - cz[lo]; sw = cw[hi] - cw[lo]
            good = sw >= 0.6 * (hi - lo)[:, None]
            return np.where(good, sz / np.maximum(sw, 1.0), np.nan)
        out = {}
        for wd in self.widths:
            m0 = mean_at(t0s, wd)
            ml, mr = mean_at(t0s - (wd + 2), wd), mean_at(t0s + (wd + 2), wd)
            fl = np.where(np.isfinite(ml) & np.isfinite(mr), 0.5 * (ml + mr), np.where(np.isfinite(ml), ml, mr))
            U = np.where(np.isfinite(m0) & np.isfinite(fl), m0 - fl, 0.0)
            out[wd] = (t0s, U)
        return s, out

    def centre_samples(self, theta, rho):
        r = self.strips(theta, rho)
        if r is None:
            return None
        s, d = r
        return {wd: d[wd][1][d[wd][0].size // 2] for wd in d}

    def best(self, theta, rho):
        """maximum over the neighbourhood (theta +- , rho +- drho): dict (width, kappa) -> (score, theta', rho', s0, s1)"""
        res = {}
        for dth in self.dth:
            r = self.strips(theta + dth, rho)
            if r is None:
                continue
            s, d = r
            for wd in self.widths:
                t0s, U = d[wd]
                Un = U / self.sig[wd]
                for kp in self.kappas:
                    sc, i0, i1 = best_runs(Un, self.min_len, kp)
                    j = int(np.argmax(sc))
                    key = (wd, kp)
                    if key not in res or sc[j] > res[key][0]:
                        res[key] = (float(sc[j]), theta + float(dth), rho + float(t0s[j]), float(s[i0[j]]), float(s[i1[j]]))
        return res


def random_lines(nx, ny, n, min_chord, rng):
    out, tries = [], 0
    R = 0.5 * math.hypot(nx, ny)
    while len(out) < n and tries < 20000:
        tries += 1
        th = rng.uniform(0, 180); rho = rng.uniform(-R, R)
        ch = chord(th, rho, nx, ny)
        if ch is not None and ch[1] - ch[0] >= min_chord:
            out.append((th, rho))
    return out


# ------------------------------------------------------------------ profile fit
def box_blur_profile(t, A, t0, a, s, c):
    from scipy.special import erf
    s = max(s, 0.3)
    return c + A * 0.5 * (erf((t - t0 + a) / (math.sqrt(2) * s)) - erf((t - t0 - a) / (math.sqrt(2) * s)))


def measure_profile(img, valid, theta, rho, s0, s1, sigma_pix, wmax=60.0, a_guess=3.0):
    """perpendicular profile (median along the run, robust against crossing galaxies) fitted with a box (x) Gaussian"""
    from scipy.optimize import least_squares
    s = np.arange(s0, s1 + 1.0, 1.0)
    t = np.arange(-wmax, wmax + 1.0, 1.0)
    Z = sample_strip(img, theta, rho, s, t)
    if not valid.all():
        Z = np.where(sample_strip(valid.astype(np.float32), theta, rho, s, t) > 0.5, Z, np.nan)
    n = np.isfinite(Z).sum(1)
    with np.errstate(all='ignore'):
        P = np.nanmedian(Z, axis=1)
    err = 1.2533 * sigma_pix / np.sqrt(np.maximum(n, 1))
    ok = np.isfinite(P) & (n > 5)
    if ok.sum() < 20:
        return None
    far = ok & (np.abs(t) > wmax * 0.6)
    off = float(np.median(P[far])) if far.sum() > 5 else float(np.median(P[ok]))
    sm = np.convolve(np.nan_to_num(P - off), np.ones(5) / 5, 'same')
    k = int(np.argmax(np.where(ok & (np.abs(t) < wmax * 0.5), sm, -1e30)))
    A0 = max(float(sm[k]), 0.2 * sigma_pix)
    best = None
    for a0 in (a_guess, 1.0, 6.0, 12.0):
        p0 = [A0, float(t[k]), min(max(a0, 0.4), wmax / 2 - 1), 1.2, off]
        try:
            r = least_squares(lambda p: (box_blur_profile(t[ok], *p) - P[ok]) / err[ok], p0, x_scale='jac',
                              bounds=([0, -wmax / 2, 0.3, 0.3, off - 5 * sigma_pix], [np.inf, wmax / 2, wmax / 2, 8.0, off + 5 * sigma_pix]))
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        return None
    A, t0, a, sb, c = best.x
    tf = np.arange(-wmax, wmax, 0.05)
    mf = box_blur_profile(tf, A, t0, a, sb, 0.0)
    above = tf[mf >= 0.5 * mf.max()]
    fwhm = float(above.max() - above.min()) if above.size else 2 * a
    return dict(amp=float(A), t0=float(t0), half_width=float(a), blur=float(sb), offset=float(c), fwhm=fwhm,
                chi2_red=float(2 * best.cost / max(1, ok.sum() - 5)), profile_snr=float(A / np.median(err[ok])),
                t=t, profile=P - c, model=box_blur_profile(t, A, t0, a, sb, 0.0))


# ------------------------------------------------------------------ main detector
def bleed_like(img, theta, rho, s0, s1, amp, level, axis_tol=3.0, fac=300.0, spike_ratio=3.0):
    """True when the run is a feature of a very bright compact source rather than a trail: (a) a run within axis_tol degrees of the pixel axes
    through a core brighter than fac x the run amplitude (CCD bleed column, ACS spikes ~2 deg off the axes), or (b) at any angle, a run through
    such a core whose brightness falls off away from the core (median in 8 px..15 % of the run from the core > spike_ratio x the median beyond
    30 %): the long diffraction spike of a bright star.  A trail has a constant amplitude along its length."""
    th = theta % 180.0
    L0 = s1 - s0
    sv = np.arange(s0 - 1.5 * L0, s1 + 1.5 * L0 + 1.0, 1.0)         # the core may lie beyond the part of the spike that was found
    import warnings
    with np.errstate(all='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        v = np.nanmax(sample_strip(img, theta, rho, sv, np.arange(-3.0, 4.0)), axis=0) - level
    if not np.isfinite(v).any() or np.nanmax(v) <= fac * max(amp, 1e-9):
        return False
    if min(abs(th - 90.0), th, 180.0 - th) <= axis_tol:
        return True
    d = np.abs(sv - sv[int(np.nanargmax(v))])
    ok = np.isfinite(v)
    inner, outer = (d > 8) & (d <= 0.15 * L0) & ok, (d > 0.3 * L0) & ok
    if inner.sum() < 5 or outer.sum() < 5:
        return False
    vi, vo = np.nanmedian(v[inner]), np.nanmedian(v[outer])
    return bool(vi > spike_ratio * max(vo, 0.2 * amp))


def _binning(shape, maxdim):
    return max(1, int(math.ceil(max(shape) / float(maxdim))))


def detect_trails(img, valid=None, threshold=8.0, min_length=None, max_trails=8, smooth=1.0, clip=3.0, cand=4.0, maxdim=768,
                  margin=2.0, edge_sigma=2.0, catalog=None, catalog_rescue=None, seed=1, n_null=100, mesh=32, min_aspect=12.0, max_fwhm=40.0, reject_bleeds=True):
    """-> dict(trails=[...], sigma_pix, binning, n_candidates, ...).  Each trail: id, x1,y1,x2,y2 (1-based FITS pixels), theta_deg, rho, length,
    halfwidth (mask half-width = box half-width + edge_sigma * blur + margin), box_halfwidth, blur, amp, amp_snr, zscore, run_score, s0, s1."""
    img = np.asarray(img, np.float32)
    ny, nx = img.shape
    valid = (np.isfinite(img) & (img != 0)) if valid is None else (valid & np.isfinite(img))
    diag = math.hypot(nx, ny)
    if min_length is None:
        min_length = max(30.0, 0.12 * diag)
    z, sig_pix, bkg = standardise(img, valid, smooth, clip, mesh)
    w = valid.astype(np.float32)
    b = _binning(img.shape, maxdim)
    zb, wb = bin_image(z, w, b)
    if b > 1:
        zb = np.clip(zb / max(robust_sigma(zb[wb > 0.9]), 1e-6), -clip, clip)
        zb[wb < 0.5] = 0.0
    cands, _ = candidates(zb, (wb > 0.5).astype(np.float32), b, nmax=60, cand=cand, min_count=max(20, min_length / b * 0.6))
    out = dict(trails=[], sigma_pix=float(sig_pix), binning=b, n_candidates=len(cands), min_length=float(min_length), threshold=float(threshold), _profiles=[])
    if not cands:
        return out
    sc = Scanner(z, w, int(min_length))
    rng = np.random.default_rng(seed)
    nulls = random_lines(nx, ny, n_null, 1.3 * min_length, rng)
    pool = {wd: [] for wd in sc.widths}
    for th, rho in nulls:
        r = sc.centre_samples(th, rho)
        if r is None:
            continue
        for wd, u in r.items():
            pool[wd].append(u[u != 0])
    for wd in sc.widths:
        v = np.concatenate(pool[wd]) if pool[wd] else np.array([1.0])
        sc.sig[wd] = max(robust_sigma(v), 1e-3)
    ns = {}
    for th, rho in nulls:
        for key, val in sc.best(th, rho).items():
            ns.setdefault(key, []).append(val[0])
    nstat = {}
    for k, v in ns.items():
        v = np.array(v)
        nstat[k] = (float(np.median(v)), max(1.4826 * float(np.median(np.abs(v - np.median(v)))), 0.05 * max(float(np.median(v)), 1.0)))
    out['null'] = {'%g/%g' % k: v for k, v in nstat.items()}
    def evaluate(c):
        res = sc.best(c['theta'], c['rho'])
        bestk, bestz = None, -1e9
        for key, val in res.items():
            if key in nstat:
                zsc = (val[0] - nstat[key][0]) / nstat[key][1]
                if zsc > bestz:
                    bestz, bestk = zsc, key
        if bestk is None:
            return None
        sc_, th, rho, s0, s1 = res[bestk]
        return dict(theta=th, rho=rho, s0=s0, s1=s1, zscore=float(bestz), run_score=float(sc_), width_strip=bestk[0], kappa=bestk[1], cand_snr=c['snr'])
    cl, cur = {}, {}
    for i, c in enumerate(cands):
        e = evaluate(c)
        if e is not None:
            cl[i], cur[i] = c, e
    out['best_zscore'] = max([f['zscore'] for f in cur.values()]) if cur else float('nan')
    # greedy acceptance with deflation: after a trail is accepted its band is blanked and the remaining candidates are re-verified, so that lines
    # that merely cross an accepted trail (or run along its wings) do not count as further trails
    accepted = []
    level = float(np.median(img[valid][::7])) if valid.any() else 0.0
    def deflate(f, sa, sb):
        """blank the band of a run in the scan arrays so that pieces/wings of the same line are not found again"""
        p1 = to_xy(f['theta'], f['rho'], sa, 0, nx, ny); p2 = to_xy(f['theta'], f['rho'], sb, 0, nx, ny)
        band = trail_mask(z.shape, [dict(x1=p1[0] + 1, y1=p1[1] + 1, x2=p2[0] + 1, y2=p2[1] + 1, halfwidth=f['width_strip'] / 2.0 + 6.0)])
        sc.z = np.where(band, 0.0, sc.z).astype(np.float32); sc.w = np.where(band, 0.0, sc.w).astype(np.float32); sc.allvalid = False

    def reverify():
        for k in list(cur):
            if cur[k]['zscore'] >= 0.5 * threshold:
                e = evaluate(cl[k])
                if e is None:
                    cur.pop(k)
                else:
                    cur[k] = e
            else:
                cur.pop(k)

    while len(accepted) < max_trails and cur:
        i = max(cur, key=lambda k: cur[k]['zscore'])
        f = cur.pop(i)
        if f['zscore'] < threshold:
            # catalogue check: a weaker candidate is accepted when >= 2 elongated, aligned catalogue objects lie on it
            if catalog is None or not catalog_rescue or f['zscore'] < catalog_rescue * threshold:
                break
            x1_, y1_ = to_xy(f['theta'], f['rho'], f['s0'], 0, nx, ny); x2_, y2_ = to_xy(f['theta'], f['rho'], f['s1'], 0, nx, ny)
            prov = dict(x1=x1_ + 1, y1=y1_ + 1, x2=x2_ + 1, y2=y2_ + 1, theta_deg=f['theta'] % 180.0, halfwidth=f['width_strip'] / 2.0 + 3.0)
            if catalog_support(catalog, prov) < 2:
                continue
            f['rescued'] = True
        if (f['s1'] - f['s0']) < min_length * 0.9:
            continue
        # shape check: a trail is long and thin; a galaxy (or a chain of them) gives a broad, short 'ridge'
        prof = measure_profile(img, valid, f['theta'], f['rho'], f['s0'], f['s1'], sig_pix, a_guess=max(1.0, f['width_strip'] / 2.0))
        fw = prof['fwhm'] if prof is not None else float(f['width_strip'])
        if fw > max_fwhm or (f['s1'] - f['s0']) < min_aspect * fw:
            out.setdefault('rejected', []).append(dict(theta=f['theta'], rho=f['rho'], length=f['s1'] - f['s0'], fwhm=fw, zscore=f['zscore'], reason='shape'))
            continue
        if reject_bleeds and prof is not None and bleed_like(img, f['theta'], f['rho'], f['s0'], f['s1'], prof['amp'], level):
            out.setdefault('rejected', []).append(dict(theta=f['theta'], rho=f['rho'], length=f['s1'] - f['s0'], fwhm=fw, zscore=f['zscore'], reason='bleed'))
            deflate(f, -diag, diag)
            reverify()
            continue
        f['prof'] = prof
        accepted.append(f)
        deflate(f, f['s0'], f['s1'])
        reverify()
    for i, f in enumerate(accepted, 1):
        th, rho, s0, s1 = f['theta'], f['rho'], f['s0'], f['s1']
        prof = f.get('prof')
        if prof is not None and abs(prof['t0']) < 30:
            rho += prof['t0']
        else:
            prof = None
        a = prof['half_width'] if prof else f['width_strip'] / 2.0
        sb_ = prof['blur'] if prof else 1.0
        hw = a + edge_sigma * sb_ + margin
        x1, y1 = to_xy(th, rho, s0, 0, nx, ny); x2, y2 = to_xy(th, rho, s1, 0, nx, ny)
        out['trails'].append(dict(rescued=bool(f.get('rescued', False)), id=i, x1=float(x1 + 1), y1=float(y1 + 1), x2=float(x2 + 1), y2=float(y2 + 1), theta_deg=float(th % 180.0), rho=float(rho),
                                  length=float(s1 - s0), halfwidth=float(hw), box_halfwidth=float(a), blur=float(sb_),
                                  fwhm=float(prof['fwhm']) if prof else float(f['width_strip']), amp=float(prof['amp']) if prof else float('nan'),
                                  amp_snr=float(prof['profile_snr']) if prof else float('nan'), zscore=f['zscore'], run_score=f['run_score'],
                                  width_strip=f['width_strip'], s0=float(s0), s1=float(s1)))
        out['_profiles'].append(prof)
    if catalog is not None:
        for t in out['trails']:
            t['catalog_support'] = catalog_support(catalog, t)
    return out


# ------------------------------------------------------------------ masks, fill, flags
def seg_distance(X, Y, t):
    """distance of pixel centres (0-based X, Y) to the trail segment (t['x1'..'y2'] are 1-based)"""
    ax, ay, bx, by = t['x1'] - 1, t['y1'] - 1, t['x2'] - 1, t['y2'] - 1
    dx, dy = bx - ax, by - ay
    L2 = max(dx * dx + dy * dy, 1e-12)
    u = np.clip(((X - ax) * dx + (Y - ay) * dy) / L2, 0.0, 1.0)
    return np.hypot(X - (ax + u * dx), Y - (ay + u * dy))


def trail_mask(shape, trails, halfwidth_scale=1.0, extra_margin=0.0, end_extend=0.0):
    """boolean mask (ny, nx): distance to the segment <= halfwidth (segments extended by end_extend px at both ends)"""
    ny, nx = shape
    m = np.zeros(shape, bool)
    for t in trails:
        t = dict(t)
        if end_extend:
            ux, uy = t['x2'] - t['x1'], t['y2'] - t['y1']
            L = math.hypot(ux, uy) or 1.0
            t['x1'] -= ux / L * end_extend; t['y1'] -= uy / L * end_extend; t['x2'] += ux / L * end_extend; t['y2'] += uy / L * end_extend
        hw = t['halfwidth'] * halfwidth_scale + extra_margin
        x0 = int(max(0, math.floor(min(t['x1'], t['x2']) - 1 - hw - 1))); x1 = int(min(nx, math.ceil(max(t['x1'], t['x2']) + hw + 1)))
        y0 = int(max(0, math.floor(min(t['y1'], t['y2']) - 1 - hw - 1))); y1 = int(min(ny, math.ceil(max(t['y1'], t['y2']) + hw + 1)))
        if x1 <= x0 or y1 <= y0:
            continue
        Y, X = np.mgrid[y0:y1, x0:x1]
        m[y0:y1, x0:x1] |= seg_distance(X, Y, t) <= hw
    return m


def fill_interpolate(img, trails, mask=None, noise=False, side=3, seed=0, sigma=None):
    """replace the trail band by a linear interpolation across the trail between its two sides (medians of `side`-pixel wide strips just outside
    the band, median-filtered along the track); optional Gaussian noise of the sky sigma so that the noise level is preserved."""
    from scipy.ndimage import distance_transform_edt, median_filter
    img = np.asarray(img, np.float32)
    out = img.copy()
    ny, nx = img.shape
    full = np.zeros(img.shape, bool) if mask is None else mask.copy()
    for t in trails:
        th, rho, s0, s1, hw = t['theta_deg'], t['rho'], t['s0'], t['s1'], t['halfwidth']
        s = np.arange(s0 - hw - 2, s1 + hw + 3, 1.0)
        ts = np.arange(-hw - side - 1, hw + side + 2, 1.0)
        Z = sample_strip(img, th, rho, s, ts)
        with np.errstate(all='ignore'):
            R = np.nanmedian(Z[(ts >= hw + 0.5) & (ts <= hw + 0.5 + side)], axis=0)
            Lf = np.nanmedian(Z[(ts <= -hw - 0.5) & (ts >= -hw - 0.5 - side)], axis=0)

        def sm(v):
            good = np.isfinite(v)
            if good.sum() < 2:
                return np.full_like(v, np.nan)
            v = v.copy()
            v[~good] = np.interp(np.flatnonzero(~good), np.flatnonzero(good), v[good])
            return median_filter(v, size=9, mode='nearest')
        R, Lf = sm(R), sm(Lf)
        x0 = int(max(0, math.floor(min(t['x1'], t['x2']) - hw - 3))); x1 = int(min(nx, math.ceil(max(t['x1'], t['x2']) + hw + 3)))
        y0 = int(max(0, math.floor(min(t['y1'], t['y2']) - hw - 3))); y1 = int(min(ny, math.ceil(max(t['y1'], t['y2']) + hw + 3)))
        Y, X = np.mgrid[y0:y1, x0:x1]
        band = seg_distance(X, Y, t) <= hw
        if not band.any():
            continue
        thr = math.radians(th)
        xr, yr = X - (nx - 1) / 2.0, Y - (ny - 1) / 2.0
        sp = xr * math.cos(thr) + yr * math.sin(thr)
        tp = -xr * math.sin(thr) + yr * math.cos(thr) - rho
        si = np.clip(sp - s[0], 0, len(s) - 1)
        Rv = np.interp(si, np.arange(len(s)), R); Lv = np.interp(si, np.arange(len(s)), Lf)
        val = Lv + (Rv - Lv) * np.clip((tp + hw) / (2 * hw), 0, 1)
        sel = band & np.isfinite(val)
        sub = out[y0:y1, x0:x1]
        sub[sel] = val[sel]
        full[y0:y1, x0:x1] |= band
    if noise and full.any():
        sg = sigma if sigma is not None else robust_sigma(img[~full & np.isfinite(img) & (img != 0)])
        out[full] += np.random.default_rng(seed).normal(0.0, sg, size=int(full.sum())).astype(np.float32)
    bad = full & ~np.isfinite(out)
    if bad.any():
        idx = distance_transform_edt(~(np.isfinite(out) & ~full), return_distances=False, return_indices=True)
        out[bad] = out[tuple(i[bad] for i in idx)]
    return out


def catalog_support(cat, t, min_elong=3.0, angle_tol=12.0, width_extra=3.0):
    """number of elongated catalogue objects (A/B >= min_elong) lying on the trail and aligned with it (within angle_tol degrees)"""
    if not cat or 'A_IMAGE' not in cat or 'B_IMAGE' not in cat:
        return 0
    x, y = np.asarray(cat['X_IMAGE'], float) - 1, np.asarray(cat['Y_IMAGE'], float) - 1
    A, B = np.asarray(cat['A_IMAGE'], float), np.asarray(cat['B_IMAGE'], float)
    th = np.asarray(cat['THETA_IMAGE'], float) if 'THETA_IMAGE' in cat else np.zeros_like(A)
    sel = (A / np.maximum(B, 1e-3) >= min_elong) & (seg_distance(x, y, t) <= t['halfwidth'] + width_extra) & \
          (np.abs((th - t['theta_deg'] + 90) % 180 - 90) <= angle_tol)
    return int(sel.sum())


def flag_catalog(cat, trails, touch_k=2.5):
    """-> (TRAIL_FLAG, TRAIL_ID, TRAIL_DIST) arrays.  TRAIL_FLAG bits: 1 = footprint (ellipse touch_k * A x B along the trail normal) reaches the
    mask band, 2 = centre inside the band, 4 = elongated and aligned (looks like a fragment of the trail itself).  TRAIL_DIST = distance of the
    centre to the band edge in pixels (negative inside the band), -1 without trails."""
    n = len(cat['X_IMAGE'])
    x, y = np.asarray(cat['X_IMAGE'], float) - 1, np.asarray(cat['Y_IMAGE'], float) - 1
    A = np.asarray(cat['A_IMAGE'], float) if 'A_IMAGE' in cat else np.full(n, 1.5)
    B = np.asarray(cat['B_IMAGE'], float) if 'B_IMAGE' in cat else A
    th = np.asarray(cat['THETA_IMAGE'], float) if 'THETA_IMAGE' in cat else np.zeros(n)
    flag = np.zeros(n, int); tid = np.zeros(n, int); dist = np.full(n, 1e9)
    for t in trails:
        edge = seg_distance(x, y, t) - t['halfwidth']
        phi = np.radians(th) - math.radians(t['theta_deg'] + 90.0)
        ext = touch_k * np.sqrt((A * np.cos(phi)) ** 2 + (B * np.sin(phi)) ** 2)
        inside = edge <= 0
        elong = (A / np.maximum(B, 1e-3) >= 3.0) & (np.abs((th - t['theta_deg'] + 90) % 180 - 90) <= 12.0) & inside
        f = (edge <= ext) * 1 + inside * 2 + elong * 4
        upd = (f > 0) & (edge < dist)
        flag = np.where(upd, f, flag); tid = np.where(upd, t['id'], tid)
        dist = np.minimum(dist, edge)
    return flag, tid, np.where(dist > 1e8, -1.0, dist)


# ------------------------------------------------------------------ multi-frame stacking with trail exclusion
def stack_frames(frames, masks=None, method='sigclip', nsig=3.0, iters=5):
    """combine registered frames with per-frame boolean masks (True = exclude): 'median', 'mean' or 'sigclip' (iterative sigma-clipped mean).
    -> combined image (NaN where no valid input), number of frames used per pixel"""
    F = np.array([np.asarray(f, np.float64) for f in frames])
    M = np.zeros(F.shape, bool) if masks is None else np.array([np.zeros(F.shape[1:], bool) if m is None else m for m in masks])
    F = np.where(M | ~np.isfinite(F), np.nan, F)
    with np.errstate(all='ignore'):
        if method == 'median':
            out = np.nanmedian(F, axis=0)
        elif method == 'mean':
            out = np.nanmean(F, axis=0)
        else:
            for _ in range(iters):
                mu = np.nanmedian(F, axis=0)
                sd = 1.4826 * np.nanmedian(np.abs(F - mu), axis=0)
                sd = np.where(sd > 0, sd, np.nanstd(F, axis=0))
                F = np.where(np.abs(F - mu) > nsig * np.maximum(sd, 1e-30), np.nan, F)
            out = np.nanmean(F, axis=0)
    return out, np.isfinite(F).sum(0)
