"""Extensions to ogfkit.trails (the straight-line detector is not changed here):
  along_profile    amplitude along the track -> duty cycle of flickering trails
  detect_curved    gap-tolerant detection of curved / segmented trails by tile-wise detection + chain linking (segments become ordinary trail dicts
                   that share a `group` id, so masks, fills, catalogue flags and the undo all work unchanged)
  trail_flux       per-object contaminating flux: the fitted trail model integrated over each catalogue object's elliptical aperture
  wcs_resample     resample a frame (image or mask) onto another frame's pixel grid through the two WCS (stacking of dithered frames)
"""
import math
import warnings
import numpy as np
from scipy.special import erf
from . import trails as T


# ------------------------------------------------------------------ flicker
def along_profile(img, trail, nbin=12, bkg=None):
    """mean amplitude (box-averaged across the trail, flank-subtracted) in bins of `nbin` px along the track.
    -> dict(s, amp, duty, n_gaps, cv): duty = fraction of bins above half of the median bin amplitude of the brightest half; flicker when duty < 0.8."""
    th, rho, s0, s1 = trail['theta_deg'], trail['rho'], trail['s0'], trail['s1']
    a = max(trail.get('box_halfwidth', 2.0), 0.7)
    sv = np.arange(s0, s1 + 1.0, 1.0)
    tc = np.arange(-a, a + 0.5, 1.0)
    fl = np.concatenate([np.arange(-a - 12, -a - 5, 1.0), np.arange(a + 6, a + 13, 1.0)])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Zc = T.sample_strip(img, th, rho, sv, tc)
        Zf = T.sample_strip(img, th, rho, sv, fl)
        c = np.nanmean(Zc, axis=0) - np.nanmedian(Zf, axis=0)
    nb = max(1, len(sv) // nbin)
    cb = np.array([np.nanmedian(c[i * nbin:(i + 1) * nbin]) for i in range(nb)])
    ok = np.isfinite(cb)
    if ok.sum() < 3:
        return dict(s=sv[::nbin], amp=cb, duty=1.0, n_gaps=0, cv=0.0)
    ref = float(np.median(np.sort(cb[ok])[len(cb[ok]) // 2:]))
    on = cb > 0.5 * ref
    gaps = int(np.sum(np.diff(on[ok].astype(int)) == -1))
    return dict(s=sv[:nb * nbin:nbin], amp=cb, duty=float(np.mean(on[ok])), n_gaps=gaps, cv=float(np.std(cb[ok]) / max(abs(ref), 1e-9)))


# ------------------------------------------------------------------ curved / segmented trails
def _seg_from_tile(t, x0, y0, nx, ny):
    d = dict(t)
    for k in ('x1', 'x2'):
        d[k] = t[k] + x0
    for k in ('y1', 'y2'):
        d[k] = t[k] + y0
    return d


def _line_params(t, nx, ny):
    """theta (deg, 0..180), rho, s0, s1 of a segment from its end points (0-based centre convention of ogfkit.trails)"""
    ax, ay, bx, by = t['x1'] - 1, t['y1'] - 1, t['x2'] - 1, t['y2'] - 1
    th = math.degrees(math.atan2(by - ay, bx - ax))
    if th < 0 or th >= 180:
        ax, ay, bx, by = bx, by, ax, ay
        th = math.degrees(math.atan2(by - ay, bx - ax))
    thr = math.radians(th)
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    rho = -(ax - cx) * math.sin(thr) + (ay - cy) * math.cos(thr)
    s0 = (ax - cx) * math.cos(thr) + (ay - cy) * math.sin(thr)
    s1 = (bx - cx) * math.cos(thr) + (by - cy) * math.sin(thr)
    return th, rho, s0, s1


def _merge_dups(segs, nx, ny):
    segs = sorted(segs, key=lambda t: -t['zscore'])
    out = []
    for t in segs:
        th, rho, s0, s1 = _line_params(t, nx, ny)
        dup = False
        for u in out:
            if abs((u['_th'] - th + 90) % 180 - 90) < 3.0 and abs(u['_rho'] - rho) < max(4.0, u['halfwidth']) and \
                    min(u['_s1'], s1) - max(u['_s0'], s0) > 0.3 * min(u['_s1'] - u['_s0'], s1 - s0):
                dup = True
                break
        if not dup:
            t = dict(t, _th=th, _rho=rho, _s0=s0, _s1=s1)
            out.append(t)
    return out


def detect_curved(img, straight=None, valid=None, tile=None, tile_threshold=6.0, chain_threshold=8.0, max_gap=None, max_turn=30.0, max_chain=8, **kw):
    """-> list of trail dicts (fields as ogfkit.trails.detect_trails + group, curved, chain_z, n_seg) for trails that are curved or broken into several
    segments and were not found as straight trails.  Tile-wise detection (tile x tile px windows, 50 % overlap, threshold `tile_threshold`, minimum
    length 0.45 tile) -> segments; segments are linked when the gap between their ends is <= max_gap (default 0.6 tile) and the direction changes by
    <= max_turn degrees; a chain is accepted when the Stouffer combination sum(z)/sqrt(n) >= chain_threshold and it has >= 2 segments."""
    img = np.asarray(img, np.float32)
    ny, nx = img.shape
    valid = (np.isfinite(img) & (img != 0)) if valid is None else valid
    tile = int(tile or max(200, min(nx, ny) // 2))
    tile = min(tile, nx, ny)
    max_gap = 0.6 * tile if max_gap is None else max_gap
    kw.setdefault('n_null', 60)
    stride = max(tile // 2, 1)
    ys = list(range(0, max(ny - tile, 0) + 1, stride)); xs = list(range(0, max(nx - tile, 0) + 1, stride))
    if ys[-1] != ny - tile:
        ys.append(max(ny - tile, 0))
    if xs[-1] != nx - tile:
        xs.append(max(nx - tile, 0))
    segs = []
    for y0 in ys:
        for x0 in xs:
            r = T.detect_trails(img[y0:y0 + tile, x0:x0 + tile], valid=valid[y0:y0 + tile, x0:x0 + tile], threshold=tile_threshold,
                                min_length=0.45 * tile, max_trails=3, **kw)
            for t in r['trails']:
                segs.append(_seg_from_tile(t, x0, y0, nx, ny))
    # drop what the straight detector already covers
    if straight:
        keep = []
        for t in segs:
            ends = [(t['x1'], t['y1']), (t['x2'], t['y2']), ((t['x1'] + t['x2']) / 2, (t['y1'] + t['y2']) / 2)]
            cov = 0
            for (x, y) in ends:
                if any(T.seg_distance(np.array([x - 1.0]), np.array([y - 1.0]), s)[0] <= s['halfwidth'] + 4 for s in straight):
                    cov += 1
            if cov < 3:
                keep.append(t)
        segs = keep
    segs = _merge_dups(segs, nx, ny)
    n = len(segs)
    if n < 2:
        return []
    # graph
    def ends(t):
        return np.array([[t['x1'], t['y1']], [t['x2'], t['y2']]], float)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            Ei, Ej = ends(segs[i]), ends(segs[j])
            best = None
            for a in (0, 1):
                for b in (0, 1):
                    g = np.hypot(*(Ei[a] - Ej[b]))
                    if best is None or g < best[0]:
                        best = (g, a, b)
            g, a, b = best
            if g > max_gap:
                continue
            ui = Ei[a] - Ei[1 - a]; uj = Ej[b] - Ej[1 - b]            # outward directions at the joint
            ui /= max(np.hypot(*ui), 1e-9); uj /= max(np.hypot(*uj), 1e-9)
            turn = math.degrees(math.acos(float(np.clip(np.dot(ui, -uj), -1, 1))))
            ok = turn <= max_turn
            if g > 6:
                gv = Ej[b] - Ei[a]; gv /= max(np.hypot(*gv), 1e-9)
                ok = ok and math.degrees(math.acos(float(np.clip(np.dot(ui, gv), -1, 1)))) <= max_turn and \
                    math.degrees(math.acos(float(np.clip(np.dot(uj, -gv), -1, 1)))) <= max_turn
            if ok:
                adj[i].append((j, turn)); adj[j].append((i, turn))
    seen, chains = set(), []
    for i in range(n):
        if i in seen:
            continue
        comp, stack = [], [i]
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k); comp.append(k)
            stack.extend(j for j, _ in adj[k] if j not in seen)
        if len(comp) >= 2:
            chains.append(comp)
    out = []
    gid = 0
    for comp in sorted(chains, key=lambda c: -sum(segs[k]['zscore'] for k in c)):
        z = sum(segs[k]['zscore'] for k in comp) / math.sqrt(len(comp))
        if z < chain_threshold or len(comp) > max_chain:
            continue
        gid += 1
        ths = [segs[k]['_th'] for k in comp]
        turn = max(abs((a - b + 90) % 180 - 90) for a in ths for b in ths)
        for k in comp:
            t = {kk: v for kk, v in segs[k].items() if not kk.startswith('_')}
            t.update(theta_deg=float(segs[k]['_th']), rho=float(segs[k]['_rho']), s0=float(segs[k]['_s0']), s1=float(segs[k]['_s1']),
                     group=gid, curved=bool(turn >= 3.0), chain_z=float(z), n_seg=len(comp), max_turn_deg=float(turn))
            out.append(t)
    return out


# ------------------------------------------------------------------ contaminating flux
def trail_model_at(t, X, Y):
    """fitted trail model (peak amp x box (x) Gaussian profile x soft ends) at 0-based pixel coordinates"""
    ax, ay, bx, by = t['x1'] - 1, t['y1'] - 1, t['x2'] - 1, t['y2'] - 1
    L = max(math.hypot(bx - ax, by - ay), 1e-9)
    ux, uy = (bx - ax) / L, (by - ay) / L
    s = (X - ax) * ux + (Y - ay) * uy
    u = -(X - ax) * uy + (Y - ay) * ux
    a, sb = t.get('box_halfwidth', 2.0), max(t.get('blur', 1.2), 0.3)
    prof = 0.5 * (erf((u + a) / (math.sqrt(2) * sb)) - erf((u - a) / (math.sqrt(2) * sb)))
    ends = 0.5 * (1 + erf(s / (math.sqrt(2) * sb))) * 0.5 * (1 + erf((L - s) / (math.sqrt(2) * sb)))
    return t.get('amp', 0.0) * prof * ends


def trail_flux(cat, trails, k=2.5):
    """-> (TRAIL_FLUX, TRAIL_FRAC): sum of the trail model over each object's ellipse (semi-axes k*A, k*B, angle THETA_IMAGE, SExtractor convention), in
    image units, and that divided by the object's FLUX_AUTO (or FLUX_ISO) when the catalogue has it, else -1.  This is the light the trail adds inside the
    aperture before any local-sky subtraction (the catalogue's background mesh removes part of a wide faint trail, so it is an upper bound there)."""
    n = len(cat['X_IMAGE'])
    x, y = np.asarray(cat['X_IMAGE'], float) - 1, np.asarray(cat['Y_IMAGE'], float) - 1
    A = np.asarray(cat['A_IMAGE'], float) if 'A_IMAGE' in cat else np.full(n, 1.5)
    B = np.asarray(cat['B_IMAGE'], float) if 'B_IMAGE' in cat else A
    th = np.radians(np.asarray(cat['THETA_IMAGE'], float)) if 'THETA_IMAGE' in cat else np.zeros(n)
    flux = np.zeros(n)
    for t in trails:
        reach = t['halfwidth'] + 3
        d = T.seg_distance(x, y, t)
        for i in np.where(d <= k * np.maximum(A, B) + reach)[0]:
            r = int(math.ceil(k * max(A[i], B[i]))) + 1
            X, Y = np.meshgrid(np.arange(int(x[i]) - r, int(x[i]) + r + 1), np.arange(int(y[i]) - r, int(y[i]) + r + 1))
            c, s_ = math.cos(th[i]), math.sin(th[i])
            dx, dy = X - x[i], Y - y[i]
            u = dx * c + dy * s_; v = -dx * s_ + dy * c
            inside = (u / (k * max(A[i], 0.5))) ** 2 + (v / (k * max(B[i], 0.5))) ** 2 <= 1.0
            flux[i] += float(trail_model_at(t, X, Y)[inside].sum())
    ref = None
    for c in ('FLUX_AUTO', 'FLUX_ISO'):
        if c in cat:
            ref = np.asarray(cat[c], float); break
    with np.errstate(all='ignore'):
        frac = flux / np.where(ref > 0, ref, np.nan) if ref is not None else np.full(n, -1.0)
    frac = np.where(np.isfinite(frac), frac, -1.0)
    return flux, frac


# ------------------------------------------------------------------ WCS registration of dithered frames
def wcs_resample(data, wcs_in, wcs_out, shape_out, order=1, cval=np.nan, min_weight=0.5):
    """resample `data` (on wcs_in) to the pixel grid of wcs_out with scipy.ndimage.map_coordinates (order 1 = bilinear; masks: use order 1 on a float mask and
    threshold).  Output pixels outside the input footprint get `cval`.  NaN input pixels (bad / masked) are excluded by normalised interpolation:
    out = I(data * valid) / I(valid), so a good pixel next to a bad one keeps its value (before, the bad pixel entered the bilinear sum as 0 and pulled
    its neighbours down); output pixels whose valid interpolation weight is < min_weight are `cval`."""
    from scipy.ndimage import map_coordinates
    ny, nx = shape_out
    Y, X = np.mgrid[0:ny, 0:nx]
    ra, dec = wcs_out.celestial.all_pix2world(X.ravel(), Y.ravel(), 0)
    xi, yi = wcs_in.celestial.all_world2pix(ra, dec, 0)
    coords = np.vstack([yi, xi])
    d = np.asarray(data, np.float64)
    good = np.isfinite(d)
    inside = (xi >= -0.5) & (xi <= d.shape[1] - 0.5) & (yi >= -0.5) & (yi <= d.shape[0] - 0.5)
    if good.all():
        out = map_coordinates(d, coords, order=order, mode='nearest')
    else:
        num = map_coordinates(np.where(good, d, 0.0), coords, order=order, mode='nearest')
        wt = map_coordinates(good.astype(np.float64), coords, order=order, mode='nearest')
        if order > 1:                                   # spline weights can overshoot: fall back to the nearest-pixel validity
            wt = np.clip(wt, 0.0, 1.0)
        with np.errstate(all='ignore'):
            out = num / wt
        inside &= wt >= min_weight
    out = np.where(inside, out, cval)
    return out.reshape(ny, nx)


def register_masks(mask, wcs_in, wcs_out, shape_out, grow=1.0):
    """boolean trail mask on the output grid (conservative: any output pixel whose bilinear mask value > 0.05; areas outside the input footprint are
    reported as 'no data' in the second return value)"""
    f = wcs_resample(mask.astype(np.float32), wcs_in, wcs_out, shape_out, order=1, cval=-1.0)
    nodata = f < 0
    m = f > 0.05
    if grow:
        from scipy.ndimage import binary_dilation
        m = binary_dilation(m, iterations=int(round(grow)))
    return m, nodata
