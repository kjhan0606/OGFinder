#!/usr/bin/env python3
"""Validation of ogfkit.trails on real HUDF F160W and M51 pixels with injected synthetic trails.
Experiments (python3 trails_validate.py --exp detect,fp,cross,phot,stack --out DIR [--quick] [--workers N]):
  detect  detection rate / position / width accuracy vs trail peak amplitude (in sky sigma, and mag/arcsec^2 for HUDF), widths, full+partial
  fp      false trails on trail-free HUDF windows, M51 variants and synthetic galaxy fields (best z-score distribution)
  cross   trails through the brightest galaxies (HUDF objects, M51 core)
  phot    aperture-photometry bias of objects near the trail: unmasked / masked / interpolated / interpolated+noise
  stack   5-frame stacks (mean/median/sigma-clip) with a different trail per frame, with and without trail masks
Writes DIR/results.json and DIR/trails_validation.md."""
import argparse, json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
from ogfkit import trails as T, trailsynth as S, trails_ext as X

HUDF_FULL = os.environ.get('TRAILS_HUDF', '/workspace/fits/hudf_f160w.fits')
M51 = os.environ.get('TRAILS_M51', '/workspace/fits/m51.fits')
ZP, PIX = 25.94, 0.06          # HUDF12 F160W zero point (AB) and pixel scale
_cache = {}


def fits_data(path):
    if path not in _cache:
        from astropy.io import fits
        _cache[path] = np.asarray(fits.getdata(path), np.float32)
    return _cache[path]


def hudf_windows(size, n, seed=11):
    d = fits_data(HUDF_FULL)
    v = (d != 0).astype(np.float32)
    ii = np.cumsum(np.cumsum(np.pad(v, ((1, 0), (1, 0))), 0), 1)
    pos = []
    for y in range(0, d.shape[0] - size, 100):
        for x in range(0, d.shape[1] - size, 100):
            if (ii[y + size, x + size] - ii[y, x + size] - ii[y + size, x] + ii[y, x]) > 0.9995 * size * size:
                pos.append((y, x))
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(pos))[:n]
    return [pos[i] for i in idx]


def get_field(spec):
    kind = spec['field']
    if kind == 'hudf':
        y, x = spec['win']; s = spec['size']
        return fits_data(HUDF_FULL)[y:y + s, x:x + s].copy()
    if kind == 'm51':
        im = fits_data(M51).copy()
        v = spec.get('variant', 0)
        if v & 1: im = im[::-1]
        if v & 2: im = im[:, ::-1]
        if v & 4: im = im.T
        return np.ascontiguousarray(im)
    if kind == 'synth':
        return S.test_field(spec['size'], spec['seed'], n_gal=spec.get('n_gal', 14), corr=spec.get('corr', 0.0))[0].astype(np.float32)
    raise ValueError(kind)


def sigma_of(img):
    key = ('sig', id(img), img.shape, float(img[0, 0]), float(img[-1, -1]))
    z, sp, b = T.standardise(img, np.isfinite(img) & (img != 0))
    return sp


def model_fwhm(w, blur):
    t = np.linspace(-60, 60, 24001)
    from scipy.special import erf
    p = 0.5 * (erf((t + w / 2) / (math.sqrt(2) * blur)) - erf((t - w / 2) / (math.sqrt(2) * blur)))
    h = t[p >= 0.5]
    return float(h.max() - h.min())


def line_extent(theta, rho, nx, ny):
    th = math.radians(theta)
    s = np.linspace(-math.hypot(nx, ny), math.hypot(nx, ny), 8001)
    x = (nx - 1) / 2 + s * math.cos(th) + rho * (-math.sin(th))
    y = (ny - 1) / 2 + s * math.sin(th) + rho * math.cos(th)
    ok = (x >= 0) & (x <= nx - 1) & (y >= 0) & (y <= ny - 1)
    return (float(s[ok].min()), float(s[ok].max())) if ok.any() else (0.0, 0.0)


def rho_through(theta, x, y, nx, ny):
    th = math.radians(theta)
    return -(x - (nx - 1) / 2) * math.sin(th) + (y - (ny - 1) / 2) * math.cos(th)


def match(trails, theta, rho, w):
    best = None
    for t in trails:
        dth = (t['theta_deg'] - theta + 90) % 180 - 90
        r = t['rho']
        if abs(dth) > 90 - 1e-9:
            pass
        # 0 and 180 deg are the same line with rho -> -rho: use the end points instead
        d = abs(dth)
        rr = r if abs((t['theta_deg'] - theta) % 360) < 90 or abs((t['theta_deg'] - theta) % 360) > 270 else -r
        if d < 1.5 and abs(rr - rho) < max(3.0, w):
            if best is None or t['zscore'] > best['zscore']:
                best = dict(t, dtheta=(t['theta_deg'] - theta + 90) % 180 - 90, drho=rr - rho)
    return best


def inject(img, sig, theta, rho, w, k, partial, rng, blur=1.2):
    ny, nx = img.shape
    s0 = s1 = None
    lo, hi = line_extent(theta, rho, nx, ny)
    if partial:
        L = hi - lo
        ln = rng.uniform(0.4, 0.6) * L
        c = rng.uniform(lo + ln / 2, hi - ln / 2)
        s0, s1 = c - ln / 2, c + ln / 2
    else:
        s0, s1 = lo - 5, hi + 5
    tr = S.trail_image(img.shape, theta, rho, w, k * sig, s0, s1, blur)
    return img + tr, tr, (max(s0, lo), min(s1, hi))


def mask_stats(trails, tr, sig, shape):
    m = T.trail_mask(shape, trails)
    vis = tr > 0.25 * sig
    cov = float((m & vis).sum() / max(vis.sum(), 1))
    return m, cov, float(m.sum() / max(vis.sum(), 1))


def run_detect(spec, cross=False):
    rng = np.random.default_rng(spec['seed'])
    img = get_field(spec)
    sig = sigma_of(img)
    ny, nx = img.shape
    theta = rng.uniform(0, 180)
    if spec.get('through') is not None:
        x0, y0 = spec['through']
        rho = rho_through(theta, x0, y0, nx, ny)
    else:
        rho = rng.uniform(-0.3, 0.3) * min(nx, ny)
    w, k = spec['w'], spec['k']
    im2, tr, (a0, a1) = inject(img, sig, theta, rho, w, k, spec.get('partial', False), rng)
    t0 = time.time()
    res = T.detect_trails(im2)
    m = match(res['trails'], theta, rho, w)
    out = dict(spec=spec, theta=theta, rho=rho, sigma=float(sig), k=k, w=w, n_det=len(res['trails']), found=m is not None,
               best_z=res.get('best_zscore'), sec=time.time() - t0)
    if m is not None:
        dth = m['dtheta']; fw = model_fwhm(w, 1.2)
        ov = max(0.0, min(m['s1'], a1) - max(m['s0'], a0)); un = max(m['s1'], a1) - min(m['s0'], a0)
        _, cov, ratio = mask_stats([m], tr, sig, img.shape)
        out.update(dtheta=dth, drho=m['drho'], fwhm_ratio=m['fwhm'] / fw, iou=ov / max(un, 1e-9), z=m['zscore'], cover=cov, mask_ratio=ratio,
                   len_ratio=(m['s1'] - m['s0']) / max(a1 - a0, 1e-9))
    out['extra'] = len(res['trails']) - (1 if m is not None else 0)
    return out


def run_fp(spec):
    img = get_field(spec)
    t0 = time.time()
    res = T.detect_trails(img)
    return dict(spec=spec, n_det=len(res['trails']), best_z=res.get('best_zscore'), z=[t['zscore'] for t in res['trails']], sec=time.time() - t0,
                rejected=len(res.get('rejected', [])))


# ---- photometry ---------------------------------------------------------------------------------
def objects(img, sig, snr_min):
    import sep
    d = np.ascontiguousarray(img, np.float32)
    bk = sep.Background(d, bw=32, bh=32)
    o = sep.extract(d - bk.back(), 4.0, err=sig, minarea=8)
    ap = 5.0
    F, _, _ = sep.sum_circle(d - bk.back(), o['x'], o['y'], ap)
    ok = F / (sig * math.sqrt(math.pi * ap * ap)) > snr_min
    edge = (o['x'] > 25) & (o['x'] < img.shape[1] - 26) & (o['y'] > 25) & (o['y'] < img.shape[0] - 26)
    sel = ok & edge
    return np.c_[o['x'][sel], o['y'][sel]], F[sel]


def aper(img, x, y, mask=None, R=5.0, r1=8.0, r2=14.0, min_frac=0.4):
    x0, y0 = int(round(x)), int(round(y)); h = int(r2) + 1
    sl = (slice(max(0, y0 - h), y0 + h + 1), slice(max(0, x0 - h), x0 + h + 1))
    Y, X = np.mgrid[sl]
    r = np.hypot(X - x, Y - y)
    sub = img[sl].astype(np.float64)
    ok = np.ones(sub.shape, bool) if mask is None else ~mask[sl]
    ann = (r >= r1) & (r <= r2) & ok & np.isfinite(sub)
    ap = (r <= R)
    if ann.sum() < 30:
        return np.nan
    sky = np.median(sub[ann])
    n_ok = (ap & ok).sum()
    if n_ok < min_frac * ap.sum():
        return np.nan
    return float(((sub - sky)[ap & ok]).sum() * ap.sum() / n_ok)


def run_phot(spec):
    import warnings; warnings.simplefilter('ignore')
    rng = np.random.default_rng(spec['seed'])
    img = get_field(spec); sig = sigma_of(img); ny, nx = img.shape
    pos, F = objects(img, sig, 8.0)
    if len(pos) < 5:
        return dict(spec=spec, rows=[], found=False)
    big = pos[F > np.percentile(F, 60)]
    c = big[rng.integers(len(big))]
    theta = rng.uniform(0, 180); rho = rho_through(theta, c[0], c[1], nx, ny)
    w, k = spec['w'], spec['k']
    im2, tr, ext = inject(img, sig, theta, rho, w, k, False, rng)
    res = T.detect_trails(im2)
    m = match(res['trails'], theta, rho, w)
    mask = T.trail_mask(img.shape, res['trails']) if res['trails'] else np.zeros(img.shape, bool)
    fi = T.fill_interpolate(im2, res['trails'], mask, noise=False) if res['trails'] else im2
    fn = T.fill_interpolate(im2, res['trails'], mask, noise=True, seed=spec['seed']) if res['trails'] else im2
    fc = T.fill_interpolate(img, res['trails'], mask, noise=False) if res['trails'] else img      # controls: same treatment on the trail-free image
    th = math.radians(theta)
    hw = float(m.get('halfwidth', w / 2 + 4.4)) if m else w / 2 + 4.4
    rows = []
    sap = sig * math.sqrt(math.pi * 25)
    for (x, y), f0 in zip(pos, F):
        d = abs(-(x - (nx - 1) / 2) * math.sin(th) + (y - (ny - 1) / 2) * math.cos(th) - rho)
        if d > 45:
            continue
        base = aper(img, x, y)
        row = dict(d=float(d), f0=float(base), snr=float(base / sap), unmasked=aper(im2, x, y) - base,
                   masked=aper(im2, x, y, mask) - base, interp=aper(fi, x, y) - base, interp_noise=aper(fn, x, y) - base, masked_ctrl=aper(img, x, y, mask) - base, interp_ctrl=aper(fc, x, y) - base,
                   grp='on' if d < hw else ('edge' if d < hw + 5.0 else 'near'), sap=float(sap))
        rows.append(row)
    return dict(spec=spec, rows=rows, found=m is not None, n_det=len(res['trails']))


def run_stack(spec):
    rng = np.random.default_rng(spec['seed'])
    img = get_field(spec); sig = sigma_of(img); ny, nx = img.shape
    pos, F = objects(img, sig, 15.0)
    nf, sf, w, k = 5, spec['sf'] * sig, spec['w'], spec['k']
    noise = [rng.normal(0, sf, img.shape).astype(np.float32) for _ in range(nf)]
    clean = [img + n for n in noise]
    frames, trs, masks, found = [], [], [], 0
    for i in range(nf):
        theta = rng.uniform(0, 180); rho = rng.uniform(-0.3, 0.3) * min(nx, ny)
        if len(pos) and i < 3:       # make sure some trails cross objects
            c = pos[rng.integers(len(pos))]; rho = rho_through(theta, c[0], c[1], nx, ny)
        tr = S.trail_image(img.shape, theta, w, 0, 0)  if False else S.trail_image(img.shape, theta, rho, w, k * sf, None, None, 1.2)
        fr = clean[i] + tr
        res = T.detect_trails(fr)
        mk = T.trail_mask(img.shape, res['trails']) if res['trails'] else np.zeros(img.shape, bool)
        found += 1 if match(res['trails'], theta, rho, w) is not None else 0
        frames.append(fr); masks.append(mk); trs.append(tr)
    ref = {'mean': T.stack_frames(clean, None, 'mean')[0], 'median': T.stack_frames(clean, None, 'median')[0],
           'sigclip': T.stack_frames(clean, None, 'sigclip')[0]}
    meth = {'mean_plain': ('mean', False), 'median_plain': ('median', False), 'sigclip_plain': ('sigclip', False),
            'mean_masked': ('mean', True), 'median_masked': ('median', True), 'sigclip_masked': ('sigclip', True)}
    stk = {n: T.stack_frames(frames, masks if mk else None, m)[0] for n, (m, mk) in meth.items()}
    tracks = np.any([t > 0.25 * sf for t in trs], axis=0)
    rows = []
    sst = sf / math.sqrt(nf) * math.sqrt(math.pi * 25)
    for (x, y), f0 in zip(pos, F):
        ix, iy = int(round(x)), int(round(y))
        near = tracks[max(0, iy - 5):iy + 6, max(0, ix - 5):ix + 6].any()
        if not near:
            continue
        row = dict(snr=float(f0 / (sig * math.sqrt(math.pi * 25))), sst=float(sst), f0=float(f0))
        for n, (m, mk) in meth.items():
            row[n] = aper(stk[n], x, y) - aper(ref[m], x, y)
        rows.append(row)
    pix = {}
    for n, (m, mk) in meth.items():
        d = (stk[n] - ref[m])[tracks]
        d = d[np.isfinite(d)]
        pix[n] = [float(np.mean(d) / (sf / math.sqrt(nf))) if len(d) else None, float(np.median(d) / (sf / math.sqrt(nf))) if len(d) else None]
    return dict(spec=spec, rows=rows, found=found, nf=nf, pix=pix)



# ---- extensions: curved / flickering trails, contaminating flux ----------------------------------
def run_curved(spec):
    import warnings; warnings.simplefilter('ignore')
    rng = np.random.default_rng(spec['seed'])
    img = get_field(spec); sig = sigma_of(img); ny, nx = img.shape
    th0 = rng.uniform(-35, 35) + (0 if rng.random() < 0.5 else 0)
    y0 = rng.uniform(0.25, 0.75) * ny
    flick = spec.get('duty')
    tr, (cx, cy, cs) = S.curved_trail_image(img.shape, (-5.0, y0), th0, spec['curv'], 1.15 * nx / max(math.cos(math.radians(th0)), 0.5), spec['w'], spec['k'] * sig,
                                           duty=flick, period=spec.get('period', 60.0))
    im2 = img + tr
    res = T.detect_trails(im2)
    cur = X.detect_curved(im2, res['trails']) if spec.get('curved', True) else []
    vis = tr > 0.25 * sig
    # truth footprint (all points of the curve, incl. gaps of a flickering trail)
    foot = np.zeros(img.shape, bool)
    ix, iy = np.round(cx).astype(int), np.round(cy).astype(int)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    foot[iy[ok], ix[ok]] = True
    from scipy.ndimage import binary_dilation, distance_transform_edt
    dist = distance_transform_edt(~foot)
    out = dict(spec=spec, n_straight=len(res['trails']), n_curved=len(cur), sigma=float(sig))
    for name, tl in (('straight', res['trails']), ('both', res['trails'] + cur)):
        m = T.trail_mask(img.shape, tl) if tl else np.zeros(img.shape, bool)
        out['cover_' + name] = float((m & vis).sum() / max(vis.sum(), 1))
        out['footcover_' + name] = float((m & (dist <= spec['w'] / 2)).sum() / max((dist <= spec['w'] / 2).sum(), 1))
        out['stray_' + name] = int((m & (dist > spec['w'] / 2 + 12)).sum())
        out['maskarea_' + name] = int(m.sum())
    out['vis'] = int(vis.sum())
    out['found_any'] = out['cover_both'] > 0.5
    if flick and res['trails']:
        best = max(res['trails'], key=lambda t: t['length'])
        ap = X.along_profile(im2, best)
        out['duty_meas'] = ap['duty']; out['gaps_meas'] = ap['n_gaps']
    return out


def run_flux(spec):
    import warnings; warnings.simplefilter('ignore')
    import sep
    rng = np.random.default_rng(spec['seed'])
    img = get_field(spec); sig = sigma_of(img); ny, nx = img.shape
    d = np.ascontiguousarray(img, np.float32)
    bk = sep.Background(d, bw=32, bh=32)
    o = sep.extract(d - bk.back(), 4.0, err=sig, minarea=8)
    if len(o) < 5:
        return dict(spec=spec, rows=[])
    big = o[np.argsort(-o['flux'])[:max(5, len(o) // 4)]]
    c = big[rng.integers(len(big))]
    theta = rng.uniform(0, 180); rho = rho_through(theta, c['x'], c['y'], nx, ny)
    w, k = spec['w'], spec['k']
    im2, tr, ext = inject(img, sig, theta, rho, w, k, False, rng)
    res = T.detect_trails(im2)
    cat = dict(NUMBER=list(range(len(o))), X_IMAGE=o['x'] + 1, Y_IMAGE=o['y'] + 1, A_IMAGE=o['a'], B_IMAGE=o['b'], THETA_IMAGE=np.degrees(o['theta']),
               FLUX_AUTO=o['flux'])
    rows = []
    if res['trails']:
        fl, fr = X.trail_flux(cat, res['trails'], k=2.5)
        Y, Xg = np.mgrid[0:ny, 0:nx]
        for i in np.where(fl != 0)[0]:
            r = int(math.ceil(2.5 * max(o['a'][i], o['b'][i]))) + 1
            x0, y0 = int(o['x'][i]), int(o['y'][i])
            sl = (slice(max(y0 - r, 0), min(y0 + r + 1, ny)), slice(max(x0 - r, 0), min(x0 + r + 1, nx)))
            yy, xx = np.mgrid[sl]
            cth, sth = math.cos(o['theta'][i]), math.sin(o['theta'][i])
            u = (xx - o['x'][i]) * cth + (yy - o['y'][i]) * sth; v = -(xx - o['x'][i]) * sth + (yy - o['y'][i]) * cth
            ins = (u / (2.5 * max(o['a'][i], 0.5))) ** 2 + (v / (2.5 * max(o['b'][i], 0.5))) ** 2 <= 1
            truth = float(tr[sl][ins].sum())
            rows.append(dict(model=float(fl[i]), truth=truth, frac=float(fr[i]), flux=float(o['flux'][i]), sap=float(sig * math.sqrt(ins.sum()))))
    return dict(spec=spec, rows=rows, found=bool(match(res['trails'], theta, rho, w)))


def dispatch(spec):
    import warnings; warnings.simplefilter('ignore')
    try:
        e = spec['exp']
        if e in ('detect', 'cross'):
            return run_detect(spec)
        if e == 'fp':
            return run_fp(spec)
        if e == 'phot':
            return run_phot(spec)
        if e == 'stack':
            return run_stack(spec)
        if e == 'curved':
            return run_curved(spec)
        if e == 'flux':
            return run_flux(spec)
    except Exception as ex:      # keep the run going, count it
        import traceback
        return dict(spec=spec, error=traceback.format_exc()[-400:])


# ---- campaigns -----------------------------------------------------------------------------------
def campaign(exp, quick):
    sp = []
    nrep = 2 if quick else 10
    wins = hudf_windows(700, 12 if quick else 40)
    if exp == 'detect':
        ks = [0.15, 0.3, 0.6, 1.2, 2.5] if quick else [0.1, 0.15, 0.2, 0.3, 0.45, 0.7, 1.0, 2.0, 4.0]
        i = 0
        for field in ('hudf', 'm51'):
            for w in ((6,) if quick else (3, 6, 12, 24)):
                for k in ks:
                    for r in range(nrep):
                        i += 1
                        s = dict(exp='detect', field=field, w=w, k=k, seed=1000 + i, partial=bool(r % 2), size=700)
                        if field == 'hudf':
                            s['win'] = wins[i % len(wins)]
                        sp.append(s)
    elif exp == 'fp':
        for i, wn in enumerate(wins):
            sp.append(dict(exp='fp', field='hudf', win=wn, size=700, tag='hudf'))
        for v in range(2 if quick else 8):
            sp.append(dict(exp='fp', field='m51', variant=v, tag='m51'))
        for i in range(4 if quick else 40):
            sp.append(dict(exp='fp', field='synth', size=700, seed=500 + i, corr=[0.0, 1.5][i % 2], n_gal=14 if i % 4 < 2 else 40, tag='synth'))
    elif exp == 'cross':
        from ogfkit import trails as _t
        n = 4 if quick else 24
        for i, wn in enumerate(wins[:n]):
            img = fits_data(HUDF_FULL)[wn[0]:wn[0] + 700, wn[1]:wn[1] + 700]
            pos, F = objects(img, sigma_of(img), 20.0)
            if len(pos) == 0:
                continue
            c = pos[int(np.argmax(F))]
            for k in (0.3, 0.6, 1.2):
                sp.append(dict(exp='cross', field='hudf', win=wn, size=700, w=6, k=k, seed=3000 + i * 7 + int(k * 10), through=[float(c[0]), float(c[1])], tag='hudf'))
        for i in range(4 if quick else 24):
            for k in (0.3, 0.6, 1.2):
                sp.append(dict(exp='cross', field='m51', w=6, k=k, seed=4000 + i * 7 + int(k * 10), through=[300.0 + (i % 5) * 4, 300.0 + (i % 3) * 4], variant=0, tag='m51'))
    elif exp == 'phot':
        for i, wn in enumerate(wins[: (6 if quick else 40)]):
            for k in (0.3, 1.0, 3.0):
                sp.append(dict(exp='phot', field='hudf', win=wn, size=700, w=6, k=k, seed=5000 + i * 5 + int(k * 3)))
    elif exp == 'curved':
        i = 0
        for field in ('hudf',):
            for curv in (0.0, 4.0, 8.0):
                for k in (1.5, 3.0, 6.0):
                    for r in range(3 if quick else 8):
                        i += 1
                        sp.append(dict(exp='curved', field=field, win=wins[i % len(wins)], size=700, w=6, k=k, curv=curv, seed=7000 + i, tag='arc'))
        for duty in (0.3, 0.5, 0.7):
            for k in (1.5, 3.0, 6.0):
                for r in range(3 if quick else 8):
                    i += 1
                    sp.append(dict(exp='curved', field='hudf', win=wins[i % len(wins)], size=700, w=6, k=k, curv=0.0, duty=duty, period=60.0, seed=7000 + i, tag='flicker'))
    elif exp == 'flux':
        for i, wn in enumerate(wins[: (6 if quick else 40)]):
            for k in (1.0, 3.0):
                sp.append(dict(exp='flux', field='hudf', win=wn, size=700, w=6, k=k, seed=8000 + i * 3 + int(k)))
    elif exp == 'stack':
        for i, wn in enumerate(wins[: (4 if quick else 24)]):
            for k in (1.0, 3.0):
                sp.append(dict(exp='stack', field='hudf', win=wn, size=700, w=6, k=k, sf=4.0, seed=6000 + i * 3 + int(k)))
    return sp


# ---- summaries -----------------------------------------------------------------------------------
def pct(a, p):
    a = np.asarray([x for x in a if x is not None and np.isfinite(x)], float)
    return float(np.percentile(a, p)) if len(a) else float('nan')


def rms(a):
    a = np.asarray([x for x in a if np.isfinite(x)], float)
    return float(np.sqrt(np.mean(a * a))) if len(a) else float('nan')


def sb_hudf(k, sig):       # peak surface brightness mag/arcsec^2 for an amplitude of k sky sigma (e/s per pixel)
    return ZP - 2.5 * math.log10(k * sig / (PIX * PIX))


def summarize(exp, R, md):
    R = [r for r in R if r and 'error' not in r]
    if exp in ('detect', 'cross'):
        md.append('\n### %s: detection rate (found/total) by peak amplitude in sky sigma\n' % ('trails through bright galaxies' if exp == 'cross' else 'trails injected into real pixels'))
        for field in ('hudf', 'm51'):
            sub = [r for r in R if r['spec']['field'] == field]
            if not sub:
                continue
            ws = sorted(set(r['w'] for r in sub)); ks = sorted(set(r['k'] for r in sub))
            sg = np.median([r['sigma'] for r in sub])
            md.append('\n**%s** (sky sigma %.4g/px%s); widths = box FWHM px; partial+full mixed\n' % (field.upper(), sg, ', peak mu = ' + ', '.join('%.1f' % sb_hudf(k, sg) for k in ks) + ' mag/arcsec^2 (HUDF F160W AB, 0.06"/px) for the columns in order' if field == 'hudf' else ''))
            md.append('| width \\ amp[sigma] | ' + ' | '.join('%g' % k for k in ks) + ' |')
            md.append('|---|' + '---|' * len(ks))
            for w in ws:
                cells = []
                for k in ks:
                    c = [r for r in sub if r['w'] == w and abs(r['k'] - k) < 1e-9]
                    cells.append('%d/%d' % (sum(r['found'] for r in c), len(c)) if c else '-')
                md.append('| %g | %s |' % (w, ' | '.join(cells)))
        det = [r for r in R if r['found']]
        md.append('\nAccuracy over %d detected trails: dtheta median %.3f deg, 68%% half-range %.3f; drho median %.2f px, RMS %.2f px; FWHM(meas)/FWHM(true) median %.2f (16-84%%: %.2f-%.2f); '
                  'length IoU median %.2f; mask covers %.3f of the visible trail (>0.25 sigma) pixels (median), mask area / visible area = %.2f; extra (unmatched) detections in injected images: %d.'
                  % (len(det), pct([r['dtheta'] for r in det], 50), 0.5 * (pct([abs(r['dtheta']) for r in det], 84)), pct([r['drho'] for r in det], 50), rms([r['drho'] for r in det]),
                     pct([r['fwhm_ratio'] for r in det], 50), pct([r['fwhm_ratio'] for r in det], 16), pct([r['fwhm_ratio'] for r in det], 84), pct([r['iou'] for r in det], 50),
                     pct([r['cover'] for r in det], 50), pct([r['mask_ratio'] for r in det], 50), sum(r['extra'] for r in R)))
        for w in sorted(set(r['w'] for r in det)):
            d = [r for r in det if r['w'] == w and r['k'] >= 0.5]
            if d:
                md.append('- w=%g (amp>=0.5 sigma, n=%d): dtheta RMS %.3f deg, drho RMS %.2f px, FWHM ratio median %.2f, coverage median %.3f' % (w, len(d), rms([r['dtheta'] for r in d]), rms([r['drho'] for r in d]), pct([r['fwhm_ratio'] for r in d], 50), pct([r['cover'] for r in d], 50)))
    if exp == 'fp':
        md.append('\n### false trails on trail-free images (default threshold)\n')
        md.append('| set | images | with detection | detections | median best z | 95% best z | max best z | shape-rejected candidates |')
        md.append('|---|---|---|---|---|---|---|---|')
        for tag in ('hudf', 'm51', 'synth'):
            sub = [r for r in R if r['spec'].get('tag') == tag]
            if sub:
                bz = [r['best_z'] for r in sub]
                md.append('| %s | %d | %d | %d | %.1f | %.1f | %.1f | %d |' % (tag, len(sub), sum(r['n_det'] > 0 for r in sub), sum(r['n_det'] for r in sub), pct(bz, 50), pct(bz, 95), pct(bz, 100), sum(r['rejected'] for r in sub)))
        allz = [r['best_z'] for r in R if r['best_z'] is not None]
        md.append('\nFraction of all %d trail-free images whose best z exceeds t: ' % len(allz) + ', '.join('t=%g: %.3f' % (t, np.mean(np.asarray(allz) >= t)) for t in (5, 6, 7, 8, 10)))
        md.append('Mean run time %.1f s/image.' % np.mean([r['sec'] for r in R]))
    if exp == 'phot':
        rows = {}
        for r in R:
            for row in r['rows']:
                rows.setdefault(r['spec']['k'], []).append(dict(row, found=r['found']))
        md.append('\n### aperture-photometry bias of objects near an injected trail (r=5 px aperture, local sky annulus, targets S/N>8, trail w=6 px through a bright object)\n')
        md.append('Bias = F(with trail, treatment) - F(original); the trail-free control applies the same mask / interpolation to the image without the trail (cost of the treatment itself: lost sky-annulus pixels, interpolated-over object light); units of the nominal aperture noise sigma_ap (sigma_pix*sqrt(N_ap)). on = object centre inside the masked band; edge = aperture reaches the band; near = band >5 px away, within 45 px.\n')
        md.append('| amp [sigma] | group | n objs | unmasked med / RMS | masked med / RMS (n lost) | masked, trail-free control | interp med / RMS | interp, trail-free control | interp+noise | trails found |')
        md.append('|---|---|---|---|---|---|---|---|---|---|')
        for k in sorted(rows):
            nt = len([r for r in R if r['spec']['k'] == k]); nf = sum(r['found'] for r in R if r['spec']['k'] == k)
            for grp in ('on', 'edge', 'near'):
                g = [x for x in rows[k] if x['grp'] == grp]
                if not g:
                    continue
                def cell(name):
                    v = np.array([x[name] / x['sap'] for x in g], float)
                    return '%+.2f / %.2f' % (np.nanmedian(v), rms(v)) if np.isfinite(v).any() else 'n/a'
                lost = sum(not np.isfinite(x['masked']) for x in g)
                md.append('| %g | %s | %d | %s | %s (%d) | %s | %s | %s | %s | %d/%d |' % (k, grp, len(g), cell('unmasked'), cell('masked'), lost, cell('masked_ctrl'), cell('interp'), cell('interp_ctrl'), cell('interp_noise'), nf, nt))
        md.append('')
        for k in sorted(rows):
            g = [x for x in rows[k] if x['grp'] != 'near' and x['f0'] > 0]
            if g:
                md.append('- amp %g sigma, on-trail+edge objects: fractional flux bias median (unmasked %+.3f, masked %+.3f, interpolated %+.3f)' % (
                    k, np.nanmedian([x['unmasked'] / x['f0'] for x in g]), np.nanmedian([x['masked'] / x['f0'] for x in g]), np.nanmedian([x['interp'] / x['f0'] for x in g])))
    if exp == 'curved':
        md.append('\n### curved and flickering trails injected into HUDF windows (default detection + `detect_curved`)\n')
        md.append('Coverage = fraction of the pixels where the trail exceeds 0.25 sigma that lie inside the mask; "straight" = ordinary detector only, "both" = plus the tile-chain curved detector; stray = masked pixels farther than w/2+12 px from the true path.\n')
        md.append('| set | curvature [deg/100px] / duty | amp [sigma] | n | cover straight (median) | cover both (median) | found (cover>0.5) straight / both | mean stray px both | duty measured (median) |')
        md.append('|---|---|---|---|---|---|---|---|---|')
        for tag in ('arc', 'flicker'):
            sub = [r for r in R if r['spec'].get('tag') == tag]
            keys = sorted(set((r['spec']['curv'] if tag == 'arc' else r['spec']['duty'], r['spec']['k']) for r in sub))
            for kk in keys:
                g = [r for r in sub if (r['spec']['curv'] if tag == 'arc' else r['spec']['duty'], r['spec']['k']) == kk]
                dm = [r['duty_meas'] for r in g if 'duty_meas' in r]
                md.append('| %s | %g | %g | %d | %.2f | %.2f | %d / %d | %.0f | %s |' % (tag, kk[0], kk[1], len(g), np.median([r['cover_straight'] for r in g]), np.median([r['cover_both'] for r in g]),
                          sum(r['cover_straight'] > 0.5 for r in g), sum(r['cover_both'] > 0.5 for r in g), np.mean([r['stray_both'] for r in g]), ('%.2f' % np.median(dm)) if dm else '-'))
    if exp == 'flux':
        rows = [dict(x, k=r['spec']['k']) for r in R for x in r['rows']]
        md.append('\n### contaminating flux per object (`trail_flux`, ellipse 2.5 A x 2.5 B) vs the true injected light in the same ellipse\n')
        if rows:
            m = np.array([x['model'] for x in rows]); t = np.array([x['truth'] for x in rows]); sap = np.array([x['sap'] for x in rows])
            sel = t > 3 * sap
            md.append('%d objects touched by a detected trail; %d with true contamination > 3 sigma_ap: model/truth median %.3f (16-84%%: %.3f-%.3f), RMS of (model - truth)/sigma_ap %.2f; for the rest RMS %.2f.' % (
                len(rows), sel.sum(), np.median(m[sel] / t[sel]), np.percentile(m[sel] / t[sel], 16), np.percentile(m[sel] / t[sel], 84), rms((m - t)[sel] / sap[sel]), rms((m - t)[~sel] / sap[~sel])))
            fr = np.array([x['frac'] for x in rows]); fl = np.array([x['flux'] for x in rows])
            md.append('Median TRAIL_FRAC (contamination / FLUX_AUTO) of the objects with > 3 sigma_ap contamination: %.3f; 90th percentile %.3f.' % (np.median(fr[sel]), np.percentile(fr[sel], 90)))
    if exp == 'stack':
        md.append('\n### 5-frame stacks, different trail per frame (w=6, frame noise 4x the HUDF pixel noise), objects within 5 px of a trail\n')
        names = ['mean_plain', 'median_plain', 'sigclip_plain', 'mean_masked', 'median_masked', 'sigclip_masked']
        md.append('Bias = F(stack with trails) - F(stack of the same frames without trails), in units of the stack aperture noise.\n')
        md.append('| trail amp [frame sigma] | n objs | trails found | ' + ' | '.join(n + ' med / RMS' for n in names) + ' |')
        md.append('|---|---|---|' + '---|' * len(names))
        for k in sorted(set(r['spec']['k'] for r in R)):
            rr = [r for r in R if r['spec']['k'] == k]; rows = [x for r in rr for x in r['rows']]
            cells = []
            for n in names:
                v = np.array([x[n] / x['sst'] for x in rows], float)
                cells.append('%+.2f / %.2f' % (np.nanmedian(v), rms(v)) if np.isfinite(v).any() else 'n/a')
            md.append('| %g | %d | %d/%d | %s |' % (k, len(rows), sum(r['found'] for r in rr), sum(r['nf'] for r in rr), ' | '.join(cells)))
        md.append('\nMean residual on the trail pixels (stack - reference, units of stack sigma, mean over trail pixels): ' + '; '.join(
            '%s %.2f' % (n, np.nanmean([r['pix'][n][0] for r in R if r['pix'][n][0] is not None])) for n in names))
    errs = [r for r in R if False]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp', default='detect,fp,cross,phot,stack')
    ap.add_argument('--out', required=True)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from multiprocessing import Pool
    allres, md = {}, ['# ogfkit.trails validation (generated by trails_validate.py%s)' % (' --quick' if a.quick else '')]
    for exp in a.exp.split(','):
        specs = campaign(exp, a.quick)
        t0 = time.time()
        with Pool(a.workers) as p:
            R = p.map(dispatch, specs, chunksize=1)
        errs = [r for r in R if r and 'error' in r]
        print('%s: %d runs in %.0f s, %d errors' % (exp, len(R), time.time() - t0, len(errs)), flush=True)
        if errs:
            print(errs[0]['error'])
        allres[exp] = R
        summarize(exp, R, md)
        json.dump(allres, open(os.path.join(a.out, 'results.json'), 'w'), default=float)
        open(os.path.join(a.out, 'trails_validation.md'), 'w').write('\n'.join(md) + '\n')
    print('\n'.join(md))


if __name__ == '__main__':
    main()
