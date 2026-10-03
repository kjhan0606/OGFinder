#!/usr/bin/env python3
"""Real-galaxy cross-check of multifit against GALFIT (same data, PSF, sigma, mask, start values, constraints).

    python galfit_real_compare.py prep  WORK [--n 30] [--hudf-dir /workspace/fits] [--m51 FILE]   # catalog, stacked PSF per band, cutouts, mask, sigma, start feedmes
    python galfit_real_compare.py run   WORK [--workers 8] [--galfit BIN] [--ld DIR]              # GALFIT and multifit on every cutout
    python galfit_real_compare.py report WORK OUT.json                                             # statistics (parameters, chi2 with a neutral arbiter, failures)

Data: HUDF F105W / F125W / F160W mosaics (same pixel grid, 0.06"/px; AB zero point from PHOTFLAM / PHOTPLAM) and the M51 image (zero point 25, arbitrary units).
Objects are chosen on F160W (extended: half-light radius > 4 px, 19.5 < mag < 23.5) and fitted in every band at the same position.
Every cutout: background-subtracted with sep (64 px mesh), constant sigma = local sep rms of the cutout (identical for both programs), PSF = median stack of isolated stars of that band
(sub-pixel recentred), mask = other detected objects (segmentation map dilated by 2 px; the target's own segment is kept), start values from the sep moments.
Models: `sersic` (free n) on every cutout; `devexp` (de Vaucouleurs + exponential disc with tied centre) on the same cutouts.  Both programs fit the same feedme, and each result is
re-rendered with GALFIT (P=1) and its chi2 evaluated in numpy with the same sigma and mask (neutral arbiter).  There is no truth on real galaxies: only the agreement of the two programs
and the chi2 reached can be assessed (not which one is closer to the physical parameters).
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path[:0] = [ROOT, HERE]
from astropy.io import fits  # noqa: E402
from ogfkit import galfitio as GI  # noqa: E402
import galfit_compare as GC  # noqa: E402

MULTIFIT = os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py')
BANDS = [('f160w', 'hudf_f160w.fits'), ('f125w', 'hudf_f125w.fits'), ('f105w', 'hudf_f105w.fits')]


def ab_zp(h):
    try:
        return -2.5 * math.log10(h['PHOTFLAM']) - 21.10 - 5 * math.log10(h['PHOTPLAM']) + 18.6921
    except Exception:
        return 25.0


def load(path):
    with fits.open(path) as hl:
        for x in hl:
            if x.data is not None and x.data.ndim == 2:
                return np.array(x.data, float), x.header
    raise ValueError(path)


def detect(d):
    import sep
    d = np.ascontiguousarray(np.nan_to_num(d))
    bk = sep.Background(d, bw=64, bh=64)
    ds = d - bk.back()
    rms = bk.rms()
    o, seg = sep.extract(ds, 2.0, err=rms, minarea=12, deblend_nthresh=32, deblend_cont=0.005, segmentation_map=True)
    fr, _ = sep.flux_radius(ds, o['x'], o['y'], 6 * np.maximum(o['a'], 1.0), 0.5)
    return ds, rms, o, seg, fr


def stack_psf(ds, rms, o, fr, size=31, nmax=40, snr_min=40):
    import sep
    from scipy.ndimage import shift as ndshift
    ny, nx = ds.shape
    flux = o['flux']
    good = (fr > 0) & (flux > 0)
    mag = np.where(good, -2.5 * np.log10(np.where(good, flux, 1)), 99)
    # stellar locus: smallest half-light radii among objects with decent S/N
    snr = flux / np.sqrt(np.maximum(o['npix'], 1)) / np.median(rms)
    cand = np.where(good & (snr > snr_min))[0]
    med = np.median(fr[cand])
    lo = np.percentile(fr[cand], 10)
    stars = [i for i in cand if abs(fr[i] - lo) < 0.12 * lo and o['a'][i] / max(o['b'][i], 0.1) < 1.25 and 40 < o["x"][i] < nx - 40 and 40 < o["y"][i] < ny - 40]
    stars = sorted(stars, key=lambda i: -flux[i])
    # isolation: no other detection within 20 px
    xy = np.c_[o['x'], o['y']]
    keep = []
    for i in stars:
        dd = np.hypot(xy[:, 0] - o['x'][i], xy[:, 1] - o['y'][i])
        dd[i] = 1e9
        if dd.min() > 20:
            keep.append(i)
    keep = keep[:nmax]
    h = size // 2 + 4
    cut = []
    for i in keep:
        x, y = o['x'][i], o['y'][i]
        xw, yw, _ = sep.winpos(ds, [x], [y], [fr[i] * 1.2 / 1.18], subpix=5)
        x, y = float(xw[0]), float(yw[0])
        ix, iy = int(round(x)), int(round(y))
        st = ds[iy - h:iy + h + 1, ix - h:ix + h + 1].copy()
        if st.shape != (2 * h + 1, 2 * h + 1):
            continue
        st = ndshift(st, (-(y - iy), -(x - ix)), order=3, mode='nearest')
        yy, xx = np.mgrid[:2 * h + 1, :2 * h + 1] - h
        rr = np.hypot(xx, yy)
        st = st - np.median(st[rr > h - 1])
        s = st[rr < 6].sum()
        if s <= 0:
            continue
        cut.append(st / s)
    if len(cut) < 4:
        raise RuntimeError('too few stars for the PSF (%d)' % len(cut))
    p = np.median(np.array(cut), axis=0)[4:-4, 4:-4]
    p = denoise_wings(p, core=4.0, rmax=min(size // 2, 12))
    return p.astype(np.float32), len(cut), float(np.median(fr[keep]) if keep else med)


def denoise_wings(p, core=4.0, rmax=12):
    """Stacked-star PSF: keep the measured pixels inside `core` px, replace the wings by the azimuthal median profile (the wing pixels of a stack of ~30 faint stars are pure noise,
    and a PSF with a negative sum of -0.5 spoils the fits), taper to zero at rmax, clip negatives, unit sum."""
    n = p.shape[0]
    yy, xx = np.mgrid[:n, :n] - n // 2
    rr = np.hypot(xx, yy)
    prof = np.array([np.median(p[(rr >= k - 0.5) & (rr < k + 0.5)]) if ((rr >= k - 0.5) & (rr < k + 0.5)).any() else 0.0 for k in range(n)])
    for k in range(2, len(prof)):                                       # monotonic decline outside the core
        prof[k] = min(prof[k], prof[k - 1])
    wing = np.interp(rr, np.arange(n), np.clip(prof, 0, None))
    w = np.clip((rr - core) / 2.0, 0.0, 1.0)                           # 0 inside the core, 1 beyond core + 2 px
    q = (1 - w) * p + w * wing
    q = q * np.clip((rmax - rr) / 3.0, 0.0, 1.0)                       # taper over the last 3 px
    q = np.clip(q, 0, None)
    return q / q.sum()


def sigma_clip_rms(a, m):
    v = a[~m]
    for _ in range(4):
        s = 1.4826 * np.median(np.abs(v - np.median(v)))
        v = v[np.abs(v - np.median(v)) < 3 * s]
    return float(np.std(v))


def prep(a):
    os.makedirs(a.work, exist_ok=True)
    sets = [(b, os.path.join(a.hudf_dir, f)) for b, f in BANDS if os.path.exists(os.path.join(a.hudf_dir, f))]
    data = {}
    for b, p in sets:
        d, h = load(p)
        data[b] = (d, h, ab_zp(h))
        print(b, d.shape, 'zp %.3f' % data[b][2], flush=True)
    det = {b: detect(data[b][0]) for b in data}
    ds, rms, o, seg, fr = det['f160w']
    zp = data['f160w'][2]
    mag = zp - 2.5 * np.log10(np.maximum(o['flux'], 1e-9))
    h = np.clip(np.round(7 * fr), 36, 80)
    ny, nx = ds.shape
    xy = np.c_[o['x'], o['y']]
    sel = [i for i in range(len(o)) if 19.5 < mag[i] < 23.5 and fr[i] > 4.0 and o['b'][i] / o['a'][i] > 0.15 and h[i] + 3 < o['x'][i] < nx - h[i] - 3 and h[i] + 3 < o['y'][i] < ny - h[i] - 3]
    print(len(o), 'detected,', len(sel), 'extended candidates', flush=True)
    rng = np.random.default_rng(7)
    sel = sorted(rng.permutation(sel)[:a.n], key=lambda i: mag[i])
    items = []
    for b in data:
        d, hd, zpb = data[b]
        dsb, rmsb, ob, segb, frb = det[b]
        psf, nst, frs = stack_psf(dsb, rmsb, ob, frb)
        pdir = os.path.join(a.work, b)
        os.makedirs(pdir, exist_ok=True)
        fits.PrimaryHDU(psf).writeto(os.path.join(pdir, 'psf.fits'), overwrite=True)
        print(b, 'PSF from %d stars, star half-light radius %.2f px, psf peak %.4f' % (nst, frs, psf.max()), flush=True)
        for k, i in enumerate(sel):
            items.append(dict(band=b, idx=int(i), k=k, x=float(o['x'][i]), y=float(o['y'][i]), zp=zpb, scale=0.06, psf_stars=nst))
    bandsel = {b: (det[b], data[b]) for b in data}
    if a.m51 and os.path.exists(a.m51):
        d, hd = load(a.m51)
        data['m51'] = (d, hd, 25.0)
        det['m51'] = detect(d)
        dsb, rmsb, ob, segb, frb = det['m51']
        pdir = os.path.join(a.work, 'm51')
        os.makedirs(pdir, exist_ok=True)
        try:
            psf, nst, frs = stack_psf(dsb, rmsb, ob, frb, size=15, nmax=30, snr_min=10)
            fits.PrimaryHDU(psf).writeto(os.path.join(pdir, 'psf.fits'), overwrite=True)
            print('m51 PSF from %d stars, star r_half %.2f px' % (nst, frs), flush=True)
            magm = 25.0 - 2.5 * np.log10(np.maximum(ob['flux'], 1e-9))
            hm = np.clip(np.round(7 * frb), 24, 60)
            nym, nxm = dsb.shape
            selm = [i for i in range(len(ob)) if frb[i] > 1.8 * frs and ob['flux'][i] > 0 and magm[i] < 19 and ob['b'][i] / ob['a'][i] > 0.15 and hm[i] + 3 < ob['x'][i] < nxm - hm[i] - 3 and hm[i] + 3 < ob['y'][i] < nym - hm[i] - 3]
            selm = sorted(selm, key=lambda i: magm[i])[:a.n_m51]
            for k, i in enumerate(selm):
                items.append(dict(band='m51', idx=int(i), k=k, x=float(ob['x'][i]), y=float(ob['y'][i]), zp=25.0, scale=1.8, psf_stars=nst))
            print('m51: %d extended objects' % len(selm), flush=True)
        except RuntimeError as e:
            print('m51 skipped:', e)
    mf = []
    for it in items:
        b = it['band']
        dsb, rmsb, ob, segb, frb = det[b]
        i = it['idx']
        if b == 'm51':
            hh = int(np.clip(round(7 * frb[i]), 24, 60))
            idx_obj = i
            ob_ = ob
        else:
            # match the F160W object in this band's catalog by position (same pixel grid)
            dd = np.hypot(ob['x'] - it['x'], ob['y'] - it['y'])
            idx_obj = int(np.argmin(dd))
            if dd[idx_obj] > 3:
                idx_obj = None
            ob_ = ob
            hh = int(h[i])
        x, y = it['x'], it['y']
        ix, iy = int(round(x)), int(round(y))
        cut = dsb[iy - hh:iy + hh + 1, ix - hh:ix + hh + 1].copy()
        sg = segb[iy - hh:iy + hh + 1, ix - hh:ix + hh + 1]
        own = sg[hh, hh] if sg[hh, hh] > 0 else (segb == (idx_obj + 1))[iy - hh:iy + hh + 1, ix - hh:ix + hh + 1].astype(int) * 0
        lab = (idx_obj + 1) if idx_obj is not None else 0
        from scipy.ndimage import binary_dilation
        others = (sg > 0) & (sg != lab)
        mask = binary_dilation(others, iterations=2) & ~((sg == lab) & (lab > 0))
        mask |= ~np.isfinite(cut)
        if (~mask).sum() < 0.6 * mask.size:
            continue
        sig = sigma_clip_rms(cut, mask | (sg > 0))
        if idx_obj is None or ob_['flux'][idx_obj] <= 0:
            continue
        zpb = it['zp']
        d0 = os.path.join(a.work, b, 'o%03d' % it['k'])
        os.makedirs(d0, exist_ok=True)
        fits.PrimaryHDU(cut.astype(np.float32)).writeto(os.path.join(d0, 'data.fits'), overwrite=True)
        fits.PrimaryHDU(np.full(cut.shape, sig, np.float32)).writeto(os.path.join(d0, 'sigma.fits'), overwrite=True)
        fits.PrimaryHDU(mask.astype(np.int16)).writeto(os.path.join(d0, 'mask.fits'), overwrite=True)
        shutil.copy(os.path.join(a.work, b, 'psf.fits'), os.path.join(d0, 'psf.fits'))
        m0 = zpb - 2.5 * math.log10(ob_['flux'][idx_obj])
        re0 = float(max(frb[idx_obj], 1.5))
        q0 = float(np.clip(ob_['b'][idx_obj] / ob_['a'][idx_obj], 0.2, 1.0))
        th = math.degrees(ob_['theta'][idx_obj])
        pa0 = ((th - 90 + 90) % 180) - 90
        cx, cy = x - ix + hh + 1, y - iy + hh + 1
        sky0 = float(np.median(cut[~mask & (sg == 0)])) if (~mask & (sg == 0)).sum() > 50 else 0.0
        n = cut.shape[0]
        for model in ('sersic', 'devexp'):
            if model == 'sersic':
                objs = [dict(type='sersic', x=cx, y=cy, mag=m0, re=re0, n=2.0, q=q0, pa=pa0)]
                cons = '1 n 0.3 to 8\n1 re 0.3 to %d\n' % (3 * hh)
            else:
                objs = [dict(type='devauc', x=cx, y=cy, mag=m0 + 0.87, re=max(0.5 * re0, 0.8), q=min(q0 + 0.2, 1.0), pa=pa0),
                        dict(type='expdisk', x=cx, y=cy, mag=m0 + 0.87, re=0.7 * re0, q=q0, pa=pa0)]
                cons = '1_2 x offset\n1_2 y offset\n1 re 0.3 to %d\n2 re 0.3 to %d\n' % (3 * hh, 3 * hh)
            md = os.path.join(d0, model)
            os.makedirs(md, exist_ok=True)
            for f in ('data.fits', 'sigma.fits', 'mask.fits', 'psf.fits'):
                shutil.copy(os.path.join(d0, f), md)
            open(os.path.join(md, 'cons.txt'), 'w').write(cons)
            GC.ZP = zpb
            txt = GC.feed_text(objs, sky0, (0.0, 0.0), mode=0, constraints='cons.txt', region=(1, n, 1, n))
            txt = txt.replace('F) none', 'F) mask.fits').replace('I) 61 61', 'I) 41 41')
            open(os.path.join(md, 'start.feedme'), 'w').write(txt)
            mf.append(dict(band=b, k=it['k'], model=model, dir=os.path.relpath(md, a.work), zp=zpb, sigma=sig, size=n, mag0=m0, re0=re0, x0=x, y0=y, cat_index=int(idx_obj), psf_stars=it['psf_stars'], nmask=int(mask.sum())))
    json.dump(mf, open(os.path.join(a.work, 'manifest.json'), 'w'))
    print('prepared', len(mf), 'fits in', len(set((m['band'], m['k']) for m in mf)), 'cutouts')


def render_chi2(a, md, feed_txt, tag, sigma):
    ed = os.path.join(md, 'e_' + tag)
    os.makedirs(ed, exist_ok=True)
    for f in ('psf.fits', 'data.fits', 'mask.fits'):
        shutil.copy(os.path.join(md, f), ed)
    txt = feed_txt
    for key, val in (('A', 'data.fits'), ('B', 'out.fits'), ('C', 'none'), ('D', 'psf.fits'), ('F', 'none'), ('G', 'none'), ('P', '1')):
        txt = re.sub(r'(?m)^%s\).*$' % key, '%s) %s' % (key, val), txt)
    open(os.path.join(ed, 'e.feedme'), 'w').write(txt)
    GC.run_galfit(a.galfit, a.ld, 'e.feedme', ed, timeout=120)
    o = os.path.join(ed, 'out.fits')
    if not os.path.exists(o):
        return None
    data = fits.getdata(os.path.join(ed, 'data.fits')).astype(float)
    mod = fits.getdata(o)
    mod = mod[0] if mod.ndim == 3 else mod
    bad = fits.getdata(os.path.join(ed, 'mask.fits')) > 0
    r = ((data - mod.astype(float)) / sigma)[~bad]
    return float(np.sum(r ** 2)), int((~bad).sum())


def comps_of(cfg):
    out = []
    for c in cfg['components']:
        d = dict(type=c['galfit_type'], x=c.get('x'), y=c.get('y'), mag=c.get('mag'), re=c.get('re'), n=c.get('n'), q=c.get('q'), pa=c.get('pa', 90) - 90 if c.get('pa') is not None else None)
        if c['galfit_type'] == 'expdisk' and c.get('re') is not None:
            d['re'] = c['re'] / GI.RS_TO_RE
        out.append(d)
    return out


def run_one(task):
    m, a = task
    md = os.path.join(a.work, m['dir'])
    res = dict(m)
    # GALFIT
    gwd = os.path.join(md, 'g')
    shutil.rmtree(gwd, ignore_errors=True)
    shutil.copytree(md, gwd, ignore=shutil.ignore_patterns('g', 'm', 'e_*'))
    try:
        t, r = GC.run_galfit(a.galfit, a.ld, 'start.feedme', gwd, timeout=a.timeout)
        res['galfit_time'] = t
        gr = GC.galfit_result(gwd)
        if gr is None:
            res['galfit'] = None
            res['galfit_err'] = (r.stdout or '')[-160:]
        else:
            cfg, chi, _ = gr
            fs = sorted(f for f in os.listdir(gwd) if re.match(r'galfit\.\d+$', f))
            txt = open(os.path.join(gwd, fs[-1])).read()
            res['galfit'] = dict(comps=comps_of(cfg), chi2nu=chi, flagged=len(re.findall(r'\*[^*\n]+\*', txt.split('# Chi')[0] if False else txt)), sky=cfg.get('sky_value_galfit'))
            res['galfit']['feed'] = txt
    except subprocess.TimeoutExpired:
        res['galfit'] = None
        res['galfit_err'] = 'timeout'
        res['galfit_time'] = a.timeout
    # multifit
    mwd = os.path.join(md, 'm')
    shutil.rmtree(mwd, ignore_errors=True)
    shutil.copytree(md, mwd, ignore=shutil.ignore_patterns('g', 'm', 'e_*'))
    t = time.time()
    r = subprocess.run([sys.executable, MULTIFIT, '-', '--config', os.path.join(mwd, 'start.feedme'), '--work', os.path.join(mwd, 'w'), '--export-feedme', os.path.join(mwd, 'w', 'gf'), '--max-nfev', str(a.max_nfev)],
                       capture_output=True, text=True, cwd=mwd, timeout=a.timeout * 3)
    res['multifit_time'] = time.time() - t
    ex = os.path.join(mwd, 'w', 'gf', 'galfit.feedme')
    if r.returncode != 0 or not os.path.exists(ex):
        res['multifit'] = None
        res['multifit_err'] = (r.stderr or r.stdout)[-200:]
    else:
        cfg = GI.parse_feedme(ex, strict=False, base_dir=mwd)
        mm = re.search(r'chi2/dof = ([0-9.]+)', r.stdout)
        res['multifit'] = dict(comps=comps_of(cfg), chi2nu=float(mm.group(1)) if mm else None, sky=cfg.get('sky_value_galfit'), feed=open(ex).read())
    try:
        if res.get('galfit'):
            res['galfit']['chi2'] = render_chi2(a, md, res['galfit']['feed'], 'g', m['sigma'])
        if res.get('multifit'):
            res['multifit']['chi2'] = render_chi2(a, md, res['multifit']['feed'], 'm', m['sigma'])
    except Exception as e:
        res['chi2_error'] = str(e)
    for w in ('galfit', 'multifit'):
        if res.get(w):
            res[w].pop('feed', None)
    return res


def run(a):
    man = json.load(open(os.path.join(a.work, 'manifest.json')))
    t0 = time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        results = list(ex.map(run_one, [(m, a) for m in man]))
    json.dump(dict(elapsed_s=time.time() - t0, results=results), open(os.path.join(a.work, 'results.json'), 'w'), default=lambda o: None)
    print('done', len(results), 'in %.0f s' % (time.time() - t0))


def ogfkit_chi2(md, feed_path, sigma):
    """chi2 of a result feedme rendered with the ogfkit (multifit) renderer, same data / sigma / mask / PSF."""
    from ogfkit import multifit as MF
    cfg = GI.parse_feedme(feed_path, strict=False, base_dir=md)
    data = fits.getdata(os.path.join(md, 'data.fits')).astype(float)
    bad = fits.getdata(os.path.join(md, 'mask.fits')) > 0
    psf = fits.getdata(os.path.join(md, 'psf.fits')).astype(float)
    psf = np.clip(psf, 0, None)
    psf /= psf.sum()
    cn = [MF._normalise(dict(c), cfg['zp']) for c in cfg['components']]
    for c in cn:
        c['x'] -= 1
        c['y'] -= 1
    mod = MF.render_model(cn, data.shape, psf=psf, sky=cfg.get('sky_value_galfit', 0.0))
    return float((((data - mod) / sigma)[~bad] ** 2).sum()), int((~bad).sum())


def stats(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return None
    return dict(n=int(len(v)), median=float(np.median(v)), nmad=float(1.4826 * np.median(np.abs(v - np.median(v)))), p90_abs=float(np.percentile(np.abs(v), 90)), max_abs=float(np.max(np.abs(v))))


def report(a):
    R = json.load(open(os.path.join(a.work, 'results.json')))
    res = R['results']
    out = dict(elapsed_s=R['elapsed_s'], n_fits=len(res), models={})
    for r in res:                                                                   # cross-evaluation of both solutions in the ogfkit renderer
        md = os.path.join(a.work, r['dir'])
        try:
            gf = sorted(f for f in os.listdir(os.path.join(md, 'g')) if re.match(r'galfit\.\d+$', f))
            r['x_gal_in_ogf'] = ogfkit_chi2(md, os.path.join(md, 'g', gf[-1]), r['sigma'])[0] if r.get('galfit') and gf else None
            r['x_mf_in_ogf'] = ogfkit_chi2(md, os.path.join(md, 'm', 'w', 'gf', 'galfit.feedme'), r['sigma'])[0] if r.get('multifit') else None
        except Exception as e:
            r['x_error'] = str(e)[:80]
    for model in ('sersic', 'devexp'):
        rs = [r for r in res if r['model'] == model]
        both = [r for r in rs if r.get('galfit') and r.get('multifit')]
        s = dict(n=len(rs), galfit_failed=sum(1 for r in rs if not r.get('galfit')), multifit_failed=sum(1 for r in rs if not r.get('multifit')), both=len(both),
                 galfit_flagged_params=sum(1 for r in both if r['galfit'].get('flagged')), galfit_time_median=float(np.median([r['galfit_time'] for r in rs])), multifit_time_median=float(np.median([r['multifit_time'] for r in rs])))
        s['galfit_failures'] = [dict(band=r['band'], k=r['k'], err=r.get('galfit_err')) for r in rs if not r.get('galfit')][:10]
        s['multifit_failures'] = [dict(band=r['band'], k=r['k'], err=(r.get('multifit_err') or '')[-120:]) for r in rs if not r.get('multifit')][:10]
        acc = {}
        for r in both:
            for ci, (g, m) in enumerate(zip(r['galfit']['comps'], r['multifit']['comps'])):
                d = dict(dx=m['x'] - g['x'], dy=m['y'] - g['y'], dmag=m['mag'] - g['mag'])
                if g['type'] != 'psf':
                    d['re_rel'] = m['re'] / g['re'] - 1
                    d['dq'] = m['q'] - g['q']
                    if g['q'] < 0.85 and m['q'] < 0.85:
                        d['dpa'] = GC.wrap_pa(m['pa'] - g['pa'])
                    if g['type'] == 'sersic':
                        d['n_rel'] = m['n'] / g['n'] - 1
                for p, v in d.items():
                    acc.setdefault('c%d.%s' % (ci, p), []).append(v)
        s['multifit_minus_galfit'] = {k: stats(v) for k, v in acc.items()}
        c = [(r['galfit']['chi2'], r['multifit']['chi2']) for r in both if r['galfit'].get('chi2') and r['multifit'].get('chi2')]
        if c:
            gg = np.array([x[0][0] for x in c]); mm = np.array([x[1][0] for x in c]); npx = np.array([x[0][1] for x in c])
            d = mm - gg
            s['chi2_arbiter'] = dict(n=len(c), chi2red_galfit_median=float(np.median(gg / npx)), chi2red_multifit_median=float(np.median(mm / npx)), diff_median=float(np.median(d)),
                                     rel_diff_median=float(np.median(d / gg)), within_0p5pct=float(np.mean(np.abs(d / gg) <= 0.005)), within_2pct=float(np.mean(np.abs(d / gg) <= 0.02)),
                                     multifit_lower_by_gt_0p5pct=float(np.mean(d / gg < -0.005)), galfit_lower_by_gt_0p5pct=float(np.mean(d / gg > 0.005)), worst_rel=float(np.max(d / gg)), best_rel=float(np.min(d / gg)))
        xs = [r for r in both if r.get('x_gal_in_ogf') and r.get('x_mf_in_ogf') and r['galfit'].get('chi2') and r['multifit'].get('chi2')]
        if xs:
            gm = np.array([r['x_gal_in_ogf'] for r in xs]); mm = np.array([r['x_mf_in_ogf'] for r in xs])
            gg = np.array([r['galfit']['chi2'][0] for r in xs]); mg = np.array([r['multifit']['chi2'][0] for r in xs])
            s['cross_chi2'] = dict(n=len(xs), renderer_agreement_galfit_solution_median_rel=float(np.median(gm / gg - 1)), renderer_agreement_multifit_solution_median_rel=float(np.median(mm / mg - 1)),
                                   renderer_agreement_p90_abs=float(np.percentile(np.abs(np.r_[gm / gg - 1, mm / mg - 1]), 90)),
                                   in_ogfkit_multifit_lower_frac=float(np.mean(mm < gm * 0.995)), in_ogfkit_galfit_lower_frac=float(np.mean(gm < mm * 0.995)), in_ogfkit_within_0p5pct_frac=float(np.mean(np.abs(mm / gm - 1) <= 0.005)),
                                   in_galfit_multifit_lower_frac=float(np.mean(mg < gg * 0.995)), in_galfit_galfit_lower_frac=float(np.mean(gg < mg * 0.995)))
        agree = [r for r in both if abs(r['multifit']['comps'][0]['mag'] - r['galfit']['comps'][0]['mag']) < 0.05]
        s['mag_within_0p05_frac'] = len(agree) / max(len(both), 1)
        s['mag_within_0p1_frac'] = float(np.mean([abs(r['multifit']['comps'][0]['mag'] - r['galfit']['comps'][0]['mag']) < 0.1 for r in both])) if both else None
        # per band
        s['by_band'] = {}
        for b in sorted({r['band'] for r in rs}):
            bb = [r for r in both if r['band'] == b]
            if bb:
                s['by_band'][b] = dict(n=len(bb), dmag=stats([r['multifit']['comps'][0]['mag'] - r['galfit']['comps'][0]['mag'] for r in bb]))
        out['models'][model] = s
    json.dump(out, open(a.out, 'w'), indent=1)
    for model, s in out['models'].items():
        print(model, 'n', s['n'], 'galfit fail', s['galfit_failed'], 'multifit fail', s['multifit_failed'], 'both', s['both'], 'mag<0.05: %.2f' % s['mag_within_0p05_frac'])
        ca = s.get('chi2_arbiter')
        if ca:
            print('  chi2 rel diff median %.4f, within 0.5%% %.2f, 2%% %.2f, multifit lower %.2f, galfit lower %.2f' % (ca['rel_diff_median'], ca['within_0p5pct'], ca['within_2pct'], ca['multifit_lower_by_gt_0p5pct'], ca['galfit_lower_by_gt_0p5pct']))
        for k in ('c0.dmag', 'c0.re_rel', 'c0.n_rel', 'c0.dq', 'c0.dpa'):
            v = s['multifit_minus_galfit'].get(k)
            if v:
                print('  %s median %+.4f nmad %.4f p90 %.4f max %.3f' % (k, v['median'], v['nmad'], v['p90_abs'], v['max_abs']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prep', 'run', 'report'])
    ap.add_argument('work')
    ap.add_argument('out', nargs='?')
    ap.add_argument('--n', type=int, default=30)
    ap.add_argument('--n-m51', type=int, default=12, dest='n_m51')
    ap.add_argument('--hudf-dir', default='/workspace/fits')
    ap.add_argument('--m51', default='/workspace/fits/m51.fits')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--galfit', default=os.environ.get('GALFIT_BIN') or shutil.which('galfit') or '')
    ap.add_argument('--ld', default=os.environ.get('GALFIT_LD') or '')
    ap.add_argument('--timeout', type=int, default=300)
    ap.add_argument('--max-nfev', type=int, default=400, dest='max_nfev')
    a = ap.parse_args()
    a.work = os.path.abspath(a.work)
    if a.cmd == 'prep':
        prep(a)
    elif a.cmd == 'run':
        run(a)
    else:
        a.out = a.out or os.path.join(a.work, 'report.json')
        report(a)


if __name__ == '__main__':
    main()
