#!/usr/bin/env python3
"""Validation of ogfkit/stacking.py against synthetic truth and real data.

    stacking_validate.py faint|align|halo|mask|real [--out report.json] [--real-img FITS --real-cat TSV]

faint  stacks of sources far below the detection limit: unbiased mean flux, bootstrap error calibration, S/N ~ sqrt(N)
align  sub-pixel vs integer alignment: width of the stack of Gaussians at random sub-pixel positions
halo   low surface brightness halo (power law, 0.5 sigma/pixel at 10 px) around bright cores in a crowded noisy field
mask   random bad regions / NaN: recovered flux unbiased
real   null-stack calibration on real HUDF/M51 pixels and a faint-source stack across HUDF bands
"""
import argparse
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import stacking as st, tsvio  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402


def gauss_stamp(img, x, y, F, sg, r=None):
    r = r or int(math.ceil(6 * sg))
    x0, y0 = int(round(x)) - r, int(round(y)) - r
    if x0 < 0 or y0 < 0 or x0 + 2 * r + 1 > img.shape[1] or y0 + 2 * r + 1 > img.shape[0]:
        return
    yy, xx = np.mgrid[y0:y0 + 2 * r + 1, x0:x0 + 2 * r + 1]
    # integrate over the pixel (centre + 3x3 subsampling) for an unbiased pixel-integrated profile
    acc = 0
    for dy in (-1 / 3, 0, 1 / 3):
        for dx in (-1 / 3, 0, 1 / 3):
            acc = acc + np.exp(-((xx + dx - x) ** 2 + (yy + dy - y) ** 2) / (2 * sg * sg))
    img[y0:y0 + 2 * r + 1, x0:x0 + 2 * r + 1] += F * acc / 9.0 / (2 * math.pi * sg * sg)


def faint(args):
    out = {}
    sh = (1400, 1400)
    sg, ap = 2.0, 8.0
    F0 = 3.0                          # per-source aperture S/N ~ F / (sqrt(pi ap^2) sigma) = 0.21
    for corr in (0.0, 1.0):
        for dist in ('const', 'lognormal'):
            for N in (100, 400, 1600):
                R = 30
                est, err, med = [], [], []
                for k in range(R):
                    rng = np.random.default_rng(1000 + k)
                    img = rng.normal(0, 1, sh)
                    if corr:
                        img = ndi.gaussian_filter(img, corr) * (2 * math.sqrt(math.pi) * corr)   # correlated noise, sigma_pix stays 1
                    xs = rng.uniform(30, sh[1] - 30, N)
                    ys = rng.uniform(30, sh[0] - 30, N)
                    F = np.full(N, F0) if dist == 'const' else rng.lognormal(math.log(F0) - 0.5 * 0.8 ** 2, 0.8, N)
                    for x, y, f in zip(xs, ys, F):
                        gauss_stamp(img, x, y, f, sg)
                    cfg = st.StackConfig(half=16, method='mean', n_boot=60, ap_r=ap, seed=k)
                    r = st.run_stack(img, None, xs, ys, cfg)
                    est.append((r['aper'], r['aper_err'], F.mean()))
                    cm = st.StackConfig(half=16, method='median', n_boot=1, ap_r=ap)
                    c = st.make_cutouts(img, None, xs, ys, cm)
                    med.append(st.aperture(st.combine(c['cube'][c['ok']], 'median'), ap))
                e = np.array(est)
                pull = (e[:, 0] - e[:, 2]) / e[:, 1]
                key = 'noise_corr%g_%s_N%d' % (corr, dist, N)
                out[key] = dict(truth=float(e[:, 2].mean()), mean_est=float(e[:, 0].mean()), bias_sigma=float((e[:, 0] - e[:, 2]).mean() / (np.std(e[:, 0] - e[:, 2], ddof=1) / math.sqrt(R))),
                                scatter=float(np.std(e[:, 0] - e[:, 2], ddof=1)), boot_err=float(e[:, 1].mean()), pull_std=float(pull.std(ddof=1)),
                                snr=float(e[:, 0].mean() / e[:, 1].mean()), median_est=float(np.mean(med)))
                print(key, {k: round(v, 3) for k, v in out[key].items()}, flush=True)
    return out


def align(args):
    out = {}
    rng = np.random.default_rng(5)
    sh = (700, 700)
    sg = 0.8
    N = 300
    img = rng.normal(0, 1, sh) * 0.3
    xs = rng.uniform(30, 670, N)
    ys = rng.uniform(30, 670, N)
    for x, y in zip(xs, ys):
        gauss_stamp(img, x, y, 3000.0, sg)
    for al in ('sub', 'int'):
        cfg = st.StackConfig(half=10, align=al, method='mean', n_boot=1, bkg='annulus')
        r = st.run_stack(img, None, xs, ys, cfg)
        s = np.nan_to_num(r['stack'])
        y, x = np.indices(s.shape) - 10
        w = np.where(np.hypot(x, y) <= 5, s, 0.0)
        sx = math.sqrt((w * x * x).sum() / w.sum())
        out[al] = dict(sigma_x=sx, peak_over_flux=float(s[10, 10] / w.sum()), centroid=[float((w * x).sum() / w.sum()), float((w * y).sum() / w.sum())])
    ideal = np.zeros((21, 21))
    gauss_stamp(ideal, 10, 10, 1.0, sg, r=10)
    y, x = np.indices(ideal.shape) - 10
    w = np.where(np.hypot(x, y) <= 5, ideal, 0.0)
    out['ideal_pixel_integrated'] = dict(sigma_x=math.sqrt((w * x * x).sum() / w.sum()), peak_over_flux=float(ideal[10, 10] / w.sum()))
    print(json.dumps(out, indent=1))
    return out


def halo_field(rng, sh, N, A=0.5, r0=10.0, slope=2.0, nneigh=900, grad=True):
    img = rng.normal(0, 1, sh)
    yy, xx = np.indices(sh)
    sky = 3.0 + 0.002 * (xx - sh[1] / 2) + 0.0015 * (yy - sh[0] / 2) if grad else 0.0
    img += sky
    xs = rng.uniform(75, sh[1] - 75, N)
    ys = rng.uniform(75, sh[0] - 75, N)
    # halo shape cut at 70 px
    rr = np.arange(-70, 71)
    RX, RY = np.meshgrid(rr, rr)
    rad = np.hypot(RX, RY)
    halo = A * (np.maximum(rad, 1.0) / r0) ** (-slope)
    halo[rad > 70] = 0
    for x, y in zip(xs, ys):
        xi, yi = int(round(x)), int(round(y))
        img[yi - 70:yi + 71, xi - 70:xi + 71] += halo      # integer position: halo centre rounded (the core carries the subpixel offset)
        gauss_stamp(img, x, y, 2e4, 1.6)
    nx = rng.uniform(10, sh[1] - 10, nneigh)
    ny = rng.uniform(10, sh[0] - 10, nneigh)
    for x, y in zip(nx, ny):
        gauss_stamp(img, x, y, rng.lognormal(math.log(60), 1.0), rng.uniform(1.5, 3.0))
    cat = np.column_stack([np.concatenate([xs, nx]), np.concatenate([ys, ny])])
    return img, xs, ys, cat, (rr, halo), sky


def halo(args):
    out = {}
    sh = (1800, 1800)
    N = 60
    modes = {
        'annulus_bkg_no_mask': dict(bkg='annulus', bkg_in=0.85),
        'none_nomask': dict(bkg='none'),
        'none_catmask': dict(bkg='none', mask_cat_scale=1.0, mask_cat_min=5.0),
        'none_detmask': dict(bkg='none', mask_detect=3.0, protect_r=4.0),
    }
    R = 12
    edges = st.radial_bins(60, 6.0, 1.3)
    sel = [4, 8, 14, 22, 32, 45]
    acc = {m: [] for m in modes}
    accn = {m: [] for m in modes}
    truth = None
    for k in range(R):
        rng = np.random.default_rng(300 + k)
        img, xs, ys, cat, (rr, h), skym = halo_field(rng, sh, N)
        sub = img - skym                                   # exact sky: tests the stacking, not the sky subtraction (a mesh sky eats the halo)
        neigh = np.column_stack([cat, np.maximum(5.0, 1.0 * 3.0 * np.full(len(cat), 2.0))])
        for m, kw in modes.items():
            cfg = st.StackConfig(half=60, n_boot=1, method='mean', **kw)
            prof = st.Profiler(121, edges)
            data = img if m.startswith('annulus') else sub
            cut = st.make_cutouts(data, None, xs, ys, cfg, neigh=neigh, self_index=np.arange(N))
            res = st.stack_cutouts(cut['cube'], cut['ok'], cfg, profiler=prof)
            nx_, ny_ = st.null_positions(data, None, 150, cfg, avoid=cat)
            nc = st.make_cutouts(data, None, nx_, ny_, cfg, neigh=neigh)
            nres = st.stack_cutouts(nc['cube'], nc['ok'], cfg, profiler=prof)
            acc[m].append(res['profile'])
            accn[m].append(nres['profile'])
        truth_prof = np.zeros(len(edges) - 1)
        # truth of the halo profile (+ core PSF) in the bins: average of the model over the pixels of each bin
        yy, xx = np.indices((121, 121)) - 60
        r = np.hypot(xx, yy)
        mod = 0.5 * (np.maximum(r, 1.0) / 10.0) ** -2.0
        ib = np.digitize(r.ravel(), edges) - 1
        for b in range(len(truth_prof)):
            truth_prof[b] = mod.ravel()[ib == b].mean() if (ib == b).any() else np.nan
        rmid = prof.rmid
    for m in modes:
        P = np.array(acc[m])
        Pn = np.array(accn[m])
        mean = P.mean(0)
        sd = P.std(0, ddof=1) / math.sqrt(R)
        meann = Pn.mean(0)
        sdn = Pn.std(0, ddof=1) / math.sqrt(R)
        corr = P - Pn
        cm = corr.mean(0)
        csd = corr.std(0, ddof=1) / math.sqrt(R)
        rows = []
        for b in range(len(rmid)):
            if rmid[b] < 8:
                continue
            rows.append([float(rmid[b]), float(truth_prof[b]), float(mean[b]), float(meann[b]), float(cm[b]), float(csd[b])])
        out[m] = rows
    out['columns'] = ['r', 'truth', 'raw_mean', 'null_mean', 'corrected_mean', 'corrected_err_of_mean']
    for m in modes:
        print(m)
        for row in out[m][::3]:
            print('  r=%5.1f truth=%.4f raw=%.4f null=%.4f corr=%.4f+-%.4f' % tuple(row))
    return out


def mask(args):
    out = {}
    sh = (1200, 1200)
    N = 600
    res = {}
    for frac in (0.0, 0.2, 0.5):
        est = []
        for k in range(12):
            rng = np.random.default_rng(900 + k)
            img = rng.normal(0, 1, sh)
            xs = rng.uniform(40, sh[1] - 40, N)
            ys = rng.uniform(40, sh[0] - 40, N)
            for x, y in zip(xs, ys):
                gauss_stamp(img, x, y, 6.0, 2.0)
            bad = np.zeros(sh, bool)
            if frac:
                # random blobs of 6 px
                lab = rng.random(sh) < 0.01
                bad = ndi.binary_dilation(lab, iterations=3)
                # scale to the requested masked fraction by thresholding a smooth random field
                f = ndi.gaussian_filter(rng.normal(size=sh), 6)
                bad = f > np.quantile(f, 1 - frac)
                img = np.where(bad, 1e6, img)                            # junk in the masked pixels
            cfg = st.StackConfig(half=16, method='mean', n_boot=1, ap_r=8, min_valid=0.3, core_r=0.0)
            r = st.run_stack(img, bad if frac else None, xs, ys, cfg)
            est.append((r['aper'], r['n']))
        e = np.array(est)
        res['masked_%.1f' % frac] = dict(mean_est=float(e[:, 0].mean()), sem=float(e[:, 0].std(ddof=1) / math.sqrt(len(e))), mean_n=float(e[:, 1].mean()), truth=6.0)
        print('masked', frac, res['masked_%.1f' % frac], flush=True)
    return res


def real(args):
    from ogfkit import imageio
    out = {}
    # (1) null stack calibration: split blank-sky positions into sets; z = aperture / bootstrap error
    for name, img_path, cat_path, half in (('hudf_f160w', args.real_img, args.real_cat, 24), ('m51', args.m51_img, args.m51_cat, 16)):
        if not (img_path and os.path.exists(img_path)):
            continue
        data, _ = imageio.load_image(img_path)
        cols, rows = tsvio.read_catalog(cat_path)
        cat = np.array([[tsvio.fnum(r['X_IMAGE']) - 1, tsvio.fnum(r['Y_IMAGE']) - 1] for r in rows])
        cfg = st.StackConfig(half=half, method='mean', n_boot=100, ap_r=6.0, bkg='annulus', mask_cat_scale=1.0, mask_cat_min=4.0, min_valid=0.5)
        neigh = np.column_stack([cat, [max(4.0, 1.0 * tsvio.fnum(r.get('A_IMAGE'), 2) * tsvio.fnum(r.get('KRON_RADIUS'), 3)) for r in rows]])
        nx_, ny_ = st.null_positions(data, None, 3000 if name != 'm51' else 1500, cfg, avoid=cat)
        call = st.make_cutouts(data, None, nx_, ny_, cfg, neigh=neigh)
        good = np.where(call['ok'])[0]
        print(name, 'null cutouts ok', len(good), 'of', len(nx_), flush=True)
        for N in (30, 100, 300):
            zs, ap, ee = [], [], []
            nsets = min(len(good) // N, 24)
            for s_ in range(nsets):
                sl = good[s_ * N:(s_ + 1) * N]
                ok = np.zeros(len(call['ok']), bool)
                ok[sl] = True
                r = st.stack_cutouts(call['cube'], ok, cfg)
                zs.append(r['aper'] / r['aper_err'])
                ap.append(r['aper'])
                ee.append(r['aper_err'])
            if not zs:
                continue
            zs = np.array(zs)
            out['%s_null_N%d' % (name, N)] = dict(nsets=len(zs), z_mean=float(zs.mean()), z_std=float(zs.std(ddof=1)), scatter_ap=float(np.std(ap, ddof=1)), mean_boot_err=float(np.mean(ee)))
            print(name, N, out['%s_null_N%d' % (name, N)], flush=True)
    # (2) faint-source stack across the HUDF bands: positions of F160W catalog sources in the faintest bin
    if args.real_cat and os.path.exists(args.real_cat) and args.bands:
        cols, rows = tsvio.read_catalog(args.real_cat)
        mag = np.array([tsvio.fnum(r['MAG_AUTO']) for r in rows])
        cat = np.array([[tsvio.fnum(r['X_IMAGE']) - 1, tsvio.fnum(r['Y_IMAGE']) - 1] for r in rows])
        neigh = np.column_stack([cat, [max(4.0, 1.5 * tsvio.fnum(r.get('A_IMAGE'), 2) * tsvio.fnum(r.get('KRON_RADIUS'), 3)) for r in rows]])
        for lo, hi in ((22, 25), (25, 26.5), (26.5, 28)):
            idx = np.where((mag >= lo) & (mag < hi))[0]
            for bname, bpath in args.bands:
                data, _ = imageio.load_image(bpath)
                cfg = st.StackConfig(half=20, method='mean', n_boot=150, ap_r=6.0, mask_cat_scale=1.0, mask_cat_min=4.0, min_valid=0.5)
                r = st.run_stack(data, None, cat[idx, 0], cat[idx, 1], cfg, null=300, avoid=cat, neigh=neigh, self_index=idx)
                nl = r['null']
                ap, ae = r['aper'] - nl['aper'], math.hypot(r['aper_err'], nl['aper_err'])
                out['bin%g-%g_%s' % (lo, hi, bname)] = dict(n_sel=int(len(idx)), n=int(r['n']), ap_mean_per_source=float(ap), err=float(ae), snr=float(ap / ae), null_ap=float(nl['aper']), null_err=float(nl['aper_err']))
                print(lo, hi, bname, out['bin%g-%g_%s' % (lo, hi, bname)], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['faint', 'align', 'halo', 'mask', 'real'])
    ap.add_argument('--out')
    ap.add_argument('--real-img')
    ap.add_argument('--real-cat')
    ap.add_argument('--m51-img')
    ap.add_argument('--m51-cat')
    ap.add_argument('--band', action='append', default=[], help='NAME=FITS (real mode)')
    a = ap.parse_args()
    a.bands = [b.split('=', 1) for b in a.band]
    res = globals()[a.mode](a)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
