#!/usr/bin/env python3
"""Validation of the depth / completeness maps (ogfkit/depth.py, plugins/completeness/depthmap.py).

    depth_validate.py synth|compmap|real [--out report.json] [--hudf FITS --m51 FITS]

synth    depth maps of noise images with known truth: white noise (analytic), correlated noise (brute force over 20000 blank apertures),
         spatially varying noise (gradient, four patches), with injected sources and masked regions; weight map input
compmap  completeness by region for an image with four noise patches: per-region lim50 versus the lim50 measured on stand-alone images of the
         same noise level (truth), collapse of the curves in (m - local depth)
real     HUDF F160W (1800^2 crop) and M51: depth, tile check, area vs depth, regional completeness
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'plugins', 'completeness'))
import warnings  # noqa: E402
warnings.filterwarnings('ignore')
from ogfkit import depth as dp, noise as nz  # noqa: E402
import completeness as cp  # noqa: E402
import depthmap as dm  # noqa: E402


def noise_image(shape, sigma_map, corr, rng):
    n = rng.normal(size=shape)
    if corr > 0:
        n = ndi.gaussian_filter(n, corr) * (2 * math.sqrt(math.pi) * corr)           # sigma_pix stays 1
    return n * sigma_map


def brute_sigma(noise, r, n, rng, box=None):
    pos = dp.blank_circle_positions(np.zeros(noise.shape, bool), r, n, rng, box=box)
    return nz.clipped_std(dp.circle_sums(noise, pos, r), 4.0)


def depth_of(img, bad=None, r=3.0, zp=25.0, nsig=5.0, rms=None, bw=64):
    P = dp.prepare(img, bad, bw=bw)
    valid = ~P['bad']
    sp = nz.robust_std(P['sub'][~P['mask']][:2000000])
    rm = P['rms'] if rms is None else dp.rms_from_file(rms[0], valid, P['rms'], rms[1])
    f, _ = dp.circle_factor(P['sub'], P['mask'], rm, r, 3000, 1)
    return dp.mag_limit_map(rm, f, zp, nsig), P, f, rm


def synth(args):
    out = {}
    sh = (1200, 1200)
    r, zp, ps = 3.0, 25.0, 0.1
    yy, xx = np.indices(sh)
    # 1. white noise, analytic
    for sig in (0.5, 2.0):
        img = np.random.default_rng(1).normal(0, sig, sh).astype('f4')
        mag, P, f_, rm_ = depth_of(img, r=r, zp=zp)
        truth = zp - 2.5 * math.log10(5 * brute_sigma(img, r, 20000, np.random.default_rng(3)))          # sep.sum_circle fractional weights: sum(w^2) < N
        out['white_sigma%g' % sig] = dict(truth=truth, median=float(np.median(mag)), diff=float(np.median(mag) - truth), factor=f_)
        # surface brightness
        L = 50
        bl = dp.box_law(P['sub'], P['mask'], P['rms'], [8, 14, 22, 36])
        fb = bl['gamma'] * float(L * L) ** bl['delta']          # power law through the smaller boxes (the L=50 box has few independent samples)
        sb = dp.sb_limit_map(P['rms'], fb, L, ps, zp, 3.0)
        tsb = zp - 2.5 * math.log10(3 * sig * L / (L * ps) ** 2)
        out['white_sigma%g' % sig].update(sb_truth=tsb, sb_median=float(np.median(sb)), sb_diff=float(np.median(sb) - tsb))
    # 2. correlated noise, brute-force truth
    for corr in (1.0, 1.5):
        rng = np.random.default_rng(2)
        nse = noise_image(sh, 1.0, corr, rng).astype('f4')
        s_true = brute_sigma(nse, r, 20000, np.random.default_rng(3))
        mag, P, f_, _ = depth_of(nse, r=r, zp=zp)
        truth = zp - 2.5 * math.log10(5 * s_true)
        out['corr%g' % corr] = dict(truth=truth, median=float(np.median(mag)), diff=float(np.median(mag) - truth), factor=f_)
    # 3. varying noise, sources, masked region
    cases = {}
    sig_grad = 1.0 + 1.5 * xx / sh[1]
    sig_patch = np.ones(sh)
    sig_patch[:, 600:] *= 1.5
    sig_patch[600:, :] *= 1.4
    for name, smap in (('gradient', sig_grad), ('patches', sig_patch)):
        rng = np.random.default_rng(4)
        nse = noise_image(sh, smap, 1.0, rng)
        img = nse.copy()
        # crowded sources: 700 galaxies/stars of varied flux
        for k in range(700):
            x, y = rng.uniform(15, sh[1] - 15), rng.uniform(15, sh[0] - 15)
            f = rng.lognormal(math.log(50), 1.0)
            s = rng.uniform(1.5, 3.5)
            x0, y0 = int(x) - 12, int(y) - 12
            Y, X = np.mgrid[y0:y0 + 25, x0:x0 + 25]
            img[y0:y0 + 25, x0:x0 + 25] += f * np.exp(-((X - x) ** 2 + (Y - y) ** 2) / (2 * s * s)) / (2 * math.pi * s * s)
        bad = np.zeros(sh, bool)
        bad[100:300, 100:400] = True
        img2 = np.where(bad, np.nan, img).astype('f4')
        mag, P, f_, rm = depth_of(img2, bad=None, r=r, zp=zp)
        # compare in 3x3 blocks
        rows = []
        for j in range(3):
            for i in range(3):
                box = (i * 400, j * 400, (i + 1) * 400 - 1, (j + 1) * 400 - 1)
                if bad[box[1]:box[3], box[0]:box[2]].mean() > 0.3:
                    continue
                st = brute_sigma(nse, r, 3000, np.random.default_rng(10 + i + 3 * j), box=box)
                truth = zp - 2.5 * math.log10(5 * st)
                mm = float(np.nanmedian(mag[box[1]:box[3], box[0]:box[2]]))
                rows.append(dict(block=[i, j], truth=truth, model=mm, diff=mm - truth))
        d = np.array([r_['diff'] for r_ in rows])
        cases[name] = dict(rows=rows, diff_mean=float(d.mean()), diff_std=float(d.std()), diff_max_abs=float(np.abs(d).max()), n=len(rows),
                           depth_range=[float(np.nanmin(mag)), float(np.nanmax(mag))])
        # weight-map input: weight = 1/sigma^2 exact
        mag_w, *_x = depth_of(img2, r=r, zp=zp, rms=(1.0 / smap ** 2, 'weight'))
        rows_w = []
        for j in range(3):
            for i in range(3):
                box = (i * 400, j * 400, (i + 1) * 400 - 1, (j + 1) * 400 - 1)
                if bad[box[1]:box[3], box[0]:box[2]].mean() > 0.3:
                    continue
                st = brute_sigma(nse, r, 3000, np.random.default_rng(10 + i + 3 * j), box=box)
                rows_w.append(float(np.nanmedian(mag_w[box[1]:box[3], box[0]:box[2]]) - (zp - 2.5 * math.log10(5 * st))))
        cases[name]['weight_map_diff_mean'] = float(np.mean(rows_w))
        cases[name]['weight_map_diff_std'] = float(np.std(rows_w))
    out['varying'] = cases
    for k, v in out.items():
        print(k, json.dumps(v)[:300], flush=True)
    return out


def compmap(args):
    out = {}
    sh = (1000, 1000)
    zp = 25.0
    sig_levels = (1.0, 1.4, 2.0, 2.8)
    r, psf_fwhm = 3.0, 3.0
    det = cp.SepDetector(zp=zp, thresh=1.5, minarea=5, local_rms=True)
    kw = dict(kind='star', mag_min=19.0, mag_max=24.0, n_bins=14, per_image=25, zp=zp, psf_fwhm=psf_fwhm, match_radius=3.0, seed=3)
    # truth: stand-alone images
    truth = {}
    for s in sig_levels:
        img = noise_image((700, 700), s, 1.0, np.random.default_rng(int(s * 10))).astype('f4')
        res = cp.run_completeness(img, det, per_bin=80, **kw)
        mag = depth_of(img, r=r, zp=zp)[0]
        truth[s] = dict(lim50=res['lim50'], lim90=res['lim90'], depth=float(np.median(mag)))
        print('standalone', s, truth[s], flush=True)
    out['standalone'] = {str(k): v for k, v in truth.items()}
    # combined image: four vertical bands
    sg = np.ones(sh)
    for i, s in enumerate(sig_levels):
        sg[:, i * 250:(i + 1) * 250] = s
    img = noise_image(sh, sg, 1.0, np.random.default_rng(77)).astype('f4')
    mag, P, f_, rm = depth_of(img, r=r, zp=zp)
    lab = np.zeros(sh, np.int32)
    for i in range(4):
        lab[:, i * 250:(i + 1) * 250] = i + 1
    res = cp.run_completeness(img, det, per_bin=300, return_records=True, **kw)
    recs = res.pop('records')
    edges = np.linspace(19.0, 24.0, 15)
    rows = []
    for k in range(1, 5):
        sel = [x for x in recs if lab[int(round(x['y'])), int(round(x['x']))] == k]
        sm = cp.summarize(sel, edges, 0, 1, 0, 1)
        dk = float(np.nanmedian(mag[lab == k]))
        t = truth[sig_levels[k - 1]]
        rows.append(dict(region=k, sigma=sig_levels[k - 1], n_inj=len(sel), lim50=sm['lim50'], lim90=sm['lim90'], truth_lim50=t['lim50'], truth_lim90=t['lim90'],
                         d50=sm['lim50'] - t['lim50'], d90=sm['lim90'] - t['lim90'], depth=dk))
    out['regions'] = rows
    off = [r_['lim50'] - r_['depth'] for r_ in rows]
    out['offset_lim50_minus_depth'] = dict(values=off, std=float(np.std(off)))
    # the same through the CLI engine (rms-class regions)
    class A:
        pass
    print(json.dumps(out)[:600])
    return out


def real(args):
    import subprocess
    out = {}
    specs = (('hudf_f160w', args.hudf, 25.94, 0.06, 5.8, 10, ['--mag-min', '26', '--mag-max', '31', '--n-bins', '14', '--per-bin', '250', '--psf-fwhm', '3']),
             ('m51', args.m51, 25.0, 1.8, 2.0, 120, ['--mag-min', '13', '--mag-max', '19', '--n-bins', '12', '--per-bin', '80', '--psf-fwhm', '3']))
    for name, img, zp, ps, rad, box, kw in specs:
        if not (img and os.path.exists(img)):
            continue
        for reg in ('rms', 'grid'):
            wd = os.path.join(args.workdir, name + '_' + reg)
            cmd = [sys.executable, os.path.join(ROOT, 'plugins', 'completeness', 'depthmap.py'), img, '--work', wd, '--mode', 'compmap', '--mag-zeropoint', str(zp),
                   '--pixel-scale', str(ps), '--aper-radius', str(rad), '--regions', reg, '--rms-classes', '4', '--grid', '3x3', '--tile', '200',
                   '--box-arcsec', str(box), '--n-workers', '4'] + kw
            subprocess.run(cmd, check=True, capture_output=True)
            d = json.load(open(os.path.join(wd, 'compmap_summary.json')))
            key = name + '_' + reg
            out[key] = dict(mag_limit=d['mag_limit'], sb_limit=d.get('sb_limit'), factor=d['factor'], tile_check=d['tile_check'], valid_fraction=d['valid_fraction'],
                            compmap={k: d['compmap'][k] for k in ('global_lim50', 'global_lim90', 'regions', 'collapse', 'n_injected')})
            if reg == 'rms':
                out[key]['area_vs_depth'] = [x for x in d['area_vs_depth'] if abs(x['mag'] * 2 - round(x['mag'] * 2)) < 1e-6]
            print(key, json.dumps(out[key])[:600], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['synth', 'compmap', 'real'])
    ap.add_argument('--out')
    ap.add_argument('--hudf')
    ap.add_argument('--m51')
    ap.add_argument('--workdir', default='/tmp/depth_real')
    a = ap.parse_args()
    res = globals()[a.mode](a)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(cp.json_clean(res), f, indent=1)


if __name__ == '__main__':
    main()
