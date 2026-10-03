#!/usr/bin/env python3
"""Validation of the bar / ring / spiral components (ogfkit/structfit.py, multifit presets bar, ring, spiral, bulge+disk+bar, bulge+disk+ring).

    structure_validate.py synth OUT.json [--n 24] [--snr hi|lo] [--workers 4]
Synthetic galaxies (101 x 101, Gaussian PSF FWHM 3 px, sky 10, rms 0.5, catalog-like start errors): control (bulge + disc), barred (Ferrers bar, PA offset 20-90 deg from the disc),
ringed (Gaussian ring in the disc plane), spiral (disc with two logarithmic arms, GALFIT rotation + m = 2 mode).  Reported: classification (selected model vs truth), false-positive rate on the
controls, parameter recovery of the selected feature, and the bias of the bulge-to-total ratio when the feature is ignored (base bulge + disc fit)."""
import json
import math
import multiprocessing as mp
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('OMP_NUM_THREADS', '1')
from ogfkit import multifit as MF, models as M, structfit as SF  # noqa: E402

ZP = 25.0
N = 101
SKY = 10.0
RMS = 0.5


def draw(kind, rng, snr):
    tot = rng.uniform(9000, 20000) if snr == 'hi' else rng.uniform(2500, 5000)
    cx, cy = 50 + rng.uniform(-1, 1), 50 + rng.uniform(-1, 1)
    pa = rng.uniform(-80, 80)
    qd = rng.uniform(0.55, 0.9) if kind != 'spiral' else rng.uniform(0.75, 0.95)
    fb = rng.uniform(0.1, 0.25)
    red = rng.uniform(9.0, 14.0)
    comps = [dict(kind='sersic', x=cx, y=cy, flux=fb * tot, re=rng.uniform(2.0, 4.0), n=rng.uniform(1.5, 3.5), q=rng.uniform(0.7, 0.95), pa=pa + rng.uniform(-20, 20))]
    info = dict(kind=kind, tot=tot)
    if kind in ('control', 'bar', 'ring'):
        fx = {'control': 0.0, 'bar': rng.uniform(0.15, 0.3), 'ring': rng.uniform(0.15, 0.3)}[kind]
        comps.append(dict(kind='exp', x=cx, y=cy, flux=(1 - fb - fx) * tot, re=red, q=qd, pa=pa))
        if kind == 'bar':
            L = rng.uniform(0.8, 1.1) * red
            pab = pa + rng.choice([-1, 1]) * rng.uniform(25, 85)
            comps.append(MF.bar_component(cx, cy, fx * tot, L, pab, q=rng.uniform(0.2, 0.32), c0=rng.uniform(0.3, 1.5)))
            info.update(bar_pa=((pab + 90) % 180) - 90, bar_len=L, bar_frac=fx)
        if kind == 'ring':
            rr = rng.uniform(0.7, 1.1) * red
            comps.append(dict(kind='gring', x=cx, y=cy, flux=fx * tot, rring=rr, sring=rng.uniform(1.5, 3.0), q=qd, pa=pa))
            info.update(ring_r=rr, ring_s=comps[-1]['sring'], ring_frac=fx)
    else:
        th = rng.choice([-1, 1]) * rng.uniform(180, 320)
        sp = MF.spiral_disc(cx, cy, (1 - fb) * tot, red, qd, pa, theta=th, ampl=rng.uniform(0.12, 0.3))
        sp['rot_out'] = 3.0 * red
        comps.append(sp)
        info.update(spiral_theta=th, spiral_amp=sp['f2a'])
    info['disc_re'] = red
    info['bulge_frac'] = fb
    return comps, info, (cx, cy, pa, qd)


def one(args):
    kind, seed, snr = args
    rng = np.random.default_rng(97 * ord(kind[0]) + seed + (5000 if snr == 'lo' else 0))
    comps, info, (cx, cy, pa, qd) = draw(kind, rng, snr)
    psf = M.gaussian_psf(21, 3.0)
    cn = [MF._normalise(c, ZP) for c in comps]
    img = MF.render_model(cn, (N, N), psf=psf) + SKY + rng.normal(0, RMS, (N, N))
    tot = info['tot']
    t0 = time.time()
    try:
        r = SF.fit_structure(img, cx + rng.uniform(-1, 1), cy + rng.uniform(-1, 1), tot * rng.uniform(0.85, 1.15), info['disc_re'] * rng.uniform(0.6, 1.0) * 1.0, float(np.clip(qd + rng.uniform(-0.15, 0.15), 0.3, 1.0)),
                             pa + rng.uniform(-15, 15), psf=psf, rms=RMS, zp=ZP)
    except Exception as e:
        return dict(info, error=str(e)[:100])
    out = dict(info, seed=seed, snr=snr, selected=r['name'], notes=r['notes'], bic={k: float(v) for k, v in r['bic'].items()}, chi2red={k: float(v['chi2_red']) for k, v in r['all'].items()}, time=time.time() - t0)
    base = r['all'].get('base')
    if base:
        cs = base['components']
        b = [c for c in cs if c['kind'] == 'sersic'][0]
        d = [c for c in cs if c['kind'] == 'exp'][0]
        fb_, fd_ = 10 ** (-0.4 * b['mag']), 10 ** (-0.4 * d['mag'])
        out['base_BT'] = fb_ / (fb_ + fd_)
    if r['res'] is not None:
        cs = r['res']['components']
        tt = sum(10 ** (-0.4 * c['mag']) for c in cs)
        out['sel_BT'] = 10 ** (-0.4 * cs[0]['mag']) / tt
        if r['name'] == 'bar':
            s = SF.bar_summary(r['res'])
            out['fit_bar'] = s
        if r['name'] == 'ring':
            out['fit_ring'] = dict(r=cs[2]['rring'], s=cs[2]['sring'], frac=10 ** (-0.4 * cs[2]['mag']) / tt)
        if r['name'] == 'spiral':
            out['fit_spiral'] = dict(theta=cs[1]['rot_theta'], amp=cs[1]['f2a'])
    # truth B/T of the bulge among all light
    out['true_BT'] = out['bulge_frac']
    return out


def wrap180(d):
    return ((d + 90.0) % 180.0) - 90.0


def summarise(res):
    s = {}
    ok = [r for r in res if 'error' not in r]
    s['errors'] = len(res) - len(ok)
    expect = dict(control='base', bar='bar', ring='ring', spiral='spiral')
    for k in expect:
        rs = [r for r in ok if r['kind'] == k]
        d = dict(n=len(rs), selected={m: sum(1 for r in rs if r['selected'] == m) for m in ('base', 'bar', 'ring', 'spiral')}, correct=sum(1 for r in rs if r['selected'] == expect[k]) / max(len(rs), 1),
                 time_median=float(np.median([r['time'] for r in rs])) if rs else None)
        bt0 = [r['base_BT'] - r['true_BT'] for r in rs if 'base_BT' in r]
        bt1 = [r['sel_BT'] - r['true_BT'] for r in rs if 'sel_BT' in r and r['selected'] == expect[k]]
        d['BT_err_base_fit_median'] = float(np.median(bt0)) if bt0 else None
        d['BT_err_selected_median'] = float(np.median(bt1)) if bt1 else None
        d['BT_abs_err_base_median'] = float(np.median(np.abs(bt0))) if bt0 else None
        d['BT_abs_err_selected_median'] = float(np.median(np.abs(bt1))) if bt1 else None
        if k == 'bar':
            f = [r for r in rs if r['selected'] == 'bar']
            d['bar_pa_err_deg_median'] = float(np.median([abs(wrap180(r['fit_bar']['pa'] - r['bar_pa'])) for r in f])) if f else None
            d['bar_pa_err_deg_p90'] = float(np.percentile([abs(wrap180(r['fit_bar']['pa'] - r['bar_pa'])) for r in f], 90)) if f else None
            d['bar_len_ratio_median'] = float(np.median([r['fit_bar']['length'] / r['bar_len'] for r in f])) if f else None
            d['bar_len_ratio_nmad'] = float(1.4826 * np.median(np.abs(np.array([r['fit_bar']['length'] / r['bar_len'] for r in f]) - d['bar_len_ratio_median']))) if f else None
            d['bar_frac_diff_median'] = float(np.median([r['fit_bar']['bar_to_total'] - r['bar_frac'] for r in f])) if f else None
            d['bar_frac_diff_p90_abs'] = float(np.percentile([abs(r['fit_bar']['bar_to_total'] - r['bar_frac']) for r in f], 90)) if f else None
        if k == 'ring':
            f = [r for r in rs if r['selected'] == 'ring']
            d['ring_r_ratio_median'] = float(np.median([r['fit_ring']['r'] / r['ring_r'] for r in f])) if f else None
            d['ring_s_ratio_median'] = float(np.median([r['fit_ring']['s'] / r['ring_s'] for r in f])) if f else None
            d['ring_frac_diff_median'] = float(np.median([r['fit_ring']['frac'] - r['ring_frac'] for r in f])) if f else None
        if k == 'spiral':
            f = [r for r in rs if r['selected'] == 'spiral']
            d['spiral_sign_correct'] = float(np.mean([np.sign(r['fit_spiral']['theta']) == np.sign(r['spiral_theta']) for r in f])) if f else None
            d['spiral_theta_ratio_median'] = float(np.median([r['fit_spiral']['theta'] / r['spiral_theta'] for r in f])) if f else None
            d['spiral_amp_diff_median'] = float(np.median([r['fit_spiral']['amp'] - r['spiral_amp'] for r in f])) if f else None
        s[k] = d
    s['false_positive_rate_controls'] = 1.0 - s['control']['correct'] if s['control']['n'] else None
    return s


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['synth'])
    ap.add_argument('out')
    ap.add_argument('--n', type=int, default=24)
    ap.add_argument('--snr', default='hi', choices=['hi', 'lo'])
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--kinds', default='control,bar,ring,spiral')
    a = ap.parse_args()
    tasks = [(k, s, a.snr) for k in a.kinds.split(',') for s in range(a.n)]
    t0 = time.time()
    with mp.Pool(a.workers) as p:
        res = p.map(one, tasks, chunksize=1)
    rep = dict(snr=a.snr, n_per_class=a.n, elapsed_s=time.time() - t0, summary=summarise(res), results=res)
    json.dump(rep, open(a.out, 'w'), default=lambda o: None)
    sm = rep['summary']
    print('elapsed %.0f s, errors %d, false positives on controls %.2f' % (rep['elapsed_s'], sm['errors'], sm['false_positive_rate_controls'] or 0))
    for k in ('control', 'bar', 'ring', 'spiral'):
        d = sm[k]
        print(k, d['n'], 'correct %.2f' % d['correct'], d['selected'], 'BT|err| base %.3f sel %s' % (d['BT_abs_err_base_median'] or 0, d['BT_abs_err_selected_median'] and round(d['BT_abs_err_selected_median'], 3)))


if __name__ == '__main__':
    main()
