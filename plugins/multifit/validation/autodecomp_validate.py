#!/usr/bin/env python3
"""Validation of the automatic structure decomposition (ogfkit/autodecomp.py).

    autodecomp_validate.py synth [OUT] [--n 40] [--snr hi|lo]   synthetic galaxies of known structure (single Sersic, exponential disc, bulge+disc,
                                                                  nucleus+bulge+disc), Gaussian PSF FWHM 3 px: type confusion, B/T and parameter recovery,
                                                                  1-D guess quality, comparison with the fixed bulge+disk preset fit and auto_select
    autodecomp_validate.py real OUT                              real HUDF F160W / M51 objects through the multifit CLI (--model decomp)
Results: JSON next to this file (or in OUT)."""
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import multifit as mf, autodecomp as ad  # noqa: E402

ZP = 25.0
FWHM = 3.0
N = 91


def mag(f):
    return ZP - 2.5 * math.log10(f)


def make_truth(kind, rng, flux):
    """Component list (multifit dicts, mag-based) of a galaxy of the given class and the truth values."""
    q_d = rng.uniform(0.35, 0.95)
    pa = rng.uniform(0, 180)
    cx, cy = 45 + rng.uniform(-0.5, 0.5), 45 + rng.uniform(-0.5, 0.5)
    t = dict(kind=kind, flux=flux, pa=pa)
    if kind == 'sersic':
        n = rng.choice([1.5, 2.5, 3.5, 4.5])
        re = rng.uniform(5, 11)
        comps = [dict(kind='sersic', x=cx, y=cy, mag=mag(flux), re=re, n=n, q=rng.uniform(0.5, 0.95), pa=pa)]
        t.update(re=re, n=n)
    elif kind == 'disc':
        h = rng.uniform(5, 11)
        comps = [dict(kind='exp', x=cx, y=cy, mag=mag(flux), re=1.678 * h, q=q_d, pa=pa)]
        t.update(re=1.678 * h, n=1.0)
    else:
        bt = rng.choice([0.15, 0.3, 0.5, 0.7])
        nb = rng.choice([1.5, 2.5, 4.0])
        reb = rng.uniform(1.8, 4.0)
        h = rng.uniform(5.5, 9.0)
        qb = rng.uniform(0.7, 1.0)
        nuc = float(rng.choice([0.06, 0.12, 0.20])) * flux if kind == 'nucleus+bulge+disk' else 0.0
        fb, fd = bt * (flux - nuc), (1 - bt) * (flux - nuc)
        comps = ([dict(kind='psf', x=cx, y=cy, mag=mag(nuc))] if nuc else []) + [
            dict(kind='sersic', x=cx, y=cy, mag=mag(fb), re=reb, n=nb, q=qb, pa=pa + rng.uniform(-20, 20)),
            dict(kind='exp', x=cx, y=cy, mag=mag(fd), re=1.678 * h, q=q_d, pa=pa)]
        t.update(bt=bt, nb=nb, reb=reb, red=1.678 * h, qb=qb, qd=q_d, nuc=nuc / flux)
    t['comps'] = comps
    return t


def one(args):
    kind, seed, flux, rms = args
    rng = np.random.default_rng(seed)
    t = make_truth(kind, rng, flux)
    psf = mf.gaussian_psf_array(FWHM)
    model = mf.render_model([mf._normalise(c, ZP) for c in t['comps']], (N, N), psf=psf, sky=0.0)
    data = model + rng.normal(0, rms, model.shape)
    # catalog-like start values with realistic errors
    yy, xx = np.mgrid[:N, :N]
    cum = np.sort(model.ravel())[::-1].cumsum()
    r50 = float(np.sqrt((cum < 0.5 * cum[-1]).sum() / math.pi))                         # circular half-light radius of the true model
    q0 = float(np.clip(t['comps'][-1].get('q', 0.8) + rng.uniform(-0.15, 0.15), 0.2, 1.0))
    pa0 = float(t['pa'] + rng.uniform(-15, 15))
    f0 = flux * rng.uniform(0.8, 1.2)
    re0 = max(1.5, r50 * rng.uniform(0.7, 1.4))
    x0, y0 = t['comps'][0]['x'], t['comps'][0]['y']
    out = dict(kind=kind, seed=seed, flux=flux, truth={k: v for k, v in t.items() if k != 'comps'})
    t0 = time.time()
    d = ad.decompose(data, x0 + rng.uniform(-0.3, 0.3), y0 + rng.uniform(-0.3, 0.3), psf=psf, rms=rms, zp=ZP, flux0=f0, re0=re0, q0=q0, pa0=pa0, fwhm=FWHM, nucleus='auto', restarts=1, max_nfev=150)
    out['t_auto'] = time.time() - t0
    out['sel'] = d['name']
    g = d.get('guess', {})
    out['guess_bt'] = g.get('bt_guess')
    out['guess_fallback'] = g.get('guess_fallback')
    out['guess_ok'] = bool(g.get('ok'))
    out['nuc_ratio'] = g.get('nucleus_ratio')
    out['flag'] = d['flag']
    out['notes'] = d['notes']
    out['bic'] = d['bic']
    out['scale'] = d.get('scale')
    if d['res'] is not None:
        cs = d['res']['components'][:d['ntarget']]
        out['chi2_red'] = d['res']['chi2_red']
        out['comps'] = [{k: c[k] for k in ('kind', 'mag', 're', 'n', 'q', 'pa') if k in c} for c in cs]
        tot = sum(c['flux'] for c in cs)
        out['mag_tot_err'] = mag(tot) - mag(flux)
        if d['name'] != 'sersic':
            b, dd = cs[-2], cs[-1]
            out['bt'] = b['flux'] / (b['flux'] + dd['flux'])
    # baselines on the same data and starts: fixed preset bulge+disk (dev + exp), auto_select among psf / sersic / bulge+disk
    t0 = time.time()
    kw = dict(psf=psf, rms=rms, zp=ZP, max_nfev=150)
    try:
        rb, _ = mf.fit_preset(data, 'bulge+disk', x0, y0, f0, re0, q0, pa0, **kw)
        cb = rb['components']
        out['base_bt'] = cb[0]['flux'] / (cb[0]['flux'] + cb[1]['flux'])
        out['base_chi2_red'] = rb['chi2_red']
        name, rsel, _ = mf.auto_select(data, x0, y0, f0, re0, q0, pa0, **kw)
        out['base_sel'] = name
    except Exception as e:
        out['base_err'] = str(e)[:60]
    out['t_base'] = time.time() - t0
    return out


def summarise(res):
    rep = {}
    kinds = ['sersic', 'disc', 'bulge+disk', 'nucleus+bulge+disk']
    truth_class = {'sersic': 'sersic', 'disc': 'sersic', 'bulge+disk': 'bulge+disk', 'nucleus+bulge+disk': 'nucleus+bulge+disk'}
    conf = {k: {} for k in kinds}
    for r in res:
        conf[r['kind']][r['sel']] = conf[r['kind']].get(r['sel'], 0) + 1
    rep['confusion (rows: truth, columns: selected)'] = conf
    ok = [r for r in res if r['sel'] == truth_class[r['kind']]]
    rep['accuracy_exact'] = len(ok) / len(res)
    simple = [r for r in res if r['kind'] in ('sersic', 'disc')]
    rep['single_correct_frac'] = float(np.mean([r['sel'] == 'sersic' for r in simple])) if simple else None
    bd = [r for r in res if r['kind'] in ('bulge+disk', 'nucleus+bulge+disk')]
    rep['two_comp_found_frac'] = float(np.mean([r['sel'] in ('bulge+disk', 'nucleus+bulge+disk') for r in bd])) if bd else None
    # base_sel: 2-comp or not
    rep['baseline_auto_select'] = dict(single_correct=float(np.mean([r.get('base_sel') == 'sersic' for r in simple])) if simple else None,
                                        two_comp_found=float(np.mean([r.get('base_sel') == 'bulge+disk' for r in bd])) if bd else None)
    fnd = [r for r in bd if 'bt' in r and r['sel'] != 'sersic']
    if fnd:
        d = np.array([r['bt'] - r['truth']['bt'] for r in fnd])
        rep['bt_auto'] = dict(n=len(d), bias=float(d.mean()), scatter=float(d.std()), nmad=float(1.4826 * np.median(np.abs(d - np.median(d)))), within_0p1=float(np.mean(np.abs(d) < 0.1)))
        rep['bt_by_truth'] = {}
        for v in sorted(set(r['truth']['bt'] for r in fnd)):
            dd = np.array([r['bt'] - v for r in fnd if r['truth']['bt'] == v])
            rep['bt_by_truth'][str(v)] = dict(n=len(dd), bias=float(dd.mean()), scatter=float(dd.std()))
        for key, f in (('reb_ratio', lambda r: r['comps'][-2]['re'] / r['truth']['reb']), ('red_ratio', lambda r: r['comps'][-1]['re'] / r['truth']['red']),
                       ('nb_diff', lambda r: r['comps'][-2]['n'] - r['truth']['nb']), ('qd_diff', lambda r: r['comps'][-1]['q'] - r['truth']['qd']),
                       ('magtot_err', lambda r: r['mag_tot_err'])):
            v = np.array([f(r) for r in fnd])
            rep[key] = dict(median=float(np.median(v)), nmad=float(1.4826 * np.median(np.abs(v - np.median(v)))))
    both = [r for r in bd if 'bt' in r and r['sel'] != 'sersic' and 'base_bt' in r]
    if both:
        da = np.array([abs(r['bt'] - r['truth']['bt']) for r in both])
        db = np.array([abs(r['base_bt'] - r['truth']['bt']) for r in both])
        rep['bt_vs_baseline_fixed_preset (same objects)'] = dict(n=len(both), auto_median_abs=float(np.median(da)), baseline_median_abs=float(np.median(db)), auto_within_0p1=float(np.mean(da < 0.1)),
                                                                  baseline_within_0p1=float(np.mean(db < 0.1)))
    g = [r for r in bd if r.get('guess_bt') is not None and not r.get('guess_fallback')]
    if g:
        d = np.array([r['guess_bt'] - r['truth']['bt'] for r in g])
        rep['bt_1d_guess'] = dict(n=len(g), bias=float(d.mean()), scatter=float(d.std()), fallback_frac=float(np.mean([bool(r.get('guess_fallback')) for r in bd])))
    rep['nucleus_detected_frac (nucleus truth)'] = float(np.mean([r['sel'] == 'nucleus+bulge+disk' for r in res if r['kind'] == 'nucleus+bulge+disk'])) if any(r['kind'] == 'nucleus+bulge+disk' for r in res) else None
    rep['nucleus_false_frac (no nucleus truth)'] = float(np.mean([r['sel'] == 'nucleus+bulge+disk' for r in res if r['kind'] != 'nucleus+bulge+disk']))
    rep['chi2_red_median'] = float(np.median([r['chi2_red'] for r in res if 'chi2_red' in r]))
    rep['time_auto_median_s'] = float(np.median([r['t_auto'] for r in res]))
    rep['time_baseline_median_s'] = float(np.median([r['t_base'] for r in res]))
    return rep


def synth(out, n=40, snr='hi'):
    flux_rng = (1.5e4, 6e4) if snr == 'hi' else (3e3, 1e4)
    jobs = []
    for ki, kind in enumerate(['sersic', 'disc', 'bulge+disk', 'nucleus+bulge+disk']):
        for i in range(n):
            rng = np.random.default_rng(777 + 1000 * ki + i)
            jobs.append((kind, 5000 * ki + i, float(np.exp(rng.uniform(math.log(flux_rng[0]), math.log(flux_rng[1])))), 1.0))
    t0 = time.time()
    with mp.get_context('fork').Pool(8) as pool:
        res = pool.map(one, jobs, chunksize=1)
    rep = summarise(res)
    rep['n_per_class'] = n
    rep['snr'] = snr
    rep['wall_s'] = time.time() - t0
    json.dump(dict(summary=rep, objects=res), open(os.path.join(out, 'autodecomp_synth_%s.json' % snr), 'w'), indent=1, default=float)
    json.dump(rep, open(os.path.join(out, 'autodecomp_synth_%s_summary.json' % snr), 'w'), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float)[:3500])


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'synth'
    out = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else HERE
    os.makedirs(out, exist_ok=True)
    n = int(sys.argv[sys.argv.index('--n') + 1]) if '--n' in sys.argv else 40
    snr = sys.argv[sys.argv.index('--snr') + 1] if '--snr' in sys.argv else 'hi'
    if mode == 'synth':
        synth(out, n, snr)
