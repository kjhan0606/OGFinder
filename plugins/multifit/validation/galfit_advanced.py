#!/usr/bin/env python3
"""Validation of the advanced GALFIT components of multifit (Fourier F*, bending B*, rotation R*, truncation, C0, Moffat / Ferrer / King / Nuker / edge-on disc).

    python galfit_advanced.py parity OUT.json [--galfit BIN] [--ld DIR]     rendered models vs the real GALFIT (P=1), with and without PSF
    python galfit_advanced.py h2h    OUT.json [--n 6] [--workers 8]         head-to-head fits on GALFIT-rendered data (multifit CLI vs GALFIT, same start feedme / constraints)
    python galfit_advanced.py recovery OUT.json [--n 12]                    synthetic truth recovery with ogfkit's own renderer (no GALFIT needed)

All model definitions are config dicts of ogfkit.multifit (1-based image coordinates); GALFIT is driven through ogfkit.galfitio.write_feedme, so the export is tested as well."""
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

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
from astropy.io import fits  # noqa: E402
from ogfkit import galfitio as GI, multifit as MF  # noqa: E402

ZP = 25.0
NY, NX = 101, 101
MULTIFIT = os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py')
SIGMA = 2.0
SKY = 50.0


def moffat_psf(fwhm=3.2, beta=3.0, size=31):
    a = fwhm / (2 * math.sqrt(2 ** (1 / beta) - 1))
    y, x = np.mgrid[:size, :size] - size // 2
    p = (1 + (x ** 2 + y ** 2) / a ** 2) ** (-beta)
    return (p / p.sum()).astype(np.float32)


def run_galfit(galfit, ld, feed, cwd, timeout=300):
    env = dict(os.environ)
    if ld:
        env['LD_LIBRARY_PATH'] = ld + ':' + env.get('LD_LIBRARY_PATH', '')
    t = time.time()
    r = subprocess.run([galfit, feed], cwd=cwd, env=env, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=timeout)
    return time.time() - t, r


S = lambda **k: dict(dict(kind='sersic', x=51.3, y=50.6, mag=16.5, re=9.0, n=2.0, q=0.7, pa=120.0), **k)

# name -> components (GALFIT-convention config dicts, 1-based coordinates)
CASES = {
    'sersic': [S()],
    'sersic_n0.8': [S(n=0.8)],
    'sersic_n4_small': [dict(kind='sersic', x=51.3, y=50.6, mag=16.5, re=4.0, n=4.0, q=0.8, pa=70.0)],
    'c0_boxy': [S(c0=0.8)],
    'c0_disky': [S(c0=-0.5)],
    'fourier_f1': [S(f1a=0.1, f1p=30.0)],
    'fourier_f3': [S(f3a=0.08, f3p=-20.0)],
    'fourier_f1f3f5': [S(f1a=0.05, f1p=10.0, f3a=0.06, f3p=-40.0, f5a=0.03, f5p=70.0)],
    'fourier_f2f4_q0.4': [S(q=0.4, f2a=0.1, f2p=15.0, f4a=0.05, f4p=-60.0)],
    'bending_b1': [S(b1=0.15)],
    'bending_b2': [S(b2=-0.1)],
    'bending_b1b2': [S(b1=0.1, b2=0.08)],
    'rot_power': [S(rot_func='power', rot_in=0.0, rot_out=25.0, rot_theta=60.0, rot_alpha=0.5, rot_incl=0.0, rot_pa=0.0)],
    'rot_log': [S(rot_func='log', rot_in=0.0, rot_out=25.0, rot_theta=-80.0, rot_ws=3.0, rot_incl=0.0, rot_pa=0.0)],
    'rot_power_incl': [S(rot_func='power', rot_in=1.0, rot_out=30.0, rot_theta=90.0, rot_alpha=-0.3, rot_incl=50.0, rot_pa=20.0)],
    'trunc_outer': [S(norm='re', i0=None, trunc_out=[1]), dict(kind='trunc', rbreak=22.0, dsoft=8.0)],
    'trunc_inner': [S(norm='re', i0=None, trunc_in=[1]), dict(kind='trunc', rbreak=16.0, dsoft=6.0)],
    'trunc_both': [S(norm='re', i0=None, trunc_in=[1], trunc_out=[2]), dict(kind='trunc', rbreak=10.0, dsoft=5.0), dict(kind='trunc', rbreak=30.0, dsoft=8.0)],
    'sersic_centre_norm': [S(norm='center', i0=None)],
    'moffat': [dict(kind='moffat', x=51.3, y=50.6, mag=17.0, fwhm=6.0, beta=2.5, q=0.8, pa=30.0)],
    'ferrer': [dict(kind='ferrer', x=51.3, y=50.6, mu=20.0, rout=30.0, alpha=3.0, beta=1.0, q=0.7, pa=60.0)],
    'king': [dict(kind='king', x=51.3, y=50.6, mu=20.0, rc=5.0, rt=35.0, alpha=2.0, q=0.8, pa=-30.0)],
    'nuker': [dict(kind='nuker', x=51.3, y=50.6, mu=20.0, rb=10.0, alpha=1.5, beta=1.8, gamma=0.5, q=0.8, pa=10.0)],
    'edgedisk': [dict(kind='edgedisk', x=51.3, y=50.6, mu=19.5, hs=3.0, rs=14.0, pa=45.0)],
    'combo': [S(f2a=0.08, f2p=20.0, b1=0.08, rot_func='power', rot_in=0.0, rot_out=25.0, rot_theta=40.0, rot_alpha=0.5, rot_incl=0.0, rot_pa=0.0)],
    'quirk_c0_plus_fourier': [S(f2a=0.08, f2p=20.0, c0=0.3)],            # GALFIT normalises C0 + Fourier differently (see docs/multifit.md): reported, not counted
}


def prepare(comps):
    out = []
    for c in comps:
        c = {k: v for k, v in c.items() if v is not None}
        if 'mu' in c:
            c['i0' if c['kind'] != 'nuker' else 'ib'] = 10 ** (-0.4 * (c.pop('mu') - ZP))
        if c.get('norm'):
            c.pop('mag', None)
        out.append(c)
    return out


def fill_sb(comps):
    """norm != total components: choose i0 so the model is bright (counts per pixel)."""
    for c in comps:
        if c.get('norm') and 'i0' not in c:
            c['i0'] = {'re': 6.0, 'center': 40.0, 'break': 10.0}[c['norm']]
    return comps


def render_ours(comps, psf, sky=0.0):
    cs = [MF._normalise(json.loads(json.dumps(c)), ZP) for c in comps]
    for c in cs:
        if 'x' in c:
            c['x'] -= 1.0
            c['y'] -= 1.0
    return MF.render_model(cs, (NY, NX), psf=psf, sky=sky)


def galfit_render(args, comps, psf, wd, sky=0.0):
    os.makedirs(wd, exist_ok=True)
    fits.PrimaryHDU(np.zeros((NY, NX), np.float32)).writeto(os.path.join(wd, 'data.fits'), overwrite=True)
    if psf is not None:
        fits.PrimaryHDU(psf.astype(np.float32)).writeto(os.path.join(wd, 'psf.fits'), overwrite=True)
    cfg = dict(components=comps, zp=ZP, exptime=1.0, sky='fixed', sky_value=sky, plate_scale=[1.0, 1.0])
    txt = GI.write_feedme(cfg, image='data.fits', output='out.fits', psf='psf.fits' if psf is not None else 'none', shape=(NY, NX), mode=1, plate_scale=(1.0, 1.0))
    open(os.path.join(wd, 'm.feedme'), 'w').write(txt)
    if os.path.exists(os.path.join(wd, 'out.fits')):
        os.remove(os.path.join(wd, 'out.fits'))
    run_galfit(args.galfit, args.ld, 'm.feedme', wd)
    o = os.path.join(wd, 'out.fits')
    return fits.getdata(o).astype(float) if os.path.exists(o) else None


def parity(args):
    psf = moffat_psf()
    res = {}
    for with_psf in (False, True):
        for name, comps in CASES.items():
            comps = fill_sb(prepare(json.loads(json.dumps(comps))))
            key = name + ('+psf' if with_psf else '')
            wd = os.path.join(args.work, key.replace('+', '_'))
            shutil.rmtree(wd, ignore_errors=True)
            try:
                g = galfit_render(args, comps, psf if with_psf else None, wd)
                m = render_ours(comps, psf if with_psf else None)
            except Exception as e:
                res[key] = dict(error=str(e)[:200])
                continue
            if g is None or not np.isfinite(g).all() or g.max() <= 0:
                res[key] = dict(error='GALFIT produced no usable model')
                continue
            d = m - g
            res[key] = dict(max_over_peak=float(np.abs(d).max() / g.max()), rms_over_peak=float(np.sqrt((d ** 2).mean()) / g.max()), flux_diff=float((m.sum() - g.sum()) / g.sum()), peak=float(g.max()))
            print('%-24s max|d|/peak %.4f  rms/peak %.5f  flux %+.4f' % (key, res[key]['max_over_peak'], res[key]['rms_over_peak'], res[key]['flux_diff']))
    ok = [v for k, v in res.items() if 'error' not in v and not k.startswith('quirk')]
    summ = dict(n=len(ok), n_error=sum('error' in v for v in res.values()), worst_max_over_peak=max(v['max_over_peak'] for v in ok), median_max_over_peak=float(np.median([v['max_over_peak'] for v in ok])),
                worst_flux_diff=max(abs(v['flux_diff']) for v in ok), within_1pct=sum(v['max_over_peak'] < 0.01 for v in ok), within_2pct=sum(v['max_over_peak'] < 0.02 for v in ok))
    json.dump(dict(summary=summ, cases=res), open(args.out, 'w'), indent=1)
    print(summ)
    return 0


# ---------------------------------------------------------------------------------------------------------------------------------- head-to-head fits / truth recovery
REL = {'re', 'n', 'fwhm', 'rbreak', 'rout', 'rc', 'rt', 'rb', 'hs', 'rs', 'h1', 'h2', 'i0', 'ib', 'dsoft', 'rot_out'}
ANG = {'pa', 'rot_pa'}
BOUNDS = {'q': (0.1, 1.0), 'n': (0.3, 8.0), 're': (0.5, 60.0), 'beta': (0.8, 10.0), 'c0': (-2.0, 2.5), 'rbreak': (4.0, 60.0), 'dsoft': (1.0, 30.0), 'fwhm': (1.0, 40.0)}


def _bounds_for(c):
    b = {}
    for k in MF.param_names(c):
        if k in BOUNDS:
            b[k] = list(BOUNDS[k])
        elif re.match(r'^f\d+a$', k):
            b[k] = [-0.6, 0.6]
    return b


def case_defs():
    """name -> dict(truth=fn(rng) -> comps, free=lambda comp -> set of free parameter names or None (= all but defaults), mask=bool)."""
    def sers(rng, **k):
        return dict(dict(kind='sersic', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), mag=rng.uniform(16.3, 17.2), re=rng.uniform(5, 9), n=rng.uniform(1.2, 3.0), q=rng.uniform(0.5, 0.9), pa=rng.uniform(0, 180)), **k)
    D = {}
    D['F_fourier'] = lambda rng: [sers(rng, f3a=rng.uniform(0.04, 0.09), f3p=rng.uniform(-40, 40), f1a=rng.uniform(0.03, 0.07), f1p=rng.uniform(-60, 60))]
    D['G_bending'] = lambda rng: [sers(rng, b1=rng.uniform(0.05, 0.12), b2=rng.uniform(-0.06, 0.06))]
    D['H_rotation'] = lambda rng: [sers(rng, rot_func='power', rot_in=0.0, rot_out=rng.uniform(22, 30), rot_theta=rng.uniform(40, 80) * rng.choice([-1, 1]), rot_alpha=0.5, rot_incl=0.0, rot_pa=0.0,
                                         fixed=['rot_in', 'rot_out', 'rot_alpha', 'rot_incl', 'rot_pa'])]
    D['I_truncation'] = lambda rng: [sers(rng, norm='re', i0=rng.uniform(5.0, 9.0), trunc_out=[1], fixed=['n']), dict(kind='trunc', rbreak=rng.uniform(18, 26), dsoft=7.0, fixed=['dsoft'])]
    D['J_boxy_moffat'] = lambda rng: [sers(rng, c0=rng.uniform(0.3, 0.9)), dict(kind='moffat', x=51 + rng.uniform(15, 25) * rng.choice([-1, 1]), y=51 + rng.uniform(-20, 20), mag=rng.uniform(17.5, 18.5), fwhm=rng.uniform(4.5, 7.0),
                                                                           beta=2.5, q=rng.uniform(0.6, 0.95), pa=rng.uniform(0, 180), fixed=['beta'])]
    D['K_tied'] = lambda rng: _tied(rng, sers)
    D['L_masked_neighbour'] = lambda rng: [sers(rng, f2a=0.05, f2p=rng.uniform(-30, 30))]
    return D


def _tied(rng, sers):
    a = sers(rng, n=rng.uniform(1.5, 3.0))
    b = dict(kind='exp', x=a['x'], y=a['y'], mag=a['mag'] - 0.3, re=rng.uniform(12, 18) * 1.678, q=rng.uniform(0.4, 0.8), pa=a['pa'] + rng.uniform(-20, 20))
    return [dict(a, mag=a['mag'] + 0.6, re=rng.uniform(2.5, 4.0), n=rng.uniform(2.5, 4.0), fixed=[]), b]


TIES_K = [['1.x', '0.x'], ['1.y', '0.y'], ['1.flux', '0.flux', 'ratio']]      # exp disc centred on the bulge; the disc is a fixed multiple of the bulge flux (GALFIT: x / y offset, magnitude offset)


def _mag_key(c):
    return 'mag'


def _perturb(c, rng):
    p = json.loads(json.dumps(c))
    for k in MF.param_names(c):
        if k in c.get('fixed', ()) or k == 'flux':
            continue
        v = c[k]
        if k in ('x', 'y'):
            p[k] = v + rng.uniform(-0.8, 0.8)
        elif k in ('q',):
            p[k] = float(np.clip(v + rng.uniform(-0.15, 0.15), 0.2, 1.0))
        elif k in ANG:
            p[k] = v + rng.uniform(-25, 25)
        elif re.match(r'^f\d+a$', k) or re.match(r'^b\d+$', k):
            p[k] = v * rng.uniform(0.3, 0.7) + rng.uniform(-0.01, 0.01)
        elif re.match(r'^f\d+p$', k):
            p[k] = v + rng.uniform(-25, 25)
        elif k == 'c0':
            p[k] = v + rng.uniform(-0.25, 0.25)
        elif k == 'rot_theta':
            p[k] = v * rng.uniform(0.75, 1.25)
        elif k == 'n':
            p[k] = float(np.clip(v * rng.uniform(0.7, 1.4), 0.5, 6))
        elif k in ('beta',):
            p[k] = v + rng.uniform(-0.3, 0.3)
        else:
            p[k] = v * rng.uniform(0.8, 1.25)
    if 'mag' in c and 'mag' not in c.get('fixed', ()):
        p['mag'] = c['mag'] + rng.uniform(-0.3, 0.3)
    return p


def _num_cmp(k, a, b):
    m = re.match(r'^f(\d+)p$', k)
    if m:
        per = 360.0 / int(m.group(1))
        return (a - b + per / 2) % per - per / 2
    if k in ANG:
        return (a - b + 90.0) % 180.0 - 90.0
    if k == 'rot_theta':
        return (a - b + 180.0) % 360.0 - 180.0
    if k in REL:
        return a / b - 1.0 if b else float('nan')
    return a - b


def _truth_with_flux(comps):
    out = []
    for c in comps:
        c = json.loads(json.dumps(c))
        out.append(c)
    return out


def _write_pair(wd, name, comps, tie, bounds_on=True):
    cfg = dict(components=comps, tie=tie, zp=ZP, exptime=1.0, sky='const', sky_value=SKY, plate_scale=[1.0, 1.0])
    return cfg


def _fit_feed(comps, tie, sky, mask='none'):
    comps = json.loads(json.dumps(comps))
    for c in comps:
        c['bounds'] = _bounds_for(c)
    cfg = dict(components=comps, tie=tie, zp=ZP, exptime=1.0, sky='const', sky_value=sky, plate_scale=[1.0, 1.0])
    cons = GI.write_constraints(cfg)
    has = len(cons.strip().splitlines()) > 1
    txt = GI.write_feedme(cfg, image='data.fits', output='out.fits', sigma='sigma.fits', psf='psf.fits', mask=mask, constraints='cons.txt' if has else 'none', shape=(NY, NX), mode=0, plate_scale=(1.0, 1.0))
    return txt, cons if has else None


def _result_comps(path, wd):
    cfg = GI.parse_feedme(path, strict=False, base_dir=wd, exptime=1.0)
    return cfg


def _flatten(cfg):
    """component list (config dicts) -> {(i, key): value} for the comparison."""
    out = {}
    for i, c in enumerate(cfg['components']):
        keys = MF.param_names(MF._normalise(json.loads(json.dumps(c)), ZP) if c['kind'] != 'trunc' else c) if c['kind'] != 'trunc' else ['rbreak', 'dsoft']
        for k in keys:
            if k == 'flux':
                out['c%d.mag' % i] = c['mag']
            elif k in c:
                out['c%d.%s' % (i, k)] = c[k]
    for key in [k for k in out if re.match(r'^c\d+\.f\d+a$', k)]:                 # canonical Fourier mode: amplitude >= 0 (a negative amplitude = phase shifted by 180 deg / m)
        if out[key] < 0:
            pk = key[:-1] + 'p'
            out[key] = -out[key]
            out[pk] = out.get(pk, 0.0) + 180.0 / int(re.search(r'f(\d+)a$', key).group(1))
    return out


def one(task):
    name, seed, args = task
    rng = np.random.default_rng(7919 * (ord(name[0]) - 60) + seed)
    truth = [MF._normalise(json.loads(json.dumps(c)), ZP) if False else c for c in case_defs()[name](rng)]
    tie = TIES_K if name == 'K_tied' else []
    wd = os.path.join(args.work, '%s_%02d' % (name, seed))
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    psf = moffat_psf()
    fits.PrimaryHDU(psf).writeto(os.path.join(wd, 'psf.fits'))
    fits.PrimaryHDU(np.zeros((NY, NX), np.float32)).writeto(os.path.join(wd, 'data.fits'))
    tcomps = json.loads(json.dumps(truth))
    if name == 'K_tied':                                                       # the truth must satisfy the ties
        tcomps[1]['x'], tcomps[1]['y'] = tcomps[0]['x'], tcomps[0]['y']
        tcomps[1]['mag'] = tcomps[0]['mag'] - 2.5 * math.log10(2.5)
    tcfg = dict(components=tcomps, zp=ZP, exptime=1.0, sky='fixed', sky_value=SKY, plate_scale=[1.0, 1.0])
    mask = 'none'
    all_truth = tcomps
    if name == 'L_masked_neighbour':                                           # bright un-fitted neighbour, masked out by the mask file given to both programs
        nb = dict(kind='sersic', x=51 + 24.0, y=51 - 6.0, mag=16.0, re=5.0, n=1.0, q=0.8, pa=40.0)
        all_truth = tcomps + [nb]
    open(os.path.join(wd, 'truth.feedme'), 'w').write(GI.write_feedme(dict(tcfg, components=all_truth), image='data.fits', output='out.fits', psf='psf.fits', shape=(NY, NX), mode=1, plate_scale=(1.0, 1.0)))
    t, r = run_galfit(args.galfit, args.ld, 'truth.feedme', wd)
    out = os.path.join(wd, 'out.fits')
    if not os.path.exists(out):
        return dict(case=name, seed=seed, error='GALFIT render failed: ' + r.stdout[-200:])
    model = fits.getdata(out).astype(float)
    noisy = model + rng.normal(0, SIGMA, model.shape)
    fits.PrimaryHDU(noisy.astype(np.float32)).writeto(os.path.join(wd, 'data.fits'), overwrite=True)
    fits.PrimaryHDU(np.full((NY, NX), SIGMA, np.float32)).writeto(os.path.join(wd, 'sigma.fits'))
    if name == 'L_masked_neighbour':
        yy, xx = np.mgrid[:NY, :NX]
        m = ((xx - (nb['x'] - 1)) ** 2 + (yy - (nb['y'] - 1)) ** 2 < 17.0 ** 2).astype(np.float32)
        fits.PrimaryHDU(m).writeto(os.path.join(wd, 'mask.fits'))
        mask = 'mask.fits'
    start = [_perturb(c, rng) for c in tcomps]
    if name == 'K_tied':
        start[1]['x'], start[1]['y'] = start[0]['x'], start[0]['y']
        start[1]['mag'] = start[0]['mag'] - 2.5 * math.log10(2.5)
    csky = float(np.median(noisy[:12, :12]))
    txt, cons = _fit_feed(start, tie, csky, mask)
    open(os.path.join(wd, 'start.feedme'), 'w').write(txt)
    if cons:
        open(os.path.join(wd, 'cons.txt'), 'w').write(cons)
    res = dict(case=name, seed=seed, truth=tcomps)
    gwd = os.path.join(wd, 'g'); mwd = os.path.join(wd, 'm')
    for d in (gwd, mwd):
        shutil.copytree(wd, d, ignore=shutil.ignore_patterns('g', 'm'))
    t, r = run_galfit(args.galfit, args.ld, 'start.feedme', gwd)
    res['galfit_time'] = t
    fs = sorted(f for f in os.listdir(gwd) if re.match(r'galfit\.\d+$', f))
    gtxt = None
    if fs:
        gtxt = open(os.path.join(gwd, fs[-1])).read()
        res['galfit'] = _flatten(_result_comps(gtxt, gwd))
        mm = re.search(r'Chi\^2/nu\s*=\s*([0-9.eE+-]+)', gtxt)
        res['galfit_chi2nu'] = float(mm.group(1)) if mm else None
    t = time.time()
    r = subprocess.run([sys.executable, MULTIFIT, '-', '--config', os.path.join(mwd, 'start.feedme'), '--work', os.path.join(mwd, 'w'), '--export-feedme', os.path.join(mwd, 'w', 'gf'), '--max-nfev', '600'],
                       capture_output=True, text=True, cwd=mwd)
    res['multifit_time'] = time.time() - t
    ex = os.path.join(mwd, 'w', 'gf', 'galfit.feedme')
    mtxt = None
    if r.returncode == 0 and os.path.exists(ex):
        mtxt = open(ex).read()
        res['multifit'] = _flatten(_result_comps(mtxt, mwd))
        mm = re.search(r'chi2/dof = ([0-9.]+)', r.stdout)
        res['multifit_chi2nu'] = float(mm.group(1)) if mm else None
    else:
        res['multifit_err'] = (r.stderr or r.stdout)[-300:]

    def arbiter(result_txt, tag):
        ed = os.path.join(wd, 'e_' + tag)
        os.makedirs(ed, exist_ok=True)
        for f in ('psf.fits', 'data.fits', 'sigma.fits', 'mask.fits'):
            if os.path.exists(os.path.join(wd, f)):
                shutil.copy(os.path.join(wd, f), ed)
        t_ = result_txt
        for key, val in (('A', 'data.fits'), ('B', 'out.fits'), ('C', 'none'), ('D', 'psf.fits'), ('F', 'none'), ('G', 'none'), ('P', '1')):
            t_ = re.sub(r'(?m)^%s\).*$' % key, '%s) %s' % (key, val), t_)
        open(os.path.join(ed, 'e.feedme'), 'w').write(t_)
        run_galfit(args.galfit, args.ld, 'e.feedme', ed)
        o = os.path.join(ed, 'out.fits')
        if not os.path.exists(o):
            return None
        good = np.ones((NY, NX), bool)
        if os.path.exists(os.path.join(ed, 'mask.fits')):
            good = fits.getdata(os.path.join(ed, 'mask.fits')) == 0
        return float((((noisy - fits.getdata(o).astype(float)) / SIGMA) ** 2)[good].sum())
    try:
        res['chi2_truth'] = arbiter(open(os.path.join(wd, 'truth.feedme')).read().replace('1) 12.0', '1) 12.0') if False else GI.write_feedme(dict(tcfg, components=tcomps, sky_value=SKY), image='data.fits', output='out.fits', psf='psf.fits', shape=(NY, NX), mode=1, plate_scale=(1.0, 1.0)), 't')
        if gtxt:
            res['chi2_galfit'] = arbiter(gtxt, 'g')
        if mtxt:
            res['chi2_multifit'] = arbiter(mtxt, 'm')
    except Exception as e:
        res['chi2_error'] = str(e)
    return res


def _err_table(results, who):
    acc = {}
    for r in results:
        if r.get(who) is None or 'truth' not in r:
            continue
        truth = _flatten(dict(components=r['truth']))
        for key, v in r[who].items():
            if key in truth:
                acc.setdefault(key, []).append(_num_cmp(key.split('.', 1)[1], v, truth[key]))
    return {k: dict(median=float(np.median(v)), rms=float(np.sqrt(np.mean(np.square(v)))), p90_abs=float(np.percentile(np.abs(v), 90))) for k, v in acc.items()}


def h2h(args):
    names = args.cases.split(',') if args.cases else list(case_defs())
    tasks = [(n, s, args) for n in names for s in range(args.n)]
    t0 = time.time()
    with ThreadPoolExecutor(args.workers) as ex:
        results = list(ex.map(one, tasks))
    summ = {}
    for n in names:
        rs = [r for r in results if r['case'] == n and 'truth' in r]
        both = [r for r in rs if r.get('galfit') and r.get('multifit')]
        s = dict(n=len(rs), galfit_ok=sum('galfit' in r for r in rs), multifit_ok=sum('multifit' in r for r in rs), both=len(both),
                 galfit_time_median=float(np.median([r['galfit_time'] for r in rs])) if rs else None, multifit_time_median=float(np.median([r['multifit_time'] for r in rs])) if rs else None)
        s['galfit_vs_truth'] = _err_table(both, 'galfit')
        s['multifit_vs_truth'] = _err_table(both, 'multifit')
        d = {}
        for r in both:
            for key, v in r['multifit'].items():
                if key in r['galfit']:
                    d.setdefault(key, []).append(_num_cmp(key.split('.', 1)[1], v, r['galfit'][key]))
        s['multifit_minus_galfit'] = {k: dict(median=float(np.median(v)), p90_abs=float(np.percentile(np.abs(v), 90)), max_abs=float(np.max(np.abs(v)))) for k, v in d.items()}
        c2 = [(r['chi2_galfit'], r['chi2_multifit'], r['chi2_truth']) for r in both if r.get('chi2_galfit') and r.get('chi2_multifit') and r.get('chi2_truth')]
        if c2:
            dd = np.array([c[1] - c[0] for c in c2])
            s['chi2_arbiter'] = dict(n=len(c2), multifit_minus_galfit_median=float(np.median(dd)), min=float(dd.min()), max=float(dd.max()), multifit_lower_frac=float(np.mean(dd < -0.5)), galfit_lower_frac=float(np.mean(dd > 0.5)),
                                     within_half_frac=float(np.mean(np.abs(dd) <= 0.5)), galfit_minus_truth_median=float(np.median([c[0] - c[2] for c in c2])), multifit_minus_truth_median=float(np.median([c[1] - c[2] for c in c2])))
        summ[n] = s
    json.dump(dict(n_per_case=args.n, sigma=SIGMA, elapsed_s=time.time() - t0, summary=summ, results=results), open(args.out, 'w'), default=lambda o: None)
    for n, s in summ.items():
        print(n, 'n', s['n'], 'galfit ok', s['galfit_ok'], 'multifit ok', s['multifit_ok'], 'times %.1f/%.1f s' % (s['galfit_time_median'] or 0, s['multifit_time_median'] or 0))
        if 'chi2_arbiter' in s:
            a = s['chi2_arbiter']
            print('   chi2 mf-galfit median %.3f [%.3f, %.3f]; mf lower %.2f, galfit lower %.2f, within 0.5: %.2f; (galfit-truth %.1f, mf-truth %.1f)' % (a['multifit_minus_galfit_median'], a['min'], a['max'], a['multifit_lower_frac'], a['galfit_lower_frac'], a['within_half_frac'], a['galfit_minus_truth_median'], a['multifit_minus_truth_median']))
        for k, v in s['multifit_minus_galfit'].items():
            print('   mf-galfit %-10s median %9.4f p90 %8.4f max %8.4f' % (k, v['median'], v['p90_abs'], v['max_abs']))
    return 0


def recovery_defs():
    D = dict(case_defs())
    cen = lambda rng: (51 + rng.uniform(-2, 2), 51 + rng.uniform(-2, 2))
    def mk(kind, **k):
        def f(rng):
            x, y = cen(rng)
            return [dict(dict(kind=kind, x=x, y=y, q=rng.uniform(0.5, 0.9), pa=rng.uniform(0, 180)), **k)]
        return f
    D['M_nuker'] = lambda rng: [dict(kind='nuker', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), ib=rng.uniform(15, 25), rb=rng.uniform(7, 10), alpha=1.5, beta=rng.uniform(1.6, 2.2), gamma=rng.uniform(0.2, 0.7),
                                  q=rng.uniform(0.6, 0.9), pa=rng.uniform(0, 180), fixed=['alpha'])]
    D['N_king'] = lambda rng: [dict(kind='king', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), i0=rng.uniform(20, 35), rc=rng.uniform(3, 5), rt=rng.uniform(25, 35), alpha=2.0, q=rng.uniform(0.6, 0.9), pa=rng.uniform(0, 180), fixed=['alpha'])]
    D['O_edgedisk'] = lambda rng: [dict(kind='edgedisk', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), i0=rng.uniform(8, 14), hs=rng.uniform(2.5, 4.0), rs=rng.uniform(12, 18), pa=rng.uniform(0, 180))]
    D['P_brokenexp'] = lambda rng: [dict(kind='brokenexp', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), i0=rng.uniform(6, 10), h1=rng.uniform(9, 12), h2=rng.uniform(4, 6), rbreak=rng.uniform(14, 18), alpha=0.5,
                                       q=rng.uniform(0.5, 0.9), pa=rng.uniform(0, 180))]
    D['Q_ferrer'] = lambda rng: [dict(kind='ferrer', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), i0=rng.uniform(25, 35), rout=rng.uniform(25, 32), alpha=3.0, beta=rng.uniform(0.8, 1.4), q=rng.uniform(0.6, 0.9), pa=rng.uniform(0, 180), fixed=['alpha'])]
    D['R_inner_truncation'] = lambda rng: [dict(kind='sersic', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), norm='re', i0=rng.uniform(5, 9), re=rng.uniform(8, 12), n=rng.uniform(1.0, 2.0), q=rng.uniform(0.6, 0.9), pa=rng.uniform(0, 180),
                                            trunc_in=[1], fixed=['n']), dict(kind='trunc', rbreak=rng.uniform(12, 16), dsoft=6.0, fixed=['dsoft'])]
    D['S_rot_log_free'] = lambda rng: [dict(kind='sersic', x=51 + rng.uniform(-2, 2), y=51 + rng.uniform(-2, 2), mag=rng.uniform(16.3, 17.2), re=rng.uniform(6, 9), n=rng.uniform(1.2, 2.5), q=rng.uniform(0.5, 0.9), pa=rng.uniform(0, 180),
                                         rot_func='log', rot_in=0.0, rot_out=rng.uniform(22, 30), rot_theta=rng.uniform(50, 90), rot_ws=3.0, rot_incl=0.0, rot_pa=0.0, fixed=['rot_in', 'rot_ws', 'rot_incl', 'rot_pa'])]
    return D


def _rec_one(task):
    name, seed, args = task
    from ogfkit import psfmodel as pm
    rng = np.random.default_rng(104729 * (ord(name[0]) - 60) + seed)
    truth = recovery_defs()[name](rng)
    tie = TIES_K if name == 'K_tied' else []
    if name == 'K_tied':
        truth[1]['x'], truth[1]['y'] = truth[0]['x'], truth[0]['y']
        truth[1]['mag'] = truth[0]['mag'] - 2.5 * math.log10(2.5)
    psf = moffat_psf().astype(float)
    cs = [MF._normalise(json.loads(json.dumps(c)), ZP) for c in truth]
    for c in cs:
        if 'x' in c:
            c['x'] -= 1.0; c['y'] -= 1.0
    model = MF.render_model(cs, (NY, NX), psf=psf, sky=SKY)
    nb_mask = None
    if name == 'L_masked_neighbour':
        nb = MF._normalise(dict(kind='sersic', x=24.0 + 50.0, y=44.0, mag=16.0, re=5.0, n=1.0, q=0.8, pa=130.0), ZP)
        model = model + MF.render_model([nb], (NY, NX), psf=psf, sky=0.0)
        yy, xx = np.mgrid[:NY, :NX]
        nb_mask = (xx - nb['x']) ** 2 + (yy - nb['y']) ** 2 < 17.0 ** 2
    data = model + rng.normal(0, SIGMA, model.shape)
    start = []
    for c, c0 in zip(truth, cs):
        p = _perturb(dict(c0, fixed=c.get('fixed', [])), rng)
        p.pop('mag', None) if 'mag' not in c else None
        start.append(p)
    if name == 'K_tied':
        start[1]['x'], start[1]['y'] = start[0]['x'], start[0]['y']
    for c in start:
        if 'x' in c:
            c['bounds'] = dict(c.get('bounds', {}))
        c['bounds'] = dict(_bounds_for(c), **c.get('bounds', {}))
    t0 = time.time()
    res = MF.fit_multistart(data, start, restarts=2, psf=psf, rms=SIGMA, mask=nb_mask, sky='const', zp=ZP, tie=[list(t) for t in tie], max_nfev=600)
    out = dict(case=name, seed=seed, truth=cs, time=time.time() - t0, chi2_red=res['chi2_red'], flags=res['flags'], converged=bool(res['converged']), errors={})
    for i, (c, f) in enumerate(zip(cs, res['components'])):
        for k in MF.param_names(c):
            if k in ('flux',) and 'mag' in f:
                out['errors']['c%d.mag' % i] = f['mag'] - (ZP - 2.5 * math.log10(c['flux']))
            elif k in f and k in c and k not in c.get('fixed', ()):
                v, tv = f[k], c[k]
                if re.match(r'^f\d+p$', k):
                    per = 360.0 / int(k[1:-1]); out['errors']['c%d.%s' % (i, k)] = (v - tv + per / 2) % per - per / 2
                else:
                    out['errors']['c%d.%s' % (i, k)] = _num_cmp(k, v, tv)
    out['errors']['sky'] = res['sky'] - SKY
    return out


def recovery(args):
    names = args.cases.split(',') if args.cases else list(recovery_defs())
    tasks = [(n, s, args) for n in names for s in range(args.n)]
    with ThreadPoolExecutor(args.workers) as ex:
        results = list(ex.map(_rec_one, tasks))
    summ = {}
    for n in names:
        rs = [r for r in results if r['case'] == n]
        acc = {}
        for r in rs:
            for k, v in r['errors'].items():
                acc.setdefault(k, []).append(v)
        summ[n] = dict(n=len(rs), chi2_red_median=float(np.median([r['chi2_red'] for r in rs])), converged=sum(r['converged'] for r in rs), time_median=float(np.median([r['time'] for r in rs])),
                       errors={k: dict(median=float(np.median(v)), rms=float(np.sqrt(np.mean(np.square(v)))), p90_abs=float(np.percentile(np.abs(v), 90))) for k, v in acc.items()})
    json.dump(dict(n_per_case=args.n, sigma=SIGMA, summary=summ, results=[{k: v for k, v in r.items() if k != 'truth'} for r in results]), open(args.out, 'w'), indent=1, default=lambda o: None)
    for n, s in summ.items():
        print('%-20s n %d chi2_red %.3f converged %d/%d  %.1f s' % (n, s['n'], s['chi2_red_median'], s['converged'], s['n'], s['time_median']))
        print('    ' + '  '.join('%s %.3g' % (k, v['rms']) for k, v in s['errors'].items()))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['parity', 'h2h', 'recovery'])
    ap.add_argument('out')
    ap.add_argument('--n', type=int, default=6)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--galfit', default=os.environ.get('GALFIT_BIN') or shutil.which('galfit') or '')
    ap.add_argument('--ld', default=os.environ.get('GALFIT_LD') or '')
    ap.add_argument('--work', default='/tmp/ogf_galfit_adv')
    ap.add_argument('--cases', default='')
    args = ap.parse_args()
    args.work = os.path.abspath(args.work)
    if args.mode != 'recovery' and (not args.galfit or not os.path.exists(args.galfit)):
        sys.exit('no GALFIT binary (--galfit / GALFIT_BIN)')
    return globals()[args.mode](args)


if __name__ == '__main__':
    sys.exit(main())
