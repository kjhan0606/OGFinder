#!/usr/bin/env python3
"""Head-to-head test of multifit and GALFIT on synthetic images with known truth.

    python galfit_compare.py OUT.json [--n 24] [--workers 8] [--galfit /path/galfit] [--ld DIR_WITH_libncurses5] [--cases A,B,C,D,E]

Images are RENDERED BY GALFIT itself (feedme P=1, so the generator is not ogfkit's renderer), noise is constant-sigma Gaussian, the same sigma image (C) and PSF (D) are given to
both programs, and both fit the same start feedme (with the same constraints file), so the two fits see identical data, model family, start values, bounds and ties.
multifit is run through its own CLI (plugins/multifit/multifit.py --config FEEDME); GALFIT's result file galfit.NN is parsed with ogfkit.galfitio.
Cases: A single Sersic (n 0.8-4); B devauc+expdisk with tied centre (constraints file); C psf + Sersic host; D Sersic with a Sersic neighbour; E Sersic + free sky gradient."""
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
from ogfkit import galfitio as GI  # noqa: E402

ZP = 25.0
SIGMA = 2.0
SKY = 50.0
NY, NX = 101, 101
MULTIFIT = os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py')


def moffat_psf(fwhm=3.2, beta=3.0, size=31):
    a = fwhm / (2 * math.sqrt(2 ** (1 / beta) - 1))
    y, x = np.mgrid[:size, :size] - size // 2
    p = (1 + (x ** 2 + y ** 2) / a ** 2) ** (-beta)
    return (p / p.sum()).astype(np.float32)


def feed_text(objs, sky, grad, data='data.fits', mode=1, free_sky_grad=False, constraints='none', region=None):
    region = region or (1, NX, 1, NY)
    t = ['A) %s' % data, 'B) out.fits', 'C) %s' % ('sigma.fits' if mode == 0 else 'none'), 'D) psf.fits', 'E) 1', 'F) none', 'G) %s' % constraints, 'H) %d %d %d %d' % region, 'I) 61 61',
         'J) %.3f' % ZP, 'K) 0.1 0.1', 'O) regular', 'P) %d' % mode, '']
    for o in objs:
        t.append(' 0) %s' % o['type'])
        t.append(' 1) %.4f %.4f %d %d' % (o['x'], o['y'], o.get('fx', 1), o.get('fy', 1)))
        t.append(' 3) %.4f 1' % o['mag'])
        if o['type'] != 'psf':
            t.append(' 4) %.4f 1' % o['re'])
            if o['type'] == 'sersic':
                t.append(' 5) %.4f 1' % o['n'])
            t.append(' 9) %.4f 1' % o['q'])
            t.append('10) %.4f 1' % o['pa'])
        t.append(' Z) 0')
    t += [' 0) sky', ' 1) %.4f 1' % sky, ' 2) %.6f %d' % (grad[0], 1 if free_sky_grad else 0), ' 3) %.6f %d' % (grad[1], 1 if free_sky_grad else 0), ' Z) 0']
    return '\n'.join(t) + '\n'


def draw_case(kind, rng):
    cx, cy = 51 + rng.uniform(-3, 3), 51 + rng.uniform(-3, 3)
    pa = rng.uniform(-90, 90)
    sers = lambda **k: dict(type='sersic', x=cx, y=cy, mag=rng.uniform(16.5, 18.0), re=rng.uniform(3.5, 9.0), n=rng.uniform(0.8, 4.0), q=rng.uniform(0.4, 0.95), pa=pa, **k)
    grad = (0.0, 0.0)
    cons = None
    if kind == 'A' or kind == 'E':
        objs = [sers()]
        if kind == 'E':
            grad = (rng.uniform(-0.04, 0.04), rng.uniform(-0.04, 0.04))
    elif kind == 'B':
        objs = [dict(type='devauc', x=cx, y=cy, mag=rng.uniform(17.8, 18.8), re=rng.uniform(1.5, 3.0), q=rng.uniform(0.7, 0.95), pa=pa + rng.uniform(-20, 20)),
                dict(type='expdisk', x=cx, y=cy, mag=rng.uniform(16.8, 17.6), re=rng.uniform(3.0, 5.0), q=rng.uniform(0.4, 0.8), pa=pa)]
        cons = '1_2 x offset\n1_2 y offset\n'
    elif kind == 'C':
        objs = [dict(type='sersic', x=cx, y=cy, mag=rng.uniform(17.0, 17.8), re=rng.uniform(4.0, 8.0), n=rng.uniform(1.0, 3.0), q=rng.uniform(0.6, 0.95), pa=pa),
                dict(type='psf', x=cx, y=cy, mag=rng.uniform(18.5, 19.8))]
        cons = '1_2 x offset\n1_2 y offset\n'
    elif kind == 'D':
        d = rng.uniform(2.5, 3.5) * 6.0
        th = rng.uniform(0, 2 * math.pi)
        objs = [sers(), dict(type='sersic', x=cx + d * math.cos(th), y=cy + d * math.sin(th), mag=rng.uniform(17.4, 18.4), re=rng.uniform(3.0, 5.0), n=rng.uniform(0.8, 2.5), q=rng.uniform(0.5, 0.95),
                             pa=rng.uniform(-90, 90))]
        objs[0]['re'] = rng.uniform(3.5, 6.0)
    else:
        raise ValueError(kind)
    return objs, grad, cons


def perturb(objs, rng, kind):
    out = []
    for o in objs:
        p = dict(o)
        p['x'] += rng.uniform(-1.0, 1.0); p['y'] += rng.uniform(-1.0, 1.0)
        p['mag'] += rng.uniform(-0.4, 0.4)
        if o['type'] != 'psf':
            p['re'] *= rng.uniform(0.7, 1.3)
            p['q'] = float(np.clip(o['q'] + rng.uniform(-0.2, 0.2), 0.2, 1.0))
            p['pa'] = o['pa'] + rng.uniform(-30, 30)
            if o['type'] == 'sersic':
                p['n'] = float(np.clip(o['n'] * rng.uniform(0.6, 1.5), 0.5, 6))
        out.append(p)
    if kind in ('B', 'C'):                                                    # tied centres must start equal
        out[1]['x'], out[1]['y'] = out[0]['x'], out[0]['y']
    return out


def run_galfit(galfit, ld, feed, cwd, timeout=300):
    env = dict(os.environ)
    if ld:
        env['LD_LIBRARY_PATH'] = ld + ':' + env.get('LD_LIBRARY_PATH', '')
    t = time.time()
    r = subprocess.run([galfit, feed], cwd=cwd, env=env, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=timeout)
    return time.time() - t, r


def galfit_result(cwd):
    fs = sorted(f for f in os.listdir(cwd) if re.match(r'galfit\.\d+$', f))
    if not fs:
        return None
    path = os.path.join(cwd, fs[-1])
    txt = open(path).read()
    cfg = GI.parse_feedme(txt, strict=False, base_dir=cwd)
    chi = re.search(r'Chi\^2/nu\s*=\s*([0-9.eE+-]+)', txt)
    flag = 'Chi^2' in txt
    return cfg, (float(chi.group(1)) if chi else None), ('*' in txt.split('# Chi^2')[0] if False else None)


def eval_chi2(args, wd, feed_text_result, tag):
    """chi2 of a result feedme (any program's), evaluated with GALFIT's own renderer (P=1) against the same data and constant sigma -> the neutral arbiter."""
    ed = os.path.join(wd, 'e_' + tag)
    os.makedirs(ed, exist_ok=True)
    for f in ('psf.fits', 'data.fits', 'sigma.fits'):
        shutil.copy(os.path.join(wd, f), ed)
    txt = re.sub(r'(?m)^A\).*$', 'A) data.fits', feed_text_result)
    txt = re.sub(r'(?m)^B\).*$', 'B) out.fits', txt)
    txt = re.sub(r'(?m)^C\).*$', 'C) none', txt)
    txt = re.sub(r'(?m)^D\).*$', 'D) psf.fits', txt)
    txt = re.sub(r'(?m)^G\).*$', 'G) none', txt)
    txt = re.sub(r'(?m)^P\).*$', 'P) 1', txt)
    open(os.path.join(ed, 'e.feedme'), 'w').write(txt)
    run_galfit(args.galfit, args.ld, 'e.feedme', ed)
    o = os.path.join(ed, 'out.fits')
    if not os.path.exists(o):
        return None
    data = fits.getdata(os.path.join(ed, 'data.fits')).astype(float)
    mod = fits.getdata(o).astype(float)
    return float(np.sum(((data - mod) / SIGMA) ** 2))


def one(task):
    kind, seed, args = task
    rng = np.random.default_rng(1000 * ord(kind) + seed)
    objs, grad, cons = draw_case(kind, rng)
    wd = os.path.join(args.work, '%s%02d' % (kind, seed))
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    fits.PrimaryHDU(moffat_psf()).writeto(os.path.join(wd, 'psf.fits'))
    fits.PrimaryHDU(np.zeros((NY, NX), np.float32)).writeto(os.path.join(wd, 'data.fits'))
    open(os.path.join(wd, 'truth.feedme'), 'w').write(feed_text(objs, SKY, grad, mode=1))
    t, r = run_galfit(args.galfit, args.ld, 'truth.feedme', wd)
    out = os.path.join(wd, 'out.fits')
    if not os.path.exists(out):
        return dict(kind=kind, seed=seed, error='galfit model render failed: ' + r.stdout[-200:])
    model = fits.getdata(out).astype(float)
    noisy = model + rng.normal(0, SIGMA, model.shape)
    fits.PrimaryHDU(noisy.astype(np.float32)).writeto(os.path.join(wd, 'data.fits'), overwrite=True)
    fits.PrimaryHDU(np.full((NY, NX), SIGMA, np.float32)).writeto(os.path.join(wd, 'sigma.fits'))
    start = perturb(objs, rng, kind)
    csky = float(np.median(noisy[:12, :12]))
    if cons:
        open(os.path.join(wd, 'cons.txt'), 'w').write(cons)
    st = feed_text(start, csky, (0.0, 0.0), mode=0, free_sky_grad=(kind == 'E'), constraints='cons.txt' if cons else 'none')
    open(os.path.join(wd, 'start.feedme'), 'w').write(st)
    res = dict(kind=kind, seed=seed, truth=objs, truth_sky=SKY, truth_grad=list(grad))
    # GALFIT
    gwd = os.path.join(wd, 'g')
    shutil.copytree(wd, gwd, ignore=shutil.ignore_patterns('g', 'm'))
    t, r = run_galfit(args.galfit, args.ld, 'start.feedme', gwd)
    gr = galfit_result(gwd)
    res['galfit_time'] = t
    if gr is None:
        res['galfit'] = None
    else:
        cfg, chi, _ = gr
        res['galfit'] = dict(comps=[dict(type=c['galfit_type'], x=c['x'], y=c['y'], mag=c['mag'], re=c.get('re'), n=c.get('n'), q=c.get('q'), pa=c.get('pa', 90) - 90, galfit_re=None) for c in cfg['components']],
                             sky=cfg.get('sky_value_galfit'), grad=cfg.get('sky_grad'), chi2nu=chi)
        for c, rc in zip(res['galfit']['comps'], cfg['components']):                 # report GALFIT-native units (R_s for expdisk)
            if rc['galfit_type'] == 'expdisk':
                c['re'] = rc['re'] / GI.RS_TO_RE
    # multifit
    mwd = os.path.join(wd, 'm')
    shutil.copytree(wd, mwd, ignore=shutil.ignore_patterns('g', 'm'))
    t = time.time()
    r = subprocess.run([sys.executable, MULTIFIT, '-', '--config', os.path.join(mwd, 'start.feedme'), '--work', os.path.join(mwd, 'w'), '--export-feedme', os.path.join(mwd, 'w', 'gf'), '--max-nfev', '400'],
                       capture_output=True, text=True, cwd=mwd)
    res['multifit_time'] = time.time() - t
    ex = os.path.join(mwd, 'w', 'gf', 'galfit.feedme')
    if r.returncode != 0 or not os.path.exists(ex):
        res['multifit'] = None
        res['multifit_err'] = (r.stderr or r.stdout)[-300:]
    else:
        cfg = GI.parse_feedme(ex, strict=False, base_dir=mwd)
        m = re.search(r'chi2/dof = ([0-9.]+)', r.stdout)
        res['multifit'] = dict(comps=[dict(type=c['galfit_type'], x=c['x'], y=c['y'], mag=c['mag'], re=c.get('re'), n=c.get('n'), q=c.get('q'), pa=c.get('pa', 90) - 90) for c in cfg['components']],
                               sky=cfg.get('sky_value_galfit'), grad=cfg.get('sky_grad'), chi2nu=float(m.group(1)) if m else None)
        for c, rc in zip(res['multifit']['comps'], cfg['components']):
            if rc['galfit_type'] == 'expdisk':
                c['re'] = rc['re'] / GI.RS_TO_RE
    try:
        res['dof_pixels'] = NX * NY
        res['chi2_truth'] = eval_chi2(args, wd, open(os.path.join(wd, 'truth.feedme')).read(), 'truth')
        if res.get('galfit'):
            fs = sorted(f for f in os.listdir(gwd) if re.match(r'galfit\.\d+$', f))
            res['galfit']['chi2'] = eval_chi2(args, wd, open(os.path.join(gwd, fs[-1])).read(), 'g')
        if res.get('multifit'):
            res['multifit']['chi2'] = eval_chi2(args, wd, open(ex).read(), 'm')
    except Exception as e:
        res['chi2_error'] = str(e)
    return res


def wrap_pa(d):
    return (d + 90.0) % 180.0 - 90.0


def param_errors(fit, truth):
    out = []
    for f, t in zip(fit['comps'], truth):
        d = dict(x=f['x'] - t['x'], y=f['y'] - t['y'], mag=f['mag'] - t['mag'])
        if t['type'] != 'psf':
            d['re_rel'] = f['re'] / t['re'] - 1
            d['q'] = f['q'] - t['q']
            d['pa'] = wrap_pa(f['pa'] - t['pa'])
            if t['type'] == 'sersic':
                d['n_rel'] = f['n'] / t['n'] - 1
        out.append(d)
    return out


def summarise(results):
    summ = {}
    for kind in sorted({r['kind'] for r in results}):
        rs = [r for r in results if r['kind'] == kind and 'truth' in r]
        s = dict(n=len(rs), galfit_ok=sum(1 for r in rs if r.get('galfit')), multifit_ok=sum(1 for r in rs if r.get('multifit')),
                 galfit_time_median=float(np.median([r['galfit_time'] for r in rs])), multifit_time_median=float(np.median([r['multifit_time'] for r in rs])))
        both = [r for r in rs if r.get('galfit') and r.get('multifit')]
        s['both_ok'] = len(both)
        for who in ('galfit', 'multifit'):
            acc = {}
            for r in both:
                for k, e in enumerate(param_errors(r[who], r['truth'])):
                    for p, v in e.items():
                        acc.setdefault('c%d.%s' % (k, p), []).append(v)
            s[who + '_vs_truth'] = {k: dict(median=float(np.median(v)), std=float(np.std(v)), p90_abs=float(np.percentile(np.abs(v), 90))) for k, v in acc.items()}
        acc = {}
        for r in both:
            for k, (g, m) in enumerate(zip(r['galfit']['comps'], r['multifit']['comps'])):
                d = dict(x=m['x'] - g['x'], y=m['y'] - g['y'], mag=m['mag'] - g['mag'])
                if g['type'] != 'psf':
                    d['re_rel'] = m['re'] / g['re'] - 1; d['q'] = m['q'] - g['q']; d['pa'] = wrap_pa(m['pa'] - g['pa'])
                    if g['type'] == 'sersic':
                        d['n_rel'] = m['n'] / g['n'] - 1
                for p, v in d.items():
                    acc.setdefault('c%d.%s' % (k, p), []).append(v)
        s['multifit_minus_galfit'] = {k: dict(median=float(np.median(v)), std=float(np.std(v)), p90_abs=float(np.percentile(np.abs(v), 90)), max_abs=float(np.max(np.abs(v)))) for k, v in acc.items()}
        c1 = [(r['galfit']['chi2nu'], r['multifit']['chi2nu']) for r in both if r['galfit'].get('chi2nu') and r['multifit'].get('chi2nu')]
        if c1:
            s['chi2nu_galfit_median'] = float(np.median([c[0] for c in c1])); s['chi2nu_multifit_median'] = float(np.median([c[1] for c in c1]))
            s['chi2nu_multifit_lower_or_equal_frac'] = float(np.mean([c[1] <= c[0] + 0.002 for c in c1]))
        c2 = [(r['galfit']['chi2'], r['multifit']['chi2'], r['chi2_truth']) for r in both if r['galfit'].get('chi2') and r['multifit'].get('chi2') and r.get('chi2_truth')]
        if c2:
            d = np.array([c[1] - c[0] for c in c2])
            s['chi2_arbiter'] = dict(n=len(c2), multifit_minus_galfit_median=float(np.median(d)), p10=float(np.percentile(d, 10)), p90=float(np.percentile(d, 90)), min=float(d.min()), max=float(d.max()),
                                     multifit_lower_frac=float(np.mean(d < -0.5)), galfit_lower_frac=float(np.mean(d > 0.5)), within_half_frac=float(np.mean(np.abs(d) <= 0.5)),
                                     galfit_minus_truth_median=float(np.median([c[0] - c[2] for c in c2])), multifit_minus_truth_median=float(np.median([c[1] - c[2] for c in c2])))
        summ[kind] = s
    return summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--n', type=int, default=24)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--galfit', default=os.environ.get('GALFIT_BIN') or shutil.which('galfit') or '')
    ap.add_argument('--ld', default=os.environ.get('GALFIT_LD') or '')
    ap.add_argument('--cases', default='A,B,C,D,E')
    ap.add_argument('--work', default='/tmp/ogf_galfit_cmp')
    args = ap.parse_args()
    if not args.galfit or not os.path.exists(args.galfit):
        sys.exit('no GALFIT binary (--galfit / GALFIT_BIN)')
    args.work = os.path.abspath(args.work)
    tasks = [(k, s, args) for k in args.cases.split(',') for s in range(args.n)]
    t0 = time.time()
    with ThreadPoolExecutor(args.workers) as ex:
        results = list(ex.map(one, tasks))
    rep = dict(n_per_case=args.n, sigma=SIGMA, zp=ZP, elapsed_s=time.time() - t0, summary=summarise(results), results=results)
    json.dump(rep, open(args.out, 'w'), default=lambda o: None)
    for k, s in rep['summary'].items():
        print(k, 'n', s['n'], 'galfit ok', s['galfit_ok'], 'multifit ok', s['multifit_ok'], 'times %.2f/%.2f s' % (s['galfit_time_median'], s['multifit_time_median']))
        if 'chi2_arbiter' in s:
            print('   chi2 (GALFIT-rendered) multifit-galfit: median %.3f p10 %.3f p90 %.3f min %.3f max %.3f; multifit lower (>0.5) %.2f, galfit lower %.2f, within 0.5: %.2f' % tuple(
                s['chi2_arbiter'][k] for k in ('multifit_minus_galfit_median', 'p10', 'p90', 'min', 'max', 'multifit_lower_frac', 'galfit_lower_frac', 'within_half_frac')))
        for p, v in s['multifit_minus_galfit'].items():
            print('   mf-galfit %-8s median %9.4f std %8.4f p90 %8.4f max %8.4f' % (p, v['median'], v['std'], v['p90_abs'], v['max_abs']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
