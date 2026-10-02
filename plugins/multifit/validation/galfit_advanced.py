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
