#!/usr/bin/env python3
"""Synthetic-truth validation of ogfkit/kinematics.py   (python3 run_kin_validation.py OUT.json [n_per_config])

Cubes: random inclined discs (arctan rotation curve, exponential flux, constant intrinsic dispersion) at several noise levels, with / without spatial seeing, per-spaxel
fits vs accretion binning, free vs fixed inclination.  Slits: the same discs through the centre at slit/major-axis angles psi = 0 and 30 deg.
Everything is compared with the generating parameters; the generator and the fitter share the disc model, which is stated in docs/spectra.md."""
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import kinematics as K, spectrasynth as SS  # noqa: E402


def stats(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return dict(n=int(v.size), median=float(np.median(v)), mean=float(np.mean(v)), std=float(np.std(v)), max_abs=float(np.max(np.abs(v)))) if v.size else dict(n=0)


def draw(rng):
    return dict(vc=rng.uniform(120, 300), rt=rng.uniform(2.0, 6.0), inc=rng.uniform(30, 75), pa=rng.uniform(0, 360), sigma0=rng.uniform(25, 70), rd=rng.uniform(3.0, 7.0),
                vsys=rng.uniform(-30, 30))


def cube_case(rng, seed, noise, seeing, bin_snr, inc_fixed):
    p = draw(rng)
    x0, y0 = 20 + rng.uniform(-1.5, 1.5), 20 + rng.uniform(-1.5, 1.5)
    w, c, t = SS.disc_cube(seed=seed, noise=noise, seeing_fwhm_pix=seeing, x0=x0, y0=y0, flux0=30.0, **p)
    m = K.fit_cube_line(c, w, t['lam_sys'], window_kms=600, snr_pix=3.0, bin_snr=bin_snr, inst_fwhm_A=2.0)
    vf = K.fit_velocity_field(m['vel'], m['vel_err'], inc_fixed=p['inc'] if inc_fixed else None)
    if not vf.get('ok'):
        return None
    sd = K.dispersion_summary(m['sigma'], m['sigma_err'], m['flux'])
    return dict(dpa=K.angle_diff(vf['pa'], t['pa']), dinc=vf['inc'] - t['inc'], dvc_rel=(vf['vc'] - t['vc']) / t['vc'], dvsini_rel=(vf['vsini'] - t['vsini']) / t['vsini'],
                drt_rel=(vf['rt'] - t['rt']) / t['rt'], dvsys=vf['vsys'] - t['vsys'], dcen=math.hypot(vf['x0'] - t['x0'], vf['y0'] - t['y0']), sigma_ratio=sd['mean'] / t['sigma0'],
                chi2r=vf['chi2r'], nspax=vf['n'], corr_vc_inc=vf['corr_vc_inc'])


def slit_case(rng, seed, noise, psi, seeing):
    p = draw(rng)
    w, d, t = SS.disc_slit(seed=seed, noise=noise, psi=psi, seeing_fwhm_pix=seeing, flux0=30.0, **p)
    prof = K.fit_slit_line(w, d, None, t['lam_sys'], window_kms=600, bin_snr=6.0, inst_fwhm_A=2.0)
    rc = K.fit_rotation_curve(prof['y'], prof['v'], prof['v_err'], inc_deg=p['inc'], slit_offset_deg=psi)
    if not rc.get('ok'):
        return None
    # the intrinsic curve along the slit: asymptote = vc sin(i) cos(psi)
    expect = t['vc'] * math.sin(math.radians(t['inc'])) * math.cos(math.radians(psi))
    sd = K.dispersion_summary(prof['sigma'], prof['sigma_err'], prof['flux'])
    # r_t along the slit in pixels: the deprojected radius along a slit at angle psi is s*sqrt(cos^2 psi + sin^2 psi / cos^2 i); report the raw ratio
    return dict(dvobs_rel=(rc['vflat_obs'] - expect) / expect, dvc_rel=(rc['vc'] - t['vc'] * math.cos(math.radians(psi)) / math.cos(math.radians(psi))) / t['vc'] if rc.get('vc') else float('nan'),
                dvsys=rc['vsys'] - t['vsys'], dy0=rc['y0'] - t['y0'], drt_rel=(rc['rt'] - t['rt']) / t['rt'] if psi == 0 else float('nan'), sigma_ratio=sd['mean'] / t['sigma0'], chi2r=rc['chi2r'],
                nbin=len(prof['y']))


def summarize(rows):
    rows = [r for r in rows if r]
    keys = [k for k in rows[0] if k not in ('nspax', 'nbin')] if rows else []
    d = dict(n_ok=len(rows))
    for k in keys:
        d[k] = stats([r[k] for r in rows])
    return d


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    R = {}
    t0 = time.time()
    cfgs = [('cube_highSN_free_inc', dict(noise=0.05, seeing=0.0, bin_snr=0.0, inc_fixed=False)),
            ('cube_highSN_fixed_inc', dict(noise=0.05, seeing=0.0, bin_snr=0.0, inc_fixed=True)),
            ('cube_lowSN_unbinned', dict(noise=0.5, seeing=0.0, bin_snr=0.0, inc_fixed=False)),
            ('cube_lowSN_binned', dict(noise=0.5, seeing=0.0, bin_snr=8.0, inc_fixed=False)),
            ('cube_seeing_2.5px', dict(noise=0.05, seeing=2.5, bin_snr=0.0, inc_fixed=False))]
    for name, kw in cfgs:
        rng = np.random.default_rng(100 + len(name))
        R[name] = summarize([cube_case(rng, 1000 + i, **kw) for i in range(n)])
        R[name]['config'] = kw
    for name, psi, noise, seeing in (('slit_psi0', 0.0, 0.05, 0.0), ('slit_psi30', 30.0, 0.05, 0.0), ('slit_psi0_noisy', 0.0, 0.5, 0.0), ('slit_psi0_seeing_2.5px', 0.0, 0.05, 2.5)):
        rng = np.random.default_rng(200 + len(name))
        R[name] = summarize([slit_case(rng, 2000 + i, noise, psi, seeing) for i in range(n)])
        R[name]['config'] = dict(psi=psi, noise=noise, seeing=seeing)
    R['elapsed_s'] = time.time() - t0
    json.dump(R, open(out, 'w'), indent=1)
    for k, v in R.items():
        if isinstance(v, dict) and 'n_ok' in v:
            print(k, v['n_ok'], {m: (round(v[m]['median'], 4), round(v[m]['std'], 4)) for m in v if m not in ('n_ok', 'config') and v[m].get('n')})


if __name__ == '__main__':
    main()
