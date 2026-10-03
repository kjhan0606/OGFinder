"""Regression cases: one or two per tool, each returns a dict of named metrics that run_regression.py compares with baselines.json.

Every case works on public data (datasets.py) and calls the same entry point as the plugin (CLI or ogfkit function)."""
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

from . import datasets as D

ROOT = D.ROOT
PY = os.environ.get('OGFINDER_PYTHON', sys.executable)
for p in (ROOT, os.path.join(ROOT, 'plugins', 'noisemodel', 'validation'), os.path.join(ROOT, 'plugins', 'psfex', 'validation'),
          os.path.join(ROOT, 'plugins', 'completeness', 'validation'), os.path.join(ROOT, 'plugins', 'stacking', 'validation')):
    if p not in sys.path:
        sys.path.insert(0, p)

CASES = {}


def case(tool, name):
    def deco(f):
        CASES[name] = dict(tool=tool, fn=f)
        return f
    return deco


def _run(args, **kw):
    p = subprocess.run([PY] + [str(a) for a in args], capture_output=True, text=True, timeout=kw.pop('timeout', 900), **kw)
    if p.returncode != 0:
        raise RuntimeError('%s failed: %s' % (args[0], p.stderr[-400:]))
    return p


def _catalog(img, work, thresh='3'):
    sex = os.path.join(ROOT, 'bin', 'ds9_sextract')
    if not os.path.exists(sex):
        raise D.Unavailable('bin/ds9_sextract missing')
    out = subprocess.run([sex, img, '--detect-thresh', thresh], capture_output=True, text=True, timeout=600, check=True).stdout
    path = os.path.join(work, 'cat.tsv')
    open(path, 'w').write(out)
    return path


def _read_tsv(path):
    L = [l for l in open(path).read().split('\n') if l and not l.startswith('#')]
    cols = L[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t'))) for l in L[1:]]


def _f(v, d=float('nan')):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


@case('extract', 'extract_hudf')
def extract_hudf(w):
    img = D.hudf_crop()
    cols, rows = _read_tsv(_catalog(img, w))
    mag = np.array([_f(r.get('MAG_AUTO')) for r in rows])
    mag = mag[np.isfinite(mag) & (mag < 90)]
    return dict(n_sources=len(rows), n_mag_lt_24=int((mag < 24).sum()), brightest_mag=float(mag.min()), median_mag=float(np.median(mag)))


@case('noisemodel', 'noise_blank_apertures')
def noise_blank(w):
    import photerr_validate as PV
    from astropy.io import fits
    img = np.nan_to_num(fits.getdata(D.hudf_crop()).astype(float))
    o = PV.blank_test(img, 'hudf_crop', radii=(3.0, 5.0), n=300)
    r3, r5 = o['rows']
    return dict(pull_std_local_r3=r3['pull_local_law']['std'], pull_std_local_r5=r5['pull_local_law']['std'], pull_std_naive_r3=r3['pull_naive']['std'],
                sigma_pix=o['sigma_pix'])


@case('noisemodel', 'photerr_injection')
def photerr_injection(w):
    import photerr_validate as PV
    from astropy.io import fits
    img = np.nan_to_num(fits.getdata(D.hudf_crop()).astype(float))[:900, :900]
    r = PV.inject_test(img, 'hudf_crop', 40, 5)
    return dict(pull_std_isolated=r['isolated/unified_full']['std'] if 'isolated/unified_full' in r else r['isolated/unified_no_contam']['std'],
                pull_std_naive_isolated=r['isolated/naive']['std'], frac_measured=r['frac_measured'])


@case('completeness', 'depth_hudf')
def depth_hudf(w):
    import depth_validate as DV
    from astropy.io import fits
    img = np.nan_to_num(fits.getdata(D.hudf_crop()).astype(float))
    mlim, P, f, rm = DV.depth_of(img, r=3.0, zp=25.94, nsig=5.0)
    ok = np.isfinite(mlim) & ~P['bad']
    return dict(depth_5sigma_r3_median=float(np.median(mlim[ok])), depth_p10=float(np.percentile(mlim[ok], 10)), depth_p90=float(np.percentile(mlim[ok], 90)),
                aperture_noise_factor=float(np.median(np.atleast_1d(f))))


@case('psfex', 'psf_heldout_residual')
def psf_heldout(w):
    import psf_real as PR
    from astropy.io import fits
    d = np.array(fits.getdata(D.hudf_crop()), float)[:1300, :1300]
    r = PR.one_field('HUDF crop', d, 3.0, 40.0, 30.0, 31)
    m = {k: v['median'] for k, v in r.get('models', {}).items()}
    out = dict(n_stars=r.get('n_used_all', 0))
    for k, v in m.items():
        out['resid_' + k] = float(v)
    return out


@case('stacking', 'stack_null_calibration')
def stack_null(w):
    from ogfkit import imageio, tsvio, stacking as st
    img = D.hudf_crop()
    cols, rows = _read_tsv(_catalog(img, w))
    data, _ = imageio.load_image(img)
    cat = np.array([[_f(r['X_IMAGE']) - 1, _f(r['Y_IMAGE']) - 1] for r in rows])
    cfg = st.StackConfig(half=24, method='mean', n_boot=100, ap_r=6.0, bkg='annulus', mask_cat_scale=1.0, mask_cat_min=4.0, min_valid=0.5)
    neigh = np.column_stack([cat, [max(4.0, _f(r.get('A_IMAGE'), 2) * _f(r.get('KRON_RADIUS'), 3)) for r in rows]])
    nx_, ny_ = st.null_positions(data, None, 2400, cfg, avoid=cat)
    call = st.make_cutouts(data, None, nx_, ny_, cfg, neigh=neigh)
    good = np.where(call['ok'])[0]
    zs = []
    N = 100
    for s_ in range(min(len(good) // N, 20)):
        ok = np.zeros(len(call['ok']), bool)
        ok[good[s_ * N:(s_ + 1) * N]] = True
        r = st.stack_cutouts(call['cube'], ok, cfg)
        zs.append(r['aper'] / r['aper_err'])
    zs = np.array(zs)
    return dict(n_sets=len(zs), z_std=float(zs.std(ddof=1)), z_mean=float(zs.mean()))


@case('isophote', 'isophote_bright_galaxy')
def isophote_gal(w):
    img = D.hudf_crop()
    cat = _catalog(img, w)
    cols, rows = _read_tsv(cat)
    ext = [r for r in rows if 'MAG_AUTO' in r and 20 < _f(r['MAG_AUTO']) < 23.5 and _f(r.get('A_IMAGE'), 0) > 5 and _f(r.get('B_IMAGE'), 0) / max(_f(r.get('A_IMAGE'), 1), 1e-6) < 0.8]
    ext.sort(key=lambda r: _f(r['MAG_AUTO']))
    nums = ','.join(r['NUMBER'] for r in ext[:3])
    tab = os.path.join(w, 'iso.tsv')
    p = _run([os.path.join(ROOT, 'plugins', 'isophote', 'isophote.py'), img, '--catalog', cat, '--numbers', nums, '--mag-zeropoint', '25.94', '--pixel-scale', '0.06',
              '--json-out', os.path.join(w, 'iso.json'), '--table-out', tab, '--n-workers', '1'])
    cols2, rows2 = _read_tsv_stdout(p.stdout)
    out = dict(n_fitted=len(rows2))
    numeric = [c for c in cols2 if c != 'NUMBER']
    for c in numeric[:12]:
        vals = np.array([_f(r[c]) for r in rows2])
        if np.isfinite(vals).any():
            out[c.lower()] = float(np.nanmedian(vals))
    return out


def _read_tsv_stdout(text):
    L = [l for l in text.split('\n') if l and not l.startswith('#')]
    cols = L[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t'))) for l in L[1:]]


@case('multifit', 'multifit_sersic')
def multifit_sersic(w):
    img = D.hudf_crop()
    cat = _catalog(img, w)
    cols, rows = _read_tsv(cat)
    ext = [r for r in rows if 21 < _f(r['MAG_AUTO']) < 24 and _f(r.get('A_IMAGE'), 0) > 4]
    ext.sort(key=lambda r: _f(r['MAG_AUTO']))
    nums = ','.join(r['NUMBER'] for r in ext[:6])
    p = _run([os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py'), img, '--catalog', cat, '--work', os.path.join(w, 'mf'), '--objects', nums, '--model', 'sersic',
              '--psf-fwhm', '2.8', '--mag-zeropoint', '25.0', '--neighbours', 'fit', '--n-workers', '1'], timeout=1500)
    cols2, rows2 = _read_tsv_stdout(p.stdout)
    sel = set(nums.split(','))
    rows2 = [r for r in rows2 if r['NUMBER'] in sel and np.isfinite(_f(r.get('GF_MAG')))]
    out = dict(n_fitted=len(rows2))
    mag0 = {r['NUMBER']: _f(r['MAG_AUTO']) for r in ext}
    out['median_abs_dmag_vs_mag_auto'] = float(np.median([abs(_f(r['GF_MAG']) - mag0[r['NUMBER']]) for r in rows2])) if rows2 else float('nan')
    for c in ('GF_RE', 'GF_N', 'GF_CHI2NU'):
        if c in cols2 and rows2:
            out['median_' + c.lower()] = float(np.median([_f(r[c]) for r in rows2]))
    return out


@case('moving', 'moving_sdss_run94')
def moving_sdss(w):
    D.sdss_run94()
    out = os.path.join(w, 'ska.json')
    _run([os.path.join(ROOT, 'moving', 'validation', 'sdss_known_asteroids.py'), out, '--camcol', 3, '--fraction', 0.99], timeout=1800)
    d = json.load(open(out))
    b, o = d['baseline'], d['orbit_stage_0.99']
    return dict(ceiling=d['ceiling'], baseline_true=b['n_true'], baseline_tracklets=b['n_tracklets'], orbit_true=o['n_true'], orbit_tracklets=o['n_tracklets'],
                orbit_precision_lb=o['precision_lower_bound'], orbit_recall=o['recall'])


@case('photoz_sed', 'photoz_closure_sdss')
def photoz_closure(w):
    h5 = os.path.join(ROOT, 'photo_z', 'data', 'sdss_specphoto.h5')
    if not os.path.isfile(h5):
        raise D.Unavailable('photo_z/data/sdss_specphoto.h5 missing')
    out = os.path.join(w, 'closure.json')
    env = dict(os.environ, OGF_EAZY_PYTHON='/nonexistent')
    try:
        _run([os.path.join(ROOT, 'plugins', 'photoz_sed', 'validation', 'closure_validate.py'), out], env=env, timeout=1200)
    except RuntimeError as e:
        if 'torch' in str(e) or 'No module' in str(e):
            raise D.Unavailable('torch / MDN checkpoint not available')
        raise
    d = json.load(open(out))
    flat = {}
    for tag, key in (('before', d['cal_to_test']['before']), ('after', d['cal_to_test']['after']), ('crossfit_after', d['crossfit_all']['after_oof'])):
        flat['ks_p_' + tag] = key['ks_p']
        flat['cov68_' + tag] = key['coverage']['0.68']
        flat['crps_' + tag] = key['crps']
    flat['x2_injected_cov68_after'] = d['injected_scale_x2']['after_oof']['coverage']['0.68']
    flat['x05_injected_cov68_after'] = d['injected_scale_x0.5']['after_oof']['coverage']['0.68']
    return flat


@case('lensmodel', 'lens_j0946_evidence')
def lens_j0946(w):
    sys.path.insert(0, os.path.join(ROOT, 'plugins', 'lensmodel', 'validation'))
    import lens_real_j0946 as LR
    return LR.regression_metrics()
