import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import multifit as mf, autodecomp as ad

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'multifit.py')
ZP = 25.0
FWHM = 3.0


def gal(comps, shape=(91, 91), rms=1.0, seed=0):
    psf = mf.gaussian_psf_array(FWHM)
    m = mf.render_model([mf._normalise(c, ZP) for c in comps], shape, psf=psf)
    return m + np.random.default_rng(seed).normal(0, rms, shape), psf


BD = [dict(kind='sersic', x=45, y=45, mag=16.9, re=3.0, n=2.5, q=0.85, pa=40), dict(kind='exp', x=45, y=45, mag=16.4, re=12.0, q=0.6, pa=30)]
SER = [dict(kind='sersic', x=45, y=45, mag=16.4, re=8.0, n=2.0, q=0.7, pa=60)]


def test_profile_formulas():
    # total flux of a Sersic / exponential from (Ie, re, n, q) against the rendered image
    n, re, q, f = 3.0, 6.0, 0.7, 5000.0
    img = mf.render_component(mf._normalise(dict(kind='sersic', x=100, y=100, flux=f, re=re, n=n, q=q, pa=0), 25.0), (201, 201), None)
    assert abs(img.sum() / f - 1) < 0.03
    ie = ad.sersic_ie_from_flux(f, re, n, q)
    assert abs(ad.sersic_total(ie, re, n, q) / f - 1) < 1e-9
    assert abs(ad.sersic_profile(np.array([re]), ie, re, n)[0] / ie - 1) < 1e-9
    assert abs(ad.exp_total(2.0, 5.0, 0.5) - 2 * math.pi * 25 * 0.5 * 2.0) < 1e-9
    assert abs(ad.bn(1.0) - 1.678) < 0.002 and abs(ad.bn(4.0) - 7.669) < 0.01
    assert abs(ad.circ_median_pa([179, 1, 0]) % 180 - 0.0) < 1.0 or abs(ad.circ_median_pa([179, 1, 0]) - 180) < 1.0
    assert ad.pa_diff(179, 1) == pytest.approx(2.0)


def synthetic_profile(bt, reb, nb, hd, qb=0.85, qd=0.6, n=60, rmax=40.0, noise=0.02, seed=1):
    rng = np.random.default_rng(seed)
    s = np.linspace(1.0, rmax, n)
    fd = 40000.0
    fb = fd * bt / (1 - bt)
    I = ad.sersic_profile(s, ad.sersic_ie_from_flux(fb, reb, nb, qb), reb, nb) + ad.exp_profile(s, fd / (2 * math.pi * hd * hd * qd), hd)
    dI = noise * I + 1e-3
    I = I + rng.normal(0, dI)
    eps = np.where(s < 2.0 * reb, 1 - qb, 1 - qd)
    pa = np.full(n, 40.0)
    gf = np.cumsum(2 * math.pi * s * np.maximum(I, 0) * (1 - eps) * np.gradient(s))
    return dict(sma=s.tolist(), intens=I.tolist(), intens_err=dI.tolist(), ellipticity=eps.tolist(), pa=pa.tolist(), growth_flux=gf.tolist())


def test_guess_from_profile_bulge_disc_and_single():
    tab = synthetic_profile(bt=0.3, reb=3.0, nb=2.5, hd=7.0)
    g = ad.guess_from_profile(ad.clean_profile(tab, FWHM), FWHM, ZP)
    assert g['ok'] and not g['guess_fallback']
    assert abs(g['bt_guess'] - 0.3) < 0.12
    b, d = g['comps']['bulge+disk']
    assert abs(d['re'] / (1.678 * 7.0) - 1) < 0.2 and abs(b['re'] / 3.0 - 1) < 0.5
    assert abs(d['q'] - 0.6) < 0.05 and abs(b['q'] - 0.85) < 0.1 and abs(d['pa'] - 40) < 1
    assert g['bic1d']['double'] < g['bic1d']['single'] - 10 and g['guess_type'] == 'bulge+disk'
    # a pure Sersic profile: the single-component guess carries the right size and index
    s = np.linspace(1, 40, 60)
    I = ad.sersic_profile(s, 20.0, 8.0, 2.0)
    tab = dict(sma=s.tolist(), intens=I.tolist(), intens_err=(0.01 * I).tolist(), ellipticity=[0.3] * 60, pa=[10.0] * 60,
               growth_flux=np.cumsum(2 * math.pi * s * I * 0.7 * np.gradient(s)).tolist())
    g = ad.guess_from_profile(ad.clean_profile(tab, FWHM), FWHM, ZP)
    assert abs(g['single']['n'] - 2.0) < 0.15 and abs(g['single']['re'] / 8.0 - 1) < 0.05 and g['guess_type'] == 'sersic'


def test_bar_feature():
    s = np.linspace(1, 40, 60)
    I = ad.sersic_profile(s, 20.0, 8.0, 1.0)
    base = dict(sma=s.tolist(), intens=I.tolist(), intens_err=(0.01 * I).tolist(), growth_flux=np.cumsum(2 * math.pi * s * I * 0.5 * np.gradient(s)).tolist())
    eps = 0.25 + 0.3 * np.exp(-0.5 * ((s - 10) / 2.5) ** 2)
    pa = np.where(s < 13, 30.0, 70.0)
    f = ad.structure_features(ad.clean_profile(dict(base, ellipticity=eps.tolist(), pa=pa.tolist()), FWHM), FWHM)
    assert f['bar'] and f['bar_twist'] > 30
    f = ad.structure_features(ad.clean_profile(dict(base, ellipticity=[0.25] * 60, pa=[30.0] * 60), FWHM), FWHM)
    assert not f['bar'] and not f['edgeon'] and abs(f['eps_out'] - 0.25) < 1e-9


def test_decompose_bulge_disc_recovers_truth():
    d, psf = gal(BD)
    o = ad.decompose(d, 45.0, 45.0, psf=psf, rms=1.0, zp=ZP, flux0=1.1 * 10 ** (-0.4 * (15.9 - ZP)), re0=7.0, q0=0.7, pa0=40, fwhm=FWHM)
    assert o['name'] == 'bulge+disk', (o['name'], o['notes'], o['bic'])
    b, dd = o['res']['components'][:2]
    bt = b['flux'] / (b['flux'] + dd['flux'])
    t_bt = 1 / (1 + 10 ** (-0.4 * (16.4 - 16.9)) ** -1) if False else (10 ** (-0.4 * 16.9)) / (10 ** (-0.4 * 16.9) + 10 ** (-0.4 * 16.4))
    assert abs(bt - t_bt) < 0.05
    assert abs(dd['re'] / 12.0 - 1) < 0.05 and abs(b['re'] / 3.0 - 1) < 0.25 and abs(dd['q'] - 0.6) < 0.03
    assert abs(o['res']['chi2_red'] - 1) < 0.15 and o['flag'] & ad.FLAGS['NOPROFILE'] == 0
    assert o['bic']['sersic'] > o['bic']['bulge+disk'] + 10


def test_decompose_single_sersic_stays_simple():
    d, psf = gal(SER, seed=3)
    o = ad.decompose(d, 45.0, 45.0, psf=psf, rms=1.0, zp=ZP, flux0=1.0e-0 * 10 ** (-0.4 * (16.4 - ZP)), re0=7.0, q0=0.7, pa0=60, fwhm=FWHM)
    assert o['name'] == 'sersic', (o['name'], o['notes'])
    c = o['res']['components'][0]
    assert abs(c['n'] - 2.0) < 0.15 and abs(c['re'] / 8.0 - 1) < 0.05


def test_nucleus_found_when_present():
    comps = [dict(kind='psf', x=45, y=45, mag=18.0)] + BD
    d, psf = gal(comps, seed=5)
    o = ad.decompose(d, 45.0, 45.0, psf=psf, rms=1.0, zp=ZP, flux0=1.0e-0 * 10 ** (-0.4 * (15.9 - ZP)), re0=7.0, q0=0.7, pa0=40, fwhm=FWHM)
    assert o['name'] == 'nucleus+bulge+disk', (o['name'], o['notes'], o['bic'])
    assert abs(o['res']['components'][0]['mag'] - 18.0) < 0.3
    o2 = ad.decompose(d, 45.0, 45.0, psf=psf, rms=1.0, zp=ZP, flux0=10 ** (-0.4 * (15.9 - ZP)), re0=7.0, q0=0.7, pa0=40, fwhm=FWHM, models=('sersic', 'bulge+disk'))
    assert o2['name'] in ('sersic', 'bulge+disk')


def test_no_profile_fallback_and_noise_only():
    rng = np.random.default_rng(2)
    d = rng.normal(0, 1.0, (61, 61))
    o = ad.decompose(d, 30.0, 30.0, psf=mf.gaussian_psf_array(FWHM), rms=1.0, zp=ZP, flux0=500.0, re0=4.0, q0=0.8, pa0=0, fwhm=FWHM, max_nfev=40)
    assert o['flag'] & ad.FLAGS['NOPROFILE'] or o['name'] is None or o['name'] == 'sersic'
    assert o['name'] in (None, 'sersic')            # nothing to decompose in noise: never a 2-component structure


def write_field(tmp):
    psf = mf.gaussian_psf_array(FWHM)
    comps = [dict(BD[0], x=60, y=60), dict(BD[1], x=60, y=60),
             dict(kind='sersic', x=170, y=60, mag=16.4, re=8.0, n=2.0, q=0.7, pa=60),
             dict(kind='psf', x=115, y=140, mag=17.0)]
    img = mf.render_model([mf._normalise(c, ZP) for c in comps], (200, 220), psf=psf) + np.random.default_rng(4).normal(0, 1.0, (200, 220))
    fits.writeto(os.path.join(tmp, 'f.fits'), img.astype(np.float32), overwrite=True)
    with open(os.path.join(tmp, 'c.tsv'), 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tKRON_RADIUS\tFLUX_AUTO\tFLUX_RADIUS\tMAG_AUTO\tCLASS_STAR\n')
        f.write('1\t61\t61\t9\t6\t30\t3.5\t%g\t7.0\t15.9\t0.0\n' % 10 ** (-0.4 * (15.9 - ZP)))
        f.write('2\t171\t61\t10\t7\t60\t3.5\t%g\t8.0\t16.4\t0.0\n' % 10 ** (-0.4 * (16.4 - ZP)))
        f.write('3\t116\t141\t2\t2\t0\t3.5\t%g\t1.8\t17.0\t0.95\n' % 10 ** (-0.4 * (17.0 - ZP)))


def read_tsv(text):
    L = [l for l in text.strip().split('\n') if l and not l.startswith('#')]
    cols = L[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t') + [''] * len(cols))) for l in L[1:]]


def test_cli_batch_decomp(tmp_path):
    write_field(str(tmp_path))
    work = tmp_path / 'w'
    cmd = [sys.executable, CLI, str(tmp_path / 'f.fits'), '--catalog', str(tmp_path / 'c.tsv'), '--work', str(work), '--model', 'decomp', '--columns', 'ad', '--psf-fwhm', str(FWHM),
           '--mag-zeropoint', str(ZP), '--rms', '1.0', '--neighbours', 'mask', '--n-workers', '2', '--max-objects', '2', '--sky', 'const']
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    cols, rows = read_tsv(p.stdout)
    import multifit as mfcli
    assert cols == ['NUMBER'] + mfcli.AD_COLUMNS
    r1, r2, r3 = rows
    assert r1['AD_TYPE'] == '2' and abs(float(r1['AD_BT']) - 0.387) < 0.15 and abs(float(r1['AD_RED']) / 12 - 1) < 0.1 and float(r1['AD_DBIC']) > 10
    assert r2['AD_TYPE'] == '1' and abs(float(r2['AD_N1']) - 2.0) < 0.3 and r2['AD_BT'] == ''
    assert r3['AD_TYPE'] == ''                                               # not fitted (max-objects 2 = the two brightest)
    assert abs(float(r1['AD_MAG']) - 15.9) < 0.05
    for f in ('multifit_decomp.tsv', 'multifit_decomp_plot.png', 'multifit_results.tsv', 'multifit_model.fits', 'multifit_residual.fits'):
        assert (work / f).stat().st_size > 100
    dc, dr = read_tsv(open(work / 'multifit_decomp.tsv').read())
    assert len(dr) == 2 and {r['DEC_MODEL'] for r in dr} == {'bulge+disk', 'sersic'}
    assert all(float(r['DEC_PROFILE']) == 1 for r in dr)
    # --model decomp with the GF_ columns, nucleus off
    p = subprocess.run(cmd[:cmd.index('--columns')] + cmd[cmd.index('--columns') + 2:] + ['--nucleus', 'off', '--objects', '1'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-500:]
    cols, rows = read_tsv(p.stdout)
    assert 'GF_BT' in cols and abs(float(rows[0]['GF_BT']) - 0.387) < 0.15 and rows[0]['GF_NCOMP'] == '2'
