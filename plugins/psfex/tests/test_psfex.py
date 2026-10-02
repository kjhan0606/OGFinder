"""PSF model plugin: recovery of a spatially varying truth, fallbacks, CLI/catalog columns, JSON round trip, and the consumers (PSF photometry, crowded photometry, tiled ZOGY)."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import psfmodel as PM, synth, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(PLUGIN, '..', '..'))
PY = sys.executable
SHAPE = (600, 600)
KW = dict(fwhm0=2.8, fwhm1=3.8, e0=0.04, e1=0.14)
PTS = [(x, y) for x in (8, 150, 300, 450, 590) for y in (8, 300, 590)]


@pytest.fixture(scope='module')
def field():
    img, t = synth.star_field(SHAPE, n=1200, mag_range=(15.5, 23), zp=25, sky=50, noise=0.5, seed=3, **KW)
    return img, t


@pytest.fixture(scope='module')
def model(field):
    return PM.build_psf_model(field[0].astype(float), fwhm_prior=3.3, snr_min=6, max_stars=400)


def truth_shapes():
    return [PM.psf_shape(synth.star_stamp(x, y, SHAPE, size=31, **KW)[0]) for x, y in PTS]


def test_recovery_of_spatially_varying_psf(model):
    m, info = model
    tr = truth_shapes()
    ss = [PM.psf_shape(m.stamp(x, y)) for x, y in PTS]
    d = np.array([s['fwhm'] / t['fwhm'] - 1 for s, t in zip(ss, tr)])
    de = np.array([s['e'] - t['e'] for s, t in zip(ss, tr)])
    print('PSFEX recovery: n_used %s order %d oversample %d | FWHM rel err rms %.3f bias %+.3f max %.3f | e err rms %.3f max %.3f' % (
        info['n_used'], m.degree, m.oversample, d.std(), d.mean(), abs(d).max(), de.std(), abs(de).max()))
    assert m.degree == 2 and m.meta['mode'] == 'empirical'
    assert d.std() < 0.03 and abs(d).max() < 0.07 and abs(d.mean()) < 0.03
    assert abs(de).max() < 0.035
    # the model genuinely varies: FWHM range of the model matches the truth range
    f = [s['fwhm'] for s in ss]
    assert max(f) - min(f) > 0.7 * (max(t['fwhm'] for t in tr) - min(t['fwhm'] for t in tr))


def test_unit_sum_and_symmetry_of_stamps(model):
    m, _ = model
    for x, y in PTS[::4]:
        st = m.stamp(x, y, 0.3, -0.2)
        assert abs(st.sum() - 1) < 1e-9 and st.shape == (m.size, m.size)
    assert m.at(300, 300).shape == (m.size, m.size) and m.at(10, 10, 11).shape == (11, 11)


def test_order_selection_and_constant_psf_is_worse(field):
    img = field[0].astype(float)
    m0, _ = PM.build_psf_model(img, fwhm_prior=3.3, snr_min=6, degree=0)
    m2, _ = PM.build_psf_model(img, fwhm_prior=3.3, snr_min=6, degree=2)
    tr = truth_shapes()
    r0 = np.std([PM.psf_shape(m0.stamp(x, y))['fwhm'] / t['fwhm'] - 1 for (x, y), t in zip(PTS, tr)])
    r2 = np.std([PM.psf_shape(m2.stamp(x, y))['fwhm'] / t['fwhm'] - 1 for (x, y), t in zip(PTS, tr)])
    print('PSFEX FWHM error rms: order 0 %.3f, order 2 %.3f' % (r0, r2))
    assert r2 < 0.5 * r0 and not m0.varies and m2.varies


def test_fallbacks_with_few_stars():
    # no stars at all -> Gaussian prior
    rng = np.random.default_rng(1)
    empty = rng.normal(0, 1, (200, 200))
    m, info = PM.build_psf_model(empty, fwhm_prior=3.0)
    assert info['mode'] == 'gaussian_prior' and abs(m.stamp(100, 100).sum() - 1) < 1e-9
    assert abs(PM.psf_shape(m.stamp(100, 100))['fwhm'] - 3.0) < 0.3
    # a handful of stars -> constant / Moffat fit, valid model
    img, t = synth.star_field((200, 200), n=6, mag_range=(16, 17), zp=25, sky=20, noise=0.5, seed=4, min_sep=30, **KW)
    m, info = PM.build_psf_model(img.astype(float), fwhm_prior=3.3, snr_min=10)
    print('PSFEX few stars: mode %s n_stars %s degree %d' % (info['mode'], info.get('n_stars'), m.degree))
    assert m.degree == 0 and info['mode'] in ('empirical', 'moffat_fit')
    f = PM.psf_shape(m.stamp(100, 100))['fwhm']
    assert 2.6 < f < 4.4


def test_json_and_fits_round_trip(model, tmp_path):
    m, _ = model
    p = str(tmp_path / 'm.json')
    m.save(p)
    m2 = PM.load_model(p)
    assert np.allclose(m.stamp(123, 456, 0.1, 0.2), m2.stamp(123, 456, 0.1, 0.2))
    json.dumps(m.summary())
    m.save(str(tmp_path / 'm.fits'))
    assert os.path.getsize(str(tmp_path / 'm.fits')) > 1000


def test_cli_with_catalog_and_tables(field, tmp_path):
    img, t = field
    f = str(tmp_path / 'f.fits')
    fits.writeto(f, img, overwrite=True)
    cat = str(tmp_path / 'c.tsv')
    tsvio.write_table(cat, ['NUMBER', 'X_IMAGE', 'Y_IMAGE'], [dict(NUMBER=i + 1, X_IMAGE=x + 1, Y_IMAGE=y + 1) for i, (x, y) in enumerate(PTS)])
    wd = str(tmp_path / 'w')
    r = subprocess.run([PY, os.path.join(PLUGIN, 'psfex.py'), f, '--work', wd, '--catalog', cat, '--snr-min', '6', '--fwhm', '3.3', '--max-stars', '400', '--grid', '5'],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    lines = r.stdout.strip().splitlines()
    head = lines[0].split('\t')
    assert head == ['NUMBER', 'PSFM_FWHM', 'PSFM_E', 'PSFM_PA', 'PSFM_NSTAR'] and len(lines) == len(PTS) + 1
    fw = [float(l.split('\t')[1]) for l in lines[1:]]
    tr = [s['fwhm'] for s in truth_shapes()]
    assert np.max(np.abs(np.array(fw) / np.array(tr) - 1)) < 0.07
    for n in ('model.json', 'model.fits', 'center.fits', 'stars.tsv', 'maps.tsv', 'info.json', 'diag.png'):
        assert os.path.getsize(os.path.join(wd, 'psfex_' + n)) > 100, n
    info = json.load(open(os.path.join(wd, 'psfex_info.json')))
    assert len(info['diag']['fwhm']) == 5 and info['summary']['degree'] == 2
    cols, rows = tsvio.read_catalog(os.path.join(wd, 'psfex_stars.tsv'))
    assert len(rows) >= 40 and 'FWHM' in cols
    # user supplied star list
    sl = str(tmp_path / 'stars.tsv')
    sel = rows[:30]
    tsvio.write_table(sl, ['ID', 'X_IMAGE', 'Y_IMAGE'], [dict(ID=i, X_IMAGE=r['X_IMAGE'], Y_IMAGE=r['Y_IMAGE']) for i, r in enumerate(sel)])
    r = subprocess.run([PY, os.path.join(PLUGIN, 'psfex.py'), f, '--work', str(tmp_path / 'w2'), '--star-list', sl, '--fwhm', '3.3'], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0 and 'PSF model' in r.stdout, r.stderr


def _isolated_truth():
    img, t = synth.star_field(SHAPE, n=350, mag_range=(16.0, 19.5), zp=25, sky=50, noise=0.5, seed=8, min_sep=16, **KW)
    return img.astype(float), t


def test_psf_photometry_uses_model(model):
    """PSF photometry with the spatially varying model vs the single central PSF: flux bias / scatter across the field."""
    from psf_phot.config import PSFPhotConfig
    from psf_phot.photometry import do_psf_photometry
    m, _ = model
    img, t = _isolated_truth()
    rng = np.random.default_rng(5)
    src = [dict(number=i + 1, x=float(x + rng.normal(0, 0.2)), y=float(y + rng.normal(0, 0.2))) for i, (x, y) in enumerate(zip(t['x'], t['y']))]
    data = img - 50.0
    cfg = PSFPhotConfig(fit_radius=9, mag_zeropoint=25.0)
    const = np.asarray(m.stamp(300, 300))
    res = {}
    for name, kw in (('constant', {}), ('model', dict(psf_model=m))):
        out = do_psf_photometry(data, const, src, cfg, n_workers=0, **kw)
        fl = np.array([o['FLUX_PSF'] for o in out])
        res[name] = fl / t['flux'] - 1
    for k, v in res.items():
        ok = np.isfinite(v)
        print('PSFEX psf_phot %-8s bias %+.4f rms %.4f (n=%d)' % (k, np.median(v[ok]), 1.4826 * np.median(np.abs(v[ok] - np.median(v[ok]))), ok.sum()))
    rc = np.nanstd(res['constant'][np.abs(res['constant']) < 0.5]); rm = np.nanstd(res['model'][np.abs(res['model']) < 0.5])
    assert rm < 0.75 * rc
    x = t['x']
    # bias is flat across the field with the model (left vs right edge), not with the constant PSF
    lo, hi = x < 150, x > 450
    assert abs(np.nanmedian(res['model'][lo]) - np.nanmedian(res['model'][hi])) < 0.03
    assert abs(np.nanmedian(res['constant'][lo]) - np.nanmedian(res['constant'][hi])) > abs(np.nanmedian(res['model'][lo]) - np.nanmedian(res['model'][hi]))


def test_crowded_photometry_uses_model(model):
    from crowded_phot.config import CrowdedPhotConfig
    from crowded_phot.iterate import crowded_photometry
    m, _ = model
    img, t = synth.star_field((300, 300), n=450, mag_range=(16.0, 21.0), zp=25, sky=50, noise=0.5, seed=9, **KW)
    data = img.astype(np.float64) - 50.0
    sel = t['mag'] < 19.5
    src = [dict(number=i + 1, x=float(x), y=float(y)) for i, (x, y) in enumerate(zip(t['x'][sel], t['y'][sel]))]
    fl_true = t['flux'][sel]
    const = np.asarray(m.stamp(150, 150))
    out = {}
    for name, kw in (('constant', {}), ('model', dict(psf_model=m))):
        res = crowded_photometry(data, const, [dict(q) for q in src], CrowdedPhotConfig(), n_workers=0, **kw)   # (the input list is extended in place)
        by = {r['NUMBER']: r for r in res}
        fl = np.array([by[i + 1]['FLUX_CROWD'] if i + 1 in by else np.nan for i in range(len(src))])
        out[name] = fl / fl_true - 1
    for k, v in out.items():
        ok = np.isfinite(v) & (np.abs(v) < 1)
        print('PSFEX crowded %-8s median %+.4f rms %.4f (n=%d of %d)' % (k, np.median(v[ok]), 1.4826 * np.median(np.abs(v[ok] - np.median(v[ok]))), ok.sum(), len(v)))
    ok = np.isfinite(out['model']) & (np.abs(out['model']) < 1)
    assert ok.sum() > 0.8 * len(src)
    # (the legacy NSTAR-style fitter is biased low by ~8 % in this crowded field - stamp truncation, linear sub-pixel shifts, free local sky - with either PSF; unchanged here)
    assert abs(np.median(out['model'][ok])) < 0.12
    rc = 1.4826 * np.median(np.abs(out['constant'] - np.median(out['constant'])))
    rm = 1.4826 * np.median(np.abs(out['model'][ok] - np.median(out['model'][ok])))
    assert rm < 0.7 * rc


def test_zogy_tiled_with_model_field():
    """Tiled ZOGY with the PSF model of each image vs one constant PSF: artifacts at bright stars and transient flux recovery (spatially varying seeing)."""
    from moving import zogy as Z
    rng = np.random.default_rng(11)
    n = 400
    x = rng.uniform(10, 590, n); y = rng.uniform(10, 590, n)
    flux = 10 ** (-0.4 * (rng.uniform(15.5, 21, n) - 25))
    sig = 0.5
    kN = dict(fwhm0=2.4, fwhm1=4.2, e0=0.02, e1=0.18)
    kR = dict(fwhm0=2.2, fwhm1=2.6, e0=0.03, e1=0.06)
    tx, ty = [], []
    while len(tx) < 20:
        a, b = rng.uniform(40, 560, 2)
        if np.min(np.hypot(x - a, y - b)) > 14:
            tx.append(a); ty.append(b)
    tx, ty = np.array(tx), np.array(ty)
    tf = sig * 25 / 0.22
    N = synth.render_stars(SHAPE, np.r_[x, tx], np.r_[y, ty], np.r_[flux, np.full(20, tf)], sky=0, noise=sig, seed=21, **kN).astype(float)
    R = synth.render_stars(SHAPE, x, y, flux, sky=0, noise=sig * 0.7, seed=22, **kR).astype(float)
    mN, _ = PM.build_psf_model(N, fwhm_prior=3.3, snr_min=6, max_stars=400)
    mR, _ = PM.build_psf_model(R, fwhm_prior=2.4, snr_min=6, max_stars=400)
    bright = (flux > 10 ** (-0.4 * (18 - 25))) & (x > 20) & (x < 580) & (y > 20) & (y < 580)
    stats = {}
    for name, (pn, pr, tile) in {'constant': (mN.at(300, 300), mR.at(300, 300), (600, 600)), 'model': (mN, mR, (150, 150))}.items():
        o = Z.zogy_tiled(N, R, pn, pr, sig, sig * 0.7, tile=tile, margin=32)
        a = o['alpha']
        rec = np.array([a[int(round(yy)) - 1:int(round(yy)) + 2, int(round(xx)) - 1:int(round(xx)) + 2].max() for xx, yy in zip(tx, ty)]) / tf
        art = np.array([a[int(round(yy)), int(round(xx))] for xx, yy in zip(x[bright], y[bright])]) / flux[bright]
        stats[name] = (np.median(rec), 1.4826 * np.median(np.abs(rec - np.median(rec))), np.sqrt(np.mean(art ** 2)))
        print('PSFEX zogy %-8s transient flux ratio median %.3f scatter %.3f | residual at bright stars (fraction of flux) rms %.3f' % ((name,) + stats[name]))
    assert stats['model'][2] < 0.4 * stats['constant'][2]
    assert stats['model'][2] < 0.06
    assert abs(stats['model'][0] - 1) < 0.15 and stats['model'][1] <= stats['constant'][1] * 1.1
