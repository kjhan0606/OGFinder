"""PSF extensions (ogfkit/psfext.py): reduced rank, wings / extended PSF, error budget, stacked-star PSF, psfex and stacking CLI options."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import models as M, psfmodel as PM, psfext as PX, synth, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(PLUGIN, '..', '..'))
PY = sys.executable
SHAPE = (600, 600)
KW = dict(fwhm0=2.8, fwhm1=3.8, e0=0.04, e1=0.14)


@pytest.fixture(scope='module')
def field():
    img, t = synth.star_field(SHAPE, n=1200, mag_range=(15.5, 23), zp=25, sky=50, noise=0.5, seed=3, **KW)
    return img.astype(float), t


@pytest.fixture(scope='module')
def model(field):
    return PM.build_psf_model(field[0], fwhm_prior=3.3, snr_min=6, max_stars=400, degree=3)


def test_reduce_rank_identity_and_truncation(model):
    mdl, _ = model
    full = PX.reduce_rank(mdl, 99)
    assert np.allclose(full.coeffs, mdl.coeffs)
    sv = PX.singular_values(mdl)
    assert len(sv) == mdl.coeffs.shape[0] - 1 and sv[0] >= sv[-1]
    r1 = PX.reduce_rank(mdl, 1)
    assert r1.meta['rank'] == 1 and 0.0 < r1.meta['rank_variance_kept'] <= 1.0
    s = np.linalg.svd(r1.coeffs[1:].reshape(len(r1.coeffs) - 1, -1), compute_uv=False)
    assert s[1] < 1e-8 * s[0] and np.allclose(r1.coeffs[0], mdl.coeffs[0])         # rank one variation, same mean
    for x, y in ((50, 50), (300, 300), (550, 100)):
        a, b = mdl.stamp(x, y), r1.stamp(x, y)
        r3 = PX.reduce_rank(mdl, 3).stamp(x, y)
        assert abs(b.sum() - 1) < 1e-9
        e1, e3 = np.abs(a - b).max() / a.max(), np.abs(a - r3).max() / a.max()
        assert e3 <= e1 + 1e-9 and e3 < 0.12                                           # more singular images -> closer to the full polynomial model
    r0 = PX.reduce_rank(mdl, 0)
    assert np.allclose(r0.stamp(10, 10), r0.stamp(500, 500))                           # no variation left


def halo_stars(seed=1, n=12, halo=0.1, size=101):
    rng = np.random.default_rng(seed)
    img = rng.normal(0, 1.0, (500, 500)) + 100.0
    core = M.moffat_psf(size, 3.0, beta=2.5)
    wide = M.moffat_psf(size, 10.0, beta=1.7, nsub=1)
    p = (1 - halo) * core + halo * wide
    xy = []
    for i in range(n):
        x, y = 60 + (i % 4) * 110 + rng.uniform(-3, 3), 60 + (i // 4) * 150 + rng.uniform(-3, 3)
        ix, iy = int(round(x)), int(round(y))
        h = size // 2
        img[iy - h:iy + h + 1, ix - h:ix + h + 1] += 3e6 * p * 1.0
        xy.append((ix, iy))
    return img, xy, p


def test_wings_measure_fit_and_extended_psf():
    img, xy, p = halo_stars()
    r, prof, ns = PX.measure_wings(img, xy, bkg=100.0, rcore=6.0, rmax=40.0, min_dist=30)
    assert ns >= 8
    w = PX.fit_wings(r, prof, rmin=10.0, rmax=40.0)
    assert w['valid'] and 2.6 < w['gamma'] < 3.8 and w['rms_dex'] < 0.1                  # Moffat beta 1.7 wing: index between 2 beta and the local slope
    yy, xx = np.mgrid[-50:51, -50:51]
    rr = np.hypot(xx, yy)
    fc = p[rr < 6].sum()
    true = p[(rr > 10) & (rr <= 40)].sum() / fc
    assert abs(PX.wing_flux(w, 10.0, 40.0) / true - 1) < 0.25
    # flat profile -> not valid
    assert PX.fit_wings(np.arange(1.0, 30.0), np.full(29, 1e-3), rmin=5.0)['valid'] is False
    mdl = PM.analytic(31, 3.0, 'moffat')
    ext, info = PX.extended_psf(mdl, w, size=101)
    assert ext.shape == (101, 101) and abs(ext.sum() + 0 - (1 - info['tail_fraction'])) < 1e-6 and ext.min() >= 0
    assert ext[50, 50] > ext[50, 80] > ext[50, 99]


def test_error_budget_zero_and_signs():
    a = M.moffat_psf(31, 3.0, beta=2.5)
    z = PX.error_budget(a, a, profiles=((1.0, 4.0),), fit=True)[0]
    assert z['max_abs_over_peak'] < 1e-12 and abs(z['dmag']) < 1e-3 and abs(z['dre_rel']) < 1e-2
    b = M.moffat_psf(31, 3.0 * 1.15, beta=2.5)                           # used PSF too wide for the data -> model fitted smaller / brighter core
    e = PX.error_budget(a, b, profiles=((1.0, 3.0),), fit=True)[0]
    assert e['max_abs_over_peak'] > 0.01 and e['dre_rel'] < -0.02          # a wider PSF makes the intrinsic size come out smaller (n=1; for n=4 the n-R_e degeneracy flips it)
    c = PX.error_budget(a, a[8:23, 8:23] / a[8:23, 8:23].sum(), profiles=((4.0, 8.0),), fit=False)[0]
    assert abs(c['flux_rel']) < 0.02


def test_stack_star_psf_positive_and_close_to_truth(field):
    img, t = field
    bkg, rms = PM.background(img)
    cand = PM.find_stars(img, bkg, rms, thresh=10.0)
    ok = np.where((cand['snr'] > 60) & (cand['fr50'] < 2.4))[0]
    xy = np.c_[cand['x'][ok], cand['y'][ok]]
    xy = xy[(xy[:, 0] > 150) & (xy[:, 0] < 450) & (xy[:, 1] > 150) & (xy[:, 1] < 450)]
    assert len(xy) >= 3
    psf, n = PX.stack_star_psf(img, xy, size=25, bkg=bkg)
    assert n >= 3 and psf.min() >= 0 and abs(psf.sum() - 1) < 1e-9
    sh = PM.psf_shape(psf)
    mdl, _ = PM.build_psf_model(img, fwhm_prior=3.3, snr_min=6, max_stars=400)
    ref = PM.psf_shape(mdl.stamp(300, 300, 0, 0, 25))
    assert abs(sh['fwhm'] / ref['fwhm'] - 1) < 0.15


def test_cli_psfex_rank_wings_assess_and_stacking_psf_stamp(field, tmp_path):
    img, t = field
    f = str(tmp_path / 'f.fits')
    fits.writeto(f, img.astype(np.float32), overwrite=True)
    wd = str(tmp_path / 'w')
    r = subprocess.run([PY, os.path.join(PLUGIN, 'psfex.py'), f, '--work', wd, '--snr-min', '6', '--fwhm', '3.3', '--max-stars', '400', '--order', '3', '--rank', '1', '--wings', '--assess', '--ext-size', '81'],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-500:]
    assert 'PCA rank 1' in r.stdout and 'wings:' in r.stdout
    mdl = PM.load_model(os.path.join(wd, 'psfex_model.json'))
    assert mdl.meta['rank'] == 1
    b = json.load(open(os.path.join(wd, 'psfex_budget.json')))
    assert 'centre_psf_at_corner' in b and len(b['centre_psf_at_corner']) == 2
    assert os.path.getsize(os.path.join(wd, 'psfex_wings.json')) > 50
    # stacking --psf-stamp on the PSF stars
    cols, rows = tsvio.read_catalog(os.path.join(wd, 'psfex_stars.tsv'))
    pos = str(tmp_path / 'pos.tsv')
    sel = [q for q in rows if int(float(q['USED'])) == 1][:25]
    tsvio.write_table(pos, ['NUMBER', 'X_IMAGE', 'Y_IMAGE'], [dict(NUMBER=i + 1, X_IMAGE=q['X_IMAGE'], Y_IMAGE=q['Y_IMAGE']) for i, q in enumerate(sel)])
    sd = str(tmp_path / 's')
    r = subprocess.run([PY, os.path.join(ROOT, 'plugins', 'stacking', 'stacking.py'), f, '--work', sd, '--catalog', pos, '--half', '12', '--n-boot', '5', '--psf-stamp'], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-400:]
    p = fits.getdata(os.path.join(sd, 'stacking_psf.fits'))
    assert p.shape == (25, 25) and abs(p.sum() - 1) < 1e-5 and p.min() >= 0
