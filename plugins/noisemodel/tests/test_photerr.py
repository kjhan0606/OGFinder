"""Unified photometric error model (ogfkit/photerr.py) and its wiring (noise.flux_error, depth, daophot.phot, noisemodel --sky-annulus / --psf)."""
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from ogfkit import noise as nz, photerr as pe, models as M, tsvio, imageio, depth as dp, daophot as dao

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def cnoise(shape, sigma, sk, seed):
    w = np.random.default_rng(seed).normal(size=shape)
    c = gaussian_filter(w, sk, mode='wrap') if sk > 0 else w
    return c * (sigma / c.std())


def test_white_noise_limit_is_legacy_formula():
    N, sp, nsky, f, g = 50.0, 3.0, 400.0, 1000.0, 2.0
    e = pe.aperture_error(f, N, sp, None, g, nsky=nsky)
    assert abs(e['ap'] ** 2 - (N * sp ** 2 * (1 + N / nsky) + f / g)) < 1e-9          # the old daophot.phot variance
    assert abs(float(nz.flux_error(f, N, sp, None, g, False)) ** 2 - (N * sp ** 2 + f / g)) < 1e-9
    law = dict(alpha=1.3, beta=0.7, ok=True)
    assert abs(float(nz.flux_error(f, N, sp, law, None, True)) - sp * 1.3 * N ** 0.7) < 1e-9
    assert abs(pe.NoiseModel(sp, 1.3, 0.7, True).sigma_sum(N) - sp * 1.3 * N ** 0.7) < 1e-9


def test_apcorr_error_term_and_mag_error():
    e = pe.aperture_error(1000.0, 30.0, 2.0, None, None, frac=0.8, frac_err=0.02)
    assert abs(e['flux_total'] - 1250.0) < 1e-9
    assert abs(e['total'] ** 2 - (e['ap'] ** 2 / 0.64 + (1250.0 * 0.02 / 0.8) ** 2)) < 1e-6
    assert abs(float(pe.mag_error(100.0, 5.0)) - 1.0857 * 0.05) < 1e-3


def test_blank_positions_shared_with_depth():
    mask = np.zeros((200, 200), bool)
    mask[:, :60] = True
    a = dp.blank_circle_positions(mask, 4, 50, np.random.default_rng(3), box=(0, 0, 199, 199))
    b = nz.blank_positions(mask, 4, 50, np.random.default_rng(3), box=(0, 0, 199, 199), max_tries=40, kmin=100)
    assert np.array_equal(a, b) and (a[:, 0] > 60).all()
    c = dp.blank_circle_positions(mask, 4, 20, np.random.default_rng(1), box=(100, 100, 150, 150))
    assert ((c >= 100) & (c <= 150)).all()


def test_blank_aperture_local_sky_pulls_correlated_noise():
    img = 100.0 + cnoise((700, 700), 5.0, 1.0, 5)
    mask = nz.source_mask(img)
    sub = img - np.median(img)
    sp = nz.robust_std(sub)
    model = pe.NoiseModel.measure(sub, mask, seed=1).add_local(img, mask, (2, 3, 5, 8), (12.0, 18.0), 'median', n=400, seed=2)
    assert model.ok and model.local['ok']
    pos = nz.blank_positions(mask, 19, 500, np.random.default_rng(9))
    for r in (2.0, 4.0, 8.0):
        ph = pe.aperture_photometry(img, pos, r, (12.0, 18.0), 'median', mask=mask)
        f, N, ns = ph['flux'], ph['N'], ph['nsky']
        loc = f / pe.aperture_error(f, N, sp, model, local=True)['ap']
        naive = f / (sp * np.sqrt(N))
        legacy = f / np.sqrt(N * sp ** 2 * (1 + N / ns))
        assert 0.85 < loc.std() < 1.2, (r, loc.std())
        assert naive.std() > 1.8 and legacy.std() > 1.7 and naive.std() > 1.5 * loc.std()
    # sky-estimation term alone for white noise: N^2 sigma^2 / nsky within the median-estimator inflation
    w = 100.0 + cnoise((600, 600), 5.0, 0.0, 6)
    mk = nz.source_mask(w)
    kap, nsky, nu = pe.annulus_factor(w - 100, mk, 12.0, 18.0, 'median', n=400, seed=1)
    assert 1.05 < kap < 1.5 and nu > 100                                                 # median of white noise: sqrt(pi/2) = 1.25


def test_neighbour_contamination_matches_direct_measurement():
    size = 61
    h = size // 2
    img = np.zeros((140, 200))
    xy = np.array([[60.0, 70.0], [72.0, 70.0], [140.0, 70.0]])
    fl = np.array([1000.0, 4000.0, 2000.0])
    for (x, y), f in zip(xy, fl):
        ix, iy = int(round(x)), int(round(y))
        img[iy - h:iy + h + 1, ix - h:ix + h + 1] += f * M.moffat_psf(size, 3.0, 2.5, dx=x - ix, dy=y - iy)
    r, ann = 4.0, (8.0, 14.0)
    ph = pe.aperture_photometry(img, xy, r, ann, 'mean')
    psf = M.moffat_psf(size, 3.0, 2.5)
    frac = pe.aperture_fraction(psf, r, ann)
    c = pe.neighbour_contamination(xy, fl, r, psf, sky_annulus=ann)
    # sky mean used above, the model assumes the annulus *mean* too -> exact bookkeeping; the median differs, so compare with mean
    expect = ph['flux'] - frac * fl                      # measured excess over the source's own net flux
    assert expect[0] < -100                               # the 4000-count neighbour sits in the sky annulus: the sky is over-subtracted
    assert abs(c['contam'][0] - expect[0]) < 0.15 * abs(expect[0]) + 1.0
    assert abs(c['contam'][2]) < 1.0                     # far star: no neighbour within the stamp
    assert c['n_neigh'][0] >= 1 and c['n_neigh'][2] == 0


def test_daophot_phot_default_unchanged_and_noise_option(tmp_path):
    img = 100.0 + cnoise((400, 400), 5.0, 1.0, 8)
    mask = nz.source_mask(img)
    xy = nz.blank_positions(mask, 19, 20, np.random.default_rng(2))
    res = dao.phot(img, xy, radii=(3.0, 5.0), sky_inner=12.0, sky_outer=18.0, sky_mode='median', gain=None)
    for p in res:
        for k, r in enumerate((3.0, 5.0)):
            area = math.pi * r * r
            assert abs(p['fluxerr'][k] - math.sqrt(area * p['skysig'] ** 2 * (1 + area / p['nsky']))) / p['fluxerr'][k] < 0.02       # exact overlap area vs pi r^2
    sub = img - np.median(img)
    model = pe.NoiseModel.measure(sub, mask, seed=1).add_local(img, mask, (2, 3, 5, 8), (12.0, 18.0), 'median', n=300, seed=2)
    res2 = dao.phot(img, xy, radii=(3.0, 5.0), sky_inner=12.0, sky_outer=18.0, sky_mode='median', noise=model)
    ratio = np.array([p2['fluxerr'][1] / p1['fluxerr'][1] for p1, p2 in zip(res, res2)])
    assert 1.8 < np.median(ratio) < 4.0                                                   # correlated noise: errors grow vs the white-noise formula


def test_noisemodel_cli_unified_columns(tmp_path):
    rng = np.random.default_rng(4)
    shape = (700, 700)
    img = 100.0 + cnoise(shape, 5.0, 1.0, 4)
    size = 61
    h = size // 2
    x = rng.uniform(40, 660, 120); y = rng.uniform(40, 660, 120)
    f = np.full(120, 600.0)
    x[1], y[1] = x[0] + 9.0, y[0]; f[1] = 6000.0                       # bright neighbour near star 0
    add = np.zeros(shape)
    for xi, yi, fi in zip(x, y, f):
        ix, iy = int(round(xi)), int(round(yi))
        add[iy - h:iy + h + 1, ix - h:ix + h + 1] += fi * M.moffat_psf(size, 3.0, 2.5, dx=xi - ix, dy=yi - iy)
    imageio.save_fits(str(tmp_path / 'im.fits'), (img + add).astype(np.float32))
    imageio.save_fits(str(tmp_path / 'psf.fits'), M.moffat_psf(size, 3.0, 2.5).astype(np.float32))
    with open(tmp_path / 'cat.tsv', 'w') as fh:
        fh.write('NUMBER\tX_IMAGE\tY_IMAGE\tFLUX_AUTO\n')
        for i in range(120):
            fh.write('%d\t%.4f\t%.4f\t%.2f\n' % (i + 1, x[i] + 1, y[i] + 1, f[i]))
    args = [PY, os.path.join(os.path.dirname(HERE), 'noisemodel.py'), str(tmp_path / 'im.fits'), '--work', str(tmp_path / 'w'), '--catalog', str(tmp_path / 'cat.tsv'), '--aperture', 'fixed',
            '--aper-radius', '4', '--n-aper', '300', '--radii', '2,3,4,6,8', '--meta-out', str(tmp_path / 'meta.json')]
    base = subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip().split('\n')
    assert 'NM_FLUXERR_LOC' not in base[0] or all(l.split('\t')[base[0].split('\t').index('NM_FLUXERR_LOC')] in ('', 'nan', 'NaN', '-') for l in base[1:3])
    out = subprocess.run(args + ['--sky-annulus', '12,18', '--psf', str(tmp_path / 'psf.fits')], capture_output=True, text=True, check=True).stdout.strip().split('\n')
    cols = out[0].split('\t')
    rows = [dict(zip(cols, l.split('\t'))) for l in out[1:]]
    c0 = tsvio.fnum(rows[0]['NM_CONTAM'])
    assert abs(c0) > 8 and abs(tsvio.fnum(rows[5]["NM_CONTAM"])) < abs(c0) / 2                                  # star 0 sees the bright companion
    tot = tsvio.fnum(rows[0]['NM_FLUXERR_TOT']); loc = tsvio.fnum(rows[0]['NM_FLUXERR_LOC'])
    assert tot > loc > tsvio.fnum(rows[0]['NM_FLUXERR']) * 0.6 and tsvio.fnum(rows[0]['NM_SKYERR']) > 0
    import json
    assert json.load(open(tmp_path / 'meta.json'))['noisemodel']['local']['ok']
