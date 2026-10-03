import json
import math
import os
import subprocess
import sys

import numpy as np

from ogfkit import lensmodel as L, lenssynth as S, lensextra as X

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, '..', 'lensmodel.py')

P = dict(theta_E=1.3, q=0.8, phi=30.0, gamma=0.04, phi_g=100.0, x0=0.0, y0=0.0)


def _quad(model, b, half=4.0):
    ims = model.solve_images(b[0], b[1], half=half, n=301)
    assert len(ims) == 4
    return np.array([[i['x'], i['y']] for i in ims])


def test_weights_and_scaled_lens():
    assert abs(X.dist_weight(0.2, 0.6, 0.6) - 1.0) < 1e-12
    w = [X.dist_weight(0.2, z, 0.6) for z in (0.4, 0.6, 1.0, 2.0, 5.0)]
    assert all(np.diff(w) > 0) and w[0] < 1 < w[-1]
    m = L.build_sie_shear(P)
    s1 = X.ScaledLens(m, 1.0)
    a = m.alpha(0.7, -0.4)
    b = s1.alpha(0.7, -0.4)
    assert abs(a[0] - b[0]) < 1e-15 and abs(a[1] - b[1]) < 1e-15
    s2 = X.ScaledLens(m, 1.5)
    assert abs(s2.alpha(0.7, -0.4)[0] - 1.5 * a[0]) < 1e-14


def test_fit_multiplane_noise_free_and_free_weight():
    m = L.build_sie_shear(P)
    w2 = X.dist_weight(0.25, 2.0, 0.6)
    g1 = _quad(m, (0.04, 0.03))
    g2 = _quad(X.ScaledLens(m, w2), (-0.05, 0.04), half=6.0)
    r = X.fit_multiplane([dict(xy=g1, z=0.6), dict(xy=g2, z=2.0)], (0, 0), sigma=0.005, z_l=0.25, n_starts=16)
    assert r['chi2'] < 1e-3
    assert abs(r['params']['theta_E'] - P['theta_E']) < 1e-3 and abs(r['params']['q'] - P['q']) < 2e-3
    r2 = X.fit_multiplane([dict(xy=g1, z=0.6), dict(xy=g2, z=None)], (0, 0), sigma=0.005, n_starts=16)
    assert abs(r2['weights'][1] / w2 - 1) < 5e-3 and r2['chi2'] < 1e-2


def test_wrong_single_plane_is_rejected():
    m = L.build_sie_shear(P)
    w2 = X.dist_weight(0.25, 2.0, 0.6)
    g1 = _quad(m, (0.04, 0.03))
    g2 = _quad(X.ScaledLens(m, w2), (-0.05, 0.04), half=6.0)
    r = X.fit_multiplane([dict(xy=np.vstack([g1, g2]), z=0.6)], (0, 0), sigma=0.005, n_starts=10)
    assert r['chi2'] > 1e3                                  # two planes cannot be one: image-plane chi2 shows it


def test_reg_and_psf_matrices():
    R = X.reg_matrix(9, 'gradient')
    assert abs(R @ np.ones(81)).max() < 1e-12
    R2 = X.reg_matrix(9, 'curvature')
    gx = np.tile(np.arange(9.0), 9)
    assert abs(R2 @ gx).max() < 1e-12                       # linear ramp has no curvature
    B = X.psf_matrix(X.gaussian_kernel(1.2), (20, 20))
    img = np.zeros((20, 20)); img[10, 10] = 1.0
    out = (B @ img.ravel()).reshape(20, 20)
    assert abs(out.sum() - 1.0) < 1e-12 and out.argmax() == 10 * 20 + 10


def _case(noise, seed=2):
    rng = np.random.default_rng(seed)
    npix, ps = 90, 0.05
    org = ((npix - 1) / 2.0,) * 2
    m = L.build_sie_shear(P)
    src = dict(x0=0.03, y0=0.02, reff=0.12, n=1.0, q=0.7, phi=40.0, amp=1.0)
    im = S.render_arcs(m, (npix, npix), ps, org, (0, 0), src, psf_sigma_pix=1.5, noise=0.0)
    img = im + rng.normal(0, noise, im.shape)
    return m, img, org, ps, src


def test_source_inversion_recovers_source_and_evidence_prefers_truth():
    noise = 0.01
    m, img, org, ps, src = _case(noise)
    mask = X.arc_mask(img, noise, org, 2.4 * P['theta_E'] / ps)
    k = X.gaussian_kernel(1.5)
    SI = X.SourceInversion(m, img, noise, mask, ps, org, (0, 0), kernel=k, n=30)
    r = SI.best()
    neff = SI.n_eff(r)
    assert 0.8 < r['chi2'] / (SI.ndata - neff) < 1.4
    g = SI.grid
    xs = g['x0'] + g['pix'] * np.arange(g['n'])
    BX, BY = np.meshgrid(xs, xs)
    ts = S.sersic_source(BX, BY, src['x0'], src['y0'], src['reff'], 1.0, src['q'], src['phi'], 1.0)
    s, mod = SI.images(r)
    cover = np.asarray((X.lens_operator(m, img.shape, ps, org, (0, 0), g)[np.flatnonzero(mask.ravel())] > 0).sum(axis=0)).reshape(g['n'], g['n']) >= 3
    sel = cover & (ts > 0.05)
    assert np.corrcoef(s[sel], ts[sel])[0, 1] > 0.95
    for dth in (1.03, 0.97):
        p2 = dict(P, theta_E=P['theta_E'] * dth)
        S2 = X.SourceInversion(L.build_sie_shear(p2), img, noise, mask, ps, org, (0, 0), kernel=k, n=30)
        assert S2.best()['evidence'] < r['evidence'] - 20


def test_fit_lens_pixels_recovers_parameters():
    noise = 0.01
    m, img, org, ps, src = _case(noise)
    mask = X.arc_mask(img, noise, org, 2.4 * P['theta_E'] / ps)
    st = dict(P, theta_E=1.33, q=0.85)
    f = X.fit_lens_pixels(img, noise, mask, ps, org, (0, 0), st, free=['theta_E', 'q'], kernel=X.gaussian_kernel(1.5), n=24, maxiter=60)
    assert abs(f['params']['theta_E'] - P['theta_E']) < 0.01 and abs(f['params']['q'] - P['q']) < 0.05


def _cli(*args):
    p = subprocess.run([sys.executable, CLI] + [str(a) for a in args], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    return p


def test_cli_source_inversion(tmp_path):
    w = str(tmp_path)
    _cli('--task', 'simulate', '--work', w, '--seed', 3, '--sim-noise', 0.01)
    truth = json.load(open(os.path.join(w, 'sim_truth.json')))
    _cli('--task', 'fit', '--work', w, '--image', os.path.join(w, 'sim_lens.fits'), '--catalog', os.path.join(w, 'sim_catalog.tsv'),
         '--image-numbers', ','.join(str(n) for n in truth['image_numbers']), '--lens-number', 1, '--sigma-pos', 0.005)
    p = _cli('--task', 'source', '--work', w, '--image', os.path.join(w, 'sim_lens.fits'), '--source-method', 'inversion', '--psf-sigma', 1.5, '--source-n', 30)
    assert 'source inversion' in p.stdout
    rep = json.load(open(os.path.join(w, 'lens_model.json')))['source_reconstruction']
    assert rep['method'] == 'inversion' and rep['n_eff'] > 10 and np.isfinite(rep['chi2_red']) and rep['lam'] > 0
    for f in ('lens_source.fits', 'lens_model_image.fits', 'lens_residual.fits'):
        assert os.path.getsize(os.path.join(w, f)) > 1000


def test_cli_multiplane(tmp_path):
    m = L.build_sie_shear(P)
    w2 = X.dist_weight(0.25, 2.0, 0.6)
    g1 = _quad(m, (0.04, 0.03))
    g2 = _quad(X.ScaledLens(m, w2), (-0.05, 0.04), half=6.0)
    sc = 0.05
    c = 100.0                                                   # 1-based lens pixel
    fmt = lambda g: ';'.join('%.5f,%.5f' % (c + x / sc, c + y / sc) for x, y in g)
    p = _cli('--task', 'multiplane', '--work', str(tmp_path), '--pixscale', sc, '--lens-x', c, '--lens-y', c, '--mp-groups', fmt(g1) + '|' + fmt(g2),
             '--mp-z', '0.6,2.0', '--z-lens', 0.25, '--sigma-pos', 0.005)
    r = json.load(open(os.path.join(str(tmp_path), 'lens_multiplane.json')))
    assert abs(r['params']['theta_E'] - P['theta_E']) < 5e-3 and abs(r['weights'][1] - w2) < 1e-9
    assert r['chi2'] < 1.0
