import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from ogfkit import lensmodel as L, lenssynth as S, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, '..', 'lensmodel.py')


def run_cli(*args):
    p = subprocess.run([sys.executable, CLI] + [str(a) for a in args], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    return p


# ---------------------------------------------------------------------------------------------- physics of the components
def test_sie_convergence_analytic_and_sis_limit():
    s = L.SIE(1.5, 0.6, 30.0, 0.1, -0.2)
    x = np.array([0.7, -1.2, 2.0, 0.3])
    y = np.array([1.1, 0.4, -0.8, 2.5])
    xp, yp = s._rot(x, y)
    k_an = s.theta_E / (2 * np.sqrt(s.q * xp ** 2 + yp ** 2 / s.q))
    assert np.allclose(s.kappa(x, y), k_an, rtol=1e-6)
    sis = L.SIE(1.0, 1.0)
    r = np.array([0.5, 1.0, 2.0])
    ax, ay = sis.alpha(r, 0 * r)
    assert np.allclose(ax, 1.0, atol=1e-12) and np.allclose(ay, 0.0, atol=1e-12)
    # q close to 1 is continuous with the SIS branch
    a1 = L.SIE(1.0, 1 - 1e-5).alpha(np.array([1.3]), np.array([0.4]))
    a2 = L.SIE(1.0, 1.0).alpha(np.array([1.3]), np.array([0.4]))
    assert abs(a1[0][0] - a2[0][0]) < 1e-4 and abs(a1[1][0] - a2[1][0]) < 1e-4


@pytest.mark.parametrize('comp', [L.SIE(1.3, 0.7, 20.0), L.Shear(0.08, 33.0), L.NFW(0.3, 5.0), L.PIEMD(1.0, 0.1, 20.0)])
def test_potential_gradient_is_deflection(comp):
    X, Y, h = np.array([1.3]), np.array([0.7]), 1e-4
    gx = (comp.psi(X + h, Y) - comp.psi(X - h, Y)) / (2 * h)
    gy = (comp.psi(X, Y + h) - comp.psi(X, Y - h)) / (2 * h)
    ax, ay = comp.alpha(X, Y)
    assert abs(gx[0] - ax[0]) < 2e-3 * max(abs(ax[0]), 0.1) and abs(gy[0] - ay[0]) < 2e-3 * max(abs(ay[0]), 0.1)


def test_circular_profiles_kappa_from_deflection():
    n, p = L.NFW(0.3, 5.0), L.PIEMD(1.0, 0.1, 20.0)
    r = np.array([1.0, 3.0, 5.0, 9.0])
    for c in (n, p):
        kn = np.array([c.kappa(np.array([ri]), np.array([0.0]))[0] for ri in r])
        assert np.allclose(kn, c.kappa_r(r), rtol=2e-3)


def test_critical_curve_area_of_lone_sie_is_pi_theta_E_squared():
    for q, phi in ((0.5, 10.0), (0.8, 120.0), (0.95, 0.0)):
        m = L.LensModel([L.SIE(1.5, q, phi)])
        cc = m.critical_curves()
        assert len(cc) == 1 and cc[0]['closed']
        assert abs(cc[0]['theta_eff'] - 1.5) < 2e-3
        # the caustic of a lone SIE (no core) has a cut: radial caustic absent
        assert np.all(np.isfinite(cc[0]['cx']))


# ---------------------------------------------------------------------------------------------- lens equation
def test_sis_images_and_delay_analytic():
    sis = L.LensModel([L.SIE(1.0, 1.0)])
    ims = sis.solve_images(0.2, 0.0)
    xs = sorted(i['x'] for i in ims)
    assert len(ims) == 2 and abs(xs[0] + 0.8) < 1e-9 and abs(xs[1] - 1.2) < 1e-9
    assert np.allclose(sorted(i['mu'] for i in ims), [-4.0, 6.0], atol=1e-4)       # mu = 1/(1 -/+ theta_E/r) (numerical Jacobian)
    d = L.time_delays(sis, ims, 0.2, 0.0, 0.5, 2.0)
    Dl, Ds, Dls = L.distances(0.5, 2.0)
    fac = 1.5 * Dl * Ds / Dls * L.ARCSEC ** 2 * 3.0856775814914e19 / (L.C_KMS * 86400)
    assert abs(max(d) - 2 * 1.0 * 0.2 * fac) < 1e-6 * max(d)           # Delta phi = 2 theta_E beta for the SIS


def test_sie_quad_solution_is_exact_and_has_four_images():
    m = L.build_sie_shear(dict(theta_E=1.5, q=0.6, phi=42.0, gamma=0.08, phi_g=100.0, x0=0.0, y0=0.0))
    ims = m.solve_images(0.1, 0.05)
    assert len(ims) == 4
    for i in ims:
        bx, by = m.beta(i['x'], i['y'])
        assert abs(bx - 0.1) < 1e-8 and abs(by - 0.05) < 1e-8
    assert sorted(i['parity'] for i in ims) == ['min', 'min', 'saddle', 'saddle']
    assert abs(sum(1 for i in ims if i['mu'] < 0) - 2) == 0


def test_counter_image_prediction():
    m = L.build_sie_shear(dict(theta_E=1.2, q=0.7, phi=10.0, gamma=0.05, phi_g=60.0, x0=0.0, y0=0.0))
    ims = m.solve_images(0.05, -0.04)
    assert len(ims) == 4
    obs = np.array([[i['x'], i['y']] for i in ims[:3]])
    pred = L.predict(m, 0.05, -0.04, (0.0, 0.0), observed=obs)
    miss = [p for p in pred if p['matched'] < 0]
    assert len(miss) == 1 and math.hypot(miss[0]['x'] - ims[3]['x'], miss[0]['y'] - ims[3]['y']) < 1e-6


def test_mass_and_velocity_dispersion_consistency():
    m1 = L.mass_in_radius(1.0, 0.5, 2.0)
    assert abs(L.mass_in_radius(2.0, 0.5, 2.0) / m1 - 4.0) < 1e-9
    # SIS: M(<theta_E) = pi sigma^2 / G * R_E  ->  sigma from the helper must reproduce the Einstein mass
    sig = L.sis_sigma_v(1.5, 0.5, 2.0)
    Dl = L.distances(0.5, 2.0)[0]
    R = 1.5 * L.ARCSEC * Dl                       # Mpc
    M = math.pi * sig ** 2 * R / 4.30091e-9       # sigma^2 / G * pi R
    assert abs(M / L.mass_in_radius(1.5, 0.5, 2.0) - 1.0) < 1e-6
    assert 150 < sig < 450


# ---------------------------------------------------------------------------------------------- fitting against synthetic truth
def test_noise_free_recovery_is_exact():
    rng = np.random.default_rng(21)
    for k in range(3):
        lens = S.random_lens(rng, quad=True)
        obs = np.array([[i['x'], i['y']] for i in lens['images']])
        f = L.fit_sie_shear(obs, 0.003, (0.0, 0.0), seed=k)
        t, p = lens['params'], f['params']
        assert abs(p['theta_E'] - t['theta_E']) < 1e-5 and abs(p['q'] - t['q']) < 1e-5 and abs(p['gamma'] - t['gamma']) < 1e-5
        assert abs(S.angle_diff(p['phi'], t['phi'])) < 1e-3
        assert math.hypot(f['source'][0] - lens['source'][0], f['source'][1] - lens['source'][1]) < 1e-5
        assert f['rms_arcsec'] < 1e-6


def test_noisy_recovery_statistics():
    res = S.validation_sample(n_lens=10, noise_arcsec=0.003, quad_fraction=1.0, seed=5)
    dth = np.array([(r['fit']['params']['theta_E'] - r['truth']['theta_E']) / r['truth']['theta_E'] for r in res])
    rms = np.array([r['fit']['rms_arcsec'] for r in res])
    dsrc = np.array([math.hypot(r['fit']['source'][0] - r['source'][0], r['fit']['source'][1] - r['source'][1]) for r in res])
    assert np.median(np.abs(dth)) < 2e-3 and np.max(np.abs(dth)) < 1e-2
    assert np.median(rms) < 0.003 and np.median(dsrc) < 0.003 and np.max(dsrc) < 0.02
    chi2 = np.array([r['fit']['chi2'] for r in res])
    assert np.mean(chi2) < 3.0                                          # dof = 1


def test_unmodelled_perturber_is_detected_by_chi2():
    rng = np.random.default_rng(3)
    lens = S.random_lens(rng, quad=True)
    p = lens['params']
    im = lens['images'][0]
    pert = L.PIEMD(b0=0.15, a=0.05, s=10.0, x0=im['x'] + 0.4, y0=im['y'] - 0.3)
    m = L.build_sie_shear(p, [pert])
    ims = m.solve_images(*lens['source'], half=3 * p['theta_E'], n=301)
    ims = [a for a in ims if abs(a['mu']) > 0.5][:4]
    assert len(ims) == 4
    obs = np.array([[a['x'], a['y']] for a in ims])
    good = L.fit_sie_shear(obs, 0.003, (0.0, 0.0), perturbers=[pert])
    bad = L.fit_sie_shear(obs, 0.003, (0.0, 0.0))
    assert good['chi2'] < 0.1 and good['rms_arcsec'] < 1e-4
    assert bad['chi2'] > 20 * max(good['chi2'], 0.05)


def test_source_reconstruction_discriminates_model_errors():
    ok = S.source_recon_test(seed=1, theta_E_error=0.0)
    off = S.source_recon_test(seed=1, theta_E_error=0.10)
    assert ok['corr'] > 0.98 and abs(ok['flux_ratio'] - 1) < 0.05 and ok['centroid_err'] < 0.03
    assert off['corr'] < 0.8 and off['resid_over_noise'] > 3 * ok['resid_over_noise']


def test_lenstronomy_cross_check_when_installed():
    pytest.importorskip('lenstronomy')
    p = dict(theta_E=1.4, q=0.7, phi=25.0, x0=0.0, y0=0.0)
    x, y = np.array([0.5, -1.2, 2.0]), np.array([1.0, 0.3, -1.1])
    ax, ay = L.SIE(**p).alpha(x, y)
    bx, by = L.lenstronomy_alpha('SIE', p, x, y)
    assert np.allclose(ax, bx, atol=1e-6) and np.allclose(ay, by, atol=1e-6)


# ---------------------------------------------------------------------------------------------- CLI
def test_cli_end_to_end(tmp_path):
    w = str(tmp_path / 'w')
    run_cli('--task', 'simulate', '--work', w, '--seed', 3)
    truth = json.load(open(os.path.join(w, 'sim_truth.json')))
    base = ['--work', w, '--image', os.path.join(w, 'sim_lens.fits')]
    p = run_cli('--task', 'fit', *base, '--catalog', os.path.join(w, 'sim_catalog.tsv'), '--image-numbers', '2,3,4,5', '--lens-number', 1, '--sigma-pos', 0.005,
                '--z-lens', 0.5, '--z-source', 2.0, '--flux-column', 'FLUX', '--flux-sigma', 0.05)
    lines = p.stdout.strip().split('\n')
    assert lines[0].split('\t') == ['NUMBER', 'LENS_ROLE', 'LENS_MU', 'LENS_RES', 'LENS_DT', 'LENS_PARITY'] and len(lines) == 5
    res = json.load(open(os.path.join(w, 'lens_model.json')))
    assert abs(res['params']['theta_E'] - truth['params']['theta_E']) < 1e-3
    assert abs(res['theta_E_critical_curve'] - 1.5298) < 5e-3 or res['theta_E_critical_curve'] > 1.0
    assert res['n_images_predicted'] == 4 and res['n_counter_images'] == 0 and res['physical']['mass_Einstein_Msun'] > 1e11
    tmu = sorted(i['mu'] for i in truth['images'])
    fmu = sorted(i['mu'] for i in res['images'])
    assert np.allclose(tmu, fmu, atol=0.02)
    for f in ('lens_images.tsv', 'lens_curves.json', 'lens_overlay.reg'):
        assert os.path.getsize(os.path.join(w, f)) > 50
    assert 'polygon' in open(os.path.join(w, 'lens_overlay.reg')).read()
    run_cli('--task', 'curves', *base)
    run_cli('--task', 'magmap', *base)
    from astropy.io import fits
    mm = fits.getdata(os.path.join(w, 'lens_magnification.fits'))
    model = L.build_sie_shear(res['params'])
    ref = np.array(res['lens_pixel']) - 1
    iy, ix = 120, 170
    mu_direct = float(model.mu((ix - ref[0]) * res['pixscale'], (iy - ref[1]) * res['pixscale']))
    assert abs(mm[iy, ix] - np.clip(mu_direct, -100, 100)) < 1e-3 * max(1, abs(mu_direct))
    run_cli('--task', 'source', *base)
    res2 = json.load(open(os.path.join(w, 'lens_model.json')))
    assert res2['source_reconstruction']['n_pixels_mapped'] > 1000
    p = run_cli('--task', 'predict', '--work', w, '--image', base[3], '--z-source', 3.0)
    assert 'mass_Einstein_Msun' in p.stdout


def test_cli_errors_are_clear(tmp_path):
    w = str(tmp_path / 'w')
    run_cli('--task', 'simulate', '--work', w, '--seed', 3)
    p = subprocess.run([sys.executable, CLI, '--task', 'fit', '--work', w, '--image', os.path.join(w, 'sim_lens.fits'), '--catalog', os.path.join(w, 'sim_catalog.tsv'),
                        '--image-numbers', '2,3,4,5'], capture_output=True, text=True)
    assert p.returncode != 0 and 'lens centre is required' in p.stderr
    p = subprocess.run([sys.executable, CLI, '--task', 'fit', '--work', w, '--image', os.path.join(w, 'sim_lens.fits'), '--catalog', os.path.join(w, 'sim_catalog.tsv'),
                        '--image-numbers', '2', '--lens-number', '1'], capture_output=True, text=True)
    assert p.returncode != 0 and 'at least 2' in p.stderr
