"""Synthetic checks for the cluster weak-lensing engine and its command line.

The catalogues are drawn from this module's own NFW model. A recovery here is
that the fitter returns the mass it was given. It is not a survey calibration.
"""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from ogfkit import wlcluster as wl

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'wlcluster.py')


@pytest.fixture(scope='module')
def cosmo():
    return wl.Cosmology(H0=70.0, Om=0.3, ns=0.96, sigma8=0.8)


def test_sigma_crit_matches_astropy_and_growth_falls(cosmo):
    from astropy.cosmology import FlatLambdaCDM
    ref = FlatLambdaCDM(H0=70.0, Om0=0.3)
    pref = (wl.C_KMS ** 2) / (4.0 * math.pi * wl.G_MPC)
    for zl, zs in ((0.3, 1.0), (0.5, 1.5), (0.2, 0.8)):
        Dl = ref.angular_diameter_distance(zl).value
        Ds = ref.angular_diameter_distance(zs).value
        try:
            Dls = ref.angular_diameter_distance(zl, zs).value
        except TypeError:
            Dls = ref.angular_diameter_distance_z1z2(zl, zs).value
        ours = float(cosmo.sigma_crit(zl, zs))
        other = pref * Ds / (Dl * Dls)
        assert abs(ours - other) / other < 5e-4
    assert math.isnan(float(np.asarray(cosmo.sigma_crit(0.3, 0.3))))
    assert cosmo.growth(0.0) == pytest.approx(1.0)
    assert cosmo.growth(0.0) > cosmo.growth(0.3) > cosmo.growth(1.0)
    # Duffy Table 1 at the pivot mass and z=0.
    assert cosmo.duffy_c(wl.DUFFY_PIVOT / cosmo.h, 0.0) == pytest.approx(wl.DUFFY_A)


def test_nfw_mass_and_lensmodel_kappa(cosmo):
    M, c, z = 1.0e15, 4.0, 0.3
    rho_s, rs, R200 = wl.nfw_scales(M, c, z, cosmo)
    assert float(wl.mass_enclosed(R200, rho_s, rs)) / M == pytest.approx(1.0, rel=1e-6)
    R = np.array([0.3, 0.7, 1.0, 1.8, 2.5])
    Sig, dSig = wl.nfw_sigma(R, rho_s, rs)
    sc = float(cosmo.sigma_crit(z, 1.0))
    from ogfkit.lensmodel import NFW
    nfw = NFW(rho_s * rs / sc, rs)
    kap = nfw.kappa_r(R)
    kbar = nfw.alpha_r(R) / R
    assert np.allclose(Sig / sc, kap, rtol=1e-5)
    assert np.allclose(dSig / sc, kbar - kap, rtol=1e-5)


def test_linear_spectrum_and_bias(cosmo):
    k = np.array([1e-8, 1e-4])
    T = wl.transfer_eh98(k, 70.0, 0.3)
    assert float(T[0]) > 0.999
    field = cosmo.linear_field()
    assert field.sigma8() == pytest.approx(0.8, rel=1e-8)
    assert np.all(np.isfinite(field.xi(np.array([1.0, 8.0, 30.0]))))
    b_low = wl.halo_bias(cosmo, 1.0e12, 0.3)
    b_high = wl.halo_bias(cosmo, 1.0e15, 0.3)
    assert b_high > 1.0
    assert abs(b_low - 1.0) < abs(b_high - 1.0)


def test_two_halo_is_smaller_inside_and_grows_outside(cosmo):
    R = np.geomspace(0.3, 10.0, 16)
    rho_s, rs, _ = wl.nfw_scales(1.0e15, 4.0, 0.3, cosmo)
    _S, d1 = wl.nfw_sigma(R, rho_s, rs)
    model = wl.HaloModel(cosmo, 0.3, R, two_halo=True)
    _St, dT = model.profiles(1.0e15, 4.0)
    d2 = dT - d1
    i1 = int(np.argmin(np.abs(R - 1.0)))
    i9 = int(np.argmin(np.abs(R - 9.0)))
    assert d2[i1] > 0 and d1[i1] > d2[i1]
    frac = d2 / dT
    assert frac[i9] > frac[i1]


def test_tangential_sign_and_inverse():
    gt, gx = wl.tangential(np.array([10.0]), np.array([0.0]), np.array([-0.1]), np.array([0.0]))
    assert float(gt[0]) == pytest.approx(0.1)
    assert float(gx[0]) == pytest.approx(0.0, abs=1e-12)
    phi = np.array([0.3, 1.1, 2.4])
    g1, g2 = wl.components_from_tangential(phi, np.array([0.07, -0.02, 0.04]), np.array([0.01, 0.0, -0.03]))
    gt2, gx2 = wl.tangential(np.cos(phi), np.sin(phi), g1, g2)
    assert np.allclose(gt2, [0.07, -0.02, 0.04])
    assert np.allclose(gx2, [0.01, 0.0, -0.03])


def test_foreground_source_is_dropped(cosmo):
    sel = wl.select_sources(
        np.array([10.0, 20.0]), np.array([0.0, 0.0]), np.zeros(2), np.zeros(2),
        np.array([0.05, 1.0]), np.ones(2), 0.3, cosmo, z_gap=0.1, rmin=0.01, rmax=5.0)
    assert list(sel['index']) == [1]


def _draw(cosmo, n, rmin, rmax, M, c, z_l, z_s, sigma, seed, reduced=True, m_obs=0.0):
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    R = np.sqrt(rmin ** 2 + u * (rmax ** 2 - rmin ** 2))
    phi = rng.random(n) * 2.0 * math.pi
    grid = np.geomspace(rmin * 0.8, rmax * 1.2, 80)
    Sig, dSig = wl.HaloModel(cosmo, z_l, grid, two_halo=False).profiles(M, c)
    sc = np.full(n, float(cosmo.sigma_crit(z_l, z_s)))
    g = wl.reduced_shear(np.interp(R, grid, Sig), np.interp(R, grid, dSig), sc, reduced=reduced)
    g1, g2 = wl.components_from_tangential(phi, g, np.zeros(n))
    if sigma:
        g1 = g1 + rng.normal(0.0, sigma, n)
        g2 = g2 + rng.normal(0.0, sigma, n)
    g1 = g1 * (1.0 + m_obs)
    g2 = g2 * (1.0 + m_obs)
    Dl = float(cosmo.angular_diameter(z_l))
    ang = R / Dl
    x = ang * np.cos(phi) / wl.ARCSEC
    y = ang * np.sin(phi) / wl.ARCSEC
    return x, y, g1, g2


def _fit(cosmo, x, y, g1, g2, z_l, z_s, rmin, rmax, m=0.0, **kw):
    n = len(x)
    sel = wl.select_sources(x, y, g1, g2, np.full(n, z_s), np.ones(n), z_l, cosmo,
                            z_gap=0.05, rmin=rmin, rmax=rmax, m=m)
    edges = wl.bin_edges(rmin, rmax, kw.pop('nbins', 8))
    kw.setdefault('two_halo', False)
    kw.setdefault('samples', 0)
    kw.setdefault('sigma_e', 1e-3)
    kw.setdefault('concentration', 'fixed')
    kw.setdefault('c_fixed', 4.0)
    kw.setdefault('reduced', True)
    return wl.fit_mass(sel, edges, cosmo, z_l, **kw)


def test_noiseless_mass_recovers_the_input(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 700, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.0, 2)
    summary, rows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=20, seed=2)
    assert summary['success']
    assert abs(math.log10(summary['M200c_Msun'] / M)) < 0.005
    assert summary['c200'] == pytest.approx(4.0, rel=1e-3)
    assert rows[0]['boost'] == pytest.approx(1.0)


def test_free_concentration_stays_near_the_input(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 900, 0.25, 3.0, M, 4.0, 0.3, 1.0, 0.0, 3)
    summary, _rows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.25, 3.0, concentration='free', sigma_e=1e-4)
    assert summary['success']
    assert abs(math.log10(summary['M200c_Msun'] / M)) < 0.03
    assert abs(summary['c200'] - 4.0) < 0.4


def test_noisy_mass_is_in_a_wide_window_and_cross_is_small(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 1600, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.25, 7)
    summary, rows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, sigma_e=0.25, samples=40, seed=7)
    assert summary['success']
    assert abs(math.log10(summary['M200c_Msun'] / M)) < 0.30
    assert abs(float(np.mean([r['gx'] for r in rows]))) < 0.04


def test_multiplicative_bias_shifts_the_mass(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 600, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.0, 4, m_obs=0.05)
    raw, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, m=0.0)
    fixed, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, m=0.05)
    assert math.log10(raw['M200c_Msun']) > math.log10(fixed['M200c_Msun']) + 0.01
    assert abs(math.log10(fixed['M200c_Msun'] / M)) < 0.015


def test_ignoring_reduced_shear_raises_the_mass(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 500, 0.15, 1.2, M, 4.0, 0.3, 1.0, 0.0, 5, reduced=True)
    on, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.15, 1.2, reduced=True, nbins=6)
    off, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.15, 1.2, reduced=False, nbins=6)
    assert abs(math.log10(on['M200c_Msun'] / M)) < 0.02
    assert math.log10(off['M200c_Msun']) > math.log10(on['M200c_Msun']) + 0.02


def test_miscentring_lowers_the_inner_excess(cosmo):
    R = np.geomspace(0.12, 3.0, 18)
    plain = wl.HaloModel(cosmo, 0.3, R, two_halo=False, miscenter=False)
    off = wl.HaloModel(cosmo, 0.3, R, two_halo=False, miscenter=True, f_cen=0.4, sigma_off=0.4)
    _s0, d0 = plain.profiles(1.0e15, 4.0)
    _s1, d1 = off.profiles(1.0e15, 4.0)
    assert d1[0] < d0[0]


def test_boost_raises_a_diluted_inner_bin(cosmo):
    z_l, z_s = 0.3, 1.0
    x, y, g1, g2 = _draw(cosmo, 400, 0.3, 3.0, 1.0e15, 4.0, z_l, z_s, 0.0, 6)
    # Extra zero-shear objects in the inner radial range, as cluster members would.
    rng = np.random.default_rng(8)
    u = rng.random(250)
    Rm = np.sqrt(0.3 ** 2 + u * (0.55 ** 2 - 0.3 ** 2))
    phi = rng.random(250) * 2.0 * math.pi
    Dl = float(cosmo.angular_diameter(z_l))
    ang = Rm / Dl
    xm = ang * np.cos(phi) / wl.ARCSEC
    ym = ang * np.sin(phi) / wl.ARCSEC
    x = np.concatenate([x, xm])
    y = np.concatenate([y, ym])
    g1 = np.concatenate([g1, np.zeros(250)])
    g2 = np.concatenate([g2, np.zeros(250)])
    sel = wl.select_sources(x, y, g1, g2, np.full(len(x), z_s), np.ones(len(x)), z_l, cosmo,
                            z_gap=0.05, rmin=0.3, rmax=3.0)
    ur = rng.random(2000)
    random_R = np.sqrt(0.3 ** 2 + ur * (3.0 ** 2 - 0.3 ** 2))
    edges = wl.bin_edges(0.3, 3.0, 6)
    bare, _b, _B = wl._profile_from_selection(sel, edges, 0.2, None)
    boosted, _b2, _B2 = wl._profile_from_selection(sel, edges, 0.2, random_R)
    assert boosted[0]['boost'] > 1.2
    assert boosted[0]['deltasigma'] > bare[0]['deltasigma']


def test_effective_redshift_of_a_single_slice(cosmo):
    z = np.array([0.4, 1.0, 1.6])
    n = np.array([0.0, 5.0, 0.0])
    assert wl.effective_source_z(cosmo, 0.3, z, n, 0.1) == pytest.approx(1.0, abs=0.01)


def _gaussian_scene(g1, g2, with_psf):
    n = 49
    cx = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    gal = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cx) ** 2) / 4.0 ** 2)
    sheared = wl.shear_image(gal, g1, g2)
    if not with_psf:
        return sheared, None
    psf = wl.gaussian_psf(n, 1.8)
    from scipy.signal import fftconvolve
    return fftconvolve(sheared, psf, mode='same'), psf


def test_metacal_recovers_a_noiseless_shear():
    obs, psf = _gaussian_scene(0.05, -0.02, True)
    rec = wl.metacal_stamp(obs, psf)
    assert rec['flag'] == 0
    assert rec['g1'] == pytest.approx(0.05, rel=0.08, abs=0.004)
    assert rec['g2'] == pytest.approx(-0.02, rel=0.08, abs=0.004)
    rnd, psf = _gaussian_scene(0.0, 0.0, True)
    rec_r = wl.metacal_stamp(rnd, psf)
    assert abs(rec_r['g1']) < 1e-3 and abs(rec_r['g2']) < 1e-3
    bare, _ = _gaussian_scene(0.05, 0.0, False)
    rec_b = wl.metacal_stamp(bare, None)
    assert rec_b['flag'] == 0
    assert rec_b['g1'] == pytest.approx(0.05, rel=0.08, abs=0.004)


def test_unresolved_and_edge_stamps_are_flagged():
    n = 49
    cx = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    tiny = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cx) ** 2) / 0.35 ** 2)
    psf = wl.gaussian_psf(n, 1.8)
    from scipy.signal import fftconvolve
    rec = wl.metacal_stamp(fftconvolve(tiny, psf, mode='same'), psf)
    assert rec['flag'] != 0
    image = np.zeros((40, 40))
    out = wl.measure_stamps(image, [1.0], [1.0], psf_fwhm_pix=0.0, stamp=32)
    assert out[0]['flag'] == 1


def _run(args):
    return subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)


def test_cli_missing_catalog_writes_nothing(tmp_path):
    work = tmp_path / 'wlcluster'
    p = _run(['--task', 'mass', '--catalog', str(tmp_path / 'absent.tsv'), '--work', str(work), '--z-lens', '0.3'])
    assert p.returncode == 1
    assert 'was not found' in p.stderr
    assert p.stdout.strip() == ''
    assert not (work / 'wl_mass.json').exists()
    assert not (work / 'wl_profile.tsv').exists()


def _write_catalog(path, x, y, g1, g2, z=1.0):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tWL_G1\tWL_G2\tZ\n')
        for i in range(len(x)):
            # 1-based pixels. Centre will be (1000, 1000) and the scale 1 arcsec/px,
            # so a catalogue pixel is 1000 + arcsec.
            f.write('%d\t%.6f\t%.6f\t%.8f\t%.8f\t%.4f\n' % (i + 1, 1000.0 + x[i], 1000.0 + y[i], g1[i], g2[i], z))


def test_cli_recovers_a_noiseless_mass_and_keeps_a_failed_fit(tmp_path, cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 500, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.0, 9)
    cat = tmp_path / 'src.tsv'
    _write_catalog(cat, x, y, g1, g2)
    work = tmp_path / 'run'
    p = _run([
        '--task', 'mass', '--catalog', str(cat), '--work', str(work),
        '--coord', 'xy', '--center-x', '1000', '--center-y', '1000', '--pixel-scale', '1',
        '--z-lens', '0.3', '--z-source', '1', '--z-column', '',
        '--concentration', 'fixed', '--c-fixed', '4', '--no-two-halo',
        '--power', 'linear', '--no-lss-cov',
        '--sigma-e', '0.001', '--samples', '20', '--seed', '9', '--min-per-bin', '8'])
    assert p.returncode == 0, p.stderr
    assert 'WL_GT' in p.stdout.splitlines()[0]
    doc = json.loads((work / 'wl_mass.json').read_text())
    assert doc['success']
    assert abs(math.log10(doc['M200c_Msun'] / M)) < 0.02
    assert doc['boost'] == 'not supplied'
    header = (work / 'wl_profile.tsv').read_text().splitlines()[0].split('\t')
    assert 'DELTASIGMA' in header and 'GT_MODEL' in header
    # Too few sources: the profile is still written and the process exits 0.
    small = tmp_path / 'small.tsv'
    _write_catalog(small, x[:5], y[:5], g1[:5], g2[:5])
    work2 = tmp_path / 'thin'
    p2 = _run([
        '--task', 'mass', '--catalog', str(small), '--work', str(work2),
        '--coord', 'xy', '--center-x', '1000', '--center-y', '1000', '--pixel-scale', '1',
        '--z-lens', '0.3', '--z-source', '1', '--concentration', 'fixed', '--c-fixed', '4',
        '--no-two-halo', '--power', 'linear', '--no-lss-cov', '--samples', '0'])
    assert p2.returncode == 0, p2.stderr
    doc2 = json.loads((work2 / 'wl_mass.json').read_text())
    assert doc2['success'] is False


def test_cli_shape_on_one_stamp(tmp_path):
    n = 49
    cx = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    gal = np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cx) ** 2) / 4.0 ** 2)
    psf = wl.gaussian_psf(n, 1.6)
    from scipy.signal import fftconvolve
    cut = fftconvolve(wl.shear_image(gal, 0.04, 0.0), psf, mode='same')
    image = np.zeros((80, 80), np.float32)
    image[20:20 + n, 30:30 + n] = cut
    # Object pixel is the centre of the pasted cut.
    ix = 30 + int(round(cx))
    iy = 20 + int(round(cx))
    fits = tmp_path / 'stamp.fits'
    from ogfkit import imageio
    imageio.save_fits(str(fits), image)
    cat = tmp_path / 'one.tsv'
    # 1-based catalogue position.
    cat.write_text('NUMBER\tX_IMAGE\tY_IMAGE\n1\t%d\t%d\n' % (ix + 1, iy + 1))
    fwhm = 1.6 * 2.0 * math.sqrt(2.0 * math.log(2.0))
    p = _run([
        '--task', 'shapes', '--catalog', str(cat), '--image', str(fits),
        '--work', str(tmp_path / 'sh'), '--psf-fwhm', '%.6f' % fwhm, '--stamp', '48'])
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert lines[0].split('\t')[:5] == ['NUMBER', 'WL_E1', 'WL_E2', 'WL_G1', 'WL_G2']
    cols = lines[0].split('\t')
    vals = lines[1].split('\t')
    g1 = float(vals[cols.index('WL_G1')])
    flag = int(float(vals[cols.index('WL_FLAG')]))
    assert flag == 0
    assert g1 == pytest.approx(0.04, rel=0.15, abs=0.008)


def test_cli_sky_east_is_left():
    import importlib.util
    spec = importlib.util.spec_from_file_location('wlcluster_cli', CLI)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    class A:
        pass
    a = A()
    a.coord = 'sky'
    a.ra = 10.0
    a.dec = 20.0
    a.ra_column = 'RA'
    a.dec_column = 'DEC'
    rows = [{'RA': '10.01', 'DEC': '20'}]
    x, y = mod.plane_arcsec(a, rows, ['RA', 'DEC'])
    assert float(x[0]) < 0
    assert float(y[0]) == pytest.approx(0.0, abs=1e-6)
    expected = -0.01 * math.cos(math.radians(20.0)) * 3600.0
    assert float(x[0]) == pytest.approx(expected, rel=1e-6)


def test_halofit_parameters_and_two_halo_under_the_nfw(cosmo):
    field = cosmo.linear_field()
    k_sigma, n_eff, curv = field.halofit_params(0.0)
    assert 0.05 < k_sigma < 5.0
    assert -2.8 < n_eff < -0.8
    assert abs(curv) < 5.0
    k, pnl, _info = field.halofit_pk(0.0)
    pl = float(field.linear_pk(np.array([1.0]), 0.0)[0])
    assert float(np.interp(1.0, k, pnl)) > pl
    R = np.geomspace(0.3, 10.0, 16)
    rho_s, rs, _ = wl.nfw_scales(1.0e15, 4.0, 0.3, cosmo)
    _S, d1 = wl.nfw_sigma(R, rho_s, rs)
    model = wl.HaloModel(cosmo, 0.3, R, two_halo=True, power='halofit')
    _St, dT = model.profiles(1.0e15, 4.0)
    d2 = dT - d1
    i1 = int(np.argmin(np.abs(R - 1.0)))
    assert d2[i1] > 0
    assert d1[i1] > d2[i1]
    x, y, g1, g2 = _draw(cosmo, 80, 0.3, 3.0, 1.0e15, 4.0, 0.3, 1.0, 0.0, 11)
    summary, _rows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=0,
                          two_halo=True, power='halofit', sigma_e=0.05, min_per_bin=2)
    text = ' '.join(summary['limits'])
    assert 'Halofit' in text
    assert 'not Halofit' not in text
    assert 'Takahashi' in text


def test_lss_covariance_scale_and_noiseless_recovery(cosmo):
    Dl = float(cosmo.angular_diameter(0.3))
    theta1 = np.array([0.7, 1.4]) / Dl
    theta2 = np.array([1.4, 2.8]) / Dl
    C = wl.lss_bin_covariance(theta1, theta2, cosmo, 1.0, power='linear')
    assert C.shape == (2, 2)
    sig = np.sqrt(np.diag(C))
    assert np.all(sig > 1e-5) and np.all(sig < 0.2)
    assert np.all(np.isfinite(C))
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 700, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.0, 12)
    summary, rows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=0,
                         lss_cov=True, power='linear', sigma_e=0.02)
    assert summary['covariance'] == 'shape+lss'
    assert summary['success']
    assert abs(math.log10(summary['M200c_Msun'] / M)) < 0.02
    assert 'large-scale-structure' in ' '.join(summary['limits']).lower()
    jk, jrows = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=0,
                     n_jackknife=8, sigma_e=0.02)
    assert jk['covariance'] == 'jackknife'
    assert all(r['gt_err'] > 0 for r in jrows)
    assert abs(math.log10(jk['M200c_Msun'] / M)) < 0.03


def test_count_slope_null_matches_and_a_lower_slope_raises_the_mass(cosmo):
    M = 1.0e15
    x, y, g1, g2 = _draw(cosmo, 500, 0.3, 3.0, M, 4.0, 0.3, 1.0, 0.0, 13)
    base, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0)
    null, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, count_slope=0.4)
    low, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, count_slope=0.2)
    assert abs(math.log10(null['M200c_Msun'] / base['M200c_Msun'])) < 1e-6
    assert math.log10(low['M200c_Msun']) > math.log10(base['M200c_Msun']) + 0.005
    assert '0.4' in ' '.join(null['limits'])


def test_tight_miscentring_prior_stays_put(cosmo):
    x, y, g1, g2 = _draw(cosmo, 400, 0.3, 3.0, 1.0e15, 4.0, 0.3, 1.0, 0.0, 14)
    wide, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=0, sigma_e=0.05,
                   fit_miscenter=True, f_cen=0.55, sigma_off=0.25,
                   f_cen_sigma=0.4, sigma_off_sigma=0.3)
    tight, _ = _fit(cosmo, x, y, g1, g2, 0.3, 1.0, 0.3, 3.0, samples=0, sigma_e=0.05,
                    fit_miscenter=True, f_cen=0.55, sigma_off=0.25,
                    f_cen_sigma=0.01, sigma_off_sigma=0.005)
    assert wide['f_cen'] > tight['f_cen'] + 0.05
    assert tight['f_cen'] == pytest.approx(0.55, abs=0.04)
    assert tight['sigma_off'] == pytest.approx(0.25, abs=0.02)
    assert 'prior' in ' '.join(tight['limits']).lower()


def test_selection_response_zero_on_a_kept_stamp_and_nonzero_on_a_cut():
    obs, psf = _gaussian_scene(0.04, 0.0, True)
    rec = wl.metacal_stamp(obs, psf)
    info = wl.selection_response([rec, rec], 0.01, min_resolution=0.05, emax=0.8)
    assert abs(info['selection_R1']) < 1e-8
    assert abs(info['selection_R2']) < 1e-8
    assert info['mean_response'] > 0.2

    def side(e1, resolved=True):
        return dict(e1=e1, e2=0.0, T=12.0, T_psf=2.0, resolved=resolved)

    recs = [
        dict(e1=0.30, e2=0.0, R1=0.7, R2=0.7, flag=0,
             sides={'+1': side(0.2), '-1': side(0.9), '+2': side(0.1), '-2': side(0.1)}),
        dict(e1=0.10, e2=0.0, R1=0.7, R2=0.7, flag=0,
             sides={'+1': side(0.2), '-1': side(0.2), '+2': side(0.1), '-2': side(0.1)}),
    ]
    cut = wl.selection_response(recs, 0.01, emax=0.5)
    # Plus keeps both (mean 0.20). Minus keeps only the second (mean 0.10).
    assert cut['selection_R1'] == pytest.approx((0.20 - 0.10) / 0.02)
    g1, g2, f1, f2 = wl.apply_selection_response(np.array([0.04]), np.array([0.01]), 0.8, 0.2, 0.0)
    assert f1 == pytest.approx(0.8 / 1.0)
    assert f2 == pytest.approx(1.0)
    assert float(g1[0]) == pytest.approx(0.04 * f1)


def test_internal_mcal_on_a_noiseless_gaussian():
    doc = wl.calibrate_shear(kind='gaussian', noise=0.0, psf_sigma=0.0)
    assert doc['success']
    assert abs(doc['m']) < 0.05
    assert abs(doc['c1']) < 0.005
    assert abs(doc['c2']) < 0.005
    assert 'not a survey' in doc['note']
    assert 'Gaussian PSF' in doc['metacal_scope']
    disc = wl.calibrate_shear(kind='sersic', n_sersic=1.0, noise=0.0, psf_sigma=0.0)
    assert disc['success']
    assert math.isfinite(disc['m'])
    assert abs(disc['c1']) < 0.02
    assert abs(disc['c2']) < 0.02


def test_cli_mcal_writes_json_without_a_catalog(tmp_path):
    work = tmp_path / 'mcal'
    p = _run(['--task', 'mcal', '--work', str(work), '--galaxy', 'gaussian',
              '--psf-fwhm', '0', '--stamp', '49', '--seed', '3'])
    assert p.returncode == 0, p.stderr
    doc = json.loads((work / 'wl_mcal.json').read_text())
    assert doc['success']
    assert 'not a survey' in doc['note']
    assert 'Gaussian PSF' in doc['metacal_scope']
    assert abs(doc['m']) < 0.05


def _catalogue(cosmo, n, rmin, rmax, mass, z_lens, z_source, seed, **extra):
    x, y, g1, g2 = _draw(cosmo, n, rmin, rmax, mass, 4.0, z_lens, z_source, 0.0, seed)
    record = dict(id='a', z_lens=z_lens, mass=mass, x_arcsec=x, y_arcsec=y, g1=g1, g2=g2,
                  z=np.full(n, z_source), weight=np.ones(n))
    record.update(extra)
    return record


def _stack(cosmo, records, **kw):
    kw.setdefault('two_halo', False)
    kw.setdefault('samples', 0)
    kw.setdefault('sigma_e', 1e-3)
    kw.setdefault('concentration', 'fixed')
    kw.setdefault('c_fixed', 4.0)
    kw.setdefault('power', 'linear')
    kw.setdefault('lss_cov', False)
    kw.setdefault('min_per_bin', 8)
    return wl.stack_shear_profiles(records, cosmo=cosmo, **kw)


def test_stacked_copies_keep_one_shear_profile(cosmo):
    mass = 1.0e15
    one = _catalogue(cosmo, 360, 0.4, 2.5, mass, 0.3, 1.0, 11)
    other = dict(one, id='b')
    alone = _stack(cosmo, [one], rmin=0.4, rmax=2.5, nbins=4, fit=False, per_cluster_mass=False)
    pair = _stack(cosmo, [one, other], rmin=0.4, rmax=2.5, nbins=4, fit=True, per_cluster_mass=True)
    assert len(pair['profile']) >= 2
    assert len(pair['profile']) == len(alone['profile'])
    for left, right in zip(alone['profile'], pair['profile']):
        assert right['gt'] == pytest.approx(left['gt'], rel=1e-6, abs=1e-8)
        assert right['n_clusters'] == 2
        assert right['deltasigma'] == pytest.approx(right['gt'] / right['sigma_crit_inv'] / 1e12, rel=1e-6)
    assert pair['success']
    assert abs(math.log10(pair['M200c_Msun'] / mass)) < 0.03
    assert pair['per_cluster'][0]['success'] and pair['per_cluster'][1]['success']
    assert pair['science_measurement'] is False
    assert 'mean tangential shear' in pair['deltasigma_note']
    assert 'member contamination' in pair['density_note']
    assert pair['covariance'] == 'shape'


def test_zero_shear_averages_near_zero(cosmo):
    record = _catalogue(cosmo, 300, 0.4, 2.5, 1.0e15, 0.3, 1.0, 12)
    record['g1'] = np.zeros_like(record['g1'])
    record['g2'] = np.zeros_like(record['g2'])
    other = dict(record, id='b')
    stacked = _stack(cosmo, [record, other], rmin=0.4, rmax=2.5, nbins=4, fit=False, per_cluster_mass=False)
    assert stacked['profile']
    assert max(abs(row['gt']) for row in stacked['profile']) < 1e-8


def test_a_stacked_image_is_not_a_shear_catalogue(cosmo):
    with pytest.raises(ValueError):
        wl.stack_shear_profiles(np.zeros((3, 12, 12)), cosmo=cosmo)
    with pytest.raises(ValueError):
        wl.stack_shear_profiles([np.zeros((12, 12)), np.zeros((12, 12))], cosmo=cosmo)
    with pytest.raises(ValueError):
        wl.stack_shear_profiles([dict(image=np.zeros((8, 8)), z_lens=0.3)], cosmo=cosmo)


def test_inverse_critical_density_is_not_the_value_at_the_mean_redshift(cosmo):
    z = np.array([0.55, 2.4])
    weight = np.array([1.0, 1.0])
    got = wl.photometric_geometry(
        cosmo, 0.25, z, weight, z_gap=0.05,
        cluster_window=(0.2, 0.3), background_window=(0.6, 1.4))
    assert got['mean_z'] == pytest.approx(1.475)
    assert got['label_window']['background'][0] == pytest.approx(0.6)
    ratio = got['mean_sigma_crit_inv'] / got['sigma_crit_inv_of_mean_z']
    assert abs(ratio - 1.0) > 1e-3
    assert 'spectroscopic redshift column' in got['note']
    with pytest.raises(ValueError):
        wl.stack_shear_profiles([dict(z_lens=0.3, z_phot=np.ones(3))], cosmo=cosmo)


def test_empty_spectroscopic_intersection_is_not_a_fit(cosmo):
    record = _catalogue(cosmo, 40, 0.4, 2.5, 1.0e15, 0.3, 1.0, 13)
    record['spec_quality'] = np.zeros(len(record['z']))
    stacked = _stack(cosmo, [record], rmin=0.4, rmax=2.5, nbins=4)
    assert stacked['success'] is False
    assert stacked['status'] == 'insufficient_sample'
    assert stacked['n_spec_intersection'] == 0
    assert stacked['science_measurement'] is False
    assert 'M200c_Msun' not in stacked
    unmatched = _catalogue(cosmo, 20, 0.4, 2.5, 1.0e15, 0.3, 1.0, 14)
    unmatched['crossmatch_ok'] = np.zeros(len(unmatched['z']), bool)
    lost = _stack(cosmo, [unmatched], rmin=0.4, rmax=2.5, nbins=4, fit=False, per_cluster_mass=False)
    assert lost['status'] == 'insufficient_sample'
    assert lost['success'] is False


def test_survey_weights_drop_the_shape_noise_formula(cosmo):
    record = _catalogue(cosmo, 300, 0.4, 2.5, 1.0e15, 0.3, 1.0, 15)
    record['weight'] = np.full(len(record['z']), 0.4)
    other = dict(record, id='b')
    omitted = _stack(cosmo, [record, other], rmin=0.4, rmax=2.5, nbins=4, fit=False, per_cluster_mass=False)
    assert omitted['covariance'] == 'omitted'
    assert omitted['profile'][0]['gt_err'] is None
    assert 'omitted' in omitted['covariance_note']
    kept = _stack(cosmo, [record, other], rmin=0.4, rmax=2.5, nbins=4, fit=False,
                  per_cluster_mass=False, inverse_variance_weights=True)
    assert kept['profile'][0]['gt_err'] > 0
    record['icl_weight'] = 3.0
    noted = _stack(cosmo, [record], rmin=0.4, rmax=2.5, nbins=4, fit=False, per_cluster_mass=False)
    assert 'same population' in noted['population_note']


def test_boost_multiplies_the_shear_and_its_error(cosmo):
    record = _catalogue(cosmo, 400, 0.4, 2.5, 1.0e15, 0.3, 1.0, 16)
    bare = _stack(cosmo, [record], rmin=0.4, rmax=2.5, nbins=4, sigma_e=0.2, fit=False, per_cluster_mass=False)
    rng = np.random.default_rng(16)
    inner = np.sqrt(0.4 ** 2 + rng.random(20) * (0.55 ** 2 - 0.4 ** 2))
    outer = np.sqrt(1.7 ** 2 + rng.random(900) * (2.5 ** 2 - 1.7 ** 2))
    random_R = np.concatenate([inner, outer])
    boosted_record = dict(record, random_R=random_R)
    boosted = _stack(cosmo, [boosted_record], rmin=0.4, rmax=2.5, nbins=4, sigma_e=0.2,
                     fit=False, per_cluster_mass=False)
    assert boosted['profile'][0]['boost'] > 1.05
    assert boosted['profile'][0]['gt_err'] == pytest.approx(
        bare['profile'][0]['gt_err'] * boosted['profile'][0]['boost'], rel=1e-6)
    wider = _stack(cosmo, [boosted_record], rmin=0.4, rmax=2.5, nbins=4, sigma_e=0.2,
                   boost_uncertainty=0.05, fit=False, per_cluster_mass=False)
    assert wider['profile'][0]['gt_err'] > boosted['profile'][0]['gt_err']


def test_two_source_planes_do_not_gain_a_large_scale_structure_term(cosmo):
    near = _catalogue(cosmo, 200, 0.4, 2.5, 1.0e15, 0.3, 0.8, 17, id='near')
    far = _catalogue(cosmo, 200, 0.4, 2.5, 1.0e15, 0.3, 1.6, 18, id='far')
    stacked = _stack(cosmo, [near, far], rmin=0.4, rmax=2.5, nbins=4, lss_cov=True,
                     fit=False, per_cluster_mass=False)
    assert stacked['covariance'] == 'shape'
    assert 'one source plane' in stacked['covariance_note']
    assert 'survey covariance' in stacked['covariance_note']


def test_r200c_stack_uses_the_caller_mass(cosmo):
    record = _catalogue(cosmo, 280, 0.5, 2.2, 1.0e15, 0.3, 1.0, 19)
    other = dict(record, id='b')
    radius = wl.r200c_mpc(record['mass'], record['z_lens'], cosmo)
    stacked = _stack(cosmo, [record, other], radius='r200c', rmin=0.25, rmax=1.2, nbins=4,
                     fit=False, per_cluster_mass=False)
    assert stacked['radius'] == 'r200c'
    assert stacked['profile']
    assert stacked['profile'][0]['R'] < 1.2
    assert 'caller mass' in stacked['radius_note']
    assert radius > 0


def test_science_measurement_waits_for_the_sensitivity_gates(cosmo):
    record = _catalogue(cosmo, 360, 0.4, 2.5, 1.0e15, 0.3, 1.0, 21)
    record['spec_quality'] = np.ones(len(record['z']))
    open_gate = _stack(
        cosmo, [record], rmin=0.4, rmax=2.5, nbins=4, per_cluster_mass=False, sensitivity=True,
        bcg_mass_offset_mpc=0.02,
        sensitivity_supplied={
            'non_sphericity': {'note': 'caller triaxial row for this methods test'},
            'projection': {'note': 'caller projection row for this methods test'}})
    assert open_gate['success']
    assert open_gate['n_spec_intersection'] > 0
    assert open_gate['science_measurement'] is True
    assert open_gate['sensitivity']['bcg_mass_offset_mpc'] == pytest.approx(0.02)
    assert open_gate['sensitivity']['non_sphericity']['source'] == 'caller'
    closed = _stack(cosmo, [record], rmin=0.4, rmax=2.5, nbins=4, per_cluster_mass=False, sensitivity=True)
    assert closed['science_measurement'] is False
    assert 'non_sphericity' in closed['science_gaps']


def test_cli_stack_writes_the_average_profile(tmp_path, cosmo):
    mass = 1.0e15
    record = _catalogue(cosmo, 280, 0.4, 2.5, mass, 0.3, 1.0, 22)
    other = dict(record, id='b')

    def _dump(item):
        out = {}
        for key, value in item.items():
            out[key] = value.tolist() if isinstance(value, np.ndarray) else value
        return out

    path = tmp_path / 'clusters.json'
    path.write_text(json.dumps({'clusters': [_dump(record), _dump(other)]}))
    work = tmp_path / 'stack'
    proc = _run([
        '--task', 'stack', '--clusters', str(path), '--work', str(work),
        '--radius', 'mpc', '--rmin', '0.4', '--rmax', '2.5', '--nbins', '4',
        '--concentration', 'fixed', '--c-fixed', '4', '--no-two-halo', '--power', 'linear',
        '--no-lss-cov', '--sigma-e', '0.001', '--samples', '0', '--min-per-bin', '8'])
    assert proc.returncode == 0, proc.stderr
    doc = json.loads((work / 'wl_stack.json').read_text())
    assert doc['success']
    assert doc['status'] == 'stacked'
    assert abs(math.log10(doc['M200c_Msun'] / mass)) < 0.03
    assert 'Gaussian PSF' in doc['metacal_scope']
    header = (work / 'wl_stack_profile.tsv').read_text().splitlines()[0].split('\t')
    assert 'GT' in header and 'DELTASIGMA' in header


def _tangential_stamps(cosmo, phis, radii, gt, z_lens=0.3, z_source=1.0):
    """Centred stamps with one tangential shear. Positions sit on those proper radii."""
    distance = float(cosmo.angular_diameter(z_lens))
    stamps, x, y, g1, g2 = [], [], [], [], []
    for phi, radius in zip(phis, radii):
        comp1, comp2 = wl.components_from_tangential(float(phi), float(gt), 0.0)
        comp1 = float(np.asarray(comp1))
        comp2 = float(np.asarray(comp2))
        image, _psf = wl.draw_galaxy(49, comp1, comp2, psf_sigma=1.2, noise=0.0, re=4.0, seed=1)
        theta = (float(radius) / distance) / wl.ARCSEC
        stamps.append(image)
        x.append(theta * math.cos(float(phi)))
        y.append(theta * math.sin(float(phi)))
        g1.append(comp1)
        g2.append(comp2)
    n = len(phis)
    return stamps, x, y, g1, g2, np.full(n, z_source), np.ones(n), z_lens


def test_tangential_stamp_stack_keeps_the_sign_and_the_sky_average_cancels(cosmo):
    angles = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    phis = np.concatenate([angles, angles])
    radii = np.concatenate([np.full(8, 0.7), np.full(8, 1.2)])
    stamps, x, y, g1, g2, z, w, z_lens = _tangential_stamps(cosmo, phis, radii, 0.15)
    doc = wl.stack_tangential_stamps(
        stamps, x, y, g1, g2, z, w, z_lens, cosmo, rmin=0.5, rmax=1.5, nbins=2, min_per_bin=4)
    assert doc['success']
    assert doc['status'] == 'diagnostic'
    assert doc['mass'] is None
    assert 'null image' in doc['note']
    assert 'DES Y3' in doc['note'] and 'lensfit' in doc['note'] and 'S19A' in doc['note']
    assert 'Gaussian PSF' in doc['metacal_scope']
    assert 'Lewis 2009' in doc['reference']
    assert len(doc['profile']) == 2
    assert len(doc['images']) == 2
    for row in doc['profile']:
        assert row['n'] == 8
        assert row['gt'] == pytest.approx(0.15)
        assert row['e_t'] > 0.05
        assert row['e_t_stack'] > 0.05
        assert row['e_t'] * row['gt'] > 0
        assert row['n_same_sign'] == row['n']
        assert row['e_abs_sky'] < 1e-6
        assert row['e_abs_sky'] < 0.02 * row['e_t_stack']
        assert abs(row['e_x_stack']) < 0.02 * row['e_t_stack']
        assert row['position_angle_span'] > 0.75 * 2.0 * math.pi
    turned = [row for row in doc['sources'] if row['e_t'] is not None]
    assert len(turned) == 16
    assert all(row['e_t'] > 0 and row['gt'] > 0 for row in turned)

    quiet, x0, y0, g10, g20, z0, w0, z0_lens = _tangential_stamps(cosmo, angles, np.full(8, 1.0), 0.0)
    zero = wl.stack_tangential_stamps(
        quiet, x0, y0, g10, g20, z0, w0, z0_lens, cosmo, rmin=0.5, rmax=1.5, nbins=1)
    assert abs(zero['profile'][0]['e_t_stack']) < 1e-8
    assert zero['profile'][0]['e_abs_sky'] < 1e-8

    flipped, xf, yf, g1f, g2f, zf, wf, zf_lens = _tangential_stamps(cosmo, angles, np.full(8, 1.0), -0.15)
    neg = wl.stack_tangential_stamps(
        flipped, xf, yf, g1f, g2f, zf, wf, zf_lens, cosmo, rmin=0.5, rmax=1.5, nbins=1)
    row = neg['profile'][0]
    assert row['gt'] == pytest.approx(-0.15)
    assert row['e_t_stack'] < -0.05
    assert row['n_same_sign'] == row['n']
    assert row['e_abs_sky'] < 1e-6


def test_stamps_at_one_angle_keep_a_sky_signal(cosmo):
    phis = np.zeros(8)
    stamps, x, y, g1, g2, z, w, z_lens = _tangential_stamps(cosmo, phis, np.full(8, 1.0), 0.15)
    doc = wl.stack_tangential_stamps(
        stamps, x, y, g1, g2, z, w, z_lens, cosmo, rmin=0.5, rmax=1.5, nbins=1)
    row = doc['profile'][0]
    assert row['position_angle_span'] == pytest.approx(0.0)
    assert row['e_t_stack'] > 0.05
    assert row['e_abs_sky'] == pytest.approx(row['e_t_stack'], rel=1e-5, abs=1e-8)
    assert row['n_same_sign'] == row['n']


def test_a_cluster_image_stack_is_not_a_postage_stamp_diagnostic(cosmo):
    with pytest.raises(ValueError, match='a stack of cluster images is not a tangential postage-stamp diagnostic'):
        wl.stack_tangential_stamps(np.zeros((2, 9, 9)), [1.0], [1.0], [0.0], [0.0], [1.0], [1.0], 0.3, cosmo)
    with pytest.raises(ValueError, match='a cluster image is not a list of postage stamps'):
        wl.stack_tangential_stamps(np.zeros((9, 9)), [1.0], [1.0], [0.0], [0.0], [1.0], [1.0], 0.3, cosmo)
    with pytest.raises(ValueError, match='postage stamps in one stack have the same shape'):
        wl.stack_tangential_stamps(
            [np.zeros((9, 9)), np.zeros((11, 11))], [1.0, 1.0], [1.0, 1.0],
            [0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [1.0, 1.0], 0.3, cosmo)

    angles = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    stamps, x, y, g1, g2, z, w, z_lens = _tangential_stamps(cosmo, angles, np.full(8, 1.0), 0.15)
    # The extra rows are a foreground galaxy and a satellite. select_sources drops both.
    blank = np.zeros((49, 49))
    doc = wl.stack_tangential_stamps(
        stamps + [blank, blank],
        list(x) + [x[0], x[1]], list(y) + [y[0], y[1]],
        list(g1) + [0.4, -0.4], list(g2) + [0.0, 0.2],
        list(z) + [0.05, 1.0], list(w) + [1.0, 1.0], z_lens, cosmo,
        rmin=0.5, rmax=1.5, nbins=1,
        member=['background'] * 8 + ['background', 'satellite'])
    assert doc['n_sources'] == 8
    assert doc['profile'][0]['n'] == 8
    assert doc['profile'][0]['gt'] == pytest.approx(0.15)

    empty = wl.stack_tangential_stamps(
        stamps, x, y, g1, g2, np.full(8, z_lens), w, z_lens, cosmo, rmin=0.5, rmax=1.5, nbins=1)
    assert empty['success'] is False
    assert empty['status'] == 'insufficient_sample'
    assert empty['profile'] == []
    assert empty['mass'] is None
