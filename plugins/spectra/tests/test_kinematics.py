import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from ogfkit import kinematics as K, spectrasynth as SS

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, '..', 'spectra.py')


def test_arctan_curve_and_disc_conventions():
    assert abs(K.arctan_rc(1e6, 200.0, 3.0) - 200.0) < 1e-3
    assert abs(K.arctan_rc(3.0, 200.0, 3.0) - 100.0) < 1e-9                # arctan(1) = pi/4 -> v = vc/2
    # PA = direction of the receding side; at 5 px from the centre along it the velocity is vsys + vc sin(i) * (2/pi) atan(5/rt)
    v = K.disc_velocity(20 + 5 * math.cos(math.radians(130)), 20 + 5 * math.sin(math.radians(130)), 20, 20, 10.0, 200.0, 4.0, 130.0, 60.0)
    assert abs(v - (10.0 + K.arctan_rc(5.0, 200.0, 4.0) * math.sin(math.radians(60)))) < 1e-9
    v_opp = K.disc_velocity(20 - 5 * math.cos(math.radians(130)), 20 - 5 * math.sin(math.radians(130)), 20, 20, 10.0, 200.0, 4.0, 130.0, 60.0)
    assert v_opp < 10.0 < v
    # along the minor axis the velocity is vsys
    assert abs(K.disc_velocity(20 + 5 * math.cos(math.radians(40)), 20 + 5 * math.sin(math.radians(40)), 20, 20, 10.0, 200.0, 4.0, 130.0, 60.0) - 10.0) < 1e-9


def test_adaptive_binning_reaches_target():
    rng = np.random.default_rng(0)
    sig = np.exp(-0.5 * ((np.arange(60) - 30) / 6.0) ** 2) * 6 + rng.normal(0, 0.3, 60)
    bins = K.adaptive_bins_1d(sig, np.ones(60), 8.0)
    assert bins and all(sig[b].sum() / math.sqrt(len(b)) >= 8.0 for b in bins)
    ids, b2 = K.accretion_bins_2d(np.outer(sig, sig) / 6, np.ones((60, 60)), 8.0)
    assert b2 and all(sum(np.outer(sig, sig)[j, i] / 6 for j, i in b) / math.sqrt(len(b)) >= 8.0 for b in b2)
    assert (ids >= 0).sum() == sum(len(b) for b in b2)


def test_slit_line_pull_and_rotation_curve():
    pulls, dsig = [], []
    for seed in range(8):
        w, d, t = SS.disc_slit(seed=seed, noise=0.15, vc=220.0, rt=3.5, inc=60.0, pa=40.0, sigma0=45.0)
        f = K.fit_slit_line(w, d, None, t['lam_sys'], window_kms=600, bin_snr=0, snr_row=6.0, inst_fwhm_A=2.0)
        assert len(f['y']) > 15
        rows = np.round(f['y']).astype(int)
        pulls += list((f['v'] - t['v_slit'][rows]) / f['v_err'])
        dsig += list(f['sigma'] / t['sigma0'])
    pulls = np.array(pulls)
    assert abs(np.mean(pulls)) < 0.25 and 0.7 < np.std(pulls) < 1.4
    assert abs(np.median(dsig) - 1.0) < 0.06
    w, d, t = SS.disc_slit(seed=11, noise=0.05, vc=220.0, rt=3.5, inc=60.0, pa=40.0, sigma0=45.0)
    f = K.fit_slit_line(w, d, None, t['lam_sys'], window_kms=600, bin_snr=6.0, inst_fwhm_A=2.0)
    rc = K.fit_rotation_curve(f['y'], f['v'], f['v_err'], inc_deg=60.0)
    assert rc['ok'] and abs(rc['vflat_obs'] - 220 * math.sin(math.radians(60))) < 4.0
    assert abs(rc['rt'] - 3.5) < 0.25 and abs(rc['y0'] - t['y0']) < 0.15 and abs(rc['vc'] - 220.0) < 5.0 and rc['chi2r'] < 2.5


def test_instrument_width_removal():
    w, d, t = SS.disc_cube(seed=2, noise=0.05, sigma0=30.0, vc=150.0, inc=50.0, inst_fwhm_A=3.0)
    ok = K.fit_cube_line(d, w, t['lam_sys'], window_kms=600, inst_fwhm_A=3.0)
    raw = K.fit_cube_line(d, w, t['lam_sys'], window_kms=600, inst_fwhm_A=0.0)
    s_ok = K.dispersion_summary(ok['sigma'], ok['sigma_err'], ok['flux'])['mean']
    s_raw = K.dispersion_summary(raw['sigma'], raw['sigma_err'], raw['flux'])['mean']
    assert abs(s_ok / 30.0 - 1) < 0.08 and s_raw > 1.8 * s_ok            # 3 A FWHM = 1.27 A sigma = 58 km/s -> sqrt(30^2 + 58^2) = 65


def test_cube_kinematics_recovery_high_sn():
    for seed, (pa, inc, vc, rt) in enumerate([(130.0, 55.0, 200.0, 4.0), (20.0, 40.0, 260.0, 3.0), (300.0, 68.0, 150.0, 5.0)]):
        w, c, t = SS.disc_cube(seed=seed + 5, noise=0.05, pa=pa, inc=inc, vc=vc, rt=rt, x0=19.4, y0=20.7, vsys=15.0, sigma0=50.0)
        m = K.fit_cube_line(c, w, t['lam_sys'], window_kms=600, snr_pix=3, inst_fwhm_A=2.0)
        vf = K.fit_velocity_field(m['vel'], m['vel_err'])
        assert vf['ok']
        assert abs(K.angle_diff(vf['pa'], pa)) < 1.0 and abs(vf['inc'] - inc) < 3.0
        assert abs(vf['vc'] - vc) / vc < 0.03 and abs(vf['rt'] - rt) / rt < 0.05 and abs(vf['vsys'] - 15.0) < 3.0
        assert math.hypot(vf['x0'] - 19.4, vf['y0'] - 20.7) < 0.1
        sd = K.dispersion_summary(m['sigma'], m['sigma_err'], m['flux'])
        assert abs(sd['mean'] / 50.0 - 1) < 0.08
        vf2 = K.fit_velocity_field(m['vel'], m['vel_err'], inc_fixed=inc)
        assert abs(vf2['vc'] - vc) / vc < 0.02 and vf2['errors']['vc'] < vf['errors']['vc'] * 1.05


def test_beam_smearing_biases_are_in_the_expected_direction():
    p = dict(pa=130.0, inc=55.0, vc=200.0, rt=3.0, sigma0=40.0, noise=0.05, seed=9)
    w, c0, t = SS.disc_cube(seeing_fwhm_pix=0.0, **p)
    w, c1, t1 = SS.disc_cube(seeing_fwhm_pix=3.0, **p)
    res = []
    for c in (c0, c1):
        m = K.fit_cube_line(c, w, t['lam_sys'], window_kms=600, inst_fwhm_A=2.0)
        res.append((K.fit_velocity_field(m['vel'], m['vel_err']), K.dispersion_summary(m['sigma'], m['sigma_err'], m['flux'])['mean']))
    assert res[1][1] > res[0][1] + 2.5                                   # the dispersion is inflated by the velocity gradient within the beam
    assert res[1][0]['vc'] < res[0][0]['vc'] or res[1][0]['rt'] > res[0][0]['rt']   # rotation curve looks flatter / slower rising


def test_binning_recovers_low_sn_cube():
    w, c, t = SS.disc_cube(seed=21, noise=0.6, pa=75.0, inc=50.0, vc=220.0, rt=3.0, sigma0=45.0)
    un = K.fit_cube_line(c, w, t['lam_sys'], window_kms=600, snr_pix=3, inst_fwhm_A=2.0)
    bn = K.fit_cube_line(c, w, t['lam_sys'], window_kms=600, bin_snr=8.0, inst_fwhm_A=2.0)
    assert bn['n_spaxels'] > un['n_spaxels']
    for m in (un, bn):
        ok = np.isfinite(m['vel'])
        m['rms'] = np.sqrt(np.mean((m['vel'] - t['velocity'])[ok] ** 2))
    assert bn['rms'] < un['rms']
    vf = K.fit_velocity_field(bn['vel'], bn['vel_err'])
    assert vf['ok'] and abs(K.angle_diff(vf['pa'], 75.0)) < 4.0 and abs(vf['vsini'] - 220 * math.sin(math.radians(50))) < 20


def test_velocity_field_refuses_too_little_data():
    v = np.full((20, 20), np.nan)
    v[5:8, 5:8] = 10.0
    assert not K.fit_velocity_field(v, np.full_like(v, 5.0))['ok']
    assert not K.fit_rotation_curve([1, 2, 3], [0, 1, 2], [1, 1, 1])['ok']


def _strict_json(path):
    def bad(c):
        raise ValueError('non-finite constant %s' % c)
    return json.load(open(path), parse_constant=bad)


def test_cli_kinematics_for_2d_and_cube(tmp_path):
    from astropy.io import fits
    d = str(tmp_path)
    w, s2, t2 = SS.disc_slit(seed=3, noise=0.05, vc=220.0, rt=3.5, inc=60.0, pa=40.0, sigma0=45.0, z=0.12)
    h = fits.PrimaryHDU(s2.astype(np.float32)); h.header['CRVAL1'] = w[0]; h.header['CDELT1'] = w[1] - w[0]; h.header['CRPIX1'] = 1
    h.writeto(os.path.join(d, 'spec_1.fits'))
    wc, cc, tc = SS.disc_cube(seed=4, noise=0.05, pa=130.0, inc=55.0, vc=200.0, rt=4.0, sigma0=40.0, z=0.12)
    h = fits.PrimaryHDU(cc.astype(np.float32)); h.header['CRVAL3'] = wc[0]; h.header['CDELT3'] = wc[1] - wc[0]; h.header['CRPIX3'] = 1
    h.writeto(os.path.join(d, 'spec_2.fits'))
    with open(os.path.join(d, 'c.tsv'), 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tZ\n1\t5\t5\t0.12\n2\t21\t21\t0.12\n3\t9\t9\t0.12\n')
    out = str(tmp_path / 'w')
    base = [sys.executable, CLI, '--catalog', os.path.join(d, 'c.tsv'), '--work', out, '--spec-dir', d, '--cube-xy', '1', '--z-column', 'Z', '--kin-inst-fwhm', '2.0', '--kin-window-kms', '600',
            '--extract-halfwidth', '3']
    p = subprocess.run(base + ['--task', 'kin', '--kin-inc', '60.0'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-600:]
    lines = [l.split('\t') for l in p.stdout.splitlines()]
    cols = lines[0]
    assert cols[0] == 'NUMBER' and 'SP_KIN_VSINI' in cols, (p.stdout[:300], p.stderr[-500:])
    assert cols[0] == 'NUMBER' and 'SP_KIN_VSINI' in cols and 'SP_KIN_PA' in cols
    rows = {l[0]: dict(zip(cols, l)) for l in lines[1:]}
    assert rows['3']['SP_KIN_VSINI'] == ''                                # object 3 has no spectrum
    # object 1: slit, inclination given -> deprojected v_c
    assert abs(float(rows['1']['SP_KIN_VC']) - 220.0) < 6.0 and abs(float(rows['1']['SP_KIN_RT']) - 3.5) < 0.3
    assert abs(float(rows['1']['SP_KIN_SIGMA']) - 45.0) < 6.0 and rows['1']['SP_KIN_PA'] == ''
    p2 = subprocess.run(base + ['--task', 'kin'], capture_output=True, text=True)
    lines = [l.split('\t') for l in p2.stdout.splitlines()]
    rows = {l[0]: dict(zip(lines[0], l)) for l in lines[1:]}
    assert abs(K.angle_diff(float(rows['2']['SP_KIN_PA']), 130.0)) < 1.5 and abs(float(rows['2']['SP_KIN_INC']) - 55.0) < 4.0 and abs(float(rows['2']['SP_KIN_VC']) - 200.0) / 200.0 < 0.04
    for f in ('kin_1.png', 'kin_2.png', 'kin_2_maps.fits', 'kin_2_vel.fits', 'kin_2_sigma.fits', 'kin_2_flux.fits', 'kin_1.tsv', 'spectra_kin.json'):
        assert os.path.getsize(os.path.join(out, f)) > 100, f
    res = _strict_json(os.path.join(out, 'spectra_kin.json'))
    assert res['2']['kind'] == 'cube' and res['1']['kind'] == '2d' and 'velocity_field' in res['2']
    hl = fits.open(os.path.join(out, 'kin_2_maps.fits'))
    assert [h.name for h in hl][1:7] == ['FLUX', 'VEL', 'VELERR', 'SIGMA', 'SIGERR', 'SNR'] and 'MODEL' in [h.name for h in hl] and hl['VEL'].data.shape == (41, 41)
    p3 = subprocess.run([a for a in base if a != '--z-column' and a != 'Z'] + ['--task', 'kin'], capture_output=True, text=True)
    assert p3.returncode == 0                                              # blind redshift fallback must not crash even when it finds nothing
