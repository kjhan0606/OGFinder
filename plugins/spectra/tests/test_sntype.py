"""Supernova types on synthetic P-Cygni spectra.  Not a survey validation and not a template match."""
import json
import os
import subprocess
import sys

import numpy as np

from ogfkit import sntype as sn
from ogfkit import specscience as sci

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'spectra.py')
TOL = 1500.0


def run_cli(args):
    p = subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    lines = [l for l in p.stdout.splitlines() if l]
    cols = lines[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t'))) for l in lines[1:]], p.stderr


def _pcygni(depths, z=0.05, v=10000.0, sigma_v=2500.0, noise=0.015, seed=1, frame='vacuum', emission=None):
    """Flat continuum minus a Gaussian absorption at rest*(1+z)*(1-v/c)."""
    w = np.arange(4500.0, 9200.0, 2.0)
    f = np.ones(w.shape, float)
    for name, depth in depths.items():
        lam = sn.FEAT[name][0]
        cen = lam * (1.0 + z) * (1.0 - v / sn.C_KMS)
        sig = max(cen * sigma_v / sn.C_KMS, 1.0)
        f = f - depth * np.exp(-0.5 * ((w - cen) / sig) ** 2)
    if emission:
        for name, sig_v, amp in emission:
            lam = sn.FEAT[name][0]
            cen = lam * (1.0 + z)
            sig = max(cen * sig_v / sn.C_KMS, 1.0)
            f = f + amp * np.exp(-0.5 * ((w - cen) / sig) ** 2)
    rng = np.random.default_rng(seed)
    f = f + rng.normal(0.0, noise, size=w.shape)
    err = np.full(w.shape, noise)
    if frame == 'air':
        w = np.asarray(sci.vacuum_to_air(w), float)
    return w, f, err


def test_si_trough_is_ia_not_ii():
    # The Si II 6355 minimum at 10000 km/s lands near 6150 A, where H-alpha at
    # ~19000 km/s would also fall.  H-beta is absent, so this is not Type II.
    w, f, e = _pcygni({'SiII6355': 0.55, 'SII5454': 0.25, 'SII5640': 0.25, 'OI7774': 0.15}, v=10000.0, seed=1)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ia', an
    assert an['a_si'] > 0.35
    assert an['v_line'] == 'SiII6355'
    assert abs(an['v'] - 10000.0) < TOL
    assert an['hydrogen'] is False
    assert an['features'].get('Hb', {}).get('detected') is not True


def test_type_ii():
    w, f, e = _pcygni({'Ha': 0.45, 'Hb': 0.30}, v=8000.0, seed=2)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'II', an
    assert an['v_line'] == 'Ha'
    assert abs(an['v'] - 8000.0) < TOL
    assert an['quality'] == 1


def test_type_ib():
    w, f, e = _pcygni({'HeI5876': 0.40, 'HeI6678': 0.30, 'HeI7065': 0.28}, v=9000.0, seed=3)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ib', an
    assert an['v_line'] in sn.HE_NAMES
    assert abs(an['v'] - 9000.0) < TOL
    assert an['helium_n'] >= 2
    assert an['hydrogen'] is False


def test_type_ic():
    w, f, e = _pcygni({'OI7774': 0.45, 'SiII6355': 0.20}, v=10000.0, seed=4)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ic', an
    assert an['a_si'] is not None and an['a_si'] < 0.35
    assert an['a_oi'] is not None and an['a_si'] / an['a_oi'] < 1.0
    assert abs(an['v'] - 10000.0) < TOL


def test_type_iib():
    w, f, e = _pcygni({'Ha': 0.30, 'Hb': 0.20, 'HeI5876': 0.35, 'HeI6678': 0.35, 'HeI7065': 0.35}, v=8000.0, seed=5)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'IIb', an
    assert an['v_line'] == 'Ha'
    assert abs(an['v'] - 8000.0) < TOL
    assert an['hydrogen'] is True and an['helium_n'] >= 2


def test_type_iin():
    w, f, e = _pcygni({}, emission=[('Ha', 400.0, 2.0)], seed=6)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'IIn', an
    assert an['v_line'] == 'Ha'
    assert abs(an['v']) < 500.0
    assert an['narrow']['Ha']['fwhm'] < 2000.0


def test_type_ibn():
    w, f, e = _pcygni({}, emission=[('HeI5876', 400.0, 1.5), ('HeI6678', 400.0, 1.2), ('HeI7065', 400.0, 1.1)], seed=16)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ibn', an
    assert an['v_line'] in sn.HE_NAMES
    assert abs(an['v']) < 500.0


def test_type_icbl():
    w, f, e = _pcygni({'OI7774': 0.25, 'SiII6355': 0.25}, v=32000.0, sigma_v=8000.0, seed=7)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ic-BL', an
    assert an['v'] >= 0.1 * sn.C_KMS
    assert abs(an['v'] - 32000.0) < TOL


def test_fast_strong_si_stays_ia():
    w, f, e = _pcygni({'SiII6355': 0.55, 'SII5454': 0.25, 'SII5640': 0.25}, v=32000.0, seed=8)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'Ia', an
    assert an['a_si'] > 0.35
    assert abs(an['v'] - 32000.0) < TOL


def test_weak_si_is_uncertain():
    w, f, e = _pcygni({'SiII6355': 0.15}, v=10000.0, seed=9)
    an = sn.sn_analysis(w, f, e, z_host=0.05)
    assert an['type'] == 'uncertain', an
    assert an['quality'] == 0
    assert an['a_si'] is not None and an['a_si'] < 0.35


def test_spike_is_not_a_line():
    w = np.arange(4500.0, 9200.0, 2.0)
    f = np.ones(w.shape, float)
    err = np.full(w.shape, 0.015)
    lam = sn.FEAT['SiII6355'][0]
    cen = lam * 1.05 * (1.0 - 10000.0 / sn.C_KMS)
    f[int(np.argmin(np.abs(w - cen)))] = 0.2
    an = sn.sn_analysis(w, f, err, z_host=0.05)
    assert an['type'] == 'uncertain', an
    assert an['features'].get('SiII6355', {}).get('detected') is not True


def test_missing_redshift_is_unknown():
    w, f, e = _pcygni({'Ha': 0.45, 'Hb': 0.30}, v=8000.0, seed=10)
    for z in (None, 0.0, float('nan'), ''):
        an = sn.sn_analysis(w, f, e, z_host=z)
        assert an['type'] == 'unknown' and an['quality'] == 0, an
        assert 'host redshift' in an['note']
        assert an['z'] is None


def test_air_frame_ii():
    w, f, e = _pcygni({'Ha': 0.45, 'Hb': 0.30}, v=8000.0, seed=11, frame='air')
    an = sn.sn_analysis(w, f, e, z_host=0.05, frame='air')
    assert an['type'] == 'II', an
    assert abs(an['v'] - 8000.0) < TOL


def test_too_few_pixels():
    w = np.linspace(5000.0, 7000.0, 10)
    an = sn.sn_analysis(w, np.ones(10), np.full(10, 0.01), z_host=0.05)
    assert an['type'] == 'unknown' and an['quality'] == 0
    assert an['z'] == 0.05


def test_cli_sntype_columns(tmp_path):
    ia = _pcygni({'SiII6355': 0.55, 'SII5454': 0.25, 'SII5640': 0.25, 'OI7774': 0.15}, v=10000.0, seed=12)
    ii = _pcygni({'Ha': 0.45, 'Hb': 0.30}, v=8000.0, seed=13)
    np.savetxt(tmp_path / 'spec_1.txt', np.c_[ia[0], ia[1], ia[2]])
    np.savetxt(tmp_path / 'spec_2.txt', np.c_[ii[0], ii[1], ii[2]])
    (tmp_path / 'c.tsv').write_text('NUMBER\tZ\n1\t0.05\n2\t0.05\n')
    cols, rows, err = run_cli(['--task', 'sntype', '--catalog', str(tmp_path / 'c.tsv'), '--work', str(tmp_path / 'w'),
                               '--spec-dir', str(tmp_path), '--spec-pattern', 'spec_{NUMBER}.txt', '--z-column', 'Z'])
    assert cols == ['NUMBER', 'SN_TYPE', 'SN_Z', 'SN_ZQ', 'SN_V', 'SN_VLINE', 'SN_A_SI', 'SN_A_HA', 'SN_A_HE', 'SN_A_OI']
    assert rows[0]['SN_TYPE'] == 'Ia'
    assert rows[1]['SN_TYPE'] == 'II'
    assert abs(float(rows[0]['SN_V']) - 10000.0) < TOL
    assert abs(float(rows[1]['SN_V']) - 8000.0) < TOL
    assert float(rows[0]['SN_A_SI']) > 0.35
    res = json.load(open(tmp_path / 'w' / 'spectra_sntype.json'))
    assert '401' in res['1']['calibration']
    assert '0.35' in res['1']['calibration']
    assert 'host redshift' in res['1']['calibration']
    assert 'spectra sntype' in err
