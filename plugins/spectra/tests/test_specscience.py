"""Science measurements on synthetic 1D spectra: calibrations, BPT type, absorption redshift, dispersion, H-alpha SFR."""
import json
import math
import os
import subprocess
import sys

import numpy as np

from ogfkit import spectra as sp
from ogfkit import specscience as sci
from ogfkit import spectrasynth as ss

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'spectra.py')
FWHM = sp.FWHM


def run_cli(args):
    p = subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    lines = [l for l in p.stdout.splitlines() if l]
    cols = lines[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t'))) for l in lines[1:]], p.stderr


def test_air_vacuum_and_ccm_and_calzetti():
    vac = sci.air_to_vacuum(6562.801)
    assert abs(vac - 6564.61) < 0.05
    assert abs(sci.vacuum_to_air(vac) - 6562.801) < 1e-3
    # V band: A(5500) / E(B-V) is Rv for the Cardelli optical curve.
    assert abs(sci.ccm_extinction(5500.0, 1.0) - 3.1) < 0.02
    assert sci.ccm_extinction(5500.0, 0.0) == 0.0
    from sed_fit.grid.generate import calzetti_extinction
    wave = np.array([4000.0, 4862.68, 5500.0, 6564.61, 9000.0])
    assert np.allclose(sci.calzetti_a(wave, 1.0 / sci.RV_CALZETTI), calzetti_extinction(wave, 1.0))


def test_prepare_sorts_drops_and_dereddens():
    w = np.array([5.0, 1.0, np.nan, 3.0, 4.0])
    f = np.array([50.0, 10.0, 10.0, 30.0, 40.0])
    e = np.array([1.0, 1.0, 1.0, 1.0, np.nan])
    prep = sci.prepare_spectrum(w, f, e, frame='vacuum', ebv_mw=0.0)
    assert prep['n_dropped'] == 2
    assert np.allclose(prep['wave'], [1.0, 3.0, 5.0])
    assert np.allclose(prep['flux'], [10.0, 30.0, 50.0])
    ebv = 0.1
    raw = np.array([5000.0, 6000.0, 7000.0])
    flux = np.array([2.0, 3.0, 4.0])
    err = np.array([0.2, 0.2, 0.2])
    A = sci.ccm_extinction(raw, ebv)
    scale = 10.0 ** (-0.4 * A)
    back = sci.prepare_spectrum(raw, flux * scale, err * scale, ebv_mw=ebv)
    assert np.allclose(back['flux'], flux)
    assert np.allclose(back['err'], err)
    same = sci.prepare_spectrum(raw, flux, err, ebv_mw=0.0)
    assert np.allclose(same['flux'], flux)


def test_detect_absorption_sign_does_not_change_emission():
    w = np.arange(6000.0, 7000.0, 1.0)
    f = np.ones_like(w)
    err = np.full_like(w, 0.02)
    f = f + 0.8 * np.exp(-0.5 * ((w - 6564.0) / 2.0) ** 2)
    f = f - 0.8 * np.exp(-0.5 * ((w - 6300.0) / 2.0) ** 2)
    em = sp.detect_lines(w, f, err, snr_min=5.0)
    ab = sp.detect_lines(w, f, err, snr_min=5.0, sign=-1.0)
    assert any(abs(d['wave'] - 6564.0) < 3 for d in em)
    assert not any(abs(d['wave'] - 6300.0) < 3 for d in em)
    assert any(abs(d['wave'] - 6300.0) < 3 for d in ab)
    assert all(d['snr'] > 0 for d in ab)


def test_bpt_and_sfr_formulas():
    x = math.log10(0.1)
    y_sf = math.log10(10.0 / (100.0 / 2.86))
    assert y_sf < sci.kauffmann2003(x) < sci.kewley2001(x)
    assert sci.classify_bpt(x, y_sf)['class_'] == 'star-forming'
    x_agn = math.log10(0.8)
    y_agn = math.log10(200.0 / 35.0)
    assert y_agn > sci.kewley2001(x_agn)
    assert sci.classify_bpt(x_agn, y_agn)['class_'] == 'AGN'
    assert sci.classify_bpt(x_agn, 1.0, log_sii_ha=0.0)['class_'] == 'Seyfert'
    assert sci.classify_bpt(x_agn, 1.0, log_sii_ha=0.2)['class_'] == 'LINER'
    assert sci.classify_bpt(x_agn, 1.0, log_oi_ha=0.0)['class_'] == 'LINER'
    assert sci.classify_bpt(x_agn, 1.5, log_oi_ha=0.0)['class_'] == 'Seyfert'
    below = sci.balmer_ebv(2.0, 1.0)
    assert below['ebv'] == 0.0 and below['used'] is False
    # A known extra reddening must come back positive.
    k_ha = float(np.asarray(sci.calzetti_k(6564.61)))
    k_hb = float(np.asarray(sci.calzetti_k(4862.68)))
    ebv = 0.2
    ratio = 2.86 * 10.0 ** (0.4 * ebv * (k_hb - k_ha))
    got = sci.balmer_ebv(ratio * 35.0, 35.0)
    assert got['used'] and abs(got['ebv'] - ebv) < 1e-6
    from astropy.cosmology import FlatLambdaCDM
    z, flux, ferr = 0.2, 1.0e-15, 1.0e-16
    dl = float(FlatLambdaCDM(H0=70.0, Om0=0.3).luminosity_distance(z).to('cm').value)
    expect = sci.SFR_HA * 4.0 * math.pi * dl * dl * flux
    got = sci.halpha_sfr(flux, z, flux_err=ferr)
    assert abs(got['sfr'] - expect) / expect < 1e-6
    assert abs(got['sfr_err'] - expect * ferr / flux) / got['sfr_err'] < 1e-6
    assert got['error_includes_balmer'] is False
    assert got['calibration'] == 'Murphy2011/KennicuttEvans2012'
    assert sci.halpha_sfr(flux, 0.0)['sfr'] is None


def _emission(z, lines, seed, fwhm_kms=200.0):
    return ss.spectrum_1d(z=z, seed=seed, wave=(5400.0, 8400.0), dlam=1.0, cont=1.0, slope=0.0,
                          noise=0.03, lines=lines, fwhm_kms=fwhm_kms, res_fwhm_A=2.0)


def test_star_forming_redshift_dispersion_and_sfr():
    lines = {'Ha': 100.0, 'Hb': 100.0 / 2.86, '[OIII]5008': 10.0, '[NII]6585': 10.0}
    w, f, e, _t = _emission(0.2, lines, seed=4)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.6, snr_min=4.0)
    assert an['type'] == 'star-forming', an['type']
    assert an['redshift']['source'] == 'emission'
    assert abs(an['redshift']['z'] - 0.2) < 0.002
    assert an['dispersion']['kind'] == 'gas' and an['dispersion']['corrected'] is True
    assert abs(an['dispersion']['sigma'] - 200.0 / FWHM) < 25.0
    assert abs(an['ha_flux'] - 100.0) / 100.0 < 0.08
    used = an['sfr']['ha_flux_dustcorr'] if an['sfr'].get('ha_flux_dustcorr') else an['ha_flux']
    expect = sci.halpha_sfr(used, an['redshift']['z'])
    assert an['sfr']['applies'] is True
    assert abs(an['sfr']['sfr'] - expect['sfr']) / expect['sfr'] < 1e-6
    # The same wavelengths written as air must come back to the same redshift.
    air = sci.vacuum_to_air(w)
    an_air = sci.galaxy_analysis(air, f, e, frame='air', inst_fwhm_a=2.0, zmax=0.6)
    assert abs(an_air['redshift']['z'] - an['redshift']['z']) < 5e-4


def test_agn_and_broad_line():
    lines = {'Ha': 100.0, 'Hb': 35.0, '[OIII]5008': 200.0, '[NII]6585': 80.0}
    w, f, e, _t = _emission(0.2, lines, seed=6)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.6)
    assert an['type'] == 'AGN', an['type']
    assert an['bpt']['class'] == 'AGN'
    assert an['sfr']['applies'] is False
    # Narrow Hb and [OIII] fix the redshift. Ha is broader than the detection kernel.
    w = np.arange(5400.0, 8400.0, 1.0)
    f = np.ones_like(w)
    z = 0.2
    for rest, flux, fwhm_kms in ((4862.68, 80.0, 250.0), (5008.24, 80.0, 250.0), (6564.61, 400.0, 3000.0)):
        lam = rest * (1.0 + z)
        sig = math.hypot(fwhm_kms / sci.C_KMS * lam, 2.0) / FWHM
        f = f + flux / (sig * math.sqrt(2.0 * math.pi)) * np.exp(-0.5 * ((w - lam) / sig) ** 2)
    e = np.full_like(w, 0.03)
    f = f + np.random.default_rng(8).normal(0.0, 0.03, len(w))
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.6, snr_min=5.0)
    assert an['type'] == 'broad-line AGN', an['type']
    assert an['broad']['flag'] is True and an['broad']['fwhm_kms'] > 1200.0


def _absorption(z, seed, wave=(3800.0, 9200.0)):
    rng = np.random.default_rng(seed)
    w = np.arange(wave[0], wave[1], 1.0)
    cont = 40.0
    f = np.full_like(w, cont)
    for _name, lam0, _wt in sci.ABS_LINES:
        lam = lam0 * (1.0 + z)
        if w[0] + 8 < lam < w[-1] - 8:
            sig = math.hypot(180.0 / sci.C_KMS * lam, 2.0) / FWHM
            f = f - cont * 0.55 * np.exp(-0.5 * ((w - lam) / sig) ** 2)
    err = np.full_like(w, 0.04)
    f = f + rng.normal(0.0, 0.04, len(w))
    return w, f, err


def test_one_false_emission_line_does_not_override_a_star():
    w, f, e = _absorption(0.0002, seed=3)
    lam = 4500.0
    sig = 2.0
    f = f + 30.0 / (sig * math.sqrt(2.0 * math.pi)) * np.exp(-0.5 * ((w - lam) / sig) ** 2)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.5, snr_min=5.0)
    assert an['type'] == 'star', (an['type'], an['redshift'], an['emission'])
    assert an['redshift']['source'] == 'absorption'
    assert abs(an['redshift']['z'] - 0.0002) < 0.001
    assert an['emission']['quality'] <= 1 or abs((an['emission']['z'] or 0.0) - an['redshift']['z']) > 0.01


def test_ultraviolet_broad_lines_set_redshift():
    z = 1.6
    lines = {'CIV': 80.0, 'CIII]': 40.0, 'MgII': 60.0}
    # Noise high enough that the 2-pixel kernel misses the broad lines. The 8-pixel
    # ultraviolet pass still fits CIV, C III] and Mg II.
    w, f, e, _t = ss.spectrum_1d(z=z, seed=9, wave=(3700.0, 9000.0), dlam=1.0, cont=1.0, slope=0.0,
                                 noise=0.05, lines=lines, fwhm_kms=4000.0, res_fwhm_A=2.0)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=2.0, snr_min=4.0)
    assert an['redshift']['source'] == 'uv', (an['redshift'], an['emission'], an['absorption'], an.get('uv'))
    assert abs(an['redshift']['z'] - z) < 0.005
    assert an['uv']['n_lines'] >= 2
    assert an['type'] == 'broad-line AGN', (an['type'], an['broad'])
    assert an['broad']['flag'] is True and an['broad']['fwhm_kms'] > 1200.0
    assert an['broad']['line'] in ('CIV', 'CIII]', 'MgII')


def test_absorption_galaxy_and_star():
    w, f, e = _absorption(0.15, seed=2)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.4, snr_min=5.0)
    assert an['redshift']['source'] == 'absorption', an['redshift']
    assert abs(an['redshift']['z'] - 0.15) < 0.002
    assert an['type'] == 'absorption galaxy', an['type']
    assert an['emission']['quality'] == 0
    w, f, e = _absorption(0.0002, seed=3)
    an = sci.galaxy_analysis(w, f, e, inst_fwhm_a=2.0, zmax=0.05, snr_min=5.0)
    assert an['type'] == 'star', (an['type'], an['redshift'])
    assert an['redshift']['source'] == 'absorption'
    assert an['redshift']['z'] < 0.002
    assert abs(an['redshift']['z'] - 0.0002) < 0.001


def test_cli_science_columns_and_json(tmp_path):
    lines = {'Ha': 100.0, 'Hb': 100.0 / 2.86, '[OIII]5008': 10.0, '[NII]6585': 10.0}
    w, f, e, _t = _emission(0.2, lines, seed=11)
    np.savetxt(tmp_path / 'spec_1.txt', np.c_[w, f, e])
    (tmp_path / 'c.tsv').write_text('NUMBER\tX_IMAGE\tY_IMAGE\n1\t1\t1\n2\t1\t1\n')
    cols, rows, err = run_cli(['--task', 'science', '--catalog', str(tmp_path / 'c.tsv'), '--work', str(tmp_path / 'w'),
                               '--spec-dir', str(tmp_path), '--spec-pattern', 'spec_{NUMBER}.txt', '--zmax', '0.6',
                               '--inst-fwhm', '2', '--h0', '70', '--omega-m', '0.3'])
    assert cols == ['NUMBER'] + sci_columns()
    assert rows[0]['SP_TYPE'] == 'star-forming'
    assert rows[0]['SP_SCI_ZSRC'] == 'emission'
    assert abs(float(rows[0]['SP_SCI_Z']) - 0.2) < 0.002
    assert rows[0]['SP_SIGMA_KIND'] == 'gas' and rows[0]['SP_SFR'] != ''
    assert rows[1]['SP_TYPE'] == '' and rows[1]['SP_SCI_Z'] == ''
    res = json.load(open(tmp_path / 'w' / 'spectra_science.json'))
    assert res['1']['calibrations']['sfr'].startswith('Murphy et al. 2011')
    assert res['1']['calibrations']['flux_unit'] == 'file'
    assert res['1']['sfr']['error_includes_balmer'] is False
    assert '2' not in res
    assert 'spectra science' in err


def sci_columns():
    return ['SP_SCI_Z', 'SP_SCI_ZERR', 'SP_SCI_ZQ', 'SP_SCI_ZSRC', 'SP_SIGMA', 'SP_SIGMA_ERR', 'SP_SIGMA_KIND',
            'SP_TYPE', 'SP_SFR', 'SP_SFR_EBV', 'SP_HA_FLUX']
