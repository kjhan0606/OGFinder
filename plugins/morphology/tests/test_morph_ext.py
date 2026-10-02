"""Extended non-parametric morphology (Petrosian radii, Kron, growth curves, smoothness, Petrosian-segmented Gini/M20) vs analytic Sersic references and noisy realisations."""
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import models as M, tsvio
from morphometry import extended as X

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
PY = sys.executable
SHAPE = (401, 401)
XC = YC = 200.0
FTOT = 1e5
CASES = [(1.0, 15, 1.0, 0), (1.0, 15, 0.6, 40), (4.0, 12, 1.0, 0), (2.0, 12, 0.7, 100)]


def gal(n, re_, q, pa, flux=FTOT, shape=SHAPE, xc=XC, yc=YC):
    return M.render_sersic(shape, xc, yc, M.sersic_Ie_from_flux(flux, re_, n, q), re_, n, q, pa)


def true_rp(n, re_, level):
    rg = np.geomspace(0.05, 20, 6000)
    eta = X.sersic_eta_ring(rg, n)
    return float(rg[np.argmax(eta < level)] * re_)


@pytest.mark.parametrize('n,re_,q,pa', CASES)
def test_petrosian_radii_flux_fractions_and_concentration(n, re_, q, pa):
    img = gal(n, re_, q, pa)
    r = X.measure_object(img, XC, YC, q, math.radians(pa), rms=1e-3, rmax=150)
    rp = {e: true_rp(n, re_, e) for e in (0.1, 0.2, 0.3)}
    errs = [r['rp'][e] / rp[e] - 1 for e in (0.1, 0.2, 0.3)]
    Fp = float(X.sersic_fraction(2 * rp[0.2] / re_, n))
    ana = [X.sersic_radius_fraction(f * Fp, n) * re_ for f in (0.2, 0.5, 0.8, 0.9)]
    meas = [r['r_frac'][f] for f in (0.2, 0.5, 0.8, 0.9)]
    cerr = r['conc'] - 5 * math.log10(ana[2] / ana[0])
    print('MORPHEXT n=%.0f q=%.1f: R_P(0.1/0.2/0.3) rel err %s | R20/50/80/90 rel err %s | C err %+.3f | F_P/F_tot %.4f (true %.4f)' % (
        n, q, ['%+.4f' % e for e in errs], ['%+.4f' % (m / a - 1) for m, a in zip(meas, ana)], cerr, r['flux_p'] / FTOT, Fp))
    assert max(abs(e) for e in errs) < 0.01
    assert max(abs(m / a - 1) for m, a in zip(meas, ana)) < 0.015 and abs(cerr) < 0.03
    assert abs(r['flux_p'] / FTOT - Fp) < 0.005 and r['flag'] == 0


@pytest.mark.parametrize('n,re_,q,pa', CASES)
def test_kron_flux_matches_analytic_aperture(n, re_, q, pa):
    img = gal(n, re_, q, pa)
    r = X.measure_object(img, XC, YC, q, math.radians(pa), rms=1e-3, rmax=150)
    frac = float(X.sersic_fraction(2.5 * r['kron_r'] / re_, n))
    print('MORPHEXT Kron n=%.0f q=%.1f: r_k %.2f px (re %d), flux in 2.5 r_k / total %.4f, analytic aperture fraction %.4f' % (n, q, r['kron_r'], re_, r['kron_flux'] / FTOT, frac))
    assert abs(r['kron_flux'] / FTOT - frac) < 0.01
    assert 0.6 * re_ < r['kron_r'] < 1.4 * re_


def test_noise_robustness_and_smoothness_noise_correction():
    n, re_, q, pa = 1.0, 15, 0.7, 30
    clean = gal(n, re_, q, pa, flux=FTOT)
    truth = true_rp(n, re_, 0.2)
    rows = {}
    for label, rms in (('S/N high (rms 0.5)', 0.5), ('S/N low (rms 8)', 8.0)):
        rp, S, S_raw = [], [], []
        for seed in range(12):
            img = clean + np.random.default_rng(seed).normal(0, rms, clean.shape)
            r = X.measure_object(img, XC, YC, q, math.radians(pa), rms=rms, rmax=150, seed=seed)
            rp.append(r['rp'][0.2]); S.append(r['smooth'])
            # the same without the noise correction
            h = int(1.5 * r['rp'][0.2]) + 3
            cut = img[int(YC) - h:int(YC) + h + 1, int(XC) - h:int(XC) + h + 1]
            yy, xx = np.indices(cut.shape)
            u = (xx - h) * math.cos(math.radians(pa)) + (yy - h) * math.sin(math.radians(pa)); v = -(xx - h) * math.sin(math.radians(pa)) + (yy - h) * math.cos(math.radians(pa))
            ap = np.sqrt(u ** 2 + (v / q) ** 2) <= 1.5 * r['rp'][0.2]
            S_raw.append(X.smoothness(cut, ap, max(0.25 * r['rp'][0.2], 3.0), None, centre=(h, h)))
        rows[label] = (np.mean(rp) / truth - 1, np.std(rp) / truth, np.mean(S), np.std(S), np.mean(S_raw))
        print('MORPHEXT noise %-20s R_P bias %+.3f scatter %.3f | S corrected %.3f +- %.3f (uncorrected %.3f)' % ((label,) + rows[label]))
    assert abs(rows['S/N high (rms 0.5)'][0]) < 0.02 and abs(rows['S/N low (rms 8)'][0]) < 0.05 and rows['S/N low (rms 8)'][1] < 0.06
    hi, lo = rows['S/N high (rms 0.5)'], rows['S/N low (rms 8)']
    assert abs(lo[2] - hi[2]) < 0.05 and lo[4] > 3 * lo[2] + 0.1           # without the correction the noise dominates S at low S/N
    assert hi[2] < 0.1 and lo[2] < 0.12


def test_smoothness_separates_clumpy_from_smooth():
    n, re_, q, pa = 1.0, 15, 0.8, 0
    smooth = gal(n, re_, q, pa, flux=0.8 * FTOT)
    rng = np.random.default_rng(5)
    clumpy = smooth.copy()
    yy, xx = np.indices(SHAPE)
    for k in range(10):
        r = rng.uniform(5, 28); t = rng.uniform(0, 2 * np.pi)
        cx, cy = XC + r * np.cos(t), YC + r * np.sin(t)
        clumpy += 0.2 * FTOT / 10 / (2 * np.pi * 1.5 ** 2) * np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / 1.5 ** 2)
    out = {}
    for name, img in (('smooth', smooth), ('clumpy', clumpy)):
        img = img + np.random.default_rng(1).normal(0, 1.0, img.shape)
        r = X.measure_object(img, XC, YC, q, 0.0, rms=1.0, rmax=150, seed=1)
        out[name] = r['smooth']
    print('MORPHEXT smoothness: smooth disk %.3f, same disk with 10 clumps (20 %% of the flux) %.3f' % (out['smooth'], out['clumpy']))
    assert out['smooth'] < 0.08 and out['clumpy'] > 0.15 and out['clumpy'] > 3 * out['smooth']


def test_petrosian_segmented_gini_m20_ordering_and_stability():
    r1 = X.measure_object(gal(1.0, 15, 1.0, 0), XC, YC, 1.0, 0.0, rms=1e-3, rmax=150)
    r4 = X.measure_object(gal(4.0, 12, 1.0, 0), XC, YC, 1.0, 0.0, rms=1e-3, rmax=150)
    img = gal(4.0, 12, 1.0, 0) + np.random.default_rng(2).normal(0, 3.0, SHAPE)
    r4n = X.measure_object(img, XC, YC, 1.0, 0.0, rms=3.0, rmax=150)
    print('MORPHEXT Gini_P / M20_P: n=1 %.3f / %.2f, n=4 %.3f / %.2f, n=4 + noise (peak S/N %.0f) %.3f / %.2f' % (r1['gini_p'], r1['m20_p'], r4['gini_p'], r4['m20_p'], gal(4.0, 12, 1.0, 0).max() / 3.0, r4n['gini_p'], r4n['m20_p']))
    assert r4['gini_p'] > r1['gini_p'] + 0.08 and r4['m20_p'] < r1['m20_p'] - 0.3
    assert abs(r4n['gini_p'] - r4['gini_p']) < 0.04 and abs(r4n['m20_p'] - r4['m20_p']) < 0.15


def test_masked_companion_and_edge_flags():
    n, re_, q, pa = 1.0, 15, 1.0, 0
    img = gal(n, re_, q, pa) + 0.0
    img += gal(1.0, 3, 1.0, 0, flux=0.3 * FTOT, xc=XC + 45, yc=YC + 10)      # bright compact companion in the outskirts
    mask = np.hypot(*np.meshgrid(np.arange(401) - (XC + 45), np.arange(401) - (YC + 10))) < 14
    truth = true_rp(n, re_, 0.2)
    r0 = X.measure_object(img, XC, YC, q, 0.0, rms=1e-3, rmax=150)
    r1 = X.measure_object(np.where(mask, 0, img), XC, YC, q, 0.0, rms=1e-3, rmax=150, mask=mask)
    print('MORPHEXT masked companion: R_P error unmasked %+.3f, masked %+.3f; F_P/F_tot unmasked %.3f masked %.3f (true 0.993)' % (r0['rp'][0.2] / truth - 1, r1['rp'][0.2] / truth - 1, r0['flux_p'] / FTOT, r1['flux_p'] / FTOT))
    assert abs(r1['rp'][0.2] / truth - 1) < 0.03 and abs(r0['rp'][0.2] / truth - 1) > 2 * abs(r1['rp'][0.2] / truth - 1)
    assert abs(r1['flux_p'] / FTOT - 0.993) < 0.02
    # near the image edge: no crash, nan Petrosian quantities are flagged
    small = gal(2.0, 8, 0.8, 30, shape=(60, 60), xc=10.0, yc=30.0)
    r = X.measure_object(small, 10.0, 30.0, 0.8, 0.5, rms=1e-3, rmax=150)
    assert r['flag'] != 0 or np.isfinite(r['rp'][0.2])
    r = X.measure_object(small, 2.0, 2.0, 0.8, 0.5, rms=1e-3)
    assert r['flag'] & 1


# --------------------------------------------------------------------------------------------------------------------- driver
def test_driver_on_synthetic_field(tmp_path):
    shape = (400, 400)
    rng = np.random.default_rng(3)
    specs = [(80, 80, 1.0, 10, 0.8, 0), (250, 90, 4.0, 8, 1.0, 0), (90, 270, 1.0, 12, 0.5, 60), (260, 280, 2.0, 9, 0.7, 120)]
    img = np.zeros(shape)
    flux = 4e4
    for x, y, n, re_, q, pa in specs:
        img += M.render_sersic(shape, x, y, M.sersic_Ie_from_flux(flux, re_, n, q), re_, n, q, pa)
    img += 50.0 + rng.normal(0, 1.0, shape)
    f = str(tmp_path / 'f.fits')
    fits.writeto(f, img.astype(np.float32), overwrite=True)
    rows = []
    for k, (x, y, n, re_, q, pa) in enumerate(specs):
        A = 0.6 * re_
        rows.append(dict(NUMBER=k + 1, X_IMAGE=x + 1, Y_IMAGE=y + 1, A_IMAGE=A, B_IMAGE=A * q, THETA_IMAGE=pa, KRON_RADIUS=3.5, FLUX_AUTO=flux, FLUX_RADIUS=re_ * 0.9))
    cat = str(tmp_path / 'c.tsv')
    tsvio.write_table(cat, ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'KRON_RADIUS', 'FLUX_AUTO', 'FLUX_RADIUS'], rows)
    wd = str(tmp_path / 'w')
    r = subprocess.run([PY, os.path.join(ROOT, 'ds9', 'library', 'ds9_morph_ext.py'), f, '--catalog', cat, '--work', wd, '--mag-zeropoint', '25', '--n-workers', '2', '--curves', '3'],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-1500:]
    L = r.stdout.strip().splitlines()
    head = L[0].split('\t')
    recs = [dict(zip(head, l.split('\t'))) for l in L[1:]]
    assert head[0] == 'NUMBER' and 'MX_RP' in head and 'MX_SMOOTH' in head and len(recs) == 4
    out = []
    for k, (x, y, n, re_, q, pa) in enumerate(specs):
        t = true_rp(n, re_, 0.2)
        m = float(recs[k]['MX_RP'])
        fp = 10 ** (-0.4 * (float(recs[k]['MX_MAG_P']) - 25.0)) / flux
        tf = float(X.sersic_fraction(2 * t / re_, n))
        out.append((m / t - 1, fp / tf - 1, float(recs[k]['MX_R50']) / (X.sersic_radius_fraction(0.5 * tf, n) * re_) - 1, float(recs[k]['MX_SMOOTH']), int(recs[k]['MX_FLAG'])))
    print('MORPHEXT driver (4 galaxies, noise 1, F=4e4): R_P rel err %s, F_P rel err %s, R50 rel err %s, S %s, flags %s' % (
        ['%+.3f' % o[0] for o in out], ['%+.3f' % o[1] for o in out], ['%+.3f' % o[2] for o in out], ['%.2f' % o[3] for o in out], [o[4] for o in out]))
    for o in out:
        assert abs(o[0]) < 0.06 and abs(o[1]) < 0.06 and abs(o[2]) < 0.06 and o[4] == 0 and o[3] < 0.2
    for nme in ('morph_ext_growth.tsv', 'morph_ext_curves.png'):
        assert os.path.getsize(os.path.join(wd, nme)) > 500
