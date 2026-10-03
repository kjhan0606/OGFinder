"""Bar (Ferrers), ring (Gaussian), spiral (logarithmic arms) components, the presets and ogfkit.structfit.fit_structure."""
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import multifit as mf, profiles as PF, structfit as sf

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'multifit.py')
ZP = 25.0
FWHM = 3.0
RMS = 0.5


def total(c, n=301):
    c = mf._normalise(c, ZP)
    return float(mf.render_component(c, (n, n), None).sum())


def test_gring_flux_peak_radius_and_shape_modifiers():
    c = dict(kind='gring', x=150, y=150, flux=5000.0, rring=40.0, sring=5.0, q=0.6, pa=30)
    assert abs(total(c) / 5000.0 - 1) < 0.01
    img = mf.render_component(mf._normalise(c, ZP), (301, 301), None)
    yy, xx = np.mgrid[:301, :301]
    t = math.radians(30)
    xr = (xx - 150) * math.cos(t) + (yy - 150) * math.sin(t)
    yr = -(xx - 150) * math.sin(t) + (yy - 150) * math.cos(t)
    r = np.hypot(xr, yr / 0.6)
    prof = [img[(r >= k) & (r < k + 1)].mean() for k in range(0, 70)]
    assert abs(int(np.argmax(prof)) - 40) <= 1                      # peak at the ring radius in the generalised radius
    boxy = dict(c, c0=1.0)
    assert abs(total(boxy) / 5000.0 - 1) < 0.015                    # boxy ring: flux is conserved through flux_factor


def test_bar_and_spiral_conserve_flux():
    for c0 in (0.0, 1.0, 1.8):
        b = mf.bar_component(150, 150, 3000.0, 30.0, 40.0, q=0.25, c0=c0)
        assert abs(total(b) / 3000.0 - 1) < 0.02, c0
    sp = mf.spiral_disc(150, 150, 4000.0, 20.0, 0.9, 10.0, theta=-250.0, ampl=0.2)
    disc = dict(kind='exp', x=150, y=150, flux=4000.0, re=20.0, q=0.9, pa=10.0)
    assert abs(total(sp) / 4000.0 - 1) < 0.03
    a = mf.render_component(mf._normalise(sp, ZP), (301, 301), None)
    d = mf.render_component(mf._normalise(disc, ZP), (301, 301), None)
    assert np.abs(a - d).max() > 0.05 * d.max()                     # arms change the image
    assert np.allclose(a, a[::-1, ::-1], atol=1e-3 * a.max())      # two arms: symmetric under 180 deg rotation about the centre


def test_presets_are_valid_and_rendering_matches_flux():
    for name in ('bar', 'ring', 'spiral', 'bulge+disk+bar', 'bulge+disk+ring'):
        cs = [mf._normalise(c, ZP) for c in mf.preset(name, 100, 100, 5000.0, 12.0, 0.6, 30.0)]
        img = mf.render_model(cs, (201, 201), psf=mf.gaussian_psf_array(FWHM))
        assert 0.85 < img.sum() / 5000.0 < 1.05, name
        for t in mf.TIES.get(name, ()):
            assert int(t[0].split('.')[0]) < len(cs) and int(t[1].split('.')[0]) < len(cs)
    with pytest.raises(ValueError):
        mf.preset('nonsense', 0, 0, 1, 1)


def galaxy(kind, seed=3):
    rng = np.random.default_rng(seed)
    psf = mf.gaussian_psf_array(FWHM)
    bulge = dict(kind='sersic', x=50, y=50, flux=1800.0, re=3.0, n=2.0, q=0.85, pa=20)
    disc = dict(kind='exp', x=50, y=50, flux=7000.0, re=12.0, q=0.7, pa=20)
    comps = [bulge, disc]
    if kind == 'bar':
        disc['flux'] = 5500.0
        comps.append(mf.bar_component(50, 50, 2500.0, 13.0, 70.0, q=0.25, c0=1.0))
    elif kind == 'ring':
        disc['flux'] = 5500.0
        comps.append(dict(kind='gring', x=50, y=50, flux=2500.0, rring=12.0, sring=2.0, q=0.7, pa=20))
    elif kind == 'spiral':
        comps = [bulge, mf.spiral_disc(50, 50, 7000.0, 12.0, 0.85, 20.0, theta=-250.0, ampl=0.25)]
    img = mf.render_model([mf._normalise(c, ZP) for c in comps], (101, 101), psf=psf) + 10.0 + rng.normal(0, RMS, (101, 101))
    return img, psf


def run_struct(kind, features, seed=3):
    img, psf = galaxy(kind, seed)
    return sf.fit_structure(img, 50.6, 49.5, 9500.0, 9.0, 0.65, 14.0, psf=psf, rms=RMS, zp=ZP, features=features, max_nfev=250)


def test_fit_structure_recovers_bar():
    r = run_struct('bar', ('bar',))
    assert r['name'] == 'bar', (r['bic'], r['notes'])
    s = sf.bar_summary(r['res'])
    assert abs(((s['pa'] - 70.0 + 90) % 180) - 90) < 8.0
    assert abs(s['length'] / 13.0 - 1) < 0.15 and abs(s['bar_to_total'] - 0.2632) < 0.06
    assert r['bic']['base'] - r['bic']['bar'] > 100


def test_fit_structure_controls_stay_simple_and_ring_found():
    r = run_struct('control', ('bar', 'ring'))
    assert r['name'] == 'base', (r['bic'], r['notes'])
    r = run_struct('ring', ('ring',))
    assert r['name'] == 'ring', (r['bic'], r['notes'])
    c = r['res']['components'][2]
    assert abs(c['rring'] / 12.0 - 1) < 0.1 and abs(c['sring'] / 2.0 - 1) < 0.3


def test_cli_structure_and_presets(tmp_path):
    img, psf = galaxy('bar')
    fits.writeto(str(tmp_path / 'f.fits'), img.astype(np.float32), overwrite=True)
    with open(tmp_path / 'c.tsv', 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tKRON_RADIUS\tFLUX_AUTO\tFLUX_RADIUS\tMAG_AUTO\tCLASS_STAR\n')
        f.write('1\t51\t51\t9\t6\t14\t3.5\t9500\t9.0\t%g\t0.0\n' % (ZP - 2.5 * math.log10(9500.0)))
    base = [sys.executable, CLI, str(tmp_path / 'f.fits'), '--catalog', str(tmp_path / 'c.tsv'), '--psf-fwhm', str(FWHM), '--mag-zeropoint', str(ZP), '--rms', str(RMS), '--neighbours', 'mask',
            '--n-workers', '1', '--sky', 'const']
    p = subprocess.run(base + ['--work', str(tmp_path / 'w'), '--model', 'structure', '--columns', 'st', '--st-features', 'bar'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-600:]
    L = [l for l in p.stdout.strip().split('\n') if l and not l.startswith('#')]
    cols = L[0].split('\t')
    row = dict(zip(cols, L[1].split('\t')))
    import multifit as cli
    assert cols == ['NUMBER'] + cli.ST_COLUMNS
    assert row['ST_TYPE'] == '1' and abs(float(row['ST_BARLEN']) / 13 - 1) < 0.2 and float(row['ST_DBIC']) > 100 and 0.1 < float(row['ST_BARFRAC']) < 0.4
    p = subprocess.run(base + ['--work', str(tmp_path / 'w2'), '--model', 'bulge+disk+ring'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-400:]
    L = [l for l in p.stdout.strip().split('\n') if l and not l.startswith('#')]
    row = dict(zip(L[0].split('\t'), L[1].split('\t')))
    assert row['GF_NCOMP'] == '3'
