"""GALFIT feedme / constraints import-export (ogfkit.galfitio), the multifit CLI on feedme input, and - when a GALFIT binary is available (env GALFIT_BIN [+ GALFIT_LD]) -
pixel-level agreement of the rendered models with GALFIT itself."""
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import galfitio as GI, models as M, multifit as MF

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
PY = sys.executable
GALFIT = os.environ.get('GALFIT_BIN') or shutil.which('galfit') or ''
GALFIT_OK = bool(GALFIT) and os.path.exists(GALFIT)

FEED = """
================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) data.fits           # Input data image (FITS file)
B) imgblock.fits       # Output data image block
C) none                # Sigma image name
D) psf.fits  none      # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data
F) none                # Bad pixel mask
G) cons.txt            # File with parameter constraints
H) 11    90   6    85  # Image region to fit (xmin xmax ymin ymax)
I) 61    61            # Size of the convolution box (x y)
J) 25.0                # Magnitude photometric zeropoint
K) 0.1  0.1            # Plate scale (dx dy)    [arcsec per pixel]
O) regular             # Display type
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
 0) sersic                 #  object type
 1) 48.5  41.25  1 0       #  position x, y
 3) 17.5     1             #  Integrated magnitude
 4) 6.0      1             #  R_e
 5) 2.5      0             #  Sersic index n
 9) 0.6      1             #  axis ratio (b/a)
10) 30.0     1             #  PA
 Z) 0
 0) expdisk
 1) 48.5  41.25  1 1
 3) 18.5     1
 4) 4.0      1
 9) 0.5      1
10) -20.0    1
 Z) 0
 0) sky
 1) 100.0    1
 2) 0.02     1
 3) -0.01    1
 Z) 0
"""
CONS = """
# comment
1_2   x   offset
1     re  2 to 12
2     q   -0.1 0.3
2     mag 17 to 19.5
"""


def test_parse_feedme_fields_and_conversions(tmp_path):
    (tmp_path / 'cons.txt').write_text(CONS)
    cfg = GI.parse_feedme(FEED, exptime=10.0, base_dir=str(tmp_path))
    assert cfg['bbox'] == [10, 90, 5, 85] and cfg['region'] == [11, 90, 6, 85] and cfg['conv_box'] == [61, 61] and cfg['zp'] == 25.0
    assert cfg['psf'].endswith('psf.fits') and cfg['mask'] == '' and cfg['sigma'] == ''
    a, b = cfg['components']
    assert a['kind'] == 'sersic' and (a['x'], a['y']) == (48.5, 41.25) and a['fixed'] == ['y', 'n'] and a['n'] == 2.5
    assert abs(a['mag'] - (17.5 - 2.5)) < 1e-9                                       # exposure time 10 s: 2.5 log10(10) = 2.5
    assert a['pa'] == 120.0 and b['pa'] == 70.0                                      # GALFIT PA + 90
    assert b['kind'] == 'exp' and abs(b['re'] - 4.0 * GI.RS_TO_RE) < 1e-9
    # sky: 3 free parameters -> plane; value moved from the GALFIT region centre (50.5, 45.5) to multifit's array centre (nx/2, ny/2 -> 1-based 51.0, 46.0)
    assert cfg['sky'] == 'plane' and cfg['sky_grad'] == [0.02, -0.01]
    assert abs(cfg['sky_value'] - (100.0 + 0.02 * (51.0 - 50.5) - 0.01 * (46.0 - 45.5))) < 1e-9
    assert cfg['tie'] == [['1.x', '0.x']]
    # constraints: absolute (to) and relative bounds, magnitude bound -> flux bound
    assert a['bounds']['re'] == [2.0, 12.0]
    assert abs(b['bounds']['q'][0] - 0.4) < 1e-9 and abs(b['bounds']['q'][1] - 0.8) < 1e-9
    f_lo, f_hi = b['bounds']['flux']
    assert abs(25.0 - 2.5 * math.log10(f_hi) - (17.0 - 2.5)) < 1e-9 and abs(25.0 - 2.5 * math.log10(f_lo) - (19.5 - 2.5)) < 1e-9


def test_unsupported_objects_and_offset_ratio_ties(tmp_path):
    bad = FEED.replace(' 0) expdisk', ' 0) fits')
    with pytest.raises(GI.FeedmeError):
        GI.parse_feedme(bad, base_dir=str(tmp_path))
    cfg = GI.parse_feedme(bad, strict=False, base_dir=str(tmp_path))
    assert len(cfg['components']) == 1 and any('skipped' in w for w in cfg['warnings'])
    (tmp_path / 'cons.txt').write_text('1_2 re offset\n1_2 pa ratio\n')
    cfg = GI.parse_feedme(FEED, base_dir=str(tmp_path))
    assert cfg['tie'] == [['1.re', '0.re', 'offset'], ['1.pa', '0.pa', 'ratio']] and not cfg['warnings']          # GALFIT semantics: offset / ratio of the start values
    with pytest.raises(GI.FeedmeError):
        GI.parse_feedme(FEED.replace(' 3) 17.5     1', ' 3)'), base_dir=str(tmp_path))


def test_gaussian_psf_and_devauc_mapping():
    t = FEED.replace(' 0) expdisk', ' 0) gaussian').replace(' 0) sersic', ' 0) devauc').replace(' 5) 2.5      0             #  Sersic index n\n', '')
    cfg = GI.parse_feedme(t, base_dir='.')
    a, b = cfg['components']
    assert a['kind'] == 'dev' and b['kind'] == 'sersic' and b['n'] == 0.5 and 'n' in b['fixed'] and abs(b['re'] - 2.0) < 1e-12      # FWHM 4 -> R_e = FWHM / 2


def test_feedme_roundtrip(tmp_path):
    (tmp_path / 'cons.txt').write_text(CONS)
    cfg = GI.parse_feedme(FEED, exptime=3.0, base_dir=str(tmp_path))
    text = GI.write_feedme(cfg, image='data.fits', psf='psf.fits', constraints='cons2.txt', shape=(80, 80))
    (tmp_path / 'cons2.txt').write_text(GI.write_constraints(cfg))
    cfg2 = GI.parse_feedme(text, exptime=3.0, base_dir=str(tmp_path))
    assert cfg2['region'] == cfg['region'] and cfg2['sky'] == 'plane' and cfg2['zp'] == cfg['zp']
    for c1, c2 in zip(cfg['components'], cfg2['components']):
        for k in ('x', 'y', 'mag', 're', 'q', 'pa'):
            assert abs(c1[k] - c2[k]) < 2e-4, k
        assert set(c1['fixed']) == set(c2['fixed']) and c1['kind'] == c2['kind']
    assert cfg2['tie'] == cfg['tie']
    assert abs(cfg2['sky_value'] - cfg['sky_value']) < 1e-4 and cfg2['sky_grad'] == cfg['sky_grad']
    for k, (lo, hi) in cfg['components'][0]['bounds'].items():
        assert np.allclose(cfg2['components'][0]['bounds'][k], (lo, hi), atol=1e-3)


def _truth_image(tmp_path, noise=1.0, seed=3):
    psf = M.moffat_psf(25, 3.0, beta=3.0); psf = psf / psf.sum()
    comps = [dict(kind='sersic', x=40.2, y=38.7, mag=17.0, re=7.0, n=2.0, q=0.65, pa=35.0 + 90.0)]
    img = MF.render_model([MF._normalise(c, 25.0) for c in comps], (81, 81), psf, 100.0) + np.random.default_rng(seed).normal(0, noise, (81, 81))
    fits.PrimaryHDU(img.astype(np.float32)).writeto(tmp_path / 'data.fits')
    fits.PrimaryHDU(psf.astype(np.float32)).writeto(tmp_path / 'psf.fits')
    return comps


def _cli(tmp_path, feed, extra=()):
    (tmp_path / 'in.feedme').write_text(feed)
    return subprocess.run([PY, os.path.join(PLUGIN, 'multifit.py'), '-', '--config', str(tmp_path / 'in.feedme'), '--work', str(tmp_path / 'w')] + list(extra), capture_output=True, text=True)


FEED_CLI = """A) data.fits
B) out.fits
C) none
D) psf.fits
E) 1
F) none
G) none
H) 1 81 1 81
I) 61 61
J) 25.0
K) 0.1 0.1
P) 0
 0) sersic
 1) 39.0 40.0 1 1
 3) 17.5 1
 4) 5.0 1
 5) 1.5 1
 9) 0.8 1
10) 0.0 1
 Z) 0
 0) sky
 1) 100.0 1
 2) 0 0
 3) 0 0
 Z) 0
"""


def test_cli_fits_from_feedme_and_exports_feedme(tmp_path):
    truth = _truth_image(tmp_path)
    r = _cli(tmp_path, FEED_CLI, ['--export-feedme', str(tmp_path / 'gf')])
    assert r.returncode == 0, r.stderr[-500:]
    out = open(tmp_path / 'gf' / 'galfit.feedme').read()
    cfg = GI.parse_feedme(out, strict=True, base_dir=str(tmp_path))
    c = cfg['components'][0]
    # 1-based image coordinates in the feedme: x = 0-based truth + 1
    assert abs(c['x'] - 41.2) < 0.1 and abs(c['y'] - 39.7) < 0.1 and abs(c['mag'] - 17.0) < 0.05 and abs(c['re'] / 7.0 - 1) < 0.06 and abs(c['n'] / 2.0 - 1) < 0.1
    assert abs(c['pa'] - 125.0) < 4.0 and abs(c['q'] - 0.65) < 0.03                       # GALFIT PA -55 = multifit 125 = truth 35 (+90 in the test's image) + 0
    assert cfg['region'] == [1, 81, 1, 81] and cfg['sky'] == 'const' and abs(cfg['sky_value'] - 100.0) < 0.3
    # JSON configs keep working
    js = dict(bbox=[0, 81, 0, 81], components=[dict(kind='sersic', x=39, y=40, mag=17.5, re=5, n=1.5, q=0.8, pa=0, bounds=dict(x=[35, 45]))], sky='const', zp=25.0)
    (tmp_path / 'c.json').write_text(json.dumps(js))
    r = subprocess.run([PY, os.path.join(PLUGIN, 'multifit.py'), str(tmp_path / 'data.fits'), '--config', str(tmp_path / 'c.json'), '--work', str(tmp_path / 'w2'), '--psf', str(tmp_path / 'psf.fits')],
                       capture_output=True, text=True)
    assert r.returncode == 0 and 'mag=17.0' in r.stdout, r.stdout[-400:] + r.stderr[-400:]


def test_cli_rejects_psf_oversampling(tmp_path):
    _truth_image(tmp_path)
    r = _cli(tmp_path, FEED_CLI.replace('E) 1', 'E) 3'))
    assert r.returncode != 0 and 'fine sampling' in (r.stderr + r.stdout)


def test_catalog_mode_exports_feedme_per_object(tmp_path):
    from ogfkit import synth
    psf = M.moffat_psf(25, 3.0, beta=3.0); psf = psf / psf.sum()
    comps = [dict(kind='sersic', x=50.0, y=50.0, mag=17.0, re=6.0, n=1.5, q=0.7, pa=40.0)]
    img = MF.render_model([MF._normalise(c, 25.0) for c in comps], (101, 101), psf, 100.0) + np.random.default_rng(1).normal(0, 1.0, (101, 101))
    fits.PrimaryHDU(img.astype(np.float32)).writeto(tmp_path / 'im.fits')
    (tmp_path / 'c.tsv').write_text('#\tNUMBER\n'.replace('#\tNUMBER\n', '') + 'NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tFLUX_AUTO\tFLUX_RADIUS\tCLASS_STAR\n7\t51\t51\t5\t3.5\t40\t%f\t6\t0.1\n' % (10 ** (-0.4 * (17.0 - 25.0))))
    r = subprocess.run([PY, os.path.join(PLUGIN, 'multifit.py'), str(tmp_path / 'im.fits'), '--catalog', str(tmp_path / 'c.tsv'), '--work', str(tmp_path / 'w'), '--psf', str(tmp_path / 'psf.fits'), '--model', 'sersic',
                        '--export-feedme', str(tmp_path / 'gf'), '--mag-zeropoint', '25', '--rms', '1.0'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    fe = tmp_path / 'gf' / 'galfit_7.feedme'
    assert fe.exists() and (tmp_path / 'gf' / 'galfit_7_psf.fits').exists()
    cfg = GI.parse_feedme(str(fe))
    c = cfg['components'][0]
    rows = [l.split('\t') for l in open(tmp_path / 'w' / 'multifit_results.tsv').read().splitlines()]
    mrow = dict(zip(rows[0], rows[1]))
    assert abs(float(mrow['MAG']) - c['mag']) < 2e-4 and abs(float(mrow['X']) - c['x']) < 2e-4 and abs(float(mrow['RE']) - c['re']) < 2e-4    # the feedme holds the fitted values
    assert cfg['psf'].endswith('galfit_7_psf.fits') and abs(c['mag'] - 17.0) < 0.5 and abs(c['x'] - 51.0) < 0.2 and cfg['region'][0] >= 1 and cfg['region'][1] <= 101


# --------------------------------------------------------------------------------------------------- against GALFIT itself
@pytest.mark.skipif(not GALFIT_OK, reason='no GALFIT binary (set GALFIT_BIN and, if needed, GALFIT_LD)')
def test_rendering_matches_galfit_pixelwise(tmp_path):
    ny, nx = 81, 93
    fits.PrimaryHDU(np.zeros((ny, nx), np.float32)).writeto(tmp_path / 'data.fits')
    y, x = np.mgrid[-12:13, -12:13]
    ps = np.exp(-0.5 * (x ** 2 + y ** 2) / 1.4 ** 2); ps /= ps.sum()
    fits.PrimaryHDU(ps.astype(np.float32)).writeto(tmp_path / 'psf.fits')
    objs = [(' 0) sersic', ' 1) 40.3 35.6 1 1', ' 3) 17.5 1', ' 4) 7.2 1', ' 5) 2.7 1', ' 9) 0.6 1', '10) 33.0 1'),
            (' 0) expdisk', ' 1) 55.0 45.2 1 1', ' 3) 18.2 1', ' 4) 4.1 1', ' 9) 0.4 1', '10) -50.0 1'),
            (' 0) devauc', ' 1) 20.5 60.1 1 1', ' 3) 18.8 1', ' 4) 3.3 1', ' 9) 0.8 1', '10) 80.0 1'),
            (' 0) gaussian', ' 1) 70.0 20.0 1 1', ' 3) 19.0 1', ' 4) 6.0 1', ' 9) 0.5 1', '10) 20.0 1'),
            (' 0) psf', ' 1) 30.2 20.7 1 1', ' 3) 19.5 1')]
    sky = (' 0) sky', ' 1) 12.0 1', ' 2) 0.02 1', ' 3) -0.03 1')
    head = 'A) data.fits\nB) out.fits\nC) none\nD) psf.fits\nE) 1\nF) none\nG) none\nH) 1 93 1 81\nI) 61 61\nJ) 25.0\nK) 0.1 0.1\nO) regular\nP) 1\n'
    env = dict(os.environ)
    if os.environ.get('GALFIT_LD'):
        env['LD_LIBRARY_PATH'] = os.environ['GALFIT_LD'] + ':' + env.get('LD_LIBRARY_PATH', '')
    worst = []
    for name, blocks in [('sersic', objs[:1]), ('expdisk', objs[1:2]), ('devauc', objs[2:3]), ('gaussian', objs[3:4]), ('psf', objs[4:5]), ('all+sky', objs + [sky])]:
        feed = head + '\n'.join('\n'.join(b) + '\n Z) 0' for b in blocks) + '\n'
        (tmp_path / 'm.feedme').write_text(feed)
        subprocess.run([GALFIT, 'm.feedme'], cwd=tmp_path, env=env, capture_output=True, stdin=subprocess.DEVNULL, timeout=120)
        g = fits.getdata(tmp_path / 'out.fits').astype(float)
        cfg = GI.parse_feedme(feed, base_dir=str(tmp_path))
        comps = [MF._normalise(c, cfg['zp']) for c in cfg['components']]
        for c in comps:
            c['x'] -= 1; c['y'] -= 1
        m = MF.render_model(comps, (ny, nx), psf=ps, sky=cfg['sky_value'] if cfg['sky'] != 'fixed' or name == 'all+sky' else 0.0, sky_grad=cfg['sky_grad'])
        rel = np.abs(m - g).max() / g.max()
        worst.append((name, rel, np.abs(m.sum() - g.sum()) / g.sum()))
        sys.stderr.write('GALFIT render %-8s max|dmodel|/peak %.4f  total-flux diff %.2e\n' % (name, rel, worst[-1][2]))
    assert max(w[1] for w in worst) < 0.02 and max(w[2] for w in worst) < 0.005             # measured: 1.4 % of peak (de Vaucouleurs cusp pixel), 0.36 % flux


def test_steep_profile_centred_on_a_pixel_is_not_overweighted():
    """Regression: an odd sub-pixel grid sampled the cusp at r = 0 exactly: n = 5 centred on a pixel had +16 % peak, +4 % flux (found against GALFIT)."""
    c = MF._normalise(dict(kind='sersic', x=40.0, y=40.0, mag=15.0, re=2.9, n=5.0, q=0.85, pa=110.0), 26.0)
    shape = (81, 81)
    img = M.render_sersic(shape, 40.0, 40.0, M.sersic_Ie_from_flux(c['flux'], 2.9, 5.0, 0.85), 2.9, 5.0, 0.85, 110.0)
    ns = 120
    o = (np.arange(ns) + 0.5) / ns - 0.5
    oy, ox = np.meshgrid(o, o, indexing='ij')
    Ie = M.sersic_Ie_from_flux(c['flux'], 2.9, 5.0, 0.85)
    y, x = np.mgrid[:81, :81].astype(float)
    ref = M._sersic_eval(x, y, 40.0, 40.0, Ie, 2.9, 5.0, 0.85, 110.0, 0.0)
    for j in range(34, 47):
        for i in range(34, 47):
            ref[j, i] = M._sersic_eval(i + ox.ravel(), j + oy.ravel(), 40.0, 40.0, Ie, 2.9, 5.0, 0.85, 110.0, 0.0).mean()
    sel = (np.abs(x - 40) <= 6) & (np.abs(y - 40) <= 6)
    assert np.abs(img - ref)[sel].max() / ref.max() < 0.01 and abs(img[sel].sum() / ref[sel].sum() - 1) < 0.005
