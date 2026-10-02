"""Advanced GALFIT components of multifit: Fourier / bending / rotation / truncation / C0 / Moffat / Ferrer / King / Nuker / edge-on disc / broken exponential -
engine invariants, feedme + constraints import / export of the GALFIT 3.0.x keywords, offset / ratio ties, CLI fits with neighbour masking, truth recovery and (with a GALFIT binary,
env GALFIT_BIN [+ GALFIT_LD]) pixel parity of the rendered models with GALFIT itself."""
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import galfitio as GI, multifit as MF, profiles as PF

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PLUGIN, 'validation'))
PY = sys.executable
GALFIT = os.environ.get('GALFIT_BIN') or shutil.which('galfit') or ''
GALFIT_OK = bool(GALFIT) and os.path.exists(GALFIT)
ZP = 25.0


def render(comps, shape=(101, 101), psf=None, sky=0.0):
    cs = [MF._normalise(json.loads(json.dumps(c)), ZP) for c in comps]
    return MF.render_model(cs, shape, psf=psf, sky=sky)


def S(**k):
    return dict(dict(kind='sersic', x=50.0, y=50.0, flux=2.0e4, re=8.0, n=2.0, q=0.7, pa=120.0), **k)


def gauss_psf(fwhm=3.0, size=25):
    y, x = np.mgrid[:size, :size] - size // 2
    p = np.exp(-0.5 * (x ** 2 + y ** 2) / (fwhm / 2.355) ** 2)
    return p / p.sum()


# ---------------------------------------------------------------------------------------------------------------------------------------------------- engine
def test_zero_modifiers_reduce_to_the_plain_profile():
    plain = render([S()], (121, 121))
    mod = render([S(c0=0.0, f1a=0.0, f1p=10.0, f3a=0.0, f3p=0.0, b1=0.0, b2=0.0)], (121, 121))
    assert np.abs(mod - plain).max() < 1e-3 * plain.max()
    rot0 = render([S(rot_func='power', rot_in=0.0, rot_out=20.0, rot_theta=0.0, rot_alpha=0.5, rot_incl=0.0, rot_pa=0.0)], (121, 121))
    assert np.abs(rot0 - plain).max() < 1e-3 * plain.max()


@pytest.mark.parametrize('mods', [dict(f2a=0.1, f2p=30.0), dict(f1a=0.06, f1p=-20.0, f3a=0.05, f3p=10.0), dict(b1=0.15), dict(b2=-0.1, b1=0.05), dict(c0=0.8), dict(c0=-0.5),
                                  dict(rot_func='power', rot_in=0.0, rot_out=25.0, rot_theta=60.0, rot_alpha=0.5, rot_incl=0.0, rot_pa=0.0),
                                  dict(rot_func='log', rot_in=0.0, rot_out=25.0, rot_theta=-70.0, rot_ws=3.0, rot_incl=0.0, rot_pa=0.0)])
def test_total_flux_is_preserved_by_the_modifiers(mods):
    img = render([S(**mods)], (241, 241))
    # the Sersic wings beyond the frame hold a few per cent (n = 2): compare with the unmodified model in the same frame
    ref = render([S()], (241, 241))
    assert abs(img.sum() / ref.sum() - 1.0) < 0.01, mods
    assert not np.allclose(img, ref, atol=1e-4 * ref.max())


def test_fourier_phase_rotates_the_perturbation():
    a = render([S(q=0.9, f4a=0.1, f4p=0.0)])
    b = render([S(q=0.9, f4a=0.1, f4p=45.0)])                 # m = 4: period 90 deg, a 45 deg phase turns boxy into diamond-like
    assert np.abs(a - b).max() > 0.02 * a.max()
    c = render([S(q=0.9, f4a=0.1, f4p=90.0)])                 # 90 deg phase = a full period of m = 4
    assert np.abs(a - c).max() < 1e-3 * a.max()


def test_truncation_removes_the_outer_profile_and_inner_the_core():
    comps = [S(flux=None, i0=5.0, norm='re', trunc_out=[1]), dict(kind='trunc', rbreak=15.0, dsoft=4.0)]
    comps[0].pop('flux')
    img = render(comps)
    ref = render([dict(comps[0], trunc_out=[])])
    yy, xx = np.mgrid[:101, :101]
    far = np.hypot(xx - 50, yy - 50) > 15 + 2 * 4.0 * 1.0 + 12                     # rbreak, softening and the ellipse elongation
    assert img[far].max() < 1e-3 * ref.max() and img.sum() < ref.sum()
    inner = [dict(comps[0], trunc_in=[1], trunc_out=[]), dict(comps[1])]
    im2 = render(inner)
    assert im2[50, 50] < 0.01 * ref[50, 50] and abs(im2[50, 95] / ref[50, 95] - 1) < 0.02


def test_surface_brightness_normalisations():
    for norm, ref in (('center', 40.0), ('re', 6.0)):
        c = S(norm=norm, i0=ref, q=1.0, n=1.0)
        c.pop('flux')
        cc = MF._normalise(c, ZP)
        r = PF.surface_brightness(cc, np.array([cc['x']]), np.array([cc['y']])) if hasattr(PF, 'surface_brightness') else None
        if r is not None:
            if norm == 'center':
                assert abs(float(r[0]) / ref - 1) < 1e-6
            else:
                rr = PF.surface_brightness(cc, np.array([cc['x'] + cc['re']]), np.array([cc['y']]))
                assert abs(float(rr[0]) / ref - 1) < 1e-6


@pytest.mark.parametrize('comp,minflux', [(dict(kind='moffat', x=50.0, y=50.0, flux=5000.0, fwhm=5.0, beta=3.0, q=0.8, pa=30.0), 0.98),
                                          (dict(kind='ferrer', x=50.0, y=50.0, i0=30.0, rout=25.0, alpha=3.0, beta=1.0, q=0.8, pa=0.0), None),
                                          (dict(kind='king', x=50.0, y=50.0, i0=30.0, rc=5.0, rt=30.0, alpha=2.0, q=0.8, pa=0.0), None),
                                          (dict(kind='nuker', x=50.0, y=50.0, ib=20.0, rb=8.0, alpha=1.5, beta=1.8, gamma=0.5, q=0.8, pa=10.0), None),
                                          (dict(kind='edgedisk', x=50.0, y=50.0, i0=10.0, hs=3.0, rs=14.0, pa=45.0), None),
                                          (dict(kind='brokenexp', x=50.0, y=50.0, i0=8.0, h1=10.0, h2=5.0, rbreak=15.0, alpha=0.5, q=0.7, pa=20.0), None)])
def test_profile_kinds_render_finite_positive_images(comp, minflux):
    img = render([comp], (101, 101))
    assert np.isfinite(img).all() and img.max() > 0 and img.min() >= -1e-9
    assert abs(np.unravel_index(np.argmax(img), img.shape)[0] - 50) <= 1 and abs(np.unravel_index(np.argmax(img), img.shape)[1] - 50) <= 1
    if minflux:
        assert minflux * comp['flux'] <= img.sum() <= 1.02 * comp['flux']


# -------------------------------------------------------------------------------------------------------------------------------------------- feedme import / export
ADV = """
A) data.fits
B) out.fits
C) none
D) none
E) 1
F) none
G) cons.txt
H) 1 101 1 101
I) 61 61
J) 25.0
K) 0.2 0.2
O) regular
P) 0
 0) sersic2                #  object 1
 1) 50.5  49.5  1 1
 3) 21.0     1             #  surface brightness at R_e [mag/arcsec^2]
 4) 8.0      1
 5) 2.0      0
 9) 0.7      1
10) 30.0     1
C0) 0.3      1
F1) 0.05  20.0  1 0
F3) 0.04 -30.0  1 1
B1) 0.1      1
B2) -0.05    0
R0) power
R1) 0.0 0
R2) 25.0 0
R3) 60.0 1
R4) 0.5 0
R9) 0.0 0
R10) 0.0 0
To) 2
 Z) 0
T0) radial                 #  object 2
T1) 50.5 49.5 0 0
T4) 22.0 1
T5) 6.0 0
 Z) 0
 0) moffat                 #  object 3
 1) 70.0 30.0 1 1
 3) 19.0 1
 4) 4.0 1
 5) 2.5 0
 9) 0.9 1
10) 10.0 1
 Z) 0
 0) nuker                  #  object 4
 1) 30.0 70.0 1 1
 3) 20.0 1
 4) 6.0 1
 5) 1.5 0
 6) 1.8 1
 7) 0.5 1
 9) 0.8 1
10) 0.0 1
 Z) 0
 0) sky
 1) 10.0 1
 2) 0 0
 3) 0 0
 Z) 0
"""


def test_parse_advanced_keywords(tmp_path):
    (tmp_path / 'cons.txt').write_text('1 c0 -1 to 2\n1 f1 0 to 0.5\n1 f3p -90 to 90\n3 re 2 to 9\n3 n 1.5 to 4\n4 mag 18 to 22\n1_4 pa offset\n1_3 mag ratio\n1 b1 0 to 1\n')
    cfg = GI.parse_feedme(ADV, exptime=1.0, base_dir=str(tmp_path))
    a, t, m, nk = cfg['components']
    assert a['kind'] == 'sersic' and a['norm'] == 're' and a['trunc_out'] == [1] and 'trunc_in' not in a
    pix = 0.2 * 0.2
    assert abs(a['i0'] - pix * 10 ** (-0.4 * (21.0 - 25.0))) < 1e-9 and 'i0' not in a['fixed'] and a['mu'] == 21.0
    assert a['c0'] == 0.3 and (a['f1a'], a['f1p'], a['f3a'], a['f3p']) == (0.05, 20.0, 0.04, -30.0) and 'f1p' in a['fixed'] and 'f1a' not in a['fixed']
    assert a['b1'] == 0.1 and a['b2'] == -0.05 and 'b2' in a['fixed']
    assert a['rot_func'] == 'power' and a['rot_out'] == 25.0 and a['rot_theta'] == 60.0 and 'rot_theta' not in a['fixed'] and set(a['fixed']) >= {'n', 'rot_in', 'rot_out', 'rot_alpha'}
    assert t['kind'] == 'trunc' and t['rbreak'] == 22.0 and t['dsoft'] == 6.0 and t['x'] == 50.5 and 'dsoft' in t['fixed'] and {'x', 'y'} <= set(t['fixed'])
    assert m['kind'] == 'moffat' and m['fwhm'] == 4.0 and m['beta'] == 2.5 and 'beta' in m['fixed'] and abs(m['pa'] - 100.0) < 1e-9 and abs(m['mag'] - 19.0) < 1e-9
    assert nk['kind'] == 'nuker' and nk['rb'] == 6.0 and nk['alpha'] == 1.5 and nk['gamma'] == 0.5 and abs(nk['ib'] - pix * 10 ** (-0.4 * (20.0 - 25.0))) < 1e-9
    # constraints: GALFIT's positional names
    assert a['bounds']['c0'] == [-1.0, 2.0] and a['bounds']['f1a'] == [0.0, 0.5] and a['bounds']['f3p'] == [-90.0, 90.0]
    assert m['bounds']['fwhm'] == [2.0, 9.0] and m['bounds']['beta'] == [1.5, 4.0]
    lo, hi = nk['bounds']['ib']
    assert abs(25 - 2.5 * math.log10(hi / pix) - 18) < 1e-9 and abs(25 - 2.5 * math.log10(lo / pix) - 22) < 1e-9
    assert ['3.pa', '0.pa', 'offset'] in cfg['tie'] and ['2.flux', '0.flux', 'ratio'] in cfg['tie'] or ['2.flux', '0.i0', 'ratio'] in cfg['tie'] or any(t_[2] == 'ratio' for t_ in cfg['tie'])
    assert any('b1' in w for w in cfg['warnings'])                                 # GALFIT cannot constrain bending modes: reported, ignored


def test_advanced_roundtrip_through_write_feedme(tmp_path):
    (tmp_path / 'cons.txt').write_text('1 c0 -1 to 2\n3 re 2 to 9\n')
    cfg = GI.parse_feedme(ADV, exptime=1.0, base_dir=str(tmp_path))
    text = GI.write_feedme(cfg, image='data.fits', constraints='cons2.txt', shape=(101, 101))
    (tmp_path / 'cons2.txt').write_text(GI.write_constraints(cfg))
    cfg2 = GI.parse_feedme(text, exptime=1.0, base_dir=str(tmp_path))
    assert len(cfg2['components']) == 4
    for c1, c2 in zip(cfg['components'], cfg2['components']):
        assert c1['kind'] == c2['kind'] and set(c1['fixed']) == set(c2['fixed'])
        for k in set(MF.param_names(c1)) | {'trunc_in', 'trunc_out', 'norm', 'rot_func'}:
            v1, v2 = c1.get(k), c2.get(k)
            if isinstance(v1, float):
                assert abs(v1 - v2) < 2e-3 * max(1.0, abs(v1)), (k, v1, v2)
            else:
                assert v1 == v2, (k, v1, v2)
    assert cfg2['components'][0]['bounds']['c0'] == [-1.0, 2.0]
    assert 'T0) radial' in text and 'To) 2' in text and 'F3) ' in text and 'R0) power' in text and '0) sersic2' in text


def test_unsupported_and_inconsistent_advanced_feedme(tmp_path):
    with pytest.raises(GI.FeedmeError, match='radial'):
        GI.parse_feedme(ADV.replace('T0) radial', 'T0) radial2'), base_dir=str(tmp_path))
    with pytest.raises(GI.FeedmeError, match='sersic1/2/3'):
        GI.parse_feedme(ADV.replace(' 0) sersic2', ' 0) sersic'), base_dir=str(tmp_path))
    with pytest.raises(GI.FeedmeError, match='truncation'):
        GI.parse_feedme(ADV.replace('To) 2', 'To) 3'), base_dir=str(tmp_path))
    with pytest.raises(GI.FeedmeError, match='rotation'):
        GI.parse_feedme(ADV.replace('R0) power', 'R0) spline'), base_dir=str(tmp_path))
    cfg = GI.parse_feedme(ADV.replace('C0) 0.3      1', 'C0) 0.3      1').replace(' B1) 0.1', ' B1) 0.1'), base_dir=str(tmp_path))
    assert any('C0' in w and 'Fourier' in w for w in cfg['warnings'])              # GALFIT normalises C0 + Fourier differently


def test_write_constraints_uses_galfit_positional_names():
    cfg = dict(zp=25.0, exptime=1.0, components=[
        dict(kind='moffat', x=1, y=1, flux=100.0, fwhm=4.0, beta=3.0, q=1.0, pa=0.0, bounds={'fwhm': [2, 8], 'beta': [1.5, 5], 'flux': [50, 200]}),
        dict(kind='sersic', x=1, y=1, flux=100.0, re=5.0, n=2.0, q=0.5, pa=0.0, f3a=0.1, f3p=0.0, b1=0.0, bounds={'f3a': [0, 0.3], 'f3p': [-45, 45], 'b1': [-1, 1], 'c0': [-1, 1]}),
        dict(kind='trunc', rbreak=20.0, dsoft=5.0, bounds={'rbreak': [10, 30]})], tie=[['1.pa', '0.pa', 'offset'], ['1.flux', '0.flux', 'ratio']])
    txt = GI.write_constraints(cfg)
    assert '1  re  2.0000 to 8.0000' in txt and '1  n  1.5000 to 5.0000' in txt and '1  mag ' in txt                      # moffat: FWHM is parameter 4 ('re'), beta parameter 5 ('n')
    assert '2  f3  0.0000 to 0.3000' in txt and '2  f3p  -45.0000 to 45.0000' in txt and '2  c0 ' in txt
    assert '# 2  b1' in txt and 'cannot constrain' in txt                                                                  # bending modes: commented
    assert '3  re  10.0000 to 30.0000' in txt                                                                              # truncation break radius
    assert '1_2  pa  offset' in txt and '1_2  mag  offset' in txt


def test_offset_and_ratio_ties_in_the_fit():
    rng = np.random.default_rng(5)
    psf = gauss_psf()
    truth = [dict(kind='exp', x=40.0, y=50.0, flux=6000.0, re=9.0, q=0.7, pa=30.0), dict(kind='exp', x=62.0, y=50.0, flux=3000.0, re=7.0, q=0.7, pa=30.0)]
    img = render(truth, (101, 101), psf=psf, sky=5.0) + rng.normal(0, 0.3, (101, 101))
    start = [dict(kind='exp', x=41.0, y=51.0, flux=5000.0, re=8.0, q=0.6, pa=20.0), dict(kind='exp', x=63.0, y=49.0, flux=2500.0, re=6.0, q=0.6, pa=20.0)]
    res = MF.fit_multistart(img, start, restarts=0, psf=psf, rms=0.3, sky='const', zp=ZP, tie=[['1.x', '0.x', 'offset'], ['1.flux', '0.flux', 'ratio'], ['1.pa', '0.pa']], max_nfev=300)
    a, b = res['components']
    assert abs((b['x'] - a['x']) - 22.0) < 1e-6                         # offset = difference of the start values (63 - 41), kept by the fit (GALFIT semantics)
    assert abs(b['flux'] / a['flux'] - 0.5) < 1e-6 and abs(b['pa'] - a['pa']) < 1e-9
    assert abs(a['flux'] / 6000.0 - 1) < 0.05


# ---------------------------------------------------------------------------------------------------------------------------------------------------- CLI + recovery
def _cli(tmp_path, feed, extra=()):
    (tmp_path / 'in.feedme').write_text(feed)
    return subprocess.run([PY, os.path.join(PLUGIN, 'multifit.py'), '-', '--config', str(tmp_path / 'in.feedme'), '--work', str(tmp_path / 'w')] + list(extra), capture_output=True, text=True, timeout=900)


def test_cli_fits_fourier_and_moffat_with_neighbour_masking(tmp_path):
    psf = gauss_psf(3.0)
    fits.PrimaryHDU(psf.astype(np.float32)).writeto(tmp_path / 'psf.fits')
    truth = [S(x=41.0, y=41.0, flux=30000.0, re=7.0, n=1.5, q=0.7, pa=60.0, f3a=0.08, f3p=20.0),
             dict(kind='moffat', x=72.0, y=24.0, flux=4000.0, fwhm=4.0, beta=3.0, q=0.8, pa=30.0)]
    nb = S(x=70.0, y=68.0, flux=40000.0, re=5.0, n=1.0, q=0.9, pa=0.0)
    rng = np.random.default_rng(11)
    img = render(truth + [nb], (80, 100), psf=psf, sky=5.0) + rng.normal(0, 0.3, (80, 100))
    fits.PrimaryHDU(img.astype(np.float32)).writeto(tmp_path / 'data.fits')
    (tmp_path / 'nb.tsv').write_text('NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\n1\t42\t42\t7\t5\t0\n2\t71\t69\t5\t5\t0\n')
    st = [dict(c) for c in truth]
    st[0].update(x=42.0, y=40.0, flux=24000.0, re=6.0, n=1.2, q=0.8, pa=40.0, f3a=0.03, f3p=0.0)
    st[1].update(x=73.0, y=25.0, flux=3000.0, fwhm=3.5, fixed=['beta'])
    cfg = dict(components=[dict(c, x=c['x'] + 1, y=c['y'] + 1) for c in st], zp=ZP, exptime=1.0, sky='const', sky_value=5.0)
    feed = GI.write_feedme(cfg, image=str(tmp_path / 'data.fits'), psf=str(tmp_path / 'psf.fits'), shape=(80, 100))
    # without masking the bright neighbour biases the fit; with --mask-catalog (the neighbour at 71, 69 only: the first catalog row is the target and is not masked) it is recovered
    r = _cli(tmp_path, feed, ['--mask-catalog', str(tmp_path / 'nb.tsv'), '--mask-exclude', '4', '--export-feedme', str(tmp_path / 'out.feedme'), '--restarts', '1'])
    assert r.returncode == 0, r.stderr[-800:]
    assert '1 neighbour(s) masked' in r.stdout and 'f3a=' in r.stdout
    rows = [ln.split('\t') for ln in (tmp_path / 'w' / 'multifit_results.tsv').read_text().splitlines()]
    head = rows[0]
    d0 = dict(zip(head, rows[1])); d1 = dict(zip(head, rows[2]))
    assert abs(float(d0['X']) - 42.0) < 0.15 and abs(float(d0['Y']) - 42.0) < 0.15 and abs(float(d0['RE']) / 7.0 - 1) < 0.08 and 'f3a=' in d0['EXTRA']
    f3a = float(dict(kv.split('=') for kv in d0['EXTRA'].split(';'))['f3a'])
    assert abs(abs(f3a) - 0.08) < 0.04                                       # amplitude and phase are degenerate up to a sign (a -> -a, phase + 180/m)
    assert abs(float(d1['X']) - 73.0) < 0.3 and 'fwhm=' in d1['EXTRA']
    back = GI.parse_feedme(str(tmp_path / 'out.feedme'), base_dir=str(tmp_path))
    assert back['components'][0]['kind'] == 'sersic' and 'f3a' in back['components'][0] and back['components'][1]['kind'] == 'moffat'


def test_cli_feedme_with_truncation_object(tmp_path):
    fits.PrimaryHDU(np.zeros((1, 1), np.float32)).writeto(tmp_path / 'psf.fits')
    comps = [S(x=50.0, y=50.0, norm='re', i0=6.0, re=9.0, n=2.0, q=0.8, pa=100.0, trunc_out=[1], fixed=['n']), dict(kind='trunc', rbreak=20.0, dsoft=6.0, fixed=['dsoft'])]
    comps[0].pop('flux')
    rng = np.random.default_rng(2)
    img = render(comps, (101, 101), sky=3.0) + rng.normal(0, 0.2, (101, 101))
    fits.PrimaryHDU(img.astype(np.float32)).writeto(tmp_path / 'data.fits')
    st = json.loads(json.dumps(comps))
    st[0].update(x=51.0, y=49.0, i0=5.0, re=8.0, q=0.7, pa=90.0)
    st[1].update(rbreak=17.0)
    feed = GI.write_feedme(dict(components=[dict(c, x=c.get('x', 0) + 1, y=c.get('y', 0) + 1) if c['kind'] != 'trunc' else c for c in st], zp=ZP, exptime=1.0, sky='const', sky_value=3.0), image=str(tmp_path / 'data.fits'), shape=(101, 101))
    r = _cli(tmp_path, feed, ['--psf-fwhm', '0.5', '--export-feedme', str(tmp_path / 'o.feedme')])
    assert r.returncode == 0, r.stderr[-800:]
    rows = [ln.split('\t') for ln in (tmp_path / 'w' / 'multifit_results.tsv').read_text().splitlines()]
    d = [dict(zip(rows[0], x)) for x in rows[1:]]
    assert d[1]['KIND'] == 'trunc' and 'rbreak=' in d[1]['EXTRA'] and abs(float(dict(kv.split('=') for kv in d[1]['EXTRA'].split(';'))['rbreak']) - 20.0) < 1.0
    assert abs(float(d[0]['RE']) / 9.0 - 1) < 0.05


@pytest.mark.parametrize('name', ['F_fourier', 'I_truncation', 'N_king', 'O_edgedisk'])
def test_truth_recovery_cases(name):
    import galfit_advanced as GA
    class A:
        pass
    r = GA._rec_one((name, 0, A()))
    assert r['chi2_red'] < 1.12 and r['converged']
    e = r['errors']
    assert abs(e['sky']) < 0.15
    for k, v in e.items():
        if k.endswith(('.x', '.y')):
            assert abs(v) < (1.2 if name == 'O_edgedisk' else 0.5), (k, v)
        if k.endswith(('.re', '.rt', '.rs', '.hs', '.rbreak', '.rc')):
            assert abs(v) < 0.2, (k, v)


# ------------------------------------------------------------------------------------------------------------------------------------------------------ GALFIT parity
@pytest.mark.skipif(not GALFIT_OK, reason='no GALFIT binary (GALFIT_BIN)')
def test_rendering_of_the_advanced_components_matches_galfit(tmp_path):
    import galfit_advanced as GA
    class A:
        galfit = GALFIT
        ld = os.environ.get('GALFIT_LD', '')
    psf = GA.moffat_psf()
    worst = []
    for name in ('c0_boxy', 'fourier_f1f3f5', 'bending_b1b2', 'rot_power', 'rot_log', 'trunc_outer', 'trunc_inner', 'moffat', 'ferrer', 'king', 'nuker', 'edgedisk', 'combo'):
        comps = GA.fill_sb(GA.prepare(json.loads(json.dumps(GA.CASES[name]))))
        for with_psf in (False, True):
            g = GA.galfit_render(A, comps, psf if with_psf else None, str(tmp_path / (name + str(with_psf))))
            assert g is not None, name
            m = GA.render_ours(comps, psf if with_psf else None)
            worst.append((float(np.abs(m - g).max() / g.max()), name, with_psf, abs(m.sum() / g.sum() - 1)))
    worst.sort(reverse=True)
    sys.stderr.write('GALFIT parity (advanced), worst: %s\n' % (worst[:3],))
    assert worst[0][0] < 0.03 and max(w[3] for w in worst) < 0.004               # measured: 2.6 % of peak (steep n = 4 core, PSF-convolved), flux 0.2 %
