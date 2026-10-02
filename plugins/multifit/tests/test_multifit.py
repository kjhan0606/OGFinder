"""Multi-component fitting: parameter recovery with known truth (single / blends vs separation / bulge+disk / PSF+host / spatially varying PSF), pull calibration,
constraints, shared sky, model selection, and the catalog / config CLI."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import models as M, multifit as MF, psfmodel as PM, synth, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
PY = sys.executable
ZP = 25.0
PSF = M.moffat_psf(25, 3.0, beta=3.0)
PSF = PSF / PSF.sum()


def image(comps, shape=(81, 81), sky=100.0, noise=1.5, seed=1, psf=PSF, grad=(0.0, 0.0)):
    img = MF.render_model([MF._normalise(c, ZP) for c in comps], shape, psf, sky, grad)
    return img + np.random.default_rng(seed).normal(0, noise, shape)


def test_single_sersic_recovery_over_profiles():
    rows = []
    for n, q, pa, re in ((1.0, 0.6, 20, 6), (2.0, 0.8, 100, 4), (4.0, 0.7, 60, 5)):
        truth = dict(kind='sersic', x=40.3, y=39.6, mag=17.0, re=re, n=n, q=q, pa=pa)
        d = image([truth], noise=1.0, seed=int(n))
        start = dict(kind='sersic', x=39.0, y=41.0, mag=17.6, re=0.7 * re, n=1.5, q=0.9, pa=0)
        r = MF.fit(d, [start], psf=PSF, rms=1.0, zp=ZP)
        c = r['components'][0]
        rows.append((n, c['mag'] - 17.0, c['re'] / re - 1, c['n'] / n - 1, c['q'] - q, ((c['pa'] - pa + 90) % 180) - 90, np.hypot(c['x'] - 40.3, c['y'] - 39.6), r['chi2_red']))
        print('MULTIFIT single n=%.0f: dmag %+.3f dre/re %+.3f dn/n %+.3f dq %+.3f dpa %+.1f dpos %.3f chi2 %.2f' % ((n,) + rows[-1][1:]))
        assert r['converged'] and r['flags'] == 0
    for n, dm, dre, dn, dq, dpa, dpos, chi in rows:
        tol = (0.05, 0.06, 0.12) if n < 3 else (0.10, 0.16, 0.12)          # n = 4: the fitted sky absorbs part of the extended wings (fixed sky: exact), see docs
        assert abs(dm) < tol[0] and abs(dre) < tol[1] and abs(dn) < tol[2] and abs(dq) < 0.03 and abs(dpa) < 6 and dpos < 0.1 and 0.85 < chi < 1.2


def test_error_calibration_pulls():
    truth = dict(kind='sersic', x=40, y=40, mag=18.0, re=5, n=2.0, q=0.7, pa=30)
    pulls = {k: [] for k in ('mag', 're', 'n', 'q', 'x')}
    tr = dict(mag=18.0, re=5.0, n=2.0, q=0.7, x=40.0)
    for seed in range(40):
        d = image([truth], noise=2.0, seed=100 + seed)
        r = MF.fit(d, [dict(kind='sersic', x=40.5, y=39.5, mag=18.3, re=4, n=1.5, q=0.8, pa=0)], psf=PSF, rms=2.0, zp=ZP)
        c = r['components'][0]
        for k in pulls:
            pulls[k].append((c[k] - tr[k]) / c['errors'][k])
    s = {k: (np.mean(v), np.std(v)) for k, v in pulls.items()}
    print('MULTIFIT pulls (mean, std): ' + ', '.join('%s %+.2f %.2f' % (k, *v) for k, v in s.items()))
    for k, (m, sd) in s.items():
        assert abs(m) < 0.5 and 0.7 < sd < 1.4, (k, m, sd)


def test_blend_recovery_vs_separation_and_baselines():
    """Simultaneous fit vs masking / ignoring the neighbour, as a function of separation (in units of re of the target)."""
    A = dict(kind='sersic', x=40, y=40, mag=17.0, re=5, n=2.0, q=0.7, pa=30)
    res = {}
    for sep in (5, 8, 14, 24):
        out = {'sim': [], 'mask': [], 'ignore': []}
        for seed in range(5):
            bx, by = 40 + sep * np.cos(0.6), 40 + sep * np.sin(0.6)
            B = dict(kind='sersic', x=bx, y=by, mag=17.6, re=4, n=1.0, q=0.8, pa=100)
            d = image([A, B], seed=seed)
            sa = dict(kind='sersic', x=40.5, y=39.5, mag=17.2, re=4, n=1.5, q=0.9, pa=0)
            sb = dict(kind='sersic', x=bx + 0.7, y=by - 0.7, mag=18.0, re=3, n=1.5, q=0.9, pa=0, bounds=dict(x=(bx - 3, bx + 3), y=(by - 3, by + 3)))
            r = MF.fit(d, [sa, sb], psf=PSF, rms=1.5, zp=ZP)
            ca = min(r['components'], key=lambda k: np.hypot(k['x'] - 40, k['y'] - 40))
            out['sim'].append((ca['mag'] - 17.0, ca['re'] / 5 - 1))
            yy, xx = np.mgrid[:81, :81]
            m = np.hypot(xx - bx, yy - by) < 1.8 * 4
            r = MF.fit(d, [sa], psf=PSF, rms=1.5, zp=ZP, mask=m)
            out['mask'].append((r['components'][0]['mag'] - 17.0, r['components'][0]['re'] / 5 - 1))
            r = MF.fit(d, [sa], psf=PSF, rms=1.5, zp=ZP)
            out['ignore'].append((r['components'][0]['mag'] - 17.0, r['components'][0]['re'] / 5 - 1))
        res[sep] = {k: (np.mean(np.abs(np.array(v)[:, 0])), np.mean(np.abs(np.array(v)[:, 1]))) for k, v in out.items()}
        print('MULTIFIT blend sep %2d px (%.1f re): |dmag| simultaneous %.3f masked %.3f ignored %.3f ; |dre/re| %.3f %.3f %.3f' % (
            sep, sep / 5, res[sep]['sim'][0], res[sep]['mask'][0], res[sep]['ignore'][0], res[sep]['sim'][1], res[sep]['mask'][1], res[sep]['ignore'][1]))
    for sep in (5, 8, 14, 24):
        assert res[sep]['sim'][0] < 0.15 and res[sep]['sim'][1] < (0.2 if sep < 8 else 0.1)
    assert res[8]['sim'][0] < 0.5 * res[8]['ignore'][0] and res[14]['sim'][0] < res[14]['ignore'][0]
    assert res[8]['sim'][0] < res[8]['mask'][0] + 0.02


def test_bulge_disk_and_psf_host():
    comps = [dict(kind='dev', x=40, y=40, mag=18.2, re=2.5, q=0.9, pa=0), dict(kind='exp', x=40, y=40, mag=17.5, re=7.0, q=0.6, pa=50)]
    bt_true = 10 ** (-0.4 * 18.2) / (10 ** (-0.4 * 18.2) + 10 ** (-0.4 * 17.5))
    d = image(comps, noise=0.6, seed=3)
    r, _ = MF.fit_preset(d, 'bulge+disk', 41, 39, 10 ** (-0.4 * (17.0 - ZP)), 5.0, 0.7, 30.0, psf=PSF, rms=0.6, zp=ZP)
    b, dk = r['components']
    bt = b['flux'] / (b['flux'] + dk['flux'])
    print('MULTIFIT bulge+disk: B/T %.3f (true %.3f), disk re %.2f (7.0) q %.2f (0.6), bulge re %.2f (2.5), centre tied dx=%.4f chi2 %.2f' % (bt, bt_true, dk['re'], dk['q'], b['re'], dk['x'] - b['x'], r['chi2_red']))
    assert abs(bt - bt_true) < 0.07 and abs(dk['re'] / 7 - 1) < 0.1 and abs(dk['x'] - b['x']) < 1e-12 and abs(dk['q'] - 0.6) < 0.05
    # AGN + host
    comps = [dict(kind='psf', x=40.2, y=39.8, mag=19.0), dict(kind='sersic', x=40.2, y=39.8, mag=17.5, re=6.0, n=2.0, q=0.75, pa=70)]
    d = image(comps, noise=0.8, seed=4)
    r, _ = MF.fit_preset(d, 'psf+sersic', 40, 40, 10 ** (-0.4 * (17.3 - ZP)), 4.0, 0.9, 0.0, psf=PSF, rms=0.8, zp=ZP)
    p, h = r['components']
    print('MULTIFIT psf+host: psf mag %.3f (19.0) host mag %.3f (17.5) re %.2f (6.0) n %.2f (2.0)' % (p['mag'], h['mag'], h['re'], h['n']))
    assert abs(p['mag'] - 19.0) < 0.15 and abs(h['mag'] - 17.5) < 0.1 and abs(h['re'] / 6 - 1) < 0.12


def test_constraints_fixed_bounds_tie_and_sky_plane():
    truth = dict(kind='sersic', x=40, y=40, mag=17.0, re=5, n=2.0, q=0.7, pa=30)
    d = image([truth], noise=1.0, seed=7, grad=(0.05, -0.03))
    # fixed n: stays exactly at the start, errors absent, other parameters still recover
    r = MF.fit(d, [dict(kind='sersic', x=40, y=40, mag=17.2, re=4, n=2.0, q=0.8, pa=0, fixed=['n'])], psf=PSF, rms=1.0, zp=ZP, sky='plane')
    c = r['components'][0]
    assert c['n'] == 2.0 and 'n' not in c['errors'] and abs(c['mag'] - 17.0) < 0.1 and abs(c['re'] / 5 - 1) < 0.1
    # the sky gradient is recovered
    print('MULTIFIT sky plane: gradient %.4f %.4f (truth 0.05 -0.03), sky %.2f (100)' % (r['sky_grad'][0], r['sky_grad'][1], r['sky']))
    assert abs(r['sky_grad'][0] - 0.05) < 0.01 and abs(r['sky_grad'][1] + 0.03) < 0.01 and abs(r['sky'] - 100) < 1.0
    # a constant sky on a gradient image leaves structure (chi2 worse): the plane is needed
    rc = MF.fit(d, [dict(kind='sersic', x=40, y=40, mag=17.2, re=4, n=2.0, q=0.8, pa=0, fixed=['n'])], psf=PSF, rms=1.0, zp=ZP, sky='const')
    assert rc['chi2_red'] > r['chi2_red'] + 0.2
    # bounds: n limited to < 1.2 is hit and flagged
    r = MF.fit(d, [dict(kind='sersic', x=40, y=40, mag=17.2, re=4, n=1.0, q=0.8, pa=0, bounds=dict(n=(0.5, 1.2)))], psf=PSF, rms=1.0, zp=ZP, sky='plane')
    assert r['flags'] & MF.FLAGS['BOUND'] and r['components'][0]['n'] <= 1.2 + 1e-9 and 'n' in r['components'][0]['at_bound']
    # tie: two components share the centre exactly
    r = MF.fit(d, [dict(kind='exp', x=40, y=40, mag=18, re=6, q=0.7, pa=30), dict(kind='dev', x=41, y=39, mag=19, re=2, q=0.9, pa=0)], psf=PSF, rms=1.0, zp=ZP, sky='plane', tie=[('1.x', '0.x'), ('1.y', '0.y')])
    assert r['components'][0]['x'] == r['components'][1]['x'] and r['components'][0]['y'] == r['components'][1]['y']
    # masked pixels are ignored: a bright artifact does not matter when masked
    d2 = d.copy(); d2[20:26, 20:26] += 500
    m = np.zeros(d.shape, bool); m[18:28, 18:28] = True
    r1 = MF.fit(d2, [dict(kind='sersic', x=40, y=40, mag=17.2, re=4, n=1.5, q=0.8, pa=0)], psf=PSF, rms=1.0, zp=ZP, sky='plane', mask=m)
    assert abs(r1['components'][0]['mag'] - 17.0) < 0.1 and r1['npix'] == d.size - m.sum()


def test_spatially_varying_psf_model_is_used():
    """A galaxy in the corner of a field with a strongly varying PSF: the PSF model at that position recovers the size; the central PSF biases it."""
    shape = (300, 300)
    kw = dict(fwhm0=2.4, fwhm1=4.4, e0=0.02, e1=0.02)
    img, t = synth.star_field(shape, n=600, mag_range=(15.5, 20), zp=25, sky=50, noise=0.5, seed=12, **kw)
    mdl, _ = PM.build_psf_model(img.astype(float), fwhm_prior=3.3, snr_min=6, degree=1)
    gx, gy = 270.0, 150.0
    f_here = synth.psf_params(gx, gy, shape, **kw)['fwhm']
    psf_true = M.moffat_psf(31, f_here, beta=3.0, nsub=5)
    gal = MF.render_component(MF._normalise(dict(kind='sersic', x=30, y=30, mag=18.0, re=3.0, n=1.0, q=0.9, pa=0), ZP), (61, 61), psf_true / psf_true.sum())
    rng = np.random.default_rng(3)
    d = gal + 50 + rng.normal(0, 0.5, gal.shape)
    out = {}
    for name, p in (('model at the galaxy', _OffsetModel(mdl, gx, gy)), ('central PSF', mdl.stamp(150, 150))):
        r = MF.fit(d, [dict(kind='sersic', x=30, y=30, mag=18.3, re=2.0, n=1.2, q=0.8, pa=0)], psf=p, rms=0.5, zp=ZP)
        out[name] = r['components'][0]['re'] / 3.0 - 1
    print('MULTIFIT PSF at the field edge (true FWHM %.2f px): re error with the PSF model %+.3f, with the central PSF %+.3f' % (f_here, out['model at the galaxy'], out['central PSF']))
    assert abs(out['model at the galaxy']) < 0.12 and abs(out['central PSF']) > 1.5 * abs(out['model at the galaxy'])


class _OffsetModel:
    """PSFModel evaluated at fixed image position (the cutout lives elsewhere in the image)."""
    def __init__(self, m, x, y):
        self.m, self.x, self.y = m, x, y
        self.size = m.size

    def stamp(self, x, y, dx=0.0, dy=0.0, size=None):
        return self.m.stamp(self.x, self.y, dx, dy, size)


def test_auto_selection_prefers_the_right_model():
    star = image([dict(kind='psf', x=40, y=40, mag=17.0)], noise=1.0, seed=5)
    gal = image([dict(kind='sersic', x=40, y=40, mag=17.0, re=6, n=2.0, q=0.7, pa=30)], noise=1.0, seed=6)
    fl = 10 ** (-0.4 * (17.0 - ZP))
    n1, r1, _ = MF.auto_select(star, 40, 40, fl, 2.0, 0.9, 0.0, psf=PSF, rms=1.0, zp=ZP)
    n2, r2, _ = MF.auto_select(gal, 40, 40, fl, 5.0, 0.7, 30.0, psf=PSF, rms=1.0, zp=ZP)
    print('MULTIFIT auto: star -> %s, galaxy -> %s' % (n1, n2))
    assert n1 == 'psf' and n2 in ('sersic', 'bulge+disk')


def test_gain_poisson_weights_and_errors_scale():
    truth = dict(kind='sersic', x=40, y=40, mag=15.5, re=5, n=2.0, q=0.7, pa=30)
    clean = MF.render_model([MF._normalise(truth, ZP)], (81, 81), PSF, 100.0)
    gain = 4.0
    d = np.random.default_rng(9).poisson(clean * gain) / gain
    r = MF.fit(d.astype(float), [dict(kind='sersic', x=40.4, y=39.7, mag=15.8, re=4, n=1.5, q=0.8, pa=0)], psf=PSF, rms=5.0, gain=gain, zp=ZP)         # rms = sky Poisson noise sqrt(100 / 4)
    c = r['components'][0]
    print('MULTIFIT poisson: dmag %+.3f dre/re %+.3f chi2 %.2f' % (c['mag'] - 15.5, c['re'] / 5 - 1, r['chi2_red']))
    assert abs(c['mag'] - 15.5) < 0.03 and abs(c['re'] / 5 - 1) < 0.05 and 0.8 < r['chi2_red'] < 1.25


# ------------------------------------------------------------------------------------------------------------------ CLI
@pytest.fixture(scope='module')
def field(tmp_path_factory):
    d = tmp_path_factory.mktemp('mf')
    shape = (260, 260)
    rng = np.random.default_rng(21)
    psf = M.moffat_psf(31, 3.0, beta=3.0); psf /= psf.sum()
    gals = [  # x, y (0-based), mag, re, n, q, pa
        (50, 50, 17.0, 5, 1.0, 0.6, 20), (130, 60, 16.5, 4, 3.0, 0.8, 80), (200, 70, 17.5, 6, 1.5, 0.5, 130),
        (60, 160, 16.8, 5, 2.0, 0.7, 40), (72, 170, 17.4, 3, 1.0, 0.9, 0),                       # blend pair, 16 px apart
        (170, 170, 17.2, 4, 4.0, 0.7, 60), (210, 210, 16.9, 5, 2.0, 0.6, 100)]
    comps = [MF._normalise(dict(kind='sersic', x=g[0], y=g[1], mag=g[2], re=g[3], n=g[4], q=g[5], pa=g[6]), ZP) for g in gals]
    img = np.zeros(shape)
    for c in comps:
        h = 45
        x0, y0 = int(c['x']) - h, int(c['y']) - h
        sl = (slice(max(y0, 0), min(y0 + 2 * h + 1, shape[0])), slice(max(x0, 0), min(x0 + 2 * h + 1, shape[1])))
        cc = dict(c, x=c['x'] - sl[1].start, y=c['y'] - sl[0].start)
        img[sl] += MF.render_component(cc, (sl[0].stop - sl[0].start, sl[1].stop - sl[1].start), psf)
    # a star
    st = (110, 200, 17.5)
    cs = dict(kind='psf', x=st[0], y=st[1], flux=10 ** (-0.4 * (st[2] - ZP)))
    img += MF.render_component(cs, shape, psf)
    img += 100.0 + rng.normal(0, 1.0, shape)
    fits.writeto(str(d / 'f.fits'), img.astype(np.float32), overwrite=True)
    rows = []
    for k, g in enumerate(gals):
        A = g[3] * 0.9
        rows.append(dict(NUMBER=k + 1, X_IMAGE=g[0] + 1 + rng.normal(0, 0.5), Y_IMAGE=g[1] + 1 + rng.normal(0, 0.5), FLUX_AUTO=10 ** (-0.4 * (g[2] - ZP)) * 0.85, MAG_AUTO=g[2] + 0.2,
                         A_IMAGE=A, B_IMAGE=A * g[5], THETA_IMAGE=g[6] + rng.normal(0, 10), FLUX_RADIUS=g[3] * 0.8, CLASS_STAR=0.05, KRON_RADIUS=3.0))
    rows.append(dict(NUMBER=len(gals) + 1, X_IMAGE=st[0] + 1.2, Y_IMAGE=st[1] + 0.8, FLUX_AUTO=10 ** (-0.4 * (st[2] - ZP)), MAG_AUTO=st[2], A_IMAGE=1.3, B_IMAGE=1.3, THETA_IMAGE=0.0,
                     FLUX_RADIUS=1.8, CLASS_STAR=0.98, KRON_RADIUS=3.0))
    cols = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'MAG_AUTO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'FLUX_RADIUS', 'CLASS_STAR', 'KRON_RADIUS']
    tsvio.write_table(str(d / 'cat.tsv'), cols, rows)
    fits.writeto(str(d / 'psf.fits'), psf.astype(np.float32), overwrite=True)
    return dict(dir=str(d), image=str(d / 'f.fits'), cat=str(d / 'cat.tsv'), psf=str(d / 'psf.fits'), gals=gals, star=st)


def run_cli(field, name, *extra, cat=True):
    wd = os.path.join(field['dir'], name)
    cmd = [PY, os.path.join(PLUGIN, 'multifit.py'), field['image'], '--work', wd, '--psf', field['psf'], '--mag-zeropoint', '25'] + (['--catalog', field['cat']] if cat else []) + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    return wd, r


def parse(stdout):
    L = stdout.strip().splitlines()
    head = L[0].split('\t')
    return [dict(zip(head, l.split('\t'))) for l in L[1:]]


def test_cli_catalog_fit_neighbours_vs_ignore(field):
    out = {}
    for mode in ('fit', 'ignore'):
        wd, r = run_cli(field, mode, '--neighbours', mode, '--model', 'sersic', '--max-objects', '20', '--n-workers', '4')
        rows = parse(r.stdout)
        assert len(rows) == len(field['gals']) + 1
        out[mode] = (wd, rows)
    wd, rows = out['fit']
    cols = rows[0].keys()
    assert 'NUMBER' in cols and all(c in cols for c in ('GF_MAG', 'GF_RE', 'GF_N', 'GF_Q', 'GF_PA', 'GF_CHI2', 'GF_FLAG', 'GF_NNEIGH', 'GF_X'))
    dm, dre, chis = [], [], []
    for k, g in enumerate(field['gals']):
        r = rows[k]
        dm.append(float(r['GF_MAG']) - g[2]); dre.append(float(r['GF_RE']) / g[3] - 1); chis.append(float(r['GF_CHI2']))
    print('MULTIFIT catalog fit (7 galaxies incl. a blend pair): median dmag %+.3f, max |dmag| %.3f, median dre/re %+.3f, max |dre/re| %.3f, chi2 median %.2f max %.2f' % (
        np.median(dm), np.max(np.abs(dm)), np.median(dre), np.max(np.abs(dre)), np.median(chis), np.max(chis)))
    assert np.max(np.abs(dm)) < 0.12 and np.max(np.abs(dre)) < 0.15 and np.median(chis) < 1.3
    # blend pair: objects 4 and 5 are 16 px apart; fitting the neighbour simultaneously must be at least as good as ignoring it
    ig = out['ignore'][1]
    e_fit = abs(float(rows[3]['GF_MAG']) - field['gals'][3][2]) + abs(float(rows[4]['GF_MAG']) - field['gals'][4][2])
    e_ign = abs(float(ig[3]['GF_MAG']) - field['gals'][3][2]) + abs(float(ig[4]['GF_MAG']) - field['gals'][4][2])
    print('MULTIFIT blend pair |dmag| sum: simultaneous %.3f, ignoring the neighbour %.3f' % (e_fit, e_ign))
    assert e_fit < e_ign and int(rows[3]['GF_NNEIGH']) >= 1
    # files
    for f in ('results.tsv', 'model.fits', 'residual.fits', 'montage.png', 'psf.json'):
        assert os.path.getsize(os.path.join(wd, 'multifit_' + f)) > 100, f
    res = fits.getdata(os.path.join(wd, 'multifit_residual.fits'))
    data = fits.getdata(field['image'])
    # residual image: the galaxies are removed (rms inside the fitted objects' areas back to the noise level)
    gx, gy = int(field['gals'][1][0]), int(field['gals'][1][1])
    core = (slice(gy - 20, gy + 21), slice(gx - 20, gx + 21))
    print('MULTIFIT residual rms in the n=3 galaxy core region: data %.2f -> residual %.2f (noise 1.0)' % (np.std(data[core]), np.std(res[core])))
    assert np.std(res[core]) < 1.3 and np.std(data[core]) > 3
    cols_, recs = tsvio.read_catalog(os.path.join(wd, 'multifit_results.tsv'))
    assert len(recs) == len(field['gals']) + 1 and 'KIND' in cols_


def test_cli_models_psf_star_and_bulge_disk(field):
    wd, r = run_cli(field, 'star', '--model', 'psf', '--objects', str(len(field['gals']) + 1), '--neighbours', 'mask')
    rows = parse(r.stdout)
    star = rows[-1]
    print('MULTIFIT CLI star (psf model): dmag %+.3f chi2 %.2f dpos (%.3f, %.3f)' % (float(star['GF_MAG']) - field['star'][2], float(star['GF_CHI2']), float(star['GF_X']) - 1 - field['star'][0], float(star['GF_Y']) - 1 - field['star'][1]))
    assert abs(float(star['GF_MAG']) - field['star'][2]) < 0.06 and float(star['GF_CHI2']) < 1.5
    assert all(rw['GF_MAG'] == '' for rw in rows[:-1])        # only the requested object is fitted
    wd, r = run_cli(field, 'auto', '--model', 'auto', '--objects', '2;%d' % (len(field['gals']) + 1), '--neighbours', 'mask')
    rows = parse(r.stdout)
    print('MULTIFIT CLI auto: galaxy ncomp %s chi2 %s; star ncomp %s' % (rows[1]['GF_NCOMP'], rows[1]['GF_CHI2'], rows[-1]['GF_NCOMP']))
    assert rows[-1]['GF_NCOMP'] == '1' and rows[-1]['GF_RE'] == ''           # the star is a PSF component (no size)
    wd, r = run_cli(field, 'bd', '--model', 'bulge+disk', '--objects', '2', '--neighbours', 'mask')
    rows = parse(r.stdout)
    assert rows[1]['GF_BT'] != '' and 0 <= float(rows[1]['GF_BT']) <= 1 and rows[1]['GF_NCOMP'] == '2'


def test_cli_config_mode(field, tmp_path):
    g = field['gals'][1]
    cfg = dict(bbox=[g[0] - 30, g[0] + 31, g[1] - 30, g[1] + 31], sky='const', tie=[], components=[
        dict(kind='sersic', x=g[0] + 1.5, y=g[1] + 0.5, mag=g[2] + 0.3, re=3.0, n=2.0, q=0.7, pa=60, fixed=['n'], bounds=dict(q=[0.3, 1.0]))])
    p = str(tmp_path / 'cfg.json')
    json.dump(cfg, open(p, 'w'))
    wd = str(tmp_path / 'w')
    r = subprocess.run([PY, os.path.join(PLUGIN, 'multifit.py'), field['image'], '--work', wd, '--psf', field['psf'], '--config', p, '--mag-zeropoint', '25'], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-1500:]
    assert 'chi2/dof' in r.stdout and 'sersic' in r.stdout
    cols, rows = tsvio.read_catalog(os.path.join(wd, 'multifit_results.tsv'))
    assert abs(float(rows[0]['MAG']) - g[2]) < 0.2 and float(rows[0]['N']) == 2.0     # n fixed at 2 (truth 3): mag still close
    assert os.path.getsize(os.path.join(wd, 'multifit_residual.fits')) > 1000
