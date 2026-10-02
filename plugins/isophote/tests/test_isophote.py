import json, os, subprocess, sys
import numpy as np
import pytest
import isophote as I
from ogfkit import models, imageio, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')
SHAPE = (301, 301)
SIG = 0.05


def synth(n=2.0, e=0.3, pa=40.0, re=15.0, Ie=10.0, c=0.0, seed=1, xc=150.0, yc=150.0, shape=SHAPE):
    img = models.render_sersic(shape, xc, yc, Ie, re, n, 1.0 - e, pa, c=c)
    return img + np.random.default_rng(seed).normal(0, SIG, shape), img


def profile_at(res, r):
    t = res['table']
    s = np.array(t['sma'])
    i = int(np.argmin(np.abs(s - r)))
    return {k: t[k][i] for k in t}


@pytest.mark.parametrize('n', [1.0, 4.0])
def test_sersic_recovers_e_pa_intensity(n):
    e_true, pa_true, re, Ie = 0.3, 40.0, 15.0, 10.0
    img, clean = synth(n=n, e=e_true, pa=pa_true, re=re, Ie=Ie)
    res, model = I.fit_galaxy(img, 150, 150, eps0=0.15, pa0=25, maxsma=75, step=0.12, zp=25.0)
    assert res['status'] == 'ok'
    t = res['table']
    s = np.array(t['sma'])
    sel = (s > 0.4 * re) & (s < 3.0 * re)
    eps = np.array(t['ellipticity'])[sel]
    pa = np.array(t['pa'])[sel]
    assert abs(np.median(eps) - e_true) < 0.01, np.median(eps)
    assert np.max(np.abs(eps - e_true)) < 0.03
    assert abs(((np.median(pa) - pa_true + 90) % 180) - 90) < 1.0, np.median(pa)
    # intensity vs the analytic profile I(a) at the semi-major axis
    bn = models.sersic_bn(n)
    for r in (0.3 * re, 1.0 * re, 2.0 * re):
        p = profile_at(res, r)
        ana = Ie * np.exp(-bn * ((p['sma'] / re) ** (1 / n) - 1))
        # mean intensity on an ellipse of fixed e: equals the profile at the same semi-major axis
        assert abs(p['intens'] / ana - 1) < 0.03, (r, p['intens'], ana)
    # centre
    assert abs(np.median(np.array(t['x0'])[sel]) - 150) < 0.1 and abs(np.median(np.array(t['y0'])[sel]) - 150) < 0.1
    # model and residual
    resid = img - model
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    rr = np.hypot(xx - 150, yy - 150)
    ann = (rr > 4) & (rr < 60)          # the n=4 cusp (r < 4 px) is not reproduced by a sampled isophote model: stated limitation
    assert abs(np.mean(resid[ann])) < 0.02
    assert np.std(resid[ann]) < 3 * SIG * 1.2


def test_growth_curve_flux_and_radii():
    n, e, re, Ie = 1.0, 0.3, 15.0, 10.0
    img, clean = synth(n=n, e=e, re=re, Ie=Ie)
    res, _ = I.fit_galaxy(img, 150, 150, eps0=0.2, pa0=40, maxsma=75, step=0.12, zp=25.0)
    t = res['table']
    p = profile_at(res, 3 * re)
    # analytic flux inside the isophote of semi-major axis a for the exponential profile (axis ratio q)
    q = 1 - e
    from scipy.special import gammainc
    bn = models.sersic_bn(n)
    Ftot = models.sersic_flux_total(Ie, re, n, q)
    a = p['sma']
    ana = Ftot * gammainc(2 * n, bn * (a / re) ** (1 / n))
    assert abs(p['growth_flux'] / ana - 1) < 0.02, (p['growth_flux'], ana)
    assert abs(p['tflux_e'] / ana - 1) < 0.02, (p['tflux_e'], ana)       # photutils' own sum agrees when nothing is masked
    # the half-light radius is relative to the flux inside the last isophote (5 re): r50 within 5 % of the truncated-curve value
    smas = np.linspace(1, 75, 400)
    f = Ftot * gammainc(2 * n, bn * (smas / re) ** (1 / n))
    r50_trunc = np.interp(0.5 * np.interp(res['summary']['rmax'], smas, f), f, smas)
    assert abs(res['summary']['r50'] / r50_trunc - 1) < 0.05, (res['summary']['r50'], r50_trunc)
    assert abs(res['summary']['mag_total'] - (25 - 2.5 * np.log10(np.interp(res['summary']['rmax'], smas, f)))) < 0.03


def test_boxy_disky_sign():
    out = {}
    for c, name in ((0.5, 'boxy'), (-0.4, 'disky'), (0.0, 'elliptical')):
        img, _ = synth(n=2.0, c=c)
        res, _ = I.fit_galaxy(img, 150, 150, eps0=0.2, pa0=40, maxsma=60, step=0.15, with_model=False)
        out[name] = res['summary']
    assert out['boxy']['shape'] == 'boxy' and out['boxy']['b4_hl'] < -I.DISKY_LIMIT
    assert out['disky']['shape'] == 'disky' and out['disky']['b4_hl'] > I.DISKY_LIMIT
    assert out['elliptical']['shape'] == 'elliptical' and abs(out['elliptical']['b4_hl']) < I.DISKY_LIMIT


def test_mask_aware():
    img, clean = synth(n=1.0, e=0.3, pa=40, re=15)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    bx, by = 150 + 38 * np.cos(np.deg2rad(40)), 150 + 38 * np.sin(np.deg2rad(40))
    blob = 1.5 * np.exp(-0.5 * (((xx - bx) ** 2 + (yy - by) ** 2) / 5.0 ** 2))      # extended companion on the major axis
    bad = img + blob
    mask = ((xx - bx) ** 2 + (yy - by) ** 2) < 16 ** 2
    kw = dict(eps0=0.2, pa0=40, maxsma=60, step=0.15, with_model=False)
    r_nomask, _ = I.fit_galaxy(bad, 150, 150, **kw)
    r_mask, _ = I.fit_galaxy(bad, 150, 150, mask=mask, **kw)

    def dev(r):
        t = r['table']; s = np.array(t['sma']); sel = (s > 25) & (s < 50)
        return np.max(np.abs(np.array(t['ellipticity'])[sel] - 0.3)), np.max(np.abs(np.array(t['intens'])[sel] / (10 * np.exp(-models.sersic_bn(1.0) * (s[sel] / 15 - 1))) - 1))
    (e_m, i_m), (e_n, i_n) = dev(r_mask), dev(r_nomask)
    print('masked  max|d eps| %.4f max|d I/I| %.3f ; unmasked %.4f %.3f' % (e_m, i_m, e_n, i_n))
    assert e_m < 0.02 and i_m < 0.05
    # mask-aware growth curve: with the companion masked the flux inside 2 re is still the analytic one (tflux_e would lose the masked pixels)
    q = 0.7
    bn = models.sersic_bn(1.0)
    from scipy.special import gammainc
    Ftot = models.sersic_flux_total(10.0, 15.0, 1.0, q)
    pm = profile_at(r_mask, 30.0)
    ana = Ftot * gammainc(2.0, bn * pm['sma'] / 15.0)
    assert abs(pm['growth_flux'] / ana - 1) < 0.03, (pm['growth_flux'], ana)
    assert pm['tflux_e'] < 0.99 * ana          # photutils' sum is biased low by the masked pixels, documented
    assert i_n > 0.10 or e_n > 0.04, (e_n, i_n)       # the companion does corrupt the unmasked fit


def write_cat(path, rows):
    cols = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'KRON_RADIUS', 'MAG_AUTO']
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(v) for v in r) + '\n')


def test_cli_two_galaxies_mask_and_outputs(tmp_path):
    shape = (301, 401)
    g1 = models.render_sersic(shape, 100, 150, 10, 12, 1.0, 0.7, 40)
    g2 = models.render_sersic(shape, 290, 150, 8, 10, 2.0, 0.5, 120)
    img = g1 + g2 + np.random.default_rng(2).normal(0, SIG, shape)
    imageio.save_fits(str(tmp_path / 'im.fits'), img)
    write_cat(str(tmp_path / 'cat.tsv'), [(1, 101, 151, 14, 10, 40, 3, 15.0), (2, 291, 151, 11, 6, 120, 3, 16.0)])
    outs = {k: str(tmp_path / k) for k in ('m.fits', 'r.fits', 't.tsv', 'j.json', 'p.png')}
    cmd = [PY, os.path.join(os.path.dirname(HERE), 'isophote.py'), str(tmp_path / 'im.fits'), '--catalog', str(tmp_path / 'cat.tsv'),
           '--mag-zeropoint', '25', '--pixel-scale', '0.1', '--maxsma-scale', '2.5', '--step', '0.15',
           '--model-out', outs['m.fits'], '--resid-out', outs['r.fits'], '--table-out', outs['t.tsv'],
           '--json-out', outs['j.json'], '--plot-out', outs['p.png']]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    lines = p.stdout.strip().split('\n')
    cols = lines[0].split('\t')
    assert cols == ['NUMBER'] + I.COLUMNS
    rows = {int(l.split('\t')[0]): dict(zip(cols, l.split('\t'))) for l in lines[1:]}
    assert set(rows) == {1, 2}
    assert abs(float(rows[1]['ISO_EPS_HL']) - 0.3) < 0.02 and abs(float(rows[2]['ISO_EPS_HL']) - 0.5) < 0.02
    assert abs(float(rows[1]['ISO_PA_HL']) - 40) < 1.5 and abs(float(rows[2]['ISO_PA_HL']) - 120) < 1.5
    assert float(rows[1]['ISO_FLAG']) == 0
    m, _ = imageio.load_image(outs['m.fits'])
    r, _ = imageio.load_image(outs['r.fits'])
    assert m.shape == shape and np.allclose(m + r, img, atol=1e-4)
    assert np.std(r[100:200, 50:150]) < 0.08
    assert os.path.getsize(outs['p.png']) > 5000
    tab = tsvio.read_catalog(outs['t.tsv'])[1]
    assert {r_['NUMBER'] for r_ in tab} == {'1', '2'}
    # a bit-flag mask file is understood (bit 2 = star), the target itself is never masked
    flags = np.zeros(shape, np.uint8)
    flags[140:160, 90:110] = 1       # "detected source" bit over the target 1: must not break the fit
    flags[100:110, 120:130] = 2
    imageio.save_fits(str(tmp_path / 'mask.fits'), flags)
    p2 = subprocess.run(cmd[:-10] + ['--mask', str(tmp_path / 'mask.fits'), '--numbers', '1'], capture_output=True, text=True)
    assert p2.returncode == 0, p2.stderr
    l2 = p2.stdout.strip().split('\n')
    assert len(l2) == 2 and 'masked' in p2.stderr


def test_effective_mask_rule():
    f = np.array([0, 1, 2, 4, 8, 16, 9, 12], np.uint8)
    assert list(imageio.effective_mask(f)) == [False, True, True, True, False, True, False, False]


@pytest.mark.skipif(not os.path.exists(FITS + '/m51.fits'), reason='m51 not available')
def test_m51_smoke():
    data, _ = imageio.load_image(FITS + '/m51.fits')
    iy, ix = np.unravel_index(np.argmax(data[250:350, 250:350]), (100, 100))
    res, model = I.fit_galaxy(data, 250 + ix, 250 + iy, eps0=0.2, pa0=0, maxsma=120, step=0.2, sma0=6, with_model=True)
    assert res['status'] in ('ok', 'few')
    assert len(res['table']['sma']) >= 5
    json.dumps(res)         # server-ready: JSON serialisable
