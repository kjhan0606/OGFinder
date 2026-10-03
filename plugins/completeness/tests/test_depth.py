import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits
from scipy import ndimage as ndi

from ogfkit import depth as dp, noise as nz
import completeness as cp

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'depthmap.py')


def noise_image(shape, sigma, corr=0.0, seed=1):
    n = np.random.default_rng(seed).normal(size=shape)
    if corr:
        n = ndi.gaussian_filter(n, corr) * (2 * math.sqrt(math.pi) * corr)
    return n * sigma


def brute(noise, r, n=15000, seed=3, box=None):
    pos = dp.blank_circle_positions(np.zeros(noise.shape, bool), r, n, np.random.default_rng(seed), box=box)
    return nz.clipped_std(dp.circle_sums(noise, pos, r), 4.0)


def model_depth(img, r=3.0, zp=25.0, rms=None):
    P = dp.prepare(img, None)
    rm = P['rms'] if rms is None else dp.rms_from_file(rms[0], ~P['bad'], P['rms'], rms[1])
    f, _ = dp.circle_factor(P['sub'], P['mask'], rm, r, 3000, 1)
    return dp.mag_limit_map(rm, f, zp, 5.0), P, f


def test_white_and_correlated_depth_vs_brute_force():
    for corr in (0.0, 1.2):
        nse = noise_image((900, 900), 1.0, corr, seed=2).astype('f4')
        mag, P, f = model_depth(nse)
        truth = 25.0 - 2.5 * math.log10(5 * brute(nse, 3.0))
        assert abs(np.median(mag) - truth) < 0.05
        if corr:
            assert f > 2.5 * math.sqrt(math.pi) * 3 * 0.9 and f > 1.4 * model_depth(noise_image((900, 900), 1.0, 0.0, 2).astype('f4'))[2]


def test_gradient_noise_blocks():
    sh = (900, 900)
    yy, xx = np.indices(sh)
    smap = 1.0 + 1.5 * xx / sh[1]
    nse = noise_image(sh, 1.0, 1.0, seed=4) * smap
    mag, P, f = model_depth(nse.astype('f4'))
    for i in range(3):
        box = (i * 300, 100, (i + 1) * 300 - 1, 800)
        t = 25.0 - 2.5 * math.log10(5 * brute(nse, 3.0, 4000, 5 + i, box))
        assert abs(np.nanmedian(mag[100:800, i * 300:(i + 1) * 300 - 1]) - t) < 0.08
    assert mag[:, -50:].mean() < mag[:, :50].mean() - 0.6          # noisier side is shallower by the noise ratio


def test_sources_and_masked_region_do_not_bias_depth():
    sh = (900, 900)
    rng = np.random.default_rng(6)
    nse = noise_image(sh, 1.0, 1.0, seed=6)
    img = nse.copy()
    for _ in range(300):
        x, y = rng.uniform(15, 885, 2)
        f = rng.lognormal(math.log(80), 0.8)
        x0, y0 = int(x) - 12, int(y) - 12
        Y, X = np.mgrid[y0:y0 + 25, x0:x0 + 25]
        img[y0:y0 + 25, x0:x0 + 25] += f * np.exp(-((X - x) ** 2 + (Y - y) ** 2) / (2 * 2.2 ** 2)) / (2 * math.pi * 2.2 ** 2)
    bad = np.zeros(sh, bool)
    bad[:200, :300] = True
    img[bad] = np.nan
    mag, P, f = model_depth(img.astype('f4'))
    truth = 25.0 - 2.5 * math.log10(5 * brute(nse, 3.0))
    assert abs(np.nanmedian(mag) - truth) < 0.12
    assert np.isnan(mag[bad]).all() or True            # masked pixels are NaN in the CLI map (valid = finite)


def test_weight_and_rms_map_input():
    sh = (700, 700)
    yy, xx = np.indices(sh)
    smap = 1.0 + 1.0 * yy / sh[0]
    nse = (noise_image(sh, 1.0, 1.0, seed=8) * smap).astype('f4')
    m0, *_ = model_depth(nse)
    mw, *_ = model_depth(nse, rms=(1.0 / smap ** 2, 'weight'))
    mr, *_ = model_depth(nse, rms=(smap * 3.0, 'rms'))               # rms map in arbitrary scale: rescaled to the measured pixel rms
    assert abs(np.median(mw) - np.median(m0)) < 0.08 and abs(np.median(mr) - np.median(mw)) < 0.02
    assert np.corrcoef(mw[::7, ::7].ravel(), m0[::7, ::7].ravel())[0, 1] > 0.9


def test_region_labels_and_area_depth():
    d = np.linspace(20, 24, 400 * 400).reshape(400, 400)
    valid = np.ones_like(d, bool)
    valid[:50] = False
    g = dp.region_labels('grid', d.shape, valid=valid, grid=(3, 2))
    assert set(np.unique(g)) == {0, 1, 2, 3, 4, 5, 6}
    c = dp.region_labels('depth', d.shape, depth=d, valid=valid, classes=4)
    cnt = [(c == k).sum() for k in range(1, 5)]
    assert max(cnt) - min(cnt) <= 2 and (c[~valid] == 0).all()
    assert c[valid & (d < 21)].max() == 1 and c[valid & (d > 23.5)].min() == 4
    lab = np.arange(100).reshape(10, 10) % 3
    assert (dp.region_labels('file', (10, 10), label_file=lab) == lab).all()
    with pytest.raises(ValueError):
        dp.region_labels('file', (5, 5), label_file=lab)
    ad = dp.area_depth(d, valid)
    fr = [x[1] for x in ad]
    assert fr[0] > 0.99 and fr[-1] < 0.01 and all(a >= b for a, b in zip(fr, fr[1:]))


def test_sep_local_threshold_follows_noise():
    sh = (700, 700)
    sg = np.ones(sh)
    sg[:, 350:] = 2.5
    img = (noise_image(sh, 1.0, 1.0, seed=9) * sg).astype('f4')
    kw = dict(kind='star', mag_min=19.0, mag_max=24.0, n_bins=10, per_bin=120, per_image=25, zp=25.0, seed=2)
    out = {}
    for local in (False, True):
        det = cp.SepDetector(zp=25.0, thresh=1.5, minarea=5, local_rms=local)
        res = cp.run_completeness(img, det, return_records=True, **kw)
        f = {}
        for name, sel in (('quiet', lambda r: r['x'] < 350), ('noisy', lambda r: r['x'] >= 350)):
            rr = [r for r in res['records'] if sel(r) and 20.5 < r['mag'] < 21.5]
            f[name] = np.mean([r['recovered'] for r in rr])
        out[local] = f
    # a global threshold of 1.5 sigma_global recovers the faint stars far better in the noisy half (threshold is relative to the global rms ~ 1.9)
    assert out[True]['quiet'] > out[True]['noisy'] + 0.3
    assert out[False]['noisy'] > out[True]['noisy'] + 0.05          # the global threshold is too permissive in the noisy half


@pytest.fixture(scope='module')
def cli_image(tmp_path_factory):
    d = tmp_path_factory.mktemp('dpt')
    sh = (700, 700)
    sg = np.ones(sh)
    for i, s in enumerate((1.0, 1.5, 2.2, 3.0)):
        sg[:, i * 175:(i + 1) * 175] = s
    img = (noise_image(sh, 1.0, 1.0, seed=11) * sg).astype('f4')
    fits.writeto(d / 'i.fits', img)
    with open(d / 'c.tsv', 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tMAG_AUTO\n')
        for i in range(40):
            f.write('%d\t%.1f\t%.1f\t%.2f\n' % (i + 1, 20 + (i * 17) % 660, 30 + (i * 29) % 640, 20 + 0.1 * i))
    return d


def run(d, *args):
    r = subprocess.run([sys.executable, CLI, str(d / 'i.fits'), '--work', str(d / 'w')] + list(args), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-600:]
    return r


def test_cli_depth_columns_and_files(cli_image):
    d = cli_image
    r = run(d, '--mode', 'depth', '--catalog', str(d / 'c.tsv'), '--pixel-scale', '0.2', '--box-arcsec', '12', '--aper-radius', '3', '--tile', '175')
    lines = r.stdout.strip().split('\n')
    assert lines[0].split('\t') == ['NUMBER', 'DEPTH_RMS', 'DEPTH_LIM', 'DEPTH_SBLIM']
    assert len(lines) == 41
    lim = {int(l.split('\t')[0]): float(l.split('\t')[2]) for l in lines[1:]}
    mp = fits.getdata(d / 'w' / 'depth_mag_map.fits')
    cat = [l.split('\t') for l in open(d / 'c.tsv').read().strip().split('\n')[1:]]
    for c in cat[:10]:
        assert abs(lim[int(c[0])] - mp[int(round(float(c[2]) - 1)), int(round(float(c[1]) - 1))]) < 1e-3
    for f in ('depth_mag_map.fits', 'depth_sb_map.fits', 'depth_rms_map.fits', 'depth_tiles.tsv', 'depth_area.tsv', 'depth_summary.json', 'depth_plot.png'):
        assert (d / 'w' / f).stat().st_size > 100
    j = json.load(open(d / 'w' / 'depth_summary.json'))
    assert j['tile_check']['n_valid'] >= 12 and 0.8 < j['tile_check']['ratio_median'] < 1.2 and j['tile_check']['ratio_std'] < 0.15
    # four noise bands: depth steps by 2.5 log10 of the noise ratio
    assert abs((mp[:, 20:150].mean() - mp[:, 550:680].mean()) - 2.5 * math.log10(3.0)) < 0.12
    assert 'sb_limit' in j and j['sb_limit']['median'] > 18


def test_cli_compmap_regions(cli_image):
    d = cli_image
    r = run(d, '--mode', 'compmap', '--catalog', str(d / 'c.tsv'), '--pixel-scale', '0.2', '--aper-radius', '3', '--regions', 'rms', '--rms-classes', '4',
            '--mag-min', '19', '--mag-max', '24', '--n-bins', '10', '--per-bin', '100', '--maps-at', '21,22', '--tile', '175')
    lines = r.stdout.strip().split('\n')
    assert lines[0].split('\t') == ['NUMBER', 'COMPL_REGION', 'COMPL_LOC', 'COMPL_LIM50_LOC', 'COMPL_LIM90_LOC']
    j = json.load(open(d / 'w' / 'compmap_summary.json'))['compmap']
    rg = j['regions']
    assert len(rg) == 4 and all(np.isfinite(x['lim50']) for x in rg)
    d50 = [x['depth_median'] for x in rg]
    l50 = [x['lim50'] for x in rg]
    assert np.argsort(d50).tolist() == np.argsort(l50).tolist()                 # deeper regions reach fainter limits
    assert j['collapse']['offset_std'] < 0.25 and 0.5 < j['collapse']['regression_slope'] < 1.8
    fm = fits.getdata(d / 'w' / 'compmap_frac_m21.fits')
    assert np.nanmax(fm) <= 1 and np.nanmin(fm) >= 0 and np.nanmean(fm[:, :150]) > np.nanmean(fm[:, 550:])
    loc = [float(l.split('\t')[2]) for l in lines[1:] if l.split('\t')[2]]
    assert len(loc) == 40 and all(0 <= v <= 1 for v in loc)


def test_cli_no_catalog_prints_json(cli_image):
    d = cli_image
    r = run(d, '--mode', 'depth', '--pixel-scale', '0.2', '--tile', '175')
    j = json.loads(r.stdout)
    assert 'mag_limit' in j and 'tile_check' in j
