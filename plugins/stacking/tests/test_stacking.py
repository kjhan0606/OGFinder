import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import stacking as st

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'stacking.py')


def stamp(img, x, y, F, sg):
    r = int(math.ceil(6 * sg))
    x0, y0 = int(round(x)) - r, int(round(y)) - r
    yy, xx = np.mgrid[y0:y0 + 2 * r + 1, x0:x0 + 2 * r + 1]
    img[y0:y0 + 2 * r + 1, x0:x0 + 2 * r + 1] += F * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sg * sg)) / (2 * math.pi * sg * sg)


def field(seed=1, n=300, F=3.0, sg=2.0, sh=(900, 900), noise=1.0):
    rng = np.random.default_rng(seed)
    img = rng.normal(0, noise, sh)
    xs = rng.uniform(30, sh[1] - 30, n)
    ys = rng.uniform(30, sh[0] - 30, n)
    for x, y in zip(xs, ys):
        stamp(img, x, y, F, sg)
    return img, xs, ys


def test_subpixel_alignment_centroid():
    rng = np.random.default_rng(2)
    img = np.zeros((300, 300))
    gx, gy = np.meshgrid(np.arange(5) * 50 + 40, np.arange(4) * 60 + 40)
    xs = gx.ravel() + rng.uniform(-0.5, 0.5, 20)
    ys = gy.ravel() + rng.uniform(-0.5, 0.5, 20)
    for x, y in zip(xs, ys):
        stamp(img, x, y, 1000.0, 1.5)
    cfg = st.StackConfig(half=10, bkg='none')
    c = st.make_cutouts(img, None, xs, ys, cfg)
    assert c['ok'].all()
    yy, xx = np.indices((21, 21)) - 10
    for cu in c['cube'][:20]:
        cx = (cu * xx).sum() / cu.sum()
        cy = (cu * yy).sum() / cu.sum()
        assert abs(cx) < 0.03 and abs(cy) < 0.03
    ci = st.make_cutouts(img, None, xs, ys, st.StackConfig(half=10, bkg='none', align='int'))
    off = [abs((cu * xx).sum() / cu.sum()) for cu in ci['cube']]
    assert max(off) > 0.1          # the integer alignment keeps the sub-pixel offset


def test_combine_nan_and_methods():
    rng = np.random.default_rng(0)
    cube = rng.normal(5, 1, (50, 9, 9)).astype('f4')
    cube[:10, 0, 0] = np.nan
    cube[0, 4, 4] = 1000.0
    m = st.combine(cube, 'mean')
    med = st.combine(cube, 'median')
    cl = st.combine(cube, 'clipmean')
    assert np.isfinite(m).all()
    assert m[4, 4] > 20 and abs(med[4, 4] - 5) < 0.5 and abs(cl[4, 4] - 5) < 0.5
    w = np.ones(50)
    w[0] = 0.0
    assert abs(st.combine(cube, 'wmean', w)[4, 4] - 5) < 0.5
    cube[:, 1, 1] = np.nan
    assert np.isnan(st.combine(cube, 'mean')[1, 1])


def test_faint_stack_unbiased_and_bootstrap_calibrated():
    img, xs, ys = field(seed=3, n=400, F=3.0)
    cfg = st.StackConfig(half=15, ap_r=8.0, n_boot=120, seed=4)
    r = st.run_stack(img, None, xs, ys, cfg)
    # analytic error: sigma sqrt(Npix) / sqrt(N)  (+ background term)
    an = math.sqrt(math.pi * 8.0 ** 2) / math.sqrt(400)
    assert abs(r['aper'] - 3.0) < 4 * r['aper_err']
    assert 0.7 * an < r['aper_err'] < 1.5 * an
    # S/N scales like sqrt(N)
    r4 = st.run_stack(img, None, xs[:100], ys[:100], cfg)
    assert 1.4 < r4['aper_err'] / r['aper_err'] < 2.8


def test_bootstrap_error_calibration_many_realisations():
    pulls = []
    for k in range(25):
        img, xs, ys = field(seed=50 + k, n=200, F=3.0, sh=(600, 600))
        r = st.run_stack(img, None, xs, ys, st.StackConfig(half=14, ap_r=8.0, n_boot=60, seed=k))
        pulls.append((r['aper'] - 3.0) / r['aper_err'])
    pulls = np.array(pulls)
    assert abs(pulls.mean()) < 0.6 and 0.6 < pulls.std(ddof=1) < 1.5


def test_mask_junk_ignored_and_rejection():
    img, xs, ys = field(seed=5, n=200, F=6.0, sh=(700, 700))
    bad = np.zeros(img.shape, bool)
    bad[:, 300:400] = True
    img[bad] = 1e7
    cfg = st.StackConfig(half=14, ap_r=8.0, n_boot=30, min_valid=0.3, core_r=0.0)
    r = st.run_stack(img, bad, xs, ys, cfg)
    assert abs(r['aper'] - 6.0) < 4 * r['aper_err']
    # a source whose core is masked is rejected with the default core check
    cfg2 = st.StackConfig(half=14, n_boot=1)
    sel = np.where((xs > 305) & (xs < 395))[0]
    c = st.make_cutouts(img, bad, xs[sel], ys[sel], cfg2)
    assert not c['ok'].any() and all(x for x in c['reason'])
    # outside the image
    c = st.make_cutouts(img, bad, np.array([-50.0, 3.0]), np.array([10.0, 3.0]), cfg2)
    assert not c['ok'].any()


def test_flux_normalisation():
    rng = np.random.default_rng(7)
    img = np.zeros((500, 500))
    gx, gy = np.meshgrid(np.arange(8) * 55 + 40, np.arange(5) * 90 + 40)
    xs = gx.ravel() + rng.uniform(-0.5, 0.5, 40)
    ys = gy.ravel() + rng.uniform(-0.5, 0.5, 40)
    F = rng.uniform(100, 5000, 40)
    for x, y, f in zip(xs, ys, F):
        stamp(img, x, y, f, 1.5)
    cfg = st.StackConfig(half=12, bkg='none', norm='flux', norm_r=8.0, n_boot=10, ap_r=8.0)
    r = st.run_stack(img, None, xs, ys, cfg)
    assert abs(r['aper'] - 1.0) < 0.01 and r['aper_err'] < 0.02
    assert np.allclose(r['cutouts']['norm'] * F, 1.0, atol=0.01)


def test_rotation_and_scale_alignment():
    rng = np.random.default_rng(11)
    img = np.zeros((900, 900))
    n = 60
    gx, gy = np.meshgrid(np.arange(10) * 80 + 60, np.arange(6) * 140 + 60)
    xs = gx.ravel() + rng.uniform(-1, 1, n)
    ys = gy.ravel() + rng.uniform(-1, 1, n)
    th = rng.uniform(0, 180, n)
    size = rng.uniform(3.0, 6.0, n)          # sigma along the major axis
    q = 0.4
    yy, xx = np.mgrid[0:900, 0:900]
    for x, y, t, s in zip(xs, ys, th, size):
        x0, y0 = int(x) - 30, int(y) - 30
        Y, X = np.mgrid[y0:y0 + 61, x0:x0 + 61]
        c, sn = math.cos(math.radians(t)), math.sin(math.radians(t))
        u = (X - x) * c + (Y - y) * sn
        v = -(X - x) * sn + (Y - y) * c
        img[y0:y0 + 61, x0:x0 + 61] += 100 * np.exp(-0.5 * ((u / s) ** 2 + (v / (q * s)) ** 2)) / (s * s)
    # rotated (major axis -> +x) and rescaled to sigma = 4 output px
    cfg = st.StackConfig(half=12, bkg='none', n_boot=1, ap_r=10.0)
    c = st.make_cutouts(img, None, xs, ys, cfg, scales=size / 4.0, angles=th)
    assert c['ok'].all()
    s = st.combine(c['cube'], 'mean')
    y, x = np.indices(s.shape) - 12
    sx = math.sqrt((s * x * x).sum() / s.sum())
    sy = math.sqrt((s * y * y).sum() / s.sum())
    assert abs(sx - 4.0) < 0.25 and abs(sy / sx - q) < 0.04
    # without alignment the stack is round and wider
    c0 = st.make_cutouts(img, None, xs, ys, cfg)
    s0 = st.combine(c0['cube'], 'mean')
    assert math.sqrt((s0 * y * y).sum() / s0.sum()) / math.sqrt((s0 * x * x).sum() / s0.sum()) > 0.8


def test_radial_profile_and_surface_brightness():
    S = 61
    y, x = np.indices((S, S)) - 30
    r = np.hypot(x, y)
    img = 10 * np.exp(-r / 5.0)
    prof = st.Profiler(S, st.radial_bins(30, 6.0, 1.25))
    p, n = prof(img)
    ref = 10 * np.exp(-prof.rmid / 5.0)
    assert np.nanmax(np.abs(p / ref - 1)[1:]) < 0.06
    q = st.Profiler(S, st.radial_bins(30), q=0.5, pa=0.0)
    e = 10 * np.exp(-np.hypot(x, y / 0.5) / 5.0)
    p2, _ = q(e)
    assert np.nanmax(np.abs(p2 / (10 * np.exp(-q.rmid / 5.0)) - 1)[1:]) < 0.06
    mu, mue, lim = st.surface_brightness(np.array([1.0, 0.1]), np.array([0.1, 0.1]), 25.0, 0.06)
    assert abs(mu[0] - (25 - 2.5 * math.log10(1.0) + 5 * math.log10(0.06))) < 1e-9
    assert np.isnan(mu[1]) and abs(mue[0] - 1.0857 * 0.1) < 1e-6 and abs(lim[0] - (25 + 5 * math.log10(0.06) - 2.5 * math.log10(0.3))) < 1e-9


def test_catalog_neighbour_and_detection_masking_recover_halo():
    rng = np.random.default_rng(21)
    sh = (900, 900)
    img = rng.normal(0, 1, sh)
    n = 120
    xs = rng.uniform(60, 840, n)
    ys = rng.uniform(60, 840, n)
    nb_x = rng.uniform(20, 880, 500)
    nb_y = rng.uniform(20, 880, 500)
    for x, y in zip(nb_x, nb_y):
        stamp(img, x, y, 400.0, 2.0)
    yy, xx = np.mgrid[-40:41, -40:41]
    halo = 0.6 * np.maximum(np.hypot(xx, yy), 1) ** -1.0 * 8      # 0.6 at r=8 ... smooth
    for x, y in zip(xs, ys):
        xi, yi = int(round(x)), int(round(y))
        img[yi - 40:yi + 41, xi - 40:xi + 41] += halo
        stamp(img, x, y, 3000.0, 1.5)
    cat = np.column_stack([np.concatenate([xs, nb_x]), np.concatenate([ys, nb_y])])
    neigh = np.column_stack([cat, np.full(len(cat), 8.0)])
    out = {}
    for name, kw in (('none', {}), ('cat', dict(mask_cat_scale=1.0, mask_cat_min=8.0)), ('det', dict(mask_detect=3.0))):
        cfg = st.StackConfig(half=30, bkg='none', n_boot=20, ap_r=3.0, **kw)
        r = st.run_stack(img, None, xs, ys, cfg, null=150, avoid=cat, neigh=neigh, self_index=np.arange(n))
        prof = r['profile'] - r['null']['profile']
        out[name] = (prof, r['null']['profile'], r['profile_err'])
    rm = st.Profiler(61, st.radial_bins(30)).rmid
    k = np.where((rm > 15) & (rm < 24))[0]
    truth = 0.6 * 8 / rm[k]
    # the unmasked stack is dominated by neighbours (500 sources of flux 400 in 800x800: ~0.5/px mean) - the masked stacks recover the halo
    for name in ('cat', 'det'):
        assert np.all(np.abs(out[name][0][k] - truth) < 4 * out[name][2][k] + 0.03), name
    assert np.abs(out['none'][1][k]).mean() > 5 * np.abs(out['cat'][1][k]).mean() or np.abs(out['none'][1][k]).mean() > 0.05


def test_null_positions_avoid_sources():
    rng = np.random.default_rng(3)
    img = rng.normal(0, 1, (500, 500))
    cat = np.array([[250.0, 250.0]])
    cfg = st.StackConfig(half=20)
    x, y = st.null_positions(img, None, 100, cfg, avoid=cat)
    assert len(x) == 100 and np.hypot(x - 250, y - 250).min() > cfg.core_r + 3
    assert x.min() >= 22 and x.max() <= 500 - 23


# ---------------------------------------------------------------------------------------------------------------- CLI
@pytest.fixture(scope='module')
def clidata(tmp_path_factory):
    d = tmp_path_factory.mktemp('stk')
    img, xs, ys = field(seed=9, n=150, F=20.0, sg=2.0, sh=(600, 600))
    hdr = fits.Header()
    hdr['CTYPE1'], hdr['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    hdr['CRPIX1'] = hdr['CRPIX2'] = 300.5
    hdr['CRVAL1'], hdr['CRVAL2'] = 150.0, 2.0
    hdr['CD1_1'], hdr['CD2_2'] = -1e-4, 1e-4
    fits.PrimaryHDU(img.astype('f4'), header=hdr).writeto(d / 'i.fits')
    with open(d / 'c.tsv', 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tKRON_RADIUS\tMAG_AUTO\n')
        for i, (x, y) in enumerate(zip(xs, ys)):
            f.write('%d\t%.3f\t%.3f\t2\t3.5\t%.2f\n' % (i + 1, x + 1, y + 1, 20 + i * 0.02))
    with open(d / 'pos.txt', 'w') as f:
        f.write('# x y\n')
        for x, y in zip(xs[:80], ys[:80]):
            f.write('%.3f %.3f\n' % (x + 1, y + 1))
    from astropy.wcs import WCS
    ra, dec = WCS(hdr).all_world2pix(xs[:60], ys[:60], 0), None
    w = WCS(hdr)
    ra, dec = w.all_pix2world(xs[:60], ys[:60], 0)
    with open(d / 'radec.csv', 'w') as f:
        f.write('RA,DEC\n')
        for a, b in zip(ra, dec):
            f.write('%.7f,%.7f\n' % (a, b))
    return d, F_TRUE if False else 20.0


F_TRUE = 20.0


def run_cli(d, *args):
    r = subprocess.run([sys.executable, CLI, str(d / 'i.fits'), '--work', str(d / 'w')] + list(args), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    return r


def test_cli_catalog_columns_and_files(clidata):
    d, F = clidata
    r = run_cli(d, '--catalog', str(d / 'c.tsv'), '--sel-col', 'MAG_AUTO', '--sel-min=20', '--sel-max=22', '--n-boot', '40', '--null', '50',
                '--subtract-null', '--ap-r', '8', '--half', '16', '--mag-zeropoint', '25', '--pixel-scale', '0.06')
    lines = r.stdout.strip().split('\n')
    assert lines[0].split('\t') == ['NUMBER', 'ST_USED', 'ST_BKG', 'ST_NORM', 'ST_APFLUX', 'ST_MASKFRAC']
    assert 90 <= len(lines) - 1 <= 105                  # MAG_AUTO 20..22 -> 101 rows
    ap = [float(l.split('\t')[4]) for l in lines[1:] if l.split('\t')[4]]
    assert abs(np.mean(ap) - F) < 3 * np.std(ap) / math.sqrt(len(ap)) + 0.5
    w = d / 'w'
    for f in ('stack.fits', 'err.fits', 'profile.tsv', 'summary.json', 'plot.png'):
        assert (w / ('stacking_' + f)).stat().st_size > 100
    j = json.load(open(w / 'stacking_summary.json'))
    assert j['n_stacked'] > 90 and abs(j['aperture_sum'] - F) < 4 * j['aperture_err'] and 'aperture_sum_corrected' in j
    hdr = open(w / 'stacking_profile.tsv').readline().split()
    assert {'R', 'I', 'I_ERR', 'NULL', 'I_CORR', 'MU', 'MU_LIM3', 'SNR'} <= set(hdr)


def test_cli_positions_and_radec(clidata):
    d, F = clidata
    run_cli(d, '--positions', str(d / 'pos.txt'), '--n-boot', '20', '--ap-r', '8', '--half', '16')
    j1 = json.load(open(d / 'w' / 'stacking_summary.json'))
    run_cli(d, '--positions', str(d / 'radec.csv'), '--n-boot', '20', '--ap-r', '8', '--half', '16')
    j2 = json.load(open(d / 'w' / 'stacking_summary.json'))
    assert j1['n_stacked'] == 80 and j2['n_stacked'] == 60
    assert abs(j1['aperture_sum'] - F) < 4 * j1['aperture_err'] and abs(j2['aperture_sum'] - F) < 4 * j2['aperture_err']


def test_cli_numbers_isolate_and_errors(clidata):
    d, F = clidata
    r = run_cli(d, '--catalog', str(d / 'c.tsv'), '--numbers', '1-10', '--n-boot', '5', '--half', '12')
    assert len(r.stdout.strip().split('\n')) == 11
    r = subprocess.run([sys.executable, CLI, str(d / 'i.fits'), '--work', str(d / 'w'), '--catalog', str(d / 'c.tsv'), '--isolate', '2000', '--half', '12'], capture_output=True, text=True)
    assert r.returncode != 0 and 'no valid cutouts' in r.stderr
    r = run_cli(d, '--catalog', str(d / 'c.tsv'), '--isolate', '25', '--n-boot', '5', '--half', '12')
    used = [l.split('\t')[1] for l in r.stdout.strip().split('\n')[1:]]
    assert 0 < used.count('1') < len(used)         # only the isolated sources are used
    p = subprocess.run([sys.executable, CLI, str(d / 'i.fits'), '--work', str(d / 'w2')], capture_output=True, text=True)
    assert p.returncode != 0


def test_cli_rescale_rotate_columns(clidata):
    d, F = clidata
    with open(d / 'c2.tsv', 'w') as f:
        f.write(open(d / 'c.tsv').readline().rstrip('\n') + '\tFLUX_RADIUS\tTHETA_IMAGE\n')
        for i, l in enumerate(open(d / 'c.tsv').readlines()[1:]):
            f.write(l.rstrip('\n') + '\t2.0\t%g\n' % (i * 7 % 180))
    run_cli(d, '--catalog', str(d / 'c2.tsv'), '--scale-col', 'FLUX_RADIUS', '--scale-unit', '4', '--angle-col', 'THETA_IMAGE', '--n-boot', '10', '--half', '16')
    h = open(d / 'w' / 'stacking_profile.tsv').readline().split()
    assert 'R_UNIT' in h
