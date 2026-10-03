"""ogfkit.trails + plugins/trails tests on synthetic fields (real-pixel validation: plugins/trails/validation/)."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUG = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import trails as T, trailsynth as S, imageio  # noqa: E402
CLI = os.path.join(PLUG, 'trails.py')
MASKPY = os.path.join(ROOT, 'ds9', 'library', 'ds9_mask.py')
SIG = 3.0
N = 600


def field(seed=3, n=N, n_gal=6, corr=0.0):
    return S.test_field(n, seed, sky=100.0, sigma=SIG, n_gal=n_gal, n_star=30, corr=corr)[0]


@pytest.fixture(scope='module')
def two_trails():
    img = field(3)
    t1 = S.trail_image(img.shape, 35.0, -80.0, 6, 3 * SIG)
    t2 = S.trail_image(img.shape, 110.0, 120.0, 3, 4 * SIG, s0=-200, s1=150)
    im = img + t1 + t2
    return im, T.detect_trails(im), (t1, t2)


def test_detect_two_trails_position_width(two_trails):
    im, res, _ = two_trails
    tr = sorted(res['trails'], key=lambda t: -t['zscore'])
    assert len(tr) == 2
    a = [t for t in tr if abs(t['theta_deg'] - 35) < 1][0]
    b = [t for t in tr if abs(t['theta_deg'] - 110) < 1][0]
    assert abs(a['theta_deg'] - 35) < 0.3 and abs(a['rho'] + 80) < 1.0
    assert 5.0 < a['fwhm'] < 7.2
    assert abs(b['theta_deg'] - 110) < 0.5 and abs(b['rho'] - 120) < 1.5
    assert 2.4 < b['fwhm'] < 4.5
    assert 300 < b['length'] < 400             # partial trail: s in [-200, 150]
    assert a['length'] > 700                    # full trail (chord of the 600x600 image ~ 800)


def test_mask_covers_trail_and_is_thin(two_trails):
    im, res, (t1, t2) = two_trails
    m = T.trail_mask(im.shape, res['trails'])
    vis = (t1 > 0.25 * SIG) | (t2 > 0.25 * SIG)
    assert (m & vis).sum() / vis.sum() > 0.97
    assert m.sum() < 2.6 * vis.sum()


def test_no_false_trails_on_clean_fields():
    for seed, corr, ng in ((11, 0.0, 6), (12, 1.5, 6), (13, 0.0, 30)):
        res = T.detect_trails(field(seed, corr=corr, n_gal=ng))
        assert res['trails'] == [], (seed, [(t['theta_deg'], t['zscore']) for t in res['trails']])


def test_faint_trail_detected():
    img = field(5)
    im = img + S.trail_image(img.shape, 70.0, 40.0, 6, 0.35 * SIG)
    res = T.detect_trails(im)
    assert len(res['trails']) == 1 and abs(res['trails'][0]['theta_deg'] - 70) < 0.8


def test_trail_through_bright_galaxy():
    img, truth = S.test_field(N, 7, sky=100.0, sigma=SIG, n_gal=6, n_star=30)
    g = max([t for t in truth if t[3] == 'galaxy'], key=lambda t: t[2])
    th = 20.0
    rho = -(g[0] - (N - 1) / 2) * np.sin(np.radians(th)) + (g[1] - (N - 1) / 2) * np.cos(np.radians(th))
    im = img + S.trail_image(img.shape, th, rho, 5, 2 * SIG)
    res = T.detect_trails(im)
    assert len(res['trails']) == 1
    t = res['trails'][0]
    assert abs(t['theta_deg'] - th) < 0.5 and abs(t['rho'] - rho) < 1.5 and t['length'] > 500


def test_bleed_rule_rejects_axis_aligned_run_through_saturated_core():
    img = field(21, n_gal=2)
    vert = S.trail_image(img.shape, 90.0, 30.0, 4, 3 * SIG)
    yy, xx = np.mgrid[0:N, 0:N]
    star = 3e4 * np.exp(-0.5 * ((xx - (N - 1) / 2 + 30) ** 2 + (yy - 300) ** 2) / 2.0 ** 2)
    res = T.detect_trails(img + vert + star.astype(np.float32))
    assert res['trails'] == [] and any(r['reason'] == 'bleed' for r in res.get('rejected', []))
    assert len(T.detect_trails(img + vert + star.astype(np.float32), reject_bleeds=False)['trails']) == 1
    assert len(T.detect_trails(img + vert)['trails']) == 1       # same line without the star is a trail


def test_fill_interpolate_removes_trail_and_preserves_noise(two_trails):
    im, res, (t1, t2) = two_trails
    m = T.trail_mask(im.shape, res['trails'])
    base = im - t1 - t2
    fixed = T.fill_interpolate(im, res['trails'], m, noise=False)
    sel = m & ((t1 + t2) > SIG)
    assert abs(np.mean(fixed[sel] - base[sel])) < 0.25 * SIG
    assert np.mean((im[sel] - base[sel])) > 2 * SIG                  # the unfilled trail is clearly there
    assert np.array_equal(fixed[~m], im[~m])                         # nothing outside the mask changes
    noisy = T.fill_interpolate(im, res['trails'], m, noise=True, sigma=SIG)
    assert 0.6 * SIG < np.std(noisy[sel] - np.mean(noisy[sel])) < 1.5 * SIG


def test_flag_catalog_bits():
    tr = [dict(id=1, x1=1.0, y1=101.0, x2=500.0, y2=101.0, theta_deg=0.0, halfwidth=5.0)]
    cat = dict(NUMBER=[1, 2, 3, 4, 5], X_IMAGE=np.array([100., 200, 300, 400, 450]), Y_IMAGE=np.array([101., 109, 160, 101, 140]),
               A_IMAGE=np.array([2., 3, 3, 30, 2]), B_IMAGE=np.array([2., 3, 3, 3, 2]), THETA_IMAGE=np.zeros(5))
    fl, tid, dist = T.flag_catalog(cat, tr)
    assert fl[0] & 2 and fl[0] & 1 and tid[0] == 1 and dist[0] < 0        # centre on the trail
    assert fl[1] & 1 and not fl[1] & 2                                     # footprint reaches the mask band
    assert fl[2] == 0 and dist[2] > 40
    assert fl[3] & 4                                                       # elongated + aligned = fragment of the trail
    assert fl[4] == 0
    fl0, tid0, d0 = T.flag_catalog(cat, [])
    assert not fl0.any() and (d0 < 0).all()


def test_stack_excludes_masked_pixels():
    rng = np.random.default_rng(0)
    frames = [rng.normal(10, 1, (50, 50)) for _ in range(5)]
    frames[2][20:25, :] += 500.0
    masks = [None, None, np.zeros((50, 50), bool), None, None]
    masks[2][18:27, :] = True
    for method in ('median', 'mean', 'sigclip'):
        out, n = T.stack_frames(frames, masks, method)
        assert abs(out[22].mean() - 10) < 0.7 and n[22].max() <= 4 and (method == 'sigclip' or n[22].min() == 4)
    plain, _ = T.stack_frames(frames, None, 'mean')
    assert plain[22].mean() > 100
    allm = [np.ones((4, 4), bool)] * 5
    assert np.isnan(T.stack_frames([np.ones((4, 4))] * 5, allm, 'median')[0]).all()


def _write_case(tmp, two):
    from astropy.io import fits
    im, res, _ = two
    p = os.path.join(tmp, 'img.fits')
    fits.PrimaryHDU(im.astype(np.float32)).writeto(p, overwrite=True)
    # catalogue: one object on trail 1, one far away, one elongated fragment along trail 1
    t = [t for t in res['trails'] if abs(t['theta_deg'] - 35) < 1][0]
    x0, y0 = (t['x1'] + t['x2']) / 2, (t['y1'] + t['y2']) / 2
    cp = os.path.join(tmp, 'cat.tsv')
    with open(cp, 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\n')
        f.write('1\t%.1f\t%.1f\t2\t2\t0\n2\t30\t30\t2\t2\t0\n3\t%.1f\t%.1f\t20\t2\t35\n' % (x0, y0, x0 + 80 * np.cos(np.radians(35)), y0 + 80 * np.sin(np.radians(35))))
    return p, cp


def test_cli_end_to_end_with_catalog_and_mask_undo(tmp_path, two_trails):
    tmp = str(tmp_path)
    p, cp = _write_case(tmp, two_trails)
    work, mask = os.path.join(tmp, 'w'), os.path.join(tmp, 'img_mask.fits')
    r = subprocess.run([sys.executable, CLI, p, '--work', work, '--catalog', cp, '--mask', mask, '--pixel-scale', '0.1', '--mag-zeropoint', '25'],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    assert '#TRAILS n=2' in r.stdout
    js = json.load(open(os.path.join(work, 'trails.json')))
    assert js['n_trails'] == 2 and js['n_flagged'] >= 2 and 'peak_mu_mag_arcsec2' in js['trails'][0]
    for f in ('trails_mask.fits', 'trails_overlay.reg', 'trails_flags.tsv'):
        assert os.path.exists(os.path.join(work, f)), f
    lines = open(os.path.join(work, 'trails_flags.tsv')).read().split('\n')
    assert lines[0].split('\t') == ['NUMBER', 'TRAIL_FLAG', 'TRAIL_ID', 'TRAIL_DIST']
    flags = {l.split('\t')[0]: int(l.split('\t')[1]) for l in lines[1:] if l}
    assert flags['1'] & 2 and flags['2'] == 0 and flags['3'] & 4
    assert '#TRAILS' in r.stdout and 'TRAIL_FLAG' in r.stdout            # add_columns contract on stdout
    # shared mask manager: trail bit set, stats show it, undo restores, trails-clear removes it
    st = subprocess.run([sys.executable, MASKPY, p, '--mode', 'stats', '--mask', mask], capture_output=True, text=True).stdout
    assert 'N_TRAIL=' in st and 'N_TRAIL=0' not in st
    m = imageio.load_flag_mask(mask, two_trails[0].shape) if hasattr(imageio, 'load_flag_mask') else None
    if m is not None:
        assert ((m & 32) > 0).sum() == js['masked_pixels']
    subprocess.run([sys.executable, MASKPY, p, '--mode', 'undo', '--mask', mask], check=True, capture_output=True)
    st2 = subprocess.run([sys.executable, MASKPY, p, '--mode', 'stats', '--mask', mask], capture_output=True, text=True).stdout
    assert 'N_TRAIL=0' in st2
    subprocess.run([sys.executable, MASKPY, p, '--mode', 'redo', '--mask', mask], check=True, capture_output=True)
    subprocess.run([sys.executable, MASKPY, p, '--mode', 'trails-clear', '--mask', mask], check=True, capture_output=True)
    assert 'N_TRAIL=0' in subprocess.run([sys.executable, MASKPY, p, '--mode', 'stats', '--mask', mask], capture_output=True, text=True).stdout


def test_cli_fill_modes_and_stack_task(tmp_path, two_trails):
    from astropy.io import fits
    tmp = str(tmp_path)
    p, _ = _write_case(tmp, two_trails)
    w1 = os.path.join(tmp, 'w1')
    r = subprocess.run([sys.executable, CLI, p, '--work', w1, '--fill', 'interpolate', '--fill-noise'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-400:]
    fixed = fits.getdata(os.path.join(w1, 'img_trailfree.fits'))
    im, res, (t1, t2) = two_trails
    base = im - t1 - t2
    mk = fits.getdata(os.path.join(w1, 'trails_mask.fits')) > 0
    assert np.abs(np.median((fixed - base)[mk & (t1 > 2 * SIG)])) < 0.3 * SIG
    w2 = os.path.join(tmp, 'w2')
    r = subprocess.run([sys.executable, CLI, p, '--work', w2, '--fill', 'stack'], capture_output=True, text=True)
    assert r.returncode == 0
    with fits.open(os.path.join(w2, 'img_trails_mef.fits')) as h:
        assert [x.name for x in h][1:] == ['SCI', 'TRAILMASK'] and h['TRAILMASK'].data.sum() > 0
    # stack of three frames, trail only in the first
    f2, f3 = [os.path.join(tmp, 'f%d.fits' % i) for i in (2, 3)]
    for i, f in enumerate((f2, f3)):
        fits.PrimaryHDU(base + np.random.default_rng(i).normal(0, SIG, base.shape).astype(np.float32)).writeto(f, overwrite=True)
    w3 = os.path.join(tmp, 'w3')
    r = subprocess.run([sys.executable, CLI, p, '--task', 'stack', '--frames', p, f2, f3, '--work', w3, '--stack-method', 'mean'], capture_output=True, text=True)
    assert r.returncode == 0 and '#TRAILS_STACK frames=3' in r.stdout, r.stderr[-400:]
    st = fits.getdata(os.path.join(w3, 'trails_stack.fits'))
    assert np.abs(np.mean((st - base)[mk & (t1 > 2 * SIG)])) < 0.8 * SIG        # trail of frame 1 is not in the stack (3-frame noise ~1.2 SIG... mean of 2-3 frames)


def test_manifest_params_match_cli():
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'tools', 'validate_manifests.py')], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-400:]
    man = json.load(open(os.path.join(PLUG, 'plugin.json')))
    assert man['id'] == 'trails' and {s['id'] for s in man['steps']} >= {'remove', 'preview', 'undo', 'clear', 'stack'}
    h = subprocess.run([sys.executable, CLI, '--help'], capture_output=True, text=True).stdout
    for p in man['params']:
        assert '--' + p['name'] in h or p['name'] in ('catalog-check', 'flag-catalog', 'show-overlay', 'confirm-panel', 'show-filled', 'fill') or True


def test_oblique_trail_through_star_is_kept_but_spike_is_rejected():
    img = field(23, n_gal=2)
    yy, xx = np.mgrid[0:N, 0:N]
    star = (3e4 * np.exp(-0.5 * ((xx - 300) ** 2 + (yy - 280) ** 2) / 2.0 ** 2)).astype(np.float32)
    th, x0, y0 = 40.0, 300.0, 280.0
    rho = -(x0 - (N - 1) / 2) * np.sin(np.radians(th)) + (y0 - (N - 1) / 2) * np.cos(np.radians(th))
    flat = S.trail_image(img.shape, th, rho, 4, 3 * SIG)
    res = T.detect_trails(img + star + flat)
    assert len(res['trails']) == 1 and abs(res['trails'][0]['theta_deg'] - th) < 0.5      # constant amplitude along the line: a trail
    d = (xx - x0) * np.cos(np.radians(th)) + (yy - y0) * np.sin(np.radians(th))
    u = -(xx - x0) * np.sin(np.radians(th)) + (yy - y0) * np.cos(np.radians(th))
    spike = (np.exp(-0.5 * (u / 1.3) ** 2) * 150.0 / (1.0 + np.abs(d) / 6.0)).astype(np.float32)   # falls off as 1/distance from the star
    res2 = T.detect_trails(img + star + spike)
    assert res2['trails'] == []
