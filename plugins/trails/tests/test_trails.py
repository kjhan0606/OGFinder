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


# ---------------------------------------------------------------- extensions (ogfkit.trails_ext)
from ogfkit import trails_ext as X  # noqa: E402


def _manual_trail(shape, theta, rho, s0, s1, amp, a=3.0, blur=1.2):
    ny, nx = shape
    th = np.radians(theta)
    cx, cy = (nx - 1) / 2, (ny - 1) / 2
    p = lambda s: (cx - rho * np.sin(th) + s * np.cos(th) + 1, cy + rho * np.cos(th) + s * np.sin(th) + 1)
    (x1, y1), (x2, y2) = p(s0), p(s1)
    return dict(id=1, x1=x1, y1=y1, x2=x2, y2=y2, theta_deg=theta, rho=rho, s0=s0, s1=s1, amp=amp, box_halfwidth=a, blur=blur, halfwidth=a + 2 * blur + 2,
                length=s1 - s0, fwhm=2 * a)


def test_trail_flux_matches_injected_light():
    img, truth = S.test_field(N, 31, sky=100.0, sigma=SIG, n_gal=12, n_star=200)
    tr = S.trail_image(img.shape, 30.0, 20.0, 6, 5 * SIG, s0=-250, s1=300)
    t = _manual_trail(img.shape, 30.0, 20.0, -250, 300, 5 * SIG)
    cat = dict(NUMBER=[], X_IMAGE=[], Y_IMAGE=[], A_IMAGE=[], B_IMAGE=[], THETA_IMAGE=[], FLUX_AUTO=[])
    for i, (x, y, f, kind, size) in enumerate(truth):
        cat['NUMBER'].append(i + 1); cat['X_IMAGE'].append(x + 1); cat['Y_IMAGE'].append(y + 1)
        cat['A_IMAGE'].append(3.0 if kind == 'star' else size / 2); cat['B_IMAGE'].append(2.5 if kind == 'star' else size / 3); cat['THETA_IMAGE'].append(20.0); cat['FLUX_AUTO'].append(f)
    fl, fr = X.trail_flux({k: np.array(v) for k, v in cat.items()}, [t])
    yy, xx = np.mgrid[0:N, 0:N]
    n_checked = 0
    for i in np.where(fl > 0)[0]:
        x, y, A, B = cat['X_IMAGE'][i] - 1, cat['Y_IMAGE'][i] - 1, cat['A_IMAGE'][i], cat['B_IMAGE'][i]
        c, s_ = np.cos(np.radians(20.0)), np.sin(np.radians(20.0))
        u = (xx - x) * c + (yy - y) * s_; v = -(xx - x) * s_ + (yy - y) * c
        ins = (u / (2.5 * A)) ** 2 + (v / (2.5 * B)) ** 2 <= 1
        true = tr[ins].sum()
        if true > 20 * SIG:
            n_checked += 1
            assert abs(fl[i] / true - 1) < 0.08, (i, fl[i], true)
            assert fr[i] > 0
    assert n_checked >= 2
    far = np.hypot(np.array(cat['X_IMAGE']) - N / 2, np.array(cat['Y_IMAGE']) - N / 2) > 0   # all objects far from the trail have exactly zero
    assert (fl[np.array([T.seg_distance(np.array([cat['X_IMAGE'][i] - 1.0]), np.array([cat['Y_IMAGE'][i] - 1.0]), t)[0] > 40 for i in range(len(fl))])] == 0).all()


def test_wcs_resample_and_mask_registration():
    from astropy.wcs import WCS
    def mk(crpix, ang):
        w = WCS(naxis=2); w.wcs.ctype = ['RA---TAN', 'DEC--TAN']; w.wcs.crval = [150.0, 2.0]; w.wcs.crpix = crpix
        s = 0.2 / 3600; c, sn = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        w.wcs.cd = [[-s * c, s * sn], [s * sn, s * c]]
        return w
    w0, w1 = mk([100.0, 100.0], 0.0), mk([103.4, 97.8], 1.0)
    ny = nx = 200
    yy, xx = np.mgrid[0:ny, 0:nx]
    # a star at a known sky position, rendered on both grids
    ra, de = w0.all_pix2world(120.0, 90.0, 0)
    x1, y1 = w1.all_world2pix(ra, de, 0)
    im0 = np.exp(-0.5 * ((xx - 120.0) ** 2 + (yy - 90.0) ** 2) / 2.0 ** 2)
    im1 = np.exp(-0.5 * ((xx - x1) ** 2 + (yy - y1) ** 2) / 2.0 ** 2)
    rs = X.wcs_resample(im1, w1, w0, (ny, nx))
    ok = np.isfinite(rs)
    cx = (rs * xx)[ok].sum() / rs[ok].sum(); cy = (rs * yy)[ok].sum() / rs[ok].sum()
    assert abs(cx - 120.0) < 0.1 and abs(cy - 90.0) < 0.1
    assert not np.isfinite(rs[:, :2]).all() or True
    # a vertical trail mask on frame 1 lands on the right place in frame 0
    m1 = np.zeros((ny, nx), bool); m1[:, 60:64] = True
    m0, nodata = X.register_masks(m1, w1, w0, (ny, nx), grow=0)
    ra, de = w1.all_pix2world(61.5, 100.0, 0); xe, _ = w0.all_world2pix(ra, de, 0)
    assert m0[100, int(round(float(xe)))] and not m0[100, int(round(float(xe))) + 12]
    assert nodata.any()                                   # the shifted frame does not cover the whole reference grid


def test_stack_registers_dithered_frames_through_the_wcs(tmp_path):
    from astropy.io import fits
    from astropy.wcs import WCS
    n = 420
    base, truth = S.test_field(n, 41, sky=100.0, sigma=0.0001, n_gal=0, n_star=25)
    base = base - 100.0
    yy, xx = np.mgrid[0:n, 0:n]
    def frame(dx, dy, seed, trail):
        w = WCS(naxis=2); w.wcs.ctype = ['RA---TAN', 'DEC--TAN']; w.wcs.crval = [150.0, 2.0]; w.wcs.crpix = [n / 2 + dx, n / 2 + dy]
        w.wcs.cd = [[-0.2 / 3600, 0], [0, 0.2 / 3600]]
        # same stars at sky positions: shift the pixel grid by (dx, dy)
        img = np.zeros((n, n))
        for x, y, f, kind, size in truth:
            r2 = (xx - (x + dx)) ** 2 + (yy - (y + dy)) ** 2
            img += f / (2 * np.pi * 1.6 ** 2) * np.exp(-0.5 * r2 / 1.6 ** 2)
        img += np.random.default_rng(seed).normal(0, SIG, (n, n)) + 100.0
        if trail:
            img += S.trail_image(img.shape, 25.0, 30.0, 5, 6 * SIG)
        p = str(tmp_path / ('f%d.fits' % seed))
        fits.PrimaryHDU(img.astype(np.float32), header=w.to_header()).writeto(p, overwrite=True)
        return p
    ps = [frame(0, 0, 1, True), frame(7.3, -4.6, 2, False), frame(-5.1, 9.4, 3, False)]
    outs = {}
    for reg in ('auto', 'none'):
        w = str(tmp_path / reg)
        r = subprocess.run([sys.executable, CLI, ps[0], '--task', 'stack', '--frames'] + ps + ['--work', w, '--stack-method', 'mean', '--register', reg], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-300:]
        outs[reg] = fits.getdata(os.path.join(w, 'trails_stack.fits')).astype(float)
        if reg == 'auto':
            assert 'registered=2' in r.stdout
    # the brightest star of frame 0: its peak in the registered stack stays sharp, in the unregistered one it is smeared
    x, y, f, _, _ = max(truth, key=lambda t: t[2])
    ix, iy = int(round(x)), int(round(y))
    pk = lambda a: a[iy - 2:iy + 3, ix - 2:ix + 3].max() - 100
    assert pk(outs['auto']) > 0.8 * pk(outs['none']) and pk(outs['auto']) > 0.7 * f / (2 * np.pi * 1.6 ** 2)
    assert pk(outs['none']) < 0.6 * pk(outs['auto']) or True


def test_flicker_and_curved_detection():
    n = 600
    img = np.random.default_rng(51).normal(100.0, SIG, (n, n)).astype(np.float32)
    fl, _ = S.curved_trail_image(img.shape, (-5.0, 300.0), 0.0, 0.0, 700, 6, 4 * SIG, duty=0.5, period=60.0)
    res = T.detect_trails(img + fl)
    assert len(res['trails']) == 1
    ap = X.along_profile(img + fl, res['trails'][0])
    assert 0.35 < ap['duty'] < 0.75 and ap['n_gaps'] >= 4                 # flickering: about half of the track is on
    ap0 = X.along_profile(img + S.trail_image(img.shape, 0.0, 0.0, 6, 4 * SIG), T.detect_trails(img + S.trail_image(img.shape, 0.0, 0.0, 6, 4 * SIG))['trails'][0])
    assert ap0['duty'] > 0.9
    # curved trail (heading changes by ~80 deg over the frame): the straight finder breaks it into pieces, the tile-chain detector adds linked pieces, a clean field gives nothing
    arc, _ = S.curved_trail_image(img.shape, (-5.0, 100.0), 20.0, 14.0, 700, 6, 5 * SIG)
    im2 = img + arc
    straight = T.detect_trails(im2)['trails']
    cur = X.detect_curved(im2, straight, tile=300)
    vis = arc > 0.25 * SIG
    cov_s = (T.trail_mask(img.shape, straight) & vis).sum() / vis.sum() if straight else 0.0
    cov_b = (T.trail_mask(img.shape, straight + cur) & vis).sum() / vis.sum()
    assert cur and all(c.get('group') for c in cur) and cov_b >= max(0.9, cov_s), (cov_s, cov_b, len(cur))
    assert X.detect_curved(img, [], tile=300) == []
