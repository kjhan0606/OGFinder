"""DAOPHOT plugin tests: synthetic crowded fields with known truth (ogfkit.synth)."""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
CLI = os.path.join(ROOT, 'plugins', 'daophot', 'daophot.py')
sys.path.insert(0, ROOT)
from ogfkit import synth, psfmodel as pm, daophot as dp, models  # noqa: E402

NOISE = 0.5


@pytest.fixture(scope='module')
def field():
    img, t = synth.star_field((600, 600), n=1200, mag_range=(15.5, 23.0), slope=0.3, noise=NOISE, seed=3)
    return img, t


@pytest.fixture(scope='module')
def run(field):
    img, t = field
    res, resid, model = dp.run_pipeline(img, fwhm=3.3, n_workers=4)
    return res, resid, model


def match_truth(res, t, r=1.0):
    S = res['stars']
    tr = cKDTree(np.c_[[s['x'] for s in S], [s['y'] for s in S]])
    d, j = tr.query(np.c_[t['x'], t['y']])
    return d, j


def true_psf_stamp(x, y, shape, size=21):
    p = synth.psf_params(x, y, shape)
    return models.moffat_psf(size, p['fwhm'], beta=3.0, e=p['e'], pa_deg=30, nsub=5)


# ---------------------------------------------------------------------------------------------------------------- PSF model
def test_psf_spatial_variation_recovered(field, run):
    img, t = field
    _, _, model = run
    dw, de = [], []
    for x in (60, 300, 540):
        for y in (60, 300, 540):
            a = pm.psf_shape(model.stamp(x, y, 0, 0, 21))
            b = pm.psf_shape(true_psf_stamp(x, y, img.shape))
            dw.append(a['fwhm'] / b['fwhm'] - 1)
            de.append(a['e'] - b['e'])
    # FWHM varies 2.8 -> 3.8 px (+36 %) and e 0.04 -> 0.14 across the field
    assert np.max(np.abs(dw)) < 0.05, dw
    assert np.median(np.abs(de)) < 0.02 and np.max(np.abs(de)) < 0.07, de      # the field corners are polynomial extrapolations
    d = pm.diagnostics(model, 5)
    assert d['fwhm'][2][4] > d['fwhm'][2][0] * 1.15          # the gradient in x is seen
    assert d['e'][4][2] > d['e'][0][2]                         # and the one in y


@pytest.mark.parametrize('kind', ['gaussian', 'moffat', 'lorentz', 'penny'])
def test_analytic_shapes(kind):
    # constant Moffat PSF, FWHM 3.2, beta 3: the moffat fit must recover it, the others must give a sensible FWHM
    psf = models.moffat_psf(25, 3.2, beta=3.0, e=0.0, nsub=5)
    r = pm.fit_analytic_shape(psf, kind)
    m = pm.analytic_model(r, 25, 1, (0, 1), (0, 1))
    tol = {'moffat': 0.03, 'penny': 0.06, 'gaussian': 0.25, 'lorentz': 0.5}[kind]       # a function that does not match the true shape only gets the core roughly right
    assert abs(pm.psf_shape(m.stamp(0, 0))['fwhm'] / pm.psf_shape(psf)['fwhm'] - 1) < tol
    if kind == 'moffat':
        assert abs(r['fx'] - 3.2) < 0.1 and abs(r['beta'] - 3.0) < 0.5


def test_few_stars_fallback(field):
    img, t = field
    for n, mode in ((0, 'gaussian_prior'), (2, 'moffat_fit')):
        xy = t['x'][np.argsort(t['mag'])][:n], t['y'][np.argsort(t['mag'])][:n]
        mdl, info = pm.build_psf_model(img, xy=np.c_[xy] if n else np.zeros((0, 2)), fwhm_prior=3.3) if n else pm.build_psf_model(img[:40, :40] * 0 + 50, fwhm_prior=3.3)
        assert info['mode'] == mode, info['mode']
        assert abs(mdl.stamp(0, 0, 0, 0).sum() - 1) < 1e-6
    # 8 stars -> constant empirical model; 30 -> order 1; order cap respected
    order = np.argsort(t['mag'])
    for n, deg in ((8, 0), (30, 1)):
        P = np.c_[t['x'][order][:n], t['y'][order][:n]]
        mdl, info = pm.build_psf_model(img, xy=P, fwhm_prior=3.3, degree=3)
        assert info['mode'] == 'empirical' and info['degree'] == deg, (n, info['degree'])


def test_psf_json_fits_roundtrip(run, tmp_path):
    _, _, model = run
    for name in ('m.json', 'm.fits'):
        model.save(str(tmp_path / name))
        m2 = pm.load_model(str(tmp_path / name))
        assert np.abs(m2.stamp(100, 200, 0.2, -0.3) - model.stamp(100, 200, 0.2, -0.3)).max() < 1e-6
    json.dumps(model.to_dict())


# -------------------------------------------------------------------------------------------------------------- FIND / PHOT
def test_find_sharpness_roundness_cuts():
    rng = np.random.default_rng(5)
    img = rng.normal(0, 1, (200, 200)).astype(np.float32) + 100
    yy, xx = np.mgrid[:200, :200]
    g = lambda x, y, s, a: a * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s ** 2)
    img += g(50, 50, 1.4, 60) + g(150, 50, 1.4, 60)                  # two stars
    img[100, 100] += 90                                              # cosmic-ray-like hot pixel
    img += g(50, 150, 5.0, 25)                                       # extended blob
    img += 60 * np.exp(-0.5 * (((xx - 150) / 0.9) ** 2 + ((yy - 150) / 3.2) ** 2))   # elongated
    bkg, rms = pm.background(img)
    f_all = dp.find(img, bkg, rms, fwhm=3.3, thresh=5, sharp=(-10, 10), round_=(-10, 10))
    f_cut = dp.find(img, bkg, rms, fwhm=3.3, thresh=5, sharp=(0.3, 0.95), round_=(-0.5, 0.5))
    near = lambda f, x, y: bool(np.any(np.hypot(f['x'] - x, f['y'] - y) < 2))
    assert near(f_all, 100, 100) and near(f_all, 150, 150)
    assert near(f_cut, 50, 50) and near(f_cut, 150, 50)
    assert not near(f_cut, 100, 100) and not near(f_cut, 50, 150) and not near(f_cut, 150, 150)


def test_phot_apertures_and_sky_modes():
    rng = np.random.default_rng(2)
    n = 301
    yy, xx = np.mgrid[:n, :n]
    sig = 1.8
    flux = 5000.0
    star = flux / (2 * np.pi * sig ** 2) * np.exp(-0.5 * ((xx - 150.3) ** 2 + (yy - 149.6) ** 2) / sig ** 2)
    sky = 200.0
    img = sky + star + rng.normal(0, 1.0, (n, n))
    clean = dp.phot(img, [(150.3, 149.6)], radii=(3, 6, 10), sky_inner=12, sky_outer=18, sky_mode='median', zp=25.0)[0]
    enclosed = [1 - math.exp(-r * r / (2 * sig * sig)) for r in (3, 6, 10)]
    for k, e in enumerate(enclosed):
        assert abs(clean['flux'][k] / (flux * e) - 1) < 0.02, (k, clean['flux'][k], flux * e)
    # contaminate the sky annulus with faint neighbours (skews the sky distribution upward)
    for _ in range(25):
        x, y = rng.uniform(130, 170, 2)
        if 11 < math.hypot(x - 150.3, y - 149.6) < 19:
            img += 4 * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / 1.5 ** 2)
    out = {}
    for mode in ('median', 'mode', 'mean'):
        p = dp.phot(img, [(150.3, 149.6)], radii=(3, 6, 10), sky_inner=12, sky_outer=18, sky_mode=mode, zp=25.0)[0]
        out[mode] = p
    assert abs(out['mode']['sky'] - sky) <= abs(out['mean']['sky'] - sky) + 1e-9      # mode least affected by faint contaminants
    assert abs(out['median']['sky'] - sky) < 1.0 and abs(out['mode']['sky'] - sky) < abs(out['median']['sky'] - sky)
    assert out['median']['magerr'][2] > 0 and np.isfinite(out['median']['mag'][2])


def test_pickpsf_auto_manual_and_neighbour_rejection():
    n = 6
    st = dict(x=np.array([50., 100, 104, 150, 30, 190]), y=np.array([50., 100, 100, 150, 100, 100]), flux=np.array([1000., 900, 500, 800, 700, 600]),
              peak=np.full(n, 50.), sharp=np.full(n, 0.6), round1=np.zeros(n), snr=np.full(n, 100.))
    xy, info = dp.pick_psf(st, (200, 200), n_max=10, half=10)
    pts = {(int(p[0]), int(p[1])) for p in xy}
    assert (100, 100) not in pts and (50, 50) in pts and (150, 150) in pts          # (100,100) has a close neighbour with >10 % of its flux
    assert info['reasons']['neighbour'] >= 1 and info['reasons']['edge'] >= 1       # (30,100) and (190,100) are within `half` of the edge
    xy2, info2 = dp.pick_psf(st, (200, 200), n_max=10, half=10, add=[(100.4, 100.2)], remove=[(50, 50)])
    pts2 = {(int(p[0]), int(p[1])) for p in xy2}
    assert (100, 100) in pts2 and (50, 50) not in pts2 and info2['n_added'] == 1 and info2['n_removed'] == 1


# ----------------------------------------------------------------------------------------------------------------- ALLSTAR
def test_allstar_flux_position_completeness_chi(field, run):
    img, t = field
    res, resid, model = run
    d, j = match_truth(res, t)
    S = res['stars']
    bright = (t['mag'] < 19.0) & (d < 1.0)
    allb = t['mag'] < 19.0
    comp_b = bright.sum() / allb.sum()
    assert comp_b > 0.97, comp_b
    mm = (t['mag'] < 20.0)
    comp = ((d < 1.0) & mm).sum() / mm.sum()
    assert comp > 0.85, comp
    dm = np.array([S[k]['mag_apc'] for k in j[bright]]) - t['mag'][bright]
    assert abs(np.median(dm)) < 0.03, np.median(dm)                    # photometric bias with the aperture correction applied
    assert 1.4826 * np.median(np.abs(dm - np.median(dm))) < 0.06
    pos = d[bright]
    assert np.median(pos) < 0.06, np.median(pos)                       # positional error (px) of bright stars
    chi = np.array([S[k]['chi'] for k in j[(d < 1.0) & (t['mag'] < 20.5)]])
    assert 0.7 < np.median(chi) < 1.4, np.median(chi)
    # unmatched (false) detections
    unm = len(S) - len(set(j[d < 1.5]))
    assert unm <= 0.05 * len(S) + 2, unm
    # the residual is at the noise level
    core = resid[60:540, 60:540]
    assert abs(np.std(core) / NOISE - 1) < 0.08, np.std(core)
    # error estimates are honest for faint-ish stars: pull of flux errors (the systematic PSF-flux bias of -1.8 % of the stamp-truncated wings is ~0.5 sigma at mag 20)
    sel = (d < 1.0) & (t['mag'] > 19.5) & (t['mag'] < 20.5)
    pulls = np.array([(S[k]['flux'] - t['flux'][i]) / S[k]['fluxerr'] for i, k in zip(np.where(sel)[0], j[sel])])
    print("DAOPHOT pulls: median %.3f, frac flux err median %.4f" % (np.median(pulls), np.median([(S[k]["flux"] / t["flux"][i] - 1) for i, k in zip(np.where(sel)[0], j[sel])])))
    assert abs(np.median(pulls)) < 0.7 and 0.6 < 1.4826 * np.median(np.abs(pulls - np.median(pulls))) < 1.6, (np.median(pulls), np.std(pulls))


def test_spatially_varying_psf_beats_constant(field):
    img, t = field
    bkg, rms = pm.background(img)
    scat = {}
    for order in (0, 2):
        mdl, _ = pm.build_psf_model(img, bkg=bkg, rms=rms, fwhm_prior=3.3, degree=order, snr_min=20)
        res, _ = dp.allstar(img, mdl, bkg=bkg, rms=rms, fwhm=3.3, n_iter=2, n_workers=4)
        d, j = match_truth(res, t)
        sel = (t['mag'] < 19.5) & (d < 1.0)
        S = res['stars']
        r = np.array([S[k]['flux'] for k in j[sel]]) / t['flux'][sel] - 1
        scat[order] = (np.std(r - np.median(r)), np.mean(np.abs(r)))
    assert scat[2][0] < 0.8 * scat[0][0], scat           # order-2 PSF: clearly smaller flux scatter than a constant PSF


def test_flags_and_cuts(run):
    res, _, _ = run
    S = res['stars']
    assert all(('flag' in s and 'good' in s) for s in S)
    assert all((s['good'] == ((s['flag'] & dp.BAD_MASK) == 0)) for s in S)
    assert any(s['flag'] & dp.FLAG_NEW for s in S)       # stars found on the residual iterations exist and are labelled
    bad = [s for s in S if not s['good']]
    assert len(bad) < 0.1 * len(S)


def test_substar_addstar_roundtrip(field, run):
    img, t = field
    res, resid, model = run
    S = res['stars']
    sub = dp.substar(img, model, S)
    assert abs(np.std(sub[60:540, 60:540]) / NOISE - 1) < 0.08
    keep = [S[0]['id']]
    sub2 = dp.substar(img, model, S, keep=keep)
    s0 = S[0]
    ix, iy = int(round(s0['x'])), int(round(s0['y']))
    assert sub2[iy, ix] - sub[iy, ix] > 0.3 * s0['flux'] * model.stamp(s0['x'], s0['y']).max()      # the kept star is still there
    pos = dp.random_positions(img.shape, 30, margin=15, seed=4, avoid=np.c_[t['x'], t['y']], min_dist=12)
    mags = np.full(len(pos), 19.0)
    out, truth = dp.addstar(img, model, pos, mags)
    added = (out - img).sum()
    assert abs(added / sum(tr['flux'] for tr in truth) - 1) < 0.01
    r2, _ = dp.allstar(out, model, fwhm=3.3, n_iter=2, n_workers=4)
    S2 = r2['stars']
    tr = cKDTree(np.c_[[s['x'] for s in S2], [s['y'] for s in S2]])
    d, j = tr.query(np.array(pos))
    ok = d < 0.7
    assert ok.mean() > 0.9
    fr = np.array([S2[k]['flux'] for k in j[ok]]) / np.array([tr_['flux'] for tr_, o in zip(truth, ok) if o]) - 1
    assert abs(np.median(fr)) < 0.03, np.median(fr)


def test_aperture_correction_model_vs_stars(run):
    res, _, _ = run
    ac = res['apcorr']
    assert ac['n'] >= 30
    for r, c, mc in zip(ac['radii'], ac['corr'], ac['model_corr']):
        if r >= 8:
            assert abs(c - mc) < 0.02, (r, c, mc)
    k = ac['radii'].index(3.0)
    assert 0.2 < ac['corr'][k] < 0.7                       # the 3 px aperture misses a third of the light of FWHM ~3.5 px stars


# ------------------------------------------------------------------------------------------------------------- DAOMATCH
def test_daomatch_transform_and_master():
    rng = np.random.default_rng(8)
    n = 120
    x = rng.uniform(20, 580, n); y = rng.uniform(20, 580, n); m = rng.uniform(16, 22, n)
    th, sc = math.radians(0.8), 1.01
    c, s = math.cos(th) * sc, math.sin(th) * sc
    x2 = c * (x - 300) - s * (y - 300) + 300 + 13.4 + rng.normal(0, 0.05, n)
    y2 = s * (x - 300) + c * (y - 300) + 300 - 7.2 + rng.normal(0, 0.05, n)
    m2 = m + 0.5 + rng.normal(0, 0.03, n)
    keep = rng.random(n) > 0.25                              # 25 % of the stars missing in frame 2
    extra = 15                                               # + spurious stars in frame 2
    A = dict(x=x, y=y, mag=m)
    B = dict(x=np.r_[x2[keep], rng.uniform(20, 580, extra)], y=np.r_[y2[keep], rng.uniform(20, 580, extra)], mag=np.r_[m2[keep], rng.uniform(18, 22, extra)])
    r = dp.match_lists(A, B, model='affine', radius=2.0)
    assert r['ok'] and r['n_match'] >= keep.sum() - 3
    assert r['rms'] < 0.12, r['rms']
    q = dp.apply_transform(r['transform'], np.c_[x[keep], y[keep]])
    assert np.sqrt(np.mean((q[:, 0] - x2[keep]) ** 2 + (q[:, 1] - y2[keep]) ** 2)) < 0.12
    assert abs(r['dmag_median'] - 0.5) < 0.02
    ms = dp.master_list([A, B], ref=1, model='affine')
    nfr = np.array([row['nframes'] for row in ms['master']])
    assert (nfr == 2).sum() >= keep.sum() - 3
    two = [row for row in ms['master'] if row['nframes'] == 2]
    assert abs(np.median([row['sigma'] for row in two]) - 0.03 * math.sqrt(2) / 2 * math.sqrt(2)) < 0.05
    # triangle matching without any guess must also cope with a pure shift and a 90 degree rotation
    B2 = dict(x=-(y - 300) + 300 + 5.0, y=(x - 300) + 300 - 3.0, mag=m)
    r2 = dp.match_lists(A, B2, model='similarity', radius=1.5)
    assert r2['ok'] and r2['n_match'] >= n - 3 and r2['rms'] < 0.05


# ------------------------------------------------------------------------------------------------------------ artificial stars
def test_artificial_star_test_reuses_completeness(field, run, tmp_path):
    import importlib.util
    img, t = field
    res, _, model = run
    spec = importlib.util.spec_from_file_location('ogf_completeness', os.path.join(ROOT, 'plugins', 'completeness', 'completeness.py'))
    cm = importlib.util.module_from_spec(spec)
    sys.modules['ogf_completeness'] = cm
    spec.loader.exec_module(cm)
    det = dp.DaophotDetector(model.to_dict(), thresh=4.0, n_iter=2, zp=25.0, fwhm=3.3)
    out = cm.run_completeness(img, det, kind='star', mag_min=19.0, mag_max=22.0, n_bins=4, per_bin=12, per_image=8, zp=25.0, psf=model, match_radius=1.5,
                              seed=3, avoid_detected=False, min_sep=8.0, n_workers=4)
    fr = [b['frac'] for b in out['bins']]
    assert fr[0] >= 0.9 and fr[-1] <= 0.5 and fr[0] > fr[-1]
    b0 = out['bins'][0]
    assert abs(b0['bias']) < 0.05 and b0['scatter'] < 0.08
    assert 19.5 < out['lim50'] < 21.5, out['lim50']


# ------------------------------------------------------------------------------------------------------------------- CLI
def test_cli_pipeline_and_catalog_columns(field, tmp_path):
    img, t = field
    f = tmp_path / 'f.fits'
    fits.PrimaryHDU(img).writeto(f)
    cat = tmp_path / 'cat.tsv'
    sel = np.where(t['mag'] < 19)[0][:20]
    with open(cat, 'w') as fh:
        fh.write('NUMBER\tX_IMAGE\tY_IMAGE\n')
        for k, i in enumerate(sel):
            fh.write('%d\t%.3f\t%.3f\n' % (k + 1, t['x'][i] + 1, t['y'][i] + 1))
    work = tmp_path / 'w'
    p = subprocess.run([sys.executable, CLI, str(f), '--mode', 'run', '--work', str(work), '--catalog', str(cat), '--fwhm', '3.3', '--n-workers', '4'],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    lines = p.stdout.strip().split('\n')
    hdr = lines[0].split('\t')
    assert hdr[0] == 'NUMBER' and 'DAO_MAG' in hdr and 'DAO_FLAG' in hdr and len(lines) == 21
    rows = [dict(zip(hdr, l.split('\t'))) for l in lines[1:]]
    nm = sum(1 for r in rows if r['DAO_MAG'])
    assert nm >= 18
    for fn in ('stars.tsv', 'psf.json', 'psf.fits', 'resid.fits', 'result.json', 'diag.png'):
        assert (work / ('daophot_' + fn)).exists(), fn
    # step-by-step modes
    for mode in ('substar', 'apcorr'):
        q = subprocess.run([sys.executable, CLI, str(f), '--mode', mode, '--work', str(work)], capture_output=True, text=True)
        assert q.returncode == 0, q.stderr[-400:]
    assert (work / 'daophot_sub.fits').exists()
    # psf-add/remove strings are accepted in 1-based coordinates
    q = subprocess.run([sys.executable, CLI, str(f), '--mode', 'pickpsf', '--work', str(work), '--fwhm', '3.3', '--psf-remove', '%.1f,%.1f' % (t['x'][sel[0]] + 1, t['y'][sel[0]] + 1)],
                       capture_output=True, text=True)
    assert q.returncode == 0 and 'PICKPSF' in q.stdout
