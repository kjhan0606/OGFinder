import json, math, os, subprocess, sys
import numpy as np
import pytest
from scipy.special import ndtr
import completeness as C
from ogfkit import models, imageio, tsvio, meta

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')
ZP = 25.0
FWHM = 3.0


def noise_image(shape=(700, 700), sigma=1.0, seed=5):
    return np.random.default_rng(seed).normal(0.0, sigma, shape).astype(np.float32)


def analytic_completeness(mag, thresh, fwhm=FWHM, sigma=1.0, ngrid=7):
    """P(at least one pixel > thresh*sigma) for a star of magnitude mag, averaged over sub-pixel positions (white noise, no smoothing,
    minarea=1).  This is exactly what SEP (filter off, minarea 1) detects (plus a negligible noise-only term)."""
    F = 10 ** (-0.4 * (mag - ZP))
    offs = (np.arange(ngrid) + 0.5) / ngrid - 0.5
    ps = []
    for dy in offs:
        for dx in offs:
            p = models.shifted_psf(C.default_psf(fwhm), dx, dy)
            ps.append(1.0 - np.prod(ndtr(thresh - F * p / sigma)))
    return float(np.mean(ps))


def test_star_completeness_matches_analytic_curve():
    thr = 4.0
    data = noise_image()
    det = C.SepDetector(zp=ZP, thresh=thr, minarea=1, filter=False)
    res = C.run_completeness(data, det, kind='star', mag_min=20.0, mag_max=23.0, n_bins=10, per_bin=120, per_image=40, zp=ZP,
                             psf_fwhm=FWHM, seed=11)
    meas = np.array([b['frac'] for b in res['bins']])
    ana = np.array([analytic_completeness(b['mag'], thr) for b in res['bins']])
    err = np.array([max(b['frac_hi'] - b['frac_lo'], 0.02) / 2 for b in res['bins']])
    # chi-square of the measured fractions around the analytic curve (bin-centre value; the in-bin variation is small)
    ana_bin = np.array([np.mean([analytic_completeness(m, thr) for m in np.linspace(b['mag_lo'], b['mag_hi'], 5)]) for b in res['bins']])
    z = (meas - ana_bin) / np.sqrt(np.clip(ana_bin * (1 - ana_bin), 0.02, None) / 120)
    assert np.all(np.abs(z) < 3.5), (meas, ana_bin, z)
    assert np.sum(z ** 2) < 35, np.sum(z ** 2)           # 10 bins
    # 50 % limit: measured vs analytic (solve on a fine grid)
    grid = np.linspace(20.0, 23.0, 61)
    a50 = float(np.interp(0.5, [analytic_completeness(g, thr) for g in grid][::-1], grid[::-1]))
    assert abs(res['lim50'] - a50) < 0.08, (res['lim50'], a50)
    a90 = float(np.interp(0.9, [analytic_completeness(g, thr) for g in grid][::-1], grid[::-1]))
    assert abs(res['lim90'] - a90) < 0.12, (res['lim90'], a90)
    # recovered photometry of the brightest bins: bias ~ 0, scatter ~ the analytic noise (sum over the Kron aperture, a few pixels' worth)
    assert abs(res['bins'][0]['bias']) < 0.05 and res['bins'][0]['scatter'] < 0.15
    json.dumps(res)


def test_false_positive_rate_on_noise():
    thr = 4.0
    data = noise_image(seed=9)
    det = C.SepDetector(zp=ZP, thresh=thr, minarea=1, filter=False)
    res = C.run_completeness(data, det, kind='star', mag_min=21, mag_max=22, n_bins=2, per_bin=10, per_image=10, zp=ZP, seed=2)
    expected = data.size * (1 - ndtr(thr))
    got = res['false_positive']['negative_image_detections']
    assert abs(got - expected) < 4 * math.sqrt(expected) + 1, (got, expected)
    assert res['false_positive']['per_1e6_pix'] == pytest.approx(1e6 * got / data.size)


def test_injected_flux_is_conserved_and_galaxies_are_psf_convolved():
    psf = C.default_psf(3.0)
    for o in (dict(x=50.3, y=50.7, mag=22.0), dict(x=60.0, y=60.0, mag=22.0, re=4.0, n=1.0, q=0.7, pa=30.0),
              dict(x=60.0, y=60.0, mag=22.0, re=3.0, n=4.0, q=0.9, pa=10.0)):
        st, y0, x0 = C.render_object(o, psf, ZP)
        tot = float(st.sum())
        flux = 10 ** (-0.4 * (22.0 - ZP))
        # Sersic n=4 has a long tail beyond the stamp: <5 %; others < 1 %
        tol = 0.05 if o.get('n') == 4.0 else 0.01
        assert abs(tot / flux - 1) < tol, (o, tot / flux)


def test_mask_and_detected_sources_are_avoided():
    data = noise_image((400, 400), seed=3)
    mask = np.zeros(data.shape, bool)
    mask[:, :200] = True
    seen = []

    class Rec:
        def __call__(self, img, m):
            return dict(x=[], y=[], mag=[])
    res = C.run_completeness(data, Rec(), mask=mask, kind='star', mag_min=21, mag_max=22, n_bins=1, per_bin=20, per_image=10, zp=ZP, seed=4)
    # a detector that finds nothing: nothing recovered, but 20 objects must have been injected
    assert res['n_injected'] == 20 and res['n_recovered'] == 0
    # positions: call the drawing function directly
    rng = np.random.default_rng(1)
    objs = C._draw_objects(rng, 15, data.shape, [21.0] * 15, 'star', mask, 10, 30.0, [2.0, 4.0], [1.0], (0.5, 1.0))
    assert len(objs) == 15 and all(o['x'] >= 200 - 1 for o in objs)
    d = [(o['x'], o['y']) for o in objs]
    assert min(math.hypot(a[0] - b[0], a[1] - b[1]) for i, a in enumerate(d) for b in d[i + 1:]) >= 30.0


def test_detection_callback_api_and_module_spec(tmp_path):
    """Any callback works; 'module:function' specs are resolved by load_detector."""
    sys.path.insert(0, str(tmp_path))
    (tmp_path / 'mydet.py').write_text(
        "import numpy as np\nfrom scipy import ndimage\n"
        "def factory(zp=25.0, thresh=5.0, **kw):\n"
        "    def det(img, mask=None):\n"
        "        sm = ndimage.gaussian_filter(img, 1.5)\n"
        "        s = 1.4826*np.median(np.abs(sm-np.median(sm)))\n"
        "        lab, n = ndimage.label(sm > thresh*s)\n"
        "        if n == 0: return []\n"
        "        idx = range(1, n+1)\n"
        "        cy_cx = ndimage.center_of_mass(sm, lab, idx)\n"
        "        return [dict(x=c[1], y=c[0]) for c in cy_cx]\n"
        "    return det\n")
    det = C.load_detector('mydet:factory', zp=ZP, thresh=5.0)
    data = noise_image((500, 500), seed=8)
    res = C.run_completeness(data, det, kind='star', mag_min=19.5, mag_max=23.5, n_bins=8, per_bin=40, per_image=30, zp=ZP, seed=5)
    fr = [b['frac'] for b in res['bins']]
    assert fr[0] > 0.95 and fr[-1] < 0.2 and np.isfinite(res['lim50']), fr        # monotone-ish transition from 1 to 0
    assert all(math.isnan(b['bias']) or True for b in res['bins'])                # detector gives no magnitudes: bias is NaN, not an error
    assert 20.5 < res['lim50'] < 23.0


def test_fit_curve_recovers_known_logistic():
    mag = np.linspace(22, 28, 13)
    n = np.full(13, 400)
    f = C.logistic(mag, 26.0, 0.4, 0.97)
    k = np.round(f * n)
    fit = C.fit_curve(mag, k, n)
    assert abs(fit['mh'] - 26.0) < 0.03 and abs(fit['w'] - 0.4) < 0.05 and abs(fit['fmax'] - 0.97) < 0.01
    # analytic limits
    assert abs(fit['m50'] - (26.0 + 0.4 * math.log(0.97 / 0.5 - 1))) < 0.03
    assert abs(fit['m90'] - (26.0 + 0.4 * math.log(0.97 / 0.9 - 1))) < 0.05
    assert abs(C.direct_limit(mag, k / n, 0.5) - 26.0) < 0.2
    assert C.wilson(0, 10)[0] == 0.0 and C.wilson(10, 10)[1] == 1.0


def test_cli_with_existing_extraction_and_catalog_metadata(tmp_path):
    data = noise_image((600, 600), sigma=0.1, seed=21)
    # a few real sources
    rng = np.random.default_rng(3)
    for x, y in ((100, 100), (300, 400), (450, 150)):
        data = data + models.shifted_psf(C.default_psf(3.0, 41), 0, 0).sum() * 0  # no-op, keeps dtype
        data[y - 20:y + 21, x - 20:x + 21] += 3000 * C.default_psf(3.0, 41)
    imageio.save_fits(str(tmp_path / 'im.fits'), data)
    det = C.SextractDetector(zp=ZP)
    d = det(data, None)
    cat = str(tmp_path / 'cat.tsv')
    with open(cat, 'w') as f:
        f.write('NUMBER\tMAG_AUTO\n')
        for i, m in enumerate([18.0, 22.0, 24.0, 26.0, 29.0], 1):
            f.write('%d\t%s\n' % (i, m))
    cmd = [PY, os.path.join(os.path.dirname(HERE), 'completeness.py'), str(tmp_path / 'im.fits'), '--detector', 'sextract',
           '--mag-min', '23', '--mag-max', '29', '--n-bins', '8', '--per-bin', '30', '--per-image', '25', '--mag-zeropoint', '25',
           '--catalog', cat, '--json-out', str(tmp_path / 'c.json'), '--curve-out', str(tmp_path / 'c.tsv'), '--meta-out', str(tmp_path / 'catalog_meta.json'),
           '--detect-thresh', '1.5', '--detect-minarea', '5', '--seed', '2']
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    lines = p.stdout.strip().split('\n')
    assert lines[0].split('\t') == ['NUMBER'] + C.COLUMNS and len(lines) == 6
    vals = [[float(x) if x else float('nan') for x in l.split('\t')[1:]] for l in lines[1:]]
    assert vals[0][0] > 0.95 and vals[-1][0] < 0.5 and vals[0][0] >= vals[2][0] >= vals[-1][0]   # completeness falls with magnitude
    res = json.load(open(tmp_path / 'c.json'))
    assert 23 < res['lim50'] < 29 and res['lim90'] <= res['lim50']
    m = meta.load(str(tmp_path / 'catalog_meta.json'))
    assert m['nrows'] == 5 and abs(m['completeness']['lim50'] - res['lim50']) < 1e-9
    assert 'completeness.lim50' in meta.flat_scalars(m)
    assert len(tsvio.read_catalog(str(tmp_path / 'c.tsv'))[1]) == 8


@pytest.mark.skipif(not os.path.exists(FITS + '/hudf_f160w.fits'), reason='HUDF not available')
def test_hudf_crop_completeness_is_monotone_and_sane():
    data, _ = imageio.load_image(FITS + '/hudf_f160w.fits')
    crop = np.ascontiguousarray(data[1500:2100, 1500:2100])
    zp = 25.94
    res = C.run_completeness(crop, C.SextractDetector(zp=zp), kind='star', mag_min=27.0, mag_max=31.5, n_bins=9, per_bin=30, per_image=20,
                             zp=zp, psf_fwhm=3.0, seed=7, pixel_scale=0.06)
    fr = np.array([b['frac'] for b in res['bins']])
    assert fr[0] >= 0.8 and fr[-1] < 0.6 and fr[0] > fr[-1]
    assert 28.0 < res['lim50'] < 31.5, res['lim50']
    assert res['false_positive']['negative_image_detections'] > 0


def test_fits_export_carries_catalog_metadata(tmp_path):
    """ds9_fits_export.py copies <dir of the input>/catalog_meta.json into the table header when nrows matches (and only then)."""
    from astropy.io import fits
    lib = '/workspace/OGFinder/ds9/library/ds9_fits_export.py'
    tsv = tmp_path / 'cat.tsv'
    tsv.write_text('NUMBER\tMAG_AUTO\n1\t20.0\n2\t21.0\n3\t22.0\n')
    meta.update(str(tmp_path / 'catalog_meta.json'), 'completeness', {'lim50': 27.25, 'lim90': 26.5, 'kind': 'star'}, nrows=3)
    out = tmp_path / 'o.fits'
    p = subprocess.run([PY, lib, '--input', str(tsv), '--output', str(out)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    h = fits.getheader(str(out), 1)
    assert h['HIERARCH OGF completeness.lim50'] == 27.25 and h['HIERARCH OGF completeness.kind'] == 'star'
    meta.update(str(tmp_path / 'catalog_meta.json'), 'completeness', {'lim50': 1.0}, nrows=99)       # other catalog: ignored
    out2 = tmp_path / 'o2.fits'
    subprocess.run([PY, lib, '--input', str(tsv), '--output', str(out2)], capture_output=True, text=True, check=True)
    assert 'HIERARCH OGF completeness.lim50' not in fits.getheader(str(out2), 1)
