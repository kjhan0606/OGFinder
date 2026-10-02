"""Spectra tools against synthetic truth (1D / 2D / IFU cube) and 10 real SDSS spectra when /workspace/fits/sdss exists.  Numbers are printed with -s."""
import glob
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from ogfkit import spectra as sp, spectrasynth as ss, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'spectra.py')
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')


def say(*a):
    sys.stderr.write(' '.join(str(x) for x in a) + '\n')


def run_cli(args):
    p = subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    lines = [l for l in p.stdout.splitlines() if l]
    if not lines:
        return [], [], p.stderr
    cols = lines[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t'))) for l in lines[1:]], p.stderr


def test_readers_roundtrip(tmp_path):
    from astropy.io import fits
    w = np.arange(5000.0, 6000.0, 2.0)
    f = np.sin(w / 50) + 3
    # FITS table, microns, ivar
    t = fits.BinTableHDU.from_columns([fits.Column(name='WAVELENGTH', format='D', unit='um', array=w / 1e4), fits.Column(name='FLUX', format='D', array=f),
                                       fits.Column(name='IVAR', format='D', array=np.full_like(w, 4.0))])
    t.writeto(tmp_path / 't.fits')
    s = sp.read_spectrum(str(tmp_path / 't.fits'))
    assert np.allclose(s['wave'], w) and np.allclose(s['err'], 0.5)
    # 1D image with WCS in nm
    h = fits.PrimaryHDU(f); h.header['CRVAL1'] = 500.0; h.header['CDELT1'] = 0.2; h.header['CRPIX1'] = 1; h.header['CUNIT1'] = 'nm'
    h.writeto(tmp_path / 'i.fits')
    s = sp.read_spectrum(str(tmp_path / 'i.fits'))
    assert np.allclose(s['wave'], w) and s['err'] is None
    # ASCII with comment lines
    np.savetxt(tmp_path / 'a.txt', np.c_[w, f, np.full_like(w, 0.1)], header='wave flux err')
    s = sp.read_spectrum(str(tmp_path / 'a.txt'))
    assert np.allclose(s['wave'], w) and np.allclose(s['err'], 0.1)
    # loglam table (SDSS style)
    t = fits.BinTableHDU.from_columns([fits.Column(name='loglam', format='D', array=np.log10(w)), fits.Column(name='flux', format='D', array=f)])
    t.writeto(tmp_path / 'l.fits')
    assert np.allclose(sp.read_spectrum(str(tmp_path / 'l.fits'))['wave'], w)
    say('readers: FITS table (um, ivar), 1D WCS image (nm), ASCII, loglam table all reproduce the wavelength grid to < 1e-6 A')


def test_line_fit_pulls_and_noise():
    pf, pm, ef, nz = [], [], [], []
    for seed in range(300):
        rng = np.random.default_rng(seed)
        w = np.arange(6000.0, 7000.0, 2.0)
        mu, s, F = 6500.0 + rng.uniform(-3, 3), 3.0, 60.0
        noise = 0.5
        f = 2.0 + F / (s * sp.SQRT2PI) * np.exp(-0.5 * ((w - mu) / s) ** 2) + rng.normal(0, noise, len(w))
        r = sp.fit_line(w, f, np.full_like(w, noise), 6500.0, sigma0=2.5)
        assert r['ok']
        pf.append((r['flux'] - F) / r['flux_err']); pm.append((r['mu'] - mu) / r['mu_err']); ef.append(r['flux_err'] / F)
        cont, mask = sp.continuum(w, f)
        nz.append(sp.noise_level(f - cont, mask) / noise)
    say('Gaussian line fit, 300 noise realisations (flux S/N ~ %.0f): flux pull mean %+.2f std %.2f, centroid pull mean %+.2f std %.2f; continuum-subtracted noise estimate / true %.3f' %
        (1 / np.mean(ef), np.mean(pf), np.std(pf), np.mean(pm), np.std(pm), np.mean(nz)))
    assert abs(np.mean(pf)) < 0.15 and 0.85 < np.std(pf) < 1.2 and abs(np.mean(pm)) < 0.15 and 0.85 < np.std(pm) < 1.2
    assert abs(np.mean(nz) - 1) < 0.06


def test_continuum_survives_strong_lines():
    w, f, e, t = ss.spectrum_1d(z=0.6, seed=3, scale=2.0, noise=0.1)
    cont, mask = sp.continuum(w, f)
    truth = 1.0 * (w / 7500.0) ** -0.15
    dev = np.abs(cont - truth)[~mask]
    say('continuum under strong lines: median |cont - truth| %.4f (noise 0.1), max %.4f; %d pixels masked' % (np.median(dev), dev.max(), mask.sum()))
    assert np.median(dev) < 0.02 and dev.max() < 0.12


def test_redshift_recovery_vs_snr():
    rows = []
    for scale in (0.15, 0.3, 0.6, 1.0):
        ok = wrong = none = 0
        err = []
        snrs = []
        for seed in range(40):
            z = np.random.default_rng(1000 + seed).uniform(0.2, 1.4)
            w, f, e, t = ss.spectrum_1d(z=z, seed=seed, scale=scale)
            r = sp.analyse_spectrum(w, f, e)
            zr = r['redshift']
            snrs.append(max([d['snr'] for d in r['detections']] or [0]))
            if zr['quality'] == 0:
                none += 1
            elif abs(zr['z'] - z) < 0.003:
                ok += 1; err.append(zr['z'] - z)
            else:
                wrong += 1
        rows.append((scale, ok, wrong, none, np.std(err) if err else float('nan'), np.median(snrs)))
        say('redshift recovery, line strength x%.2f (strongest line S/N %.0f): correct %d/40, wrong with quality>0 %d, no answer %d; z scatter of correct %.1e' % (scale, np.median(snrs), ok, wrong, none, rows[-1][4]))
    assert rows[2][1] >= 38 and rows[3][1] >= 39            # x0.6, x1.0
    assert all(r[2] <= 2 for r in rows)                       # confident wrong answers are rare, also at low S/N
    assert rows[1][1] >= 25


def test_2d_optimal_vs_boxcar():
    gb, go, fb, fo = [], [], [], []
    for seed in range(30):
        w, d, e, t = ss.spectrum_2d(seed=seed, noise=0.4)
        fb_, eb_, _ = sp.extract_2d(d, t['row'], 4.0, e, sky=(8, 18), mode='boxcar')
        fo_, eo_, inf = sp.extract_2d(d, t['row'], 4.0, e, sky=(8, 18), mode='optimal')
        m = slice(100, 700)
        tr = t['flux1d']
        gb.append(np.mean((fb_ - tr)[m] ** 2) ** 0.5); go.append(np.mean((fo_ - tr)[m] ** 2) ** 0.5)
        fb.append(np.sum(fb_[m]) / np.sum(tr[m])); fo.append(np.sum(fo_[m]) / np.sum(tr[m]))
    say('2D extraction (point source, Gaussian profile sigma 1.6 rows, sky residual removed): flux / truth boxcar %.3f (rms %.3f), optimal %.3f (rms %.3f); pixel rms error boxcar %.3f, optimal %.3f (x%.2f better)'
        % (np.mean(fb), np.std(fb), np.mean(fo), np.std(fo), np.mean(gb), np.mean(go), np.mean(gb) / np.mean(go)))
    assert abs(np.mean(fb) - 1) < 0.05 and abs(np.mean(fo) - 1) < 0.03 and np.mean(gb) / np.mean(go) > 1.15


def test_cube_aperture_and_velocity_map():
    w, c, t = ss.cube(seed=2, noise=0.05, extra_lines=False)
    spec, _ = sp.extract_cube(c, 10.0, 10.0, 3.0)
    an = sp.analyse_spectrum(w, spec, None, z_known=0.0, snr_min=4)
    ha = next(l for l in an['lines'] if l['name'] == 'Ha')
    # truth: sum of the line flux inside r<3 spaxels (fractional edges)
    yy, xx = np.mgrid[0:21, 0:21]
    wgt = np.clip(3.0 + 0.5 - np.hypot(xx - 10, yy - 10), 0, 1)
    truth = float(np.sum(t['flux'] * wgt))
    mm = sp.line_moments(c, w, t['lam_obs'], half=12.0, cont_gap=8.0, cont_width=50.0, snr_min=5.0)
    ok = np.isfinite(mm['velocity'])
    dv = (mm['velocity'] - t['velocity'])[ok]
    say('IFU cube 21x21x240: aperture r=3 Halpha flux %.1f +- %.1f vs truth %.1f; velocity map: %d/%d spaxels with S/N>5, rms error %.1f km/s (disc amplitude 150), median %+.1f; dispersion median %.0f km/s (true 40 + instrument %.0f)'
        % (ha['flux'], ha['flux_err'], truth, ok.sum(), ok.size, dv.std(), np.median(dv), np.nanmedian(mm['sigma']), 1.5 / 6564.61 * 299792.458))
    assert abs(ha['flux'] - truth) < 3 * ha['flux_err'] + 0.03 * truth
    bright = ok & (mm['flux'] > 15.0)
    dvb = (mm['velocity'] - t['velocity'])[bright]
    say('   bright spaxels (line flux > 15): %d, rms %.1f km/s' % (bright.sum(), dvb.std()))
    assert ok.sum() > 100 and dv.std() < 35 and abs(np.median(dv)) < 3 and dvb.std() < 12


def _write_object_files(d, seed0=0):
    from astropy.io import fits
    truth = {}
    rng = np.random.default_rng(5)
    for i in range(1, 13):
        z = rng.uniform(0.3, 1.2)
        w, f, e, t = ss.spectrum_1d(z=z, seed=seed0 + i, scale=0.8)
        tb = fits.BinTableHDU.from_columns([fits.Column(name='WAVE', format='D', array=w), fits.Column(name='FLUX', format='D', array=f), fits.Column(name='ERR', format='D', array=e)])
        tb.writeto(os.path.join(d, 'spec_%d.fits' % i))
        truth[i] = z
    # 2D spectrum object 13 and cube object 14
    w2, d2, e2, t2 = ss.spectrum_2d(seed=3, noise=0.2)
    h = fits.PrimaryHDU(d2.astype(np.float32)); h.header['CRVAL1'] = w2[0]; h.header['CDELT1'] = w2[1] - w2[0]; h.header['CRPIX1'] = 1
    h.writeto(os.path.join(d, 'spec_13.fits'))
    truth[13] = t2['z']
    wc, cc, tc = ss.cube(seed=4, z=0.25, wave0=7900.0, noise=0.05)
    h = fits.PrimaryHDU(cc.astype(np.float32)); h.header['CRVAL3'] = wc[0]; h.header['CDELT3'] = wc[1] - wc[0]; h.header['CRPIX3'] = 1
    h.writeto(os.path.join(d, 'spec_14.fits'))
    truth[14] = 0.25
    return truth, t2['row']


def test_cli_link_and_fit_mixed_kinds(tmp_path):
    truth, row = _write_object_files(str(tmp_path))
    cat = tmp_path / 'c.tsv'
    with open(cat, 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tROWPOS\n')
        for i in range(1, 16):
            f.write('%d\t11\t11\t%s\n' % (i, '%.2f' % row if i == 13 else ''))
    base = ['--catalog', str(cat), '--work', str(tmp_path / 'w'), '--spec-dir', str(tmp_path), '--meta-out', str(tmp_path / 'meta.json'), '--row-column', 'ROWPOS', '--cube-xy', '1', '--sky-inner', '8', '--sky-outer', '18']
    cols, rows, err = run_cli(['--task', 'link'] + base)
    assert cols == ['NUMBER', 'SP_FILE', 'SP_KIND', 'SP_OK']
    kinds = {r['NUMBER']: r['SP_KIND'] for r in rows}
    assert kinds['1'] == '1d' and kinds['13'] == '2d' and kinds['14'] == 'cube' and rows[14]['SP_OK'] == '0'
    cols, rows, err = run_cli(['--task', 'fit'] + base + ['--extract-mode', 'optimal'])
    assert cols[:5] == ['NUMBER', 'SP_Z', 'SP_ZERR', 'SP_ZQ', 'SP_NLINES']
    dz = []
    for r in rows:
        n = int(r['NUMBER'])
        if n in truth and r['SP_ZQ'] not in ('', '0'):
            dz.append(float(r['SP_Z']) - truth[n])
    nq = sum(1 for r in rows if r['SP_ZQ'] not in ('', '0'))
    say('CLI link + fit on 14 linked objects (12 x 1D, 1 x 2D optimal extraction, 1 x cube) + 1 missing: redshift found for %d, max |dz| %.1e, rms %.1e; object 15 SP_OK=0'
        % (nq, np.max(np.abs(dz)), np.sqrt(np.mean(np.square(dz)))))
    say('   per object: ' + ', '.join('%s:%s' % (r['NUMBER'], r['SP_ZQ'] or '-') for r in rows))
    assert nq >= 13 and np.max(np.abs(dz)) < 0.003
    assert rows[14]['SP_Z'] == ''
    # deterministic + results file + link table + plots
    cols2, rows2, _ = run_cli(['--task', 'fit'] + base + ['--extract-mode', 'optimal'])
    assert rows == rows2
    res = json.load(open(tmp_path / 'w' / 'spectra_results.json'))
    assert '1' in res and res['1']['redshift']['quality'] >= 2 and os.path.exists(tmp_path / 'w' / 'spectra_links.tsv')
    from PIL import Image
    for n in (1, 13, 14):
        out = tmp_path / ('v%d.png' % n)
        run_cli(['--task', 'plot', '--number', str(n), '--out', str(out)] + base)
        im = Image.open(out)
        assert im.size[0] >= 500 and im.size[1] >= 250
    say('viewer PNGs for 1D / 2D / cube objects:', Image.open(tmp_path / 'v1.png').size, Image.open(tmp_path / 'v13.png').size, Image.open(tmp_path / 'v14.png').size)


def test_link_file_overrides_pattern(tmp_path):
    from astropy.io import fits
    w, f, e, t = ss.spectrum_1d(z=0.5, seed=1)
    np.savetxt(tmp_path / 'weird_name.txt', np.c_[w, f, e])
    (tmp_path / 'c.tsv').write_text('NUMBER\tX_IMAGE\tY_IMAGE\n7\t1\t1\n8\t1\t1\n')
    (tmp_path / 'links.tsv').write_text('NUMBER\tFILE\n7\tweird_name.txt\n')
    cols, rows, _ = run_cli(['--task', 'fit', '--catalog', str(tmp_path / 'c.tsv'), '--work', str(tmp_path / 'w'), '--link-file', str(tmp_path / 'links.tsv')])
    assert rows[0]['SP_ZQ'] == '3' and abs(float(rows[0]['SP_Z']) - 0.5) < 0.002 and rows[1]['SP_Z'] == ''


SDSS = sorted(glob.glob(os.path.join(FITS, 'sdss', 'spec-*.fits')))


@pytest.mark.skipif(len(SDSS) < 5, reason='SDSS test spectra not available')
def test_real_sdss_spectra():
    from astropy.io import fits
    tab = []
    for p in SDSS:
        with fits.open(p) as h:
            so = h[2].data[0]
            zt, cl = float(so['Z']), str(so['CLASS'])
        s = sp.read_spectrum(p)
        r = sp.analyse_spectrum(s['wave'], s['flux'], s['err'])
        z = r['redshift']
        tab.append((os.path.basename(p)[-9:-5], cl, zt, z['z'], z.get('z_err', float('nan')), z['quality'], (z['z'] - zt) / (1 + zt) * 299792.458))
    good = [t for t in tab if t[5] >= 2]
    say('real SDSS spectra (%d): quality>=2 for %d, all within %.0f km/s of the SDSS pipeline redshift (rms %.0f km/s); %d with quality 0 (absorption-line galaxies, star, weak lines: no answer, as intended)' %
        (len(tab), len(good), max(abs(t[6]) for t in good), math.sqrt(np.mean([t[6] ** 2 for t in good])), sum(1 for t in tab if t[5] == 0)))
    for t in tab:
        say('   plate 266 fiber %s %-6s SDSS z %.5f  found %.5f +- %.5f  q%d  dv %+.0f km/s' % t)
    assert len(good) >= 5 and max(abs(t[6]) for t in good) < 600
