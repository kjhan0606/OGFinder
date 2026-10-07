"""Dither stack (trails.py --task stack) regressions found on the real HST ACS/WFC visit jc8m32010 (GO-13498; Astrafex example
examples/hst_dither_stack): lookup-table WCS, NaN leakage in the resampling, MEF chips / DQ / ERR / EXPTIME / sky handling, and the
cosmic-ray rejection of few exposures."""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLUG = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, PLUG)
from ogfkit import trails as T, trails_ext as X  # noqa: E402
import trails as TR  # noqa: E402
CLI = os.path.join(PLUG, 'trails.py')


def _wcs(crpix, shape=None):
    from astropy.wcs import WCS
    w = WCS(naxis=2); w.wcs.ctype = ['RA---TAN', 'DEC--TAN']; w.wcs.crval = [109.4, 37.75]; w.wcs.crpix = crpix
    w.wcs.cd = [[-0.05 / 3600, 0], [0, 0.05 / 3600]]
    return w


def test_resample_masked_pixel_does_not_pull_neighbours_down():
    img = np.full((40, 40), 100.0); img[20, 20] = np.nan
    w0, w1 = _wcs([20.0, 20.0]), _wcs([20.5, 20.0])          # half-pixel shift
    rs = X.wcs_resample(img, w1, w0, img.shape)
    fin = rs[np.isfinite(rs)]
    assert np.allclose(fin, 100.0)                            # before: the NaN entered as 0 -> neighbours ~50
    assert np.isfinite(rs).sum() >= 40 * 40 - 3 * 40          # only the edge column and the bad pixel are lost


def _lookup_mef(path, data, crpix):
    """an HST-like file whose WCS has a lookup-table distortion (CPDIS / WCSDVARR): astropy needs the HDUList to read it"""
    from astropy.io import fits
    h = _wcs(crpix).to_header()
    h['EXTNAME'] = 'SCI'; h['EXTVER'] = 1; h['EXPTIME'] = 100.0
    for i in (1, 2):
        h['CPDIS%d' % i] = 'LOOKUP'; h['DP%d' % i] = 'EXTVER: %d' % i
        h['DP%d' % i] = 'NAXES: 2'; h['DP%d' % i] = 'AXIS.1: 1'; h['DP%d' % i] = 'AXIS.2: 2'
    hl = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data.astype(np.float32), h)])
    for i in (1, 2):
        t = fits.ImageHDU(np.full((4, 4), 0.02, np.float32), name='WCSDVARR', ver=i)
        for k, v in (('CRPIX1', 0.0), ('CRPIX2', 0.0), ('CRVAL1', 0.0), ('CRVAL2', 0.0), ('CDELT1', 20.0), ('CDELT2', 20.0)):
            t.header[k] = v
        hl.append(t)
    hl.writeto(path, overwrite=True)


def test_lookup_table_wcs_registers_and_survives_in_the_stack(tmp_path):
    from astropy.io import fits
    from astropy.wcs import WCS
    n = 60
    yy, xx = np.mgrid[0:n, 0:n]
    ps = []
    for k, dx in enumerate((0.0, 3.4)):
        img = 10 + 500 * np.exp(-0.5 * ((xx - 30 - dx) ** 2 + (yy - 30) ** 2) / 1.2 ** 2)
        p = str(tmp_path / ('f%d.fits' % k)); _lookup_mef(p, img, [30.0 + dx, 30.0]); ps.append(p)
    hdr = fits.getheader(ps[0], 1)
    try:
        WCS(hdr); plain = True
    except ValueError:
        plain = False
    assert not plain                                          # the failure mode: header-only WCS raises
    assert TR._wcs(hdr, fits.open(ps[0])) is not None
    w = str(tmp_path / 'w')
    r = subprocess.run([sys.executable, CLI, ps[0], '--task', 'stack', '--frames'] + ps + ['--work', w, '--stack-method', 'mean'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-400:]
    assert 'registered=1' in r.stdout                         # before: registered=0 (silently unregistered)
    with fits.open(os.path.join(w, 'trails_stack.fits')) as hl:
        assert [h.name for h in hl[1:]] == ['WCSDVARR', 'WCSDVARR']
        WCS(hl[0].header, hl)                                 # the stack's WCS is readable with the tables
        s = hl[0].data.astype(float)
    assert s[30, 30] > 0.9 * s.max() and np.unravel_index(np.argmax(s), s.shape) == (30, 30)


def _noise_frames(nf=4, n=300, seed=1):
    r = np.random.default_rng(seed)
    return [r.normal(0, 1, (n, n)) for _ in range(nf)], [np.ones((n, n))] * nf


def test_crrej_rejects_single_and_double_hits_without_overclipping():
    F, V = _noise_frames()
    F[1][50, 50] += 100; F[2][60, 60] += 10
    F[0][80, 80] += 30; F[3][80, 80] += 12                    # the same pixel hit in 2 of 4 exposures
    F[0][90, 90] += 300; F[2][90, 90] += 320; F[0][90, 91] += 40; F[2][90, 91] += 30   # 2 bright hits (the jc8m32 case: median spike -> D)
    o, n = T.stack_frames(F, method='crrej', variances=V)
    assert n[50, 50] == 3 and n[60, 60] == 3 and n[80, 80] == 2 and n[90, 90] == 2 and n[90, 91] == 2
    assert abs(o[90, 90]) < 3 and abs(o[90, 91]) < 3
    assert abs(o[80, 80]) < 3 and abs(o[50, 50]) < 3
    assert (n < 4).mean() < 0.002                             # false rejections of good noise pixels
    o2, n2 = T.stack_frames(F, method='sigclip')
    assert (n2 < 4).mean() > 0.1                              # documents the MAD clipping of 4 values: ~20 % of good pixels lost
    o3, n3 = T.stack_frames(F[:2], method='crrej', variances=V[:2])   # two exposures: R = minimum
    assert n3[50, 50] == 1 and abs(o3[50, 50]) < 4 and (n3 < 2).mean() < 0.002


def test_crrej_keeps_star_cores_with_registration_differences():
    n = 41
    yy, xx = np.mgrid[0:n, 0:n]
    flux, s = 2e4, 0.9                                        # undersampled star, ~20 % peak differences between exposures
    F, V = [], []
    r = np.random.default_rng(3)
    for dx, dy in ((0, 0), (0.25, -0.2), (-0.2, 0.25), (0.15, 0.15)):
        m = flux / (2 * np.pi * s * s) * np.exp(-0.5 * ((xx - 20 - dx) ** 2 + (yy - 20 - dy) ** 2) / s ** 2)
        var = m + 25.0
        F.append(m + r.normal(0, np.sqrt(var))); V.append(var)
    o, nn = T.stack_frames(F, method='crrej', variances=V)
    assert (nn[17:24, 17:24] == 4).all()                     # nothing rejected in the core
    assert abs(o[12:29, 12:29].sum() / flux - 1) < 0.02


def _hst_like(path, chips, exptime, sky):
    """two-chip MEF (SCI/ERR/DQ per chip) in electrons"""
    from astropy.io import fits
    hl = [fits.PrimaryHDU(header=fits.Header([('EXPTIME', exptime)]))]
    for ver, (data, crpix, dq) in enumerate(chips, 1):
        h = _wcs(crpix).to_header()
        e = np.sqrt(np.maximum(data + sky * exptime, 0) + 25.0)
        hl += [fits.ImageHDU((data + sky * exptime).astype(np.float32), h, name='SCI', ver=ver),
               fits.ImageHDU(e.astype(np.float32), h, name='ERR', ver=ver), fits.ImageHDU(dq.astype(np.int16), h, name='DQ', ver=ver)]
    fits.HDUList(hl).writeto(path, overwrite=True)


def test_mef_chips_dq_exptime_and_sky(tmp_path):
    """frame 1 is dithered by 30 px so its second chip fills part of the reference chip; exposure times and sky levels differ; a DQ-flagged
    hot column in frame 0 must not reach the stack"""
    from astropy.io import fits
    ny, nx = 80, 100
    rate = 2.0                                                # source rate e-/s of a flat patch at sky x = 60..70 of the reference chip
    paths = []
    for k, (shift, texp, sky) in enumerate(((0, 100.0, 1.0), (30, 300.0, 3.0))):
        chips = []
        for c in range(2):
            xoff = c * nx - shift                            # chip c covers reference columns xoff .. xoff + nx - 1
            img = np.zeros((ny, nx)); dq = np.zeros((ny, nx), int)
            xs = np.arange(nx) + xoff
            img[:, (xs >= 60) & (xs < 70)] = rate * texp
            if k == 0 and c == 0:
                img[:, 20] += 5e4; dq[:, 20] = 16 | 2        # hot column flagged in DQ
            chips.append((img, [1.0 - xoff, 1.0], dq))
        p = str(tmp_path / ('h%d.fits' % k)); _hst_like(p, chips, texp, sky); paths.append(p)
    w = str(tmp_path / 'w')
    r = subprocess.run([sys.executable, CLI, paths[0], '--task', 'stack', '--frames'] + paths + ['--work', w, '--stack-method', 'crrej', '--all-chips',
                        '--dq-bits', '2', '--scale', 'exptime', '--sky', 'median', '--max-trails', '0'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-400:]
    s = fits.getdata(os.path.join(w, 'trails_stack.fits')).astype(float)
    nmap = fits.getdata(os.path.join(w, 'trails_stack_n.fits'))
    js = json.load(open(os.path.join(w, 'trails_stack.json')))
    assert js['registered'] == 1 and len(js['frames']) == 4
    assert np.allclose(s[:, 62:68], rate, atol=0.05)         # e-/s, both exposures agree after scaling
    assert np.abs(s[:, 30:50]).max() < 0.05                  # sky removed
    assert abs(s[:, 20]).max() < 0.05 and (nmap[:, 20] == 1).all()   # DQ column excluded, frame 1 alone there
    assert (nmap[:, 75:95] == 2).all()                       # frame 1's second chip fills the reference chip edge region it covers
