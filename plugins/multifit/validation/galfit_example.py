#!/usr/bin/env python3
"""GALFIT's own distribution example (galfit-example/EXAMPLE: gal.fits, psf.fits, galfit.feedme; HST/NICMOS-like galaxy) fitted by GALFIT and by multifit from the same feedme.
    python galfit_example.py EXAMPLE_DIR WORKDIR [--galfit BIN] [--ld DIR]
GALFIT builds its own sigma image from the header (GAIN, NCOMBINE, ...) - a different noise model - so, for a like-for-like comparison, both programs are run a second time
with the SAME constant sigma image (robust background sigma of the data), and chi^2 of every result is evaluated with GALFIT's own renderer (P=1)."""
import os
import re
import shutil
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
from astropy.io import fits  # noqa: E402
from ogfkit import galfitio as GI, imageio  # noqa: E402


def run(galfit, ld, cwd, feed):
    env = dict(os.environ)
    if ld:
        env['LD_LIBRARY_PATH'] = ld + ':' + env.get('LD_LIBRARY_PATH', '')
    subprocess.run([galfit, feed], cwd=cwd, env=env, capture_output=True, stdin=subprocess.DEVNULL, timeout=300)


def last_result(d):
    fs = sorted(f for f in os.listdir(d) if re.match(r'galfit\.\d+$', f))
    return open(os.path.join(d, fs[-1])).read() if fs else None


def chi2_const(galfit, ld, wd, res_text, sigma):
    ed = os.path.join(wd, 'eval')
    os.makedirs(ed, exist_ok=True)
    for f in ('gal.fits', 'psf.fits'):
        shutil.copy(os.path.join(wd, f), ed)
    t = re.sub(r'(?m)^P\).*$', 'P) 1', res_text)
    for k, v in (('A', 'gal.fits'), ('B', 'o.fits'), ('C', 'none'), ('D', 'psf.fits'), ('G', 'none')):
        t = re.sub(r'(?m)^%s\).*$' % k, '%s) %s' % (k, v), t)
    open(os.path.join(ed, 'e.feedme'), 'w').write(t)
    run(galfit, ld, ed, 'e.feedme')
    cfg = GI.parse_feedme(t)
    x0, x1, y0, y1 = cfg['region']
    d = fits.getdata(os.path.join(ed, 'gal.fits')).astype(float)[y0 - 1:y1, x0 - 1:x1]
    m = fits.getdata(os.path.join(ed, 'o.fits')).astype(float)[y0 - 1:y1, x0 - 1:x1]
    return float(np.sum(((d - m) / sigma) ** 2)), d.size


def show(tag, res_text):
    cfg = GI.parse_feedme(res_text, strict=False)
    c = cfg['components'][0]
    print('  %-34s mag %.3f  re %.3f  n %.3f  q %.3f  PA %.2f  x,y %.3f %.3f  sky %.3f' % (tag, c['mag'] + 0, c['re'], c['n'], c['q'], c['pa'] - 90, c['x'], c['y'], cfg['sky_value_galfit']))


def main():
    ex, wd = sys.argv[1], sys.argv[2]
    galfit = sys.argv[sys.argv.index('--galfit') + 1] if '--galfit' in sys.argv else os.environ.get('GALFIT_BIN', 'galfit')
    ld = sys.argv[sys.argv.index('--ld') + 1] if '--ld' in sys.argv else os.environ.get('GALFIT_LD', '')
    shutil.rmtree(wd, ignore_errors=True)
    shutil.copytree(ex, wd)
    data = fits.getdata(os.path.join(wd, 'gal.fits')).astype(float)
    hdr = fits.getheader(os.path.join(wd, 'gal.fits'))
    ex_t = float(hdr.get('EXPTIME', 1.0))
    print('data %s, EXPTIME %.0f; magnitudes below are GALFIT magnitudes (exposure time included)' % (data.shape, ex_t))
    # 1. GALFIT with its own sigma
    run(galfit, ld, wd, 'galfit.feedme')
    g_own = last_result(wd)
    show('GALFIT (own Poisson sigma)', g_own)
    # 2. multifit from the same feedme (constant robust sigma = its default)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py'), '-', '--config', os.path.join(wd, 'galfit.feedme'), '--work', os.path.join(wd, 'mf'),
                        '--export-feedme', os.path.join(wd, 'mf', 'gf')], capture_output=True, text=True)
    m_own = open(os.path.join(wd, 'mf', 'gf', 'galfit.feedme')).read()
    show('multifit (robust constant sigma)', m_own)
    # 3. like for like: constant sigma image for both
    bg = data[:12, :12]
    sig = float(imageio.robust_sigma(data))
    fits.PrimaryHDU(np.full(data.shape, sig, np.float32)).writeto(os.path.join(wd, 'sigma.fits'))
    feed = open(os.path.join(wd, 'galfit.feedme')).read()
    feed = re.sub(r'(?m)^C\).*$', 'C) sigma.fits', feed)
    open(os.path.join(wd, 'const.feedme'), 'w').write(feed)
    for f in [f for f in os.listdir(wd) if re.match(r'galfit\.\d+$', f)]:
        os.remove(os.path.join(wd, f))
    run(galfit, ld, wd, 'const.feedme')
    g_c = last_result(wd)
    show('GALFIT (constant sigma %.4g)' % sig, g_c)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'plugins', 'multifit', 'multifit.py'), '-', '--config', os.path.join(wd, 'const.feedme'), '--work', os.path.join(wd, 'mf2'),
                        '--export-feedme', os.path.join(wd, 'mf2', 'gf')], capture_output=True, text=True)
    m_c = open(os.path.join(wd, 'mf2', 'gf', 'galfit.feedme')).read()
    show('multifit (constant sigma %.4g)' % sig, m_c)
    for tag, t in (('GALFIT const-sigma fit', g_c), ('multifit const-sigma fit', m_c), ('GALFIT own-sigma fit', g_own), ('multifit own fit', m_own)):
        c2, n = chi2_const(galfit, ld, wd, t, sig)
        print('  chi2 under the constant sigma, %-26s %.2f  (%d pixels, chi2/N %.4f)' % (tag, c2, n, c2 / n))


if __name__ == '__main__':
    main()
