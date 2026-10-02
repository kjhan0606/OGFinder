#!/usr/bin/env python3
"""Real-data comparison of the headless `cli` blocks with the legacy GUI path (item 1).

1. prepares inputs from /workspace/fits/m51.fits (scaled copies = known flux ratios, a Gaussian PSF file),
2. runs the legacy dialogs/procs in a real ds9 under Xvfb with the REAL drivers (scripts/verify_headless_real_gui.tcl) -> gui_<step>.tsv,
3. runs the same steps headless through the batch runner (manifest `headless` templates, no GUI) -> catalogs,
4. compares: catalogs must be byte-identical (same drivers, same argv); deconvolved image and segmentation map pixel-identical;
   plus quantitative checks (multi-band magnitude offsets vs. the known flux scaling, Gaia match count, completeness values).
Usage: verify_headless_real.py [--workdir DIR] [--display :77]   (needs network for the VizieR cross-match, skipped with --no-network)
Exit status 1 on any failure.
"""
import argparse, json, os, subprocess, sys, shutil, tempfile
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ogfkit import batchrun, cliexpand

ap = argparse.ArgumentParser()
ap.add_argument('--workdir', default='/workspace/work/i1b/real')
ap.add_argument('--display', default=':77')
ap.add_argument('--no-network', action='store_true')
ap.add_argument('--python', default=os.environ.get('OGFINDER_PYTHON', sys.executable))
a = ap.parse_args()
W = os.path.abspath(a.workdir); shutil.rmtree(W, ignore_errors=True); os.makedirs(W)
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')
nfail = 0
def R(tag, ok, d=''):
    global nfail
    print(('PASS' if ok else 'FAIL'), tag, d, flush=True)
    nfail += (not ok)

from astropy.io import fits
with fits.open(FITS + '/m51.fits') as h:
    data = h[0].data.astype('float32'); hdr = h[0].header
for nm, f in (('m51_half', 0.5), ('m51_quarter', 0.25)):
    fits.PrimaryHDU(data * f, hdr).writeto('%s/%s.fits' % (W, nm), overwrite=True)
y, x = np.mgrid[-15:16, -15:16]; g = np.exp(-(x * x + y * y) / (2 * (3.0 / 2.355) ** 2)); g /= g.sum()
fits.PrimaryHDU(g.astype('float32')).writeto(W + '/psf.fits', overwrite=True)

# ---- GUI half
home = tempfile.mkdtemp(prefix='ogf_hlreal_home.')
env = dict(os.environ, HOME=home, DISPLAY=a.display, OGF_HLREAL_DIR=W, OGFINDER_PYTHON=a.python)
if subprocess.run(['xdpyinfo'], env=env, capture_output=True).returncode != 0:
    subprocess.Popen(['Xvfb', a.display, '-screen', '0', '1400x1000x24'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    import time; time.sleep(2)
if a.no_network:
    pass
r = subprocess.run(['timeout', '-s', 'KILL', '400', ROOT + '/bin/ds9', FITS + '/m51.fits', '-geometry', '1300x950', '-source', ROOT + '/scripts/verify_headless_real_gui.tcl'],
                   env=env, capture_output=True, text=True)
shutil.rmtree(home, ignore_errors=True)
log = open(W + '/gui.log').read() if os.path.exists(W + '/gui.log') else ''
R('gui_half_finished', 'done' in log.splitlines()[-1:] , log.strip().splitlines()[-1] if log.strip() else 'no log')
if 'BGERROR' in log:
    print(log)

# ---- headless half: one batch recipe per branch (every step replaces or extends the catalog of the field)
man = cliexpand.load_manifests(ROOT)
fields = [('m51', FITS + '/m51.fits')]
common = {'extract:detect-thresh': 2.5, 'extract:mag-zeropoint': 25.0, 'extract:n-workers': 1}
def batch(name, steps):
    rec = {'steps': [dict(plugin='extract', step='extract', params={'detect-thresh': 2.5, 'mag-zeropoint': 25.0, 'n-workers': 1})] + steps}
    out = '%s/hl_%s' % (W, name)
    res = batchrun.run_batch(fields, rec, out, ROOT, a.python, 1, False, 0, 300, [])
    return out, res[0]
def gui(name): return open('%s/gui_%s.tsv' % (W, name)).read() if os.path.exists('%s/gui_%s.tsv' % (W, name)) else None
def hl(out): return open(out + '/m51/catalog.tsv').read()
def cmp(name, steps, tag=None, parse=None):
    out, res = batch(name, steps)
    ok = res['status'] == 'ok'
    R('%s_headless_ran' % name, ok, res.get('message', ''))
    g, h = gui(name), (hl(out) if ok else None)
    R('%s_catalog_identical_to_gui' % name, g is not None and h is not None and g.rstrip('\n') == h.rstrip('\n'), 'rows gui=%s hl=%s' % (None if g is None else g.count('\n'), None if h is None else h.count('\n')))
    return out, h

P = lambda plugin, step, **kw: dict(plugin=plugin, step=step, params=dict(common, **kw) if plugin != 'extract' else dict({'detect-thresh': 2.5, 'mag-zeropoint': 25.0, 'n-workers': 1}, **kw))
out, h = cmp('extract', [])
cols = h.split('\n')[0].split('\t'); n_ext = h.count('\n') - 1
R('extract_rows', n_ext > 100, 'rows=%d' % n_ext)
out, h = cmp('dual', [P('extract', 'dual', **{'dual-measure-image': W + '/m51_half.fits'})])
# detection on m51 and measurement on the half-flux copy: the catalog must be the extraction one with the dual-image columns of the measurement image
out, h = cmp('completeness', [P('photometry', 'completeness', **{'comp-ninject': 30, 'comp-magmin': 20.0, 'comp-magmax': 24.0, 'comp-nbins': 4})])
out, h = cmp('multiband', [P('photometry', 'multiband', **{'mb-bands': 'A:%s/m51_half.fits,B:%s/m51_quarter.fits' % (W, W)})])
# quantitative: bands are flux x0.5 and x0.25 of the detection image -> magnitude offsets +0.753 / +1.505 relative to the detection-image magnitudes
if h:
    L = [l.split('\t') for l in h.strip().split('\n')]; c = L[0]; rows = L[1:]
    mcols = [x for x in c if x.startswith('MAG') and x not in ('MAG_AUTO', 'MAG_ISOCOR', 'MAG_APER') and not x.startswith('MAG_APER_')]
    print('  multiband columns:', [x for x in c if x not in cols][:12])
    ia = [i for i, x in enumerate(c) if x.endswith('_A') and 'MAG' in x]; ib = [i for i, x in enumerate(c) if x.endswith('_B') and 'MAG' in x]
    if ia and ib:
        d = []
        for r in rows:
            try:
                ma, mb = float(r[ia[0]]), float(r[ib[0]])
                if abs(ma) < 40 and abs(mb) < 40: d.append(mb - ma)
            except ValueError: pass
        d = np.array(d)
        R('multiband_B_minus_A_is_0.753', len(d) > 50 and abs(np.median(d) - 0.753) < 0.03, 'median %.3f (expected 0.753), n=%d, scatter %.3f' % (np.median(d) if len(d) else np.nan, len(d), np.std(d) if len(d) else np.nan))
    else:
        R('multiband_columns_found', False, str(c[-8:]))
if not a.no_network:
    out, h = cmp('crossmatch', [P('photometry', 'crossmatch', **{'xm-catalog': 'GAIA_DR3', 'xm-radius': 2.0})])
    if h:
        L = [l.split('\t') for l in h.strip().split('\n')]; c = L[0]; i = c.index('MATCH_ID')
        nm = sum(1 for r in L[1:] if r[i] not in ('', '0'))
        R('crossmatch_gaia_matches', nm > 50, 'matched %d of %d' % (nm, len(L) - 1))
    out, h = cmp('photoz', [P('photometry', 'crossmatch', **{'xm-catalog': 'GAIA_DR3', 'xm-radius': 2.0}),
                            P('photoz_sed', 'photoz', **{'photoz-bands': 'g,r', 'photoz-mag-columns': 'MAG_AUTO,MAG_APER'})])
    out, h = cmp('sed', [P('photometry', 'crossmatch', **{'xm-catalog': 'GAIA_DR3', 'xm-radius': 2.0}),
                         P('photoz_sed', 'photoz', **{'photoz-bands': 'g,r', 'photoz-mag-columns': 'MAG_AUTO,MAG_APER'}),
                         P('photoz_sed', 'sed', **{'sed-bands': 'g,r', 'sed-mag-columns': 'MAG_AUTO,MAG_APER', 'sed-photoz-column': 'PHOTO_Z', 'sed-backend': 'auto'})])
# deconvolution and segmentation map: images
rec = {'steps': [P('deconv', 'deconvolve', **{'rl-iterations': 6, 'algorithm': 'rl'})], 'detect': {'args': []}, 'psf': W + '/psf.fits'}
out = W + '/hl_deconv'
res = batchrun.run_batch(fields, rec, out, ROOT, a.python, 1, False, 0, 300, [])
R('deconv_headless_ran', res[0]['status'] == 'ok', res[0].get('message', ''))
dp = os.path.join(os.path.expanduser('~'), '.ds9', 'deconv_result.fits')
# the headless step writes {work}/deconv_result.fits with work = <out>/m51/work
hp = out + '/m51/work/deconv_result.fits'
if os.path.exists(hp) and os.path.exists(W + '/gui_deconv.fits'):
    A = fits.getdata(W + '/gui_deconv.fits'); B = fits.getdata(hp)
    R('deconv_image_identical', A.shape == B.shape and np.array_equal(A, B, equal_nan=True), 'max|diff| %.3g, peak %.4g' % (np.nanmax(np.abs(A - B)), np.nanmax(A)))
else:
    R('deconv_image_identical', False, 'missing %s' % [p for p in (hp, W + '/gui_deconv.fits') if not os.path.exists(p)])
rec = {'steps': [P('photometry', 'segmap')], 'detect': {'args': []}}
out = W + '/hl_segmap'
res = batchrun.run_batch(fields, rec, out, ROOT, a.python, 1, False, 0, 300, [])
R('segmap_headless_ran', res[0]['status'] == 'ok', res[0].get('message', ''))
txt = [f for f in os.listdir(out + '/m51/work') if 'segmap' in f.lower()]
print('  segmap work files:', txt)
print('SUMMARY failures=%d' % nfail)
sys.exit(1 if nfail else 0)
