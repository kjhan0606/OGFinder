#!/usr/bin/env python3
"""Real barred galaxy test: NGC 3351 (M95, SB(r)b) in the SDSS r-band frame 3836-4-84 (public, downloaded from data.sdss.org and cached).

    python barred_real.py OUT.json [--cache DIR] [--half 120] [--galfit BIN --ld DIR]
Bulge + disc base fit, residual-based bar / ring / spiral detection and refit (ogfkit.structfit), stars masked with sep, PSF = stacked frame stars (azimuthal wings).
Literature (for comparison only): bar position angle ~ 112 deg E of N, semi-major axis ~ 50 arcsec (Erwin 2005; Martinet & Friedli 1997).
The bar PA is converted from the multifit convention (counter-clockwise from +x) to degrees east of north with the frame WCS.  Optionally the final model is passed to GALFIT as start
(Ferrers bar + exp disc + Sersic bulge) and refit there; the parameters of both programs are listed."""
def _import_ogfmeas():
    """In-tree measurements. Finds ogfmeas from this file so a script does not need PYTHONPATH."""
    import pathlib
    import sys
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas
    return ogfmeas.measurement_library()

import argparse
import bz2
import json
import math
import os
import shutil
import subprocess
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path[:0] = [ROOT, HERE]
from astropy.io import fits  # noqa: E402
from ogfkit import multifit as MF, structfit as SF  # noqa: E402

URL = 'https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/3836/4/frame-r-003836-4-0084.fits.bz2'
RA, DEC = 160.9906, 11.7037
PIX = 0.396          # arcsec per unbinned pixel


def fetch(cache):
    os.makedirs(cache, exist_ok=True)
    f = os.path.join(cache, 'frame-r-003836-4-0084.fits')
    if not os.path.exists(f):
        with urllib.request.urlopen(URL, timeout=120) as r:
            raw = r.read()
        open(f, 'wb').write(bz2.decompress(raw))
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--cache', default=os.environ.get('OGF_DATA_CACHE', os.path.expanduser('~/.cache/ogfinder_regression')))
    ap.add_argument('--half', type=int, default=100, help='half-size of the cutout in BINNED pixels')
    ap.add_argument('--bin', type=int, default=2)
    a = ap.parse_args()
    sep = _import_ogfmeas()
    from astropy.wcs import WCS
    import warnings
    warnings.simplefilter('ignore')
    f = fetch(a.cache)
    with fits.open(f) as h:
        img = np.array(h[0].data, float)
        hdr = h[0].header
        nmgy = np.array(h[1].data, float)                     # calibration vector (nanomaggies per count) along x
    wcs = WCS(hdr)
    x0, y0 = [float(v) for v in wcs.all_world2pix([[RA, DEC]], 0)[0]]
    img = img * np.mean(nmgy)                                 # nanomaggies -> zero point 22.5
    B = a.bin
    ny_, nx_ = (img.shape[0] // B) * B, (img.shape[1] // B) * B
    img = img[:ny_, :nx_].reshape(ny_ // B, B, nx_ // B, B).sum(axis=(1, 3))     # B x B binning (sum): zero point unchanged, PSF narrower in pixels
    x0, y0 = (x0 + 0.5) / B - 0.5, (y0 + 0.5) / B - 0.5
    bk = sep.Background(np.ascontiguousarray(img), bw=64, bh=64)
    ds = img - bk.back()
    rms = bk.globalrms
    o, seg = sep.extract(ds, 3.0, err=bk.rms(), minarea=6, deblend_nthresh=32, deblend_cont=0.005, segmentation_map=True)
    from galfit_real_compare import stack_psf
    fr, _ = sep.flux_radius(ds, o['x'], o['y'], 6 * np.maximum(o['a'], 1.0), 0.5)
    psf, nst, frs = stack_psf(ds, bk.rms(), o, fr, size=19, nmax=60, snr_min=60)
    h_ = a.half
    ix, iy = int(round(x0)), int(round(y0))
    sl = (slice(iy - h_, iy + h_ + 1), slice(ix - h_, ix + h_ + 1))
    cut = ds[sl].copy()
    sg = seg[sl]
    from scipy.ndimage import binary_dilation
    obj_area = np.bincount(seg.ravel())
    lab_c = sg[h_, h_]
    mask = np.zeros(cut.shape, bool)
    for lab in np.unique(sg[sg > 0]):
        if lab == lab_c:
            continue
        m = sg == lab
        yy, xx = np.nonzero(m)
        if np.hypot(xx.mean() - h_, yy.mean() - h_) < 15 and obj_area[lab] > 300:      # bright knots inside the nucleus region are galaxy, not stars
            continue
        mask |= binary_dilation(m, iterations=3)
    mask[np.hypot(*np.mgrid[:2 * h_ + 1, :2 * h_ + 1] - h_) < 4] = False
    flux0 = float(cut[~mask & (np.hypot(*np.mgrid[:2 * h_ + 1, :2 * h_ + 1] - h_) < 0.8 * h_)].sum())
    q0, pa0 = 0.75, 0.0
    t_flux = max(flux0, 1.0)
    r = SF.fit_structure(cut, h_ + (x0 - ix), h_ + (y0 - iy), t_flux, 30.0, q0, pa0, psf=psf, rms=rms, mask=mask, zp=22.5, max_nfev=200, features=('bar', 'ring', 'spiral'), geometry='outer')
    out = dict(source=URL, galaxy='NGC 3351 (M95)', cutout=2 * h_ + 1, pixel_scale=PIX * B, psf_stars=nst, selected=r['name'], notes=r['notes'], bic={k: float(v) for k, v in r['bic'].items()}, scale=float(r['scale']),
               chi2red={k: float(v['chi2_red']) for k, v in r['all'].items()}, masked_fraction=float(mask.mean()))
    # sky-angle conversion: multifit PA is ccw from +x (pixel frame); degrees east of north via the WCS
    cd = wcs.pixel_scale_matrix if hasattr(wcs, 'pixel_scale_matrix') else None

    def pa_en(pa_xccw):
        t = math.radians(pa_xccw)
        dxy = np.array([math.cos(t), math.sin(t)])
        p0 = np.array([[x0, y0]]); p1 = p0 + 10 * dxy[None, :]
        s0 = wcs.all_pix2world(p0, 0)[0]; s1 = wcs.all_pix2world(p1, 0)[0]
        de = (s1[0] - s0[0]) * math.cos(math.radians(s0[1])); dn = s1[1] - s0[1]
        return float(math.degrees(math.atan2(de, dn)) % 180.0)
    for name, res in r['all'].items():
        cs = res['components']
        if name == 'bar+ring':
            out['ring_in_bar_ring'] = dict(r_arcsec=float(cs[3]['rring'] * PIX * B), w_arcsec=float(cs[3]['sring'] * PIX * B))
        tot = sum(10 ** (-0.4 * c['mag']) for c in cs)
        out['model_' + name] = [dict(kind=c['kind'], frac=10 ** (-0.4 * c['mag']) / tot, **{k: float(c[k]) for k in ('re', 'n', 'q', 'pa', 'rout', 'c0', 'rring', 'sring') if k in c}) for c in cs]
        if name in ('bar', 'bar+ring') and 'bar' not in out or name == r['name'] and name in ('bar', 'bar+ring'):
            b = cs[2]
            out['bar'] = dict(model=name, length_pix=float(b['rout']), length_arcsec=float(b['rout'] * PIX * B), pa_xccw=float(b['pa']), pa_east_of_north=pa_en(b['pa']), q=float(b['q']), c0=float(b['c0']), frac=out['model_bar'][2]['frac'])
            d = cs[1]
            out['disc_pa_east_of_north'] = pa_en(d['pa']); out['disc_q'] = float(d['q']); out['disc_re_arcsec'] = float(d['re'] * PIX * B)
    np.save(os.path.join(a.cache, 'm95_cut.npy'), cut)
    json.dump(out, open(a.out, 'w'), indent=1)
    print({k: out[k] for k in ('selected', 'bic', 'chi2red', 'notes')})
    if 'bar' in out:
        print('bar', out['bar'], 'disc q %.2f PA %.0f re %.0f"' % (out['disc_q'], out['disc_pa_east_of_north'], out['disc_re_arcsec']))


if __name__ == '__main__':
    main()
