"""Public data used by the regression set.  Everything is cached under $OGF_DATA_CACHE (default ~/.cache/ogfinder_regression) and fetched on demand;
a dataset that cannot be obtained (offline, archive down) makes its cases SKIP, never FAIL."""
import os
import sys

CACHE = os.path.expanduser(os.environ.get('OGF_DATA_CACHE', '~/.cache/ogfinder_regression'))
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HUDF_URL = 'https://archive.stsci.edu/hlsps/hudf12/hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits'    # HST WFC3/IR F160W, 60 mas, HUDF12 HLSP (MAST)
HUDF_BOX = (1000, 1000, 2400, 2400)                                                                     # x0, y0, x1, y1 of the crop (0-based, numpy slicing)
LOCAL_HUDF = [os.environ.get('OGF_FITS_DIR', '/workspace/fits') + '/hudf_f160w.fits']


class Unavailable(Exception):
    pass


def _download(url, dest, timeout=900):
    import requests
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            with open(dest + '.part', 'wb') as f:
                for ch in r.iter_content(1 << 20):
                    f.write(ch)
        os.replace(dest + '.part', dest)
    except Exception as e:
        raise Unavailable('download of %s failed: %s' % (url, str(e)[:100]))


def hudf_crop():
    """1400 x 1400 pixel crop of the HUDF12 WFC3/IR F160W mosaic (HUDF12 HLSP) -> path of a small FITS file"""
    dest = os.path.join(CACHE, 'hudf12_f160w_crop.fits')
    if os.path.isfile(dest):
        return dest
    from astropy.io import fits
    src = next((p for p in LOCAL_HUDF if os.path.isfile(p)), None)
    if src is None:
        src = os.path.join(CACHE, 'hudf12_f160w_full.fits')
        if not os.path.isfile(src):
            _download(HUDF_URL, src)
    with fits.open(src) as h:
        hd = h[0].header.copy()
        x0, y0, x1, y1 = HUDF_BOX
        d = h[0].data[y0:y1, x0:x1].astype('float32')
    for k in ('CRPIX1', 'CRPIX2'):
        hd[k] = hd[k] - (x0 if k == 'CRPIX1' else y0)
    hd['HISTORY'] = 'crop [%d:%d, %d:%d] of %s' % (y0, y1, x0, x1, os.path.basename(HUDF_URL))
    os.makedirs(CACHE, exist_ok=True)
    fits.writeto(dest, d, hd, overwrite=True)
    return dest


def j0946():
    """SDSS J0946+1006 (SLACS), HST ACS/WFC F814W, program 10886, 500 x 500 pixel cutout (MAST)"""
    sys.path.insert(0, os.path.join(ROOT, 'plugins', 'lensmodel', 'validation'))
    import lens_real_j0946 as LR
    try:
        return LR.get_cutout()
    except Exception as e:
        raise Unavailable('J0946 cutout: %s' % str(e)[:100])


def sdss_run94():
    d = os.path.join(CACHE, 'sdss_run94')
    if not os.path.isdir(d):
        raise Unavailable('SDSS run 94 frames are not cached and are large: run moving/validation/sdss_known_asteroids.py once (fetches them)')
    return d
