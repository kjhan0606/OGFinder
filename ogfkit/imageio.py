"""FITS image / mask helpers.  Masks follow the shared mask manager (plugins/mask, ds9/library/ds9_mask.py)."""
import numpy as np

# bits of the flag image written by ds9_mask.py
SRC, STAR, ADD, ERASE, IMPORT = 1, 2, 4, 8, 16
MASKBITS = SRC | STAR | ADD | IMPORT


def effective_mask(flags):
    """Bit-flag mask -> boolean 'masked' array: (flags & 23) != 0 and (flags & 8) == 0 (same rule as ds9_mask.py)."""
    f = np.asarray(flags).astype(np.int64)
    return ((f & MASKBITS) != 0) & ((f & ERASE) == 0)


def load_image(path, hdu=None):
    """-> (float32 2-D array, astropy header).  First HDU with 2-D data unless hdu is given."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as hl:
        if hdu is None:
            for i, h in enumerate(hl):
                if h.data is not None and h.data.ndim == 2:
                    hdu = i
                    break
            else:
                raise ValueError('no 2-D image in %s' % path)
        data = np.array(hl[hdu].data, dtype=np.float32)
        hdr = hl[hdu].header.copy()
    return data, hdr


def load_mask(path, shape=None):
    """Mask file -> bool array (True = masked).  A file with values > 1 is read as the bit-flag image of the mask
    manager (see effective_mask); a 0/1 file (the manager's *_bool.fits) is used as is.  Returns None for a falsy path."""
    if not path:
        return None
    from astropy.io import fits
    with fits.open(path, memmap=False) as hl:
        a = np.asarray(hl[0].data if hl[0].data is not None else hl[1].data)
    m = effective_mask(a) if a.max() > 1 else (a != 0)
    if shape is not None and m.shape != tuple(shape):
        raise ValueError('mask shape %s differs from image shape %s' % (m.shape, tuple(shape)))
    return m


def save_fits(path, data, header=None, extra=None):
    """Write a float32 FITS image; copy the WCS-related cards of header; extra = dict of cards to add."""
    from astropy.io import fits
    h = fits.Header()
    if header is not None:
        for k in header.keys():
            if k in ('SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND', 'BSCALE', 'BZERO', 'COMMENT', 'HISTORY', ''):
                continue
            try:
                h[k] = header[k]
            except Exception:
                pass
    for k, v in (extra or {}).items():
        h[k] = v
    fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header=h).writeto(path, overwrite=True)
    return path


def robust_sigma(a, mask=None):
    """1.4826 * MAD of the finite, unmasked pixels (iteratively clipped at 4 sigma)."""
    v = np.asarray(a, dtype=float)
    ok = np.isfinite(v)
    if mask is not None:
        ok &= ~mask
    v = v[ok]
    if v.size == 0:
        return float('nan')
    for _ in range(3):
        med = np.median(v)
        s = 1.4826 * np.median(np.abs(v - med))
        if s <= 0:
            break
        keep = np.abs(v - med) < 4 * s
        if keep.all():
            break
        v = v[keep]
    return float(1.4826 * np.median(np.abs(v - np.median(v))))
