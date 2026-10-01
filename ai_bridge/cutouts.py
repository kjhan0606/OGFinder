"""Cutout generation helper: per object, per band; pixel or arcsec size; asinh / zscale / linear / none;
FITS / PNG / npy output.  Needs numpy + astropy (imported lazily; the rest of the bridge is stdlib only).

Convention: the first image is the reference grid; catalog x,y (1-based, SExtractor) refer to it.  For other
images with a different pixel grid the object's RA/Dec is converted through that image's WCS.
"""
import base64
import hashlib
import io
import os
import struct
import zlib

NORMALIZATIONS = ('none', 'linear', 'zscale', 'asinh')
FORMATS = ('png', 'fits', 'npy')


def _np():
    import numpy as np
    return np


def write_png_gray(path_or_buf, arr8):
    """Minimal 8-bit grayscale PNG writer (stdlib zlib; no Pillow needed). arr8: 2-D uint8 array, row 0 = top."""
    h, w = arr8.shape
    raw = b''.join(b'\x00' + arr8[i].tobytes() for i in range(h))

    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    png = (b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 0, 0, 0, 0))
           + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
    if hasattr(path_or_buf, 'write'):
        path_or_buf.write(png)
    else:
        with open(path_or_buf, 'wb') as f:
            f.write(png)
    return png


def normalize(data, how='asinh'):
    """-> float32 array in [0,1] (NaN -> 0).  how: none (returns data as float32 unchanged), linear (1..99.5 %),
    zscale (astropy ZScaleInterval, linear), asinh (zscale interval, asinh stretch a=0.1)."""
    np = _np()
    d = np.asarray(data, dtype=np.float64)
    if how == 'none':
        return d.astype(np.float32)
    fin = d[np.isfinite(d)]
    if fin.size == 0:
        return np.zeros(d.shape, np.float32)
    if how == 'linear':
        lo, hi = np.percentile(fin, [1.0, 99.5])
    elif how in ('zscale', 'asinh'):
        from astropy.visualization import ZScaleInterval
        lo, hi = ZScaleInterval().get_limits(fin)
    else:
        raise ValueError('unknown normalization %r (use %s)' % (how, '|'.join(NORMALIZATIONS)))
    if not hi > lo:
        hi = lo + 1.0
    x = np.clip((np.nan_to_num(d, nan=lo) - lo) / (hi - lo), 0.0, 1.0)
    if how == 'asinh':
        a = 0.1
        x = np.arcsinh(x / a) / np.arcsinh(1.0 / a)
    return x.astype(np.float32)


def header_pixel_scale(path):
    """arcsec/pixel from the first image HDU header only (cheap), or None."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=True) as hl:
        for h in hl:
            if h.header.get('NAXIS', 0) >= 2:
                hd = h.header
                try:
                    w = WCS(hd)
                    if w.has_celestial:
                        from astropy.wcs.utils import proj_plane_pixel_scales
                        return float(proj_plane_pixel_scales(w)[0]) * 3600.0
                except Exception:
                    pass
                for k in ('PIXSCALE', 'PIXSCAL1'):
                    if k in hd:
                        return float(hd[k])
                return None
    return None


class ImageSet:
    """Lazily opened FITS images {band: path}; reads the first image HDU with 2-D data."""

    def __init__(self, images):
        from astropy.io import fits
        from astropy.wcs import WCS
        self.bands = list(images)
        self.paths = dict(images)
        self.data, self.hdr, self.wcs = {}, {}, {}
        for b, p in images.items():
            with fits.open(p, memmap=False) as hl:
                for h in hl:
                    if h.data is not None and h.data.ndim >= 2:
                        d = h.data
                        while d.ndim > 2:
                            d = d[0]
                        self.data[b] = d
                        self.hdr[b] = h.header.copy()
                        try:
                            w = WCS(h.header)
                            self.wcs[b] = w if w.has_celestial else None
                        except Exception:
                            self.wcs[b] = None
                        break
                else:
                    raise ValueError('no 2-D image in %s' % p)
        self.ref = self.bands[0]

    def pixel_scale_arcsec(self, band=None):
        band = band or self.ref
        h = self.hdr[band]
        w = self.wcs.get(band)
        if w is not None:
            try:
                from astropy.wcs.utils import proj_plane_pixel_scales
                return float(proj_plane_pixel_scales(w)[0]) * 3600.0
            except Exception:
                pass
        for k in ('PIXSCALE', 'PIXSCAL1'):
            if k in h:
                return float(h[k])
        return None

    def center_pixel(self, band, rec):
        """0-based (cx, cy) for this object in this band, or None."""
        if band == self.ref or self.data[band].shape == self.data[self.ref].shape and self._same_wcs(band):
            if rec.get('x') is None or rec.get('y') is None:
                return None
            return rec['x'] - 1.0, rec['y'] - 1.0
        w = self.wcs.get(band)
        if w is not None and rec.get('ra') is not None and rec.get('dec') is not None:
            px, py = w.all_world2pix([[rec['ra'], rec['dec']]], 0)[0]
            return float(px), float(py)
        if rec.get('x') is not None and rec.get('y') is not None and self.data[band].shape == self.data[self.ref].shape:
            return rec['x'] - 1.0, rec['y'] - 1.0
        return None

    def _same_wcs(self, band):
        a, b = self.hdr[self.ref], self.hdr[band]
        for k in ('CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD2_2', 'CDELT1', 'CDELT2'):
            if a.get(k) != b.get(k):
                return False
        return True


def cut(data, cx, cy, size):
    """size x size cutout centred on 0-based (cx, cy) (rounded to a pixel); NaN padded outside the image.
    Returns (array, x0, y0) with x0,y0 the 0-based lower-left corner in the parent image."""
    np = _np()
    size = int(size)
    ix, iy = int(round(cx)), int(round(cy))
    x0, y0 = ix - size // 2, iy - size // 2
    out = np.full((size, size), np.nan, dtype=np.float32)
    ny, nx = data.shape
    sx0, sy0 = max(x0, 0), max(y0, 0)
    sx1, sy1 = min(x0 + size, nx), min(y0 + size, ny)
    if sx1 > sx0 and sy1 > sy0:
        out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = data[sy0:sy1, sx0:sx1]
    return out, x0, y0


class CutoutMaker:
    """cfg keys: size_pix | size_arcsec, normalize, format (png|fits|npy), flip (png: row 0 = top, default True)."""

    def __init__(self, imageset, outdir, size_pix=None, size_arcsec=None, norm='asinh', fmt='png', pixel_scale=None):
        if fmt not in FORMATS:
            raise ValueError('cutout format %r not in %s' % (fmt, FORMATS))
        if norm not in NORMALIZATIONS:
            raise ValueError('normalize %r not in %s' % (norm, NORMALIZATIONS))
        self.im, self.outdir, self.norm, self.fmt = imageset, outdir, norm, fmt
        self.pixel_scale = pixel_scale
        if size_pix is None and size_arcsec is None:
            size_pix = 64
        self.size_pix, self.size_arcsec = size_pix, size_arcsec
        os.makedirs(outdir, exist_ok=True)

    def size_for(self, band):
        if self.size_pix is not None:
            return int(self.size_pix)
        ps = self.pixel_scale or self.im.pixel_scale_arcsec(band)
        if not ps:
            raise ValueError('size_arcsec needs a pixel scale (WCS or --pixel-scale)')
        return max(8, int(round(self.size_arcsec / ps)))

    def make(self, rec, band, fmt=None):
        """-> dict {path, format, size, x0, y0, sha256} or None when the object is not on this image."""
        np = _np()
        fmt = fmt or self.fmt
        if fmt not in FORMATS:
            raise ValueError('cutout format %r not in %s' % (fmt, FORMATS))
        c = self.im.center_pixel(band, rec)
        if c is None:
            return None
        size = self.size_for(band)
        arr, x0, y0 = cut(self.im.data[band], c[0], c[1], size)
        if not np.isfinite(arr).any():
            return None
        name = 'obj%s_%s.%s' % (rec['id'], band, fmt)
        path = os.path.join(self.outdir, name)
        if fmt == 'png':
            n = normalize(arr, self.norm if self.norm != 'none' else 'linear')
            a8 = np.flipud((n * 255.0 + 0.5).astype(np.uint8))
            write_png_gray(path, a8)
        elif fmt == 'npy':
            np.save(path, normalize(arr, self.norm) if self.norm != 'none' else arr.astype(np.float32))
        else:
            from astropy.io import fits
            h = fits.Header()
            src = self.im.hdr[band]
            for k in src:
                if k.startswith(('CTYPE', 'CRVAL', 'CD', 'PC', 'CDELT', 'CUNIT', 'RADESYS', 'EQUINOX', 'LONPOLE', 'LATPOLE',
                                 'PIXSCALE', 'FILTER', 'BUNIT')) and k not in h:
                    h[k] = src[k]
            h['CRPIX1'] = float(src.get('CRPIX1', 1.0)) - x0
            h['CRPIX2'] = float(src.get('CRPIX2', 1.0)) - y0
            h['OBJID'] = str(rec['id'])
            h['CUTX0'] = x0 + 1
            h['CUTY0'] = y0 + 1
            h['COMMENT'] = 'OGFinder ai_bridge cutout (raw pixel values)'
            fits.PrimaryHDU(data=arr, header=h).writeto(path, overwrite=True)
        with open(path, 'rb') as f:
            blob = f.read()
        return {'path': path, 'format': fmt, 'size': size, 'x0': int(x0), 'y0': int(y0),
                'sha256': hashlib.sha256(blob).hexdigest(), 'bytes': len(blob),
                'wcs': _wcs_keywords(self.im.hdr[band], x0, y0, size)}


def _wcs_keywords(hdr, x0, y0, size):
    out = {}
    for k in ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2', 'RADESYS', 'EQUINOX'):
        if k in hdr:
            v = hdr[k]
            out[k] = v if isinstance(v, (int, float, str)) else str(v)
    out['CRPIX1'] = float(hdr.get('CRPIX1', 1.0)) - x0
    out['CRPIX2'] = float(hdr.get('CRPIX2', 1.0)) - y0
    out['NAXIS1'] = out['NAXIS2'] = int(size)
    return out


def b64_file(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('ascii')
