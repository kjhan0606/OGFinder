#!/usr/bin/env python3
"""
OGFinder shared mask manager (one mask per image, bit-flag FITS).

Bits of the flag image (uint8):
     1  detected source (segmentation + dilation)
     2  bright star / diffraction region
     4  manual add (regions, grow, invert)
     8  manual erase / protect  (wins over every masking bit)
    16  imported / reprojected from another mask or band
    32  satellite / aircraft trail (plugins/trails, ogfkit/trails.py); survives Auto Mask, removable with --mode trails-clear
A pixel is *masked* when (flags & 55) != 0 and (flags & 8) == 0.

Files (MASK = ~/.ds9/mask_<base>.fits):
    MASK                       flag image (this program's master copy)
    MASK with _bool.fits       0/1 image of the effective mask; this is what
                               ds9_icl.py / ds9_lsbg.py read via --mask
    MASK stem + .undo/ .redo/  compressed snapshots (undo/redo stacks)

Usage:
    ds9_mask.py IMAGE --mode auto --mask MASK [--detect-thresh 5 --minarea 5 --expand-factor 1.5
                                  --max-dilate-radius 20 --bright-star-mag-limit 18
                                  --bright-star-radius-scale 10 --mag-threshold 99
                                  --mag-zeropoint 25 --catalog cat.tsv --fresh]
    ds9_mask.py IMAGE --mode add|erase --mask MASK --regions file.reg
    ds9_mask.py IMAGE --mode grow|shrink --mask MASK --pixels N
    ds9_mask.py IMAGE --mode invert|clear|undo|redo|stats --mask MASK
    ds9_mask.py IMAGE --mode trails --mask MASK --file trailmask.fits [--fresh]   (set bit 32 where the file is non-zero; --fresh first clears old trail bits)
    ds9_mask.py IMAGE --mode trails-clear --mask MASK                            (clear bit 32)
    ds9_mask.py IMAGE --mode import --mask MASK --file other.fits
    ds9_mask.py IMAGE --mode export --mask MASK --file out.fits [--boolean]
    ds9_mask.py IMAGE --mode reproject --mask MASK --target-image band.fits --output band_mask.fits
    ds9_mask.py IMAGE --mode masked --mask MASK --masked-output out.fits [--interp-method linear]
Every mode ends with one line:  #MASK_STATS key=value ...
"""
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

import sys
import os
import re
import glob
import shutil
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

SRC, STAR, ADD, ERASE, IMPORT, TRAIL = 1, 2, 4, 8, 16, 32
MASKBITS = SRC | STAR | ADD | IMPORT | TRAIL
KEEP_ON_AUTO = ADD | ERASE | IMPORT | TRAIL
MAX_UNDO = 10

WCS_KEYS = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2',
            'CD2_1', 'CD2_2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
            'CDELT1', 'CDELT2', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2',
            'CRVAL3', 'EQUINOX', 'RADESYS', 'LONPOLE', 'LATPOLE']


def log(*a):
    print(*a, file=sys.stderr)


def effective(flags):
    return ((flags & MASKBITS) != 0) & ((flags & ERASE) == 0)


# ------------------------------------------------------------------ I/O
def load_image(path):
    from astropy.io import fits
    with fits.open(path) as hdul:
        for h in hdul:
            if h.data is not None and h.data.ndim >= 2:
                data = h.data
                header = h.header.copy()
                break
        else:
            sys.exit("ERROR: no image data in %s" % path)
    while data.ndim > 2:
        data = data[0]
    return np.asarray(data, dtype=np.float64), header


def image_header(path):
    from astropy.io import fits
    with fits.open(path) as hdul:
        for h in hdul:
            if h.data is not None and h.data.ndim >= 2:
                return h.header.copy()
    sys.exit("ERROR: no image in %s" % path)


def bool_path(mask_path):
    return re.sub(r'\.fits?$', '', mask_path) + '_bool.fits'


def stack_dir(mask_path, which):
    return re.sub(r'\.fits?$', '', mask_path) + '.' + which


def write_fits(arr, header, path, extra=None):
    from astropy.io import fits
    hdu = fits.PrimaryHDU(arr)
    if header is not None:
        for k in WCS_KEYS:
            if k in header:
                hdu.header[k] = header[k]
    for k, (v, c) in (extra or {}).items():
        hdu.header[k] = (v, c)
    tmp = path + '.tmp'
    hdu.writeto(tmp, overwrite=True)
    os.replace(tmp, path)


def save_mask(flags, header, path):
    flags = flags.astype(np.uint8)
    write_fits(flags, header, path,
               {'OGFMASK': (1, 'OGFinder bit-flag mask'),
                'MSKBITS': ('1src 2star 4add 8erase 16imp 32trail', 'flag meaning')})
    write_fits(effective(flags).astype(np.uint8), header, bool_path(path),
               {'OGFBOOL': (1, 'effective boolean mask')})


def load_mask(path, shape=None):
    from astropy.io import fits
    if not os.path.exists(path):
        return None, None
    with fits.open(path) as hdul:
        d = hdul[0].data
        h = hdul[0].header.copy()
    if d is None:
        return None, None
    if 'OGFMASK' not in h:      # legacy boolean-like mask: treat nonzero as source
        d = np.where(d != 0, SRC, 0)
    d = np.asarray(d).astype(np.uint8)
    if shape is not None and d.shape != shape:
        log("WARNING: mask shape %s != image shape %s; ignoring mask" % (d.shape, shape))
        return None, None
    return d, h


# ------------------------------------------------------------ undo/redo
def _snap_files(d):
    return sorted(glob.glob(os.path.join(d, '*.npz')))


def push_snapshot(mask_path, flags, which='undo', clear_redo=True):
    d = stack_dir(mask_path, which)
    os.makedirs(d, exist_ok=True)
    files = _snap_files(d)
    n = int(os.path.basename(files[-1])[:-4]) + 1 if files else 1
    np.savez_compressed(os.path.join(d, '%06d.npz' % n),
                        flags=flags if flags is not None else np.zeros(1, np.uint8),
                        none=np.array(flags is None))
    files = _snap_files(d)
    for f in files[:-MAX_UNDO]:
        os.remove(f)
    if which == 'undo' and clear_redo:
        for f in _snap_files(stack_dir(mask_path, 'redo')):
            os.remove(f)


def pop_snapshot(mask_path, which):
    d = stack_dir(mask_path, which)
    files = _snap_files(d)
    if not files:
        return False, None
    with np.load(files[-1]) as z:
        flags = None if bool(z['none']) else z['flags']
    os.remove(files[-1])
    return True, flags


def depth(mask_path):
    return (len(_snap_files(stack_dir(mask_path, 'undo'))),
            len(_snap_files(stack_dir(mask_path, 'redo'))))


# ---------------------------------------------------------------- stats
def stats_line(flags, data=None, mask_path=None, note=''):
    if flags is None:
        return "#MASK_STATS N_MASKED=0 FRACTION=0 N_PIXELS=0 EXISTS=0 UNDO=0 REDO=0"
    eff = effective(flags)
    n = int(eff.sum())
    tot = flags.size
    s = ("#MASK_STATS N_MASKED=%d FRACTION=%.6f N_PIXELS=%d EXISTS=1 "
         "N_SRC=%d N_STAR=%d N_ADD=%d N_ERASE=%d N_IMPORT=%d N_TRAIL=%d" %
         (n, n / tot, tot, int(((flags & SRC) != 0).sum()),
          int(((flags & STAR) != 0).sum()), int(((flags & ADD) != 0).sum()),
          int(((flags & ERASE) != 0).sum()), int(((flags & IMPORT) != 0).sum()),
          int(((flags & TRAIL) != 0).sum())))
    if data is not None:
        nodata = (data == 0) | ~np.isfinite(data)
        s += " N_NODATA=%d FRAC_OF_VALID=%.6f" % (
            int(nodata.sum()),
            int((eff & ~nodata).sum()) / max(1, int((~nodata).sum())))
    if mask_path:
        u, r = depth(mask_path)
        s += " UNDO=%d REDO=%d" % (u, r)
    if note:
        s += " NOTE=%s" % note.replace(' ', '_')
    return s


# ------------------------------------------------------------- regions
_num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def parse_regions(path):
    """Return (include_shapes, exclude_shapes); each shape = (name, [floats])."""
    inc, exc = [], []
    system = 'image'
    with open(path) as f:
        for raw in f:
            line = raw.split('#')[0].strip()
            if not line or line.startswith('global'):
                continue
            low = line.lower()
            if low in ('image', 'physical', 'linear'):
                system = 'image'
                continue
            if low in ('fk5', 'fk4', 'icrs', 'galactic', 'ecliptic', 'wcs'):
                system = 'sky'
                continue
            for part in line.split(';'):
                part = part.strip()
                m = re.match(r'^([-+]?)\s*(\w+)\s*\(([^)]*)\)', part)
                if not m:
                    continue
                sign, name, args = m.groups()
                name = name.lower()
                if name not in ('circle', 'ellipse', 'box', 'polygon', 'annulus'):
                    continue
                if system == 'sky':
                    sys.exit("ERROR: regions must be in image coordinates")
                vals = [float(x) for x in re.findall(_num, args)]
                (exc if sign == '-' else inc).append((name, vals))
    return inc, exc


def rasterize(shapes, shape):
    """Boolean image of union of shapes; DS9 image coords (1-based centres)."""
    ny, nx = shape
    out = np.zeros(shape, bool)
    for name, v in shapes:
        if name in ('circle', 'annulus', 'ellipse', 'box') and len(v) < 3:
            continue
        if name == 'polygon':
            if len(v) < 6:
                continue
            px, py = np.array(v[0::2]), np.array(v[1::2])
            x0 = max(0, int(np.floor(px.min())) - 2)
            x1 = min(nx, int(np.ceil(px.max())) + 1)
            y0 = max(0, int(np.floor(py.min())) - 2)
            y1 = min(ny, int(np.ceil(py.max())) + 1)
            if x1 <= x0 or y1 <= y0:
                continue
            yy, xx = np.mgrid[y0:y1, x0:x1]
            X, Y = xx + 1.0, yy + 1.0
            inside = np.zeros(X.shape, bool)
            n = len(px)
            j = n - 1
            for i in range(n):
                xi, yi, xj, yj = px[i], py[i], px[j], py[j]
                cond = ((yi > Y) != (yj > Y))
                with np.errstate(divide='ignore', invalid='ignore'):
                    xint = (xj - xi) * (Y - yi) / (yj - yi) + xi
                inside ^= cond & (X < xint)
                j = i
            out[y0:y1, x0:x1] |= inside
            continue
        cx, cy = v[0], v[1]
        if name == 'circle':
            a = b = v[2]; ang = 0.0; boxw = None
        elif name == 'annulus':
            a = b = v[-1]; ang = 0.0; boxw = None
        elif name == 'ellipse':
            a, b = v[2], (v[3] if len(v) > 3 else v[2])
            ang = v[4] if len(v) > 4 else 0.0
            boxw = None
        else:  # box
            w, h = v[2], (v[3] if len(v) > 3 else v[2])
            ang = v[4] if len(v) > 4 else 0.0
            a, b = w / 2.0, h / 2.0; boxw = True
        rr = np.hypot(a, b) + 1
        x0 = max(0, int(np.floor(cx - rr)) - 1); x1 = min(nx, int(np.ceil(cx + rr)) + 1)
        y0 = max(0, int(np.floor(cy - rr)) - 1); y1 = min(ny, int(np.ceil(cy + rr)) + 1)
        if x1 <= x0 or y1 <= y0:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        dx, dy = xx + 1.0 - cx, yy + 1.0 - cy
        t = np.deg2rad(ang)
        u = dx * np.cos(t) + dy * np.sin(t)
        w_ = -dx * np.sin(t) + dy * np.cos(t)
        if boxw:
            ins = (np.abs(u) <= a) & (np.abs(w_) <= b)
        else:
            ins = (u / a) ** 2 + (w_ / b) ** 2 <= 1.0
        if name == 'annulus' and len(v) >= 4:
            r_in = v[2]
            ins &= (dx * dx + dy * dy) >= r_in * r_in
        out[y0:y1, x0:x1] |= ins
    return out


# ---------------------------------------------------------------- modes
def parse_tsv_catalog(path):
    with open(path) as f:
        lines = [l for l in f.read().splitlines() if l.strip()]
    hdr = lines[0].split('\t')
    idx = {h.strip(): i for i, h in enumerate(hdr)}
    if 'X_IMAGE' not in idx or 'MAG_AUTO' not in idx:
        return None
    x, y, m = [], [], []
    for l in lines[1:]:
        p = l.split('\t')
        try:
            x.append(float(p[idx['X_IMAGE']])); y.append(float(p[idx['Y_IMAGE']]))
            m.append(float(p[idx['MAG_AUTO']]))
        except (ValueError, IndexError):
            continue
    return {'x_image': np.array(x), 'y_image': np.array(y), 'mag_auto': np.array(m)}


def mode_auto(a):
    sep = _import_ogfmeas()
    from icl.masking import create_source_mask, mask_bright_stars
    data, header = load_image(a.image)
    ny, nx = data.shape
    old, _ = load_mask(a.mask, data.shape)
    log("Image %dx%d" % (nx, ny))
    dc = np.ascontiguousarray(data, dtype=np.float64)
    nodata = (dc == 0) | ~np.isfinite(dc)
    dc = np.where(np.isfinite(dc), dc, 0.0)
    bkg = sep.Background(dc, mask=nodata)
    sub = dc - bkg.back()
    objects, seg = sep.extract(sub, a.detect_thresh, err=bkg.rms(), minarea=a.minarea,
                               segmentation_map=True, mask=nodata)
    log("auto mask: detected %d sources (thresh %.2f sigma, minarea %d)" %
        (len(objects), a.detect_thresh, a.minarea))
    if a.mag_threshold < 90:
        from lsbg.masking import mask_by_magnitude
        src, labels = mask_by_magnitude(
            data, seg, objects, mag_threshold=a.mag_threshold,
            expand_factor=a.expand_factor, max_dilate_radius=a.max_dilate_radius,
            lsb_protect=a.lsb_protect, lsb_mu_threshold=a.lsb_mu_threshold,
            pixel_scale=a.pixel_scale, mag_zeropoint=a.mag_zeropoint,
            n_workers=a.n_workers)
        log("auto mask: %d sources brighter than %.1f masked" % (len(labels), a.mag_threshold))
    else:
        src = create_source_mask(data, seg, a.expand_factor,
                                 max_dilate_radius=a.max_dilate_radius,
                                 n_workers=a.n_workers)
    star = np.zeros(data.shape, bool)
    cat = None
    if a.catalog and os.path.exists(a.catalog):
        cat = parse_tsv_catalog(a.catalog)
    if cat is None and len(objects):
        cat = {'x_image': objects['x'] + 1.0, 'y_image': objects['y'] + 1.0,
               'mag_auto': -2.5 * np.log10(np.maximum(objects['flux'], 1e-30)) + a.mag_zeropoint}
    if cat is not None:
        mask_bright_stars(star, cat, a.bright_star_mag_limit, a.bright_star_radius_scale)
    flags = np.zeros(data.shape, np.uint8)
    flags[src] |= SRC
    flags[star] |= STAR
    if old is not None and not a.fresh:
        flags |= (old & KEEP_ON_AUTO)
        push_snapshot(a.mask, old)
        log("auto mask: kept manual bits (add/erase/import) from existing mask")
    else:
        push_snapshot(a.mask, old)
    save_mask(flags, header, a.mask)
    print(stats_line(flags, data, a.mask))


def _edit(a, fn, note=''):
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is None:
        flags = np.zeros(data.shape, np.uint8)
    push_snapshot(a.mask, flags)
    new = fn(flags.copy(), data)
    save_mask(new, header, a.mask)
    print(stats_line(new, data, a.mask, note))


def mode_region(a, erase):
    def fn(flags, data):
        inc, exc = parse_regions(a.regions)
        r = rasterize(inc, flags.shape)
        if exc:
            r &= ~rasterize(exc, flags.shape)
        log("regions: %d included shapes, %d excluded; %d pixels" %
            (len(inc), len(exc), int(r.sum())))
        if erase:
            flags[r] |= ERASE
            flags[r] &= ~np.uint8(ADD)
        else:
            flags[r] |= ADD
            flags[r] &= ~np.uint8(ERASE)
        return flags
    _edit(a, fn, 'erase' if erase else 'add')


def mode_grow(a, shrink):
    from scipy.ndimage import distance_transform_edt
    n = float(a.pixels)

    def fn(flags, data):
        eff = effective(flags)
        if shrink:
            keep = distance_transform_edt(np.pad(eff, 1, constant_values=True))[1:-1, 1:-1] > n
            removed = eff & ~keep
            flags[removed] |= ERASE
            flags[removed] &= ~np.uint8(ADD)
        else:
            new = distance_transform_edt(~eff) <= n
            added = new & ~eff
            flags[added] |= ADD
            flags[added] &= ~np.uint8(ERASE)
        return flags
    _edit(a, fn)


def mode_invert(a):
    def fn(flags, data):
        eff = effective(flags)
        out = np.where(eff, ERASE, ADD).astype(np.uint8)
        return out
    _edit(a, fn)


def mode_clear(a):
    _edit(a, lambda f, d: np.zeros_like(f))


def mode_undo_redo(a, redo):
    data, header = load_image(a.image)
    cur, _ = load_mask(a.mask, data.shape)
    ok, prev = pop_snapshot(a.mask, 'redo' if redo else 'undo')
    if not ok:
        print(stats_line(cur, data, a.mask, 'nothing_to_' + ('redo' if redo else 'undo')))
        return
    push_snapshot(a.mask, cur, 'undo' if redo else 'redo', clear_redo=False)
    if prev is None:
        prev = np.zeros(data.shape, np.uint8)
    save_mask(prev, header, a.mask)
    print(stats_line(prev, data, a.mask))


def mode_stats(a):
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is not None and not os.path.exists(bool_path(a.mask)):
        save_mask(flags, header, a.mask)
    print(stats_line(flags, data, a.mask))


def mode_import(a):
    from astropy.io import fits
    from icl.masking import reproject_mask
    data, header = load_image(a.image)
    with fits.open(a.file) as hdul:
        md = hdul[0].data
        mh = hdul[0].header
    md = np.asarray(md)
    while md.ndim > 2:
        md = md[0]
    flags, _ = load_mask(a.mask, data.shape)
    push_snapshot(a.mask, flags)
    if 'OGFMASK' in mh and md.shape == data.shape:
        new = md.astype(np.uint8)
        log("import: OGFinder flag mask, same grid -> exact copy")
    else:
        r = reproject_mask(md, mh, header)
        new = np.where(r, IMPORT, 0).astype(np.uint8)
    save_mask(new, header, a.mask)
    print(stats_line(new, data, a.mask))


def mode_trails(a, clear=False):
    from astropy.io import fits
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is None:
        flags = np.zeros(data.shape, np.uint8)
    push_snapshot(a.mask, flags)
    new = flags.copy()
    if clear or a.fresh:
        new &= np.uint8(~TRAIL & 0xFF)
    if not clear:
        with fits.open(a.file) as hdul:
            tm = np.asarray(hdul[0].data)
        while tm.ndim > 2:
            tm = tm[0]
        if tm.shape != data.shape:
            sys.exit("ERROR: trail mask shape %s != image shape %s" % (tm.shape, data.shape))
        new[tm != 0] |= TRAIL
        new[tm != 0] &= ~np.uint8(ERASE)
    save_mask(new, header, a.mask)
    print(stats_line(new, data, a.mask, 'trails_cleared' if clear else 'trails_set'))


def mode_export(a):
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is None:
        sys.exit("ERROR: no mask to export")
    if a.boolean:
        write_fits(effective(flags).astype(np.uint8), header, a.file,
                   {'OGFBOOL': (1, 'effective boolean mask')})
    else:
        write_fits(flags, header, a.file, {'OGFMASK': (1, 'OGFinder bit-flag mask')})
    print(stats_line(flags, data, a.mask, 'exported'))


def mode_reproject(a):
    from icl.masking import reproject_mask
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is None:
        sys.exit("ERROR: no mask to reproject")
    th = image_header(a.target_image)
    r = reproject_mask(effective(flags).astype(np.uint8), header, th)
    out = np.where(r, IMPORT, 0).astype(np.uint8)
    old, _ = load_mask(a.output, out.shape)
    if old is not None and not a.fresh:
        push_snapshot(a.output, old)
        out |= (old & (ADD | ERASE))
    save_mask(out, th, a.output)
    print(stats_line(out, None, a.output, 'reprojected'))


def mode_masked(a):
    from icl.masking import interpolate_masked
    data, header = load_image(a.image)
    flags, _ = load_mask(a.mask, data.shape)
    if flags is None:
        sys.exit("ERROR: no mask")
    m = effective(flags)
    md = interpolate_masked(data, m, method=a.interp_method)
    write_fits(md.astype(np.float32), header, a.masked_output)
    print(stats_line(flags, data, a.mask, 'masked_image_written'))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('image')
    p.add_argument('--mode', required=True, choices=['auto', 'add', 'erase', 'grow', 'shrink', 'invert',
                                    'clear', 'undo', 'redo', 'stats', 'import', 'export',
                                    'reproject', 'masked', 'trails', 'trails-clear'])
    p.add_argument('--mask', required=True)
    p.add_argument('--regions')
    p.add_argument('--pixels', type=float, default=1)
    p.add_argument('--file')
    p.add_argument('--boolean', action='store_true')
    p.add_argument('--target-image')
    p.add_argument('--output')
    p.add_argument('--masked-output')
    p.add_argument('--interp-method', default='linear')
    p.add_argument('--fresh', action='store_true', help='auto: discard manual bits')
    p.add_argument('--detect-thresh', type=float, default=5.0)
    p.add_argument('--minarea', type=int, default=5)
    p.add_argument('--expand-factor', type=float, default=1.5)
    p.add_argument('--max-dilate-radius', type=int, default=20)
    p.add_argument('--bright-star-mag-limit', type=float, default=18.0)
    p.add_argument('--bright-star-radius-scale', type=float, default=10.0)
    p.add_argument('--mag-threshold', type=float, default=99.0,
                   help='mask only sources brighter than this (>=90: all sources)')
    p.add_argument('--lsb-protect', action='store_true')
    p.add_argument('--lsb-mu-threshold', type=float, default=24.0)
    p.add_argument('--pixel-scale', type=float, default=0.06)
    p.add_argument('--mag-zeropoint', type=float, default=25.0)
    p.add_argument('--catalog')
    p.add_argument('--n-workers', type=int, default=0)
    a = p.parse_args()
    m = a.mode
    if m == 'auto': mode_auto(a)
    elif m in ('add', 'erase'):
        if not a.regions:
            sys.exit('--regions required')
        mode_region(a, m == 'erase')
    elif m == 'grow': mode_grow(a, False)
    elif m == 'shrink': mode_grow(a, True)
    elif m == 'invert': mode_invert(a)
    elif m == 'clear': mode_clear(a)
    elif m == 'trails': mode_trails(a)
    elif m == 'trails-clear': mode_trails(a, clear=True)
    elif m == 'undo': mode_undo_redo(a, False)
    elif m == 'redo': mode_undo_redo(a, True)
    elif m == 'stats': mode_stats(a)
    elif m == 'import': mode_import(a)
    elif m == 'export': mode_export(a)
    elif m == 'reproject': mode_reproject(a)
    elif m == 'masked': mode_masked(a)


if __name__ == '__main__':
    main()
