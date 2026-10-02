"""Multiwavelength cross-match: table readers (FITS / CSV / TSV / VOTable), sky or pixel matching with separation and match flags, systematic offsets,
chance-match estimate.  Pure functions on numpy arrays; no GUI, no network (the TAP retrieval lives in plugins/xmatch/xmatch.py).

Flag bits (XM_FLAG):  1 = at least one counterpart within the radius   2 = more than one candidate within the radius (ambiguous)
                      4 = mutual nearest neighbours (the reference source's nearest object is this one too)   8 = the reference source is the nearest of another object as well
"""
import io
import math

import numpy as np

RA_NAMES = ('ra', 'ra_deg', 'raj2000', '_raj2000', 'alpha_j2000', 'ra_icrs', 'radeg', 'ra_j2000', 'ra2000', 'alpha', 'right_ascension', 'ra_obj', 'ra_', 's_ra')
DEC_NAMES = ('dec', 'dec_deg', 'dej2000', '_dej2000', 'delta_j2000', 'dec_icrs', 'decdeg', 'dec_j2000', 'dec2000', 'delta', 'declination', 'dec_obj', 'dec_', 's_dec', 'de_icrs', 'decj2000')


def find_columns(names, ra=None, dec=None):
    low = {n.lower(): n for n in names}
    r = ra if ra and ra in names else next((low[k] for k in RA_NAMES if k in low), None)
    d = dec if dec and dec in names else next((low[k] for k in DEC_NAMES if k in low), None)
    return r, d


def _to_float_array(col):
    a = np.asarray(col)
    if a.dtype.kind in 'fiu':
        return a.astype(float)
    if np.ma.isMaskedArray(col):
        a = np.ma.filled(col.astype(float), np.nan)
        return np.asarray(a, float)
    out = np.full(len(a), np.nan)
    for i, v in enumerate(a):
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            pass
    return out


def parse_sky(ra_col, dec_col):
    """-> (ra_deg, dec_deg).  Numeric columns are taken as degrees; strings like '13:29:52.7' / '+47:11:43' (hms/dms) are parsed with astropy."""
    ra_a = np.asarray(ra_col)
    if ra_a.dtype.kind in 'fiu' or np.ma.isMaskedArray(ra_col):
        return _to_float_array(ra_col), _to_float_array(dec_col)
    try:
        v = np.array([float(x) for x in ra_a[:5]])
        return _to_float_array(ra_col), _to_float_array(dec_col)
    except (TypeError, ValueError):
        pass
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    sc = SkyCoord([str(x) for x in ra_a], [str(x) for x in np.asarray(dec_col)], unit=(u.hourangle, u.deg))
    return sc.ra.deg, sc.dec.deg


def read_table(path, fmt='auto'):
    """-> dict name -> numpy array (original order kept in the key 'columns' of the returned (cols, data) pair).  FITS tables, VOTable (.xml/.vot/.votable), CSV, TSV/whitespace."""
    p = str(path)
    low = p.lower()
    from astropy.table import Table
    if fmt == 'auto':
        if low.endswith(('.fits', '.fit', '.fits.gz', '.fz')):
            fmt = 'fits'
        elif low.endswith(('.xml', '.vot', '.votable')):
            fmt = 'votable'
        elif low.endswith('.csv'):
            fmt = 'csv'
        else:
            fmt = 'ascii'
    if fmt == 'fits':
        t = Table.read(p, format='fits')
    elif fmt == 'votable':
        t = Table.read(p, format='votable')
    elif fmt == 'csv':
        t = Table.read(p, format='ascii.csv')
    else:
        with open(p) as fh:
            head = fh.readline()
        t = Table.read(p, format='ascii.tab' if '\t' in head else 'ascii')
    cols = list(t.colnames)
    return cols, {c: t[c] for c in cols}


def read_text_table(text, delimiter=','):
    """CSV/TSV text (e.g. a TAP response) -> (cols, data) with numeric columns converted."""
    import csv
    rd = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    rows = list(rd)
    cols = list(rd.fieldnames or [])
    data = {}
    for c in cols:
        vals = [r.get(c, '') for r in rows]
        try:
            data[c] = np.array([float(v) if v not in ('', None) else np.nan for v in vals])
        except ValueError:
            data[c] = np.array(vals, dtype=object)
    return cols, data


def sky_match(ra1, dec1, ra2, dec2, radius_arcsec):
    """All pairs within the radius.  Returns dict(i1, i2, sep (arcsec), dra, ddec (arcsec, reference - object, dra includes cos dec)).  Uses astropy's KD-tree search."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ok1 = np.isfinite(ra1) & np.isfinite(dec1); ok2 = np.isfinite(ra2) & np.isfinite(dec2)
    idx1 = np.nonzero(ok1)[0]; idx2 = np.nonzero(ok2)[0]
    if not len(idx1) or not len(idx2):
        e = np.zeros(0)
        return dict(i1=e.astype(int), i2=e.astype(int), sep=e, dra=e, ddec=e)
    c1 = SkyCoord(ra1[idx1] * u.deg, dec1[idx1] * u.deg)
    c2 = SkyCoord(ra2[idx2] * u.deg, dec2[idx2] * u.deg)
    j1, j2, sep, _ = c2.search_around_sky(c1, radius_arcsec * u.arcsec)
    sep = sep.arcsec
    dra = ((ra2[idx2[j2]] - ra1[idx1[j1]] + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec1[idx1[j1]])) * 3600.0
    ddec = (dec2[idx2[j2]] - dec1[idx1[j1]]) * 3600.0
    return dict(i1=idx1[j1], i2=idx2[j2], sep=sep, dra=dra, ddec=ddec)


def pixel_match(x1, y1, x2, y2, radius):
    from scipy.spatial import cKDTree
    ok1 = np.isfinite(x1) & np.isfinite(y1); ok2 = np.isfinite(x2) & np.isfinite(y2)
    idx1 = np.nonzero(ok1)[0]; idx2 = np.nonzero(ok2)[0]
    if not len(idx1) or not len(idx2):
        e = np.zeros(0)
        return dict(i1=e.astype(int), i2=e.astype(int), sep=e, dra=e, ddec=e)
    t2 = cKDTree(np.c_[x2[idx2], y2[idx2]])
    lists = t2.query_ball_point(np.c_[x1[idx1], y1[idx1]], radius)
    i1, i2 = [], []
    for k, l in enumerate(lists):
        for j in l:
            i1.append(idx1[k]); i2.append(idx2[j])
    i1 = np.array(i1, int); i2 = np.array(i2, int)
    dx = x2[i2] - x1[i1]; dy = y2[i2] - y1[i1]
    return dict(i1=i1, i2=i2, sep=np.hypot(dx, dy), dra=dx, ddec=dy)


def resolve(n1, pairs, n2=None):
    """Per object of catalog 1: nearest counterpart, number of candidates, flags; per object of catalog 2 bookkeeping for the mutual test.
    Returns dict(best (index into catalog 2 or -1), sep, dra, ddec, n, second_sep, flag) arrays of length n1."""
    i1, i2, sep = pairs['i1'], pairs['i2'], pairs['sep']
    best = np.full(n1, -1, int); bsep = np.full(n1, np.nan); n = np.zeros(n1, int); second = np.full(n1, np.nan)
    bdra = np.full(n1, np.nan); bddec = np.full(n1, np.nan)
    order = np.lexsort((sep, i1))
    seen = {}
    for k in order:
        a = i1[k]
        c = seen.get(a, 0)
        if c == 0:
            best[a] = i2[k]; bsep[a] = sep[k]; bdra[a] = pairs['dra'][k]; bddec[a] = pairs['ddec'][k]
        elif c == 1:
            second[a] = sep[k]
        seen[a] = c + 1
        n[a] = c + 1
    # nearest object of catalog 1 for every catalog-2 source (over the pairs)
    nearest1 = {}
    order2 = np.lexsort((sep, i2))
    for k in order2:
        j = i2[k]
        if j not in nearest1:
            nearest1[j] = i1[k]
    flag = np.zeros(n1, int)
    claims = {}
    for a in np.nonzero(best >= 0)[0]:
        claims[best[a]] = claims.get(best[a], 0) + 1
    for a in np.nonzero(best >= 0)[0]:
        f = 1
        if n[a] > 1:
            f |= 2
        if nearest1.get(best[a]) == a:
            f |= 4
        if claims[best[a]] > 1:
            f |= 8
        flag[a] = f
    return dict(best=best, sep=bsep, dra=bdra, ddec=bddec, n=n, second_sep=second, flag=flag)


def systematic_offset(res, clip=3.0, niter=5):
    """Median offset (dra, ddec) of the unique mutual matches, sigma-clipped; returns (dra, ddec, n, scatter_dra, scatter_ddec)."""
    sel = (res['flag'] & 1 > 0) & ((res['flag'] & 2) == 0) & ((res['flag'] & 4) > 0)
    x, y = res['dra'][sel], res['ddec'][sel]
    if len(x) < 5:
        return 0.0, 0.0, len(x), float('nan'), float('nan')
    m = np.ones(len(x), bool)
    for _ in range(niter):
        mx, my = np.median(x[m]), np.median(y[m])
        sx = 1.4826 * np.median(np.abs(x[m] - mx)); sy = 1.4826 * np.median(np.abs(y[m] - my))
        r = np.hypot((x - mx) / max(sx, 1e-9), (y - my) / max(sy, 1e-9))
        mm = r < clip
        if mm.sum() < 5 or np.array_equal(mm, m):
            break
        m = mm
    return float(np.median(x[m])), float(np.median(y[m])), int(m.sum()), float(sx), float(sy)


def chance_matches(ra1, dec1, ra2, dec2, radius_arcsec, offsets_arcsec=(30.0, 45.0, 60.0, 90.0, 120.0), n_shifts=None):
    """Expected number of chance coincidences within the radius, estimated by shifting catalog 2 on the sky by offsets much larger than the radius (8 directions each).
    Returns dict(mean_matched_objects, per_shift [..]) = number of catalog-1 objects with at least one counterpart in the shifted catalog."""
    counts = []
    for off in offsets_arcsec:
        for ang in np.arange(0, 360, 45.0):
            dd = off * math.sin(math.radians(ang)) / 3600.0
            dr = off * math.cos(math.radians(ang)) / 3600.0
            cosd = np.maximum(np.cos(np.radians(dec2)), 0.05)
            r2 = (ra2 + dr / cosd) % 360.0; d2 = dec2 + dd
            p = sky_match(ra1, dec1, r2, d2, radius_arcsec)
            # edge correction: objects of catalog 1 outside the footprint of the shifted catalog 2 cannot have chance partners
            dra = (ra1 - np.nanmedian(r2) + 180.0) % 360.0 - 180.0; rr = (r2 - np.nanmedian(r2) + 180.0) % 360.0 - 180.0
            inside = (dra >= np.nanmin(rr)) & (dra <= np.nanmax(rr)) & (dec1 >= np.nanmin(d2)) & (dec1 <= np.nanmax(d2))
            f = float(np.mean(inside[np.isfinite(ra1)])) if np.isfinite(ra1).any() else 1.0
            c = len(np.unique(p['i1']))
            counts.append(c / f if f > 0.3 else c)
    return dict(mean=float(np.mean(counts)), std=float(np.std(counts)), counts=[round(float(c), 2) for c in counts])


def crossmatch(ra1, dec1, ra2, dec2, radius_arcsec, apply_shift=False, pixel=False, x1=None, y1=None, x2=None, y2=None):
    """High-level: match, resolve, optional systematic shift (matches redone after removing the median offset), chance estimate (sky only).
    Returns dict(res (resolve output), shift (dra, ddec), n_matched, chance)."""
    n1 = len(ra1 if not pixel else x1)
    pairs = pixel_match(x1, y1, x2, y2, radius_arcsec) if pixel else sky_match(ra1, dec1, ra2, dec2, radius_arcsec)
    res = resolve(n1, pairs)
    shift = systematic_offset(res)
    out_shift = (0.0, 0.0)
    if apply_shift and shift[2] >= 5 and not pixel:
        cosd = np.maximum(np.cos(np.radians(dec2)), 0.05)
        ra2s = (ra2 - shift[0] / 3600.0 / cosd) % 360.0
        dec2s = dec2 - shift[1] / 3600.0
        pairs = sky_match(ra1, dec1, ra2s, dec2s, radius_arcsec)
        res = resolve(n1, pairs)
        out_shift = (shift[0], shift[1])
    chance = None if pixel else chance_matches(ra1, dec1, ra2, dec2, radius_arcsec)
    return dict(res=res, shift=dict(dra=shift[0], ddec=shift[1], n=shift[2], scatter_dra=shift[3], scatter_ddec=shift[4], applied=out_shift), n_matched=int(np.sum(res['best'] >= 0)), chance=chance)
