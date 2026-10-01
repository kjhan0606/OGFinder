"""Catalog TSV <-> standard input records, and output TSV writing."""
import csv
import io
import json
import math
import re

from . import contracts

BAD_MAG = 90.0   # SExtractor 99 / "no measurement" convention (>= this means "none")


def read_catalog(path):
    """Returns (header, rows) where rows are lists of str.  Lines starting with '#' are ignored."""
    with open(path, 'r', encoding='utf-8', newline='') as f:
        txt = f.read()
    lines = [l for l in txt.split('\n') if l.strip() != '' and not l.startswith('#')]
    if not lines:
        raise ValueError('catalog %s is empty' % path)
    header = lines[0].split('\t')
    rows = [l.split('\t') for l in lines[1:]]
    return header, rows


def _num(s):
    try:
        v = float(str(s).strip())
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def detect_mag_columns(header):
    """MAG_<band> columns -> {band: (mag_col, err_col or None)}.  MAG_APER_n / MAG_ISOCOR etc. are skipped
    unless explicitly requested with --mag-columns."""
    skip = re.compile(r'^MAG_(APER|ISO|ISOCOR|BEST|MODEL|PETRO|PSF|CROWD)(_\d+)?$')
    out = {}
    for c in header:
        if c.startswith('MAG_') and not skip.match(c):
            b = c[4:]
            ec = 'MAGERR_' + b
            out[b] = (c, ec if ec in header else None)
    return out


def build_records(header, rows, id_col='NUMBER', x_col='X_IMAGE', y_col='Y_IMAGE', ra_col='ALPHA_J2000',
                  dec_col='DELTA_J2000', mag_columns=None, pixel_scale=None, psf_fwhm=None, numbers=None):
    """Standard input records (dicts), in catalogue order.  `numbers` restricts to those ids."""
    idx = {c: i for i, c in enumerate(header)}
    if id_col not in idx:
        raise ValueError('catalog has no %s column (needed to merge results back)' % id_col)
    if mag_columns is None:
        mc = detect_mag_columns(header)
    else:
        mc = {}
        for c in mag_columns:
            if c not in idx:
                raise ValueError('mag column %s not in catalog' % c)
            b = c[4:] if c.startswith('MAG_') else c
            ec = ('MAGERR_' + b) if ('MAGERR_' + b) in idx else None
            mc[b] = (c, ec)
    want = None if numbers is None else set(str(n) for n in numbers)
    recs = []
    seen = set()
    for r in rows:
        r = r + [''] * (len(header) - len(r))
        oid = r[idx[id_col]].strip()
        if want is not None and oid not in want:
            continue
        if oid in seen:
            raise ValueError('duplicate %s value %r in catalog' % (id_col, oid))
        seen.add(oid)
        mags, errs = {}, {}
        for b, (c, ec) in mc.items():
            m = _num(r[idx[c]])
            mags[b] = None if (m is None or m >= BAD_MAG) else m
            e = _num(r[idx[ec]]) if ec else None
            errs[b] = None if (mags[b] is None or e is None or e >= BAD_MAG) else e
        recs.append({
            'id': oid,
            'x': _num(r[idx[x_col]]) if x_col in idx else None,
            'y': _num(r[idx[y_col]]) if y_col in idx else None,
            'ra': _num(r[idx[ra_col]]) if ra_col in idx else None,
            'dec': _num(r[idx[dec_col]]) if dec_col in idx else None,
            'mags': mags, 'mag_errs': errs,
            'pixel_scale_arcsec': pixel_scale, 'psf_fwhm_arcsec': psf_fwhm,
            'catalog_row': {h: r[i] for i, h in enumerate(header)},
        })
    return recs


def fmt_value(v, typ):
    if v is None:
        return ''
    if typ == 'float':
        return repr(float(v))
    if typ == 'int':
        return str(int(v))
    if typ == 'bool':
        return 'True' if v else 'False'
    s = str(v)
    return s.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')


def output_columns(task, mapped_fields, suffix=''):
    """[(column_name, type)] in output order: the mapped task outputs, then provenance columns (all text)."""
    cols = [(n + suffix, t) for n, t in mapped_fields]
    cols += [(c, 'str') for c in contracts.prov_columns(task)]
    return cols


def write_tsv(fileobj, id_col, ids, columns, table):
    """table: {id: {column: value}}; deterministic order = ids order."""
    fileobj.write('\t'.join([id_col] + [c for c, _ in columns]) + '\n')
    for i in ids:
        row = table.get(i, {})
        fileobj.write('\t'.join([i] + [fmt_value(row.get(c), t) for c, t in columns]) + '\n')


def parse_tsv_table(text):
    """-> (header, {id: {col: str}}) with the first column as id."""
    lines = [l for l in text.split('\n') if l.strip() != '' and not l.startswith('#')]
    if not lines:
        return [], {}
    h = lines[0].split('\t')
    out = {}
    for l in lines[1:]:
        f = l.split('\t')
        f += [''] * (len(h) - len(f))
        out[f[0]] = dict(zip(h[1:], f[1:]))
    return h, out
