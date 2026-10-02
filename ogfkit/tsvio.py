"""Catalog TSV helpers (the OGFinder catalog format: header row, tab separated, NUMBER is the key column)."""
import math
import sys


def read_catalog(path):
    """-> (columns, rows) with rows = list of dict col -> str.  Comment lines (#...) and blank lines are skipped."""
    with open(path, encoding='utf-8') as f:
        lines = [l.rstrip('\n') for l in f if l.strip() and not l.startswith('#')]
    if not lines:
        return [], []
    cols = [c.strip() for c in lines[0].split('\t')]
    rows = []
    for l in lines[1:]:
        v = l.split('\t')
        v += [''] * (len(cols) - len(v))
        rows.append({c: v[i].strip() for i, c in enumerate(cols)})
    return cols, rows


def fnum(s, default=float('nan')):
    try:
        v = float(s)
        return v
    except (TypeError, ValueError):
        return default


def fmt(v, nd=6):
    """Cell text: '' for None/NaN, shortest sensible float otherwise."""
    if v is None:
        return ''
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return ''
        return ('%.' + str(nd) + 'g') % v
    return str(v)


def write_columns(cols, rows, fh=None):
    """Print 'NUMBER <cols>' TSV (the add_columns contract of a CLI step).  rows = list of (number, dict)."""
    fh = fh or sys.stdout
    fh.write('\t'.join(['NUMBER'] + list(cols)) + '\n')
    for num, d in rows:
        fh.write('\t'.join([str(num)] + [fmt(d.get(c)) for c in cols]) + '\n')


def write_table(path, cols, records):
    """Write a long-form TSV (records = list of dict)."""
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for r in records:
            f.write('\t'.join(fmt(r.get(c)) for c in cols) + '\n')


def parse_numbers(text):
    """'1,5,7-9' -> [1,5,7,8,9]; '' -> []"""
    out = []
    for t in (text or '').replace(' ', '').split(','):
        if not t:
            continue
        if '-' in t[1:]:
            a, b = t.split('-', 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(t))
    return out
