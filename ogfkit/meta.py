"""Catalog metadata: a small JSON dictionary kept next to the working catalog (`catalog_meta.json` in the work directory).

Plugins add keys with update(path, key, value) (e.g. "completeness": {...}); the GUI's Save Catalog writes `<name>.meta.json` and the
FITS table export puts the scalar values into the table header, so that limits such as the 50 % completeness magnitude travel with the catalog.
`nrows` records the catalog size the metadata belongs to; consumers ignore metadata whose nrows differs from the catalog's."""
import json
import os


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def update(path, key, value, nrows=None):
    d = load(path)
    d[key] = value
    if nrows is not None:
        d['nrows'] = int(nrows)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(d, f, indent=1, sort_keys=True)
    os.replace(tmp, path)
    return d


def flat_scalars(d, prefix=''):
    """Flatten nested dicts to {'completeness.lim50': 26.1, ...} keeping only scalar leaves (for FITS header keywords)."""
    out = {}
    for k, v in d.items():
        kk = prefix + k
        if isinstance(v, dict):
            out.update(flat_scalars(v, kk + '.'))
        elif isinstance(v, (int, float, str, bool)) and not (isinstance(v, float) and v != v):
            out[kk] = v
    return out
