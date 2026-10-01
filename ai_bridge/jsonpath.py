"""Tiny JSONPath-like accessor and type conversion for response mapping.

Path syntax:  a.b[0].c   a.items[-1]   a.rows[*].id (wildcard -> list)   optional leading "$."
"""
import math
import re

MISSING = object()
_TOK = re.compile(r'\[(-?\d+|\*)\]|([^.\[\]]+)')


def tokens(path):
    p = path.strip()
    if p.startswith('$'):
        p = p[1:].lstrip('.')
    return [(m.group(1) if m.group(1) is not None else m.group(2), m.group(1) is not None)
            for m in _TOK.finditer(p)]


def get_path(obj, path):
    """Return the value at `path` or MISSING."""
    cur = [obj]
    wild = False
    for tok, is_idx in tokens(path):
        nxt = []
        for c in cur:
            if is_idx:
                if not isinstance(c, list):
                    continue
                if tok == '*':
                    wild = True
                    nxt.extend(c)
                else:
                    i = int(tok)
                    if -len(c) <= i < len(c):
                        nxt.append(c[i])
            else:
                if isinstance(c, dict) and tok in c:
                    nxt.append(c[tok])
        cur = nxt
        if not cur:
            return MISSING
    if wild:
        return cur
    return cur[0] if cur else MISSING


def _is_missing_number(v, spec):
    return v in (spec.get('missing_values') or [])


def convert(value, spec):
    """Convert a raw response value according to a field spec.

    spec: type float|int|str|bool|json, scale, offset, missing_values, join (for lists).
    Returns None for missing/NaN.  Raises ValueError for unconvertible values.
    """
    t = spec.get('type', 'str')
    if value is MISSING or value is None:
        return None
    if isinstance(value, (list, tuple)) and t in ('str',):
        value = spec.get('join', ',').join('' if x is None else str(x) for x in value)
    if t in ('float', 'int'):
        if isinstance(value, bool):
            raise ValueError('boolean where a number is expected')
        if isinstance(value, str):
            s = value.strip()
            if s == '' or s.lower() in ('nan', 'null', 'none'):
                return None
            value = float(s)
        if not isinstance(value, (int, float)):
            raise ValueError('not a number: %r' % (value,))
        if _is_missing_number(value, spec):
            return None
        v = float(value)
        if math.isnan(v) or math.isinf(v):
            return None
        v = v * float(spec.get('scale', 1.0)) + float(spec.get('offset', 0.0))
        return int(round(v)) if t == 'int' else v
    if t == 'bool':
        if isinstance(value, str):
            s = value.strip().lower()
            if s in ('true', '1', 'yes', 't', 'y'):
                return True
            if s in ('false', '0', 'no', 'f', 'n'):
                return False
            raise ValueError('not a boolean: %r' % value)
        return bool(value)
    if t == 'json':
        import json
        return json.dumps(value, sort_keys=True, separators=(',', ':'))
    if isinstance(value, float) and math.isnan(value):
        return None
    return str(value)
