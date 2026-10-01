"""Response mapping: service JSON -> {column: value} per object (with type conversion and units)."""
from . import contracts, units
from .jsonpath import MISSING, convert, get_path


class MapError(Exception):
    pass


def default_fields(task, items):
    """Fields used when the profile has no response.fields (local_command / python_callable / mock):
    the task's standard columns that occur in the results; for 'generic' every key except reserved ones."""
    reserved = {'id', 'error', 'model', 'request_id'}
    keys = []
    for it in items:
        if isinstance(it, dict):
            for k in it:
                if k not in keys:
                    keys.append(k)
    if task in contracts.TASKS and task != 'generic':
        names = contracts.output_names(task)
        sel = [n for n in names if n in keys]
        return {n: {'path': n, 'type': contracts.output_spec(task, n)['type']} for n in sel}
    sel = sorted(k for k in keys if k not in reserved)
    out = {}
    for k in sel:
        v = next((it[k] for it in items if isinstance(it, dict) and it.get(k) is not None), None)
        t = 'float' if isinstance(v, float) else 'int' if isinstance(v, int) and not isinstance(v, bool) else \
            'bool' if isinstance(v, bool) else 'json' if isinstance(v, (dict, list)) else 'str'
        out[k] = {'path': k, 'type': t}
    return out


def apply_fields(item, fields, root=None):
    """item: parsed JSON (dict) for one object.  Paths starting with '$root.' (or '$.') are looked up in the whole
    response `root` instead.  Returns {col: value}; raises MapError for unconvertible values / unit problems."""
    out = {}
    found = 0
    for col, spec in fields.items():
        path = spec['path']
        src = item
        if path.startswith('$root.'):
            src, path = root, path[len('$root.'):]
        raw = get_path(src, path)
        if raw is MISSING:
            out[col] = None
            continue
        found += 1
        try:
            v = convert(raw, spec)
            if v is not None and spec.get('unit') and spec.get('to_unit') and spec.get('type', 'str') in ('float', 'int'):
                v = v * units.factor(spec['unit'], spec['to_unit'])
        except ValueError as e:
            raise MapError('field %s (path %s): %s' % (col, spec['path'], e))
        if v is not None and col in contracts.UNIT_INTERVAL and isinstance(v, (int, float)) and not 0.0 <= v <= 1.0:
            raise MapError('field %s = %r outside [0, 1] (contract); check the profile mapping' % (col, v))
        out[col] = v
    if fields and found == 0:
        raise MapError('none of the mapped paths exist in the response (%s)' % ', '.join(s['path'] for s in fields.values()))
    return out
