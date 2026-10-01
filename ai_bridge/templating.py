"""Request templating: {placeholders} inside JSON-like templates.

* a string that is exactly "{name}" is replaced by the (typed) value: numbers stay numbers,
  dicts/lists stay JSON structures, None becomes null.
* a placeholder inside a longer string is replaced by its text (dict/list -> compact JSON).
* "{{" and "}}" are literal braces.
* an unknown placeholder raises TemplateError (nothing is sent).
"""
import json
import re

from .errors import TemplateError

_PH = re.compile(r'\{([A-Za-z_][A-Za-z0-9_]*)\}')
_LB, _RB = '\x00LB\x00', '\x00RB\x00'


def placeholders(obj):
    """All placeholder names used in a template (any JSON-like structure)."""
    out = []

    def walk(o):
        if isinstance(o, str):
            for m in _PH.finditer(o.replace('{{', _LB).replace('}}', _RB)):
                if m.group(1) not in out:
                    out.append(m.group(1))
        elif isinstance(o, dict):
            for k, v in o.items():
                walk(k)
                walk(v)
        elif isinstance(o, (list, tuple)):
            for v in o:
                walk(v)
    walk(obj)
    return out


def _text(v):
    if v is None:
        return ''
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, (dict, list)):
        return json.dumps(v, sort_keys=True, separators=(',', ':'))
    if isinstance(v, float):
        return repr(v)
    return str(v)


def render(obj, lookup):
    """lookup(name) -> value, raises KeyError for unknown names."""
    def rstr(s):
        t = s.replace('{{', _LB).replace('}}', _RB)
        m = _PH.fullmatch(t)
        if m:
            try:
                return lookup(m.group(1))
            except KeyError:
                raise TemplateError('unknown placeholder {%s}' % m.group(1))

        def rep(mm):
            try:
                return _text(lookup(mm.group(1)))
            except KeyError:
                raise TemplateError('unknown placeholder {%s}' % mm.group(1))
        return _PH.sub(rep, t).replace(_LB, '{').replace(_RB, '}')

    def walk(o):
        if isinstance(o, str):
            r = rstr(o)
            return r.replace(_LB, '{').replace(_RB, '}') if isinstance(r, str) else r
        if isinstance(o, dict):
            return {walk(k) if isinstance(k, str) else k: walk(v) for k, v in o.items()}
        if isinstance(o, list):
            return [walk(v) for v in o]
        return o
    return walk(obj)
