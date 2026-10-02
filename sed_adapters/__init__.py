"""Adapters to external photo-z / SED-fitting codes: EAZY (photoz), CIGALE, Bagpipes, Prospector (sed_fit).  See docs/sed_codes.md.

Each adapter module has process(records, params, task) -> (results, model_string), run(records, params, context) (ai_bridge python_callable) and a
`python -m sed_adapters.<code>_adapter` entry (ai_bridge local_command)."""
import importlib

CODES = ('eazy', 'cigale', 'bagpipes', 'prospector')


class _Lazy(dict):
    def __missing__(self, key):
        if key not in CODES:
            raise KeyError(key)
        m = importlib.import_module('%s.%s_adapter' % (__name__, key))
        self[key] = m
        return m


ADAPTERS = _Lazy()
