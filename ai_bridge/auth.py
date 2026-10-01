"""Authentication from environment variables.  The profile stores only the NAME of the variable.
Values are read at send time, never printed, never cached, never written to provenance; any occurrence
in error text is replaced by ***.
"""
import os

from .errors import AuthError


class Auth:
    def __init__(self, spec):
        spec = spec or {}
        self.scheme = spec.get('scheme', 'none')
        self.env = spec.get('env')
        self.header = spec.get('header')
        self.param = spec.get('param')
        self.prefix = spec.get('prefix', 'Bearer ')

    def is_set(self):
        return self.scheme == 'none' or bool(self.env and os.environ.get(self.env))

    def _value(self):
        v = os.environ.get(self.env or '')
        if not v:
            raise AuthError('environment variable %s (auth.%s) is not set' % (self.env, self.scheme))
        return v

    def apply(self, headers, query):
        """Returns (headers, query) with credentials added (new dicts)."""
        headers, query = dict(headers), dict(query)
        if self.scheme == 'none':
            return headers, query
        v = self._value()
        if self.scheme == 'bearer_env':
            headers['Authorization'] = self.prefix + v
        elif self.scheme == 'header_env':
            headers[self.header] = v
        elif self.scheme == 'query_env':
            query[self.param] = v
        return headers, query

    def describe(self):
        """Redacted description for dry-run output: never the value."""
        if self.scheme == 'none':
            return {}, {}
        st = '<value of $%s: %s>' % (self.env, 'set' if os.environ.get(self.env or '') else 'NOT SET')
        if self.scheme == 'bearer_env':
            return {'Authorization': self.prefix + st}, {}
        if self.scheme == 'header_env':
            return {self.header: st}, {}
        return {}, {self.param: st}

    def scrub(self, text):
        v = os.environ.get(self.env or '') if self.scheme != 'none' else None
        if v and len(v) >= 4:
            text = text.replace(v, '***')
        return text
