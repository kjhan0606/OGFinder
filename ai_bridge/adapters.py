"""Built-in minimal adapters (examples / templates, not model implementations):

  https_json / https_multipart   -> "generic_rest"
  local_command                  JSON over stdin/stdout
  python_callable                module:function
  mock                           deterministic FAKE values (SERVICE=mock)
  tap_query                      IVOA TAP ADQL cone search by object position (non-ML)

Adapter protocol (used by runner.py):
    a.cutout_needs()            -> list of (fmt, band_or_None, kind) the templates use (or None = profile-driven)
    a.prepare(batch, cut)       -> Prepared   (no side effects; `cut(rec, fmt, band)` returns cutout info or raises)
    a.send(prepared, timeout)   -> (obj, meta)  raises errors.RequestFailed
    a.extract(obj, meta, batch) -> {id: {'values': {col: v}} | {'error': str}}, model, request_id
"""
import base64
import csv
import hashlib
import io
import json
import mimetypes
import os
import re
import subprocess
import sys
import urllib.parse

from . import CONTRACT_VERSION, contracts, httpio, mapping, templating
from .auth import Auth
from .errors import ProfileError, RequestFailed, TemplateError

CUT_RE = re.compile(r'^cutout(?:_(png|fits|npy))?_(base64|path|url)(?:_(.+))?$')
MOCK_TIME = '1970-01-01T00:00:00Z'


class Prepared:
    def __init__(self, ids, key_material, preview, size_bytes=0, payload=None):
        self.ids = ids
        self.key_material = key_material
        self.preview = preview
        self.size_bytes = size_bytes
        self.payload = payload


def public_record(rec, cutouts=None):
    """The standard input record as sent to local commands / callables (cutouts as file descriptors)."""
    d = {k: rec.get(k) for k in ('id', 'x', 'y', 'ra', 'dec', 'mags', 'mag_errs', 'pixel_scale_arcsec',
                                 'psf_fwhm_arcsec', 'catalog_row')}
    d['cutouts'] = cutouts or {}
    return d


def _cutout_digest(c):
    return {k: {'sha256': v.get('sha256'), 'format': v.get('format'), 'size': v.get('size')} for k, v in (c or {}).items()}


def short_base64(s, n=24):
    return '<base64 %d chars, sha256 %s>' % (len(s), hashlib.sha256(s.encode()).hexdigest()[:12])


def redact_for_display(o):
    """Shorten big base64 blobs for dry-run output."""
    if isinstance(o, dict):
        return {k: redact_for_display(v) for k, v in o.items()}
    if isinstance(o, list):
        return [redact_for_display(v) for v in o]
    if isinstance(o, str) and len(o) > 400 and re.fullmatch(r'[A-Za-z0-9+/=\r\n]+', o):
        return short_base64(o)
    return o


class Base:
    def __init__(self, profile, ctx):
        self.p = profile
        self.ctx = ctx                    # dict: params, bands, cutout_cfg, mock flags, ...
        self.name = profile['name']
        self.task = profile.get('task')
        self.fields = (profile.get('response') or {}).get('fields')

    def cutout_needs(self):
        return None

    def check(self):
        return True, 'nothing to check'


# ---------------------------------------------------------------- lookup for templates
class Lookup:
    """Placeholder resolver for one record (batch_size 1) or for item rendering."""

    def __init__(self, adapter, rec, cut, extra=None):
        self.a, self.rec, self.cut, self.extra = adapter, rec, cut, extra or {}

    def __call__(self, name):
        r, p = self.rec, self.a.p
        if name in self.extra:
            return self.extra[name]
        simple = ('id', 'x', 'y', 'ra', 'dec', 'mags', 'mag_errs', 'catalog_row', 'pixel_scale_arcsec', 'psf_fwhm_arcsec')
        if name in simple:
            return r.get(name)
        if name in ('mags_json', 'mag_errs_json', 'catalog_row_json'):
            return json.dumps(r.get(name[:-5]), sort_keys=True, separators=(',', ':'))
        if name == 'task':
            return self.a.task
        if name == 'service':
            return self.a.name
        if name == 'contract':
            return CONTRACT_VERSION
        if name in ('bands', 'bands_json'):
            b = list(self.a.ctx.get('bands') or [])
            return b if name == 'bands' else json.dumps(b)
        if name == 'radius_deg':
            ra = self.a.ctx['params'].get('radius_arcsec')
            if ra is None:
                raise KeyError(name)
            return float(ra) / 3600.0
        if name.startswith('param_'):
            k = name[6:]
            if k in self.a.ctx['params']:
                return self.a.ctx['params'][k]
            raise KeyError(name)
        if name == 'cutouts_json':
            fmt = (p.get('cutouts') or {}).get('format', 'png')
            out = {}
            for b in self.a.ctx.get('bands') or []:
                c = self.cut(r, fmt, b)
                out[b] = {'format': fmt, 'size': c['size'], 'base64': _b64(c['path'])}
            return json.dumps(out, sort_keys=True, separators=(',', ':'))
        m = CUT_RE.match(name)
        if m:
            fmt, kind, band = m.group(1), m.group(2), m.group(3)
            fmt = fmt or (p.get('cutouts') or {}).get('format', 'png')
            if kind == 'url':
                return self.a.cutout_url(r, band)
            band = band or (self.a.ctx.get('bands') or [None])[0]
            c = self.cut(r, fmt, band)
            return _b64(c['path']) if kind == 'base64' else c['path']
        raise KeyError(name)


def _b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('ascii')


def template_cutout_needs(tpl, default_fmt, default_band):
    out = []
    for ph in templating.placeholders(tpl):
        if ph == 'cutouts_json':
            out.append((default_fmt, '*', 'base64'))
            continue
        m = CUT_RE.match(ph)
        if m and m.group(2) != 'url':
            out.append((m.group(1) or default_fmt, m.group(3) or default_band, m.group(2)))
    return out


# ---------------------------------------------------------------- generic REST
class GenericREST(Base):
    """https_json (JSON body or query string) and https_multipart (image upload)."""

    def __init__(self, profile, ctx):
        super().__init__(profile, ctx)
        self.auth = Auth(profile.get('auth'))
        self.multipart = profile['transport'] == 'https_multipart'
        self.req = profile.get('request') or {}
        self.resp = profile.get('response') or {}

    def _templates(self):
        r = self.req
        if self.multipart:
            return [r.get('multipart'), r.get('query')]
        return [r.get('template'), r.get('item_template'), r.get('query')]

    def cutout_needs(self):
        fmt = (self.p.get('cutouts') or {}).get('format', 'png')
        out = []
        for t in self._templates():
            if t is not None:
                out += template_cutout_needs(t, fmt, None)
        return out

    def cutout_url(self, rec, band):
        t = (self.p.get('cutouts') or {}).get('url_template')
        if not t:
            raise KeyError('cutout_url (profile has no cutouts.url_template)')
        cu = self.p['cutouts']
        vals = {'id': rec['id'], 'ra': rec.get('ra'), 'dec': rec.get('dec'), 'band': band or '',
                'size_arcsec': cu.get('size_arcsec'), 'size_pix': cu.get('size_pix')}

        def lk(n):
            if n in vals:
                if vals[n] is None:
                    raise KeyError(n)
                return vals[n]
            raise KeyError(n)
        return templating.render(t, lk)

    def _url(self):
        base = self.p['base_url'].rstrip('/')
        ep = (self.p.get('endpoint') or '').lstrip('/')
        return base + ('/' + ep if ep else '')

    def prepare(self, batch, cut):
        p = self.p
        method = p.get('method', 'POST')
        headers = {'Accept': 'application/json'}
        headers.update(p.get('headers') or {})
        ids = [r['id'] for r in batch]
        query = {}
        if self.req.get('query'):
            query = templating.render(self.req['query'], Lookup(self, batch[0], cut)) if len(batch) == 1 else \
                templating.render(self.req['query'], Lookup(self, batch[0], cut))
            query = {k: ('' if v is None else v) for k, v in query.items()}
        body, size = None, 0
        material = {'method': method, 'url': self._url(), 'headers': dict(headers), 'query': dict(query),
                    'contract': CONTRACT_VERSION, 'task': self.task}
        preview = {'method': method, 'url': self._url(), 'query': dict(query), 'ids': ids}
        if self.multipart:
            rec = batch[0]
            lk = Lookup(self, rec, cut)
            mp = self.req['multipart']
            fields = templating.render(mp.get('fields', {}), lk)
            files = []
            fdesc = {}
            for fname, tpl in mp['files'].items():
                path = templating.render(tpl, lk)
                with open(path, 'rb') as f:
                    data = f.read()
                ct = mimetypes.guess_type(path)[0] or 'application/octet-stream'
                files.append((fname, os.path.basename(path), ct, data))
                fdesc[fname] = {'filename': os.path.basename(path), 'content_type': ct, 'bytes': len(data),
                                'sha256': hashlib.sha256(data).hexdigest()}
            fields = {k: ('' if v is None else v if isinstance(v, str) else json.dumps(v)) for k, v in fields.items()}
            ctype, body = httpio.multipart(fields, files)
            headers['Content-Type'] = ctype
            size = len(body)
            material.update(fields=fields, files={k: {x: v[x] for x in ('filename', 'sha256')} for k, v in fdesc.items()})
            preview.update(content_type='multipart/form-data', fields=fields, files=fdesc)
        elif method in ('POST', 'PUT'):
            tpl = self.req.get('template', {})
            if len(batch) == 1 and 'items' not in templating.placeholders(tpl):
                payload = templating.render(tpl, Lookup(self, batch[0], cut))
            else:
                item_t = self.req.get('item_template') or {'id': '{id}', 'ra': '{ra}', 'dec': '{dec}', 'mags': '{mags}'}
                items = [templating.render(item_t, Lookup(self, r, cut)) for r in batch]
                payload = templating.render(tpl, Lookup(self, batch[0], cut, extra={'items': items, 'n_items': len(items),
                                                                                  'batch_ids': ids}))
            body = json.dumps(payload, sort_keys=False, separators=(',', ':'), ensure_ascii=True).encode()
            headers.setdefault('Content-Type', 'application/json')
            size = len(body)
            material['body'] = payload
            preview['body'] = redact_for_display(payload)
        preview['bytes'] = size
        # headers shown in dry-run (credentials redacted); real headers are built in send()
        h2, q2 = self.auth.describe()
        preview['headers'] = dict(headers, **h2)
        if q2:
            preview['query'] = dict(preview['query'], **q2)
        return Prepared(ids, material, preview, size, payload=(method, headers, query, body))

    def send(self, prepared, timeout):
        method, headers, query, body = prepared.payload
        headers, query = self.auth.apply(headers, query)      # raises AuthError if env var is unset
        url = self._url()
        if query:
            url += ('&' if '?' in url else '?') + urllib.parse.urlencode(query)
        try:
            status, rh, data = httpio.request(method, url, headers, body, timeout)
        except RequestFailed as e:
            raise RequestFailed(self.auth.scrub(str(e)), e.retryable, e.status, e.retry_after)
        try:
            obj = json.loads(data.decode('utf-8'))
        except ValueError:
            raise RequestFailed('response is not JSON (HTTP %d, %d bytes)' % (status, len(data)), retryable=False)
        meta = {'http_status': status}
        for k, hk in (('request_id', 'request_id_header'), ('model', 'model_header')):
            if self.resp.get(hk):
                v = rh.get(self.resp[hk].lower())
                if v:
                    meta[k] = v
        return obj, meta

    def extract(self, obj, meta, batch):
        r = self.resp
        ids = [x['id'] for x in batch]
        model = meta.get('model') or _str_or_none(_first(obj, r.get('model_path')))
        rid = meta.get('request_id') or _str_or_none(_first(obj, r.get('request_id_path')))
        res = {}
        fields = self.fields
        items_path = r.get('items_path')
        if items_path:
            items = _first(obj, items_path)
            if not isinstance(items, list):
                return {i: {'error': 'response.items_path %r is not a list in the response' % items_path} for i in ids}, model, rid
            by_id = {}
            if r.get('id_path'):
                for it in items:
                    k = _first(it, r['id_path'])
                    if k is not None:
                        by_id[str(k)] = it
            else:
                if len(items) != len(batch):
                    return {i: {'error': 'response has %d items for %d objects and no id_path to match them' % (len(items), len(batch))}
                            for i in ids}, model, rid
                by_id = {i: it for i, it in zip(ids, items)}
            for i in ids:
                res[i] = self._one(by_id.get(i, mapping.MISSING if False else None), obj, fields, missing_msg='object id %s not in response' % i)
        else:
            if len(batch) != 1:
                return {i: {'error': 'batch_size > 1 needs response.items_path'} for i in ids}, model, rid
            res[ids[0]] = self._one(obj, obj, fields)
        return res, model, rid

    def _one(self, item, root, fields, missing_msg=None):
        if item is None:
            return {'error': missing_msg or 'empty response'}
        ep = self.resp.get('error_path')
        if ep:
            e = _first(item, ep)
            if e not in (None, '', False):
                return {'error': 'service error: %s' % (e if isinstance(e, str) else json.dumps(e))[:300]}
        try:
            return {'values': mapping.apply_fields(item, fields, root)}
        except mapping.MapError as e:
            return {'error': 'mapping: %s' % e}

    def check(self):
        h = self.p.get('health') or {}
        path = h.get('path')
        url = self.p['base_url'].rstrip('/') + ('/' + path.lstrip('/') if path else '')
        headers, q = {}, {}
        if h.get('send_auth', True):
            try:
                headers, q = self.auth.apply({}, {})
            except Exception as e:
                return False, str(e)
        if q:
            url += '?' + urllib.parse.urlencode(q)
        try:
            st, rh, data = httpio.request('GET', url, headers, None, float(self.p.get('timeout_s', 30)))
        except RequestFailed as e:
            return False, self.auth.scrub(str(e))
        return True, 'GET %s -> HTTP %d (%d bytes)' % (url.split('?')[0], st, len(data))


def _first(obj, path):
    if not path:
        return None
    from .jsonpath import MISSING, get_path
    v = get_path(obj, path)
    return None if v is MISSING else v


def _str_or_none(v):
    return None if v is None else str(v)


# ---------------------------------------------------------------- local command / python callable
class _RecordAdapter(Base):
    """Common part of local_command / python_callable / mock: sends a list of standard records."""

    def cutout_needs(self):
        cu = self.p.get('cutouts') or {}
        if cu.get('enabled') or cu.get('use'):
            return [(cu.get('format', 'png'), b, 'path') for b in (self.ctx.get('bands') or [None])]
        return []

    def _records(self, batch, cut):
        out = []
        needs = self.cutout_needs()
        for r in batch:
            c = {}
            for fmt, band, kind in needs:
                info = cut(r, fmt, band)
                c[band or 'main'] = {k: info[k] for k in ('path', 'format', 'size', 'sha256', 'wcs')}
            out.append(public_record(r, c))
        return out

    def prepare(self, batch, cut):
        recs = self._records(batch, cut)
        req = {'contract': CONTRACT_VERSION, 'task': self.task, 'service': self.name,
               'params': self.ctx['params'], 'records': recs}
        material = {'contract': CONTRACT_VERSION, 'task': self.task, 'params': self.ctx['params'],
                    'records': [dict(r, cutouts=_cutout_digest(r['cutouts'])) for r in recs]}
        body = json.dumps(req, sort_keys=True, separators=(',', ':'))
        return Prepared([r['id'] for r in batch], material, {'ids': [r['id'] for r in batch], 'bytes': len(body),
                                                              'request': redact_for_display(req)}, len(body), payload=req)

    def extract(self, obj, meta, batch):
        ids = [x['id'] for x in batch]
        model = _str_or_none(obj.get('model')) if isinstance(obj, dict) else None
        rid = _str_or_none(obj.get('request_id')) if isinstance(obj, dict) else None
        items = obj.get('results') if isinstance(obj, dict) else None
        if not isinstance(items, list):
            return {i: {'error': 'response has no "results" list'} for i in ids}, model, rid
        by_id = {str(it.get('id')): it for it in items if isinstance(it, dict)}
        fields = self.fields or mapping.default_fields(self.task, items)
        res = {}
        for i in ids:
            it = by_id.get(i)
            if it is None:
                res[i] = {'error': 'object id %s not in results' % i}
            elif it.get('error') not in (None, ''):
                res[i] = {'error': 'service error: %s' % str(it['error'])[:300]}
            else:
                try:
                    res[i] = {'values': mapping.apply_fields(it, fields, obj)}
                except mapping.MapError as e:
                    res[i] = {'error': 'mapping: %s' % e}
        return res, model, rid


class LocalCommand(_RecordAdapter):
    """Runs a user-provided executable (argv list, no shell). stdin: one JSON document
    {"contract","task","service","params","records":[...]}; stdout: {"results":[{"id":..,<columns>..}],"model":..,"request_id":..}.
    Only a minimal environment is passed (PATH, HOME, LANG, TMPDIR, VIRTUAL_ENV, PYTHONPATH + profile.env_passthrough)."""

    def _env(self):
        keep = ['PATH', 'HOME', 'LANG', 'LC_ALL', 'TMPDIR', 'VIRTUAL_ENV', 'PYTHONPATH', 'SYSTEMROOT'] + \
            list(self.p.get('env_passthrough') or [])
        return {k: os.environ[k] for k in keep if k in os.environ}

    def send(self, prepared, timeout):
        cmd = [self._expand(a) for a in self.p['command']]
        try:
            pr = subprocess.run(cmd, input=json.dumps(prepared.payload).encode(), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, timeout=timeout, env=self._env(), cwd=self.p.get('cwd') or None)
        except subprocess.TimeoutExpired:
            raise RequestFailed('command timed out after %ss' % timeout, retryable=True)
        except FileNotFoundError:
            raise RequestFailed('command not found: %s' % cmd[0], retryable=False)
        except OSError as e:
            raise RequestFailed('cannot run command: %s' % e, retryable=False)
        if pr.returncode != 0:
            tail = pr.stderr.decode('utf-8', 'replace').strip().splitlines()[-1:] or ['']
            raise RequestFailed('command exited with status %d: %s' % (pr.returncode, tail[0][:200]), retryable=False)
        try:
            return json.loads(pr.stdout.decode('utf-8')), {}
        except ValueError:
            raise RequestFailed('command stdout is not JSON (%d bytes)' % len(pr.stdout), retryable=False)

    def _expand(self, a):
        return a.replace('{python}', sys.executable)

    def check(self):
        import shutil
        c = self._expand(self.p['command'][0])
        w = c if os.path.isabs(c) and os.access(c, os.X_OK) else shutil.which(c)
        return (bool(w), ('executable: %s' % w) if w else 'executable not found: %s' % c)


class PythonCallable(_RecordAdapter):
    """module:function(records, params, context) -> list of result dicts, or {"results": [...], "model":..,"request_id":..}.
    Optional profile keys: "sys_path": [dirs to add]."""

    def _resolve(self):
        mod, fn = self.p['callable'].split(':')
        for d in self.p.get('sys_path') or []:
            d = os.path.expanduser(d)
            if d not in sys.path:
                sys.path.insert(0, d)
        import importlib
        try:
            m = importlib.import_module(mod)
            f = getattr(m, fn)
        except Exception as e:
            raise RequestFailed('cannot import %s: %s: %s' % (self.p['callable'], type(e).__name__, e), retryable=False)
        if not callable(f):
            raise RequestFailed('%s is not callable' % self.p['callable'], retryable=False)
        return f

    def send(self, prepared, timeout):
        f = self._resolve()
        req = prepared.payload
        try:
            out = f(req['records'], req['params'], {'task': self.task, 'service': self.name, 'contract': CONTRACT_VERSION})
        except Exception as e:
            raise RequestFailed('callable raised %s: %s' % (type(e).__name__, e), retryable=False)
        if isinstance(out, list):
            out = {'results': out}
        try:
            json.dumps(out)
        except (TypeError, ValueError) as e:
            raise RequestFailed('callable result is not JSON-serialisable: %s' % e, retryable=False)
        return out, {}

    def check(self):
        try:
            self._resolve()
            return True, 'imported %s' % self.p['callable']
        except RequestFailed as e:
            return False, str(e)


class Mock(_RecordAdapter):
    """FAKE deterministic values. Every output is tagged SERVICE=mock / MODEL=mock-model (FAKE) so it can never be
    mistaken for a real result.  Values depend only on (object id, column name), so replays are identical."""

    def cutout_needs(self):
        # exercise the cutout path whenever the task contract wants cutouts and images were given
        need = contracts.TASKS.get(self.task, {}).get('needs', [])
        if self.ctx.get('bands') and any(n.startswith('cutouts') for n in need):
            fmt = (self.p.get('cutouts') or {}).get('format', 'png')
            return [(fmt, b, 'path') for b in self.ctx['bands']]
        return []

    def prepare(self, batch, cut):
        p = super().prepare(batch, cut)
        p.key_material['mock'] = True
        return p

    def send(self, prepared, timeout):
        task = prepared.payload['task']
        if task in (None, 'any'):
            raise RequestFailed('mock needs a task', retryable=False)
        out = []
        for rec in prepared.payload['records']:
            if self.ctx.get('mock_fail_ids') and rec['id'] in self.ctx['mock_fail_ids']:
                out.append({'id': rec['id'], 'error': 'mock: simulated failure for this object'})
                continue
            it = {'id': rec['id']}
            for name in contracts.output_names(task):
                it[name] = self._fake(rec['id'], name, contracts.output_spec(task, name))
            if task == 'photoz':        # keep the fake values self-consistent: P16 <= P50 <= P84
                z, e = it['PHOTOZ'], it['PHOTOZ_ERR']
                it.update(PHOTOZ_P16=round(max(0.0, z - e), 6), PHOTOZ_P50=z, PHOTOZ_P84=round(z + e, 6))
            if task == 'generic':
                it['MOCK_VALUE'] = self._unit('%s|%s' % (rec['id'], 'MOCK_VALUE'))
            out.append(it)
        return {'results': out, 'model': 'mock-model (FAKE VALUES)', 'request_id': 'mock-%s' % hashlib.sha256(
            json.dumps(prepared.ids).encode()).hexdigest()[:10]}, {}

    @staticmethod
    def _unit(s):
        return int(hashlib.sha256(s.encode()).hexdigest()[:12], 16) / float(16 ** 12)

    def _fake(self, oid, col, spec):
        u = self._unit('%s|%s' % (oid, col))
        g = spec['mock']
        if g[0] == 'float':
            return round(g[1] + u * (g[2] - g[1]), 6)
        if g[0] == 'choice':
            c = g[1]
            v = c[int(u * len(c)) % len(c)]
            return v == 'True' if spec['type'] == 'bool' else v
        return '%s%s' % (g[1], oid)

    def extract(self, obj, meta, batch):
        res, model, rid = super().extract(obj, meta, batch)
        return res, model, rid

    def check(self):
        return True, 'mock service (FAKE values) - always available offline'


# ---------------------------------------------------------------- IVOA TAP (non-ML example)
class TapQuery(Base):
    """Synchronous TAP query (IVOA TAP 1.x, <base_url>/sync, REQUEST=doQuery, LANG=ADQL) per object.  request.adql is a
    template with {ra} {dec} {radius_deg} {id} {param_*}.  The response (CSV) is exposed to the mapping as
        {"n_rows": N, "columns": [...], "rows": [ {col: value,...}, ... ]}   e.g. path "rows[0].main_id".
    {id} is inserted with single quotes doubled (ADQL string escaping)."""

    def __init__(self, profile, ctx):
        super().__init__(profile, ctx)
        self.auth = Auth(profile.get('auth'))
        self.req = profile['request']

    def cutout_needs(self):
        return []

    def _url(self):
        return self.p['base_url'].rstrip('/') + '/sync'

    def prepare(self, batch, cut):
        rec = batch[0]
        if rec.get('ra') is None or rec.get('dec') is None:
            raise TemplateError('object has no RA/Dec (catalog columns ALPHA_J2000/DELTA_J2000)')
        base = Lookup(self, rec, cut)

        def lk(n):
            v = base(n)
            return v.replace("'", "''") if isinstance(v, str) else v
        adql = templating.render(self.req['adql'], lk)
        if not isinstance(adql, str):
            raise TemplateError('adql did not render to text')
        fmt = self.req.get('format', 'csv')
        form = {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': fmt, 'QUERY': adql}
        form.update({k: str(v) for k, v in (self.req.get('extra_form') or {}).items()})
        h = {'Accept': 'text/csv, text/plain, */*', 'Content-Type': 'application/x-www-form-urlencoded'}
        h.update(self.p.get('headers') or {})
        body = urllib.parse.urlencode(form).encode()
        material = {'url': self._url(), 'form': form, 'headers': h}
        h2, q2 = self.auth.describe()
        preview = {'method': 'POST', 'url': self._url(), 'headers': dict(h, **h2), 'form': form, 'ids': [rec['id']], 'bytes': len(body)}
        return Prepared([rec['id']], material, preview, len(body), payload=(h, body))

    def send(self, prepared, timeout):
        h, body = prepared.payload
        headers, q = self.auth.apply(h, {})
        url = self._url() + (('?' + urllib.parse.urlencode(q)) if q else '')
        try:
            st, rh, data = httpio.request('POST', url, headers, body, timeout)
        except RequestFailed as e:
            m = re.search(r'<INFO[^>]*name="QUERY_STATUS"[^>]*value="ERROR"[^>]*>(.*?)</INFO>', e.body or '', re.S)
            if m:      # TAP reports ADQL errors as HTTP 400 + VOTable
                raise RequestFailed('TAP error: %s' % m.group(1).strip().replace('\n', ' ')[:300], retryable=False)
            raise RequestFailed(self.auth.scrub(str(e)), e.retryable, e.status, e.retry_after)
        text = data.decode('utf-8', 'replace')
        if text.lstrip().startswith('<'):
            m = re.search(r'<INFO[^>]*name="QUERY_STATUS"[^>]*value="ERROR"[^>]*>(.*?)</INFO>', text, re.S)
            raise RequestFailed('TAP error: %s' % ((m.group(1).strip() if m else text[:200]).replace('\n', ' ')[:300]), retryable=False)
        delim = '\t' if self.req.get('format') == 'tsv' else ','
        rd = csv.DictReader(io.StringIO(text), delimiter=delim)
        rows = [dict(r) for r in rd]
        return {'n_rows': len(rows), 'columns': list(rd.fieldnames or []), 'rows': rows}, {'http_status': st}

    def extract(self, obj, meta, batch):
        i = batch[0]['id']
        model = self.p.get('model_label') or 'TAP ' + urllib.parse.urlsplit(self.p['base_url']).netloc
        if obj.get('n_rows') == 0:      # a cone search without a match is a valid answer: empty cells, no error
            return {i: {'values': {c: None for c in self.fields}}}, model, None
        try:
            return {i: {'values': mapping.apply_fields(obj, self.fields)}}, model, None
        except mapping.MapError as e:
            return {i: {'error': 'mapping: %s' % e}}, None, None

    def check(self):
        url = self.p['base_url'].rstrip('/') + '/availability'
        try:
            st, rh, data = httpio.request('GET', url, {}, None, float(self.p.get('timeout_s', 30)))
        except RequestFailed as e:
            return False, str(e)
        txt = data.decode('utf-8', 'replace')
        m = re.search(r'<(?:\w+:)?available>\s*(true|false)\s*<', txt)
        return (m is None or m.group(1) == 'true'), 'GET %s -> HTTP %d%s' % (
            url, st, (' available=' + m.group(1)) if m else ' (no VOSI availability document recognised)')


ADAPTERS = {'https_json': GenericREST, 'https_multipart': GenericREST, 'local_command': LocalCommand,
            'python_callable': PythonCallable, 'mock': Mock, 'tap_query': TapQuery}


def make_adapter(profile, ctx):
    if profile.get('transport') == 'agent_cli':       # lazy: agent_cli imports this module
        from .agent_cli import AgentCLI
        return AgentCLI(profile, ctx)
    cls = ADAPTERS.get(profile.get('transport'))
    if cls is None:
        raise ProfileError('unknown transport %r' % profile.get('transport'))
    return cls(profile, ctx)
