"""Service profiles: loading, defaults, validation, redacted display.

A profile file is JSON:   {"schema": "ogf-ai-bridge/1", "services": [ {profile}, ... ]}
(JSON has no comments: keys starting with "_" are ignored and may be used for notes.)
Secrets are NEVER stored: auth only names an environment variable.
"""
import copy
import hashlib
import json
import os
import re

from . import CONTRACT_VERSION, contracts, templating
from .errors import ProfileError

DEFAULT_PROFILE_FILE = os.path.join('~', '.ds9', 'ai_services.json')
TRANSPORTS = ('https_json', 'https_multipart', 'local_command', 'python_callable', 'tap_query', 'mock')
NETWORK_TRANSPORTS = ('https_json', 'https_multipart', 'tap_query')
AUTH_SCHEMES = ('none', 'bearer_env', 'header_env', 'query_env')
METHODS = ('GET', 'POST', 'PUT')
NAME_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$')
ENV_RE = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')

DEFAULTS = {
    'enabled': True, 'template': False, 'description': '', 'method': 'POST', 'headers': {},
    'auth': {'scheme': 'none'}, 'batch_size': 1, 'timeout_s': 30.0, 'retries': 2, 'backoff_s': 1.0,
    'backoff_factor': 2.0, 'backoff_max_s': 60.0, 'rate_limit_per_s': 0.0, 'max_payload_bytes': 10_000_000,
    'params': {}, 'cutouts': {}, 'request': {}, 'response': {},
}


def profile_path(path=None):
    p = path or os.environ.get('OGF_AI_SERVICES_FILE') or DEFAULT_PROFILE_FILE
    return os.path.expanduser(p)


def builtin_mock():
    return {'name': 'mock', 'task': 'any', 'transport': 'mock', 'enabled': True, 'builtin': True,
            'description': 'BUILT-IN MOCK: deterministic FAKE values for offline tests. Never real results.'}


def canonical(profile):
    """Canonical JSON text of a profile (keys sorted, no "_" notes); used for hashing."""
    def strip(o):
        if isinstance(o, dict):
            return {k: strip(v) for k, v in o.items() if not str(k).startswith('_')}
        if isinstance(o, list):
            return [strip(v) for v in o]
        return o
    return json.dumps(strip(profile), sort_keys=True, separators=(',', ':'), ensure_ascii=True)


def profile_hash(profile):
    return hashlib.sha256(canonical(profile).encode()).hexdigest()


def load_file(path=None):
    """-> (profiles list, path, error or None).  Missing file is not an error (empty list)."""
    p = profile_path(path)
    if not os.path.exists(p):
        return [], p, None
    try:
        with open(p, 'r', encoding='utf-8') as f:
            doc = json.load(f)
    except (OSError, ValueError) as e:
        return [], p, 'cannot read %s: %s' % (p, e)
    if isinstance(doc, dict) and 'services' in doc:
        lst = doc['services']
    elif isinstance(doc, list):
        lst = doc
    else:
        return [], p, '%s: expected {"services": [...]}' % p
    if isinstance(lst, dict):                       # tolerate {"services": {"name": {...}}}
        lst = [dict(v, name=k) for k, v in lst.items()]
    if not isinstance(lst, list) or not all(isinstance(x, dict) for x in lst):
        return [], p, '%s: "services" must be a list of objects' % p
    return lst, p, None


def all_services(path=None):
    """Profiles from the file plus the built-in mock (unless the file defines its own 'mock')."""
    lst, p, err = load_file(path)
    names = [x.get('name') for x in lst]
    if 'mock' not in names:
        lst = lst + [builtin_mock()]
    return lst, p, err


def find(name, path=None):
    lst, p, err = all_services(path)
    if err:
        raise ProfileError(err)
    for x in lst:
        if x.get('name') == name:
            return x
    raise ProfileError('service %r not found in %s (known: %s)' % (name, p, ', '.join(str(x.get('name')) for x in lst)))


def with_defaults(profile):
    p = copy.deepcopy(profile)
    for k, v in DEFAULTS.items():
        p.setdefault(k, copy.deepcopy(v))
    for k in ('request', 'response', 'cutouts', 'params', 'headers'):
        if p.get(k) is None:
            p[k] = {}
    return p


def env_status(profile):
    a = (profile.get('auth') or {})
    if a.get('scheme', 'none') == 'none':
        return None, '-'
    e = a.get('env')
    return e, ('set' if e and os.environ.get(e) else 'unset')


def validate(profile, task=None):
    """-> (errors, warnings) lists of strings.  Does not touch the network."""
    errs, warns = [], []
    pr = with_defaults(profile)
    nm = pr.get('name')
    if not isinstance(nm, str) or not NAME_RE.match(nm or ''):
        errs.append('name: required, letters/digits/_.- only (got %r)' % (nm,))
    t = pr.get('task')
    if t != 'any' and t not in contracts.TASKS:
        errs.append('task: %r is not one of %s' % (t, ', '.join(contracts.TASK_NAMES)))
    if t == 'any' and pr.get('transport') != 'mock':
        errs.append('task "any" is only allowed for the mock transport')
    if task and t not in ('any', task):
        errs.append('profile is for task %r, requested %r' % (t, task))
    tr = pr.get('transport')
    if tr not in TRANSPORTS:
        errs.append('transport: %r is not one of %s' % (tr, ', '.join(TRANSPORTS)))
    for k, lo in (('batch_size', 1), ('retries', 0)):
        if not isinstance(pr.get(k), int) or isinstance(pr.get(k), bool) or pr[k] < lo:
            errs.append('%s: must be an integer >= %d' % (k, lo))
    for k in ('timeout_s', 'backoff_s', 'backoff_factor', 'backoff_max_s', 'rate_limit_per_s', 'max_payload_bytes'):
        if not isinstance(pr.get(k), (int, float)) or isinstance(pr.get(k), bool) or pr[k] < 0:
            errs.append('%s: must be a number >= 0' % k)
    if isinstance(pr.get('timeout_s'), (int, float)) and pr['timeout_s'] == 0:
        errs.append('timeout_s: must be > 0')
    a = pr.get('auth') or {}
    sch = a.get('scheme', 'none')
    if sch not in AUTH_SCHEMES:
        errs.append('auth.scheme: %r is not one of %s' % (sch, ', '.join(AUTH_SCHEMES)))
    if sch != 'none':
        if not ENV_RE.match(str(a.get('env', ''))):
            errs.append('auth.env: must be the NAME of an environment variable (the secret itself is never stored)')
        if sch == 'header_env' and not a.get('header'):
            errs.append('auth.header: required for header_env')
        if sch == 'query_env' and not a.get('param'):
            errs.append('auth.param: required for query_env')
        for k in a:
            if k not in ('scheme', 'env', 'header', 'param', 'prefix') and not k.startswith('_'):
                errs.append('auth.%s: unknown key (only the env-var NAME is stored; no secrets in profiles)' % k)
    for k in list(pr['headers']):
        if re.search(r'authorization|api[-_]?key|token|secret', k, re.I):
            warns.append('headers.%s looks like a credential header: put secrets in the environment and use auth.scheme' % k)
    secret_like = re.compile(r'(Bearer\s+[A-Za-z0-9._~+/=-]{12,}|sk-[A-Za-z0-9]{16,}|[A-Za-z0-9]{32,})')
    for k, v in pr['headers'].items():
        if isinstance(v, str) and secret_like.search(v):
            errs.append('headers.%s: value looks like a literal secret; use auth.scheme + environment variable' % k)
    pp = pr.get('provenance_prefix')
    if pp is not None and not re.match(r'^[A-Z][A-Z0-9]{0,11}$', str(pp)):
        errs.append('provenance_prefix: 1-12 characters A-Z0-9, starting with a letter')
    if tr in NETWORK_TRANSPORTS:
        url = pr.get('base_url')
        if not isinstance(url, str) or not re.match(r'^https?://', url or ''):
            errs.append('base_url: required (https://...) for transport %s' % tr)
        elif url.startswith('http://') and not re.match(r'^http://(localhost|127\.0\.0\.1|\[::1\])(:\d+)?(/|$)', url) \
                and not pr.get('allow_insecure_http'):
            errs.append('base_url: plain http:// is refused except for localhost (set "allow_insecure_http": true to override, NOT recommended)')
        if isinstance(pr.get('base_url'), str) and 'example.invalid' in pr['base_url']:
            warns.append('base_url uses the placeholder host example.invalid: this profile is a TEMPLATE and cannot work until edited')
        if pr.get('method') not in METHODS:
            errs.append('method: %r is not one of %s' % (pr.get('method'), ', '.join(METHODS)))
        if pr.get('endpoint') is not None and not isinstance(pr.get('endpoint'), str):
            errs.append('endpoint: must be a string path')
    if tr == 'https_multipart':
        mp = pr['request'].get('multipart')
        if not isinstance(mp, dict) or not isinstance(mp.get('files'), dict) or not mp['files']:
            errs.append('request.multipart.files: required for https_multipart ({"form_field": "{cutout_png_path}"})')
        if pr['batch_size'] != 1:
            errs.append('batch_size must be 1 for https_multipart')
    if tr == 'local_command':
        c = pr.get('command')
        if not isinstance(c, list) or not c or not all(isinstance(x, str) for x in c):
            errs.append('command: required, a list of strings (argv, no shell)')
        for e in pr.get('env_passthrough') or []:
            if not ENV_RE.match(str(e)):
                errs.append('env_passthrough: %r is not an environment variable name' % (e,))
    if tr == 'python_callable':
        if not re.match(r'^[A-Za-z_][\w.]*:[A-Za-z_]\w*$', str(pr.get('callable', ''))):
            errs.append('callable: required, "module:function"')
    if tr == 'tap_query':
        if not isinstance(pr['request'].get('adql'), str):
            errs.append('request.adql: required for tap_query (ADQL text with {placeholders})')
        if pr['batch_size'] != 1:
            errs.append('batch_size must be 1 for tap_query')
        if pr['request'].get('format', 'csv') not in ('csv', 'tsv'):
            errs.append('request.format: csv or tsv')
    # templates
    known = set(KNOWN_PLACEHOLDERS) | {'items', 'n_items', 'batch_ids'}
    cut_pat = re.compile(r'^cutout(s)?(_(png|fits|npy))?(_(base64|path|url))?(_[A-Za-z0-9.-]+)?$')
    tpl_objs = [pr['request'].get(k) for k in ('template', 'query', 'item_template', 'adql', 'multipart')]
    for tp in tpl_objs:
        if tp is None:
            continue
        for ph in templating.placeholders(tp):
            if ph in known or ph.startswith('param_') or cut_pat.match(ph):
                continue
            errs.append('request template uses unknown placeholder {%s} (known: %s, cutout_<fmt>_<base64|path|url>[_BAND], param_<name>)'
                        % (ph, ', '.join(sorted(KNOWN_PLACEHOLDERS))))
    if pr['batch_size'] > 1 and tr in ('https_json',):
        if 'items' not in templating.placeholders(pr['request'].get('template', {})):
            errs.append('batch_size > 1 needs "{items}" in request.template (and optionally request.item_template)')
    # response mapping
    rs = pr['response']
    fields = rs.get('fields')
    default_fields = tr == 'mock' or (tr in ('local_command', 'python_callable') and not fields)
    if not default_fields:
        if not isinstance(fields, dict) or not fields:
            errs.append('response.fields: required, {"COLUMN": {"path": "a.b", "type": "float"}}')
            fields = {}
    fields = fields or {}
    for col, spec in fields.items():
        if not re.match(r'^[A-Za-z][A-Za-z0-9_]*$', col):
            errs.append('response.fields: column name %r must be letters/digits/_' % col)
            continue
        if col == 'NUMBER' or re.match(r'^AI_[A-Z0-9]+_(SERVICE|MODEL|REQID|TIME|ERROR)$', col):
            errs.append('response.fields: column name %r is reserved' % col)
        if not isinstance(spec, dict) or not isinstance(spec.get('path'), str) or not spec['path']:
            errs.append('response.fields.%s.path: required dotted path string' % col)
            continue
        ty = spec.get('type', 'str')
        if ty not in ('float', 'int', 'str', 'bool', 'json'):
            errs.append('response.fields.%s.type: %r not in float,int,str,bool,json' % (col, ty))
        if spec.get('unit') or spec.get('to_unit'):
            from . import units
            try:
                units.factor(spec.get('unit') or spec.get('to_unit'), spec.get('to_unit') or spec.get('unit'))
            except ValueError as e:
                errs.append('response.fields.%s: %s' % (col, e))
        if t in contracts.TASKS and t != 'generic':
            cs = contracts.output_spec(t, col)
            if cs is None:
                warns.append('response.fields.%s is not a standard column of task %r (kept as an extra column)' % (col, t))
            elif cs['type'] != ty:
                warns.append('response.fields.%s: type %s differs from the task contract (%s)' % (col, ty, cs['type']))
    if t in contracts.TASKS and t != 'generic' and fields and not default_fields:
        if not any(contracts.output_spec(t, c) for c in fields):
            errs.append('response.fields maps none of the standard columns of task %r (%s)' % (t, ', '.join(contracts.output_names(t))))
    if tr == 'tap_query' and t != 'generic':
        warns.append('tap_query is a non-ML example adapter; task "generic" is recommended')
    cu = pr.get('cutouts') or {}
    if cu.get('size_pix') and cu.get('size_arcsec'):
        errs.append('cutouts: give size_pix or size_arcsec, not both')
    if cu.get('normalize', 'asinh') not in ('none', 'linear', 'zscale', 'asinh'):
        errs.append('cutouts.normalize: none|linear|zscale|asinh')
    if cu.get('format', 'png') not in ('png', 'fits', 'npy'):
        errs.append('cutouts.format: png|fits|npy')
    return errs, warns


KNOWN_PLACEHOLDERS = (
    'id', 'x', 'y', 'ra', 'dec', 'mags', 'mag_errs', 'mags_json', 'mag_errs_json', 'catalog_row', 'catalog_row_json',
    'pixel_scale_arcsec', 'psf_fwhm_arcsec', 'task', 'service', 'bands', 'bands_json', 'radius_deg',
    'contract', 'cutouts_json',
)


def redacted(profile):
    """Copy that is safe to print (profiles never hold secrets; the env value is never read here)."""
    p = copy.deepcopy(profile)
    e, st = env_status(p)
    if e:
        p['auth']['_env_status'] = st
    return p
