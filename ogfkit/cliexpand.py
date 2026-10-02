"""Headless expansion of plugin.json "cli" templates (the Python twin of ::ogf::step::build_argv in ds9/library/ogf_core.tcl).

Tokens: {python} {plugin_dir} {work} {image} {catalog} {root} {mask} {PARAM} {OTHER:PARAM} {cat:KEY} (empty) {script:NAME};
elements may be {"if": "param", "argv": [...]} (included when the parameter is true / non-empty / non-zero) or {"if_file": "{token}", "argv": [...]}.
"""
import json
import os
import re

TOKEN = re.compile(r'\{([A-Za-z0-9_.:,-]+)\}')


def load_manifests(root):
    out = {}
    d = os.path.join(root, 'plugins')
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name, 'plugin.json')
        if os.path.isfile(p):
            m = json.load(open(p))
            m['_dir'] = os.path.join(d, name)
            out[m['id']] = m
    return out


def param_text(v):
    if isinstance(v, bool):
        return '1' if v else '0'
    return str(v)


def defaults(manifest):
    return {p['name']: p.get('default', '') for p in manifest.get('params', [])}


def find_step(manifest, sid):
    for s in manifest.get('steps', []):
        if s['id'] == sid:
            return s
    raise KeyError('plugin %s has no step %s' % (manifest['id'], sid))


def _truthy(v):
    if v is None:
        return False
    s = param_text(v)
    return not (s == '' or s == '0' or s.lower() in ('false', 'no', 'off'))


def expand_string(s, ctx, manifests, pid, params):
    def rep(m):
        key = m.group(1)
        if key in ctx:
            return str(ctx[key])
        mm = re.match(r'^script:(.+)$', key)
        if mm:
            for d in (os.path.join(ctx['root'], 'bin'), os.path.join(ctx['root'], 'ds9', 'library')):
                p = os.path.join(d, mm.group(1))
                if os.path.exists(p):
                    return p
            return os.path.join(ctx['root'], 'ds9', 'library', mm.group(1))
        if key.startswith('cat:'):
            return ''
        mm = re.match(r'^([A-Za-z0-9_]+):([A-Za-z0-9_.-]+)$', key)
        if mm and mm.group(1) in manifests and mm.group(2) in defaults(manifests[mm.group(1)]):
            return param_text(defaults(manifests[mm.group(1)])[mm.group(2)])
        if key in params:
            return param_text(params[key])
        raise KeyError('unknown token {%s} in argv template of plugin %s' % (key, pid))
    return TOKEN.sub(rep, s)


def build_argv(manifests, pid, sid, ctx, overrides=None):
    """-> argv list for step `sid` of plugin `pid`.  ctx needs python, plugin_dir (set automatically), work, image, catalog, root."""
    m = manifests[pid]
    step = find_step(m, sid)
    if 'cli' not in step:
        raise KeyError('step %s.%s has no cli template' % (pid, sid))
    params = defaults(m)
    params.update(overrides or {})
    ctx = dict(ctx)
    ctx.setdefault('plugin_dir', m['_dir'])
    ctx.setdefault('mask', '')
    argv = []
    for e in step['cli']:
        if isinstance(e, dict) and 'argv' in e:
            if 'if_file' in e:
                f = expand_string(e['if_file'], ctx, manifests, pid, params)
                if not f or not os.path.isfile(f):
                    continue
            else:
                if not _truthy(params.get(e['if'])):
                    continue
            argv += [expand_string(a, ctx, manifests, pid, params) for a in e['argv']]
        else:
            argv.append(expand_string(e, ctx, manifests, pid, params))
    return argv
