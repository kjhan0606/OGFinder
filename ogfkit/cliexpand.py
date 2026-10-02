"""Headless expansion of plugin.json "cli" templates (the Python twin of ::ogf::step::build_argv in ds9/library/ogf_core.tcl).

Tokens: {python} {plugin_dir} {work} {image} {catalog} {root} {mask} {sextract} {image_tail} {base} {psf} {PARAM} {OTHER:PARAM} {cat:KEY} (empty) {script:NAME};
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
    """The step object; a step whose GUI path is a Tcl proc / dialog and that has a "headless" block is returned with that block merged in
    (cli, output, record, ... of the block replace the step's; proc and variants are dropped)."""
    for s in manifest.get('steps', []):
        if s['id'] == sid:
            if 'headless' in s:
                s = dict({k: v for k, v in s.items() if k not in ('proc', 'variants', 'headless')}, **s['headless'])
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
        if ':' in key and key in params:                     # recipe override of another plugin's parameter, e.g. "extract:detect-thresh"
            return param_text(params[key])
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
    ctx.setdefault('psf', '')
    ctx.setdefault('sextract', os.path.join(ctx['root'], 'bin', 'ds9_sextract'))
    ctx.setdefault('image_tail', os.path.basename(ctx.get('image', '') or ''))
    ctx.setdefault('base', re.sub(r'\.(fits|fit|fts)$', '', re.sub(r'\.gz$', '', ctx['image_tail'])))
    argv = []
    for e in step['cli']:
        if isinstance(e, dict) and 'argv' in e:
            if 'if_eq' in e:
                ie = e['if_eq']
                if param_text(params.get(ie[0])) not in [param_text(v) for v in ie[1:]]:
                    continue
            elif 'if_file' in e:
                f = expand_string(e['if_file'], ctx, manifests, pid, params)
                if not f or not os.path.isfile(f):
                    continue
            elif 'if_not_file' in e:
                f = expand_string(e['if_not_file'], ctx, manifests, pid, params)
                if f and os.path.isfile(f):
                    continue
            elif 'if_not' in e:
                if _truthy(params.get(e['if_not'])):
                    continue
            else:
                conds = e['if'].split() if isinstance(e['if'], str) else e['if']
                if not all(_truthy(params.get(c)) for c in conds):
                    continue
            argv += [expand_string(a, ctx, manifests, pid, params) for a in e['argv']]
        else:
            argv.append(expand_string(e, ctx, manifests, pid, params))
    return argv
