"""Reproducibility bundles: what was run, with which software, on which inputs, giving which outputs - and a verify that re-runs and compares.

bundle (zip):  manifest.json (schema 1: created, tool/git, host, packages, binaries, inputs[sha256], outputs[sha256, rows, columns], kind, verify command),
               params.json, steps.json, requirements.lock, catalog_final.tsv (or per-field catalogs), README.txt, verify.sh,
               recipe.json + fields.txt (kind "batch") or session.py (kind "session"), optional inputs/ (copies, --include-inputs).
verify:        environment diff (warning), input checksums (failure), re-run (batch recipe through ogfkit.batchrun, or the session script in replay mode),
               table comparison (identical / numerically within rtol / different).
"""
import datetime
import hashlib
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import zipfile

SCHEMA = 1
KEY_PACKAGES = ['numpy', 'scipy', 'astropy', 'sep', 'matplotlib', 'scikit-learn', 'photutils', 'eazy', 'bagpipes', 'torch', 'Pillow', 'numba']


def sha256_file(path, bufsize=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def installed_packages():
    """{name: version} of every distribution visible to this interpreter (no network, no pip needed)."""
    from importlib import metadata
    out = {}
    for d in metadata.distributions():
        n = d.metadata['Name']
        if n:
            out[n] = d.version
    return dict(sorted(out.items(), key=lambda kv: kv[0].lower()))


def requirements_lock(pk=None):
    pk = pk or installed_packages()
    return '\n'.join('%s==%s' % (k, v) for k, v in pk.items()) + '\n'


def git_info(root):
    def g(*a):
        try:
            return subprocess.run(['git', '-C', root] + list(a), capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:
            return ''
    head = g('rev-parse', 'HEAD')
    if not head:
        return dict(available=False)
    st = [l for l in g('status', '--porcelain').splitlines() if l.strip()]
    tracked = [l for l in st if not l.startswith('??')]
    return dict(available=True, head=head, branch=g('rev-parse', '--abbrev-ref', 'HEAD'), describe=g('describe', '--always', '--dirty'), dirty=bool(tracked),
                modified_files=[l[3:] for l in tracked][:50], diff_sha256=sha256_bytes(g('diff', 'HEAD').encode()) if tracked else '')


def binary_info(root):
    out = {}
    for n in ('ds9_sextract', 'ds9'):
        p = os.path.join(root, 'bin', n)
        if os.path.isfile(p):
            out[n] = dict(path=p, size=os.path.getsize(p), sha256=sha256_file(p))
    return out


def file_record(path, name=None):
    r = dict(name=name or os.path.basename(path), path=os.path.abspath(path), exists=os.path.isfile(path))
    if r['exists']:
        r.update(size=os.path.getsize(path), sha256=sha256_file(path))
    return r


def read_tsv(path):
    lines = [l.rstrip('\n') for l in open(path, encoding='utf-8') if l.strip() and not l.startswith('#')]
    if not lines:
        return [], []
    cols = lines[0].split('\t')
    return cols, [l.split('\t') + [''] * (len(cols) - len(l.split('\t'))) for l in lines[1:]]


def table_info(path):
    cols, rows = read_tsv(path)
    return dict(rows=len(rows), columns=len(cols), column_names=cols)


def compare_tables(a, b, rtol=0.0, atol=0.0):
    """-> dict(status identical|numeric|different, max_rel_diff, n_cells_different, columns_different, detail)."""
    if sha256_file(a) == sha256_file(b):
        return dict(status='identical', max_rel_diff=0.0, n_cells_different=0, columns_different=[], detail='byte-identical')
    ca, ra = read_tsv(a); cb, rb = read_tsv(b)
    if ca != cb:
        return dict(status='different', max_rel_diff=None, n_cells_different=None, columns_different=sorted(set(ca) ^ set(cb)), detail='different column lists (%d vs %d)' % (len(ca), len(cb)))
    if len(ra) != len(rb):
        return dict(status='different', max_rel_diff=None, n_cells_different=None, columns_different=[], detail='different row counts (%d vs %d)' % (len(ra), len(rb)))
    nd = 0; mx = 0.0; bad = set(); nonnum = 0
    for x, y in zip(ra, rb):
        for j, (u, v) in enumerate(zip(x, y)):
            if u == v:
                continue
            try:
                fu, fv = float(u), float(v)
                rel = abs(fu - fv) / max(abs(fu), abs(fv), 1e-300)
                if abs(fu - fv) <= atol + rtol * max(abs(fu), abs(fv)):
                    mx = max(mx, rel); continue
                mx = max(mx, rel)
            except ValueError:
                nonnum += 1
            nd += 1; bad.add(ca[j])
    if nd == 0:
        return dict(status='numeric', max_rel_diff=mx, n_cells_different=0, columns_different=[], detail='all differences within rtol=%g atol=%g (max relative %.3g)' % (rtol, atol, mx))
    return dict(status='different', max_rel_diff=mx, n_cells_different=nd, columns_different=sorted(bad), detail='%d cells differ in %d columns (max relative %.3g)' % (nd, len(bad), mx))


README = """OGFinder reproducibility bundle (schema %d)
created: %s
kind: %s

Contents: manifest.json (software versions, git state, input and output checksums), params.json, steps.json, requirements.lock, %s.
Verify (re-runs the processing and compares with the bundled outputs):
    python3 plugins/repro/repro.py verify %s [--inputs-dir DIR] [--rtol 0]
Bit-identical results require the same OGFinder commit, binaries and package versions as listed in manifest.json; differences in versions are reported by verify.
"""


def make_bundle(out_zip, kind, root, inputs, outputs, params=None, steps=None, recipe=None, fields_txt=None, session_script=None, extra_files=None,
                include_inputs=False, note='', verify_args=''):
    """inputs / outputs: list of (name, path).  Returns the manifest dict."""
    pk = installed_packages()
    now = datetime.datetime.now().astimezone()
    man = dict(schema=SCHEMA, kind=kind, created=now.isoformat(timespec='seconds'), created_utc=now.astimezone(datetime.timezone.utc).isoformat(timespec='seconds'), note=note,
               tool=git_info(root), binaries=binary_info(root), host=dict(platform=platform.platform(), machine=platform.machine(), python=sys.version.split()[0], executable=sys.executable),
               packages={k: pk[k] for k in pk if k.lower() in [x.lower() for x in KEY_PACKAGES]}, n_packages_locked=len(pk))
    man['inputs'] = [file_record(p, n) for n, p in inputs]
    man['outputs'] = []
    files = {}
    for n, p in outputs:
        if not os.path.isfile(p):
            continue
        rec = file_record(p, n)
        rec.update({k: v for k, v in table_info(p).items() if k != 'column_names'}) if p.endswith('.tsv') else None
        arc = 'outputs/%s' % n
        rec['archive'] = arc
        files[arc] = p
        man['outputs'].append(rec)
    if include_inputs:
        for r in man['inputs']:
            if r['exists']:
                r['archive'] = 'inputs/%s' % r['name']
                files[r['archive']] = r['path']
    man['verify'] = dict(command='python3 plugins/repro/repro.py verify %s' % os.path.basename(out_zip) + (' ' + verify_args if verify_args else ''))
    man['has'] = dict(recipe=recipe is not None, session=session_script is not None, inputs=include_inputs)
    with zipfile.ZipFile(out_zip, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('manifest.json', json.dumps(man, indent=1, sort_keys=True))
        z.writestr('params.json', json.dumps(params or {}, indent=1, sort_keys=True))
        z.writestr('steps.json', json.dumps(steps or [], indent=1))
        z.writestr('requirements.lock', requirements_lock(pk))
        z.writestr('README.txt', README % (SCHEMA, man['created'], kind, 'recipe.json + fields.txt' if recipe is not None else ('session.py' if session_script else 'outputs/'), os.path.basename(out_zip)))
        z.writestr('verify.sh', '#!/bin/sh\n%s "$@"\n' % man['verify']['command'])
        if recipe is not None:
            z.writestr('recipe.json', json.dumps(recipe, indent=1))
        if fields_txt is not None:
            z.writestr('fields.txt', fields_txt)
        if session_script:
            z.write(session_script, 'session.py')
        for arc, p in files.items():
            z.write(p, arc)
        for arc, p in (extra_files or {}).items():
            z.write(p, arc)
    return man


def read_manifest(zpath):
    with zipfile.ZipFile(zpath) as z:
        return json.loads(z.read('manifest.json'))


def verify_bundle(zpath, root, python=None, inputs_dir=None, workdir=None, rtol=0.0, atol=0.0, strict_env=False, rerun=True, jobs=1, session_args=None):
    """-> dict(ok, checks=[{name, status PASS|WARN|FAIL|SKIP, detail}], ...).  Never raises for content problems."""
    checks = []

    def add(name, status, detail=''):
        checks.append(dict(name=name, status=status, detail=detail))
    try:
        zf = zipfile.ZipFile(zpath)
        man = json.loads(zf.read('manifest.json'))
    except Exception as e:
        return dict(ok=False, checks=[dict(name='bundle', status='FAIL', detail='cannot read bundle: %s' % e)])
    add('bundle', 'PASS' if man.get('schema') == SCHEMA else 'FAIL', 'schema %s, kind %s, created %s' % (man.get('schema'), man.get('kind'), man.get('created')))
    wd = workdir or tempfile.mkdtemp(prefix='ogf_verify_')
    os.makedirs(wd, exist_ok=True)
    zf.extractall(wd)
    # environment
    pk = installed_packages()
    diffs = []
    for k, v in (man.get('packages') or {}).items():
        have = next((pv for pn, pv in pk.items() if pn.lower() == k.lower()), None)
        if have != v:
            diffs.append('%s %s -> %s' % (k, v, have))
    pyv = (man.get('host') or {}).get('python', '')
    if pyv and pyv.split('.')[:2] != sys.version.split()[0].split('.')[:2]:
        diffs.append('python %s -> %s' % (pyv, sys.version.split()[0]))
    add('environment', ('FAIL' if strict_env else 'WARN') if diffs else 'PASS', '; '.join(diffs[:8]) if diffs else 'package versions match (%d key packages)' % len(man.get('packages') or {}))
    gi = git_info(root)
    mt = man.get('tool') or {}
    if mt.get('available') and gi.get('available'):
        same = gi['head'] == mt['head'] and gi.get('diff_sha256', '') == mt.get('diff_sha256', '')
        add('code version', 'PASS' if same else 'WARN', ('commit %s' % gi['head'][:9]) if same else 'bundle %s%s, now %s%s' % (mt['head'][:9], '+local changes' if mt.get('dirty') else '', gi['head'][:9], '+local changes' if gi.get('dirty') else ''))
    else:
        add('code version', 'SKIP', 'git information not available')
    bn = man.get('binaries') or {}
    for n, r in bn.items():
        p = os.path.join(root, 'bin', n)
        if os.path.isfile(p):
            add('binary ' + n, 'PASS' if sha256_file(p) == r['sha256'] else 'WARN', 'sha256 ' + ('matches' if sha256_file(p) == r['sha256'] else 'differs from the bundle'))
    # inputs
    paths = {}
    bad_inputs = False
    for r in man.get('inputs', []):
        cand = []
        if inputs_dir:
            cand.append(os.path.join(inputs_dir, r['name']))
        cand.append(r['path'])
        if r.get('archive'):
            cand.append(os.path.join(wd, r['archive']))
        found = next((c for c in cand if os.path.isfile(c)), None)
        if found is None:
            add('input ' + r['name'], 'FAIL', 'not found (looked at %s)' % ', '.join(cand[:3])); bad_inputs = True
            continue
        ok = r.get('sha256') == sha256_file(found)
        add('input ' + r['name'], 'PASS' if ok else 'FAIL', ('sha256 matches %s' % found) if ok else 'checksum differs from the bundle (%s)' % found)
        bad_inputs |= not ok
        paths[r['name']] = found
    # outputs of the bundle itself must match the manifest
    for r in man.get('outputs', []):
        p = os.path.join(wd, r['archive'])
        ok = os.path.isfile(p) and sha256_file(p) == r['sha256']
        add('bundled output ' + r['name'], 'PASS' if ok else 'FAIL', 'sha256 matches manifest' if ok else 'bundled file differs from the manifest checksum (bundle altered?)')
    if not rerun or bad_inputs:
        if bad_inputs:
            add('re-run', 'SKIP', 'inputs missing or changed')
        return dict(ok=not any(c['status'] == 'FAIL' for c in checks), checks=checks, workdir=wd)
    # re-run
    produced = {}
    kind = man.get('kind')
    try:
        if kind == 'batch':
            from . import batchrun
            recipe = json.load(open(os.path.join(wd, 'recipe.json')))
            fields = batchrun.read_fields(os.path.join(wd, 'fields.txt'))
            fields = [(n, paths.get(os.path.basename(im), paths.get(n, im))) for n, im in fields]
            if inputs_dir:
                fields = [(n, os.path.join(inputs_dir, os.path.basename(im)) if os.path.isfile(os.path.join(inputs_dir, os.path.basename(im))) else im) for n, im in fields]
            out = os.path.join(wd, 'rerun')
            res = batchrun.run_batch(fields, recipe, out, root, python or sys.executable, jobs=jobs)
            bad = [r for r in res if r['status'] != 'ok']
            add('re-run (batch)', 'FAIL' if bad else 'PASS', ('%d of %d fields failed: %s' % (len(bad), len(res), bad[0]['failed_step'] + ': ' + bad[0]['message'])) if bad else '%d fields, %d steps' % (len(res), sum(r['steps_ok'] for r in res)))
            for r in man['outputs']:
                n = r['name'].split('/')[0]
                produced[r['name']] = os.path.join(out, n, 'catalog.tsv')
        elif kind == 'session':
            out = os.path.join(wd, 'rerun')
            cmd = [python or sys.executable, os.path.join(wd, 'session.py'), '--mode', 'replay', '--outdir', out, '--ogfinder-root', root, '--python', python or sys.executable] + list(session_args or [])
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            open(os.path.join(wd, 'rerun.log'), 'w').write(r.stdout + r.stderr)
            add('re-run (session replay)', 'PASS' if r.returncode == 0 else 'FAIL', 'rc=%d, log %s' % (r.returncode, os.path.join(wd, 'rerun.log')))
            for o in man['outputs']:
                hits = []
                for dp, dn, fn in os.walk(out):
                    if os.path.basename(o['name']) in fn:
                        hits.append(os.path.join(dp, os.path.basename(o['name'])))
                if hits:
                    produced[o['name']] = sorted(hits)[0]
        else:
            add('re-run', 'SKIP', 'kind %r has no re-run (integrity checks only)' % kind)
    except Exception as e:
        add('re-run', 'FAIL', '%s: %s' % (type(e).__name__, e))
    for r in man.get('outputs', []):
        if r['name'] not in produced:
            if kind in ('batch', 'session'):
                add('compare ' + r['name'], 'FAIL', 'the re-run did not produce this file')
            continue
        if not r['name'].endswith('.tsv'):
            ok = sha256_file(produced[r['name']]) == r['sha256']
            add('compare ' + r['name'], 'PASS' if ok else 'FAIL', 'sha256 ' + ('identical' if ok else 'differs'))
            continue
        c = compare_tables(os.path.join(wd, r['archive']), produced[r['name']], rtol, atol)
        add('compare ' + r['name'], {'identical': 'PASS', 'numeric': 'PASS', 'different': 'FAIL'}[c['status']], '%s: %s' % (c['status'], c['detail']))
    return dict(ok=not any(c['status'] == 'FAIL' for c in checks), checks=checks, workdir=wd)
