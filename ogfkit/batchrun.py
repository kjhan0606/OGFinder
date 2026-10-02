"""Batch runner: a recipe of plugin steps over many fields, with parallel workers, resume, per-field logs and a summary.

Recipe (JSON):
  {"detect": {"args": ["--detect-thresh", "1.5"]}            # optional: bin/ds9_sextract <image> + args -> catalog   (or)
   "catalog": "/data/{name}.tsv",                            # optional: catalog file per field ({name} = field name)
   "steps": [{"plugin": "cluster", "step": "members", "params": {"dens-sigma": 150}}, ...]}
Fields file: one field per line, "NAME path/to/image.fits" (extra tokens ignored; '#' comments).
Per field the runner keeps  <out>/<name>/{catalog.tsv, state.json, field.log, work/}; <out>/summary.tsv|json and <out>/status.json (live progress).
Steps are run exactly as the GUI does (same argv from the manifest templates); add_columns outputs are merged into the field catalog by NUMBER.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback

from . import cliexpand, tsvio


def sha_file(p, n=1 << 20):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        h.update(f.read(n))
    h.update(str(os.path.getsize(p)).encode())
    return h.hexdigest()


def read_fields(path):
    out = []
    for l in open(path):
        l = l.strip()
        if not l or l.startswith('#'):
            continue
        t = l.split()
        out.append((t[0], t[1] if len(t) > 1 else ''))
    names = [n for n, _ in out]
    if len(set(names)) != len(names):
        raise ValueError('duplicate field names in %s' % path)
    return out


def read_tsv_text(text):
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith('#')]
    if not lines:
        return [], []
    cols = lines[0].split('\t')
    rows = []
    for l in lines[1:]:
        v = l.split('\t'); v += [''] * (len(cols) - len(v))
        rows.append(v)
    return cols, rows


def merge_add_columns(cat_cols, cat_rows, out_text):
    """Merge a step's stdout (NUMBER + new columns) into the catalog; same-named columns are replaced.  -> (cols, rows, added)"""
    ocols, orows = read_tsv_text(out_text)
    if not ocols or ocols[0] != 'NUMBER':
        raise ValueError('step output does not start with a NUMBER column')
    ni = cat_cols.index('NUMBER')
    byn = {r[0]: r for r in orows}
    cols = list(cat_cols); rows = [list(r) for r in cat_rows]
    added = []
    for j, c in enumerate(ocols[1:], 1):
        if c in cols:
            k = cols.index(c)
        else:
            cols.append(c); k = len(cols) - 1
            for r in rows:
                r.append('')
            added.append(c)
        for r in rows:
            o = byn.get(r[ni])
            r[k] = o[j] if o is not None else ''
    return cols, rows, added


def write_tsv(path, cols, rows):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')
    os.replace(tmp, path)


def atomic_json(path, obj):
    tmp = path + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


class Field:
    def __init__(self, name, image, recipe, outdir, root, python, manifests, resume, retries, timeout, status_cb=None):
        self.name, self.image, self.recipe = name, image, recipe
        self.dir = os.path.join(outdir, name)
        self.root, self.python, self.manifests = root, python, manifests
        self.resume, self.retries, self.timeout = resume, retries, timeout
        self.cb = status_cb or (lambda *a, **k: None)
        os.makedirs(os.path.join(self.dir, 'work'), exist_ok=True)
        self.logf = open(os.path.join(self.dir, 'field.log'), 'a')

    def log(self, msg):
        self.logf.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg)); self.logf.flush()

    def run_cmd(self, argv, tag):
        """-> (rc, stdout, stderr, seconds); timeout -> rc 124."""
        t0 = time.time()
        try:
            r = subprocess.run(argv, capture_output=True, text=True, timeout=self.timeout or None, cwd=self.dir)
            rc, out, err = r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired as e:
            rc, out, err = 124, (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or ''), 'timeout after %ss' % self.timeout
        except OSError as e:
            rc, out, err = 127, '', str(e)
        dt = time.time() - t0
        self.log('%s: rc=%d %.1fs' % (tag, rc, dt))
        if err.strip():
            for l in err.strip().splitlines()[-12:]:
                self.log('  | ' + l)
        return rc, out, err, dt

    def plan(self):
        st = []
        if self.recipe.get('detect') is not None:
            st.append(dict(kind='detect', label='detect'))
        elif self.recipe.get('catalog'):
            st.append(dict(kind='catalog', label='catalog'))
        for s in self.recipe.get('steps', []):
            st.append(dict(kind='step', label='%s.%s' % (s['plugin'], s['step']), spec=s))
        return st

    def run(self):
        t_start = time.time()
        stp = os.path.join(self.dir, 'state.json')
        old = json.load(open(stp)) if (self.resume and os.path.isfile(stp)) else {}
        state = dict(name=self.name, image=self.image, steps={})
        plan = self.plan()
        cat_path = os.path.join(self.dir, 'catalog.tsv')
        cols, rows = [], []
        chain = hashlib.sha256((self.image + (sha_file(self.image) if self.image and os.path.isfile(self.image) else '')).encode()).hexdigest()
        res = dict(name=self.name, status='ok', steps_ok=0, steps_cached=0, steps_failed=0, failed_step='', message='', n_objects=0, n_columns=0, seconds=0.0)
        broken = False
        for i, s in enumerate(plan):
            lab = s['label']
            self.cb(self.name, 'running', lab, i, len(plan))
            if broken:
                state['steps'][lab] = dict(status='skipped')
                continue
            try:
                key, argv = self._key_and_argv(s, chain, cols, rows, cat_path)
            except Exception as e:
                self.log('%s: cannot build command: %s' % (lab, e))
                state['steps'][lab] = dict(status='failed', message=str(e)); res.update(status='failed', failed_step=lab, message=str(e)); res['steps_failed'] += 1; broken = True
                continue
            prev = old.get('steps', {}).get(lab)
            snap = os.path.join(self.dir, 'work', 'snap_%d.tsv' % i)
            if self.resume and prev and prev.get('status') == 'ok' and prev.get('key') == key and os.path.isfile(snap) and prev.get('chain'):
                cols, rows = read_tsv_text(open(snap).read())
                write_tsv(cat_path, cols, rows)
                state['steps'][lab] = dict(prev, status='ok', cached=True)
                res['steps_cached'] += 1
                chain = prev['chain']
                self.log('%s: cached' % lab)
                continue
            ok, msg, secs, cols, rows = self._exec(s, argv, cols, rows, cat_path)
            if not ok and self.retries:
                for k in range(self.retries):
                    self.log('%s: retry %d' % (lab, k + 1))
                    ok, msg, secs, cols, rows = self._exec(s, argv, cols, rows, cat_path)
                    if ok:
                        break
            if ok:
                h = hashlib.sha256(open(cat_path, 'rb').read()).hexdigest() if os.path.isfile(cat_path) else ''
                chain = hashlib.sha256((key + h).encode()).hexdigest()
                state['steps'][lab] = dict(status='ok', key=key, chain=chain, seconds=round(secs, 2), message=msg)
                import shutil
                shutil.copyfile(cat_path, os.path.join(self.dir, 'work', 'snap_%d.tsv' % i))
                res['steps_ok'] += 1
            else:
                state['steps'][lab] = dict(status='failed', key=key, seconds=round(secs, 2), message=msg)
                res.update(status='failed', failed_step=lab, message=msg); res['steps_failed'] += 1; broken = True
            atomic_json(stp, state)
        atomic_json(stp, state)
        res['n_objects'] = len(rows); res['n_columns'] = len(cols); res['seconds'] = round(time.time() - t_start, 2)
        self.cb(self.name, res['status'], '', len(plan), len(plan))
        self.logf.close()
        return res

    def _key_and_argv(self, s, chain, cols, rows, cat_path):
        if s['kind'] == 'detect':
            argv = [os.path.join(self.root, 'bin', 'ds9_sextract'), self.image] + [str(a) for a in self.recipe['detect'].get('args', [])]
            return hashlib.sha256((chain + json.dumps(argv)).encode()).hexdigest(), argv
        if s['kind'] == 'catalog':
            p = self.recipe['catalog'].replace('{name}', self.name)
            if not os.path.isfile(p):
                raise FileNotFoundError('catalog file %s not found' % p)
            return hashlib.sha256((chain + p + sha_file(p)).encode()).hexdigest(), [p]
        spec = s['spec']
        ctx = dict(python=self.python, work=os.path.join(self.dir, 'work'), image=self.image, catalog=os.path.join(self.dir, 'work', 'cat_in.tsv'), root=self.root)
        argv = cliexpand.build_argv(self.manifests, spec['plugin'], spec['step'], ctx, spec.get('params'))
        return hashlib.sha256((chain + json.dumps(argv)).encode()).hexdigest(), argv

    def _exec(self, s, argv, cols, rows, cat_path):
        lab = s['label']
        if s['kind'] == 'catalog':
            text = open(argv[0]).read()
            c, r = read_tsv_text(text)
            if not c or 'NUMBER' not in c:
                return False, 'catalog has no NUMBER column', 0.0, cols, rows
            write_tsv(cat_path, c, r)
            self.log('catalog: %d rows, %d columns' % (len(r), len(c)))
            return True, '%d rows' % len(r), 0.0, c, r
        if s['kind'] == 'detect':
            if not self.image or not os.path.isfile(self.image):
                return False, 'image not found: %r' % self.image, 0.0, cols, rows
            rc, out, err, dt = self.run_cmd(argv, 'detect')
            if rc != 0:
                return False, ('detect failed rc=%d: %s' % (rc, err.strip().splitlines()[-1] if err.strip() else '')), dt, cols, rows
            c, r = read_tsv_text(out)
            if not c or 'NUMBER' not in c:
                return False, 'detection produced no catalog', dt, cols, rows
            write_tsv(cat_path, c, r)
            return True, '%d sources' % len(r), dt, c, r
        # plugin step
        spec = s['spec']
        step = cliexpand.find_step(self.manifests[spec['plugin']], spec['step'])
        if not cols:
            return False, 'no catalog yet (put detect/catalog first)', 0.0, cols, rows
        os.makedirs(os.path.join(self.dir, 'work'), exist_ok=True)
        write_tsv(os.path.join(self.dir, 'work', 'cat_in.tsv'), cols, rows)
        rc, out, err, dt = self.run_cmd(argv, lab)
        if rc != 0:
            return False, ('rc=%d: %s' % (rc, err.strip().splitlines()[-1] if err.strip() else '')), dt, cols, rows
        mode = (step.get('output') or {}).get('mode')
        try:
            if mode == 'add_columns':
                c, r, added = merge_add_columns(cols, rows, out)
                write_tsv(cat_path, c, r)
                return True, 'added %d columns' % len(added), dt, c, r
            if mode == 'set':
                c, r = read_tsv_text(out)
                write_tsv(cat_path, c, r)
                return True, 'catalog replaced (%d rows)' % len(r), dt, c, r
            if mode == 'text':
                open(os.path.join(self.dir, 'work', '%s_output.txt' % lab.replace('.', '_')), 'w').write(out)
        except Exception as e:
            return False, 'bad output: %s' % e, dt, cols, rows
        write_tsv(cat_path, cols, rows)
        return True, 'ok', dt, cols, rows


def _worker(args):
    (name, image, recipe, outdir, root, python, resume, retries, timeout, status_path) = args
    manifests = cliexpand.load_manifests(root)

    def cb(n, state, step, i, tot):
        try:
            _status_update(status_path, n, dict(state=state, step=step, index=i, total=tot, t=time.time()))
        except Exception:
            pass
    try:
        f = Field(name, image, recipe, outdir, root, python, manifests, resume, retries, timeout, cb)
        return f.run()
    except Exception as e:
        traceback.print_exc()
        return dict(name=name, status='failed', steps_ok=0, steps_cached=0, steps_failed=1, failed_step='', message='%s: %s' % (type(e).__name__, e), n_objects=0, n_columns=0, seconds=0.0)


def _status_update(path, name, d):
    import fcntl
    lock = path + '.lock'
    with open(lock, 'w') as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            cur = json.load(open(path))
        except Exception:
            cur = {}
        cur[name] = d
        atomic_json(path, cur)


SUMMARY_COLS = ['name', 'status', 'n_objects', 'n_columns', 'steps_ok', 'steps_cached', 'steps_failed', 'failed_step', 'seconds', 'message']


def run_batch(fields, recipe, outdir, root, python=None, jobs=1, resume=False, retries=0, timeout=0, only=None):
    """-> list of per-field result dicts (also written to summary.tsv/json).  Fields are processed by `jobs` worker processes."""
    from multiprocessing import get_context
    python = python or sys.executable
    os.makedirs(outdir, exist_ok=True)
    status_path = os.path.join(outdir, 'status.json')
    sel = [(n, im) for n, im in fields if not only or n in only]
    atomic_json(status_path, {n: dict(state='queued', step='', index=0, total=0, t=time.time()) for n, _ in sel})
    args = [(n, im, recipe, outdir, root, python, resume, retries, timeout, status_path) for n, im in sel]
    if jobs <= 1:
        results = [_worker(a) for a in args]
    else:
        with get_context('spawn').Pool(jobs) as p:
            results = p.map(_worker, args, chunksize=1)
    prev = {}
    sp = os.path.join(outdir, 'summary.json')
    if only and os.path.isfile(sp):
        prev = {r['name']: r for r in json.load(open(sp))}
    for r in results:
        prev[r['name']] = r
    order = [n for n, _ in fields]
    allr = [prev[n] for n in order if n in prev]
    with open(os.path.join(outdir, 'summary.tsv'), 'w') as f:
        f.write('\t'.join(SUMMARY_COLS) + '\n')
        for r in allr:
            f.write('\t'.join(str(r.get(c, '')) for c in SUMMARY_COLS) + '\n')
    atomic_json(sp, allr)
    return results
