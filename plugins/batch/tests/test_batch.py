"""Batch runner tests: headless template expansion, parallel == serial, resume (cached / changed parameter), failure isolation, timeout, throughput on many fields."""
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUG = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import cliexpand, batchrun  # noqa: E402

CLI = os.path.join(PLUG, 'batch.py')
SEXT = os.path.join(ROOT, 'ogfmeas', 'sextract.py')
need_sext = pytest.mark.skipif(not os.path.exists(SEXT), reason='ogfmeas/sextract.py missing')


def run(args):
    return subprocess.run([sys.executable, CLI, '--root', ROOT] + args, capture_output=True, text=True)


def make_fields(tmp, n, size=300):
    from astropy.io import fits
    lines = []
    for k in range(n):
        rng = np.random.RandomState(k)
        im = rng.randn(size, size) * 5 + 100
        yy, xx = np.mgrid[0:size, 0:size]
        for i in range(25):
            x, y = rng.rand(2) * (size - 40) + 20
            im += rng.uniform(150, 900) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 2.0 ** 2))
        p = os.path.join(tmp, 'f%d.fits' % k)
        fits.writeto(p, im.astype('float32'), overwrite=True)
        lines.append('F%02d %s\n' % (k, p))
    open(os.path.join(tmp, 'fields.txt'), 'w').write(''.join(lines))
    return os.path.join(tmp, 'fields.txt')


def cat_hash(out, name):
    return hashlib.sha256(open(os.path.join(out, name, 'catalog.tsv'), 'rb').read()).hexdigest()


def test_all_cli_templates_expand_headlessly():
    man = cliexpand.load_manifests(ROOT)
    n = 0
    for pid, m in man.items():
        for s in m['steps']:
            if 'cli' in s:
                ctx = dict(python='PY', work='/w', image='/i.fits', catalog='/c.tsv', root=ROOT)
                argv = cliexpand.build_argv(man, pid, s['id'], ctx)
                assert argv and argv[0] == 'PY' and all(isinstance(a, str) for a in argv)
                n += 1
    print('expanded %d cli templates' % n)
    assert n >= 30
    with pytest.raises(KeyError):
        cliexpand.find_step(man['cluster'], 'nope')

def test_overrides_and_conditions():
    man = cliexpand.load_manifests(ROOT)
    ctx = dict(python='PY', work='/w', image='/i.fits', catalog='/c.tsv', root=ROOT)
    a0 = cliexpand.build_argv(man, 'xmatch', 'match', ctx)
    a1 = cliexpand.build_argv(man, 'xmatch', 'match', ctx, {'apply-shift': True, 'radius-arcsec': 2.5})
    assert '--apply-shift' not in a0 and '--apply-shift' in a1
    assert a1[a1.index('--radius-arcsec') + 1] == '2.5'


@need_sext
def test_serial_equals_parallel_and_resume(tmp_path):
    fl = make_fields(str(tmp_path), 6)
    recipe = {"detect": {"args": ["--detect-thresh", "3"]}, "steps": [{"plugin": "noisemodel", "step": "model"}, {"plugin": "example_hello", "step": "greet"}]}
    rp = tmp_path / 'recipe.json'; rp.write_text(json.dumps(recipe))
    t0 = time.time(); r1 = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o1'), '--jobs', '1', '--python', sys.executable]); t1 = time.time() - t0
    t0 = time.time(); r3 = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o3'), '--jobs', '3', '--python', sys.executable]); t3 = time.time() - t0
    assert r1.returncode == 0 and r3.returncode == 0, r1.stdout + r1.stderr + r3.stdout + r3.stderr
    print('6 fields x 3 steps: jobs=1 %.1fs, jobs=3 %.1fs (speed-up x%.2f)' % (t1, t3, t1 / t3))
    for k in range(6):
        assert cat_hash(str(tmp_path / 'o1'), 'F%02d' % k) == cat_hash(str(tmp_path / 'o3'), 'F%02d' % k)
    s = json.load(open(tmp_path / 'o3' / 'summary.json'))
    assert len(s) == 6 and all(r['status'] == 'ok' and r['steps_ok'] == 3 and r['n_objects'] >= 20 for r in s)
    assert os.path.getsize(tmp_path / 'o3' / 'F00' / 'field.log') > 100
    # resume: nothing recomputed
    r = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o3'), '--jobs', '3', '--resume', '--python', sys.executable])
    assert r.returncode == 0
    s = json.load(open(tmp_path / 'o3' / 'summary.json'))
    print('resume: cached steps per field', [r_['steps_cached'] for r_ in s])
    assert all(r_['steps_cached'] == 3 and r_['steps_ok'] == 0 for r_ in s)
    assert cat_hash(str(tmp_path / 'o3'), 'F00') == cat_hash(str(tmp_path / 'o1'), 'F00')
    # an explicit parameter equal to its default does not invalidate anything
    names = cliexpand.defaults(cliexpand.load_manifests(ROOT)['noisemodel'])
    recipe['steps'][0]['params'] = {'bw': names['bw']}
    rp.write_text(json.dumps(recipe))
    run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o3'), '--jobs', '3', '--resume', '--python', sys.executable])
    s = json.load(open(tmp_path / 'o3' / 'summary.json'))
    assert all(r_['steps_cached'] == 3 for r_ in s)
    recipe['steps'][0]['params'] = {'bw': 32}
    rp.write_text(json.dumps(recipe))
    run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o3'), '--jobs', '3', '--resume', '--python', sys.executable])
    s = json.load(open(tmp_path / 'o3' / 'summary.json'))
    print('changed noisemodel bw: ran %s cached %s' % ([r_['steps_ok'] for r_ in s], [r_['steps_cached'] for r_ in s]))
    assert all(r_['steps_cached'] == 1 and r_['steps_ok'] == 2 for r_ in s)
    recipe['steps'][0]['params'] = {}
    # detection threshold changed -> everything downstream recomputed
    recipe['detect']['args'] = ["--detect-thresh", "4"]
    rp.write_text(json.dumps(recipe))
    r = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o3'), '--jobs', '3', '--resume', '--python', sys.executable])
    s = json.load(open(tmp_path / 'o3' / 'summary.json'))
    print('changed detection threshold: ran', [r_['steps_ok'] for r_ in s], 'cached', [r_['steps_cached'] for r_ in s])
    assert all(r_['steps_ok'] == 3 and r_['steps_cached'] == 0 for r_ in s)


@need_sext
def test_failure_isolation_and_timeout(tmp_path):
    fl = make_fields(str(tmp_path), 3)
    lines = open(fl).read().splitlines()
    lines[1] = 'F01 /nonexistent/missing.fits'
    open(fl, 'w').write('\n'.join(lines) + '\n')
    recipe = {"detect": {"args": ["--detect-thresh", "3"]}, "steps": [{"plugin": "example_hello", "step": "greet"}]}
    rp = tmp_path / 'r.json'; rp.write_text(json.dumps(recipe))
    r = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o'), '--jobs', '2', '--python', sys.executable])
    assert r.returncode == 1
    s = {x['name']: x for x in json.load(open(tmp_path / 'o' / 'summary.json'))}
    assert s['F01']['status'] == 'failed' and s['F01']['failed_step'] == 'detect' and 'not found' in s['F01']['message']
    assert s['F00']['status'] == 'ok' and s['F02']['status'] == 'ok'
    # timeout
    r = run(['--fields', fl, '--recipe', str(rp), '--outdir', str(tmp_path / 'o2'), '--jobs', '1', '--timeout', '0.01', '--only', 'F00', '--python', sys.executable])
    s = json.load(open(tmp_path / 'o2' / 'summary.json'))
    assert s[0]['status'] == 'failed' and 'rc=124' in s[0]['message'] or 'timeout' in s[0]['message']


def test_recipe_validation_and_dry_run(tmp_path):
    (tmp_path / 'f.txt').write_text('A /x.fits\nB /y.fits\n')
    def rr(steps):
        (tmp_path / 'r.json').write_text(json.dumps({"catalog": "/c.tsv", "steps": steps}))
        return run(['--fields', str(tmp_path / 'f.txt'), '--recipe', str(tmp_path / 'r.json'), '--outdir', str(tmp_path / 'o')])
    assert 'unknown plugin' in rr([{"plugin": "nope", "step": "x"}]).stderr
    assert 'no cli template' in rr([{"plugin": "cluster", "step": "plot"}]).stderr or 'no step' in rr([{"plugin": "cluster", "step": "plot"}]).stderr
    assert 'unknown parameter' in rr([{"plugin": "xmatch", "step": "match", "params": {"bogus": 1}}]).stderr
    rr([{"plugin": "example_hello", "step": "greet"}])
    r = run(['--fields', str(tmp_path / 'f.txt'), '--recipe', str(tmp_path / 'r.json'), '--outdir', str(tmp_path / 'o'), '--dry-run'])
    assert r.returncode == 0 and 'A /x.fits -> catalog' in r.stdout


def test_many_fields_throughput(tmp_path):
    """60 catalog-only fields, 6 workers: no failures; throughput and the parallel speed-up are reported."""
    n = 60
    cat = tmp_path / 'c.tsv'
    cat.write_text('NUMBER\tX_IMAGE\tY_IMAGE\n' + ''.join('%d\t%d\t%d\n' % (i, i, i) for i in range(1, 51)))
    (tmp_path / 'f.txt').write_text(''.join('S%03d /none.fits\n' % k for k in range(n)))
    (tmp_path / 'r.json').write_text(json.dumps({"catalog": str(cat), "steps": [{"plugin": "example_hello", "step": "greet"}]}))
    res = {}
    for jobs in (1, 6):
        t0 = time.time()
        r = run(['--fields', str(tmp_path / 'f.txt'), '--recipe', str(tmp_path / 'r.json'), '--outdir', str(tmp_path / ('o%d' % jobs)), '--jobs', str(jobs), '--python', sys.executable])
        res[jobs] = time.time() - t0
        assert r.returncode == 0, r.stdout[-500:]
    s = json.load(open(tmp_path / 'o6' / 'summary.json'))
    assert len(s) == n and all(x['status'] == 'ok' and x['n_objects'] == 50 for x in s)
    print('%d fields: jobs=1 %.1fs (%.1f fields/s), jobs=6 %.1fs (%.1f fields/s)' % (n, res[1], n / res[1], res[6], n / res[6]))
    assert res[6] < res[1]
