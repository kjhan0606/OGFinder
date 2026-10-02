#!/usr/bin/env python3
"""Reproducibility bundles (zip) and their verification.

    repro.py bundle  --out B.zip --kind session|catalog --catalog C.tsv [--inputs a.fits,b.fits] [--params P.json] [--steps S.json] [--session-script S.py] [--include-inputs] [--note TXT]
    repro.py bundle-batch --out B.zip --batch-out DIR --fields F.txt --recipe R.json [--include-inputs]
    repro.py verify  B.zip [--inputs-dir DIR] [--rtol 0] [--atol 0] [--strict-env] [--no-rerun] [--jobs N] [--json OUT.json] [--workdir DIR] [-- <extra args for the session script>]
    repro.py info    B.zip

verify exit status: 0 = every check PASS/WARN/SKIP, 1 = at least one FAIL.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from ogfkit import repro, cliexpand  # noqa: E402


def cmd_bundle(a):
    inputs = [(os.path.basename(p), p) for p in a.inputs.split(',') if p]
    params = json.load(open(a.params)) if a.params else {}
    steps = json.load(open(a.steps)) if a.steps else []
    outs = [('catalog_final.tsv', a.catalog)] if a.catalog else []
    man = repro.make_bundle(a.out, a.kind, a.root, inputs, outs, params, steps, session_script=a.session_script or None, include_inputs=a.include_inputs, note=a.note,
                            verify_args=a.verify_args)
    print('bundle %s: kind %s, %d inputs, %d outputs, %d key packages, git %s%s' % (a.out, man['kind'], len(man['inputs']), len(man['outputs']), len(man['packages']),
                                                                                    (man['tool'].get('head') or 'n/a')[:9], ' (dirty)' if man['tool'].get('dirty') else ''))
    return 0


def cmd_bundle_batch(a):
    recipe = json.load(open(a.recipe))
    fields_txt = open(a.fields).read()
    from ogfkit import batchrun
    fields = batchrun.read_fields(a.fields)
    man = cliexpand.load_manifests(a.root)
    params = {}
    for s in recipe.get('steps', []):
        d = cliexpand.defaults(man[s['plugin']])
        d.update(s.get('params') or {})
        params['%s.%s' % (s['plugin'], s['step'])] = d
    steps = []
    for n, im in fields:
        sp = os.path.join(a.batch_out, n, 'state.json')
        if os.path.isfile(sp):
            st = json.load(open(sp))
            steps.append(dict(field=n, steps={k: {kk: vv for kk, vv in v.items() if kk in ('status', 'seconds', 'message')} for k, v in st.get('steps', {}).items()}))
    inputs = [(os.path.basename(im), im) for n, im in fields if im and os.path.isfile(im)]
    outputs = [('%s/catalog.tsv' % n, os.path.join(a.batch_out, n, 'catalog.tsv')) for n, _ in fields]
    mo = repro.make_bundle(a.out, 'batch', a.root, inputs, outputs, params, steps, recipe=recipe, fields_txt=fields_txt, include_inputs=a.include_inputs, note=a.note)
    print('bundle %s: kind batch, %d fields, %d inputs, %d outputs' % (a.out, len(fields), len(mo['inputs']), len(mo['outputs'])))
    return 0


def cmd_verify(a, extra):
    rep = repro.verify_bundle(a.bundle, a.root, a.python, a.inputs_dir, a.workdir, a.rtol, a.atol, a.strict_env, not a.no_rerun, a.jobs, extra)
    for c in rep['checks']:
        print('%-5s %-34s %s' % (c['status'], c['name'], c['detail']))
    nf = sum(c['status'] == 'FAIL' for c in rep['checks']); nw = sum(c['status'] == 'WARN' for c in rep['checks'])
    print('VERIFY %s  checks=%d failures=%d warnings=%d' % ('PASS' if rep['ok'] else 'FAIL', len(rep['checks']), nf, nw))
    if a.json:
        json.dump(rep, open(a.json, 'w'), indent=1)
    return 0 if rep['ok'] else 1


def cmd_info(a):
    m = repro.read_manifest(a.bundle)
    print(json.dumps({k: m[k] for k in ('schema', 'kind', 'created', 'tool', 'host', 'packages', 'inputs', 'outputs', 'verify') if k in m}, indent=1))
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    extra = []
    if '--' in argv:
        i = argv.index('--'); extra = argv[i + 1:]; argv = argv[:i]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('bundle')
    b.add_argument('--out', required=True); b.add_argument('--kind', default='session', choices=['session', 'catalog'])
    b.add_argument('--catalog', default=''); b.add_argument('--inputs', default=''); b.add_argument('--params', default=''); b.add_argument('--steps', default='')
    b.add_argument('--session-script', default=''); b.add_argument('--include-inputs', action='store_true'); b.add_argument('--note', default=''); b.add_argument('--verify-args', default='')
    bb = sp.add_parser('bundle-batch')
    bb.add_argument('--out', required=True); bb.add_argument('--batch-out', required=True); bb.add_argument('--fields', required=True); bb.add_argument('--recipe', required=True)
    bb.add_argument('--include-inputs', action='store_true'); bb.add_argument('--note', default='')
    v = sp.add_parser('verify')
    v.add_argument('bundle'); v.add_argument('--inputs-dir', default=None); v.add_argument('--rtol', type=float, default=0.0); v.add_argument('--atol', type=float, default=0.0)
    v.add_argument('--strict-env', action='store_true'); v.add_argument('--no-rerun', action='store_true'); v.add_argument('--jobs', type=int, default=1)
    v.add_argument('--json', default=''); v.add_argument('--workdir', default=None); v.add_argument('--python', default=None)
    i = sp.add_parser('info'); i.add_argument('bundle')
    for p in (b, bb, v, i):
        p.add_argument('--root', default=ROOT)
    a = ap.parse_args(argv)
    if a.cmd == 'bundle':
        return cmd_bundle(a)
    if a.cmd == 'bundle-batch':
        return cmd_bundle_batch(a)
    if a.cmd == 'verify':
        return cmd_verify(a, extra)
    return cmd_info(a)


if __name__ == '__main__':
    sys.exit(main())
