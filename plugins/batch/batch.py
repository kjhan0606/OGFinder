#!/usr/bin/env python3
"""Run a recipe of plugin steps over many fields (parallel workers, resume, per-field logs, summary table).

    batch.py --fields fields.txt --recipe recipe.json --outdir out [--jobs 4] [--resume] [--retries 1] [--timeout 600] [--only A,B] [--dry-run]
    batch.py --example-recipe recipe.json

fields.txt: "NAME path/to/image.fits" per line.  recipe.json: see ogfkit/batchrun.py (detect | catalog, then steps of any plugin that has a cli template).
Outputs: <out>/<NAME>/{catalog.tsv, state.json, field.log, work/}, <out>/summary.tsv|json, <out>/status.json (live progress; the GUI reads it).
Exit status 1 when at least one field failed.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit import batchrun, cliexpand  # noqa: E402

EXAMPLE = {
    "detect": {"args": ["--detect-thresh", "3.0", "--detect-minarea", "5"]},
    "steps": [
        {"plugin": "noisemodel", "step": "model"},
        {"plugin": "cluster", "step": "density", "params": {"dens-sigma": 150}},
    ],
}


def main(argv=None):
    import signal
    signal.signal(signal.SIGTERM, lambda *a: sys.exit(143))      # Cancel in the GUI: the worker pool is terminated on exit
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fields')
    ap.add_argument('--recipe')
    ap.add_argument('--outdir')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--retries', type=int, default=0)
    ap.add_argument('--timeout', type=float, default=0.0, help='seconds per step (0 = none)')
    ap.add_argument('--only', default='')
    ap.add_argument('--python', default=sys.executable)
    ap.add_argument('--root', default=ROOT)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--example-recipe', default='')
    ap.add_argument('--expand', default='', help='PLUGIN.STEP: print the expanded argv as JSON (with --ctx-json / --params-json)')
    ap.add_argument('--ctx-json', default='{}')
    ap.add_argument('--params-json', default='{}')
    a = ap.parse_args(argv)
    if a.expand:
        man = cliexpand.load_manifests(a.root)
        pid, sid = a.expand.split('.', 1)
        ctx = dict(root=a.root, python=a.python, work='', image='', catalog='')
        ctx.update(json.loads(a.ctx_json))
        print(json.dumps(cliexpand.build_argv(man, pid, sid, ctx, json.loads(a.params_json))))
        return 0
    if a.example_recipe:
        json.dump(EXAMPLE, open(a.example_recipe, 'w'), indent=1)
        print('wrote', a.example_recipe)
        return 0
    if not (a.fields and a.recipe and a.outdir):
        ap.error('--fields, --recipe and --outdir are required')
    if not os.path.isfile(a.fields):
        raise SystemExit('batch: fields file not found: %s' % a.fields)
    if not os.path.isfile(a.recipe):
        raise SystemExit('batch: recipe not found: %s' % a.recipe)
    fields = batchrun.read_fields(a.fields)
    recipe = json.load(open(a.recipe))
    man = cliexpand.load_manifests(a.root)
    for s in recipe.get('steps', []):          # validate the recipe before anything runs
        if s.get('plugin') not in man:
            raise SystemExit('batch: unknown plugin %r in recipe' % s.get('plugin'))
        try:
            st = cliexpand.find_step(man[s['plugin']], s['step'])
        except KeyError as e:
            raise SystemExit('batch: %s' % e)
        if 'cli' not in st:
            raise SystemExit('batch: step %s.%s has no cli template and cannot run headless' % (s['plugin'], s['step']))
        unknown = [k for k in (s.get('params') or {}) if k not in cliexpand.defaults(man[s['plugin']])]
        if unknown:
            raise SystemExit('batch: unknown parameter(s) %s for %s' % (', '.join(unknown), s['plugin']))
    only = [x for x in a.only.split(',') if x]
    if a.dry_run:
        for n, im in fields:
            if only and n not in only:
                continue
            f = batchrun.Field.__new__(batchrun.Field)
            f.recipe = recipe
            print(n, im, '->', ' | '.join(s['label'] for s in batchrun.Field.plan(f)))
        return 0
    res = batchrun.run_batch(fields, recipe, a.outdir, a.root, a.python, a.jobs, a.resume, a.retries, a.timeout, only)
    print('%-16s %-8s %8s %5s %3s %3s %3s %8s  %s' % ('field', 'status', 'objects', 'cols', 'ok', 'cac', 'bad', 'seconds', 'failed step / message'))
    for r in res:
        print('%-16s %-8s %8d %5d %3d %3d %3d %8.1f  %s' % (r['name'], r['status'], r['n_objects'], r['n_columns'], r['steps_ok'], r['steps_cached'], r['steps_failed'], r['seconds'],
                                                              (r['failed_step'] + ': ' + r['message']) if r['status'] != 'ok' else ''))
    print('SUMMARY %s  ok=%d failed=%d' % (os.path.join(a.outdir, 'summary.tsv'), sum(r['status'] == 'ok' for r in res), sum(r['status'] != 'ok' for r in res)))
    return 1 if any(r['status'] != 'ok' for r in res) else 0


if __name__ == '__main__':
    sys.exit(main())
