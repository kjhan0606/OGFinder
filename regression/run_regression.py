#!/usr/bin/env python3
"""Regression-validation set: real public data, one or two cases per tool, headline numbers compared with baselines.json (lo/hi ranges).

    python3 regression/run_regression.py [--list] [--only a,b] [--tool t1,t2] [--report out.json] [--record]   # --record only prints the observed values

Exit status 0 = all run cases inside their ranges (cases whose data cannot be fetched are SKIP), 1 = at least one metric outside / case error.
Data cache: $OGF_DATA_CACHE (default ~/.cache/ogfinder_regression), fetched on demand (see regression/datasets.py, docs/testing.md)."""
import argparse
import json
import math
import os
import sys
import tempfile
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from regression import cases as C, datasets as D  # noqa: E402

BASE = os.path.join(HERE, 'baselines.json')


def compare(metrics, ranges):
    """-> (ok, [problems])"""
    bad = []
    for k, (lo, hi) in ranges.items():
        v = metrics.get(k)
        if v is None or (isinstance(v, float) and math.isnan(v)):
            bad.append('%s missing' % k)
        elif not (lo <= v <= hi):
            bad.append('%s = %.5g outside [%g, %g]' % (k, v, lo, hi))
    return not bad, bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--only', default='')
    ap.add_argument('--tool', default='')
    ap.add_argument('--report', default='')
    ap.add_argument('--record', action='store_true', help='print observed metrics (for writing baselines), no comparison')
    a = ap.parse_args(argv)
    base = json.load(open(BASE)) if os.path.exists(BASE) else {}
    names = [n for n in C.CASES if (not a.only or n in a.only.split(',')) and (not a.tool or C.CASES[n]['tool'] in a.tool.split(','))]
    if a.list:
        for n in names:
            print('%-28s %-12s %s' % (n, C.CASES[n]['tool'], 'baseline' if n in base else 'NO BASELINE'))
        return 0
    rep = {}
    nfail = 0
    for n in names:
        t0 = time.time()
        entry = dict(tool=C.CASES[n]['tool'])
        with tempfile.TemporaryDirectory(prefix='ogf_reg_') as w:
            try:
                m = C.CASES[n]['fn'](w)
                entry['metrics'] = {k: (float(v) if isinstance(v, (int, float)) else v) for k, v in m.items()}
                if a.record:
                    entry['status'] = 'RECORDED'
                elif n not in base:
                    entry['status'] = 'FAIL'; entry['problems'] = ['no baseline']
                else:
                    ok, bad = compare(entry['metrics'], base[n]['ranges'])
                    entry['status'] = 'PASS' if ok else 'FAIL'
                    entry['problems'] = bad
            except D.Unavailable as e:
                entry['status'] = 'SKIP'; entry['problems'] = [str(e)]
            except Exception as e:
                entry['status'] = 'FAIL'; entry['problems'] = ['%s: %s' % (type(e).__name__, str(e)[:300])]
                traceback.print_exc(file=sys.stderr)
        entry['seconds'] = round(time.time() - t0, 1)
        rep[n] = entry
        nfail += entry['status'] == 'FAIL'
        print('%-28s %-6s %6.0fs  %s' % (n, entry['status'], entry['seconds'], '; '.join(entry.get('problems', []))[:160] if entry['status'] != 'PASS' else
              ' '.join('%s=%.4g' % (k, v) for k, v in list(entry['metrics'].items())[:4] if isinstance(v, float))), flush=True)
        if a.record:
            print('   ', json.dumps(entry.get('metrics')), flush=True)
    if a.report:
        json.dump(rep, open(a.report, 'w'), indent=1)
    npass = sum(e['status'] == 'PASS' for e in rep.values())
    print('regression: %d pass, %d fail, %d skip of %d cases' % (npass, nfail, sum(e['status'] == 'SKIP' for e in rep.values()), len(rep)))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
