#!/usr/bin/env python3
"""Add pysersic (pysersic_map.py output) to a galfit_compare.py report for case A:  python compare_pysersic.py CMP.json PYS.json WORKDIR OUT.json [--galfit BIN --ld DIR]
chi^2 of the pysersic MAP solution is evaluated with GALFIT's renderer like the others (same data, same constant sigma)."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..'))
import galfit_compare as GC  # noqa: E402
from ogfkit import galfitio as GI  # noqa: E402


def main():
    cmp_json, pys_json, work, out = sys.argv[1:5]
    GC_args = type('A', (), {})()
    GC_args.galfit = sys.argv[sys.argv.index('--galfit') + 1] if '--galfit' in sys.argv else os.environ.get('GALFIT_BIN', 'galfit')
    GC_args.ld = sys.argv[sys.argv.index('--ld') + 1] if '--ld' in sys.argv else os.environ.get('GALFIT_LD', '')
    rep = json.load(open(cmp_json))
    pys = {p['seed']: p for p in json.load(open(pys_json)) if 'error' not in p}
    rows = []
    for r in rep['results']:
        if r['kind'] != 'A' or r['seed'] not in pys or not r.get('galfit') or not r.get('multifit'):
            continue
        p = pys[r['seed']]
        t = r['truth'][0]
        wd = os.path.join(work, 'A%02d' % r['seed'])
        cfg = dict(components=[dict(kind='sersic', x=p['x'], y=p['y'], mag=p['mag'], re=p['re'], n=p['n'], q=p['q'], pa=p['pa'] + 90.0, fixed=[], bounds={})], sky='const', sky_value=p['sky'], sky_grad=[0, 0], zp=GC.ZP, exptime=1.0)
        txt = GC.write_feedme_for_eval(cfg) if hasattr(GC, 'write_feedme_for_eval') else GI.write_feedme(cfg, image='data.fits', psf='psf.fits', shape=(GC.NY, GC.NX), region=[1, GC.NX, 1, GC.NY], mode=1)
        c2 = GC.eval_chi2(GC_args, wd, txt, 'p')
        row = dict(seed=r['seed'], pysersic_chi2=c2, galfit_chi2=r['galfit'].get('chi2'), multifit_chi2=r['multifit'].get('chi2'), truth_chi2=r.get('chi2_truth'))
        for who, f in (('pysersic', p), ('galfit', r['galfit']['comps'][0]), ('multifit', r['multifit']['comps'][0])):
            row[who] = dict(x=f['x'] - t['x'], y=f['y'] - t['y'], mag=f['mag'] - t['mag'], re_rel=f['re'] / t['re'] - 1, n_rel=f['n'] / t['n'] - 1, q=f['q'] - t['q'], pa=GC.wrap_pa(f['pa'] - t['pa']))
        row['time_pysersic'] = p['time']
        rows.append(row)
    summ = {}
    for who in ('pysersic', 'galfit', 'multifit'):
        summ[who] = {k: dict(median=float(np.median([r[who][k] for r in rows])), std=float(np.std([r[who][k] for r in rows])), p90_abs=float(np.percentile(np.abs([r[who][k] for r in rows]), 90))) for k in rows[0][who]}
    d = np.array([r['pysersic_chi2'] - r['galfit_chi2'] for r in rows if r['pysersic_chi2'] and r['galfit_chi2']])
    summ['chi2_pysersic_minus_galfit'] = dict(n=int(d.size), median=float(np.median(d)), p10=float(np.percentile(d, 10)), p90=float(np.percentile(d, 90)), max=float(d.max()), min=float(d.min()),
                                              pysersic_worse_by_gt1=float(np.mean(d > 1.0)))
    summ['time_pysersic_median'] = float(np.median([r['time_pysersic'] for r in rows]))
    json.dump(dict(rows=rows, summary=summ), open(out, 'w'))
    print('n', len(rows))
    for who in ('pysersic', 'galfit', 'multifit'):
        print(who, ' '.join('%s %+.3f/%.3f' % (k, v['median'], v['p90_abs']) for k, v in summ[who].items()))
    print(summ['chi2_pysersic_minus_galfit'], 'pysersic time median %.1f s (CPU, jax)' % summ['time_pysersic_median'])


if __name__ == '__main__':
    main()
