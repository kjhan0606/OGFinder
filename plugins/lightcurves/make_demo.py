#!/usr/bin/env python3
"""Write a demo long-form light-curve table for catalog objects: make_demo.py OUTDIR NUMBERS(comma list) [--seed 3] -> OUTDIR/lc_demo.tsv, truth.tsv.
Classes cycle through the synthetic library; `number` links each light curve to a catalog NUMBER."""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from ogfkit import lcsynth  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('numbers')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--mag', action='store_true', help='write magnitudes instead of fluxes')
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.RandomState(a.seed)
    nums = [n for n in a.numbers.split(',') if n]
    with open(os.path.join(a.outdir, 'lc_demo.tsv'), 'w') as f, open(os.path.join(a.outdir, 'truth.tsv'), 'w') as g:
        f.write('id\tnumber\tmjd\t%s\t%s\thost_offset_re\n' % (('mag', 'mag_err') if a.mag else ('flux', 'flux_err')))
        g.write('number\tclass\n')
        for i, n in enumerate(nums):
            c = lcsynth.CLASSES[i % len(lcsynth.CLASSES)]
            lc = lcsynth.make_lc(c, rng)
            g.write('%s\t%s\n' % (n, c))
            for t, fl, e in zip(lc['t'], lc['flux'], lc['err']):
                if a.mag:
                    if fl <= 0:
                        continue
                    f.write('LC%s\t%s\t%.4f\t%.4f\t%.4f\t%s\n' % (n, n, 59000 + t, lcsynth.flux2mag(fl), 1.0857 * e / fl, '' if lc['host_offset_re'] is None else '%.3f' % lc['host_offset_re']))
                else:
                    f.write('LC%s\t%s\t%.4f\t%.5g\t%.5g\t%s\n' % (n, n, 59000 + t, fl, e, '' if lc['host_offset_re'] is None else '%.3f' % lc['host_offset_re']))


if __name__ == '__main__':
    main()
