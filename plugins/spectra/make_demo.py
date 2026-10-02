#!/usr/bin/env python3
"""Write demonstration spectra for catalog objects:  make_demo.py OUTDIR NUMBER[,NUMBER...] [--seed 1]

One synthetic emission-line spectrum spec_<NUMBER>.fits (FITS table WAVE/FLUX/ERR, redshift 0.3-1.3 drawn from the seed) per object, plus truth.tsv (NUMBER Z).
Used by the GUI check and for trying the Spectra plugin without data."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from astropy.io import fits  # noqa: E402

from ogfkit import spectrasynth as ss  # noqa: E402


def main(argv):
    out = argv[1]
    nums = [int(v) for v in argv[2].split(',') if v.strip()]
    seed = int(argv[argv.index('--seed') + 1]) if '--seed' in argv else 1
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(seed)
    with open(os.path.join(out, 'truth.tsv'), 'w') as fh:
        fh.write('NUMBER\tZ\n')
        for n in nums:
            z = float(rng.uniform(0.3, 1.3))
            w, f, e, _ = ss.spectrum_1d(z=z, seed=seed * 1000 + n, scale=0.8)
            tb = fits.BinTableHDU.from_columns([fits.Column(name='WAVE', format='D', array=w), fits.Column(name='FLUX', format='D', array=f), fits.Column(name='ERR', format='D', array=e)])
            tb.writeto(os.path.join(out, 'spec_%d.fits' % n), overwrite=True)
            fh.write('%d\t%.5f\n' % (n, z))


if __name__ == '__main__':
    main(sys.argv)
