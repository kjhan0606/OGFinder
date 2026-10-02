"""Run with the interpreter that has eazy-py: photometry of EAZY's own templates (single template per object, redshifted through the FILTER.RES curves) -> JSON on stdout."""
import json
import os
import sys

import numpy as np

root = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
os.environ['EAZYCODE'] = root
import eazy.filters as F
import eazy.templates as T

res = F.FilterFile(os.path.join(root, 'filters', 'FILTER.RES.latest'))
ids = {'F435W': 233, 'F606W': 236, 'F775W': 238, 'F850LP': 240, 'F105W': 202, 'F125W': 203, 'F160W': 205}
tpl = T.read_templates_file(os.path.join(root, 'templates', 'eazy_v1.3.spectra.param'))[:7]
rng = np.random.default_rng(7)
out = []
for k in range(n):
    z = float(rng.uniform(0.2, 3.0)); j = int(rng.integers(0, len(tpl))); mag = float(rng.uniform(21.0, 23.5))
    fl = {b: float(tpl[j].integrate_filter(res.filters[i - 1], scale=1.0, z=z, include_igm=True)) for b, i in ids.items()}
    sc = 10 ** (-0.4 * (mag - 23.9)) / fl['F160W']
    mags, errs = {}, {}
    for b in ids:
        f = fl[b] * sc; e = max(0.04 * f, 0.02)
        f = f + rng.normal() * e
        mags[b] = float(23.9 - 2.5 * np.log10(f)) if f > 0 else None
        errs[b] = float(1.0857 * e / f) if f > 0 else None
    out.append(dict(id=str(k + 1), z=z, template=j, mags=mags, mag_errs=errs))
json.dump(out, sys.stdout)
