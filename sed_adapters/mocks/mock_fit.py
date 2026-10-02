#!/usr/bin/env python3
"""MOCK of a Bagpipes / Prospector fit script (test double): same JSON protocol as native/bagpipes_fit.py, TOY model grid (sed_adapters.toy) instead of the real library."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from sed_adapters import toy  # noqa: E402


def main():
    req = json.load(sys.stdin)
    rows = []
    for k in range(len(req['ids'])):
        fl = np.array(req['flux'][k]); er = np.array(req['err'][k])
        ok = (fl > -90) & (er > 0)
        bands = [b for b, o in zip(req['bands'], ok) if o]
        rows.append(toy.fit_sed(fl[ok], er[ok], bands, req['z'][k]))
    json.dump({'version': 'mock-fit 1.0 (toy)', 'rows': rows}, sys.stdout)


if __name__ == '__main__':
    main()
