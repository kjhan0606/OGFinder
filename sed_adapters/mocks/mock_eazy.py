#!/usr/bin/env python3
"""MOCK of the classic EAZY executable (test double, NOT EAZY): `mock_eazy.py zphot.param` reads the EAZY catalog / translate / param files and the FILTER.RES.info
the param file points to, fits the TOY templates (sed_adapters.toy) on a redshift grid and writes `<MAIN_OUTPUT_FILE>.zout` in EAZY's column layout."""
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from sed_adapters import filters, toy, common  # noqa: E402


def read_param(path):
    d = {}
    for ln in open(path):
        t = ln.split(None, 1)
        if len(t) == 2 and not ln.startswith('#'):
            d[t[0]] = t[1].strip()
    return d


def main():
    P = read_param(sys.argv[1])
    info = P['FILTERS_RES'] + '.info'
    names = {}
    for ln in open(info):
        m = re.match(r'\s*(\d+)\s+(\S+)', ln)
        if m:
            names[int(m.group(1))] = m.group(2)
    nl = [names[k] for k in sorted(names)]
    idx2band = {}
    for b, e in filters.B.items():
        k = filters.eazy_filter_index(e, nl)
        if k is not None:
            idx2band[k] = b
    tr = {}
    for ln in open(P.get('TRANSLATE_FILE', 'zphot.translate')):
        t = ln.split()
        if len(t) == 2:
            tr[t[0]] = t[1]
    cols = None
    rows = []
    for ln in open(P['CATALOG_FILE']):
        if ln.startswith('#'):
            cols = ln[1:].split(); continue
        if ln.strip():
            rows.append(ln.split())
    zg = np.arange(float(P['Z_MIN']), float(P['Z_MAX']) + 1e-9, float(P['Z_STEP']))
    out = ['# id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68 l95 u95 l99 u99 nfilt q_z z_peak peak_prob z_mc']
    for r in rows:
        d = dict(zip(cols, r))
        bands, fl, er = [], [], []
        for c in cols:
            if c.startswith('F_') or tr.get(c, '').startswith('F') and c in tr:
                ecol = 'E_' + c[2:]
                fk = tr.get(c); k = int(fk[1:]) if fk else None
                f, e = float(d[c]), float(d.get(ecol, -99))
                if k in idx2band and f > -90 and e > 0:
                    bands.append(idx2band[k]); fl.append(f); er.append(e)
        if len(bands) < int(P.get('N_MIN_COLORS', 3)):
            out.append(' '.join([d['id'], '-1'] + ['-99'] * 18)); continue
        fl, er = np.array(fl), np.array(er)
        chi, _ = toy.fit_photoz(fl, er, bands, zg)
        s = common.pdf_summary(zg, np.exp(-0.5 * (chi - chi.min())))
        j = int(np.argmin(chi))
        out.append(' '.join([d['id'], d.get('z_spec', '-1'), '%.4f' % zg[j], '%.4f' % s['z_mean'], '%.4f' % chi.min(), '%.4f' % s['z_peak'], '%.4f' % chi[j], '%.4f' % s['z_mean'], '1.0',
                             '%.4f' % s['p16'], '%.4f' % s['p84'], '%.4f' % s['p16'], '%.4f' % s['p84'], '%.4f' % s['p16'], '%.4f' % s['p84'], str(len(bands)), '0.0', '%.4f' % s['z_peak'], '1.0', '%.4f' % s['z_peak']]))
    open(P['MAIN_OUTPUT_FILE'] + '.zout', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
