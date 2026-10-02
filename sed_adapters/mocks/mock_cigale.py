#!/usr/bin/env python3
"""MOCK of `pcigale run` (test double, NOT CIGALE): reads pcigale.ini + the observations table from the current directory, fits the TOY grid, writes out/results.txt with
CIGALE's column naming (bayes.stellar.m_star, bayes.sfh.sfr, ..., best.reduced_chi_square)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from sed_adapters import filters, toy  # noqa: E402


def main():
    ini = {}
    for ln in open('pcigale.ini'):
        if '=' in ln and not ln.strip().startswith('['):
            k, v = ln.split('=', 1)
            ini[k.strip()] = v.strip()
    cols, rows = None, []
    for ln in open(ini['data_file']):
        if ln.startswith('#'):
            cols = ln[1:].split(); continue
        if ln.strip():
            rows.append(dict(zip(cols, ln.split())))
    name2band = {e['cigale']: b for b, e in filters.B.items()}
    os.makedirs('out', exist_ok=True)
    head = ['id', 'bayes.stellar.m_star', 'bayes.stellar.m_star_err', 'bayes.sfh.sfr', 'bayes.sfh.sfr_err', 'bayes.sfh.age_main', 'bayes.sfh.age_main_err', 'bayes.dust.e_bv_lines',
            'bayes.stellar.metallicity', 'best.reduced_chi_square']
    out = ['# ' + ' '.join(head)]
    for r in rows:
        z = float(r['redshift'])
        bands, fl, er = [], [], []
        for c in cols:
            if c in name2band and r[c] != 'nan':
                bands.append(name2band[c]); fl.append(float(r[c]) * 1e3); er.append(float(r[c + '_err']) * 1e3)
        res = toy.fit_sed(np.array(fl), np.array(er), bands, z)
        m = 10 ** res['LOG_MASS']
        out.append(' '.join([r['id'], '%.6g' % m, '%.6g' % (m * np.log(10) * res['LOG_MASS_ERR']), '%.6g' % res['SFR'], '%.6g' % (0.1 * res['SFR']),
                             '%.6g' % (10 ** res['LOG_AGE'] / 1e6), '%.6g' % (10 ** res['LOG_AGE'] / 1e6 * np.log(10) * res['LOG_AGE_ERR']), '%.6g' % (res['AV'] / 4.05),
                             '%.6g' % (0.02 * 10 ** res['LOG_Z']), '%.6g' % res['SED_CHI2']]))
    open('out/results.txt', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
