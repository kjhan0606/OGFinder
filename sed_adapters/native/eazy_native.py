#!/usr/bin/env python3
"""Runs eazy-py on the photometry in the JSON on stdin (uJy fluxes, AB 23.9); prints {"version":..., "rows":[...]}.  Needs `pip install eazy` + the eazy-photoz repository."""
import json
import os
import sys
import tempfile

import numpy as np


def main():
    req = json.load(sys.stdin)
    wd = req['wd']
    os.chdir(wd)
    if req.get('eazy_data'):
        os.environ['EAZYCODE'] = req['eazy_data']
    import eazy
    import eazy.filters as F
    import eazy.photoz as P
    from astropy.table import Table
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    from sed_adapters import filters as reg, common
    res = F.FilterFile(req['filters_res'])
    names = [f.name for f in res.filters]
    bands = req['bands']
    flux = np.array(req['flux']); err = np.array(req['err'])
    cat = Table()
    cat['id'] = np.arange(1, len(req['ids']) + 1)
    tr, idx = [], []
    for j, b in enumerate(bands):
        k = reg.eazy_filter_index(reg.B[b], names)
        if k is None:
            continue
        cat['f_' + b] = flux[:, j]; cat['e_' + b] = err[:, j]
        tr += ['f_%s F%d' % (b, k), 'e_%s E%d' % (b, k)]
        idx.append(k)
    cat.write('catalog.fits', overwrite=True)
    open('zphot.translate', 'w').write('\n'.join(tr) + '\n')
    prior_filt = idx[-1]
    params = {'CATALOG_FILE': 'catalog.fits', 'MAIN_OUTPUT_FILE': 'photz', 'FILTERS_RES': req['filters_res'], 'TEMPLATES_FILE': req['templates_file'],
              'Z_MIN': req['z_min'], 'Z_MAX': req['z_max'], 'Z_STEP': req['z_step'], 'PRIOR_ABZP': 23.9, 'PRIOR_FILTER': prior_filt,
              'APPLY_PRIOR': 'y' if req.get('prior') else 'n', 'CAT_HAS_EXTCORR': 'n', 'MW_EBV': 0.0}
    if req.get('temp_err_file'):
        params['TEMP_ERR_FILE'] = req['temp_err_file']
    elif req.get('eazy_data'):
        te = os.path.join(req['eazy_data'], 'templates', 'TEMPLATE_ERROR.eazy_v1.0')
        if os.path.exists(te):
            params['TEMP_ERR_FILE'] = te
    ez = P.PhotoZ(param_file=None, translate_file='zphot.translate', zeropoint_file=None, params=params, load_prior=bool(req.get('prior')), load_products=False)
    ez.fit_catalog(n_proc=0, get_best_fit=True)
    zg = np.asarray(ez.zgrid)
    rows = []
    zb = np.asarray(ez.zbest)
    for i in range(len(req['ids'])):
        lnp = np.asarray(ez.lnp[i])
        if not np.all(np.isfinite(lnp)) or lnp.max() <= -1e20:
            rows.append(dict(z=None, error='no valid p(z)')); continue
        s = common.pdf_summary(zg, np.exp(lnp - lnp.max()))
        n = max(int(np.sum(flux[i] > -90)) - 1, 1)
        rows.append(dict(z=float(zb[i]) if np.isfinite(zb[i]) else s['z_peak'], p16=s['p16'], p50=s['p50'], p84=s['p84'], sigma=s['z_sigma'], chi2=float(ez.chi2_best[i]) / n))
    json.dump({'version': eazy.__version__, 'rows': rows}, sys.stdout)


if __name__ == '__main__':
    main()
