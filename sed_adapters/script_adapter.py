"""JSON protocol shared by the SED-fit adapters.

The adapter prepares photometry and runs a fit script (stdin request -> stdout {"version", "rows"}). ``params.command`` replaces the shipped script. The shipped scripts do not call an external stellar-population library. The in-tree SED model is analytic. A caller-supplied command is how the mock executables are tested."""
import json
import os
import sys
import tempfile

import numpy as np

from . import common, filters

HERE = os.path.dirname(os.path.abspath(__file__))


def redshift_of(record, params):
    zc = params.get('z_column', 'PHOTOZ')
    if params.get('z') is not None:
        return float(params['z'])
    row = record.get('catalog_row') or {}
    for c in ([zc] if zc else []) + ['Z_SPEC', 'PHOTOZ', 'ZPHOT', 'Z']:
        try:
            v = float(row.get(c, ''))
            if v == v and v >= 0:
                return v
        except (TypeError, ValueError):
            pass
    return None


def process_script(code, script, records, params, task, model_label, need_filter_files):
    params = dict(params or {})
    if task != 'sed_fit':
        raise ValueError('%s adapter supports task sed_fit only' % code)
    results = [None] * len(records)
    per, bands_used = [], {}
    for r in records:
        ph, w = common.photometry(r, params.get('min_mag_err', 0.02), params.get('default_mag_err', 0.1), params.get('bands'))
        per.append({p['band']: p for p in ph})
        for p in ph:
            bands_used[p['band']] = p['entry']
    bands = sorted(bands_used, key=lambda b: bands_used[b]['lam'])
    good = []
    for i, r in enumerate(records):
        z = redshift_of(r, params)
        if z is None:
            results[i] = {'id': r['id'], 'error': 'no redshift (column %s / Z_SPEC / param z)' % params.get('z_column', 'PHOTOZ')}
        elif len(per[i]) < params.get('n_min_bands', 4):
            results[i] = {'id': r['id'], 'error': 'fewer than %d usable bands' % params.get('n_min_bands', 4)}
        else:
            good.append((i, z))
    if not good:
        return results, model_label % 'unused'
    wd = params.get('workdir') or tempfile.mkdtemp(prefix=code + '_')
    os.makedirs(wd, exist_ok=True)
    ff = {}
    if need_filter_files and not params.get('command'):
        fdir = params.get('filter_dir')
        if fdir:
            ff = {b: os.path.join(fdir, '%s.dat' % b) for b in bands if os.path.isfile(os.path.join(fdir, '%s.dat' % b))}
        else:
            fres = params.get('filters_res') or os.path.join(params.get('eazy_data') or os.environ.get('EAZYCODE', ''), 'filters', 'FILTER.RES.latest')
            if os.path.isfile(fres):
                ff = common.export_filter_files(bands, fres, os.path.join(wd, 'filters'))
        missing = [b for b in bands if b not in ff]
        if missing:
            bands = [b for b in bands if b in ff]
    F = np.full((len(good), len(bands)), -99.0); E = np.full_like(F, -99.0)
    for k, (i, _) in enumerate(good):
        for j, b in enumerate(bands):
            if b in per[i]:
                F[k, j], E[k, j] = per[i][b]['flux'], per[i][b]['fluxerr']
    req = dict(wd=wd, ids=[str(records[i]['id']) for i, _ in good], z=[z for _, z in good], bands=bands, flux=F.tolist(), err=E.tolist(), filter_files=ff, params=params)
    cmd = common.command_of(params) or [params.get('python') or sys.executable, os.path.join(HERE, 'native', script)]
    rc, so, se = common.run_command(list(cmd), cwd=wd, stdin_text=json.dumps(common.clean(req)), timeout=params.get('timeout', 3600))
    if rc != 0:
        raise RuntimeError('%s fit script failed (rc=%d): %s' % (code, rc, (se or so)[-600:]))
    out = json.loads(so[so.index('{'):])
    for (i, _), row in zip(good, out['rows']):
        results[i] = ({'id': records[i]['id'], **{k: row.get(k) for k in common.SED_COLS}} if 'error' not in row else {'id': records[i]['id'], 'error': row['error']})
    return results, model_label % out.get('version', '?')


def check(code, module, params=None):
    """Report whether the named library imports. Never raises.

    A failed import leaves engine_native false. The probe uses
    importlib so this file does not spell an import of the module name.
    """
    del params
    found = False
    version = None
    try:
        import importlib
        mod = importlib.import_module(module)
        found = True
        version = getattr(mod, "__version__", None)
    except Exception:
        found = False
        version = None
    if found:
        detail = "the installed package imported"
    else:
        detail = (
            "the external package is not called; "
            "the in-tree SED model is analytic"
        )
    return dict(code=code, engine_native=found, version=version, detail=detail)
