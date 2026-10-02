"""CIGALE adapter (SED fitting at the record's redshift, pdf_analysis).  Writes the CIGALE input files (observations table in mJy, `pcigale.ini`), runs `command`
(default ['pcigale', 'run']) in the work directory and parses `out/results.txt`.  CIGALE could not be installed where this adapter was written (its source is only
hosted on gitlab.lam.fr, unreachable from the build box): the file formats follow the CIGALE documentation and are exercised against the mock executable
`mocks/mock_cigale.py`, NOT against the real code - see docs/sed_codes.md.

params: command, workdir, z_column/z, filter_names {band: cigale filter name} (overrides sed_adapters.filters), sed_modules, module_params (dict module -> {param: values}),
cores, variables (list of CIGALE output variables), av_per_ebv, timeout, min_mag_err, default_mag_err, n_min_bands."""
import os
import re
import sys
import tempfile

import numpy as np

from . import common, filters, script_adapter

CODE = 'cigale'
TASKS = ('sed_fit',)
MODEL = 'CIGALE (%s)'
DEF_MODULES = ['sfhdelayed', 'bc03', 'nebular', 'dustatt_modified_starburst', 'redshifting']
DEF_PARAMS = {'sfhdelayed': {'tau_main': [100, 500, 2000, 8000], 'age_main': [200, 1000, 3000, 8000]},
              'bc03': {'imf': [0], 'metallicity': [0.004, 0.02]},
              'nebular': {'logU': [-2.0]},
              'dustatt_modified_starburst': {'E_BV_lines': [0.0, 0.1, 0.3, 0.6, 1.0]}}
DEF_VARS = ['stellar.m_star', 'sfh.sfr', 'sfh.age_main', 'dust.e_bv_lines']


def _fmt(v):
    return ', '.join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)


def write_inputs(wd, ids, zs, bands, F, E, names, params):
    cols = ['id', 'redshift'] + sum([[names[b], names[b] + '_err'] for b in bands], [])
    with open(os.path.join(wd, 'observations.txt'), 'w') as f:
        f.write('# ' + ' '.join(cols) + '\n')
        for i, oid in enumerate(ids):
            row = [str(oid), '%.5f' % zs[i]]
            for j in range(len(bands)):
                if F[i, j] > -90:
                    row += ['%.8g' % (F[i, j] * 1e-3), '%.8g' % (E[i, j] * 1e-3)]        # uJy -> mJy
                else:
                    row += ['nan', 'nan']
            f.write(' '.join(row) + '\n')
    mods = params.get('sed_modules') or DEF_MODULES
    mp = dict(DEF_PARAMS)
    mp.update(params.get('module_params') or {})
    vars_ = params.get('variables') or DEF_VARS
    lines = ['data_file = observations.txt', 'parameters_file = ', 'sed_modules = ' + ', '.join(mods), 'analysis_method = pdf_analysis', 'cores = %d' % params.get('cores', 4), '',
             '[sed_modules_params]']
    for m in mods:
        lines.append('  [[%s]]' % m)
        for k, v in (mp.get(m) or {}).items():
            lines.append('    %s = %s' % (k, _fmt(v)))
        if m == 'redshifting':
            lines.append('    redshift = ')
    lines += ['', '[analysis_params]', '  variables = ' + ', '.join(vars_), '  save_best_sed = False', '  save_chi2 = none', '  lim_flag = False', '  mock_flag = False',
              '  redshift_decimals = 2', '  blocks = 1']
    with open(os.path.join(wd, 'pcigale.ini'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def parse_results(path):
    with open(path) as f:
        head, out = None, {}
        for ln in f:
            if ln.startswith('#') or head is None:
                head = ln.lstrip('#').split()
                continue
            v = ln.split()
            if v:
                out[v[0]] = dict(zip(head, v))
    return out


def row_from(d, av_per_ebv=4.05):
    g = lambda k: float(d[k]) if k in d and d[k] not in ('nan', 'NaN', '') else float('nan')
    m, me = g('bayes.stellar.m_star'), g('bayes.stellar.m_star_err')
    age, agee = g('bayes.sfh.age_main'), g('bayes.sfh.age_main_err')
    ebv = g('bayes.dust.e_bv_lines')
    return dict(LOG_MASS=float(np.log10(m)) if m > 0 else None, LOG_MASS_ERR=float(me / (m * np.log(10))) if m > 0 and me == me else None,
                LOG_AGE=float(np.log10(age * 1e6)) if age > 0 else None, LOG_AGE_ERR=float(agee / (age * np.log(10))) if age > 0 and agee == agee else None,
                LOG_Z=float(np.log10(g('bayes.stellar.metallicity') / 0.02)) if g('bayes.stellar.metallicity') > 0 else None,
                AV=ebv * av_per_ebv if ebv == ebv else None, SFR=g('bayes.sfh.sfr'), SED_CHI2=g('best.reduced_chi_square'))


def process(records, params=None, task='sed_fit'):
    params = dict(params or {})
    if task != 'sed_fit':
        raise ValueError('cigale adapter supports task sed_fit only')
    results = [None] * len(records)
    per, used = [], {}
    for r in records:
        ph, _ = common.photometry(r, params.get('min_mag_err', 0.02), params.get('default_mag_err', 0.1), params.get('bands'))
        per.append({p['band']: p for p in ph})
        for p in ph:
            used[p['band']] = p['entry']
    bands = sorted(used, key=lambda b: used[b]['lam'])
    names = {b: (params.get('filter_names') or {}).get(b) or used[b]['cigale'] for b in bands}
    good, zs = [], []
    for i, r in enumerate(records):
        z = script_adapter.redshift_of(r, params)
        if z is None:
            results[i] = {'id': r['id'], 'error': 'no redshift'}
        elif len(per[i]) < params.get('n_min_bands', 4):
            results[i] = {'id': r['id'], 'error': 'fewer than %d usable bands' % params.get('n_min_bands', 4)}
        else:
            good.append(i); zs.append(z)
    if not good:
        return results, MODEL % 'unused'
    wd = params.get('workdir') or tempfile.mkdtemp(prefix='cigale_')
    os.makedirs(wd, exist_ok=True)
    F = np.full((len(good), len(bands)), -99.0); E = np.full_like(F, -99.0)
    for k, i in enumerate(good):
        for j, b in enumerate(bands):
            if b in per[i]:
                F[k, j], E[k, j] = per[i][b]['flux'], per[i][b]['fluxerr']
    write_inputs(wd, [records[i]['id'] for i in good], zs, bands, F, E, names, params)
    cmd = common.command_of(params) or ['pcigale', 'run']
    rc, so, se = common.run_command(list(cmd), cwd=wd, timeout=params.get('timeout', 3600))
    if rc != 0:
        raise RuntimeError('CIGALE command failed (rc=%d): %s' % (rc, (se or so)[-500:]))
    res = parse_results(os.path.join(wd, 'out', 'results.txt'))
    for i in good:
        d = res.get(str(records[i]['id']))
        results[i] = {'id': records[i]['id'], **row_from(d, params.get('av_per_ebv', 4.05))} if d else {'id': records[i]['id'], 'error': 'no row in results.txt'}
    return results, MODEL % ('external ' + os.path.basename(cmd[0]))


def check(params=None):
    import shutil
    cmd = common.command_of(params or {}) or ['pcigale', 'run']
    return dict(code=CODE, engine_native=bool(shutil.which(cmd[0])), version=None, detail=None if shutil.which(cmd[0]) else '%s not found on PATH' % cmd[0])


def run(records, params, context=None):
    return process(records, params, (context or {}).get('task', 'sed_fit'))[0]


def main():
    return common.adapter_main(process)


if __name__ == '__main__':
    sys.exit(main())
