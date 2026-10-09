"""Photo-z adapter.

engine=native calls ogfmeas.photoz.template_chi2 when ``params['templates']`` is a list of numeric ``[rest_wavelength, rest_flux]`` pairs. Band photometry is then sampled at registry pivot wavelengths (a coarse point sample, not a filter integral). Template-file paths are not read, and the external photo-z package is not called. Without numeric templates the native engine stops.
engine=external: the classic file interface - writes catalog / zphot.param / zphot.translate, runs `command <param file>`, parses `<MAIN_OUTPUT_FILE>.zout`.
engine=package runs the user-installed eazy-py package in a separate interpreter, only when that name is requested. A missing package is an error and is not labeled as eazy.

Entry points: `process(records, params, task)` (pure, JSON in/out), `run(records, params, context)` (ai_bridge python_callable), `python -m sed_adapters.eazy_adapter`
(ai_bridge local_command, JSON on stdin/stdout).

params (all optional): engine, python, command (argv list), templates (numeric [rest_wavelength, rest_flux] pairs; native engine only), eazy_data (external engine; dir with filters/ templates/, default $EAZYCODE), filters_res, filters_info, templates_file,
z_min z_max z_step, min_mag_err, default_mag_err, point_estimate (peak|ml), workdir, timeout, prior (bool), n_min_colors.
"""
import json
import math
import os
import re
import sys
import tempfile

import numpy as np

from . import common, filters

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = 'eazy'
TASKS = ('photoz',)
MODEL = 'EAZY (%s)'


def filter_names(info_path):
    """Names from a FILTER.RES.info file ('  1 name lambda_c= ...') in order."""
    out = {}
    with open(info_path) as f:
        for ln in f:
            m = re.match(r'\s*(\d+)\s+(\S+)', ln)
            if m:
                out[int(m.group(1))] = m.group(2)
    return [out[k] for k in sorted(out)]


def default_paths(params):
    root = params.get('eazy_data') or os.environ.get('EAZYCODE') or ''
    fres = params.get('filters_res') or (os.path.join(root, 'filters', 'FILTER.RES.latest') if root else '')
    info = params.get('filters_info') or (fres + '.info' if fres else '')
    tpl = params.get('templates_file') or (os.path.join(root, 'templates', 'eazy_v1.3.spectra.param') if root else '')
    return root, fres, info, tpl


def collect(records, params):
    """-> (bands (registry names present in any record), matrices flux/err (nrec, nband; -99 where missing), per-record warnings)."""
    per, warns, used = [], [], {}
    for r in records:
        ph, w = common.photometry(r, params.get('min_mag_err', 0.02), params.get('default_mag_err', 0.1), params.get('bands'))
        per.append({p['band']: p for p in ph})
        warns.append(w)
        for p in ph:
            used[p['band']] = p['entry']
    bands = sorted(used, key=lambda b: used[b]['lam'])
    F = np.full((len(records), len(bands)), -99.0)
    E = np.full_like(F, -99.0)
    for i, d in enumerate(per):
        for j, b in enumerate(bands):
            if b in d:
                F[i, j], E[i, j] = d[b]['flux'], d[b]['fluxerr']
    return bands, F, E, warns


def write_eazy_inputs(wd, records, bands, F, E, idx, params, fres, tpl):
    """Classic EAZY files; returns the param-file path.  Column names F_<band>/E_<band>, mapped to filter numbers in zphot.translate."""
    cols = ['id'] + sum([['F_' + b, 'E_' + b] for b in bands], []) + ['z_spec']
    with open(os.path.join(wd, 'catalog.cat'), 'w') as f:
        f.write('# ' + ' '.join(cols) + '\n')
        for i, r in enumerate(records):
            vals = [str(r['id'])] + sum([['%.8g' % F[i, j], '%.8g' % E[i, j]] for j in range(len(bands))], []) + ['-1']
            f.write(' '.join(vals) + '\n')
    with open(os.path.join(wd, 'zphot.translate'), 'w') as f:
        for b, k in zip(bands, idx):
            f.write('F_%s F%d\nE_%s E%d\n' % (b, k, b, k))
    lines = ['CATALOG_FILE catalog.cat', 'CATALOG_FORMAT ascii.commented_header', 'MAIN_OUTPUT_FILE photz', 'FILTERS_RES %s' % fres, 'TEMPLATES_FILE %s' % tpl,
             'Z_MIN %g' % params.get('z_min', 0.01), 'Z_MAX %g' % params.get('z_max', 6.0), 'Z_STEP %g' % params.get('z_step', 0.01), 'N_MIN_COLORS %d' % params.get('n_min_colors', 3),
             'APPLY_PRIOR %s' % ('y' if params.get('prior') else 'n'), 'PRIOR_ABZP 23.9', 'FIX_ZSPEC n', 'CAT_HAS_EXTCORR n', 'TRANSLATE_FILE zphot.translate', 'WAVELENGTH_FILE', 'TEMP_ERR_FILE']
    p = os.path.join(wd, 'zphot.param')
    with open(p, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return p


def parse_zout(path):
    """EAZY .zout -> {id: dict}.  Header '# id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68 ...'."""
    out, head = {}, None
    with open(path) as f:
        for ln in f:
            if ln.startswith('#'):
                t = ln[1:].split()
                if 'id' in t and head is None:
                    head = t
                continue
            if not ln.strip() or head is None:
                continue
            v = ln.split()
            out[v[0]] = dict(zip(head, v))
    return out


def from_zout(d, point='peak'):
    g = lambda k: float(d[k]) if k in d else float('nan')
    z = g('z_peak') if point == 'peak' and 'z_peak' in d else (g('z_m1') if point == 'ml' else g('z_p'))
    p16, p84 = g('l68'), g('u68')
    chi = g('chi_p') if 'chi_p' in d else g('chi_a')
    nf = g('nfilt')
    return common.photoz_row(z, p16, g('z_m1'), p84, chi / max(nf - 1, 1) if nf == nf else chi)


def process(records, params=None, task='photoz'):
    params = dict(params or {})
    if task != 'photoz':
        raise ValueError('eazy adapter supports task photoz only')
    if not records:
        label = MODEL % 'no records' if params.get('engine', 'native') == 'external' else 'ogfmeas.photoz'
        return [], label
    bands, F, E, warns = collect(records, params)
    nmin = params.get('n_min_colors', 3)
    engine = params.get('engine', 'native')
    _root, fres, info, tpl = default_paths(params)
    results = [None] * len(records)
    good = [i for i in range(len(records)) if int(np.sum(F[i] > -90)) >= nmin]
    for i in range(len(records)):
        if i not in good:
            results[i] = {'id': records[i]['id'], 'error': 'fewer than %d usable bands' % nmin}
    if not good:
        label = MODEL % 'unused' if engine == 'external' else 'ogfmeas.photoz'
        return results, label
    sub = [records[i] for i in good]
    wd = params.get('workdir') or tempfile.mkdtemp(prefix='eazy_')
    os.makedirs(wd, exist_ok=True)
    timeout = params.get('timeout', 900)
    if engine == 'external':
        cmd = common.command_of(params)
        if not cmd:
            raise ValueError('engine external needs params.command (argv list: the EAZY executable or a compatible program)')
        if not os.path.isfile(info):
            raise ValueError('FILTER.RES.info not found (%s): set filters_info / eazy_data' % info)
        names = filter_names(info)
        idx, keep = [], []
        for j, b in enumerate(bands):
            k = filters.eazy_filter_index(filters.B[b], names)
            if k is None:
                warns[0].append('band %s: no matching filter in %s' % (b, info))
            else:
                idx.append(k); keep.append(j)
        bands2 = [bands[j] for j in keep]
        zp = write_eazy_inputs(wd, sub, bands2, F[good][:, keep], E[good][:, keep], idx, params, fres, tpl)
        rc, so, se = common.run_command(list(cmd) + [os.path.basename(zp)], cwd=wd, timeout=timeout)
        if rc != 0:
            raise RuntimeError('EAZY command failed (rc=%d): %s' % (rc, (se or so)[-400:]))
        zo = parse_zout(os.path.join(wd, 'photz.zout'))
        ver = 'external %s' % os.path.basename(cmd[0])
        for i in good:
            d = zo.get(str(records[i]['id']))
            results[i] = ({'id': records[i]['id'], **from_zout(d, params.get('point_estimate', 'peak'))} if d else {'id': records[i]['id'], 'error': 'no row in .zout'})
    elif engine == 'package':
        py = params.get('python') or sys.executable
        req = dict(wd=wd, bands=bands, flux=F[good].tolist(), err=E[good].tolist(), ids=[str(records[i]['id']) for i in good],
                   filters_res=fres, templates_file=tpl, eazy_data=_root or None,
                   z_min=params.get('z_min', 0.01), z_max=params.get('z_max', 6.0), z_step=params.get('z_step', 0.01),
                   prior=bool(params.get('prior')))
        rc, so, se = common.run_command([py, os.path.join(HERE, 'native', 'eazy_package.py')], cwd=wd, stdin_text=json.dumps(common.clean(req)), timeout=timeout)
        if rc != 0:
            raise RuntimeError('eazy package helper failed (rc=%d): %s' % (rc, (se or so)[-600:]))
        out = json.loads(so[so.index('{'):])
        ver = out.get('version') or 'eazy'
        for i, row in zip(good, out['rows']):
            results[i] = ({'id': records[i]['id'], **common.photoz_row(row['z'], row['p16'], row['p50'], row['p84'], row['chi2'], row['sigma'])} if row.get('z') is not None
                          else {'id': records[i]['id'], 'error': row.get('error', 'fit failed')})
    elif engine == 'native':
        py = params.get('python') or sys.executable
        req = dict(wd=wd, bands=bands, flux=F[good].tolist(), err=E[good].tolist(), ids=[str(records[i]['id']) for i in good],
                   z_min=params.get('z_min', 0.01), z_max=params.get('z_max', 6.0), z_step=params.get('z_step', 0.01))
        if params.get('templates') is not None:
            req['templates'] = params['templates']
        if params.get('wavelength') is not None:
            req['wavelength'] = params['wavelength']
        if params.get('redshifts') is not None:
            req['redshifts'] = params['redshifts']
        rc, so, se = common.run_command([py, os.path.join(HERE, 'native', 'eazy_native.py')], cwd=wd, stdin_text=json.dumps(common.clean(req)), timeout=timeout)
        if rc != 0:
            raise RuntimeError('photo-z helper failed (rc=%d): %s' % (rc, (se or so)[-600:]))
        out = json.loads(so[so.index('{'):])
        ver = out.get('version') or 'ogfmeas.photoz'
        for i, row in zip(good, out['rows']):
            results[i] = ({'id': records[i]['id'], **common.photoz_row(row['z'], row['p16'], row['p50'], row['p84'], row['chi2'], row['sigma'])} if row.get('z') is not None
                          else {'id': records[i]['id'], 'error': row.get('error', 'fit failed')})
    else:
        raise ValueError('unknown eazy engine %r' % engine)
    for r, w in zip(results, warns):
        if w and 'error' not in r:
            r['warning'] = '; '.join(w[:3])
    if engine == 'external':
        model = MODEL % ver
    elif engine == 'package':
        model = 'eazy (%s)' % ver
    else:
        model = 'ogfmeas.photoz' if ver == 'ogfmeas.photoz' else ver
    return results, model


def check(params=None):
    """Availability report (never raises). The native engine does not call an external package."""
    params = params or {}
    _root, fres, _info, tpl = default_paths(params)
    return dict(code=CODE, engine_native=False, version='ogfmeas.photoz', filters_res=os.path.isfile(fres), templates=os.path.isfile(tpl),
                detail='the external package is not called; numeric rest-frame templates are fit by ogfmeas.photoz')


def run(records, params, context=None):
    res, _ = process(records, params, (context or {}).get('task', 'photoz'))
    return res


def main():
    return common.adapter_main(process)


if __name__ == '__main__':
    sys.exit(main())
