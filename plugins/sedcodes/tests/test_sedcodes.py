"""SED-code adapters: unit conversion, band/filter mapping, native file formats, parsing, entry points (python_callable, local_command, ai_bridge profile, plugin CLI).
The external codes are replaced by the mock executables in sed_adapters/mocks (TOY templates): this validates the plumbing against a known truth, NOT any code's accuracy."""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from sed_adapters import common, filters, toy, eazy_adapter, cigale_adapter, bagpipes_adapter, prospector_adapter
from sed_adapters import mocks as _mocks

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
MOCK = os.path.join(ROOT, 'sed_adapters', 'mocks')
PY = sys.executable
BANDS = ['F435W', 'F606W', 'F775W', 'F850LP', 'F105W', 'F125W', 'F160W']
INFO = os.path.join(MOCK, 'FILTER.RES.toy.info')
FRES = os.path.join(MOCK, 'FILTER.RES.toy')


def population(n=80, seed=5, zmax=3.2, noise=0.05, avs=(0, 0.25, 0.5, 1.0)):
    rng = np.random.default_rng(seed)
    recs, truth = [], []
    for k in range(n):
        i = int(rng.integers(0, toy.NT)); z = float(rng.uniform(0.2, zmax)); mag = float(rng.uniform(21.5, 24.5)); av = float(rng.choice(avs))
        fl, err, sc = toy.simulate(i, z, mag, BANDS, 'F160W', av=av, rng=rng, noise_frac=noise)
        mags = {b: (float(common.ujy_to_mag(f)) if f > 0 else None) for b, f in zip(BANDS, fl)}
        me = {b: float(1.0857 * e / f) if f > 0 else None for b, f, e in zip(BANDS, fl, err)}
        recs.append(dict(id=str(k + 1), mags=mags, mag_errs=me, catalog_row={'Z_SPEC': '%.4f' % z}))
        truth.append(dict(z=z, i=i, av=av, scale=sc, logm=toy.log_mass(i, sc, z)))
    return recs, truth


def nmad(dz):
    return 1.4826 * np.median(np.abs(dz - np.median(dz)))


def test_photometry_conversion():
    assert abs(float(common.mag_to_ujy(23.9)) - 1.0) < 1e-12 and abs(float(common.mag_to_ujy(22.9)) - 10 ** 0.4) < 1e-9
    assert abs(float(common.ujy_to_mag(1.0)) - 23.9) < 1e-12
    f = float(common.mag_to_ujy(25.0)); e = float(common.magerr_to_ujy_err(25.0, 0.1))
    assert abs(e / f - 0.1 * math.log(10) / 2.5) < 1e-12
    rec = {'mags': {'F160W': 24.0, 'F105W': None, 'AUTO': 22.0, 'XYZ': 23.0, 'g': 23.0}, 'mag_errs': {'F160W': 0.05}}
    ph, w = common.photometry(rec, min_mag_err=0.02)
    assert [p['band'] for p in ph] == ['F160W'] and abs(ph[0]['magerr'] - math.hypot(0.05, 0.02)) < 1e-12
    assert any('XYZ' in x for x in w) and not any('AUTO' in x for x in w)
    ph2, _ = common.photometry({'mags': {'F160W': 24.0}, 'mag_errs': {}}, default_mag_err=0.2, min_mag_err=0.0)
    assert abs(ph2[0]['magerr'] - 0.2) < 1e-12


def test_pdf_summary():
    z = np.linspace(0, 4, 4001)
    s = common.pdf_summary(z, np.exp(-0.5 * ((z - 1.5) / 0.1) ** 2))
    assert abs(s['p50'] - 1.5) < 1e-3 and abs(s['p16'] - 1.4) < 3e-3 and abs(s['p84'] - 1.6) < 3e-3 and abs(s['z_sigma'] - 0.1) < 2e-3


def test_filter_registry_and_eazy_index():
    names = eazy_adapter.filter_names(INFO)
    seen = {}
    for b, e in filters.B.items():
        assert e['lam'] > 0 and e['dl'] > 0 and e['sedpy'] and e['cigale'] and e['bagpipes'] and e['eazy']
        k = filters.eazy_filter_index(e, names)
        assert k is not None and k not in seen, (b, k)
        seen[k] = b
    assert filters.lookup('acs_f606w')['name'] == 'F606W' and filters.lookup('f160w')['name'] == 'F160W' and filters.lookup('nonsense') is None
    real = '/workspace/eazy_data/eazy-photoz/filters/FILTER.RES.latest'
    if os.path.exists(real):
        flt = common.read_filter_res(real)
        rn = [f[0] for f in flt]
        miss = [b for b, e in filters.B.items() if filters.eazy_filter_index(e, rn) is None]
        assert not miss, miss
        assert filters.eazy_filter_index(filters.B['F160W'], rn) == 205 and filters.eazy_filter_index(filters.B['F606W'], rn) == 236 and filters.eazy_filter_index(filters.B['F105W'], rn) == 202


def test_toy_lyman_break_and_flux_positive():
    assert toy.band_flux(3, 0.5, 'F160W') > 0
    f435_hi = toy.band_flux(3, 4.5, 'F435W'); f435_lo = toy.band_flux(3, 0.5, 'F435W')
    assert f435_hi < 0.02 * f435_lo


@pytest.fixture(scope='module')
def eazy_run(tmp_path_factory):
    recs, truth = population(avs=(0.0,))
    wd = str(tmp_path_factory.mktemp('eazy'))
    res, model = eazy_adapter.process(recs, dict(engine='external', command=[PY, os.path.join(MOCK, 'mock_eazy.py')], filters_res=FRES, templates_file='toy', workdir=wd, z_max=4.0, z_step=0.02), 'photoz')
    return recs, truth, res, model, wd


def test_eazy_external_recovers_truth_and_files(eazy_run):
    recs, truth, res, model, wd = eazy_run
    z = np.array([r['PHOTOZ'] for r in res]); zt = np.array([t['z'] for t in truth])
    dz = (z - zt) / (1 + zt)
    s, out = nmad(dz), float(np.mean(np.abs(dz) > 0.15))
    lo = np.array([r['PHOTOZ_P16'] for r in res]); hi = np.array([r['PHOTOZ_P84'] for r in res])
    print('EAZY(mock) sigma_NMAD %.4f outliers %.3f ; 68%% interval contains truth %.2f' % (s, out, np.mean((zt >= lo) & (zt <= hi))))
    assert s < 0.03 and out < 0.12 and all('error' not in r for r in res) and 'EAZY' in model
    # the adapter chain equals the direct TOY fit (same grid)
    zg = np.arange(0.01, 4.0 + 1e-9, 0.02)
    r0 = recs[3]
    ph, _ = common.photometry(r0)
    chi, _ = toy.fit_photoz(np.array([p['flux'] for p in ph]), np.array([p['fluxerr'] for p in ph]), [p['band'] for p in ph], zg)
    assert abs(res[3]['PHOTOZ'] - toy.fit_photoz.__globals__['np'].asarray(zg)[int(np.argmin(chi))]) <= 0.021
    # file formats written for EAZY: fluxes are micro-Jansky (AB 23.9), translate maps to FILTER.RES numbers
    cat = open(os.path.join(wd, 'catalog.cat')).read().split('\n')
    hdr = cat[0][1:].split(); row = cat[1].split(); d = dict(zip(hdr, row))
    assert hdr[0] == 'id' and 'F_F160W' in hdr and 'E_F160W' in hdr and hdr[-1] == 'z_spec'
    assert abs(float(d['F_F160W']) / float(common.mag_to_ujy(recs[0]['mags']['F160W'])) - 1) < 1e-6
    names = eazy_adapter.filter_names(INFO)
    tr = dict(l.split() for l in open(os.path.join(wd, 'zphot.translate')) if l.strip())
    assert 'f160w' in names[int(tr['F_F160W'][1:]) - 1] and 'f606w' in names[int(tr['F_F606W'][1:]) - 1]
    p = dict(l.split(None, 1) for l in open(os.path.join(wd, 'zphot.param')) if l.strip() and len(l.split(None, 1)) == 2)
    assert p['CATALOG_FILE'].strip() == 'catalog.cat' and p['Z_STEP'].strip() == '0.02' and p['APPLY_PRIOR'].strip() == 'n'


def test_eazy_missing_data_and_errors():
    recs, _ = population(4)
    recs[0]['mags'] = {'F160W': 24.0}                      # one band only
    res, _ = eazy_adapter.process(recs, dict(engine='external', command=[PY, os.path.join(MOCK, 'mock_eazy.py')], filters_res=FRES, z_step=0.05), 'photoz')
    assert 'error' in res[0] and 'PHOTOZ' in res[1]
    with pytest.raises(ValueError):
        eazy_adapter.process(recs, dict(engine='external', filters_res=FRES), 'photoz')
    with pytest.raises(RuntimeError):
        eazy_adapter.process(recs, dict(engine='external', command=[PY, '-c', 'import sys; sys.exit(3)'], filters_res=FRES), 'photoz')
    with pytest.raises(ValueError):
        eazy_adapter.process(recs, dict(engine='external', command=['x']), 'sed_fit')


def _with_z(recs, truth):
    for r, t in zip(recs, truth):
        r['catalog_row']['PHOTOZ'] = '%.4f' % t['z']
    return recs


def sed_stats(res, truth):
    dm = np.array([r['LOG_MASS'] - t['logm'] for r, t in zip(res, truth)])
    da = np.array([r['AV'] - t['av'] for r, t in zip(res, truth)])
    return float(np.sqrt(np.mean(dm ** 2))), float(np.median(dm)), float(np.sqrt(np.mean(da ** 2)))


def test_cigale_external_files_and_recovery(tmp_path):
    recs, truth = population(40, seed=8, zmax=2.5, noise=0.04)
    recs = _with_z(recs, truth)
    wd = str(tmp_path / 'cg')
    res, model = cigale_adapter.process(recs, dict(command=[PY, os.path.join(MOCK, 'mock_cigale.py')], workdir=wd, z_column='PHOTOZ'), 'sed_fit')
    rms, med, rav = sed_stats(res, truth)
    print('CIGALE(mock) logM rms %.3f median %.3f ; Av rms %.3f' % (rms, med, rav))
    assert rms < 0.35 and abs(med) < 0.15 and rav < 0.6 and all('error' not in r for r in res)
    ini = open(os.path.join(wd, 'pcigale.ini')).read()
    assert 'data_file = observations.txt' in ini and 'analysis_method = pdf_analysis' in ini and 'variables = stellar.m_star' in ini and '[[redshifting]]' in ini
    lines = open(os.path.join(wd, 'observations.txt')).read().split('\n')
    hdr = lines[0][1:].split(); d = dict(zip(hdr, lines[1].split()))
    assert hdr[:2] == ['id', 'redshift'] and 'hst.wfc3.F160W' in hdr and 'hst.wfc3.F160W_err' in hdr
    assert abs(float(d['hst.wfc3.F160W']) * 1e3 / float(common.mag_to_ujy(recs[0]['mags']['F160W'])) - 1) < 1e-5          # mJy
    assert abs(float(d['redshift']) - truth[0]['z']) < 1e-4
    # no redshift -> row error, user filter names respected
    r2 = [dict(recs[0], catalog_row={})]
    out, _ = cigale_adapter.process(r2, dict(command=[PY, os.path.join(MOCK, 'mock_cigale.py')], workdir=str(tmp_path / 'cg2')), 'sed_fit')
    assert 'error' in out[0] and 'redshift' in out[0]['error']


@pytest.mark.parametrize('mod', [bagpipes_adapter, prospector_adapter])
def test_script_protocol_with_mock(mod, tmp_path):
    recs, truth = population(30, seed=9, zmax=2.2, noise=0.04)
    recs = _with_z(recs, truth)
    res, model = mod.process(recs, dict(command=[PY, os.path.join(MOCK, 'mock_fit.py')], workdir=str(tmp_path / 'w'), z_column='PHOTOZ', filters_res='/nonexistent'), 'sed_fit') \
        if mod is prospector_adapter else \
        mod.process(recs, dict(command=[PY, os.path.join(MOCK, 'mock_fit.py')], workdir=str(tmp_path / 'w'), z_column='PHOTOZ', filters_res=FRES), 'sed_fit')
    rms, med, rav = sed_stats(res, truth)
    print('%s(mock) logM rms %.3f median %.3f Av rms %.3f' % (mod.CODE, rms, med, rav))
    assert rms < 0.35 and abs(med) < 0.15 and all('error' not in r for r in res) and 'mock-fit' in model


def test_entry_points_python_callable_and_local_command(tmp_path):
    recs, truth = population(6, seed=3, avs=(0.0,))
    params = dict(engine='external', command=[PY, os.path.join(MOCK, 'mock_eazy.py')], filters_res=FRES, z_step=0.05, workdir=str(tmp_path / 'a'))
    r1 = eazy_adapter.run(recs, params, {'task': 'photoz'})
    req = json.dumps({'contract': 'ogf-ai-bridge/1', 'task': 'photoz', 'service': 'x', 'params': dict(params, workdir=str(tmp_path / 'b')), 'records': recs})
    p = subprocess.run([PY, '-m', 'sed_adapters.eazy_adapter'], input=req, capture_output=True, text=True, cwd=ROOT, env=dict(os.environ, PYTHONPATH=ROOT))
    assert p.returncode == 0, p.stderr[-300:]
    out = json.loads(p.stdout)
    assert 'EAZY' in out['model'] and [x['id'] for x in out['results']] == [x['id'] for x in r1]
    assert max(abs(a['PHOTOZ'] - b['PHOTOZ']) for a, b in zip(out['results'], r1)) < 1e-9
    # python_callable signature: module:function(records, params, context)
    import importlib
    f = getattr(importlib.import_module('sed_adapters.bagpipes_adapter'), 'run')
    assert callable(f)


def write_cat(path, recs, truth):
    cols = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'Z_SPEC'] + sum([['MAG_' + b, 'MAGERR_' + b] for b in BANDS], [])
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for k, (r, t) in enumerate(zip(recs, truth)):
            row = [r['id'], '%d' % (10 + k), '%d' % (20 + k), '%.4f' % t['z']]
            for b in BANDS:
                row += ['%.4f' % r['mags'][b] if r['mags'][b] is not None else '99', '%.4f' % r['mag_errs'][b] if r['mag_errs'][b] is not None else '99']
            fh.write('\t'.join(row) + '\n')


def test_plugin_cli_and_ai_bridge_profile(tmp_path):
    recs, truth = population(25, seed=12, avs=(0.0,))
    cat = str(tmp_path / 'cat.tsv'); write_cat(cat, recs, truth)
    sc = os.path.join(ROOT, 'plugins', 'sedcodes', 'sedcodes.py')
    cmd = [PY, sc, '--task', 'photoz', '--code', 'eazy', '--catalog', cat, '--work', str(tmp_path / 'w'), '--engine', 'external', '--command', '%s %s' % (PY, os.path.join(MOCK, 'mock_eazy.py')),
           '--filters-res', FRES, '--z-step', '0.02', '--z-max', '4.0', '--meta-out', str(tmp_path / 'catalog_meta.json')]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-400:]
    lines = [l.split('\t') for l in p.stdout.strip().split('\n')]
    assert lines[0][:3] == ['NUMBER', 'EZ_Z', 'EZ_ZERR'] and len(lines) == 26
    z = np.array([float(l[1]) for l in lines[1:]]); zt = np.array([t['z'] for t in truth])
    assert nmad((z - zt) / (1 + zt)) < 0.03
    assert json.load(open(tmp_path / 'catalog_meta.json'))['sedcodes']['code'] == 'eazy'
    # sed fit through the CLI with the mock fit script
    cmd2 = [PY, sc, '--task', 'sed_fit', '--code', 'bagpipes', '--catalog', cat, '--work', str(tmp_path / 'w2'), '--command', '%s %s' % (PY, os.path.join(MOCK, 'mock_fit.py')),
            '--filters-res', FRES, '--z-column', 'Z_SPEC']
    p2 = subprocess.run(cmd2, capture_output=True, text=True)
    assert p2.returncode == 0, p2.stderr[-400:]
    l2 = [l.split('\t') for l in p2.stdout.strip().split('\n')]
    assert l2[0][:3] == ['NUMBER', 'SC_LOGM', 'SC_LOGM_ERR'] and all(x[1] != '' for x in l2[1:])
    # unsupported task for the code
    p3 = subprocess.run([PY, sc, '--task', 'sed_fit', '--code', 'eazy', '--catalog', cat, '--work', str(tmp_path / 'w3')], capture_output=True, text=True)
    assert p3.returncode == 2
    # ai_bridge profile file written by the plugin, run through ds9_ai_bridge.py (local_command)
    prof = str(tmp_path / 'ai.json')
    assert subprocess.run([PY, sc, '--write-profiles', prof], capture_output=True).returncode == 0
    d = json.load(open(prof))
    for s in d['services']:
        s['enabled'] = True
        s['env_passthrough'] = ['PYTHONPATH']
    json.dump(d, open(prof, 'w'))
    out_tsv = str(tmp_path / 'ai_out.tsv')
    br = subprocess.run([PY, os.path.join(ROOT, 'ds9', 'library', 'ds9_ai_bridge.py'), '--mode', 'run', '--services-file', prof, '--service', 'sedcodes_eazy', '--task', 'photoz', '--catalog', cat,
                         '--output', out_tsv, '--no-cache', '--param', 'engine=external', '--param', 'command=%s %s' % (PY, os.path.join(MOCK, 'mock_eazy.py')), '--param', 'filters_res=' + FRES,
                         '--param', 'z_step=0.05', '--param', 'workdir=' + str(tmp_path / 'brw')], capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=ROOT), cwd=ROOT)
    assert br.returncode == 0, (br.stdout + br.stderr)[-600:]
    rows = [l.split('\t') for l in open(out_tsv).read().strip().split('\n')]
    ci = rows[0].index('PHOTOZ')
    zb = np.array([float(r[ci]) for r in rows[1:]])
    assert len(zb) == 25 and nmad((zb - zt) / (1 + zt)) < 0.04


def test_check_reports():
    r = eazy_adapter.check({'python': PY})
    assert r['code'] == 'eazy' and 'engine_native' in r
