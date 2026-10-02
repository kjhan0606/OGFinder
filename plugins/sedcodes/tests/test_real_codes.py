"""Smoke / accuracy tests against the REAL codes when they are installed (skipped otherwise): eazy-py (pip install eazy + eazy-photoz data) and Bagpipes.
Interpreters: $OGF_EAZY_PYTHON (default /workspace/eazy_venv/bin/python), $EAZYCODE (default /workspace/eazy_data/eazy-photoz), $OGF_BAGPIPES_PYTHON (default /workspace/sed_venv/bin/python)."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from sed_adapters import common, eazy_adapter, bagpipes_adapter, prospector_adapter

HERE = os.path.dirname(os.path.abspath(__file__))
EZ_PY = os.environ.get('OGF_EAZY_PYTHON', '/workspace/eazy_venv/bin/python')
EZ_DATA = os.environ.get('EAZYCODE', '/workspace/eazy_data/eazy-photoz')
BP_PY = os.environ.get('OGF_BAGPIPES_PYTHON', '/workspace/sed_venv/bin/python')
FRES = os.path.join(EZ_DATA, 'filters', 'FILTER.RES.latest')


def nmad(dz):
    return 1.4826 * np.median(np.abs(dz - np.median(dz)))


@pytest.mark.skipif(not (os.path.exists(EZ_PY) and os.path.exists(FRES)), reason='eazy-py environment not available')
def test_real_eazy_py(tmp_path):
    p = subprocess.run([EZ_PY, os.path.join(HERE, 'make_eazy_truth.py'), EZ_DATA, '60'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-400:]
    objs = json.loads(p.stdout[p.stdout.index('['):])
    recs = [dict(id=o['id'], mags=o['mags'], mag_errs=o['mag_errs'], catalog_row={}) for o in objs]
    res, model = eazy_adapter.process(recs, dict(engine='native', python=EZ_PY, eazy_data=EZ_DATA, workdir=str(tmp_path / 'w'), z_max=4.0, z_step=0.01, min_mag_err=0.0), 'photoz')
    assert all('error' not in r for r in res), [r for r in res if 'error' in r][:2]
    z = np.array([r['PHOTOZ'] for r in res]); zt = np.array([o['z'] for o in objs])
    dz = (z - zt) / (1 + zt)
    lo = np.array([r['PHOTOZ_P16'] for r in res]); hi = np.array([r['PHOTOZ_P84'] for r in res])
    cov = float(np.mean((zt >= lo) & (zt <= hi)))
    print('REAL %s: sigma_NMAD %.4f outliers %.3f, 68%% interval coverage %.2f, median chi2_red %.2f' % (model, nmad(dz), np.mean(np.abs(dz) > 0.15), cov, np.median([r['PHOTOZ_CHI2'] for r in res])))
    assert nmad(dz) < 0.05 and np.mean(np.abs(dz) > 0.15) < 0.15 and 'eazy-py' in model


def _bagpipes_ok():
    if not os.path.exists(BP_PY) or not os.path.exists(FRES):
        return False
    return subprocess.run([BP_PY, '-c', 'import bagpipes'], capture_output=True).returncode == 0


@pytest.mark.skipif(not _bagpipes_ok(), reason='bagpipes environment not available')
def test_real_bagpipes(tmp_path):
    bands = ['F435W', 'F606W', 'F775W', 'F850LP', 'F105W', 'F125W', 'F160W']
    fdir = str(tmp_path / 'filters')
    common.export_filter_files(bands, FRES, fdir)
    env = dict(os.environ, MPLBACKEND='Agg')
    p = subprocess.run([BP_PY, os.path.join(HERE, 'make_bagpipes_truth.py'), fdir, '3'], capture_output=True, text=True, env=env, cwd=str(tmp_path))
    assert p.returncode == 0, p.stderr[-500:]
    objs = json.loads(p.stdout[p.stdout.index('['):])
    recs = [dict(id=o['id'], mags=o['mags'], mag_errs=o['mag_errs'], catalog_row={'Z_SPEC': '%.5f' % o['z']}) for o in objs]
    res, model = bagpipes_adapter.process(recs, dict(python=BP_PY, filter_dir=fdir, workdir=str(tmp_path / 'w'), n_live=300, min_mag_err=0.0, z_column='Z_SPEC'), 'sed_fit')
    assert all('error' not in r for r in res), res
    dm = np.array([r['LOG_MASS'] - o['logm'] for r, o in zip(res, objs)])
    pull = np.array([(r['LOG_MASS'] - o['logm']) / max(r['LOG_MASS_ERR'], 0.03) for r, o in zip(res, objs)])
    print('REAL %s: logM - truth = %s (pull %s) ; Av - truth = %s' % (model, np.round(dm, 3), np.round(pull, 2), np.round([r['AV'] - o['av'] for r, o in zip(res, objs)], 2)))
    assert np.all(np.isfinite(dm)) and np.max(np.abs(dm)) < 0.5 and 'Bagpipes' in model


def _prospector_ok():
    py = os.environ.get('OGF_PROSPECTOR_PYTHON', BP_PY)
    if not os.path.exists(py):
        return False
    code = 'import prospect, fsps; fsps.StellarPopulation(zcontinuous=1)'
    return subprocess.run([py, '-c', code], capture_output=True, env=dict(os.environ, SPS_HOME=os.environ.get('SPS_HOME', ''))).returncode == 0


@pytest.mark.skipif(not _prospector_ok(), reason='working python-fsps (SPS_HOME) + prospector not available')
def test_real_prospector(tmp_path):
    pytest.skip('not exercised: no verified FSPS installation in the development environment')
