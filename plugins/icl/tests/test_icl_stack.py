"""Menu command for the scaled BCG+diffuse stack."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit.imageio import load_image, save_fits  # noqa: E402
from ogfkit.wlcluster import Cosmology  # noqa: E402
from icl.wgl import theta_delta_arcsec  # noqa: E402

CLI = os.path.join(ROOT, 'plugins', 'icl', 'icl_stack.py')
PLUGIN = os.path.join(ROOT, 'plugins', 'icl', 'plugin.json')
FLAGS = ['--u-max', '0.5', '--n-side', '11', '--delta', '200', '--h0', '70', '--omega-m', '0.3']


def _run(args):
    return subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)


def _frame(path, counts, variance=None):
    image = np.full((41, 41), float(counts))
    save_fits(path, image)
    record = {
        'image': os.path.basename(path),
        'center': [20, 20],
        'pixel_scale': None,
        'mass': 1.0e15,
        'z': 0.3,
    }
    if variance is not None:
        var_path = path.replace('.fits', '_var.fits')
        save_fits(var_path, np.full((41, 41), float(variance)))
        record['variance'] = os.path.basename(var_path)
    return record


def _write(tmp_path, records, **extra):
    cosmo = Cosmology(H0=70.0, Om=0.3)
    scale = float(theta_delta_arcsec(1.0e15, 0.3, cosmo) / 20.0)
    for record in records:
        record['pixel_scale'] = scale
    payload = dict(extra)
    payload['clusters'] = records
    path = tmp_path / 'clusters.json'
    path.write_text(json.dumps(payload))
    return path


def test_cli_stack_writes_the_bcg_diffuse_frames(tmp_path):
    first = _frame(str(tmp_path / 'a.fits'), 4.0)
    second = _frame(str(tmp_path / 'b.fits'), 4.0)
    first['id'] = 'A'
    second['id'] = 'B'
    second['sky'] = 'global'
    second['product'] = 'deepCoadd'
    path = _write(tmp_path, [first, second])
    work = tmp_path / 'iclstack'
    proc = _run(['--clusters', str(path), '--work', str(work)] + FLAGS)
    assert proc.returncode == 0, proc.stderr
    doc = json.loads((work / 'icl_stack.json').read_text())
    assert doc['product'] == 'bcg+diffuse'
    assert doc['n'] == 2
    assert doc['science_measurement'] is False
    assert doc['named_data_mass'] is False
    assert doc['independent_evidence'] is False
    assert doc['inverse_variance'] is False
    assert 'luminosity_unit' in doc['survey_empty']
    image, _header = load_image(str(work / 'icl_stack_u.fits'))
    weight, _header = load_image(str(work / 'icl_stack_u_weight.fits'))
    kpc, _header = load_image(str(work / 'icl_stack_kpc.fits'))
    assert image.shape == (11, 11)
    assert kpc.shape == (11, 11)
    assert np.isfinite(image[5, 5])
    assert weight[5, 5] == pytest.approx(2.0)
    header = (work / 'icl_stack_profile.tsv').read_text().splitlines()[0].split('\t')
    assert header == ['U_MIN', 'U_MAX', 'N_CLUSTERS', 'IDS']
    assert (work / 'icl_stack_kpc_profile.tsv').read_text().splitlines()[0].startswith('KPC_MIN')
    assert 'science_measurement false' in proc.stdout
    assert not (work / 'icl_stack_u_ivar.fits').exists()


def test_inverse_variance_is_a_separate_file(tmp_path):
    faint = _frame(str(tmp_path / 'faint.fits'), 4.0, variance=1.0)
    bright = _frame(str(tmp_path / 'bright.fits'), 12.0, variance=3.0)
    path = _write(tmp_path, [faint, bright])
    work = tmp_path / 'ivar'
    proc = _run(['--clusters', str(path), '--work', str(work)] + FLAGS)
    assert proc.returncode == 0, proc.stderr
    mean, _header = load_image(str(work / 'icl_stack_u.fits'))
    inverse, _header = load_image(str(work / 'icl_stack_u_ivar.fits'))
    assert np.isfinite(mean[5, 5]) and np.isfinite(inverse[5, 5])
    assert inverse[5, 5] != pytest.approx(mean[5, 5])
    doc = json.loads((work / 'icl_stack.json').read_text())
    assert doc['inverse_variance'] is True
    assert doc['science_measurement'] is False


def test_local_sky_and_a_wrong_product_write_no_stack(tmp_path):
    frame = _frame(str(tmp_path / 'a.fits'), 4.0)
    local = dict(frame, sky='local-sky')
    path = _write(tmp_path, [local])
    work = tmp_path / 'rejected'
    proc = _run(['--clusters', str(path), '--work', str(work)] + FLAGS)
    assert proc.returncode == 1
    assert 'local-sky' in proc.stderr
    assert not (work / 'icl_stack.json').exists()
    other = dict(local, sky='global', product='calexp')
    path.write_text(json.dumps({'clusters': [other]}))
    # The pixel scale was filled by _write on the previous record object.
    proc = _run(['--clusters', str(path), '--work', str(work)] + FLAGS)
    assert proc.returncode == 1
    assert 'deepCoadd' in proc.stderr
    assert not (work / 'icl_stack.json').exists()


def test_a_missing_list_writes_nothing(tmp_path):
    work = tmp_path / 'missing'
    proc = _run(['--clusters', str(tmp_path / 'absent.json'), '--work', str(work)] + FLAGS)
    assert proc.returncode == 1
    assert 'was not found' in proc.stderr
    assert not (work / 'icl_stack.json').exists()


def test_menu_step_calls_the_stack_program():
    plugin = json.loads(open(PLUGIN, encoding='utf-8').read())
    step = next(item for item in plugin['steps'] if item['id'] == 'stack')
    assert step['label'] == 'Stack BCG+diffuse'
    cli = [item for item in step['cli'] if isinstance(item, str)]
    assert cli[1].endswith('icl_stack.py')
    assert '--clusters' in cli and '--work' in cli
    assert '--u-max' in cli and '--n-side' in cli and '--delta' in cli
    names = {item['name']: item['default'] for item in plugin['params']}
    assert names['u-max'] == 2.0
    assert names['n-side'] == 201
    assert names['stack-delta'] == 200.0
    assert names['h0'] == 70.0
    assert names['omega-m'] == 0.3
    assert names['clusters'] == ''
