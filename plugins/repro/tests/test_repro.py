"""Bundle round trip: batch run -> bundle -> verify (PASS), tampered input / bundle / code environment (FAIL / WARN), table comparison, bundle contents."""
import json
import os
import shutil
import subprocess
import sys
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUG = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import repro  # noqa: E402

CLI = os.path.join(PLUG, 'repro.py')
BATCH = os.path.join(ROOT, 'plugins', 'batch', 'batch.py')
need = pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, 'ogfmeas', 'sextract.py')), reason='ogfmeas/sextract.py missing')


def sh(cmd):
    return subprocess.run([sys.executable] + cmd, capture_output=True, text=True)


@pytest.fixture(scope='module')
def batch_run(tmp_path_factory):
    from astropy.io import fits
    tmp = tmp_path_factory.mktemp('rb')
    lines = []
    for k in range(3):
        rng = np.random.RandomState(k)
        im = rng.randn(250, 250) * 5 + 100
        yy, xx = np.mgrid[0:250, 0:250]
        for i in range(20):
            x, y = rng.rand(2) * 200 + 25
            im += rng.uniform(150, 900) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / 8.0)
        p = str(tmp / ('img%d.fits' % k)); fits.writeto(p, im.astype('float32'), overwrite=True)
        lines.append('R%d %s\n' % (k, p))
    (tmp / 'fields.txt').write_text(''.join(lines))
    (tmp / 'recipe.json').write_text(json.dumps({"detect": {"args": ["--detect-thresh", "3"]}, "steps": [{"plugin": "noisemodel", "step": "model"}, {"plugin": "example_hello", "step": "greet"}]}))
    r = sh([BATCH, '--root', ROOT, '--fields', str(tmp / 'fields.txt'), '--recipe', str(tmp / 'recipe.json'), '--outdir', str(tmp / 'out'), '--jobs', '3'])
    assert r.returncode == 0, r.stdout + r.stderr
    z = str(tmp / 'b.zip')
    r = sh([CLI, 'bundle-batch', '--root', ROOT, '--out', z, '--batch-out', str(tmp / 'out'), '--fields', str(tmp / 'fields.txt'), '--recipe', str(tmp / 'recipe.json')])
    assert r.returncode == 0, r.stdout + r.stderr
    return tmp, z


@need
def test_bundle_contents(batch_run):
    tmp, z = batch_run
    names = zipfile.ZipFile(z).namelist()
    for n in ('manifest.json', 'params.json', 'steps.json', 'requirements.lock', 'README.txt', 'verify.sh', 'recipe.json', 'fields.txt', 'outputs/R0/catalog.tsv'):
        assert n in names, n
    m = repro.read_manifest(z)
    print('bundle: %d files, %d bytes; git %s dirty=%s; %d key packages %s; %d locked packages; inputs %d (sha256 %s...)' % (
        len(names), os.path.getsize(z), m['tool'].get('head', '')[:9], m['tool'].get('dirty'), len(m['packages']), sorted(m['packages']), m['n_packages_locked'], len(m['inputs']), m['inputs'][0]['sha256'][:12]))
    assert m['schema'] == 1 and m['kind'] == 'batch' and len(m['inputs']) == 3 and len(m['outputs']) == 3
    assert 'numpy' in m['packages'] and m['binaries']['ogfmeas_sextract']['sha256']
    assert m['inputs'][0]['sha256'] == repro.sha256_file(m['inputs'][0]['path'])
    lock = zipfile.ZipFile(z).read('requirements.lock').decode()
    assert 'numpy==' in lock and len(lock.splitlines()) == m['n_packages_locked']
    p = json.loads(zipfile.ZipFile(z).read('params.json'))
    assert 'noisemodel.model' in p and p['noisemodel.model']['bw'] == 64


@need
def test_round_trip_verify_passes(batch_run, tmp_path):
    tmp, z = batch_run
    r = sh([CLI, 'verify', z, '--root', ROOT, '--workdir', str(tmp_path / 'w'), '--json', str(tmp_path / 'rep.json')])
    print(r.stdout)
    assert r.returncode == 0, r.stdout + r.stderr
    rep = json.load(open(tmp_path / 'rep.json'))
    st = {c['name']: c['status'] for c in rep['checks']}
    assert rep['ok'] and st['re-run (batch)'] == 'PASS'
    cmp_ = [c for c in rep['checks'] if c['name'].startswith('compare')]
    assert len(cmp_) == 3 and all(c['status'] == 'PASS' and 'identical' in c['detail'] for c in cmp_)
    assert all(st[k] == 'PASS' for k in st if k.startswith('input') or k.startswith('bundled'))


@need
def test_tampered_input_is_detected(batch_run, tmp_path):
    tmp, z = batch_run
    d = tmp_path / 'in'; d.mkdir()
    for n in ('img0.fits', 'img1.fits', 'img2.fits'):
        shutil.copy(tmp / n, d / n)
    from astropy.io import fits
    with fits.open(d / 'img1.fits', mode='update') as h:
        h[0].data[10, 10] += 1.0
    r = sh([CLI, 'verify', z, '--root', ROOT, '--inputs-dir', str(d), '--workdir', str(tmp_path / 'w')])
    assert r.returncode == 1 and 'FAIL  input img1.fits' in r.stdout and 'SKIP  re-run' in r.stdout
    shutil.copy(tmp / 'img1.fits', d / 'img1.fits')                 # restored copy in another directory verifies
    r = sh([CLI, 'verify', z, '--root', ROOT, '--inputs-dir', str(d), '--workdir', str(tmp_path / 'w2')])
    assert r.returncode == 0, r.stdout


@need
def test_altered_bundle_and_changed_processing(batch_run, tmp_path):
    tmp, z = batch_run
    # (a) a bundled catalog edited after the fact
    z2 = str(tmp_path / 'edited.zip')
    with zipfile.ZipFile(z) as zi, zipfile.ZipFile(z2, 'w') as zo:
        for n in zi.namelist():
            data = zi.read(n)
            if n == 'outputs/R1/catalog.tsv':
                data = data.replace(b'\t', b'\t', 1)[:-2] + b'X\n'
            zo.writestr(n, data)
    r = sh([CLI, 'verify', z2, '--root', ROOT, '--workdir', str(tmp_path / 'w'), '--no-rerun'])
    assert r.returncode == 1 and 'FAIL  bundled output R1/catalog.tsv' in r.stdout
    # (b) different processing: the recipe inside the bundle is changed -> the re-run no longer reproduces the catalog
    z3 = str(tmp_path / 'proc.zip')
    with zipfile.ZipFile(z) as zi, zipfile.ZipFile(z3, 'w') as zo:
        for n in zi.namelist():
            data = zi.read(n)
            if n == 'recipe.json':
                data = data.replace(b'"3"', b'"5"')
            zo.writestr(n, data)
    r = sh([CLI, 'verify', z3, '--root', ROOT, '--workdir', str(tmp_path / 'w3')])
    print(r.stdout)
    assert r.returncode == 1 and 'FAIL  compare R0/catalog.tsv' in r.stdout


def test_compare_tables_numeric_tolerance(tmp_path):
    a = tmp_path / 'a.tsv'; b = tmp_path / 'b.tsv'
    a.write_text('NUMBER\tX\tN\n1\t1.0000000\tabc\n2\t2.5000000\tdef\n')
    b.write_text('NUMBER\tX\tN\n1\t1.0000001\tabc\n2\t2.5000000\tdef\n')
    assert repro.compare_tables(str(a), str(a))['status'] == 'identical'
    assert repro.compare_tables(str(a), str(b))['status'] == 'different'
    c = repro.compare_tables(str(a), str(b), rtol=1e-6)
    assert c['status'] == 'numeric' and abs(c['max_rel_diff'] - 1e-7) < 1e-8
    b.write_text('NUMBER\tX\tN\n1\t1.0\tabc\n2\t2.5\tXYZ\n')
    c = repro.compare_tables(str(a), str(b), rtol=1e-3)
    assert c['status'] == 'different' and c['columns_different'] == ['N']
    b.write_text('NUMBER\tX\n1\t1\n')
    assert 'column' in repro.compare_tables(str(a), str(b))['detail']


def test_environment_difference_is_a_warning_or_failure(batch_run, tmp_path):
    tmp, z = batch_run
    z2 = str(tmp_path / 'env.zip')
    with zipfile.ZipFile(z) as zi, zipfile.ZipFile(z2, 'w') as zo:
        for n in zi.namelist():
            data = zi.read(n)
            if n == 'manifest.json':
                m = json.loads(data); m['packages']['numpy'] = '0.0.1'; data = json.dumps(m).encode()
            zo.writestr(n, data)
    rep = repro.verify_bundle(z2, ROOT, rerun=False)
    e = [c for c in rep['checks'] if c['name'] == 'environment'][0]
    assert e['status'] == 'WARN' and 'numpy 0.0.1' in e['detail'] and rep['ok']
    rep = repro.verify_bundle(z2, ROOT, rerun=False, strict_env=True)
    assert not rep['ok']


def test_unreadable_bundle(tmp_path):
    (tmp_path / 'x.zip').write_text('not a zip')
    r = sh([CLI, 'verify', str(tmp_path / 'x.zip'), '--root', ROOT])
    assert r.returncode == 1 and 'cannot read bundle' in r.stdout


def test_catalog_bundle_with_included_inputs(tmp_path):
    from astropy.io import fits
    img = tmp_path / 'a.fits'; fits.writeto(str(img), np.zeros((10, 10), 'float32'))
    cat = tmp_path / 'c.tsv'; cat.write_text('NUMBER\tX\n1\t2\n')
    z = str(tmp_path / 'c.zip')
    r = sh([CLI, 'bundle', '--root', ROOT, '--out', z, '--kind', 'catalog', '--catalog', str(cat), '--inputs', str(img), '--include-inputs', '--note', 'unit test'])
    assert r.returncode == 0, r.stderr
    assert 'inputs/a.fits' in zipfile.ZipFile(z).namelist()
    os.remove(img)                                             # original gone: the bundled copy is used
    rep = repro.verify_bundle(z, ROOT)
    st = {c['name']: c['status'] for c in rep['checks']}
    assert rep['ok'] and st['input a.fits'] == 'PASS'
