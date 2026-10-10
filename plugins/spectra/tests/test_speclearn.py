"""Synthetic checks for the learned spectrum steps.

Stars are two made-up families, not an MK atlas.  Galaxy names are the labels
given to the neighbour library, not a survey catalog.  The supernova match
uses a temporary template directory and does not ship a library.
"""
import os
import subprocess
import sys

import numpy as np

from ogfkit import speclearn as learn

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'spectra.py')


def _star_wave_flux(letter, seed):
    rng = np.random.default_rng(seed)
    w = np.linspace(3600.0, 9200.0, 1400)
    x = (w - 3600.0) / (9200.0 - 3600.0)
    if letter == 'A':
        f = 2.2 - 1.6 * x
        f = f - 0.55 * np.exp(-0.5 * ((w - 4340.0) / 25.0) ** 2)
    else:
        f = 0.35 + 1.8 * x
        f = f - 0.55 * np.exp(-0.5 * ((w - 7600.0) / 40.0) ** 2)
    f = np.clip(f, 0.05, None) + rng.normal(0.0, 0.01, size=w.shape)
    return w, f


def test_temperature_letter_skips_wd_and_cv():
    assert learn.temperature_letter('F3V') == 'F'
    assert learn.temperature_letter('  k5  ') == 'K'
    assert learn.temperature_letter('') is None
    assert learn.temperature_letter('WD') is None
    assert learn.temperature_letter('CV') is None
    assert learn.temperature_letter('Carbon') is None


def test_star_families_roundtrip(tmp_path):
    train_x, train_y = [], []
    for letter, seed0 in (('A', 10), ('M', 40)):
        for i in range(16):
            w, f = _star_wave_flux(letter, seed0 + i)
            train_x.append(learn.on_grid(w, f, None, 0.0, learn.STAR_WAVE))
            train_y.append(letter)
    assert all(v is not None for v in train_x)
    model = learn.fit_star(np.vstack(train_x), train_y, n_comp=6, seed=1)
    path = tmp_path / 'star.npz'
    learn.save_model(path, model)
    loaded = learn.load_model(path)
    assert loaded['kind'] == 'star'
    n_ok = 0
    for letter, seed0 in (('A', 100), ('M', 130)):
        for i in range(6):
            w, f = _star_wave_flux(letter, seed0 + i)
            an = learn.classify_star(w * 1.001, f, None, 0.001, loaded)
            n_ok += int(an['type'] == letter and an['quality'] == 1)
    assert n_ok == 12


def _galaxy_flux(kind, z, seed):
    rng = np.random.default_rng(seed)
    rest = np.linspace(4000.0, 7200.0, 1600)
    if kind == 'star-forming':
        f = np.ones(rest.shape)
        f = f + 3.0 * np.exp(-0.5 * ((rest - 6563.0) / 6.0) ** 2)
        f = f + 1.1 * np.exp(-0.5 * ((rest - 5007.0) / 6.0) ** 2)
    else:
        f = 0.55 + 0.9 * (rest - 4000.0) / 3200.0
        f = f + 2.4 * np.exp(-0.5 * ((rest - 4861.0) / 45.0) ** 2)
    f = f + rng.normal(0.0, 0.02, size=f.shape)
    return rest * (1.0 + z), f


def test_galaxy_neighbour_vote_separates_two_families():
    xs, labels = [], []
    for i in range(24):
        w, f = _galaxy_flux('star-forming', 0.08, 200 + i)
        xs.append(learn.on_grid(w, f, None, 0.08, learn.GAL_WAVE))
        labels.append('star-forming')
        w, f = _galaxy_flux('broad-line AGN', 0.12, 400 + i)
        xs.append(learn.on_grid(w, f, None, 0.12, learn.GAL_WAVE))
        labels.append('broad-line AGN')
    for i in range(8):
        for kind, z, seed in (('star-forming', 0.08, 600), ('broad-line AGN', 0.12, 700)):
            w, f = _galaxy_flux(kind, z, seed + i)
            xs.append(learn.on_grid(w, f, None, z, learn.GAL_WAVE))
            labels.append('')
    assert all(v is not None for v in xs)
    model = learn.fit_galaxy(np.vstack(xs), labels, latent=4, hidden=32, epochs=60, batch=16, lr=0.04, k=5, seed=2)
    assert np.isfinite(model['loss'])
    n_ok = 0
    for kind, z, seed in (('star-forming', 0.08, 800), ('broad-line AGN', 0.12, 900)):
        for i in range(8):
            w, f = _galaxy_flux(kind, z, seed + i)
            an = learn.classify_galaxy(w, f, None, z, model)
            n_ok += int(an['type'] == kind and an['quality'] == 1)
    assert n_ok >= 14, (n_ok, model['loss'])


def _write_template(path, kind, center):
    w = np.linspace(4300.0, 7200.0, 500)
    f = 1.0 - 0.55 * np.exp(-0.5 * ((w - center) / 70.0) ** 2)
    with open(path, 'w') as fh:
        fh.write('# type %s\n' % kind)
        for wi, fi in zip(w, f):
            fh.write('%.4f %.6f\n' % (wi, fi))


def test_supernova_match_and_flat(tmp_path):
    d = tmp_path / 'tpl'
    d.mkdir()
    _write_template(d / 'Ia_a.txt', 'Ia', 6200.0)
    _write_template(d / 'II_a.txt', 'II', 5000.0)
    templates = learn.load_templates(str(d))
    assert {t['type'] for t in templates} == {'Ia', 'II'}
    z = 0.05
    w = np.linspace(4300.0, 7200.0, 500) * (1.0 + z)
    f = 1.0 - 0.55 * np.exp(-0.5 * ((w / (1.0 + z) - 6200.0) / 70.0) ** 2)
    an = learn.match_supernova(w, f, None, z, templates)
    assert an['type'] == 'Ia', an
    assert an['quality'] == 1 and an['score'] >= learn.SN_MIN_R
    flat = learn.match_supernova(w, np.ones(w.shape), None, z, templates)
    assert flat['type'] == 'uncertain' and flat['quality'] == 0
    missing = learn.match_supernova(w, f, None, None, templates)
    assert missing['type'] == 'unknown' and missing['quality'] == 0
    try:
        learn.load_templates(str(tmp_path / 'absent'))
    except FileNotFoundError as ex:
        assert 'was not found' in str(ex)
    else:
        raise AssertionError('missing template directory did not fail')


def _run(args):
    return subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)


def test_cli_missing_model_writes_nothing(tmp_path):
    cat = tmp_path / 'cat.tsv'
    cat.write_text('NUMBER\tZ\n1\t0.01\n')
    work = tmp_path / 'work'
    p = _run(['--task', 'startype', '--catalog', str(cat), '--work', str(work), '--model', str(tmp_path / 'missing.npz')])
    assert p.returncode == 1
    assert 'was not found' in p.stderr
    assert p.stdout.strip() == ''
    assert not work.exists()
    p2 = _run(['--task', 'snmatch', '--catalog', str(cat), '--work', str(work), '--template-dir', str(tmp_path / 'no-templates')])
    assert p2.returncode == 1
    assert 'was not found' in p2.stderr
    assert p2.stdout.strip() == ''
    assert not work.exists()


def test_cli_startype_column(tmp_path):
    rows, letters = [], []
    for letter, seed0 in (('A', 10), ('M', 40)):
        for i in range(12):
            w, f = _star_wave_flux(letter, seed0 + i)
            rows.append(learn.on_grid(w, f, None, 0.0, learn.STAR_WAVE))
            letters.append(letter)
    model = learn.fit_star(np.vstack(rows), letters, n_comp=4, seed=1)
    model_path = tmp_path / 'star.npz'
    learn.save_model(model_path, model)
    spec = tmp_path / 'spec_1.txt'
    w, f = _star_wave_flux('A', 1000)
    with open(spec, 'w') as fh:
        for wi, fi in zip(w, f):
            fh.write('%.4f %.6f\n' % (wi, fi))
    cat = tmp_path / 'cat.tsv'
    cat.write_text('NUMBER\tZ\n1\t0\n')
    work = tmp_path / 'out'
    p = _run(['--task', 'startype', '--catalog', str(cat), '--work', str(work),
              '--spec-dir', str(tmp_path), '--spec-pattern', 'spec_{NUMBER}.txt',
              '--model', str(model_path), '--z-column', 'Z'])
    assert p.returncode == 0, p.stderr[-500:]
    lines = [ln for ln in p.stdout.splitlines() if ln]
    assert lines[0].split('\t')[:4] == ['NUMBER', 'ST_TYPE', 'ST_P', 'ST_Z']
    assert lines[1].split('\t')[1] == 'A'
    assert (work / 'spectra_startype.json').is_file()
