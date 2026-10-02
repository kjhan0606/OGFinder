"""Merger / interaction indicators (morphometry.features + `ds9_morph_ext.py --features`): asymmetry, peaks, Lotz class, shells, tails, pairs; controls must stay clean."""
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from ogfkit import models as M, tsvio
from morphometry import extended as X, features as FT

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
PY = sys.executable
sys.path.insert(0, os.path.join(HERE, '..', 'validation'))
import features_validate as FV  # noqa: E402


def test_lotz_classification_lines():
    # at M20 = -2: merger line G = 0.61, E/S0/Sa line G = 0.52
    assert FT.lotz_class(0.65, -2.0) == 'merger'
    assert FT.lotz_class(0.55, -2.0) == 'E/S0/Sa'
    assert FT.lotz_class(0.40, -2.0) == 'Sb-Irr'
    assert FT.lotz_class(float('nan'), 0.0) == 'unknown'


def _host(n=2.5, re_=11.0, q=0.8, pa=30.0, flux=1e5, shape=(241, 241), xc=120.0, yc=120.0):
    return M.render_sersic(shape, xc, yc, M.sersic_Ie_from_flux(flux, re_, n, q), re_, n, q, pa)


def test_asymmetry_is_background_corrected_and_grows_with_companion():
    rng = np.random.default_rng(3)
    sh = (241, 241)
    ap = FT._ell_coords(sh, 120, 120, 1.0, 0.0)[2] <= 40
    vals = []
    for ratio in (0.0, 0.1, 0.3, 1.0):
        img = _host() + (M.render_sersic(sh, 142, 125, M.sersic_Ie_from_flux(1e5 * ratio, 5, 1, 0.9), 5, 1, 0.9, 0) if ratio else 0)
        a = FT.asymmetry(img + rng.normal(0, 1.0, sh), ap, 120, 120, rng.normal(0, 1.0, sh))[0]
        vals.append(a)
    assert abs(vals[0]) < 0.02 and all(np.diff(vals) > 0) and vals[3] > 0.2


def test_peaks_find_double_nucleus():
    sh = (201, 201)
    rng = np.random.default_rng(1)
    one = _host(shape=sh, xc=100, yc=100) + rng.normal(0, 1, sh)
    two = one + M.render_sersic(sh, 114, 100, M.sersic_Ie_from_flux(3e4, 2.5, 1, 1.0), 2.5, 1, 1.0, 0)
    ap = FT._ell_coords(sh, 100, 100, 1, 0)[2] <= 40
    assert len(FT.peaks(one, ap, 1.0, 2.0, 100, 100)) == 1
    assert len(FT.peaks(two, ap, 1.0, 2.0, 100, 100)) == 2


def test_shell_tail_detection_and_clean_controls():
    det_s = det_t = fp = 0
    n = 8
    for s in range(n):
        for kind in ('control', 'shell', 'tail'):
            r = FV.trial((kind, s, 0.0 if kind == 'control' else 0.05, 2.0))
            if kind == 'control':
                fp += (r['n_shell'] + r['n_tail']) > 0
            elif kind == 'shell':
                det_s += r['n_shell'] > 0
            else:
                det_t += r['n_tail'] > 0
    assert fp == 0 and det_s >= n - 2 and det_t >= n - 2, (fp, det_s, det_t)


def test_merger_flags_on_equal_mass_merger_but_not_on_minor_companion():
    flags = {}
    for key, (ratio, sep) in {'major': (1.0, 1.6), 'minor': (0.1, 2.5)}.items():
        flags[key] = [FV.trial(('merger', s, (ratio, sep), 2.0))['flags'] for s in range(6)]
    assert all(f & 3 for f in flags['major']) and sum(bool(f & 3) for f in flags['minor']) <= 2


def test_find_pairs_classification():
    x = np.array([100.0, 130.0, 400.0, 100.0, 300.0])
    y = np.array([100.0, 100.0, 400.0, 160.0, 300.0])
    rp = np.array([20.0, 15.0, 20.0, 5.0, 10.0])
    fl = np.array([1000.0, 600.0, 1000.0, 50.0, 1000.0])
    p = FT.find_pairs(x, y, rp, fl)
    assert p['kind'][0] == 2 and p['kind'][1] == 2 and p['touching'][0] and abs(p['sep'][0] - 30.0) < 1e-9     # 1000 vs 600: major
    assert p['kind'][2] == 0 and p['kind'][4] == 0 and p['n_comp'][2] == 0                                        # isolated
    assert p['kind'][3] == 0                                                                                       # flux ratio 1:20 is not even a minor pair


def test_driver_features_only_on_synthetic_field(tmp_path):
    shape = (400, 400)
    rng = np.random.default_rng(5)
    img = np.zeros(shape)
    specs = [(100, 100, 2.5, 11, 0.8, 30, 1e5), (135, 108, 1.0, 8, 0.9, 0, 6e4), (300, 100, 2.5, 11, 0.8, 30, 1e5), (300, 300, 2.5, 11, 0.8, 30, 1e5)]
    for x, y, n, re_, q, pa, fl in specs:
        img += M.render_sersic(shape, x, y, M.sersic_Ie_from_flux(fl, re_, n, q), re_, n, q, pa)
    rr = FT._ell_coords(shape, 300, 300, 0.8, math.radians(30))[2]
    arc = np.exp(-0.5 * ((rr - 18) / 2.5) ** 2) * (np.abs(np.arctan2(*(np.mgrid[:400, :400] - 300)[::-1] if False else (np.mgrid[:400, :400][0] - 300, np.mgrid[:400, :400][1] - 300))) < 1.0)
    img += arc / arc.sum() * 6e3                                                              # a shell on object 4 (6 % of its flux)
    img += 20.0 + rng.normal(0, 1.0, shape)
    f = str(tmp_path / 'f.fits')
    fits.writeto(f, img.astype(np.float32), overwrite=True)
    rows = [dict(NUMBER=k + 1, X_IMAGE=x + 1, Y_IMAGE=y + 1, A_IMAGE=0.6 * re_, B_IMAGE=0.6 * re_ * q, THETA_IMAGE=pa, KRON_RADIUS=3.5, FLUX_AUTO=fl, FLUX_RADIUS=re_ * 0.9)
            for k, (x, y, n, re_, q, pa, fl) in enumerate(specs)]
    cat = str(tmp_path / 'c.tsv')
    tsvio.write_table(cat, list(rows[0]), rows)
    r = subprocess.run([PY, os.path.join(ROOT, 'ds9', 'library', 'ds9_morph_ext.py'), f, '--catalog', cat, '--work', str(tmp_path / 'w'), '--features-only', '--n-workers', '2', '--curves', '0'],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr[-1500:]
    L = r.stdout.strip().splitlines()
    head = L[0].split('\t')
    assert head[0] == 'NUMBER' and 'MX_ASYM' in head and 'MX_MORPH_FLAGS' in head and 'MX_RP' not in head
    rec = [dict(zip(head, l.split('\t'))) for l in L[1:]]
    flags = [int(float(x['MX_MORPH_FLAGS'])) for x in rec]
    assert flags[0] & FT.FLAGS['PAIR_MAJOR'] and flags[1] & (FT.FLAGS['PAIR_MAJOR'] | FT.FLAGS['PAIR_MINOR'])         # objects 1 and 2 are a close pair
    assert not flags[2] & (FT.FLAGS['PAIR_MAJOR'] | FT.FLAGS['PAIR_MINOR'] | FT.FLAGS['SHELL'] | FT.FLAGS['TAIL'])      # isolated smooth control
    assert int(rec[3]['MX_NSHELL']) >= 1 and float(rec[2]['MX_ASYM']) < 0.05
