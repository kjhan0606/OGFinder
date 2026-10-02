#!/usr/bin/env python3
"""Independent-truth test: images of SIE+shear lenses are produced by lenstronomy's own lens-equation solver (not by ogfkit/lensmodel.py), noise is added
and the images are fitted with ogfkit.  Needs lenstronomy (pip install lenstronomy).   python3 crosscheck_lenstronomy.py [n] [sigma_arcsec]"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from ogfkit import lensmodel as L  # noqa: E402
from lenstronomy.LensModel.lens_model import LensModel as LM  # noqa: E402
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver  # noqa: E402

n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
sig = float(sys.argv[2]) if len(sys.argv) > 2 else 0.003
rng = np.random.default_rng(2024)
rows = []
tries = 0
while len(rows) < n and tries < 500:
    tries += 1
    p = dict(theta_E=rng.uniform(0.9, 2.2), q=rng.uniform(0.55, 0.95), phi=rng.uniform(0, 180), gamma=rng.uniform(0, 0.10), phi_g=rng.uniform(0, 180))
    q = p['q']
    e1, e2 = (1 - q) / (1 + q) * math.cos(2 * math.radians(p['phi'])), (1 - q) / (1 + q) * math.sin(2 * math.radians(p['phi']))
    g1, g2 = p['gamma'] * math.cos(2 * math.radians(p['phi_g'])), p['gamma'] * math.sin(2 * math.radians(p['phi_g']))
    lm = LM(['SIE', 'SHEAR'])
    kw = [dict(theta_E=p['theta_E'], e1=e1, e2=e2, center_x=0, center_y=0), dict(gamma1=g1, gamma2=g2)]
    r = rng.uniform(0.05, 0.25) * p['theta_E']
    a = rng.uniform(0, 2 * math.pi)
    bx, by = r * math.cos(a), r * math.sin(a)
    xs, ys = LensEquationSolver(lm).image_position_from_source(bx, by, kw, min_distance=0.01, search_window=8 * p['theta_E'], precision_limit=1e-10, num_iter_max=100)
    if len(xs) != 4:
        continue
    obs = np.column_stack([xs, ys]) + rng.normal(0, sig, (4, 2))
    f = L.fit_sie_shear(obs, sig, (0.0, 0.0), seed=len(rows))
    t = p
    rows.append(dict(dth=(f['params']['theta_E'] - t['theta_E']) / t['theta_E'], dq=f['params']['q'] - t['q'],
                     dphi=((f['params']['phi'] - t['phi'] + 90) % 180) - 90 if t['q'] < 0.9 else float('nan'), dg=f['params']['gamma'] - t['gamma'],
                     dsrc=math.hypot(f['source'][0] - bx, f['source'][1] - by), rms=f['rms_arcsec'], chi2=f['chi2']))
print('lenstronomy-generated quads: %d, noise %.4f arcsec' % (len(rows), sig))
for k in ('dth', 'dq', 'dphi', 'dg', 'dsrc', 'rms', 'chi2'):
    v = np.array([r[k] for r in rows])
    v = v[np.isfinite(v)]
    print('  %-5s median %+.4g  std %.4g  max|.| %.4g' % (k, np.median(v), np.std(v), np.max(np.abs(v))))
