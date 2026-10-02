#!/usr/bin/env python3
"""Run pysersic (find_MAP, flat sky) on the case-A images written by galfit_compare.py.   Needs `pip install pysersic` (own venv):
    python pysersic_map.py WORKDIR OUT.json [--n 30]
pysersic conventions (checked against GALFIT-rendered images): xc = x_GALFIT - 1, flux in counts, ellip = 1 - q, theta = radians(PA_GALFIT)."""
import json
import math
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
import jax  # noqa: E402
sys.path.insert(0, __import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__file__), '..', '..', '..')))
from ogfkit import galfitio as GI  # noqa: E402
from pysersic import FitSingle  # noqa: E402
from pysersic.priors import PySersicSourcePrior  # noqa: E402

work, out = sys.argv[1], sys.argv[2]
n = int(sys.argv[sys.argv.index('--n') + 1]) if '--n' in sys.argv else 30
res = []
for s in range(n):
    wd = '%s/A%02d' % (work, s)
    data = fits.getdata(wd + '/data.fits').astype(np.float32)
    psf = fits.getdata(wd + '/psf.fits').astype(np.float32)
    rms = np.full_like(data, 2.0)
    g = GI.parse_feedme(wd + '/start.feedme', base_dir=wd)
    c0 = g['components'][0]
    st = dict(x=c0['x'], y=c0['y'], flux=10 ** (-0.4 * (c0['mag'] - 25.0)), re=c0['re'], sky=g['sky_value_galfit'])
    t = time.time()
    try:
        prior = PySersicSourcePrior(profile_type='sersic', sky_type='flat', sky_guess=st['sky'], sky_guess_err=2.0)       # priors centred on the SAME start values given to GALFIT / multifit
        prior.set_gaussian_prior('xc', st['x'] - 1, 2.0).set_gaussian_prior('yc', st['y'] - 1, 2.0)
        prior.set_gaussian_prior('flux', st['flux'], 0.5 * st['flux']).set_truncated_gaussian_prior('r_eff', st['re'], 0.5 * st['re'], 0.3, 50.0)
        prior.set_uniform_prior('n', 0.5, 6.0).set_uniform_prior('ellip', 0.0, 0.8).set_uniform_prior('theta', -2 * math.pi, 2 * math.pi)
        fit = FitSingle(data=data, rms=rms, psf=psf, prior=prior)
        m = fit.find_MAP(rkey=jax.random.PRNGKey(s))
        p = {k: float(np.asarray(v)) for k, v in m.items() if np.ndim(v) == 0 or np.size(v) == 1}
        res.append(dict(seed=s, time=time.time() - t, x=p['xc'] + 1, y=p['yc'] + 1, mag=25.0 - 2.5 * math.log10(p['flux']), re=p['r_eff'], n=p['n'], q=1 - p['ellip'],
                        pa=((math.degrees(p['theta']) + 90.0) % 180.0) - 90.0, sky=p.get('sky_back'), raw=p))
    except Exception as e:
        res.append(dict(seed=s, error=str(e)[:200], time=time.time() - t))
    json.dump(res, open(out, 'w'))
    print(s, res[-1].get('mag'), '%.1fs' % (time.time() - t), flush=True)
