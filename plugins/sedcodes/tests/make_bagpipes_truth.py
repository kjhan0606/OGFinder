"""Run with the interpreter that has Bagpipes: photometry of a Bagpipes model galaxy (delayed-exponential SFH, Calzetti dust) in the test bands, written as JSON."""
import json
import os
import sys

import numpy as np

fdir, nobj = sys.argv[1], int(sys.argv[2])
os.environ.setdefault('MPLBACKEND', 'Agg')
import bagpipes as pipes

bands = ['F435W', 'F606W', 'F775W', 'F850LP', 'F105W', 'F125W', 'F160W']
ff = [os.path.join(fdir, b + '.dat') for b in bands]
rng = np.random.default_rng(11)
out = []
for k in range(nobj):
    z = float(rng.uniform(0.4, 1.2))
    comp = dict(redshift=z, exponential=dict(age=float(rng.uniform(1.0, 6.0)), tau=float(rng.uniform(0.5, 4.0)), massformed=float(rng.uniform(9.8, 10.8)), metallicity=float(rng.uniform(0.6, 1.2))),
                dust=dict(type='Calzetti', Av=float(rng.uniform(0.1, 1.0))))
    m = pipes.model_galaxy(comp, filt_list=ff, spec_wavs=np.arange(4000., 9000., 50.))
    f = np.asarray(m.photometry, float)                      # erg/s/cm2/A
    lam = np.array([np.average(np.loadtxt(p)[:, 0], weights=np.loadtxt(p)[:, 1]) for p in ff])
    fnu = f * lam ** 2 / 2.99792458e18                       # erg/s/cm2/Hz
    ujy = fnu / 1e-29                                         # 1 uJy = 1e-29 erg/s/cm2/Hz
    e = np.maximum(0.05 * ujy, 0.02)
    ujy = ujy + rng.normal(size=ujy.size) * e
    mags = {b: float(23.9 - 2.5 * np.log10(u)) for b, u in zip(bands, ujy)}
    errs = {b: float(1.0857 * ee / u) for b, u, ee in zip(bands, ujy, e)}
    # true stellar mass (surviving), from the model's sfh bookkeeping
    out.append(dict(id=str(k + 1), z=z, mags=mags, mag_errs=errs, logm_formed=comp['exponential']['massformed'], av=comp['dust']['Av'], age_gyr=comp['exponential']['age'],
                    logm=float(m.sfh.stellar_mass) if hasattr(m.sfh, 'stellar_mass') else None))
json.dump(out, sys.stdout)
