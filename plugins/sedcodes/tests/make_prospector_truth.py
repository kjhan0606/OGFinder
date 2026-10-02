"""Run with the interpreter that has Prospector + python-fsps (SPS_HOME set): photometry of Prospector model galaxies (parametric delay-tau SFH, dust2) in the
test bands, with 5 % noise, written as JSON with the truth (surviving log mass, log age, Av).  usage: make_prospector_truth.py NOBJ"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
from prospect.models.templates import TemplateLibrary
from prospect.models.sedmodel import SedModel
from prospect.sources import CSPSpecBasis
from prospect.utils.obsutils import fix_obs
from sedpy.observate import load_filters
from sed_adapters import filters as reg

n = int(sys.argv[1])
bands = ['F435W', 'F606W', 'F775W', 'F850LP', 'F105W', 'F125W', 'F160W']
flt = load_filters([reg.B[b]['sedpy'] for b in bands])
sps = CSPSpecBasis(zcontinuous=1, compute_vega_mags=False)
rng = np.random.default_rng(11)
out = []
for k in range(n):
    z = float(rng.uniform(0.4, 1.2))
    mp = TemplateLibrary['parametric_sfh']
    mp['zred'].update(init=z, isfree=False)
    mp['dust_type']['init'] = 2
    model = SedModel(mp)
    th = model.theta.copy()
    lab = list(model.theta_labels())
    truth = dict(mass=10 ** rng.uniform(9.8, 10.8), logzsol=float(rng.uniform(-0.6, 0.0)), dust2=float(rng.uniform(0.1, 0.9)), tage=float(rng.uniform(1.5, 6.0)), tau=float(rng.uniform(0.6, 4.0)))
    for key, v in truth.items():
        th[lab.index(key)] = v
    obs = fix_obs(dict(wavelength=None, spectrum=None, unc=None, filters=flt, maggies=np.ones(len(flt)), maggies_unc=np.ones(len(flt)), phot_mask=np.ones(len(flt), bool),
                       phot_wave=np.array([f.wave_effective for f in flt])))
    _, ph, mfrac = model.predict(th, obs=obs, sps=sps)
    ujy = np.asarray(ph) * 3631.0 * 1e6
    e = np.maximum(0.05 * ujy, 0.02)
    ujy = ujy + rng.normal(size=ujy.size) * e
    mags = {b: float(23.9 - 2.5 * np.log10(u)) for b, u in zip(bands, ujy)}
    errs = {b: float(1.0857 * ee / u) for b, u, ee in zip(bands, ujy, e)}
    out.append(dict(id=str(k + 1), z=z, mags=mags, mag_errs=errs, logm=float(np.log10(truth['mass'] * mfrac)), logage=float(np.log10(truth['tage'] * 1e9)),
                    av=1.086 * truth['dust2'], logzsol=truth['logzsol']))
json.dump(out, sys.stdout)
