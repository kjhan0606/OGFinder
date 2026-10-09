#!/usr/bin/env python3
"""Fits each object with Prospector (https://github.com/bd-j/prospector, v1.x API: dict obs, SedModel, CSPSpecBasis, dynesty) at the record's redshift and prints the
SED-fit columns as JSON.  Needs `pip install astro-prospector astro-sedpy dynesty fsps` and SPS_HOME (FSPS data) in the interpreter that runs this script.
Model: parametric delay-tau SFH (sfh=4), free mass, logzsol, dust2, tage, tau; Calzetti-like dust2 (dust_type=2)."""
import json
import math
import os
import sys

import numpy as np


def wquant(x, w, q):
    i = np.argsort(x)
    c = np.cumsum(w[i]); c /= c[-1]
    return float(np.interp(q, c, x[i]))


def main():
    try:
        import prospect  # noqa: F401
    except Exception:
        sys.stderr.write(
            "prospect was not found. the external package is not called; "
            "the in-tree SED model is analytic.\n")
        return 2
    req = json.load(sys.stdin)
    P = req.get('params') or {}
    os.chdir(req.get('wd') or '.')
    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors
    from prospect.models.sedmodel import SedModel
    from prospect.sources import CSPSpecBasis
    from prospect.fitting import fit_model
    from prospect.utils.obsutils import fix_obs
    from sedpy.observate import load_filters
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
    from sed_adapters import filters as reg
    bands = req['bands']
    names = [reg.B[b]['sedpy'] for b in bands]
    flt = load_filters(names)
    sps = CSPSpecBasis(zcontinuous=1, compute_vega_mags=False)
    rows = []
    for k, oid in enumerate(req['ids']):
        z = req['z'][k]
        flux = np.array(req['flux'][k]); err = np.array(req['err'][k])
        good = (flux > -90) & (err > 0)
        maggies = np.where(good, flux * 1e-6 / 3631.0, 0.0)           # uJy -> maggies (AB)
        mu = np.where(good, err * 1e-6 / 3631.0, 1.0)
        obs = dict(wavelength=None, spectrum=None, unc=None, filters=flt, maggies=maggies, maggies_unc=mu, phot_mask=good, phot_wave=np.array([f.wave_effective for f in flt]))
        obs = fix_obs(obs)
        mp = TemplateLibrary['parametric_sfh']
        mp['zred'].update(init=float(z), isfree=False)
        mp['mass']['prior'] = priors.LogUniform(mini=10 ** P.get('log_mass_min', 6.0), maxi=10 ** P.get('log_mass_max', 13.0))
        mp['mass']['init'] = 1e9
        mp['logzsol']['prior'] = priors.TopHat(mini=-1.5, maxi=0.3)
        mp['dust2']['prior'] = priors.TopHat(mini=0.0, maxi=2.5)
        mp['tage']['init'] = 3.0
        mp['tage']['prior'] = priors.TopHat(mini=0.05, maxi=max(0.06, min(13.7, 13.8 - 0.0)))
        mp['tau']['prior'] = priors.LogUniform(mini=0.1, maxi=30.0)
        mp['dust_type']['init'] = 2
        model = SedModel(mp)
        out = fit_model(obs, model, sps, optimize=False, emcee=False, dynesty=True, nested_method='rwalk', nested_nlive_init=int(P.get('n_live', 200)),
                        nested_dlogz_init=float(P.get('dlogz', 0.1)), nested_target_n_effective=int(P.get('n_eff', 500)), verbose=False)
        res = out['sampling'][0]
        w = np.exp(res['logwt'] - res['logz'][-1])
        th = np.asarray(res['samples'], float)
        lab = list(model.theta_labels())
        col = lambda n: th[:, lab.index(n)]
        rng = np.random.default_rng(1)
        idx = rng.choice(len(w), size=min(120, len(w)), p=w / w.sum())
        mfrac = np.empty(len(idx)); chi = np.empty(len(idx))
        for j, ii in enumerate(idx):
            _, ph, mf = model.predict(th[ii], obs=obs, sps=sps)
            mfrac[j] = mf
            chi[j] = np.sum(((maggies - ph)[good] / mu[good]) ** 2)
        mass = col('mass')[idx]; tage = col('tage')[idx]; tau = col('tau')[idx]
        lm = np.log10(mass * mfrac)
        la = np.log10(tage * 1e9)
        a = mass / (tau ** 2 * (1 - np.exp(-tage / tau) * (1 + tage / tau)))          # delay-tau amplitude: M_formed = A tau^2 [1 - e^{-T/tau}(1+T/tau)]
        sfr = a * tage * np.exp(-tage / tau) / 1e9
        q = lambda x: (float(np.median(x)), float(0.5 * (np.percentile(x, 84.13) - np.percentile(x, 15.87))))
        rows.append(dict(LOG_MASS=q(lm)[0], LOG_MASS_ERR=q(lm)[1], LOG_AGE=q(la)[0], LOG_AGE_ERR=q(la)[1], LOG_Z=float(np.median(col('logzsol')[idx])),
                         AV=float(np.median(col('dust2')[idx])) * 1.086, SFR=float(np.median(sfr)),
                         SED_CHI2=float(chi.min() / max(int(good.sum()) - 5, 1))))
    try:
        import importlib.metadata as md
        ver = md.version('astro-prospector')
    except Exception:
        ver = 'unknown'
    json.dump({'version': ver, 'rows': rows}, sys.stdout)


if __name__ == '__main__':
    sys.exit(main() or 0)
