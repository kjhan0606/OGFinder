#!/usr/bin/env python3
"""Fits each object with Bagpipes (https://github.com/ACCarnall/bagpipes) at a FIXED redshift and prints the SED-fit
columns as JSON.  stdin: {"ids", "z", "bands", "flux" (uJy), "err" (uJy), "filter_files" {band: path}, "params"}.
The sampler is nautilus. MultiNest is not called. Model: delayed-exponential SFH + Calzetti dust + nebular emission, priors from params (defaults below).  Posterior summaries: median, half of the 16-84 interval."""
import json
import os
import sys
import tempfile

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from sed_fit.backends.bagpipes_backend import nautilus_sampler


def main():
    try:
        import bagpipes  # noqa: F401
    except Exception:
        sys.stderr.write(
            "bagpipes was not found. the external package is not called; "
            "the in-tree SED model is analytic.\n")
        return 2
    req = json.load(sys.stdin)
    P = req.get('params') or {}
    try:
        sampler = nautilus_sampler(P.get('sampler'))
    except RuntimeError as exc:
        sys.stderr.write(str(exc) + "\n")
        return 2
    os.chdir(req.get('wd') or tempfile.mkdtemp(prefix='bagpipes_'))
    os.environ.setdefault('MPLBACKEND', 'Agg')
    import bagpipes as pipes
    bands = req['bands']
    ffiles = [req['filter_files'][b] for b in bands]
    rows = []
    for k, oid in enumerate(req['ids']):
        z = req['z'][k]
        if z is None:
            rows.append(dict(error='no redshift'))
            continue
        flux = np.array(req['flux'][k]); err = np.array(req['err'][k])
        use = np.where((flux > -90) & (err > 0))[0]
        photo = np.array([[flux[j], err[j]] for j in use])
        fl = [ffiles[j] for j in use]

        def load_data(_id, photo=photo):
            return photo

        # Bagpipes phot_units 'mujy' = microjansky (default)
        dust = {'type': 'Calzetti', 'Av': tuple(P.get('av_range', (0., 3.)))}
        exp = {'age': tuple(P.get('age_range', (0.1, 13.0))), 'tau': tuple(P.get('tau_range', (0.1, 10.0))), 'massformed': tuple(P.get('mass_range', (6., 13.))),
               'metallicity': tuple(P.get('metallicity_range', (0.2, 2.5)))}
        inst = {'redshift': float(z), 'exponential': exp, 'dust': dust}
        if P.get('nebular', False):
            inst['nebular'] = {'logU': -3.0}
        gal = pipes.galaxy(str(oid), load_data, filt_list=fl, spectrum_exists=False, phot_units='mujy')
        import hashlib
        run = str(P.get('run', 'adapter')) + '_' + hashlib.sha1(json.dumps([photo.tolist(), z, sorted(inst)]).encode()).hexdigest()[:10]
        fit = pipes.fit(gal, inst, run=run)
        fit.fit(verbose=False, n_live=int(P.get('n_live', 500)), sampler=sampler)
        s = fit.posterior.samples
        def q(name):
            a = np.asarray(s[name], float)
            lo, md, hi = np.percentile(a, [15.87, 50, 84.13])
            return float(md), float(0.5 * (hi - lo))
        lm, lme = q('stellar_mass')
        age = np.asarray(s['exponential:age'], float)        # Gyr
        la = np.log10(age * 1e9)
        sfr = np.asarray(s['sfr'], float) if 'sfr' in s else np.full(age.size, np.nan)
        try:       # Bagpipes: lnlike = -chi2/2 + K,  K = -0.5 sum ln(2 pi sigma^2)  (best sample)
            kk = -0.5 * float(np.sum(np.log(2 * np.pi * np.asarray(fit.galaxy.photometry)[:, 2] ** 2)))
            chi = -2.0 * (float(np.max(np.asarray(fit.results['lnlike'], float))) - kk)
        except Exception as exc:
            sys.stderr.write('chi2 unavailable: %r\n' % (exc,))
            chi = float('nan')
        rows.append(dict(LOG_MASS=lm, LOG_MASS_ERR=lme, LOG_AGE=float(np.median(la)), LOG_AGE_ERR=float(0.5 * np.subtract(*np.percentile(la, [84.13, 15.87]))),
                         LOG_Z=float(np.log10(np.median(np.asarray(s['exponential:metallicity'], float)))), AV=float(np.median(np.asarray(s['dust:Av'], float))),
                         SFR=float(np.median(sfr)), SED_CHI2=chi / max(len(use) - 4, 1)))
    try:
        import importlib.metadata as md
        ver = md.version('bagpipes')
    except Exception:
        ver = getattr(pipes, '__version__', 'unknown')
    json.dump({'version': ver, 'rows': rows}, sys.stdout)


if __name__ == '__main__':
    sys.exit(main() or 0)
