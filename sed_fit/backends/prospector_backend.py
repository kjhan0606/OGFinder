"""Prospector backend."""

import numpy as np
from .base import SPSBackend, SEDResult


class ProspectorBackend(SPSBackend):
    """SED fitting using Prospector (prospect).

    Requires: pip install prospect dynesty
    Features native fitting via dynesty nested sampling.
    """

    @classmethod
    def is_available(cls) -> bool:
        try:
            import prospect
            return True
        except ImportError:
            return False

    @classmethod
    def name(cls) -> str:
        return "prospector"

    def _build_model(self, z=None):
        """Build a simple Prospector parametric model."""
        from prospect.models import SpecModel
        from prospect.models.templates import TemplateLibrary

        model_params = TemplateLibrary['parametric_sfh']
        if z is not None and np.isfinite(z) and z > 0:
            model_params['zred']['init'] = z
            model_params['zred']['isfree'] = False

        model_params['imf_type'] = {'N': 1, 'isfree': False, 'init': 1}
        return SpecModel(model_params)

    def _build_sps(self):
        """Build CSPSpecBasis."""
        from prospect.sources import CSPSpecBasis
        return CSPSpecBasis(zcontinuous=1)

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        from .filters import map_filters
        from .params import canonical_to_prospector
        from sedpy.observate import load_filters

        params = canonical_to_prospector(z, log_mass, log_age, log_Z,
                                          Av, log_tau)
        sedpy_filters = map_filters(bands, 'prospector')
        filters = load_filters(sedpy_filters)

        sps = self._build_sps()
        model = self._build_model(z)

        theta = model.theta.copy()
        model.params['logzsol'] = np.array([params['logzsol']])
        model.params['dust2'] = np.array([params['dust2']])
        model.params['tage'] = np.array([params['tage']])
        model.params['tau'] = np.array([params['tau']])
        model.params['mass'] = np.array([params['mass']])

        spec, phot, extras = model.predict(model.theta, sps=sps,
                                            obs={'filters': filters})
        # phot is in maggies -> AB mag
        mags = np.where(phot > 0, -2.5 * np.log10(phot), 99.0)
        return mags.astype(np.float32)

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        """Native fitting using dynesty nested sampling."""
        import sys

        try:
            return self._fit_native(mags, mag_errs, photo_z, bands, **kw)
        except Exception as e:
            print(f"WARNING: Prospector fit failed ({e}), "
                  f"falling back to chi2", file=sys.stderr)
            return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                     n_grid=kw.get('n_grid', 5000))

    def _fit_native(self, mags, mag_errs, photo_z, bands, **kw):
        """Run Prospector dynesty fitting for each source."""
        import sys
        from prospect.fitting import fit_model
        from sedpy.observate import load_filters
        from .filters import map_filters

        sedpy_filters = map_filters(bands, 'prospector')
        filters = load_filters(sedpy_filters)
        sps = self._build_sps()
        n_sources = mags.shape[0]

        log_mass_out = np.full(n_sources, np.nan)
        log_mass_err = np.full(n_sources, np.nan)
        log_age_out = np.full(n_sources, np.nan)
        log_age_err = np.full(n_sources, np.nan)
        log_Z_out = np.full(n_sources, np.nan)
        Av_out = np.full(n_sources, np.nan)
        log_tau_out = np.full(n_sources, np.nan)
        sfr_out = np.full(n_sources, np.nan)
        chi2_out = np.full(n_sources, np.nan)

        for i in range(n_sources):
            obs_m = mags[i]
            errs = mag_errs[i]
            z = photo_z[i]

            valid = np.isfinite(obs_m) & (obs_m > 0) & (obs_m < 90)
            if np.sum(valid) < 2:
                continue

            # Convert mags to maggies
            maggies = 10 ** (-0.4 * obs_m)
            err_default = np.where(np.isfinite(errs) & (errs > 0),
                                    errs, 0.1)
            maggies_unc = maggies * err_default * np.log(10) / 2.5

            obs_dict = {
                'filters': [filters[j] for j in range(len(bands)) if valid[j]],
                'maggies': maggies[valid],
                'maggies_unc': maggies_unc[valid],
            }

            model = self._build_model(z if np.isfinite(z) and z > 0
                                       else None)

            try:
                output = fit_model(obs_dict, model, sps,
                                    nlive_init=100, nested_maxiter=500)
                result = output['sampling'][0]
                chain = result.samples
                weights = np.exp(result.logwt - result.logz[-1])

                # Extract medians
                for k, pname in enumerate(model.free_params):
                    vals = chain[:, k]
                    med = np.average(vals, weights=weights)
                    if pname == 'mass':
                        log_mass_out[i] = np.log10(max(med, 1))
                        log_mass_err[i] = 0.3
                    elif pname == 'tage':
                        log_age_out[i] = np.log10(max(med, 1e-3)) + 9
                        log_age_err[i] = 0.3
                    elif pname == 'logzsol':
                        log_Z_out[i] = med
                    elif pname == 'dust2':
                        Av_out[i] = med / 0.4
                    elif pname == 'tau':
                        log_tau_out[i] = np.log10(max(med, 1e-3)) + 9

                chi2_out[i] = output.get('chi2_best', np.nan)

                if np.isfinite(log_mass_out[i]) and np.isfinite(
                        log_age_out[i]):
                    mass_val = 10 ** log_mass_out[i]
                    age_gyr = 10 ** (log_age_out[i] - 9)
                    sfr_out[i] = np.log10(
                        max(mass_val / (age_gyr * 1e9), 1e-5))
            except Exception:
                continue

            if (i + 1) % 5 == 0:
                print(f"  prospector: {i+1}/{n_sources}", file=sys.stderr)

        return SEDResult(
            log_mass=log_mass_out,
            log_mass_err=log_mass_err,
            log_age=log_age_out,
            log_age_err=log_age_err,
            log_Z=log_Z_out,
            Av=Av_out,
            log_tau=log_tau_out,
            sfr=sfr_out,
            chi2=chi2_out,
        )
