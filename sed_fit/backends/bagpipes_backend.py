"""Bagpipes backend."""

import numpy as np
from .base import SPSBackend, SEDResult


class BagpipesBackend(SPSBackend):
    """SED fitting using Bagpipes (Bayesian Analysis of Galaxies).

    Requires: pip install bagpipes
    Features native Bayesian SED fitting via MultiNest/nautilus.
    """

    @classmethod
    def is_available(cls) -> bool:
        try:
            import bagpipes
            return True
        except ImportError:
            return False

    @classmethod
    def name(cls) -> str:
        return "bagpipes"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        import bagpipes
        from .params import canonical_to_bagpipes
        from .filters import map_filters

        params = canonical_to_bagpipes(z, log_mass, log_age, log_Z,
                                        Av, log_tau)
        bp_filters = map_filters(bands, 'bagpipes')

        model = bagpipes.model_galaxy(
            params,
            filt_list=bp_filters,
            spec_wavs=np.arange(1000, 50000, 10.0),
        )

        # model.photometry is in microjanskys -> convert to AB mag
        fluxes = model.photometry  # uJy
        mags = np.where(
            fluxes > 0,
            -2.5 * np.log10(fluxes * 1e-6) + 8.90,  # AB mag from Jy
            99.0,
        )
        return mags.astype(np.float32)

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        """Native Bayesian fitting using bagpipes.fit()."""
        import sys

        try:
            return self._fit_native(mags, mag_errs, photo_z, bands, **kw)
        except Exception as e:
            print(f"WARNING: Bagpipes native fit failed ({e}), "
                  f"falling back to chi2", file=sys.stderr)
            return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                     n_grid=kw.get('n_grid', 5000))

    def _fit_native(self, mags, mag_errs, photo_z, bands, **kw):
        """Run Bagpipes Bayesian fitting for each source."""
        import sys
        import bagpipes
        from .filters import map_filters

        bp_filters = map_filters(bands, 'bagpipes')
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

        # Fit prior specification
        fit_instructions = {
            'redshift': (0.0, 10.0),
            'delayed': {
                'age': (0.001, 15.0),
                'tau': (0.01, 10.0),
                'massformed': (1e6, 1e13),
                'metallicity': (0.0001, 0.06),
            },
            'dust': {
                'type': 'Calzetti',
                'Av': (0.0, 4.0),
            },
        }

        for i in range(n_sources):
            obs = mags[i]
            errs = mag_errs[i]
            z = photo_z[i]

            valid = np.isfinite(obs) & (obs > 0) & (obs < 90)
            if np.sum(valid) < 2:
                continue

            # Fix redshift to photo-z
            fit_instr = dict(fit_instructions)
            if np.isfinite(z) and z > 0:
                fit_instr['redshift'] = (z - 0.1, z + 0.1)

            # Prepare flux array (AB mag -> uJy)
            fluxes = 10 ** (-0.4 * (obs - 8.90)) * 1e6  # uJy
            flux_errs = fluxes * errs * np.log(10) / 2.5  # propagate

            def load_func(idx):
                return np.column_stack([fluxes[valid], flux_errs[valid]])

            try:
                galaxy = bagpipes.galaxy(
                    i, load_func,
                    filt_list=[bp_filters[j] for j, v in enumerate(valid) if v],
                    spectrum_exists=False,
                )
                fit = bagpipes.fit(galaxy, fit_instr, run='temp')
                fit.fit(verbose=False)

                post = fit.posterior
                log_mass_out[i] = np.log10(np.median(
                    post.samples['delayed:massformed']))
                log_mass_err[i] = np.std(np.log10(
                    post.samples['delayed:massformed']))
                log_age_out[i] = np.log10(
                    np.median(post.samples['delayed:age'])) + 9
                log_age_err[i] = 0.3
                log_Z_out[i] = np.log10(
                    np.median(post.samples['delayed:metallicity']) / 0.02)
                Av_out[i] = np.median(post.samples['dust:Av'])
                log_tau_out[i] = np.log10(
                    np.median(post.samples['delayed:tau'])) + 9
                chi2_out[i] = fit.min_chi2

                mass_val = 10 ** log_mass_out[i]
                age_gyr = 10 ** (log_age_out[i] - 9)
                sfr_out[i] = np.log10(
                    max(mass_val / (age_gyr * 1e9), 1e-5))
            except Exception:
                continue

            if (i + 1) % 10 == 0:
                print(f"  bagpipes: {i+1}/{n_sources}", file=sys.stderr)

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
