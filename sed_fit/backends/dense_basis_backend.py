"""Dense Basis backend."""

import numpy as np
from .base import SPSBackend, SEDResult


class DenseBasisBackend(SPSBackend):
    """SED fitting using Dense Basis (GP-based SFH atlas).

    Requires: pip install dense_basis
    Very fast fitting (~7ms/source) using pre-computed atlas.
    """

    _atlas = None

    @classmethod
    def is_available(cls) -> bool:
        try:
            import dense_basis
            return True
        except ImportError:
            return False

    @classmethod
    def name(cls) -> str:
        return "dense_basis"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        import dense_basis as db
        from .filters import map_filters

        db_filters = map_filters(bands, 'dense_basis')

        # Use dense_basis to generate SED at given parameters
        # Dense Basis uses GP-based SFH, we approximate with delayed-tau
        priors = {
            'zred': z,
            'logmass': log_mass,
            'logzsol': log_Z,
            'dust2': Av * 0.4,
        }

        try:
            sed = db.generate_sed(priors, db_filters)
            # sed returns maggies -> AB mag
            mags = np.where(sed > 0, -2.5 * np.log10(sed), 99.0)
        except Exception:
            mags = np.full(len(bands), 99.0)

        return mags.astype(np.float32)

    def _get_atlas(self, bands):
        """Load or create cached atlas for fitting."""
        import dense_basis as db
        from .filters import map_filters

        if DenseBasisBackend._atlas is None:
            db_filters = map_filters(bands, 'dense_basis')
            DenseBasisBackend._atlas = db.load_atlas(db_filters)
        return DenseBasisBackend._atlas

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        """Atlas-based fitting (~7ms per source)."""
        import sys

        try:
            return self._fit_atlas(mags, mag_errs, photo_z, bands, **kw)
        except Exception as e:
            print(f"WARNING: Dense Basis fit failed ({e}), "
                  f"falling back to chi2", file=sys.stderr)
            return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                     n_grid=kw.get('n_grid', 5000))

    def _fit_atlas(self, mags, mag_errs, photo_z, bands, **kw):
        """Run dense_basis atlas fitting."""
        import sys
        import dense_basis as db
        from .filters import map_filters

        db_filters = map_filters(bands, 'dense_basis')
        atlas = self._get_atlas(bands)
        n_sources = mags.shape[0]

        log_mass_out = np.full(n_sources, np.nan)
        log_mass_err_out = np.full(n_sources, np.nan)
        log_age_out = np.full(n_sources, np.nan)
        log_age_err_out = np.full(n_sources, np.nan)
        log_Z_out = np.full(n_sources, np.nan)
        Av_out = np.full(n_sources, np.nan)
        log_tau_out = np.full(n_sources, np.nan)
        sfr_out = np.full(n_sources, np.nan)
        chi2_out = np.full(n_sources, np.nan)

        for i in range(n_sources):
            obs = mags[i]
            errs = mag_errs[i]
            z = photo_z[i]

            valid = np.isfinite(obs) & (obs > 0) & (obs < 90)
            if np.sum(valid) < 2:
                continue

            # Convert to maggies
            maggies = 10 ** (-0.4 * obs)
            err_default = np.where(np.isfinite(errs) & (errs > 0),
                                    errs, 0.1)
            maggies_unc = maggies * err_default * np.log(10) / 2.5

            try:
                result = db.fit(maggies, maggies_unc,
                                z if np.isfinite(z) and z > 0 else 0.1,
                                atlas)
                log_mass_out[i] = result.get('logmass', np.nan)
                log_mass_err_out[i] = result.get('logmass_err', 0.3)
                log_Z_out[i] = result.get('logzsol', np.nan)
                Av_out[i] = result.get('dust2', np.nan)
                if np.isfinite(Av_out[i]):
                    Av_out[i] /= 0.4
                sfr_out[i] = result.get('logsfr', np.nan)
                chi2_out[i] = result.get('chi2', np.nan)

                # Approximate age from mass-weighted age if available
                if 'mwa' in result:
                    log_age_out[i] = np.log10(
                        max(result['mwa'], 1e-3)) + 9
                    log_age_err_out[i] = 0.3
            except Exception:
                continue

            if (i + 1) % 100 == 0:
                print(f"  dense_basis: {i+1}/{n_sources}", file=sys.stderr)

        return SEDResult(
            log_mass=log_mass_out,
            log_mass_err=log_mass_err_out,
            log_age=log_age_out,
            log_age_err=log_age_err_out,
            log_Z=log_Z_out,
            Av=Av_out,
            log_tau=log_tau_out,
            sfr=sfr_out,
            chi2=chi2_out,
        )
