"""FSPS (python-fsps) backend."""

import numpy as np
from .base import SPSBackend, SEDResult


class FSPSBackend(SPSBackend):
    """Stellar population synthesis using python-fsps.

    Requires: pip install fsps (+ FSPS_DIR environment variable).
    """

    _sp = None  # cached StellarPopulation instance

    @classmethod
    def is_available(cls) -> bool:
        try:
            import fsps
            return True
        except ImportError:
            return False

    @classmethod
    def name(cls) -> str:
        return "fsps"

    def _get_sp(self):
        """Get or create cached StellarPopulation object."""
        if FSPSBackend._sp is None:
            import fsps
            FSPSBackend._sp = fsps.StellarPopulation(
                zcontinuous=1,
                imf_type=1,    # Chabrier
                dust_type=0,   # Calzetti-like
                sfh=4,         # delayed-tau
            )
        return FSPSBackend._sp

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        from .filters import map_filters
        from .params import canonical_to_fsps

        sp = self._get_sp()
        params = canonical_to_fsps(z, log_mass, log_age, log_Z, Av, log_tau)

        sp.params['logzsol'] = params['logzsol']
        sp.params['dust2'] = params['dust2']
        sp.params['tau'] = params['tau']

        fsps_filters = map_filters(bands, 'fsps')
        mags = sp.get_mags(
            tage=params['tage'],
            redshift=params['zred'],
            bands=fsps_filters,
        )

        # Scale by mass (FSPS returns for 1 Msun formed)
        mass_scale = -2.5 * np.log10(params['mass'])
        mags = mags + mass_scale

        return mags.astype(np.float32)

    def generate_grid(self, n_samples, bands, output_path, seed=42,
                      param_ranges=None):
        """Optimized batch generation reusing SP object."""
        import sys
        import h5py
        import os
        from ..config import PHYS_PARAM_NAMES
        from .base import DEFAULT_PARAM_RANGES
        from .filters import map_filters

        sp = self._get_sp()
        rng = np.random.RandomState(seed)
        ranges = dict(DEFAULT_PARAM_RANGES)
        if param_ranges:
            ranges.update(param_ranges)

        n_bands = len(bands)
        fsps_filters = map_filters(bands, 'fsps')

        params = np.zeros((n_samples, 6), dtype=np.float32)
        for i, pname in enumerate(PHYS_PARAM_NAMES):
            lo, hi = ranges[pname]
            params[:, i] = rng.uniform(lo, hi, size=n_samples)

        mags = np.zeros((n_samples, n_bands), dtype=np.float32)
        for j in range(n_samples):
            z_val = float(params[j, 0])
            log_Z_val = float(params[j, 3])
            Av_val = float(params[j, 4])
            tau_gyr = 10 ** (float(params[j, 5]) - 9)
            tage_gyr = 10 ** (float(params[j, 2]) - 9)
            mass = 10 ** float(params[j, 1])

            sp.params['logzsol'] = log_Z_val
            sp.params['dust2'] = Av_val * 0.4
            sp.params['tau'] = tau_gyr

            try:
                m = sp.get_mags(tage=tage_gyr, redshift=z_val,
                                bands=fsps_filters)
                mags[j] = m - 2.5 * np.log10(mass)
            except Exception:
                mags[j] = np.nan

            if (j + 1) % 1000 == 0:
                print(f"  fsps: {j+1}/{n_samples}", file=sys.stderr)

        # Add noise
        noise_scale = 0.01 + 0.05 * np.clip(
            (mags - 20.0) / 10.0, 0.0, 1.0)
        mags += rng.normal(0, noise_scale).astype(np.float32)

        os.makedirs(os.path.dirname(os.path.abspath(output_path)),
                     exist_ok=True)
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('params', data=params)
            f.create_dataset('mags', data=mags)
            f.attrs['param_names'] = PHYS_PARAM_NAMES
            f.attrs['band_names'] = bands
            f.attrs['n_samples'] = n_samples
            f.attrs['n_bands'] = n_bands
            f.attrs['backend'] = 'fsps'

        print(f"Generated {n_samples} SEDs [backend=fsps]", file=sys.stderr)
        print(f"Saved to {output_path}", file=sys.stderr)
        return params, mags

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        n_grid = kw.get('n_grid', 5000)
        return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                 n_grid=n_grid)
