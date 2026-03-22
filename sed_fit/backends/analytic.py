"""Analytic SPS backend (always available).

Wraps the existing simple_sed_model from grid/generate.py.
"""

import numpy as np
from .base import SPSBackend, SEDResult


class AnalyticBackend(SPSBackend):
    """Simple analytic SED model (power-law + Calzetti dust).

    Always available — no external dependencies. Uses the existing
    simple_sed_model() function from sed_fit.grid.generate.
    """

    @classmethod
    def is_available(cls) -> bool:
        return True

    @classmethod
    def name(cls) -> str:
        return "analytic"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        from ..grid.generate import simple_sed_model, FILTER_WAVELENGTHS

        wavelengths = np.array([FILTER_WAVELENGTHS[b] for b in bands])
        return simple_sed_model(z, log_mass, log_age, log_Z, Av, log_tau,
                                wavelengths)

    def generate_grid(self, n_samples, bands, output_path, seed=42,
                      param_ranges=None):
        """Use the optimized existing implementation."""
        from ..grid.generate import generate_sps_grid
        return generate_sps_grid(n_samples, bands, output_path, seed=seed)

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        n_grid = kw.get('n_grid', 10000)
        return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                 n_grid=n_grid)
