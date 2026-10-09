"""Registry name for an external stellar-population backend.

The external package is not called. SED fits use the in-tree analytic model
(``get_backend('analytic')`` or ``get_backend('auto')``).
"""

from .base import SPSBackend

_MSG = (
    "Backend 'dense_basis' does not call the external package. "
    "The in-tree SED model is analytic."
)


class DenseBasisBackend(SPSBackend):
    @classmethod
    def is_available(cls) -> bool:
        return False

    @classmethod
    def name(cls) -> str:
        return "dense_basis"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau, bands):
        raise RuntimeError(_MSG)

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw):
        raise RuntimeError(_MSG)
