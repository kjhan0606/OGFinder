"""SPS backend registry.

Provides get_backend(), list_available(), list_all() for discovering
and instantiating SPS backends.
"""

from .base import SPSBackend, SEDResult
from .analytic import AnalyticBackend
from .fsps_backend import FSPSBackend
from .bagpipes_backend import BagpipesBackend
from .prospector_backend import ProspectorBackend
from .cigale_backend import CIGALEBackend
from .dense_basis_backend import DenseBasisBackend

# Priority order (most capable first, analytic last as fallback)
_REGISTRY = [
    FSPSBackend,
    ProspectorBackend,
    BagpipesBackend,
    CIGALEBackend,
    DenseBasisBackend,
    AnalyticBackend,
]


def get_backend(name=None) -> SPSBackend:
    """Get an SPS backend instance.

    Parameters
    ----------
    name : str or None
        Backend name. If None or 'auto', returns the first available
        backend (in priority order). Use 'analytic' for the always-
        available fallback.

    Returns
    -------
    SPSBackend instance.

    Raises
    ------
    ValueError
        If the requested backend is not available.
    """
    if name is None or name == 'auto':
        for cls in _REGISTRY:
            if cls.is_available():
                return cls()
        # Should never happen (AnalyticBackend is always available)
        return AnalyticBackend()

    name_lower = name.lower()
    for cls in _REGISTRY:
        if cls.name() == name_lower:
            if not cls.is_available():
                raise ValueError(
                    f"Backend '{name}' is not available. "
                    f"Install the required package.")
            return cls()

    raise ValueError(
        f"Unknown backend '{name}'. "
        f"Available: {list_all()}")


def list_available() -> list:
    """Return names of backends whose packages are installed."""
    return [cls.name() for cls in _REGISTRY if cls.is_available()]


def list_all() -> list:
    """Return names of all registered backends."""
    return [cls.name() for cls in _REGISTRY]
