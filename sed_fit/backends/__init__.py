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

# dense_basis stays disconnected. The other names are explicit requests:
# fsps, bagpipes, and prospector import the user-installed package, and
# cigale runs the user-installed program. auto skips every one of them
# and does not call is_available() on them, so an installed package cannot
# replace the in-tree analytic model.
_NOT_CALLED = frozenset({"dense_basis"})
_EXPLICIT_IMPORT = frozenset({"fsps", "bagpipes", "prospector"})
_AUTO_SKIP = frozenset({
    "fsps", "bagpipes", "prospector", "cigale", "dense_basis",
})
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
        Backend name. None or 'auto' returns the in-tree analytic model.
        It does not select fsps, bagpipes, prospector, cigale, or
        dense_basis, even when those packages are installed. An explicit
        fsps, bagpipes, or prospector name imports that package. cigale
        runs ``pcigale`` when that name is requested and the program is
        on PATH. dense_basis is not called.

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
            if cls.name() in _AUTO_SKIP:
                continue
            if cls.is_available():
                return cls()
        # Should never happen (AnalyticBackend is always available)
        return AnalyticBackend()

    name_lower = name.lower()
    for cls in _REGISTRY:
        if cls.name() == name_lower:
            if not cls.is_available():
                if name_lower in _NOT_CALLED:
                    raise ValueError(
                        f"Backend '{name}' does not call the external package. "
                        f"The in-tree SED model is analytic.")
                if name_lower == "cigale":
                    raise ValueError(
                        "Backend 'cigale': the external program 'pcigale' "
                        "was not found on PATH. "
                        "The in-tree SED model is analytic.")
                if name_lower in _EXPLICIT_IMPORT:
                    raise ValueError(
                        f"Backend '{name}' was not found. "
                        f"The in-tree SED model is analytic.")
                raise ValueError(
                    f"Backend '{name}' is not available. "
                    f"Install the required package.")
            return cls()

    raise ValueError(
        f"Unknown backend '{name}'. "
        f"Available: {list_all()}")


def list_available() -> list:
    """Return names of backends that can be called in-process."""
    return [cls.name() for cls in _REGISTRY if cls.is_available()]


def list_all() -> list:
    """Return names of all registered backends."""
    return [cls.name() for cls in _REGISTRY]
