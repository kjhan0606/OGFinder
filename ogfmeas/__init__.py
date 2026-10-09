"""In-tree source measurements for the standalone application.

The names exported here match the calls already used on background-subtracted
images: background, extract, circular and elliptical sums, Kron radius, flux
radius, and a windowed centroid. Submodules ``sersic``, ``photoz``, and ``nbody``
are separate original fits. Nothing here is a copy of an external extractor,
orbit integrator, or population-synthesis code.
"""
from ogfmeas.aperture import flux_radius, kron_radius, sum_circle, sum_circann, sum_ellipse, winpos
from ogfmeas.background import Background
from ogfmeas.extract import extract, set_extract_pixstack, set_sub_object_limit


def measurement_library():
    """Return the measurement module the caller asked for.

    The default is this package, including when ``sep`` is installed.
    ``OGF_USE_SEP=1`` returns the installed ``sep`` module. A missing
    package raises and is not replaced by this package.
    """
    import os
    import sys
    if os.environ.get("OGF_USE_SEP") != "1":
        return sys.modules[__name__]
    try:
        import sep
    except Exception as exc:
        raise RuntimeError(
            "sep was not found. The in-tree measurements are ogfmeas."
        ) from exc
    return sep


__all__ = [
    "Background",
    "extract",
    "flux_radius",
    "kron_radius",
    "measurement_library",
    "set_extract_pixstack",
    "set_sub_object_limit",
    "sum_circle",
    "sum_circann",
    "sum_ellipse",
    "winpos",
]
