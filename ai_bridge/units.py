"""Declared-unit conversion for response fields (small, explicit table; no guessing)."""

_ANGLE = {"rad": 1.0, "deg": 3.141592653589793 / 180.0, "arcmin": 3.141592653589793 / 180.0 / 60.0,
          "arcsec": 3.141592653589793 / 180.0 / 3600.0, "mas": 3.141592653589793 / 180.0 / 3.6e6}
_FLUX = {"Jy": 1.0, "mJy": 1e-3, "uJy": 1e-6, "nJy": 1e-9}
_LEN = {"m": 1.0, "km": 1e3, "cm": 1e-2, "AA": 1e-10, "nm": 1e-9, "um": 1e-6}
_TIME = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0}
_TABLES = [_ANGLE, _FLUX, _LEN, _TIME]


def factor(unit_in, unit_out):
    """Multiplicative factor unit_in -> unit_out.  ValueError when unknown / incompatible."""
    if unit_in == unit_out:
        return 1.0
    for t in _TABLES:
        if unit_in in t and unit_out in t:
            return t[unit_in] / t[unit_out]
    raise ValueError("cannot convert %r to %r (supported: %s)" % (
        unit_in, unit_out, "; ".join(",".join(t) for t in _TABLES)))
