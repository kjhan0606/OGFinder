"""Filter name mapping between canonical names and SPS code conventions."""

# Canonical -> code-specific filter name mappings
# Canonical names are HST/JWST filter names (e.g., F435W)

_FSPS_FILTERS = {
    # ACS/WFC
    "F435W":  "acs_wfc_f435w",
    "F475W":  "acs_wfc_f475w",
    "F555W":  "acs_wfc_f555w",
    "F606W":  "acs_wfc_f606w",
    "F625W":  "acs_wfc_f625w",
    "F775W":  "acs_wfc_f775w",
    "F814W":  "acs_wfc_f814w",
    "F850LP": "acs_wfc_f850lp",
    # WFC3/IR
    "F098M": "wfc3_ir_f098m",
    "F105W": "wfc3_ir_f105w",
    "F110W": "wfc3_ir_f110w",
    "F125W": "wfc3_ir_f125w",
    "F140W": "wfc3_ir_f140w",
    "F160W": "wfc3_ir_f160w",
    # WFC3/UVIS
    "F225W": "wfc3_uvis_f225w",
    "F275W": "wfc3_uvis_f275w",
    "F336W": "wfc3_uvis_f336w",
    "F390W": "wfc3_uvis_f390w",
    # JWST/NIRCam
    "F070W": "jwst_f070w",
    "F090W": "jwst_f090w",
    "F115W": "jwst_f115w",
    "F150W": "jwst_f150w",
    "F200W": "jwst_f200w",
    "F277W": "jwst_f277w",
    "F356W": "jwst_f356w",
    "F444W": "jwst_f444w",
    # SDSS
    "u": "sdss_u0",
    "g": "sdss_g0",
    "r": "sdss_r0",
    "i": "sdss_i0",
    "z": "sdss_z0",
}

_BAGPIPES_FILTERS = {
    # ACS/WFC (SVO notation)
    "F435W":  "HST/ACS_WFC.F435W",
    "F475W":  "HST/ACS_WFC.F475W",
    "F555W":  "HST/ACS_WFC.F555W",
    "F606W":  "HST/ACS_WFC.F606W",
    "F625W":  "HST/ACS_WFC.F625W",
    "F775W":  "HST/ACS_WFC.F775W",
    "F814W":  "HST/ACS_WFC.F814W",
    "F850LP": "HST/ACS_WFC.F850LP",
    # WFC3/IR
    "F098M": "HST/WFC3_IR.F098M",
    "F105W": "HST/WFC3_IR.F105W",
    "F110W": "HST/WFC3_IR.F110W",
    "F125W": "HST/WFC3_IR.F125W",
    "F140W": "HST/WFC3_IR.F140W",
    "F160W": "HST/WFC3_IR.F160W",
    # WFC3/UVIS
    "F225W": "HST/WFC3_UVIS2.F225W",
    "F275W": "HST/WFC3_UVIS2.F275W",
    "F336W": "HST/WFC3_UVIS2.F336W",
    "F390W": "HST/WFC3_UVIS2.F390W",
    # JWST/NIRCam
    "F070W": "JWST/NIRCam.F070W",
    "F090W": "JWST/NIRCam.F090W",
    "F115W": "JWST/NIRCam.F115W",
    "F150W": "JWST/NIRCam.F150W",
    "F200W": "JWST/NIRCam.F200W",
    "F277W": "JWST/NIRCam.F277W",
    "F356W": "JWST/NIRCam.F356W",
    "F444W": "JWST/NIRCam.F444W",
    # SDSS
    "u": "SLOAN/SDSS.u",
    "g": "SLOAN/SDSS.g",
    "r": "SLOAN/SDSS.r",
    "i": "SLOAN/SDSS.i",
    "z": "SLOAN/SDSS.z",
}

_CIGALE_FILTERS = {
    # ACS/WFC
    "F435W":  "hst.acs.wfc1_f435w",
    "F475W":  "hst.acs.wfc1_f475w",
    "F555W":  "hst.acs.wfc1_f555w",
    "F606W":  "hst.acs.wfc1_f606w",
    "F625W":  "hst.acs.wfc1_f625w",
    "F775W":  "hst.acs.wfc1_f775w",
    "F814W":  "hst.acs.wfc1_f814w",
    "F850LP": "hst.acs.wfc1_f850lp",
    # WFC3/IR
    "F098M": "hst.wfc3.ir_f098m",
    "F105W": "hst.wfc3.ir_f105w",
    "F110W": "hst.wfc3.ir_f110w",
    "F125W": "hst.wfc3.ir_f125w",
    "F140W": "hst.wfc3.ir_f140w",
    "F160W": "hst.wfc3.ir_f160w",
    # WFC3/UVIS
    "F225W": "hst.wfc3.uvis2_f225w",
    "F275W": "hst.wfc3.uvis2_f275w",
    "F336W": "hst.wfc3.uvis2_f336w",
    "F390W": "hst.wfc3.uvis2_f390w",
    # JWST/NIRCam
    "F070W": "jwst.nircam.f070w",
    "F090W": "jwst.nircam.f090w",
    "F115W": "jwst.nircam.f115w",
    "F150W": "jwst.nircam.f150w",
    "F200W": "jwst.nircam.f200w",
    "F277W": "jwst.nircam.f277w",
    "F356W": "jwst.nircam.f356w",
    "F444W": "jwst.nircam.f444w",
    # SDSS
    "u": "sdss.up",
    "g": "sdss.gp",
    "r": "sdss.rp",
    "i": "sdss.ip",
    "z": "sdss.zp",
}

# Prospector uses the same filter names as FSPS (sedpy)
_PROSPECTOR_FILTERS = dict(_FSPS_FILTERS)

# Dense Basis uses the same filter names as FSPS (sedpy)
_DENSE_BASIS_FILTERS = dict(_FSPS_FILTERS)

_BACKEND_FILTER_MAPS = {
    "fsps": _FSPS_FILTERS,
    "bagpipes": _BAGPIPES_FILTERS,
    "cigale": _CIGALE_FILTERS,
    "prospector": _PROSPECTOR_FILTERS,
    "dense_basis": _DENSE_BASIS_FILTERS,
    "analytic": {},  # analytic uses canonical names directly
}


def map_filters(bands, backend_name):
    """Map canonical filter names to a specific backend's convention.

    Parameters
    ----------
    bands : list of str
        Canonical filter names (e.g., ['F435W', 'F814W']).
    backend_name : str
        Backend name (e.g., 'fsps', 'bagpipes').

    Returns
    -------
    list of str
        Backend-specific filter names.

    Raises
    ------
    KeyError
        If a filter is not found in the mapping.
    """
    fmap = _BACKEND_FILTER_MAPS.get(backend_name, {})
    if not fmap:
        return list(bands)

    result = []
    for b in bands:
        if b in fmap:
            result.append(fmap[b])
        else:
            # Try case-insensitive
            matched = False
            for key, val in fmap.items():
                if key.lower() == b.lower():
                    result.append(val)
                    matched = True
                    break
            if not matched:
                raise KeyError(
                    f"Filter '{b}' not found in {backend_name} filter map. "
                    f"Available: {list(fmap.keys())}")
    return result
