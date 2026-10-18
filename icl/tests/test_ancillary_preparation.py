"""Synthetic checks of the independent SFD lookup and catalogue request contract."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

SOURCE = Path(__file__).resolve().parents[2] / 'analysis/icl_wgl_membership/hsc_debug/prepare_ancillary.py'
spec = importlib.util.spec_from_file_location('prepare_ancillary', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def maps(tmp_path, values=(.1, .2)):
    result = {}
    for pole, sign, value in zip(('ngp', 'sgp'), (1, -1), values):
        w = WCS(naxis=2)
        w.wcs.ctype = ['GLON-ZEA', 'GLAT-ZEA']
        w.wcs.crval = [90, sign*90]
        w.wcs.crpix = [32.5, 32.5]
        w.wcs.cdelt = [-3, 3]
        path = tmp_path / (pole+'.fits')
        fits.writeto(path, np.full((64,64), value), w.to_header())
        result[pole] = path
    return result


def test_galactic_hemisphere_and_unscaled_reddening(tmp_path):
    coords = SkyCoord(l=[10, 240], b=[75,-75], unit='deg', frame='galactic').icrs
    np.testing.assert_allclose(module.query_sfd(coords, maps(tmp_path)), [.1,.2])


def test_invalid_map_values_are_rejected(tmp_path):
    coords = SkyCoord(l=[10], b=[75], unit='deg', frame='galactic').icrs
    with pytest.raises(ValueError, match='Invalid SFD'):
        module.query_sfd(coords, maps(tmp_path, values=(-1,.2)))


def test_queries_use_arcsec_and_preserve_target_id():
    q = module.catalogue_queries(dict(id='123', ra_deg=359.9, dec_deg=-20))
    assert '3.0)' in q['extinction_sql'] and 'a_i' in q['extinction_sql']
    assert '1800.0)' in q['random_candidates_sql']
    assert "'123' AS target_cluster_id" in q['random_candidates_sql']
    assert 'IS FALSE' in q['random_candidates_sql']
    assert module.COEFFICIENTS['i'] != module.COEFFICIENTS['i2']


@pytest.mark.parametrize('ident,ra', [("1'; DROP TABLE x", 10), ('1', float('nan'))])
def test_queries_refuse_unsafe_or_nonfinite_inputs(ident, ra):
    with pytest.raises(ValueError):
        module.catalogue_queries(dict(id=ident, ra_deg=ra, dec_deg=0))
