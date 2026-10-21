import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

from audit_empirical_psf_matched_control_support import center_offset_pixel, filters_match, stamp_inside


def tan_wcs(ra, dec):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crpix = [501.0, 501.0]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.cdelt = [-0.168 / 3600, 0.168 / 3600]
    return wcs


def test_same_offset_maps_to_same_pixel_offset_for_translated_tan_frames():
    target_center = SkyCoord(150 * u.deg, 2 * u.deg, frame="icrs")
    control_center = SkyCoord(152 * u.deg, 2.5 * u.deg, frame="icrs")
    target_wcs = tan_wcs(150, 2)
    control_wcs = tan_wcs(152, 2.5)

    x, y, east, north = center_offset_pixel(
        target_wcs, control_wcs, target_center, control_center, (681.25, 738.75)
    )

    assert np.allclose((x - 500, y - 500), (181.25, 238.75), atol=1e-5)
    assert east < 0
    assert north > 0


def test_stamp_inside_requires_full_requested_margin():
    assert stamp_inside((2508, 2509), 200, 200)
    assert stamp_inside((2508, 2509), 2308.9, 2307.9)
    assert not stamp_inside((2508, 2509), 199.99, 200)
    assert not stamp_inside((2508, 2509), 2309, 200)


def test_filter_matching_keeps_i_and_i2_separate():
    target = {"band": "HSC-I", "coadd_header_filter": "i2"}
    assert filters_match(target, {"band": "HSC-I", "coadd_header_filter": "i2"})
    assert not filters_match(target, {"band": "HSC-I", "coadd_header_filter": "i"})
