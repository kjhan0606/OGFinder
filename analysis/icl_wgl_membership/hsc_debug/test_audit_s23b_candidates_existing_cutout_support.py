import numpy as np
import pytest
from astropy.wcs import WCS

from audit_s23b_candidates_existing_cutout_support import contained_candidate_rows


def test_wcs_containment_includes_interior_and_excludes_outside_and_nonfinite():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [5.0, 5.0]
    wcs.wcs.cdelt = [-1.0 / 3600.0, 1.0 / 3600.0]
    wcs.wcs.crval = [150.0, 2.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    candidates = [
        {"id": "inside", "ra_deg": 150.0, "dec_deg": 2.0, "z": 0.4, "richness": 30.0},
        {"id": "outside", "ra_deg": 150.01, "dec_deg": 2.0, "z": 0.4, "richness": 30.0},
        {"id": "bad", "ra_deg": float("nan"), "dec_deg": 2.0, "z": 0.4, "richness": 30.0},
    ]

    result = contained_candidate_rows(wcs, (9, 9), candidates)

    assert [row["candidate_id"] for row in result] == ["inside"]
    assert result[0]["x_pixel"] == pytest.approx(4.0)
    assert result[0]["y_pixel"] == pytest.approx(4.0)
