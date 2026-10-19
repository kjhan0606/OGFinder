import pytest

from audit_wl_bmode_spatial_distribution import (
    EXPECTED_COLUMNS, grid_cell, normalize_shape_header, passes_published_source_cuts,
)


def test_normalize_shape_header_accepts_export_comment_prefix():
    assert normalize_shape_header(["# " + EXPECTED_COLUMNS[0], *EXPECTED_COLUMNS[1:]]) == EXPECTED_COLUMNS


def test_normalize_shape_header_leaves_unmarked_or_other_fields_unchanged():
    assert normalize_shape_header(EXPECTED_COLUMNS) == EXPECTED_COLUMNS
    assert normalize_shape_header(["#other", *EXPECTED_COLUMNS[1:]])[0] == "#other"


def test_grid_cell_handles_ra_wrap_domain_and_dec_edges():
    assert grid_cell(359.99, -90.0) == (359, 0)
    assert grid_cell(0.0, 90.0) == (0, 180)
    assert grid_cell(12.2, -0.2, 0.5) == (24, 179)


@pytest.mark.parametrize("ra, dec", [(-0.01, 0), (360, 0), (1, -90.01), (1, 90.01)])
def test_grid_cell_rejects_invalid_sky_position(ra, dec):
    with pytest.raises(ValueError):
        grid_cell(ra, dec)


def test_published_shape_and_tomographic_cuts_are_applied():
    row = {"e1": "1.2", "e2": "0", "sigma_e": "0.4", "resolution": "0.3",
           "mag_a10": "25.5", "source_zbin": "4"}
    assert passes_published_source_cuts(row)
    row["e1"] = "2.0"
    assert not passes_published_source_cuts(row)
    row["e1"] = "0.1"
    row["source_zbin"] = "0"
    assert not passes_published_source_cuts(row)
