import pytest

from run_core_wing_hybrid_grid import (
    configuration_key, configurations, grid_spec,
)


def test_predeclared_hybrid_grid_has_expected_unique_configurations():
    rows = configurations()
    keys = [configuration_key(row) for row in rows]

    assert len(rows) == 12
    assert len(set(keys)) == 12
    assert {row["splice_radius_arcsec"] for row in rows} == {1.8, 2.4, 3.0}
    assert {row["wing_amplitude"] for row in rows} == {0.50, 0.65, 0.80, 1.00}
    assert {row["blend_width_arcsec"] for row in rows} == {0.6}


def test_grid_spec_pins_filter_cluster_and_false_approval_flags():
    result = grid_spec("1056", "i")

    assert result["cluster_id"] == "1056"
    assert result["filter"] == "i"
    assert result["n_configurations"] == 12
    assert result["science_measurement"] is False
    assert result["psf_approved"] is False
    with pytest.raises(ValueError, match="filter_name"):
        grid_spec("1056", "r")
