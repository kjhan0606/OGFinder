import pytest

from match_psf_gaia_batched import (
    angular_separation_arcsec,
    associate,
    build_query,
)


def test_cone_union_query_contains_only_requested_positions_and_gaia_fields():
    targets = [
        {"ra_deg": 12.345678901234, "dec_deg": -4.5},
        {"ra_deg": 120.0, "dec_deg": 30.0},
    ]

    query = build_query(targets)

    assert query.count("CIRCLE(") == 2
    assert "12.345678901234" in query
    assert "120.000000000000,30.000000000000,2.0/3600.0" in query
    assert "tap_upload" not in query
    assert "g.classprob_dsc_combmod_star" in query
    assert "phot_bp_mean_mag" not in query


def test_spherical_separation_wraps_across_zero_ra():
    assert angular_separation_arcsec(359.99999, 0.0, 0.00001, 0.0) == pytest.approx(0.072, abs=1e-5)


def test_local_association_keeps_overlapping_candidate_matches():
    targets = [
        {"ra_deg": 359.99999, "dec_deg": 0.0},
        {"ra_deg": 0.00001, "dec_deg": 0.0},
    ]
    sources = [{"source_id": "123", "ra": "0.0", "dec": "0.0"}]

    matches = associate(targets, sources, offset=250)

    assert [row["candidate_id"] for row in matches] == ["250", "251"]
    assert all(float(row["separation_arcsec"]) < 0.1 for row in matches)
