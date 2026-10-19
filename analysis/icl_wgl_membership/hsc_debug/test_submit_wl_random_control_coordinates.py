import csv

import pytest

from submit_wl_random_control_coordinates import coordinate_query, read_coordinates


def test_coordinate_query_preserves_candidate_filters_without_catalog_rows():
    lens = {"lens_id": "camira-01", "ra_deg": 150.0, "dec_deg": 2.0, "radius_arcsec": 600.0}
    sql = coordinate_query(lens)
    upper = sql.upper()
    assert "SELECT 'CAMIRA-01' AS LENS_ID, R.RA AS RA_DEG, R.DEC AS DEC_DEG" in upper
    assert "FROM PDR3_WIDE.RANDOM AS R" in upper
    assert "R.ADJUST_DENSITY < 0.01" in upper
    assert "R.ISPRIMARY IS TRUE" in upper
    assert "R.I_PIXELFLAGS_EDGE IS FALSE" in upper
    assert "COUNT(" not in upper
    assert "SELECT *" not in upper


def test_read_coordinates_checks_exact_per_lens_counts_and_coordinates(tmp_path):
    path = tmp_path / "coordinates.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["lens_id", "ra_deg", "dec_deg"])
        writer.writerow(["a", "359.5", "-30"])
        writer.writerow(["a", "0.5", "30"])
        writer.writerow(["b", "180", "0"])
    assert read_coordinates(path, {"a": 2, "b": 1}) == 3


def test_read_coordinates_rejects_wrong_counts_and_invalid_coordinates(tmp_path):
    path = tmp_path / "coordinates.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["lens_id", "ra_deg", "dec_deg"])
        writer.writerow(["a", "361", "0"])
    with pytest.raises(ValueError, match="invalid RA"):
        read_coordinates(path, {"a": 1})

    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["lens_id", "ra_deg", "dec_deg"])
        writer.writerow(["a", "20", "0"])
        writer.writerow(["b", "20.1", "0"])
        writer.writerow(["b", "20.2", "0"])
    with pytest.raises(ValueError, match="per-lens candidate counts"):
        read_coordinates(path, {"a": 2, "b": 1})
