import csv
import gzip

import pytest

from prepare_icl_expansion_coverage_plan import attach_members, destination, render_query


def test_destination_is_great_circle_and_wraps_ra():
    ra, dec = destination(359.9, 0.0, 0.01, 90.0)
    assert ra < 360.0 and ra > 0.0
    assert dec == pytest.approx(0.0, abs=1e-10)
    north_ra, north_dec = destination(10.0, 20.0, 0.02, 0.0)
    assert north_dec > 20.0
    assert north_ra == pytest.approx(10.0, abs=1e-10)


def test_coverage_sql_returns_only_patch_metadata():
    sql = render_query([
        {"candidate_id": "17", "probe_id": "center", "ra_deg": 10.5, "dec_deg": -2.0}
    ]).upper()
    assert "PDR3_WIDE.MOSAIC" in sql
    assert "PUBLIC.SKYMAP" in sql
    assert "PDR3_WIDE.PATCH_QA" in sql
    assert "PATCH_CONTAINS" in sql
    assert "BAND = 'HSC-I'" in sql
    assert "OBJECT_ID" not in sql and "PHOTOZ" not in sql and "SHAPE" not in sql
    assert "SELECT *" not in sql


def test_member_attachment_checks_cluster_identity(tmp_path):
    path = tmp_path / "members.dat.gz"
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        stream.write("150.00000 2.00000 31.000 0.50000 149.90000 2.10000 22.0 0.3\n")
        stream.write("150.00000 2.00000 25.000 0.40000 150.01000 2.01000 21.0 1.0\n")
    rows = [{"id": "4", "ra_deg": 150.0, "dec_deg": 2.0, "z": 0.4, "richness": 25.0}]
    members = attach_members(rows, path)
    assert members == {"4": [(150.01, 2.01)]}

    with gzip.open(path, "wt", encoding="utf-8") as stream:
        stream.write("150.00000 2.00000 25.000 0.50000 150.01000 2.01000 21.0 1.0\n")
    with pytest.raises(ValueError, match="lack local member rows"):
        attach_members(rows, path)
