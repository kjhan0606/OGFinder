import json

import pytest

from prepare_wl_random_shape_count_plan import prepare, shape_count_query


def _centers():
    return {
        "schema": "astrafek-wl-random-centers-v1",
        "status": "selected_for_shape_count_preflight",
        "science_measurement": False,
        "h0_km_s_mpc": 70.0,
        "omega_m": 0.3,
        "centers": [
            {"random_center_id": "R-001", "matched_lens_id": "L-1",
             "ra_deg": 12.0, "dec_deg": -2.0, "z_lens": 0.31},
            {"random_center_id": "R-002", "matched_lens_id": "L-1",
             "ra_deg": 13.0, "dec_deg": -1.0, "z_lens": 0.31},
            {"random_center_id": "R-003", "matched_lens_id": "L-2",
             "ra_deg": 30.0, "dec_deg": 1.0, "z_lens": 0.42},
        ],
    }


def test_plan_is_count_only_and_preserves_random_to_lens_mapping(tmp_path):
    source = tmp_path / "centers.json"
    source.write_text(json.dumps(_centers()))
    work = tmp_path / "plan"

    manifest = prepare(source, work, batch_size=2)

    assert manifest["status"] == "prepared_not_submitted"
    assert manifest["n_random_centers"] == 3
    assert manifest["n_matched_lenses"] == 2
    assert manifest["n_query_batches"] == 2
    assert manifest["output_columns"] == ["random_center_id", "matched_lens_id", "n_shape_sources"]
    assert manifest["network_used"] is False
    assert manifest["query_submitted"] is False
    first_sql = (work / "random_shape_counts_000.sql").read_text()
    assert "'R-001' AS random_center_id" in first_sql
    assert "'L-1' AS matched_lens_id" in first_sql
    assert "COUNT(DISTINCT b.object_id) AS n_shape_sources" in first_sql
    assert "s19a_wide.weaklensing_hsm_regauss" in first_sql
    assert "e1" not in first_sql.lower()
    assert "photoz" not in first_sql.lower()
    assert manifest["batches"][0]["sha256"]


def test_query_literals_escape_ids_and_numeric_fields_are_exact():
    center = {"random_center_id": "R'1", "matched_lens_id": "L'2",
              "ra_deg": 23.125, "dec_deg": -1.5, "radius_arcsec": 900.0}
    sql = shape_count_query(center)
    assert "'R''1' AS random_center_id" in sql
    assert "'L''2' AS matched_lens_id" in sql
    assert "coneSearch(c.i_coord, 23.125, -1.5, 900)" in sql
    assert "SELECT *" not in sql


@pytest.mark.parametrize("mutate, message", [
    (lambda d: d.update(schema="wrong"), "schema"),
    (lambda d: d.update(status="provisional"), "selected"),
    (lambda d: d["centers"].append(d["centers"][0]), "unique"),
])
def test_preflight_refuses_unselected_or_ambiguous_centers(tmp_path, mutate, message):
    source = tmp_path / "centers.json"
    doc = _centers()
    mutate(doc)
    source.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match=message):
        prepare(source, tmp_path / "plan")


def test_preflight_refuses_nonfinite_radius(tmp_path):
    source = tmp_path / "centers.json"
    source.write_text(json.dumps(_centers()))
    with pytest.raises(ValueError, match="positive"):
        prepare(source, tmp_path / "plan", radius_mpc=float("nan"))


def test_preflight_requires_matched_centers_to_share_target_redshift(tmp_path):
    source = tmp_path / "centers.json"
    document = _centers()
    document["centers"][1]["z_lens"] = 0.32
    source.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="same lens redshift"):
        prepare(source, tmp_path / "plan")
