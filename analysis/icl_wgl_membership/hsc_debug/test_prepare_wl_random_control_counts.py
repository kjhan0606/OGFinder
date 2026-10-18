import json
from pathlib import Path

import pytest

import prepare_wl_random_control_counts as planner


def _lens(index, radius=650.0):
    return {
        "lens_id": str(1000 + index),
        "ra_deg": 130.0 + index * 0.001,
        "dec_deg": -1.0,
        "radius_arcsec": radius,
    }


def test_random_count_sql_returns_only_aggregate_and_excludes_overlapping_apertures():
    sql = planner.random_count_query(_lens(1, radius=650.25))

    assert sql.startswith("SELECT '1001' AS lens_id, COUNT(*) AS n_random_candidates")
    assert "FROM pdr3_wide.random AS r" in sql
    assert "coneSearch(r.coord, 130.001, -1, 1800)" in sql
    assert "AND NOT coneSearch(r.coord, 130.001, -1, 1300.5)" in sql
    assert "r.adjust_density < 0.01" in sql
    assert "r.isprimary IS TRUE" in sql
    assert "r.patch_inner" not in sql
    assert "r.tract_inner" not in sql
    assert "r.i_pixelflags_edge IS FALSE" in sql
    assert "r.i_pixelflags_bright_objectcenter IS FALSE" in sql
    assert "ra," not in sql.lower()
    assert "dec," not in sql.lower()
    assert "object_id" not in sql.lower()
    assert "SELECT r.*" not in sql


def test_random_count_sql_rejects_empty_annulus_and_unsafe_lens_id():
    with pytest.raises(ValueError, match="annulus is empty"):
        planner.random_count_query(_lens(0, radius=900.0))
    unsafe = _lens(1)
    unsafe["lens_id"] = "1'; SELECT ra"
    with pytest.raises(ValueError, match="unsafe characters"):
        planner.random_count_query(unsafe)


def test_plan_writes_32_local_only_count_queries(tmp_path, monkeypatch):
    preflight = tmp_path / "preflight"
    preflight.mkdir()
    (preflight / "preflight_manifest.json").write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(planner, "load_plan", lambda _path: {
        "lenses": [_lens(i) for i in range(800)],
    })
    output = tmp_path / "count-plan"

    manifest = planner.prepare(preflight, output)

    assert manifest["status"] == "prepared_not_submitted"
    assert manifest["network_used"] is False
    assert manifest["query_submitted"] is False
    assert manifest["random_coordinates_returned"] is False
    assert manifest["output_columns"] == ["lens_id", "n_random_candidates"]
    assert manifest["n_query_batches"] == 32
    assert len(list(output.glob("random_counts_*.sql"))) == 32
    first = (output / "random_counts_000.sql").read_text(encoding="utf-8")
    assert first.count("COUNT(*) AS n_random_candidates") == 25
    assert "SELECT *" not in first
    stored = json.loads((output / "random_control_count_plan.json").read_text())
    assert stored["execution_requires_explicit_user_approval"] is True
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        planner.prepare(preflight, output)
