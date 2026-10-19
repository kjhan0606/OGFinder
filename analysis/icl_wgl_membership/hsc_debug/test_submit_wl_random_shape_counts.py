import csv
import json

import pytest

from prepare_wl_random_shape_count_plan import prepare
from submit_wl_random_shape_counts import read_count_rows, validate_plan


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


def test_submitter_validates_the_exact_generated_plan(tmp_path):
    source = tmp_path / "centers.json"
    source.write_text(json.dumps(_centers()))
    plan_dir = tmp_path / "plan"
    prepare(source, plan_dir, batch_size=2)

    manifest, centers, batches, plan_sha = validate_plan(plan_dir)

    assert manifest["status"] == "prepared_not_submitted"
    assert len(centers) == 3
    assert len(batches) == 2
    assert len(plan_sha) == 64


def test_submitter_rejects_changed_sql_or_source_manifest(tmp_path):
    source = tmp_path / "centers.json"
    source.write_text(json.dumps(_centers()))
    plan_dir = tmp_path / "plan"
    prepare(source, plan_dir, batch_size=2)
    sql_path = plan_dir / "random_shape_counts_000.sql"
    sql_path.write_text(sql_path.read_text() + "-- changed\n")
    with pytest.raises(ValueError, match="SQL changed"):
        validate_plan(plan_dir)


def test_count_result_validation_preserves_match_and_count_fields(tmp_path):
    path = tmp_path / "counts.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["# random_center_id", "matched_lens_id", "n_shape_sources"])
        writer.writerow(["R-1", "L-1", "42"])
        writer.writerow(["R-2", "L-1", "0"])
    expected = [
        {"random_center_id": "R-1", "matched_lens_id": "L-1"},
        {"random_center_id": "R-2", "matched_lens_id": "L-1"},
    ]
    rows = read_count_rows(path, expected)
    assert rows == {
        "R-1": {"matched_lens_id": "L-1", "n_shape_sources": 42},
        "R-2": {"matched_lens_id": "L-1", "n_shape_sources": 0},
    }


@pytest.mark.parametrize("records, message", [
    ([{"random_center_id": "R-1", "matched_lens_id": "L-x", "n_shape_sources": "1"}], "association"),
    ([{"random_center_id": "R-1", "matched_lens_id": "L-1", "n_shape_sources": "-1"}], "invalid"),
    ([{"random_center_id": "R-1", "matched_lens_id": "L-1", "n_shape_sources": "1"}], "cover"),
])
def test_count_result_validation_rejects_bad_rows(tmp_path, records, message):
    path = tmp_path / "counts.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, ["random_center_id", "matched_lens_id", "n_shape_sources"])
        writer.writeheader()
        writer.writerows(records)
    expected = [
        {"random_center_id": "R-1", "matched_lens_id": "L-1"},
        {"random_center_id": "R-2", "matched_lens_id": "L-1"},
    ]
    with pytest.raises(ValueError, match=message):
        read_count_rows(path, expected)
