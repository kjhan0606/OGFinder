import csv
import hashlib

import pytest

from screen_wl_random_completed_batches import load_verified_completed_batches


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("# lens_id,ra_deg,dec_deg\n")
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_loads_only_complete_verified_coordinate_batches(tmp_path):
    rows0 = [["L1", "10.0", "-1.0"], ["L1", "11.0", "-1.0"]]
    digest = write_csv(tmp_path / "random_coordinates_000.csv", rows0)
    batches = [{"n_lenses": 1, "expected_rows": 2}, {"n_lenses": 1, "expected_rows": 1}]
    lens_ids = ["L1", "L2"]
    counts = {"L1": 2, "L2": 1}
    checkpoint = {"status": "running", "plan_sha256": "plan", "jobs": {
        "random_coordinates_000": {"status": "complete", "job_id": "100", "rows": 2,
                                    "sha256": digest},
        "random_coordinates_001": {"status": "submitted", "job_id": "101"},
    }}

    coordinates, completed = load_verified_completed_batches(
        tmp_path, batches, lens_ids, counts, 1, checkpoint, "plan")

    assert coordinates == [{"lens_id": "L1", "ra_deg": 10.0, "dec_deg": -1.0},
                           {"lens_id": "L1", "ra_deg": 11.0, "dec_deg": -1.0}]
    assert completed == [{"name": "random_coordinates_000", "job_id": "100",
                          "lens_ids": ["L1"],
                          "n_lenses": 1, "n_coordinates": 2, "sha256": digest}]


def test_rejects_changed_hash_and_mismatched_per_lens_counts(tmp_path):
    path = tmp_path / "random_coordinates_000.csv"
    digest = write_csv(path, [["L1", "10.0", "-1.0"]])
    batch = [{"n_lenses": 1, "expected_rows": 1}]
    checkpoint = {"status": "running", "plan_sha256": "plan", "jobs": {
        "random_coordinates_000": {"status": "complete", "job_id": "100", "rows": 1,
                                    "sha256": "0" * 64}}}
    with pytest.raises(ValueError, match="missing or changed"):
        load_verified_completed_batches(tmp_path, batch, ["L1"], {"L1": 1}, 1,
                                        checkpoint, "plan")

    checkpoint["jobs"]["random_coordinates_000"]["sha256"] = digest
    with pytest.raises(ValueError, match="per-lens candidate counts"):
        load_verified_completed_batches(tmp_path, batch, ["L1"], {"L1": 2}, 1,
                                        checkpoint, "plan")


def test_rejects_a_different_plan_and_no_completed_batches(tmp_path):
    checkpoint = {"status": "running", "plan_sha256": "old", "jobs": {}}
    with pytest.raises(ValueError, match="different plan"):
        load_verified_completed_batches(tmp_path, [], [], {}, 1, checkpoint, "new")

    checkpoint["plan_sha256"] = "plan"
    with pytest.raises(ValueError, match="no verified completed"):
        load_verified_completed_batches(tmp_path, [{"n_lenses": 1, "expected_rows": 1}],
                                        ["L1"], {"L1": 1}, 1, checkpoint, "plan")
