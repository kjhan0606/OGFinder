import csv

import pytest

from prepare_wl_shape_rows import OUTPUT_COLUMNS
from submit_wl_shape_rows import validate_rows


def _write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def _row(lens_id, object_id):
    return {key: "" for key in OUTPUT_COLUMNS} | {
        "lens_id": lens_id, "object_id": object_id,
    }


def test_row_export_checks_pair_counts_and_accepts_empty_lens(tmp_path):
    path = tmp_path / "rows.csv"
    _write_rows(path, [_row("lens-a", "123"), _row("lens-b", "456")])
    assert validate_rows(path, {"lens-a": 1, "lens-b": 1, "lens-c": 0}) == 2


def test_row_export_rejects_duplicate_pairs_and_count_mismatch(tmp_path):
    path = tmp_path / "rows.csv"
    _write_rows(path, [_row("lens-a", "123"), _row("lens-a", "123")])
    with pytest.raises(ValueError, match=r"duplicate \(lens_id, object_id\)"):
        validate_rows(path, {"lens-a": 2})
    _write_rows(path, [_row("lens-a", "123")])
    with pytest.raises(ValueError, match="disagree with preflight"):
        validate_rows(path, {"lens-a": 2})


def test_row_export_rejects_unrequested_lens(tmp_path):
    path = tmp_path / "rows.csv"
    _write_rows(path, [_row("unexpected", "123")])
    with pytest.raises(ValueError, match="outside the requested batch"):
        validate_rows(path, {"lens-a": 1})
