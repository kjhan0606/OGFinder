import csv

import pytest

from prepare_wl_shape_rows import OUTPUT_COLUMNS, lens_query
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


def test_lens_query_uses_catalog_table_for_each_shape_product():
    sql = lens_query({
        "lens_id": "lens-a", "ra_deg": 10, "dec_deg": -2,
        "z_lens": 0.3, "radius_arcsec": 1200,
    })
    assert "b.i_apertureflux_10_mag AS mag_a10" in sql
    assert "a.i_sdssshape_psf_shape11 AS psf_shape11" in sql
    assert "a.i_sdssshape_psf_shape22 AS psf_shape22" in sql
    assert "a.i_sdssshape_psf_shape12 AS psf_shape12" in sql
    assert "c.i_ra AS source_ra_deg" in sql
    assert "a.i_hsmshaperegauss_e1 AS e1" in sql


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
