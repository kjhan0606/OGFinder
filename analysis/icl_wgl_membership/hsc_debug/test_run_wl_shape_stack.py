import csv
import io

import pytest

from fetch_ancillary_catalogues import catalogue_reader
from prepare_wl_shape_rows import OUTPUT_COLUMNS
from run_wl_shape_stack import _arrays_from_csv


def _csv_reader(raw_masks=(1, 0), z_bins=(1, 2)):
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["# " + OUTPUT_COLUMNS[0], *OUTPUT_COLUMNS[1:]])
    for index, raw_mask in enumerate(raw_masks):
        writer.writerow([
            "lens-1", 150.0, 2.0, 0.3, str(9007199254740993 + index),
            150.01 + 0.001 * index, 2.0, 0.1, -0.04, 0.2, 0.3, 2.0,
            0.01, 0.001, -0.001, 24.7, 0.5, 1.0, 1.0, 0.0,
            raw_mask, z_bins[index],
        ])
    buffer.seek(0)
    return catalogue_reader(buffer)


def test_csv_adapter_preserves_ids_bins_and_explicit_mask_polarity():
    arrays, observed, n_rows = _arrays_from_csv(
        _csv_reader(), masked_value=1, lens_counts={"lens-1": 2},
    )
    assert n_rows == 2
    assert observed == {"lens-1": 2}
    assert arrays["object_id"].tolist() == ["9007199254740993", "9007199254740994"]
    assert arrays["source_zbin"].tolist() == [1, 2]
    assert arrays["b_mode_pass"].tolist() == [False, True]


def test_csv_adapter_rejects_nonbinary_hsc_mask_values():
    with pytest.raises(ValueError, match="not binary"):
        _arrays_from_csv(_csv_reader(raw_masks=(0, 2)), masked_value=1,
                         lens_counts={"lens-1": 2})


def test_csv_adapter_rejects_noninteger_redshift_bins():
    with pytest.raises(ValueError, match="integer values"):
        _arrays_from_csv(_csv_reader(z_bins=(1.5, 2)), masked_value=1,
                         lens_counts={"lens-1": 2})
