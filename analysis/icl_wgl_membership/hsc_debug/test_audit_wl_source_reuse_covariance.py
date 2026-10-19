import numpy as np
import pytest

from audit_wl_source_reuse_covariance import source_reuse_summary


def test_source_reuse_reports_repeated_objects_and_endpoint_cells():
    pairs = {
        "object_id": np.array(["a", "a", "a", "b", "b", "c"]),
        "lens_ra_deg": np.array([1.0, 6.0, 1.0, 6.0, 6.0, 1.0]),
        "lens_dec_deg": np.zeros(6),
        "source_ra_deg": np.array([11.0, 11.0, 11.0, 21.0, 21.0, 31.0]),
        "source_dec_deg": np.zeros(6),
    }
    summary = source_reuse_summary(pairs, np.ones(6, dtype=bool), 5.0)

    assert summary["selected_pair_rows"] == 6
    assert summary["unique_source_objects"] == 3
    assert summary["source_objects_reused_across_lens_pairs"] == 2
    assert summary["fraction_unique_sources_reused"] == pytest.approx(2 / 3)
    assert summary["pair_rows_from_reused_sources"] == 5
    assert summary["fraction_pair_rows_from_reused_sources"] == pytest.approx(5 / 6)
    assert summary["maximum_lens_pair_multiplicity_per_source"] == 3
    assert summary["reused_sources_spanning_multiple_lens_endpoint_cells"] == 1
    assert summary["reused_sources_spanning_multiple_source_endpoint_cells"] == 0
    assert summary["source_endpoint_cell_assignment_is_unique_per_object_id"] is True
