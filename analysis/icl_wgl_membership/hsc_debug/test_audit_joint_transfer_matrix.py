import numpy as np
import pytest

from audit_joint_transfer_matrix import matrix_diagnostics, summarize_filter, transfer_matrix


def test_transfer_matrix_uses_measurement_rows_and_injected_signal_columns():
    response = {
        "icl_response_to_diffuse_unit_injection": 0.2,
        "icl_response_to_all_satellite_envelopes": 0.6,
        "satellite_response_to_smooth_diffuse": 0.01,
        "satellite_response_to_satellite_unit_injection": 0.4,
    }
    np.testing.assert_allclose(transfer_matrix(response), [[0.2, 0.6], [0.01, 0.4]])
    diag = matrix_diagnostics(response)
    assert diag["determinant"] == pytest.approx(0.074)
    assert diag["condition_number_2"] > 1


def test_transfer_matrix_rejects_nonfinite_responses():
    response = {key: 0.1 for key in (
        "icl_response_to_diffuse_unit_injection",
        "icl_response_to_all_satellite_envelopes",
        "satellite_response_to_smooth_diffuse",
        "satellite_response_to_satellite_unit_injection",
    )}
    response["satellite_response_to_smooth_diffuse"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        transfer_matrix(response)


def test_filter_summary_counts_ill_conditioned_and_negative_determinants():
    rows = [
        {"determinant": -0.1, "condition_number_2": 2.0,
         "matrix": [[0.2, 0.4], [0.0, 0.5]]},
        {"determinant": 0.01, "condition_number_2": 150.0,
         "matrix": [[0.2, 0.6], [0.0, 0.01]]},
    ]
    result = summarize_filter(rows)
    assert result["n_pairs"] == 2
    assert result["n_negative_determinant"] == 1
    assert result["n_condition_gt_10"] == 1
    assert result["n_condition_gt_100"] == 1
    assert result["median_icl_satellite_cross_to_diffuse_response"] == pytest.approx(2.5)
