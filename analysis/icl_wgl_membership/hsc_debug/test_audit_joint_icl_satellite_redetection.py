import hashlib
import json

import pytest

from audit_joint_icl_satellite_redetection import audit


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_audit_verifies_partial_redetection_batch_and_keeps_science_gate(tmp_path):
    frame_path = tmp_path / "cluster_1.json"
    frame_path.write_text(json.dumps({
        "cluster_id": "1", "filter": "i",
        "schema": "hsc-joint-icl-satellite-redetection-injection-v1",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "injection": {"source_detection_and_masks_rerun": True,
                      "original_fits_modified": False},
        "measurements_by_scenario": {
            "baseline": {}, "icl_only": {}, "satellite_only": {}, "both": {},
        },
        "redetection_by_scenario": {
            "baseline": {"n_bcg_candidates": 1, "bcg_mask_policy": "preserve_unique_bcg"},
            "icl_only": {"n_bcg_candidates": 1, "bcg_mask_policy": "preserve_unique_bcg"},
            "satellite_only": {"n_bcg_candidates": 1, "bcg_mask_policy": "preserve_unique_bcg"},
            "both": {"n_bcg_candidates": 1, "bcg_mask_policy": "preserve_unique_bcg"},
        },
    }))
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps({
        "schema": "hsc-joint-icl-satellite-fixed-mask-batch-v1",
        "status": "partial", "science_measurement": False,
        "injection": {"detection_rerun": True,
                      "masks_and_segmentation_fixed": False,
                      "unresolved_bcg_policy": "mask_all",
                      "poisson_noise_added": False,
                      "psf_convolution_applied": False,
                      "original_fits_modified": False},
        "selection": {"n_selected_cluster_filter_pairs": 2,
                      "selected_pairs": [{"cluster_id": "1", "filter": "i"},
                                         {"cluster_id": "2", "filter": "i2"}]},
        "n_completed": 1, "n_errors": 1,
        "records": [{"cluster_id": "1", "filter": "i", "path": str(frame_path),
                     "sha256": sha(frame_path)}],
        "errors": [{"cluster_id": "2", "filter": "i2", "error_type": "ValueError",
                    "error": "synthetic test gate"}],
    }))

    result = audit(summary_path)

    assert result["status"] == "verified_partial"
    assert result["science_measurement"] is False
    assert result["n_selected"] == 2
    assert result["frame_sha256_verified"] == 1
    assert result["errors_by_type_and_message"] == {"ValueError: synthetic test gate": 1}
    assert result["unresolved_bcg_policy"] == "mask_all"
    assert result["bcg_mask_policy_scenarios"] == {"preserve_unique_bcg": 4}


def test_audit_rejects_missing_selected_pair(tmp_path):
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps({
        "schema": "hsc-joint-icl-satellite-fixed-mask-batch-v1",
        "status": "partial", "science_measurement": False,
        "injection": {"detection_rerun": True,
                      "masks_and_segmentation_fixed": False,
                      "poisson_noise_added": False,
                      "psf_convolution_applied": False,
                      "original_fits_modified": False},
        "selection": {"n_selected_cluster_filter_pairs": 1,
                      "selected_pairs": [{"cluster_id": "1", "filter": "i"}]},
        "n_completed": 0, "n_errors": 0, "records": [], "errors": [],
    }))
    with pytest.raises(ValueError, match="partitioned exactly"):
        audit(summary_path)
