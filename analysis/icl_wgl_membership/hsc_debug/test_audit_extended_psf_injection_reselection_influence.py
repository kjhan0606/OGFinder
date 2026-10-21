import hashlib
import json

import pytest

from audit_extended_psf_injection_reselection_influence import audit_batch


def _make_batch(tmp_path, n=3):
    records = []
    for index in range(n):
        ident = str(100 + index)
        red = {"icl_response_to_diffuse_unit_injection": 0.2 + index,
               "icl_response_to_all_satellite_envelopes": 0.3 + index,
               "satellite_response_to_smooth_diffuse": 0.01 + index,
               "satellite_response_to_satellite_unit_injection": 0.4 + index}
        fixed = {key: value - 0.05 for key, value in red.items()}
        delta = {key: red[key] - fixed[key] for key in red}
        measure = {name: {"satellite_support": {"n_satellites": 4, "n_control_profiles": 5},
                          "icl_coverage_200_300": 0.6}
                   for name in ("baseline", "icl_only", "satellite_only", "both")}
        frame = {"status": "complete", "science_measurement": False, "filter": "i",
                 "transfer_response": red, "fixed_mask_transfer_response": fixed,
                 "redetection_minus_fixed_mask_response": delta,
                 "measurements_by_scenario": measure,
                 "measurements_by_scenario_fixed_masks": measure}
        path = tmp_path / f"cluster_{ident}.json"
        path.write_text(json.dumps(frame), encoding="utf-8")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        records.append({"cluster_id": ident, "filter": "i", "path": str(path),
                        "sha256": digest, "transfer_response": red,
                        "fixed_mask_transfer_response": fixed,
                        "redetection_minus_fixed_mask_response": delta})
    (tmp_path / "summary.json").write_text(json.dumps({
        "status": "complete", "science_measurement": False,
        "frame_filter": "i", "records": records,
    }), encoding="utf-8")


def test_audit_recomputes_frame_support_and_leave_one_out(tmp_path):
    _make_batch(tmp_path)
    result = audit_batch(tmp_path, expected_frames=3)
    assert result["status"] == "validated_exploratory_sensitivity_only"
    assert result["metrics"]["icl_response_to_diffuse_unit_injection"]["redetected"]["median"] == pytest.approx(1.2)
    assert result["metrics"]["icl_response_to_diffuse_unit_injection"]["fixed_mask"]["median"] == pytest.approx(1.15)
    assert result["frames"][0]["baseline_support"]["n_satellites"] == 4


def test_audit_rejects_modified_frame_or_bad_reselection_arithmetic(tmp_path):
    _make_batch(tmp_path)
    p = tmp_path / "cluster_100.json"
    data = json.loads(p.read_text())
    data["redetection_minus_fixed_mask_response"]["icl_response_to_diffuse_unit_injection"] = 9
    p.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="hash mismatch"):
        audit_batch(tmp_path, expected_frames=3)
