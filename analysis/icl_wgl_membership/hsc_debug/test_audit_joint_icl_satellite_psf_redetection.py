import pytest

from audit_joint_icl_satellite_psf_redetection import validate_support_accounting


def frame_result(gaps=()):
    scenarios = ("baseline", "icl_only", "satellite_only", "both")
    measurements = {}
    for scenario in scenarios:
        missing = scenario in gaps
        measurements[scenario] = {
            "icl_raw_200_300_njy_arcsec2": 1.0,
            "icl_plane_corrected_200_300_njy_arcsec2": 0.9,
            "satellite_delta_50_75_njy_arcsec2": None if missing else 0.2,
            "satellite_support": {
                "status": ("no_satellite_passed_minimum_support" if missing else "supported"),
                "n_satellites": 0 if missing else 1,
                "n_control_profiles": 0 if missing else 3,
            },
        }
    diffuse_to_satellite_missing = "baseline" in gaps or "icl_only" in gaps
    satellite_to_satellite_missing = "baseline" in gaps or "satellite_only" in gaps
    return {
        "support_status": "partial_satellite_support" if gaps else "complete",
        "satellite_support_gaps": list(gaps),
        "measurements_by_scenario": measurements,
        "transfer_response": {
            "icl_response_to_diffuse_unit_injection": 0.5,
            "icl_plane_corrected_response_to_diffuse_unit_injection": 0.4,
            "icl_response_to_all_satellite_envelopes": 0.2,
            "icl_plane_corrected_response_to_all_satellite_envelopes": 0.1,
            "satellite_response_to_smooth_diffuse": None if diffuse_to_satellite_missing else 0.01,
            "satellite_response_to_satellite_unit_injection": None if satellite_to_satellite_missing else 0.03,
            "both_minus_sum_residual_satellite_njy_arcsec2": None if gaps else 0.0,
        },
    }


def test_support_gap_is_reported_without_dropping_finite_icl_metrics():
    assert validate_support_accounting(frame_result(("baseline",))) == ["baseline"]


def test_all_supported_scenarios_pass_accounting():
    assert validate_support_accounting(frame_result()) == []


def test_support_gap_with_finite_satellite_metric_is_rejected():
    result = frame_result(("baseline",))
    result["measurements_by_scenario"]["baseline"]["satellite_delta_50_75_njy_arcsec2"] = 0.3
    with pytest.raises(ValueError, match="missing satellite support"):
        validate_support_accounting(result)
