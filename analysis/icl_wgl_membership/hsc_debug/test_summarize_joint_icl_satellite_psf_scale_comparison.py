import pytest

from summarize_joint_icl_satellite_psf_scale_comparison import compare_runs


def make_run(scale, plane_ratio, satellite_diffuse):
    records = {}
    for cluster, filt in (("10", "i"), ("11", "i"), ("12", "i2")):
        records[(cluster, filt)] = {
            "transfer_response": {
                "icl_response_to_diffuse_unit_injection": 2.0,
                "icl_plane_corrected_response_to_diffuse_unit_injection": 2.0 * plane_ratio,
                "icl_response_to_all_satellite_envelopes": 0.4,
                "icl_plane_corrected_response_to_all_satellite_envelopes": 0.3,
                "satellite_response_to_smooth_diffuse": satellite_diffuse,
                "satellite_response_to_satellite_unit_injection": 0.6,
            },
            "satellite_support_gaps": [],
        }
    return {
        "code_sha256": {"tool.py": "same"},
        "summary": {"inputs": {"psf_plan": "same"},
                    "injection": {"diffuse_scale_kpc": scale}},
        "records": records,
    }


def test_comparison_uses_matched_frames_and_reports_descriptive_retention():
    result = compare_runs(make_run(100, 0.98, 0.01), make_run(500, 0.62, None))
    i = result["filters"]["i"]
    i2 = result["filters"]["i2"]

    assert i["n_matched_frames"] == 2
    assert i2["n_matched_frames"] == 1
    ratio = i["icl_plane_retention_ratio"]
    assert ratio["100kpc"]["median"] == pytest.approx(0.98)
    assert ratio["500kpc"]["median"] == pytest.approx(0.62)
    assert ratio["paired_500_minus_100"]["median"] == pytest.approx(-0.36)
    assert i["transfer_responses"]["satellite_response_to_smooth_diffuse"]["500kpc"]["n"] == 0


def test_comparison_rejects_different_runtime_code_or_inputs():
    reference = make_run(100, 1.0, 0.0)
    extended = make_run(500, 0.6, 0.0)
    extended["summary"]["inputs"]["psf_plan"] = "different"
    with pytest.raises(ValueError, match="input, PSF plan"):
        compare_runs(reference, extended)
