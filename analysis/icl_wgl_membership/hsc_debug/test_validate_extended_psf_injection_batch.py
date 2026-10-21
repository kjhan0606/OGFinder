import hashlib
import json

import pytest

from validate_extended_psf_injection_batch import (
    audit_batch, recompute_transfer_response, validate_hybrid_metadata,
    validate_point_source_profile,
)


def _record():
    scenarios = {
        "baseline": {"icl_raw_200_300_njy_arcsec2": 10.0,
                     "satellite_delta_50_75_njy_arcsec2": 20.0},
        "icl_only": {"icl_raw_200_300_njy_arcsec2": 13.0,
                     "satellite_delta_50_75_njy_arcsec2": 21.0},
        "satellite_only": {"icl_raw_200_300_njy_arcsec2": 18.0,
                           "satellite_delta_50_75_njy_arcsec2": 25.0},
        "both": {"icl_raw_200_300_njy_arcsec2": 22.0,
                 "satellite_delta_50_75_njy_arcsec2": 29.0},
    }
    return {
        "schema": "astrafek-extended-psf-joint-injection-v1",
        "status": "complete", "cluster_id": "10", "filter": "i2",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "injection": {
            "diffuse_amplitude_at_100_njy_arcsec2": 1.0,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": 2.0,
            "HSC_I_to_i2_is_proxy_sensitivity_only": True,
        },
        "measurements_by_scenario": scenarios,
        "transfer_response": {
            "icl_response_to_diffuse_unit_injection": 3.0,
            "icl_response_to_all_satellite_envelopes": 4.0,
            "satellite_response_to_smooth_diffuse": 1.0,
            "satellite_response_to_satellite_unit_injection": 2.5,
            "both_minus_sum_residual_icl_njy_arcsec2": 1.0,
            "both_minus_sum_residual_satellite_njy_arcsec2": 3.0,
        },
    }


def _add_fixed_mask_branch(record):
    fixed = json.loads(json.dumps(record["measurements_by_scenario"]))
    fixed["icl_only"]["icl_raw_200_300_njy_arcsec2"] = 12.0
    fixed["satellite_only"]["icl_raw_200_300_njy_arcsec2"] = 16.0
    fixed["icl_only"]["satellite_delta_50_75_njy_arcsec2"] = 20.5
    record["schema"] = "astrafek-extended-psf-joint-injection-v2"
    record["measurements_by_scenario_fixed_masks"] = fixed
    fixed_response = recompute_transfer_response(
        record, measurement_key="measurements_by_scenario_fixed_masks")
    record["fixed_mask_transfer_response"] = fixed_response
    record["redetection_minus_fixed_mask_response"] = {
        key: record["transfer_response"][key] - value
        for key, value in fixed_response.items()
    }
    return record


def test_recompute_transfer_arithmetic_from_saved_scenarios():
    result = recompute_transfer_response(_record())
    assert result["icl_response_to_diffuse_unit_injection"] == pytest.approx(3.0)
    assert result["icl_response_to_all_satellite_envelopes"] == pytest.approx(4.0)
    assert result["satellite_response_to_smooth_diffuse"] == pytest.approx(1.0)
    assert result["satellite_response_to_satellite_unit_injection"] == pytest.approx(2.5)
    assert result["both_minus_sum_residual_icl_njy_arcsec2"] == pytest.approx(1.0)
    assert result["both_minus_sum_residual_satellite_njy_arcsec2"] == pytest.approx(3.0)


def test_recompute_fixed_mask_and_redetection_selection_effects():
    record = _add_fixed_mask_branch(_record())
    fixed = recompute_transfer_response(
        record, measurement_key="measurements_by_scenario_fixed_masks")
    assert fixed["icl_response_to_diffuse_unit_injection"] == pytest.approx(2.0)
    assert record["redetection_minus_fixed_mask_response"][
        "icl_response_to_diffuse_unit_injection"] == pytest.approx(1.0)


def test_batch_audit_verifies_checksums_filter_flags_and_arithmetic(tmp_path):
    record = _add_fixed_mask_branch(_record())
    record_path = tmp_path / "cluster_10.json"
    record_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    summary = {
        "status": "complete", "frame_filter": "i2",
        "public_reference_transfer_mode": "exploratory_HSC_I_proxy_for_i2",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "n_selected": 1, "n_completed": 1, "n_errors": 0,
        "reference": {"psf_approved": False},
        "records": [{"path": str(record_path),
                      "sha256": hashlib.sha256(record_path.read_bytes()).hexdigest()}],
    }
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

    result = audit_batch(tmp_path, expected_filter="i2", expected_frames=1)

    assert result["status"] == "validated_exploratory_transfer_only"
    assert result["n_frames_validated"] == 1
    assert result["max_stored_vs_recomputed_response_error"] == 0.0
    assert result["frame_response_summaries"]["icl_response_to_diffuse_unit_injection"]["median"] == 3.0
    assert result["fixed_mask_frame_response_summaries"][
        "icl_response_to_diffuse_unit_injection"]["median"] == 2.0
    assert result["redetection_minus_fixed_mask_summaries"][
        "icl_response_to_diffuse_unit_injection"]["median"] == 1.0

    record["filter"] = "i"
    record_path.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        audit_batch(tmp_path, expected_filter="i2", expected_frames=1)


def _hybrid_record():
    record = _add_fixed_mask_branch(_record())
    params = {"splice_radius_arcsec": 2.4,
              "blend_width_arcsec": 0.6,
              "wing_amplitude": 0.65}
    record["hybrid_psf"] = {
        "status": "exploratory_core_wing_hybrid_sensitivity_only",
        "psf_approved": False, "science_measurement": False,
        "splice_radius_arcsec": 2.4, "blend_width_arcsec": 0.6,
        "wing_amplitude_pre_normalization": 0.65,
        "hybrid_signed_sum": 1.0,
        "hybrid_signed_sum_before_normalization": 0.998,
    }
    core_sha = "a" * 64
    satellite_sha = "b" * 64
    satellite_meta = dict(record["hybrid_psf"])
    satellite_meta.update({
        "target_index": 0,
        "local_core_filename": "satellite_psf.fits",
        "local_core_sha256": satellite_sha,
    })
    record["hybrid_satellite_psfs"] = [satellite_meta]
    record["support"] = {"n_injected_targets_in_50_500_kpc": 1}
    record["input_sha256"] = {
        "/data/center_psf.fits": core_sha,
        "/data/satellite_psf.fits": satellite_sha,
    }
    record["injection"]["hybrid_parameters"] = params
    return record, params


def test_hybrid_psf_validator_checks_parameters_normalization_and_flags():
    record, params = _hybrid_record()
    validate_hybrid_metadata(record, params)

    record["hybrid_psf"]["psf_approved"] = True
    with pytest.raises(ValueError, match="unsafe science flags"):
        validate_hybrid_metadata(record, params)


def test_hybrid_psf_validator_rejects_mismatched_amplitude():
    record, params = _hybrid_record()
    record["hybrid_psf"]["wing_amplitude_pre_normalization"] = 0.8
    with pytest.raises(ValueError, match="parameter mismatch"):
        validate_hybrid_metadata(record, params)


def test_hybrid_psf_validator_checks_each_satellite_core_and_checksum():
    record, params = _hybrid_record()
    validate_hybrid_metadata(record, params)

    record["hybrid_satellite_psfs"][0]["local_core_sha256"] = "c" * 64
    with pytest.raises(ValueError, match="checksum does not match"):
        validate_hybrid_metadata(record, params)

    record, params = _hybrid_record()
    record["hybrid_satellite_psfs"] = []
    with pytest.raises(ValueError, match="cover every injected satellite"):
        validate_hybrid_metadata(record, params)


def test_batch_audit_recomputes_hybrid_mode_provenance(tmp_path):
    record, params = _hybrid_record()
    record_path = tmp_path / "cluster_10.json"
    record_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    summary = {
        "status": "complete", "frame_filter": "i2",
        "public_reference_transfer_mode": "exploratory_HSC_I_proxy_for_i2",
        "psf_model_mode": "position_core_public_wing_hybrid_sensitivity",
        "hybrid_parameters": params,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "n_selected": 1, "n_completed": 1, "n_errors": 0,
        "reference": {"psf_approved": False},
        "records": [{"path": str(record_path),
                      "sha256": hashlib.sha256(record_path.read_bytes()).hexdigest()}],
    }
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    result = audit_batch(tmp_path, expected_filter="i2", expected_frames=1)
    assert result["n_frames_validated"] == 1

    record["hybrid_psf"]["hybrid_signed_sum"] = 0.8
    record_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    summary["records"][0]["sha256"] = hashlib.sha256(record_path.read_bytes()).hexdigest()
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(ValueError, match="not signed-unit-normalized"):
        audit_batch(tmp_path, expected_filter="i2", expected_frames=1)


def test_point_source_profile_validator_checks_physical_annuli():
    record = {"injection": {
        "point_source_added_to_redetection_scenarios": False,
        "separate_unit_point_source_flux_njy": 1.0,
    }, "point_source_wing_profile": {
        "source_flux_njy": 1.0, "pixel_scale_arcsec": 0.168,
        "kpc_per_pixel": 5.2,
        "annuli": {
            "50_100_kpc": {"n_valid_pixels": 100, "geometric_pixel_count": 900.,
                           "mask_coverage": 0.8,
                           "mean_surface_brightness_njy_arcsec2_per_njy": 0.2,
                           "masked_annulus_flux_fraction": 0.01},
            "200_300_kpc": {"n_valid_pixels": 200, "geometric_pixel_count": 1200.,
                            "mask_coverage": 0.7,
                            "mean_surface_brightness_njy_arcsec2_per_njy": 0.03,
                            "masked_annulus_flux_fraction": 0.02},
        },
    }}

    validate_point_source_profile(record)

    record["point_source_wing_profile"]["annuli"]["200_300_kpc"]["mask_coverage"] = 1.2
    with pytest.raises(ValueError, match="invalid point-source annulus"):
        validate_point_source_profile(record)
