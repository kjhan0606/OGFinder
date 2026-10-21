import hashlib
import json

import pytest

from compare_extended_psf_injection_variants import (
    compare_batches,
    summarize_paired_differences,
)
from validate_extended_psf_injection_batch import recompute_transfer_response


def _write_variant(root, *, asymmetric, delta=0.0, filter_name="i2"):
    root.mkdir()
    record = {
        "schema": "astrafek-extended-psf-joint-injection-v1",
        "status": "complete", "cluster_id": "10", "filter": filter_name,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "support": {"n_injected_targets_in_50_500_kpc": 8},
        "injection": {
            "diffuse_amplitude_at_100_njy_arcsec2": 1.0,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": 1.0,
            "diffuse_scale_kpc": 100.0,
            "satellite_scale_kpc": 15.0,
            "n_injected_satellite_centers": 8,
            "unresolved_bcg_policy": "mask_all",
            "HSC_I_to_i2_is_proxy_sensitivity_only": filter_name == "i2",
        },
        "measurements_by_scenario": {
            "baseline": {"icl_raw_200_300_njy_arcsec2": 10.0,
                         "satellite_delta_50_75_njy_arcsec2": 20.0},
            "icl_only": {"icl_raw_200_300_njy_arcsec2": 13.0 + delta,
                         "satellite_delta_50_75_njy_arcsec2": 21.0},
            "satellite_only": {"icl_raw_200_300_njy_arcsec2": 18.0,
                               "satellite_delta_50_75_njy_arcsec2": 25.0},
            "both": {"icl_raw_200_300_njy_arcsec2": 22.0 + delta,
                     "satellite_delta_50_75_njy_arcsec2": 29.0},
        },
    }
    record["transfer_response"] = recompute_transfer_response(record)
    frame_path = root / "cluster_10.json"
    frame_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    symmetric_name = "HSC_PSF_i_asymmetric.fits" if asymmetric else "HSC_PSF_i.fits"
    summary = {
        "status": "complete", "frame_filter": filter_name,
        "public_reference_transfer_mode": ("exploratory_HSC_I_proxy_for_i2" if filter_name == "i2"
                                            else "HSC_I_reference_for_i_coadd"),
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "n_selected": 1, "n_completed": 1, "n_errors": 0,
        "reference": {"name": symmetric_name,
                      "sha256": "b" * 64 if asymmetric else "a" * 64,
                      "psf_approved": False},
        "records": [{"path": str(frame_path),
                      "sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest()}],
    }
    (root / "summary.json").write_text(json.dumps(summary) + "\n", encoding="utf-8")


def test_paired_variant_comparison_requires_same_frames_and_quantifies_delta(tmp_path):
    symmetric = tmp_path / "symmetric"
    asymmetric = tmp_path / "asymmetric"
    _write_variant(symmetric, asymmetric=False)
    _write_variant(asymmetric, asymmetric=True, delta=0.25)

    result = compare_batches(symmetric, asymmetric, expected_frames=1)

    assert result["n_matched_frames"] == 1
    assert result["paired_asymmetric_minus_symmetric"][
        "icl_response_to_diffuse_unit_injection"]["median_asymmetric_minus_symmetric"] == pytest.approx(0.25)
    assert result["science_measurement"] is False


def test_difference_summary_rejects_empty_or_nonfinite_values():
    with pytest.raises(ValueError):
        summarize_paired_differences([])
    with pytest.raises(ValueError):
        summarize_paired_differences([0.0, float("nan")])


def test_paired_variant_comparison_accepts_i_filter(tmp_path):
    symmetric = tmp_path / "symmetric_i"
    asymmetric = tmp_path / "asymmetric_i"
    _write_variant(symmetric, asymmetric=False, filter_name="i")
    _write_variant(asymmetric, asymmetric=True, delta=0.25, filter_name="i")

    result = compare_batches(symmetric, asymmetric, expected_frames=1, filter_name="i")

    assert result["status"] == "complete_paired_psf_model_sensitivity"
    assert result["frame_filter"] == "i"
    assert result["n_matched_frames"] == 1
    assert result["science_measurement"] is False
