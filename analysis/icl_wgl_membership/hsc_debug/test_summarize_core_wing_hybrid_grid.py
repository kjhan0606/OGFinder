import json
import hashlib

import pytest

import summarize_core_wing_hybrid_grid as grid_summary
from summarize_core_wing_hybrid_grid import (
    finite_summary, optional_summary, summarize_grid, validate_frame_set,
)


def test_grid_metric_summaries_keep_model_values_descriptive():
    result = finite_summary([0.1, 0.2, 0.3])

    assert result["n_frames"] == 3
    assert result["median"] == pytest.approx(0.2)
    assert result["minimum"] == pytest.approx(0.1)
    assert optional_summary([])["median"] is None


def test_grid_summary_rejects_duplicate_or_mismatched_frame_sets():
    reference = validate_frame_set(["1056", "430"], 2)

    with pytest.raises(ValueError, match="duplicate frame IDs"):
        validate_frame_set(["1056", "1056"], 2)
    with pytest.raises(ValueError, match="identical frame sets"):
        validate_frame_set(["1056", "527"], 2, reference)


def test_grid_summary_rejects_running_or_incomplete_checkpoint(tmp_path):
    (tmp_path / "grid_status.json").write_text(json.dumps({
        "status": "running",
        "spec": {"n_configurations": 1,
                 "science_measurement": False,
                 "psf_approved": False},
        "records": [],
    }), encoding="utf-8")

    with pytest.raises(ValueError, match="incomplete"):
        summarize_grid(tmp_path)


def test_grid_summary_includes_fixed_and_point_source_mask_diagnostics(tmp_path,
                                                                      monkeypatch):
    key = "r2.4_a0.65"
    parameters = {"splice_radius_arcsec": 2.4,
                  "blend_width_arcsec": 0.6,
                  "wing_amplitude": 0.65}
    response = {
        "icl_response_to_diffuse_unit_injection": 0.22,
        "icl_response_to_all_satellite_envelopes": 0.5,
        "satellite_response_to_smooth_diffuse": -0.01,
        "satellite_response_to_satellite_unit_injection": 0.31,
    }
    frame = {
        "cluster_id": "1056",
        "transfer_response": response,
        "fixed_mask_transfer_response": response,
        "redetection_minus_fixed_mask_response": {key: 0.0 for key in response},
        "point_source_wing_profile": {
            "annuli": {
                "50_100_kpc": {
                    "mask_coverage": 0.46,
                    "masked_annulus_flux_fraction": 0.002,
                    "mean_surface_brightness_njy_arcsec2_per_njy": 1e-5,
                },
                "200_300_kpc": {
                    "mask_coverage": 0.55,
                    "masked_annulus_flux_fraction": 0.0008,
                    "mean_surface_brightness_njy_arcsec2_per_njy": 2e-6,
                },
            },
        },
    }
    batch_dir = tmp_path / key
    batch_dir.mkdir()
    frame_path = batch_dir / "cluster_1056.json"
    frame_path.write_text(json.dumps(frame), encoding="utf-8")
    frame_hash = hashlib.sha256(frame_path.read_bytes()).hexdigest()
    batch_summary = {
        "hybrid_parameters": parameters,
        "records": [{"path": str(frame_path), "sha256": frame_hash}],
    }
    summary_path = batch_dir / "summary.json"
    summary_path.write_text(json.dumps(batch_summary), encoding="utf-8")
    summary_hash = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    checkpoint = {
        "status": "complete",
        "spec": {"n_configurations": 1, "n_frames": 1, "filter": "i",
                 "science_measurement": False, "psf_approved": False},
        "records": [{"configuration": key, "status": "complete",
                     "parameters": parameters, "summary_sha256": summary_hash}],
    }
    (tmp_path / "grid_status.json").write_text(json.dumps(checkpoint),
                                                encoding="utf-8")
    monkeypatch.setattr(grid_summary, "audit_batch", lambda *args, **kwargs: {
        "n_frames_validated": 1,
        "max_stored_vs_recomputed_response_error": 0.0,
    })

    result = summarize_grid(tmp_path)

    assert result["n_configurations_validated"] == 1
    assert result["configurations"][0]["frame_ids"] == ["1056"]
    assert result["point_source_model_sensitivity_envelopes"]["200_300_kpc"][
        "masked_annulus_flux_fraction"]["minimum_frame_value"] == pytest.approx(0.0008)
