import csv
import hashlib
import json

import pytest

from analyze_joint_icl_satellite_batch import audit, spatial_bootstrap


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_spatial_bootstrap_is_seeded_and_keeps_blocks_together():
    first = spatial_bootstrap([1., 2., 10.], ["a", "a", "b"], 200, 19)
    second = spatial_bootstrap([1., 2., 10.], ["a", "a", "b"], 200, 19)
    assert first == second
    assert first["n_clusters"] == 3
    assert first["n_spatial_blocks_5deg"] == 2
    assert first["largest_block_clusters"] == 2


def test_audit_checks_outputs_and_builds_filter_spatial_summary(tmp_path):
    icl = tmp_path / "icl.tsv"
    sat = tmp_path / "sat.tsv"
    write_tsv(icl, ["id", "filter", "coverage_200_300"], [
        {"id": "1", "filter": "i", "coverage_200_300": ".5"},
        {"id": "2", "filter": "i", "coverage_200_300": ".6"},
    ])
    write_tsv(sat, ["cluster_id", "filter", "radial_bin_kpc", "target_coverage", "control_coverage"], [
        {"cluster_id": "1", "filter": "i", "radial_bin_kpc": "50-75", "target_coverage": ".7", "control_coverage": ".7"},
        {"cluster_id": "2", "filter": "i", "radial_bin_kpc": "50-75", "target_coverage": ".8", "control_coverage": ".8"},
    ])
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"clusters": [
        {"id": "1", "ra_deg": 7.0, "dec_deg": 1.0},
        {"id": "2", "ra_deg": 8.0, "dec_deg": 1.0},
    ]}))
    metric_keys = [
        "icl_response_to_diffuse_unit_injection",
        "icl_response_to_all_satellite_envelopes",
        "satellite_response_to_smooth_diffuse",
        "satellite_response_to_satellite_unit_injection",
    ]
    records = []
    for ident, value in (("1", 0.2), ("2", 0.4)):
        response = {k: value for k in metric_keys}
        response.update({"both_minus_sum_residual_icl_njy_arcsec2": 0.,
                         "both_minus_sum_residual_satellite_njy_arcsec2": 0.})
        frame_path = tmp_path / f"cluster_{ident}.json"
        frame_path.write_text(json.dumps({
            "cluster_id": ident, "filter": "i", "science_measurement": False,
            "observational_validation_passed": False, "flux_covariance_available": False,
            "script_sha256": "frame-code-hash", "transfer_response": response,
        }))
        records.append({"cluster_id": ident, "filter": "i", "path": str(frame_path),
                        "sha256": sha(frame_path), "transfer_response": response})
    batch = tmp_path / "batch.json"
    batch.write_text(json.dumps({
        "status": "complete", "n_errors": 0, "n_completed": 2,
        "science_measurement": False, "single_frame_script_sha256": "frame-code-hash",
        "inputs": {str(icl): sha(icl), str(sat): sha(sat)}, "records": records,
    }))
    result = audit(batch, manifest, icl, sat, tmp_path / "audit.json", replicas=200, seed=3)
    assert result["status"] == "verified_exploratory"
    assert result["n_clusters"] == 2
    assert result["n_spatial_blocks_5deg_total"] == 1
    assert result["by_filter"]["i"]["metrics"][metric_keys[0]]["mean"] == pytest.approx(.3)
