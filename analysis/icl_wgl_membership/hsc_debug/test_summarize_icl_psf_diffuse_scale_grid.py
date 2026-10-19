import pytest

from summarize_icl_psf_diffuse_scale_grid import summarize_grid


def make_run(scale, n=2):
    records = {}
    for filt in ("i", "i2"):
        for index in range(n):
            key = (str(index + 1), filt)
            raw = 0.2 + 0.01 * index
            retention = 1.0 - 0.001 * scale
            records[key] = {"transfer_response": {
                "icl_response_to_diffuse_unit_injection": raw,
                "icl_plane_corrected_response_to_diffuse_unit_injection": raw * retention,
            }}
    return {
        "scale_kpc": scale,
        "directory": f"/run/{scale}",
        "summary_sha256": f"summary-{scale}",
        "audit_path": f"/run/{scale}/audit.json",
        "audit_sha256": f"audit-{scale}",
        "records": records,
        "code_sha256": {"code.py": "code-hash"},
        "inputs": {"manifest": "input-hash"},
    }


def test_scale_grid_pairs_identical_frames_and_reports_change_from_anchor():
    runs = {scale: make_run(scale) for scale in (100.0, 200.0, 300.0, 500.0)}

    result = summarize_grid(runs)

    assert result["science_measurement"] is False
    assert result["scale_lengths_kpc"] == [100.0, 200.0, 300.0, 500.0]
    assert result["filters"]["i"]["n_common_frames_for_paired_retention"] == 2
    assert result["filters"]["i"]["paired_retention_change_from_100kpc"]["300"]["median"] == pytest.approx(-0.2)


def test_scale_grid_rejects_unequal_frame_support():
    runs = {scale: make_run(scale) for scale in (100.0, 200.0, 500.0)}
    del runs[200.0]["records"][("1", "i")]

    with pytest.raises(ValueError, match="identical cluster/filter support"):
        summarize_grid(runs)


def test_scale_grid_requires_both_preexisting_anchor_scales():
    with pytest.raises(ValueError, match="include the existing 100 and 500"):
        summarize_grid({200.0: make_run(200.0), 300.0: make_run(300.0)})
