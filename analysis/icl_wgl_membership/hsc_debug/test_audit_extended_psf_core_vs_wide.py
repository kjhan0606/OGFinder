import numpy as np
import pytest
import hashlib
import json
from pathlib import Path

from audit_extended_psf_core_vs_wide import load_batch, quantiles


def test_quantiles_are_descriptive_and_mark_zero_overlap():
    result = quantiles([-2.0, -1.0, 0.0, 1.0, 4.0])
    assert result["median"] == pytest.approx(0.0)
    assert result["q16"] == pytest.approx(-1.36)
    assert result["q84"] == pytest.approx(2.08)
    assert result["q16_q84_contains_zero"] is True


def test_quantiles_reject_empty_or_non_finite_input():
    with pytest.raises(ValueError):
        quantiles([])
    with pytest.raises(ValueError):
        quantiles([1.0, np.nan])


def test_load_batch_selects_the_requested_subfilter_and_checks_hashes(tmp_path: Path):
    frame = {
        "cluster_id": "10",
        "filter": "i2",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
    }
    frame_path = tmp_path / "cluster_10.json"
    frame_path.write_text(json.dumps(frame), encoding="utf-8")
    digest = hashlib.sha256(frame_path.read_bytes()).hexdigest()
    summary = {
        "status": "complete",
        "n_completed": 1,
        "n_errors": 0,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "records": [{"path": str(frame_path), "sha256": digest}],
    }
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

    _, selected = load_batch(tmp_path, 1, "test", filter_name="i2")

    assert list(selected) == ["10"]
    _, other_filter = load_batch(tmp_path, 1, "test", filter_name="i")
    assert other_filter == {}
