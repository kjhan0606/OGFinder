import json

import pytest

from summarize_core_wing_hybrid_grid import summarize_grid


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
