from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
import pytest

import audit_hsc_psf_picker_target_control_batch as audit


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path: Path):
    plan_dir = tmp_path / "plan"
    download_dir = tmp_path / "download"
    plan_dir.mkdir()
    download_dir.mkdir()
    mapping = [
        {"ra_deg": 10.0, "dec_deg": 2.0,
         "roles": [{"role": "target_cluster_center", "cluster_id": "A"}]},
        {"ra_deg": 11.0, "dec_deg": -3.0,
         "roles": [{"role": "matched_control_center", "control_id": "B"}]},
    ]
    mapping_path = plan_dir / "position_mapping.json"
    mapping_path.write_text(json.dumps(mapping), encoding="utf-8")
    coordinate_path = plan_dir / "psf_picker_coordinates.txt"
    coordinate_path.write_text("#? ra dec\n10.000000000 +2.000000000\n11.000000000 -3.000000000\n",
                               encoding="ascii")
    plan = {
        "schema": audit.EXPECTED_PLAN_SCHEMA,
        "status": "prepared_not_submitted",
        "n_requests_remaining": 2,
        "rerun": "pdr3_wide",
        "type": "coadd",
        "filter": "HSC-I",
        "centered": True,
        "coordinate_list_sha256": _sha(coordinate_path),
        "position_mapping_sha256": _sha(mapping_path),
    }
    (plan_dir / "psf_picker_plan.json").write_text(json.dumps(plan), encoding="utf-8")
    psf = np.zeros((3, 3), dtype=np.float64)
    psf[1, 1] = 1.0
    names = [
        "psf_123_1_2_10.000000000_+2.000000000_HSC-I.fits",
        "psf_124_2_1_11.000000000_-3.000000000_HSC-I.fits",
    ]
    for name in names:
        fits.PrimaryHDU(psf).writeto(download_dir / name)
    return plan_dir, download_dir, names


def test_audits_complete_batch_and_records_hashes(tmp_path):
    plan_dir, download_dir, names = _fixture(tmp_path)
    output = tmp_path / "audit.json"
    result = audit.audit_batch(plan_dir, download_dir, output, expected_requests=2)
    assert result["status"] == "verified"
    assert result["n_approved"] == result["n_returned"] == 2
    assert result["n_matched_coordinates"] == 2
    assert {row["file"] for row in result["files"]} == set(names)
    assert all(row["sha256"] for row in result["files"])
    assert result["science_measurement"] is False
    assert result["psf_approved"] is False


def test_rejects_missing_position(tmp_path):
    plan_dir, download_dir, names = _fixture(tmp_path)
    (download_dir / names[0]).unlink()
    with pytest.raises(ValueError, match="expected 2 FITS files"):
        audit.audit_batch(plan_dir, download_dir, tmp_path / "audit.json", expected_requests=2)


def test_refuses_to_overwrite_audit(tmp_path):
    plan_dir, download_dir, _ = _fixture(tmp_path)
    output = tmp_path / "audit.json"
    output.write_text("preserve", encoding="utf-8")
    with pytest.raises(FileExistsError):
        audit.audit_batch(plan_dir, download_dir, output, expected_requests=2)
    assert output.read_text(encoding="utf-8") == "preserve"
