import hashlib
import subprocess

import pytest

from verify_icl_target_control_2d_stack import (candidate_audit_pairs, git_version_matches,
                                               resolve_mask_factor)


@pytest.mark.parametrize("factor, extension", [
    (2.5, "MASK_K2P5"),
    (3.5, "MASK_K3P5"),
    (5.0, "MASK_K5P0"),
])
def test_resolves_mask_factor_from_stack_definition(factor, extension):
    summary = {"stack_definition": {
        "mask": f"target {extension}; matched-control SEP/Kron {factor:g}; native masks retained"
    }}
    assert resolve_mask_factor(summary) == factor


def test_prefers_and_validates_explicit_mask_factor():
    summary = {"mask_factor": 3.5, "stack_definition": {"mask": "legacy description"}}
    assert resolve_mask_factor(summary) == 3.5


def test_rejects_disagreeing_target_and_control_factors():
    summary = {"stack_definition": {
        "mask": "target MASK_K2P5; matched-control SEP/Kron 3.5"
    }}
    with pytest.raises(ValueError, match="disagree"):
        resolve_mask_factor(summary)


def test_rejects_missing_or_unsupported_mask_factor():
    with pytest.raises(ValueError, match="no numeric mask factor"):
        resolve_mask_factor({"stack_definition": {"mask": "unknown"}})
    with pytest.raises(ValueError, match="unsupported"):
        resolve_mask_factor({"mask_factor": 4.0})


def test_matches_archived_source_from_git_history(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.invalid"], check=True)
    source = repo / "tool.py"
    original = b"print('run version')\n"
    source.write_bytes(original)
    subprocess.run(["git", "-C", str(repo), "add", "tool.py"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "archive run code"], check=True)
    source.write_text("print('current version')\n")

    assert git_version_matches(source, hashlib.sha256(original).hexdigest(), repo)
    assert not git_version_matches(source, "0" * 64, repo)


def test_independently_reconstructs_candidate_control_pairs():
    audit = {
        "schema": "astrafek-s23b-candidate-existing-cutout-support-v1",
        "status": "header_wcs_audit_complete", "science_measurement": False,
        "hsc_query_submitted": False, "data_downloaded": False,
        "image_pixels_read": False, "n_candidates": 1333,
        "candidate_table_sha256": "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c",
        "target": {"n_images": 800, "n_images_with_candidate_centers": 0,
                   "n_unique_candidates_inside": 0, "matches": []},
        "matched_controls": {
            "inputs_sha256": {"manifest": "m" * 64, "download_status": "s" * 64},
            "n_images": 798, "n_images_with_candidate_centers": 1,
            "matches": [{"target_cluster_id": "t1", "frame_id": "c1", "filter": "i",
                         "n_candidate_centers_inside": 1,
                         "candidates": [{"candidate_id": "42"}]}],
            "n_unique_candidates_inside": 1,
        },
    }
    controls = {"t1": {"id": "c1", "coadd_header_filter": "i"}}
    assert candidate_audit_pairs(audit, controls, "m" * 64, "s" * 64) == {
        "t1": ("c1", "i", ["42"])
    }
    audit["matched_controls"]["inputs_sha256"]["manifest"] = "x" * 64
    with pytest.raises(ValueError, match="checksums"):
        candidate_audit_pairs(audit, controls, "m" * 64, "s" * 64)
