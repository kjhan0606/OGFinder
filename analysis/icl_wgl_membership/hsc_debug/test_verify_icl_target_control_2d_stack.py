import hashlib
import subprocess

import pytest

from verify_icl_target_control_2d_stack import git_version_matches, resolve_mask_factor


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
