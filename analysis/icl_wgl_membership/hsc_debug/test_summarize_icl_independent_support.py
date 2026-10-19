import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).with_name("summarize_icl_independent_support.py")
SPEC = importlib.util.spec_from_file_location("summarize_icl_independent_support", SCRIPT)
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_work(root, factor):
    work = root / f"K{factor}"
    work.mkdir()
    radial = work / report.RADIAL_FILE
    fields = [
        "filter", "r_min_kpc", "r_max_kpc", "n_pairs_independent_delta",
        "n_pairs_common_support", "n_pairs_added_vs_common",
        "mean_delta_independent_njy_arcsec2", "pair_sem_delta_independent_njy_arcsec2",
        "mean_independent_delta_on_common_pairs_njy_arcsec2",
        "mean_common_delta_on_same_pairs_njy_arcsec2", "median_target_coverage",
        "median_control_coverage", "median_common_coverage",
    ]
    rows = []
    for band in ("i", "i2"):
        for index in range(20):
            lo, hi = index * 25, (index + 1) * 25
            rows.append({
                "filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                "n_pairs_independent_delta": 10, "n_pairs_common_support": 8,
                "n_pairs_added_vs_common": 2,
                "mean_delta_independent_njy_arcsec2": 2.0,
                "pair_sem_delta_independent_njy_arcsec2": 0.5,
                "mean_independent_delta_on_common_pairs_njy_arcsec2": 2.0,
                "mean_common_delta_on_same_pairs_njy_arcsec2": 1.0,
                "median_target_coverage": 0.5, "median_control_coverage": 0.4,
                "median_common_coverage": 0.3,
            })
    with radial.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    summary = {
        "status": "complete_exploratory", "science_measurement": False,
        "mask_factor": factor, "n_pairs": 783, "n_by_filter": {"i": 302, "i2": 481},
        "outputs": {"radial_profiles_sha256": sha(radial)},
    }
    summary_path = work / "summary.json"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    verification = {
        "status": "verified", "science_measurement": False, "n_pairs": 783,
        "summary_sha256": sha(summary_path),
    }
    (work / "independent_verification.json").write_text(json.dumps(verification), encoding="utf-8")
    return work


def test_report_requires_three_verified_mask_runs_and_formats_same_pair_shift(tmp_path):
    works = [make_work(tmp_path, factor) for factor in (2.5, 3.5, 5.0)]
    result = report.render(works, tmp_path / "report.md")
    assert result.count("| K=") == 36
    assert "100%" in result
    assert "science_measurement=false" in result


def test_report_rejects_tampered_radial_profile(tmp_path):
    work = make_work(tmp_path, 2.5)
    with (work / report.RADIAL_FILE).open("a", encoding="utf-8") as stream:
        stream.write("tamper\n")
    with pytest.raises(ValueError, match="checksum"):
        report.verified_profile(work)
