#!/usr/bin/env python3
"""Compare the offline S23B candidate pool with ICL sample-size scenarios.

This is a catalog-only projection. It submits no HSC query and treats the
current joint-support yield as an explicit, uncertain extrapolation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data")
DEFAULT_CATALOG_AUDIT = DEFAULT_ROOT / "catalogue_s23b_camira_audit_2026-10-20_02"
DEFAULT_POWER_PLAN = DEFAULT_ROOT / "sample_800_2026-10-13/icl_sample_size_plan_2026-10-21_01/sample_size_plan.json"
DEFAULT_OUTPUT = DEFAULT_ROOT / "sample_800_2026-10-13/icl_s23b_expansion_feasibility_2026-10-21_01/summary.json"
EXPECTED_CATALOG_SUMMARY_SHA256 = "f98e9bdb0b7f5cfaf217132cafa6c8a0cf9187a0ec448cc9ee1c50826c800e20"
EXPECTED_CANDIDATE_SHA256 = "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c"
EXPECTED_POWER_PLAN_SHA256 = "df3d6d06adecaa310ab94df94bbebd1441262087520c16385927f285846f0b88"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(catalog_audit_dir: Path, power_plan_path: Path) -> dict:
    catalog_audit_dir = catalog_audit_dir.resolve()
    catalog_summary_path = catalog_audit_dir / "summary.json"
    candidate_path = catalog_audit_dir / "greedy_candidate_additions.tsv"
    for path in (catalog_summary_path, candidate_path, power_plan_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(catalog_summary_path) != EXPECTED_CATALOG_SUMMARY_SHA256:
        raise ValueError("catalog audit summary is not the reviewed pinned result")
    if sha256(candidate_path) != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("S23B candidate table differs from the reviewed result")
    if sha256(power_plan_path) != EXPECTED_POWER_PLAN_SHA256:
        raise ValueError("sample-size plan differs from the reviewed result")

    catalog = json.loads(catalog_summary_path.read_text(encoding="utf-8"))
    power = json.loads(power_plan_path.read_text(encoding="utf-8"))
    if (catalog.get("schema") != "astrafek-s23b-camira-catalog-audit-v1"
            or catalog.get("status") != "complete_catalog_only_feasibility_audit"
            or catalog.get("science_measurement") is not False
            or catalog.get("query_submitted") is not False):
        raise ValueError("S23B input is not a completed offline catalog audit")
    if (power.get("schema") != "astrafek-icl-satellite-sample-size-plan-v1"
            or power.get("status") != "exploratory_design_sensitivity"
            or power.get("science_measurement") is not False):
        raise ValueError("sample-size input is not an exploratory design plan")

    with candidate_path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "id" not in reader.fieldnames:
            raise ValueError("candidate table has no cluster ID column")
        candidate_ids = [row["id"] for row in reader]
    if len(set(candidate_ids)) != len(candidate_ids):
        raise ValueError("candidate table contains duplicate cluster IDs")
    expected_candidates = int(catalog["candidate_pool"]["greedy_additions_after_existing_800_and_6arcmin_separation"])
    if len(candidate_ids) != expected_candidates:
        raise ValueError("candidate table row count disagrees with the catalog audit")

    parent_n = int(power["input"]["current_parent_lenses"])
    yields = {key: float(value) for key, value in power["input"]["observed_support_yield"].items()}
    if set(yields) != {"i", "i2"}:
        raise ValueError("unexpected filter support-yield keys")
    projected_new = {band: len(candidate_ids) * rate for band, rate in yields.items()}
    current_support = {key: int(value) for key, value in power["input"]["current_joint_support_rows"].items()}
    projected_total_support = {band: current_support[band] + projected_new[band] for band in yields}
    scenarios = []
    for effect in (0.20, 0.293):
        for deff in (1.0, 1.5):
            row = next((item for item in power["scenarios"]
                        if item["planning_partial_r"] == effect
                        and item["target_power"] == 0.80
                        and item["multiplicity_scenario"] == "bonferroni_168"
                        and item["design_effect"] == deff), None)
            if row is None:
                raise ValueError(f"missing sample-size scenario: r={effect}, design_effect={deff}")
            required = int(row["inflated_usable_rows"])
            scenarios.append({
                "assumed_partial_r": effect,
                "design_effect": deff,
                "required_i2_usable_rows": required,
                "projected_i2_usable_rows_after_expansion": projected_total_support["i2"],
                "projected_fraction_of_required": projected_total_support["i2"] / required,
                "projected_i2_parent_lenses_at_constant_yield": int(row["projected_parent_lenses_i2"]),
                "candidate_pool_alone_meets_required_rows": projected_total_support["i2"] >= required,
            })
    return {
        "schema": "astrafek-icl-s23b-expansion-feasibility-v1",
        "status": "offline_catalog_projection_only",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "source_shape_photoz_or_image_data_received": False,
        "current_parent_lenses": parent_n,
        "offline_new_candidate_upper_bound": len(candidate_ids),
        "catalog_only_parent_upper_bound": parent_n + len(candidate_ids),
        "current_joint_support_rows": current_support,
        "observed_support_yield_assumed_for_new_candidates": yields,
        "expected_new_joint_support_rows_if_yield_transfers": projected_new,
        "projected_total_support_rows_if_yield_transfers": projected_total_support,
        "sample_size_scenarios": scenarios,
        "catalog_provenance": {
            str(catalog_summary_path): sha256(catalog_summary_path),
            str(candidate_path): sha256(candidate_path),
        },
        "sample_size_plan_sha256": sha256(power_plan_path),
        "limitations": [
            "The 1,333 additions are a catalog-only upper bound; PDR3 image, mask, PSF, depth, and matched-control support are unknown.",
            "The 3.25% i and 8.5% i2 joint-support yields are extrapolated from the selected 800-cluster sample, not measured for candidates.",
            "The r=0.20 and r=0.293 values are planning assumptions; 0.293 is a post-selection exploratory estimate.",
            "This projection does not address PSF wings, sky transfer, completeness, or spatial covariance and cannot establish a detection.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-audit-dir", type=Path, default=DEFAULT_CATALOG_AUDIT)
    parser.add_argument("--power-plan", type=Path, default=DEFAULT_POWER_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    result = run(args.catalog_audit_dir, args.power_plan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "catalog_only_parent_upper_bound": result["catalog_only_parent_upper_bound"],
        "projected_total_support_rows_if_yield_transfers": result["projected_total_support_rows_if_yield_transfers"],
        "output": str(args.output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
