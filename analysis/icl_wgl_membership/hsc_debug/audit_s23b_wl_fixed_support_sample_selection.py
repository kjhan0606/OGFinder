#!/usr/bin/env python3
"""Describe how the all-annulus S23B weak-lensing support cut selects lenses.

This offline audit reads the approved 80-lens manifest and the already
verified 74-lens per-annulus support table. It makes no network request,
does not read source-row CSVs or images, and is not a shear measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = DATA / "wl_s23b_wl_common_support_plan_2026-11-10_01"
SUPPORT_DIR = DATA / "wl_s23b_matched_radial_support_2026-11-11_02"
QA_DIR = DATA / "wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01"
DEFAULT_OUTPUT = DATA / "wl_s23b_wl_fixed_support_selection_2026-11-11_01"

EXPECTED_SHA256 = {
    "selected_lenses.tsv": "93ae6a2331928d8c593286015dc5330858fe5c63217e64297ef1a81daf4e59af",
    "support_summary.json": "c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234",
    "per_lens_radial_support.tsv": "6744ece545781faffebef65fc7a0143284da3a234425b6bed540bd1e9d19625c",
    "qa_summary.json": "eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445",
}
EXPECTED_GLOBAL_LENS_IDS_SHA256 = "78e35dcab6bc000160d391ae3b16516bc342ee9f74a550a6106bb256afef7a08"
EXPECTED_FIXED52_LENS_IDS_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"
EXPECTED_APPROVED_ROW_PLAN_SHA256 = "7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc"
ATTRIBUTES = (
    "z_lens", "richness", "n_raw_shape_sources",
    "eligible_control_candidates_6mpc", "eligible_control_candidates_8mpc",
    "old_union_sources", "old_union_area_fraction_4096",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def id_sha256(values: set[str]) -> str:
    return hashlib.sha256("\n".join(sorted(values)).encode()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def distribution(values: list[float]) -> dict:
    if not values:
        return {key: None for key in ("mean", "median", "q16", "q84", "min", "max")}
    return {
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "q16": quantile(values, 0.16),
        "q84": quantile(values, 0.84),
        "min": min(values),
        "max": max(values),
    }


def standardized_mean_difference(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2 or len(right) < 2:
        return None
    pooled = math.sqrt((statistics.variance(left) + statistics.variance(right)) / 2.0)
    if pooled == 0:
        return None
    return (statistics.mean(left) - statistics.mean(right)) / pooled


def load_inputs(plan_dir: Path, support_dir: Path) -> tuple[dict, dict[str, dict], list[dict]]:
    paths = {
        "selected_lenses.tsv": plan_dir / "selected_lenses.tsv",
        "support_summary.json": support_dir / "summary.json",
        "per_lens_radial_support.tsv": support_dir / "per_lens_radial_support.tsv",
    }
    for name, path in paths.items():
        observed = sha256(path)
        if observed != EXPECTED_SHA256[name]:
            raise ValueError(f"pinned input hash changed for {name}: {observed}")

    summary = json.loads(paths["support_summary.json"].read_text(encoding="utf-8"))
    if (summary.get("status") != "complete_local_radial_matched_support_audit"
            or summary.get("science_measurement") is not False
            or summary.get("network_used") is not False
            or summary.get("cas_query_submitted") is not False
            or summary.get("new_data_downloaded") is not False
            or summary.get("row_plan_sha256") not in (None, EXPECTED_APPROVED_ROW_PLAN_SHA256)):
        raise ValueError("radial-support summary is incomplete or has an unexpected scope")
    qa_summary_path = Path(summary.get("provenance", {}).get("qa_summary_path", ""))
    if (qa_summary_path != QA_DIR / "summary.json"
            or sha256(qa_summary_path) != EXPECTED_SHA256["qa_summary.json"]):
        raise ValueError("radial-support audit is not tied to the pinned local QA summary")
    qa_summary = json.loads(qa_summary_path.read_text(encoding="utf-8"))
    if (qa_summary.get("row_plan_sha256") != EXPECTED_APPROVED_ROW_PLAN_SHA256
            or qa_summary.get("science_measurement") is not False
            or qa_summary.get("network_used") is not False
            or qa_summary.get("additional_data_downloaded") is not False):
        raise ValueError("local QA summary differs from the approved scope or row plan")
    if (summary.get("provenance", {}).get("qa_summary_sha256")
            != EXPECTED_SHA256["qa_summary.json"]):
        raise ValueError("radial-support summary QA hash differs from the pinned input")
    if summary.get("global_sample", {}).get("selected_lens_sets") != 74:
        raise ValueError("pinned radial audit does not describe the expected 74-lens support set")

    manifest_rows = read_tsv(paths["selected_lenses.tsv"])
    manifest = {str(row["lens_id"]): row for row in manifest_rows}
    if len(manifest) != 80 or len(manifest_rows) != 80:
        raise ValueError("approved manifest must contain exactly 80 unique lenses")
    for row in manifest.values():
        for field in ATTRIBUTES:
            value = float(row[field])
            if not math.isfinite(value):
                raise ValueError(f"non-finite cluster property {field}")

    radial_rows = read_tsv(paths["per_lens_radial_support.tsv"])
    by_lens: dict[str, list[dict]] = {}
    for row in radial_rows:
        lens_id = str(row["matched_lens_id"])
        if lens_id not in manifest:
            raise ValueError(f"support table contains an unrequested lens ID: {lens_id}")
        by_lens.setdefault(lens_id, []).append(row)
    if len(by_lens) != 74 or len(radial_rows) != 74 * 7:
        raise ValueError("radial table must contain 74 globally selected lenses × 7 annuli")
    global_ids = set(by_lens)
    if id_sha256(global_ids) != EXPECTED_GLOBAL_LENS_IDS_SHA256:
        raise ValueError("global 74-lens membership differs from the pinned radial audit")
    if summary.get("global_sample", {}).get("lens_id_sha256") != EXPECTED_GLOBAL_LENS_IDS_SHA256:
        raise ValueError("global lens ID hash disagrees with the radial-support summary")

    fixed_ids = set()
    for lens_id, rows in by_lens.items():
        if len(rows) != 7 or {int(row["radial_bin"]) for row in rows} != set(range(7)):
            raise ValueError(f"lens {lens_id} does not have all seven annuli")
        if all(int(row["lens_pair_rows_z234"]) > 0
               and int(row["n_positive_controls_z234"]) >= 4 for row in rows):
            fixed_ids.add(lens_id)
    fixed_evidence = summary.get("fixed_all_radius_support", {}).get("4", {})
    if (len(fixed_ids) != 52
            or id_sha256(fixed_ids) != EXPECTED_FIXED52_LENS_IDS_SHA256
            or fixed_evidence.get("n_lens_sets_supported_in_all_seven_annuli") != 52
            or fixed_evidence.get("lens_ids_sha256") != EXPECTED_FIXED52_LENS_IDS_SHA256):
        raise ValueError("derived fixed 52-lens set differs from the pinned support result")
    return summary, manifest, radial_rows


def make_summary(plan_dir: Path, support_dir: Path, output: Path) -> dict:
    support, manifest, radial_rows = load_inputs(plan_dir, support_dir)
    global_ids = {str(row["matched_lens_id"]) for row in radial_rows}
    fixed_ids = set()
    by_lens: dict[str, list[dict]] = {}
    for row in radial_rows:
        by_lens.setdefault(str(row["matched_lens_id"]), []).append(row)
    for lens_id, rows in by_lens.items():
        if all(int(row["lens_pair_rows_z234"]) > 0
               and int(row["n_positive_controls_z234"]) >= 4 for row in rows):
            fixed_ids.add(lens_id)
    comparison_ids = global_ids - fixed_ids
    globally_excluded_ids = set(manifest) - global_ids
    groups = {
        "fixed_all_annuli_52": fixed_ids,
        "global_support_not_fixed_22": comparison_ids,
        "global_support_74": global_ids,
        "outside_global_support_6": globally_excluded_ids,
        "all_manifest_lenses_80": set(manifest),
    }
    attributes = {}
    for name, ids in groups.items():
        attributes[name] = {
            field: distribution([float(manifest[lens_id][field]) for lens_id in sorted(ids)])
            for field in ATTRIBUTES
        }
    smd = {}
    for field in ATTRIBUTES:
        fixed = [float(manifest[lens_id][field]) for lens_id in sorted(fixed_ids)]
        comparison = [float(manifest[lens_id][field]) for lens_id in sorted(comparison_ids)]
        smd[field] = standardized_mean_difference(fixed, comparison)

    script_path = Path(__file__).resolve()
    return {
        "schema": "astrafek-s23b-wl-fixed-support-selection-audit-v1",
        "status": "complete_local_fixed_support_sample_selection_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "source_csvs_read": False,
        "scope": "descriptive comparison of cluster properties for the all-annulus 4/5-control sample against the other globally supported lenses",
        "fixed_sample_rule": "positive lens z-bin 2-4 source-pair support and at least four positive matched controls in every one of seven 0.1-3 proper-Mpc annuli",
        "group_counts": {name: len(ids) for name, ids in groups.items()},
        "cluster_property_distributions": attributes,
        "fixed52_minus_other22_standardized_mean_difference": smd,
        "interpretation": [
            "This is a descriptive selection audit, not a hypothesis test or a shear measurement.",
            "n_raw_shape_sources is a pre-cut catalog-availability proxy, not a calibrated depth map.",
            "The fixed sample is not automatically representative of all S23B clusters; interpreting an ICL-WL relation beyond these 52 lenses requires selection-aware analysis.",
        ],
        "provenance": {
            "runner_path": str(script_path),
            "runner_sha256": sha256(script_path),
            "approved_row_plan_sha256": EXPECTED_APPROVED_ROW_PLAN_SHA256,
            "selected_lenses_path": str(plan_dir / "selected_lenses.tsv"),
            "selected_lenses_sha256": EXPECTED_SHA256["selected_lenses.tsv"],
            "radial_summary_path": str(support_dir / "summary.json"),
            "radial_summary_sha256": EXPECTED_SHA256["support_summary.json"],
            "per_lens_radial_support_path": str(support_dir / "per_lens_radial_support.tsv"),
            "per_lens_radial_support_sha256": EXPECTED_SHA256["per_lens_radial_support.tsv"],
            "global_74_lens_ids_sha256": EXPECTED_GLOBAL_LENS_IDS_SHA256,
            "fixed_52_lens_ids_sha256": EXPECTED_FIXED52_LENS_IDS_SHA256,
        },
        "output_directory": str(output),
    }


def write_assessment(path: Path, summary: dict) -> None:
    groups = summary["cluster_property_distributions"]
    smd = summary["fixed52_minus_other22_standardized_mean_difference"]
    lines = [
        "# S23B fixed-support weak-lensing sample selection audit (2026-11-11)",
        "",
        "This offline audit compares the 52-lens all-annulus 4/5-control sample with the other 22 lenses in the globally supported 74-lens set. It reads only the approved 80-lens manifest and the hashed radial-support summary table; it reads no source CSVs or images and makes no network/CAS request.",
        "",
        "The fixed sample requires positive lens source support and at least four positive matched controls in each of seven 0.1–3 proper-Mpc annuli, using HSC Y3 source bins 2–4. The comparison group is the remaining 22 lenses in the parent 74-lens supported sample; six manifest lenses are outside that parent sample.",
        "",
        "| Cluster property | Fixed 52 mean (median) | Other 22 mean (median) | Fixed-minus-other SMD |",
        "|---|---:|---:|---:|",
    ]
    labels = {
        "z_lens": "Lens redshift",
        "richness": "CAMIRA richness",
        "n_raw_shape_sources": "Raw S19A shape-source count",
        "eligible_control_candidates_6mpc": "Eligible control candidates within 6 Mpc",
        "eligible_control_candidates_8mpc": "Eligible control candidates within 8 Mpc",
        "old_union_sources": "Previous union-source count",
        "old_union_area_fraction_4096": "Previous union area fraction",
    }
    for field in ATTRIBUTES:
        left = groups["fixed_all_annuli_52"][field]
        right = groups["global_support_not_fixed_22"][field]
        lines.append(
            f"| {labels[field]} | {left['mean']:.3f} ({left['median']:.3f}) | "
            f"{right['mean']:.3f} ({right['median']:.3f}) | {smd[field]:.3f} |"
        )
    lines.extend([
        "",
        "SMD is the fixed-minus-comparison mean difference divided by the pooled sample standard deviation. It is descriptive; no p-values or significance claims are made.",
        "",
        f"The fixed sample has a median raw shape-source count of {groups['fixed_all_annuli_52']['n_raw_shape_sources']['median']:.0f}, versus {groups['global_support_not_fixed_22']['n_raw_shape_sources']['median']:.0f} in the other 22 lenses. Its median redshift is {groups['fixed_all_annuli_52']['z_lens']['median']:.3f}, versus {groups['global_support_not_fixed_22']['z_lens']['median']:.3f}. The richness distribution is also reported because the fixed sample may not preserve the parent cluster population.",
        "",
        "The raw shape-source count measures pre-cut catalog availability, not calibrated depth. These comparisons do not correct selection bias. A lensing or ICL–lensing interpretation beyond the selected 52 needs a selection-aware model or an explicit restriction to this supported subset.",
        "",
        "`science_measurement=false`; no shear profile or halo mass is computed here.",
        "",
        "## Provenance",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Approved row plan SHA-256: `{summary['provenance']['approved_row_plan_sha256']}`",
        f"- Selected lens manifest SHA-256: `{summary['provenance']['selected_lenses_sha256']}`",
        f"- Radial support summary SHA-256: `{summary['provenance']['radial_summary_sha256']}`",
        f"- Per-lens radial support SHA-256: `{summary['provenance']['per_lens_radial_support_sha256']}`",
        f"- Output directory: `{summary['output_directory']}`",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def run(plan_dir: Path, support_dir: Path, output: Path) -> dict:
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    summary = make_summary(Path(plan_dir), Path(support_dir), output)
    output.mkdir(parents=True, exist_ok=False)
    write_assessment(output / "assessment.md", summary)
    summary["output_sha256"] = {
        "assessment.md": sha256(output / "assessment.md"),
    }
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "group_counts": summary["group_counts"],
        "fixed52_vs_other22": {
            field: {
                "fixed_median": summary["cluster_property_distributions"]["fixed_all_annuli_52"][field]["median"],
                "other_median": summary["cluster_property_distributions"]["global_support_not_fixed_22"][field]["median"],
                "smd": summary["fixed52_minus_other22_standardized_mean_difference"][field],
            }
            for field in ("z_lens", "richness", "n_raw_shape_sources")
        },
        "output_directory": str(output),
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=PLAN_DIR)
    parser.add_argument("--support-dir", type=Path, default=SUPPORT_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan_dir, args.support_dir, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
