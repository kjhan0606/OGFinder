#!/usr/bin/env python3
"""Audit same-cluster overlap of S23B ICL/satellite and weak-lensing results.

The script joins already saved, pinned pair tables for the fixed 52-lens
weak-lensing sample. It reads no raw source CSVs or science-image pixels and
makes no network/CAS request. It reports overlap and support only; it does
not fit a joint ICL--satellite--lensing relation.
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
SUPPORT_DIR = DATA / "wl_s23b_matched_radial_support_2026-11-11_02"
WL_DIR = DATA / "wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_01"
ICL_DIR = DATA / "s23b_icl_satellite_relation_2026-11-07_03"
DEFAULT_OUTPUT = DATA / "s23b_icl_wl_fixed_support_overlap_2026-11-11_01"

EXPECTED_SHA256 = {
    "support_summary.json": "c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234",
    "per_lens_radial_support.tsv": "6744ece545781faffebef65fc7a0143284da3a234425b6bed540bd1e9d19625c",
    "wl_summary.json": "2391fede047fc2de25b34082b0c0d6635e88910d5701a6f9d16b5498ec66602f",
    "wl_per_lens_set_profile.tsv": "5bc599d97123264058011313b0ad35c3def49b08bba98a016cc5edf6ff03906c",
    "icl_summary.json": "fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39",
    "icl_pair_table.tsv": "98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994",
}
EXPECTED_ROW_PLAN_SHA256 = "7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc"
EXPECTED_FIXED52_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def id_sha256(values: set[str]) -> str:
    return hashlib.sha256("\n".join(sorted(values)).encode()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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
    return {
        "n": len(values),
        "q16": quantile(values, 0.16),
        "median": statistics.median(values) if values else None,
        "q84": quantile(values, 0.84),
        "min": min(values) if values else None,
        "max": max(values) if values else None,
    }


def load_inputs() -> tuple[dict, dict, list[dict], dict, list[dict], list[dict]]:
    paths = {
        "support_summary.json": SUPPORT_DIR / "summary.json",
        "per_lens_radial_support.tsv": SUPPORT_DIR / "per_lens_radial_support.tsv",
        "wl_summary.json": WL_DIR / "summary.json",
        "wl_per_lens_set_profile.tsv": WL_DIR / "per_lens_set_profile.tsv",
        "icl_summary.json": ICL_DIR / "summary.json",
        "icl_pair_table.tsv": ICL_DIR / "paired_icl_satellite_profiles.tsv",
    }
    for name, path in paths.items():
        actual = sha256(path)
        if actual != EXPECTED_SHA256[name]:
            raise ValueError(f"pinned input hash changed for {name}: {actual}")

    support_summary = read_json(paths["support_summary.json"])
    wl_summary = read_json(paths["wl_summary.json"])
    icl_summary = read_json(paths["icl_summary.json"])
    radial_rows = read_tsv(paths["per_lens_radial_support.tsv"])
    wl_rows = read_tsv(paths["wl_per_lens_set_profile.tsv"])
    icl_rows = read_tsv(paths["icl_pair_table.tsv"])

    if (support_summary.get("science_measurement") is not False
            or support_summary.get("network_used") is not False
            or support_summary.get("new_data_downloaded") is not False
            or support_summary.get("global_sample", {}).get("selected_lens_sets") != 74):
        raise ValueError("radial-support input has unexpected scope or sample")
    qa_path = Path(support_summary["provenance"]["qa_summary_path"])
    qa_summary = read_json(qa_path)
    if (sha256(qa_path) != support_summary["provenance"]["qa_summary_sha256"]
            or qa_summary.get("row_plan_sha256") != EXPECTED_ROW_PLAN_SHA256
            or qa_summary.get("science_measurement") is not False):
        raise ValueError("support table no longer points to the approved local shape-row audit")
    if (wl_summary.get("science_measurement") is not False
            or wl_summary.get("network_used") is not False
            or wl_summary.get("new_data_downloaded") is not False
            or wl_summary.get("sample", {}).get("n_fixed_lens_sets") != 52
            or wl_summary.get("sample", {}).get("fixed_lens_ids_sha256") != EXPECTED_FIXED52_SHA256
            or wl_summary.get("provenance", {}).get("approved_row_plan_sha256") != EXPECTED_ROW_PLAN_SHA256):
        raise ValueError("weak-lensing profile is not the pinned exploratory 52-lens sample")
    if (icl_summary.get("science_measurement") is not False
            or icl_summary.get("sample", {}).get("science_image_pixels_read") is not False
            or icl_summary.get("sample", {}).get("network_used") is not False
            or icl_summary.get("sample", {}).get("supported_pairs_by_filter") != {"i": 232, "i2": 756}):
        raise ValueError("ICL/satellite profile is incomplete or has unexpected scope")
    return support_summary, wl_summary, radial_rows, icl_summary, icl_rows, wl_rows


def run(output: Path) -> dict:
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    support_summary, wl_summary, radial_rows, icl_summary, icl_rows, wl_rows = load_inputs()

    by_lens: dict[str, list[dict]] = {}
    for row in radial_rows:
        by_lens.setdefault(str(row["matched_lens_id"]), []).append(row)
    fixed_ids = {
        lens_id for lens_id, rows in by_lens.items()
        if len(rows) == 7 and all(
            int(row["lens_pair_rows_z234"]) > 0
            and int(row["n_positive_controls_z234"]) >= 4
            for row in rows
        )
    }
    if len(fixed_ids) != 52 or id_sha256(fixed_ids) != EXPECTED_FIXED52_SHA256:
        raise ValueError("derived fixed 52-lens set differs from the pinned sample")

    wl_by_lens: dict[str, list[dict]] = {}
    for row in wl_rows:
        lens_id = str(row["matched_lens_id"])
        if lens_id not in fixed_ids:
            raise ValueError(f"weak-lensing profile contains a non-fixed lens: {lens_id}")
        wl_by_lens.setdefault(lens_id, []).append(row)
        for field in ("lens_delta_sigma_msun_mpc2", "random_delta_sigma_msun_mpc2",
                      "delta_sigma_random_subtracted_msun_mpc2", "cross_random_subtracted_msun_mpc2"):
            if not math.isfinite(float(row[field])):
                raise ValueError(f"non-finite weak-lensing value for {lens_id}, {field}")
    if (len(wl_by_lens) != 52 or len(wl_rows) != 52 * 7
            or any(len(rows) != 7 for rows in wl_by_lens.values())):
        raise ValueError("weak-lensing table does not cover all 52 lenses in seven annuli")

    icl_by_filter: dict[str, list[dict]] = {"i": [], "i2": []}
    seen: set[tuple[str, str]] = set()
    for row in icl_rows:
        target_id, band = str(row["target_id"]), str(row["filter"])
        if band not in icl_by_filter or target_id not in fixed_ids:
            continue
        key = (target_id, band)
        if key in seen:
            raise ValueError(f"duplicate ICL target/filter record: {key}")
        seen.add(key)
        if str(row["science_measurement"]).lower() != "false":
            raise ValueError("an ICL/satellite input row unexpectedly promotes a science measurement")
        icl_by_filter[band].append(row)

    by_filter = {}
    filter_lens_ids = {}
    for band, rows in icl_by_filter.items():
        ids = {str(row["target_id"]) for row in rows}
        filter_lens_ids[band] = ids
        support = [min(float(row["inner_common_support_fraction"]),
                       float(row["outer_common_support_fraction"])) for row in rows]
        if not all(math.isfinite(value) and 0 <= value <= 1 for value in support):
            raise ValueError(f"invalid ICL common support fraction in filter {band}")
        by_filter[band] = {
            "n_fixed_wl_lenses": len(fixed_ids),
            "n_icl_satellite_profiles": len(rows),
            "n_unique_lenses": len(ids),
            "n_lenses_with_min_both_annulus_support": {
                str(floor): sum(value >= floor for value in support)
                for floor in (0.05, 0.10, 0.20)
            },
            "minimum_inner_outer_common_support_fraction": distribution(support),
            "inner_common_support_fraction": distribution([
                float(row["inner_common_support_fraction"]) for row in rows]),
            "outer_common_support_fraction": distribution([
                float(row["outer_common_support_fraction"]) for row in rows]),
            "icl_relative_contrast_njy_arcsec2": distribution([
                float(row["icl_relative_contrast_njy_arcsec2"]) for row in rows]),
            "log10_inner_outer_candidate_density_ratio": distribution([
                float(row["log10_inner_outer_candidate_density_ratio"]) for row in rows]),
            "member_candidates_inner_median": statistics.median([
                int(row["camira_member_candidates_50_300_in_image"]) for row in rows]) if rows else None,
            "member_candidates_outer_median": statistics.median([
                int(row["camira_member_candidates_300_800_in_image"]) for row in rows]) if rows else None,
        }

    wl_ids = set(wl_by_lens)
    icl_ids = filter_lens_ids["i"] | filter_lens_ids["i2"]
    both_filters = filter_lens_ids["i"] & filter_lens_ids["i2"]
    if (len(icl_rows) != 988 or len(icl_by_filter["i"]) != 25
            or len(icl_by_filter["i2"]) != 27 or icl_ids != wl_ids or both_filters):
        raise ValueError("S23B ICL/WL overlap changed from the pinned 52-lens cross-match")

    script_path = Path(__file__).resolve()
    summary = {
        "schema": "astrafek-s23b-icl-wl-fixed-support-overlap-v1",
        "status": "complete_local_icl_satellite_weak_lensing_overlap_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "source_csvs_read": False,
        "joint_relation_fit_performed": False,
        "wl_annuli_mpc": [[0.1, 0.2], [0.2, 0.4], [0.4, 0.7], [0.7, 1.0], [1.0, 1.5], [1.5, 2.0], [2.0, 3.0]],
        "fixed_lens_count": len(fixed_ids),
        "fixed_lens_ids_sha256": id_sha256(fixed_ids),
        "overlap": {
            "n_fixed_lenses_with_wl_profile": len(wl_ids),
            "n_fixed_lenses_with_supported_icl_satellite_profile": len(icl_ids),
            "n_fixed_lenses_with_both_wl_and_icl_satellite_products": len(wl_ids & icl_ids),
            "n_lenses_with_both_i_and_i2_profiles": len(both_filters),
            "i_only_lenses": len(filter_lens_ids["i"] - filter_lens_ids["i2"]),
            "i2_only_lenses": len(filter_lens_ids["i2"] - filter_lens_ids["i"]),
        },
        "by_filter": by_filter,
        "interpretation": [
            "All 52 fixed weak-lensing lenses have a supported ICL/satellite profile in exactly one of i or i2; none has both filters in this overlap.",
            "The filter-specific joint samples contain 25 and 27 lenses. Requiring at least 10% common pixels in both ICL annuli reduces them to 18 and 22; a 20% floor leaves 5 and 7.",
            "This sample supports only low-dimensional exploratory joint analysis. Sparse common-pixel coverage, the filter split, and unresolved WL/ICL systematics prevent a confirmatory claim.",
        ],
        "provenance": {
            "runner_path": str(script_path),
            "runner_sha256": sha256(script_path),
            "approved_row_plan_sha256": EXPECTED_ROW_PLAN_SHA256,
            "fixed_lens_id_sha256": EXPECTED_FIXED52_SHA256,
            "support_summary_sha256": EXPECTED_SHA256["support_summary.json"],
            "per_lens_radial_support_sha256": EXPECTED_SHA256["per_lens_radial_support.tsv"],
            "wl_summary_sha256": EXPECTED_SHA256["wl_summary.json"],
            "wl_per_lens_profile_sha256": EXPECTED_SHA256["wl_per_lens_set_profile.tsv"],
            "icl_summary_sha256": EXPECTED_SHA256["icl_summary.json"],
            "icl_pair_table_sha256": EXPECTED_SHA256["icl_pair_table.tsv"],
        },
        "output_directory": str(output),
    }
    output.mkdir(parents=True, exist_ok=False)
    assessment = output / "assessment.md"
    lines = [
        "# S23B ICL–위성–weak-lensing 고정표본 교집합 점검 (2026-11-11)",
        "",
        "이미 저장된 S23B CAMIRA member-position/ICL profile과 고정 52-lens weak-lensing per-lens profile을 해시 검증해 결합 가능 표본을 셌다. 원시 source CSV와 영상 픽셀은 읽지 않았고, CAS/네트워크 요청 및 데이터 다운로드도 없었다. 공동 관계 적합은 하지 않았으며 `science_measurement=false`를 유지한다.",
        "",
        "모든 52개 weak-lensing 렌즈에 ICL–위성 결과가 있지만, 각 렌즈는 i 또는 i2 중 한 필터에만 나타난다. 두 필터 결과가 모두 있는 렌즈는 0개다.",
        "",
        "| 필터 | 고정 WL 표본 | ICL–위성 공통 지원 렌즈 | 안쪽/바깥쪽 각각 ≥10% 공통 픽셀 | 각각 ≥20% | 중앙 안쪽/바깥쪽 공통 픽셀 비율 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        row = by_filter[band]
        lines.append(
            f"| {band} | 52 | {row['n_unique_lenses']} | "
            f"{row['n_lenses_with_min_both_annulus_support']['0.1']} | "
            f"{row['n_lenses_with_min_both_annulus_support']['0.2']} | "
            f"{row['inner_common_support_fraction']['median']:.1%} / "
            f"{row['outer_common_support_fraction']['median']:.1%} |"
        )
    lines.extend([
        "",
        "The 10% and 20% thresholds require both the 50–300 kpc and 300–800 kpc ICL annuli to meet the floor. The full fixed WL sample has 52×7 per-lens radial estimates, but the ICL/filter-specific matched sample has only 25 i and 27 i2 lenses. The two band samples are disjoint.",
        "",
        "The overlap is sufficient to prepare a low-dimensional exploratory relation, not a stable multi-covariate detection. ICL common-pixel support is modest and drops sharply under a 20% floor; WL source selection also favors higher raw shape-catalog availability and somewhat lower CAMIRA richness than the other 22 globally supported lenses. Any joint result must remain scoped to this selected subset and account for sky covariance and image support.",
        "",
        "## Provenance",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Fixed 52-lens ID SHA-256: `{summary['fixed_lens_ids_sha256']}`",
        f"- ICL–satellite pair-table SHA-256: `{summary['provenance']['icl_pair_table_sha256']}`",
        f"- WL per-lens profile SHA-256: `{summary['provenance']['wl_per_lens_profile_sha256']}`",
        f"- Output directory: `{output}`",
        "",
    ])
    assessment.write_text("\n".join(lines), encoding="utf-8")
    summary["output_sha256"] = {"assessment.md": sha256(assessment)}
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "overlap": summary["overlap"],
        "by_filter": {
            band: {
                "n_lenses": stats["n_unique_lenses"],
                "n_ge10pct_both_annuli": stats["n_lenses_with_min_both_annulus_support"]["0.1"],
                "n_ge20pct_both_annuli": stats["n_lenses_with_min_both_annulus_support"]["0.2"],
            }
            for band, stats in by_filter.items()
        },
        "output_directory": str(output),
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
