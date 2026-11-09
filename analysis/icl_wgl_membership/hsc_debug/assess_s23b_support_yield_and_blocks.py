#!/usr/bin/env python3
"""Bound S23B expansion value using observed ICL/satellite support and sky cells.

This offline design analysis combines the existing S23B candidate ledger,
catalog-only expansion ceiling, observed 60% and 40% satellite support yields,
and Fisher-z spatial-block power approximation. It requests no data, selects
no clusters, and makes no science measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from plan_icl_spatial_block_power import required_units


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUTS = {
    "power": DATA / "icl_association_power_plan_2026-10-19_01/summary.json",
    "s23b_expansion": DATA / "icl_s23b_expansion_feasibility_2026-10-21_02/summary.json",
    "s23b_ledger": DATA / "icl_s23b_spatial_candidate_ledger_2026-10-22_02/summary.json",
    "support_60": DATA / "icl_satellite_support_yield_2026-11-05_01/summary.json",
    "support_40": DATA / "icl_satellite_coverage_threshold_2026-11-05_02/summary.json",
}
DEFAULT_OUTPUT = DATA / "icl_s23b_support_yield_and_blocks_2026-11-05_01"
OUTER_BIN = "50-75"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_inputs(paths: dict[str, Path]) -> dict[str, dict]:
    result = {name: json.loads(path.read_text(encoding="utf-8"))
              for name, path in paths.items()}
    expected = {
        "power": ("astrafek-icl-spatial-block-power-plan-v1", "planning_approximation"),
        "s23b_expansion": ("astrafek-icl-s23b-expansion-feasibility-v1", "offline_catalog_projection_only"),
        "s23b_ledger": ("astrafek-s23b-spatial-candidate-ledger-v1", "offline_coordinate_and_footprint_ledger_complete"),
        "support_60": ("astrafek-icl-satellite-support-yield-v1", "complete_offline_diagnostic"),
        "support_40": ("astrafek-icl-satellite-coverage-threshold-sensitivity-v1", "complete_exploratory_sensitivity"),
    }
    for name, doc in result.items():
        schema, status = expected[name]
        if doc.get("schema") != schema or doc.get("status") != status:
            raise ValueError(f"unexpected {name} input schema/status")
        if doc.get("science_measurement") is not False:
            raise ValueError(f"{name} input no longer has science_measurement=false")
    if (result["s23b_expansion"].get("hsc_query_submitted") is not False
            or result["s23b_expansion"].get("source_shape_photoz_or_image_data_received") is not False
            or result["s23b_ledger"].get("hsc_query_submitted") is not False
            or result["s23b_ledger"].get("data_downloaded") is not False
            or result["support_40"].get("data_requested") is not False
            or result["support_40"].get("hsc_network_used") is not False):
        raise ValueError("one or more inputs include an unexpected HSC request or transfer")
    return result


def assess(docs: dict[str, dict]) -> dict:
    expansion = docs["s23b_expansion"]
    ledger = docs["s23b_ledger"]
    support60 = docs["support_60"]
    support40 = docs["support_40"]
    power = docs["power"]

    n_parent = int(expansion["current_parent_lenses"])
    n_candidates = int(expansion["offline_new_candidate_upper_bound"])
    n_projected_parent = int(expansion["catalog_only_parent_upper_bound"])
    stats = ledger["statistics"]
    n_new_cells = int(stats["candidate_cells_new_vs_parent"])
    n_current_cells = int(stats["parent_unique_5deg_cells"])
    n_total_cells = int(stats["candidate_unique_5deg_cells"])
    if n_parent != 800 or n_projected_parent != n_parent + n_candidates:
        raise ValueError("S23B catalog-only parent ceiling is internally inconsistent")
    if n_total_cells != n_current_cells + n_new_cells:
        raise ValueError("S23B spatial-cell counts do not reconcile")

    p60 = support60["support_intersection_by_radial_bin_filter"]["i2"][OUTER_BIN]
    p40 = support40["results_by_filter"]["i2"]
    rows = {
        "coverage_60pct": {
            "observed_i2_pairs": int(p60["common_cluster_filter_pairs"]),
            "observed_center_eligible_pairs": int(p60["common_center_eligible_pairs"]),
            "partial_r": float(
                support60["satellite_candidate_diagnostics"].get("i2_partial_r", 0.0)
                if isinstance(support60.get("satellite_candidate_diagnostics"), dict)
                else p60.get("partial_r", 0.0)
            ),
        },
        "coverage_40pct": {
            "observed_i2_pairs": int(p40["n_40pct_support"]),
            "observed_center_eligible_pairs": int(p40["center_eligible_pairs_40pct"]),
            "partial_r": float(p40["partial_r_40pct"]),
            "partial_r_ci95": p40["partial_r_40pct_block_bootstrap_ci95"],
        },
    }
    # The 60% association estimate is stored in the paired 40/60 sensitivity;
    # require it to agree with the earlier 60% support assessment.
    rows["coverage_60pct"]["partial_r"] = float(p40["partial_r_60pct"])
    rows["coverage_60pct"]["partial_r_ci95"] = p40["partial_r_60pct_block_bootstrap_ci95"]

    for scenario in rows.values():
        scenario["observed_support_fraction"] = scenario["observed_i2_pairs"] / n_parent
        scenario["observed_center_eligible_fraction"] = scenario["observed_center_eligible_pairs"] / n_parent
        scenario["projected_i2_pairs_after_catalog_pool"] = (
            scenario["observed_support_fraction"] * n_projected_parent)
        scenario["projected_center_eligible_pairs_after_catalog_pool"] = (
            scenario["observed_center_eligible_fraction"] * n_projected_parent)
        scenario["projected_added_i2_pairs_from_candidates"] = (
            scenario["observed_support_fraction"] * n_candidates)

    comparisons = int(power["n_comparisons"])
    target_power = 0.8
    for name, scenario in rows.items():
        r = scenario["partial_r"]
        scenario["required_blocks_single_endpoint_alpha_0p05"] = required_units(
            r, 0.05, target_power, int(power["n_controls"]))
        scenario["required_blocks_bonferroni_30_tests"] = required_units(
            r, 0.05 / comparisons, target_power, int(power["n_controls"]))

    single_endpoint_r = float(power["observed_effect_r_scenario"])
    single_endpoint_req = required_units(single_endpoint_r, 0.05, target_power,
                                         int(power["n_controls"]))
    familywise_req = required_units(single_endpoint_r, 0.05 / comparisons,
                                    target_power, int(power["n_controls"]))
    return {
        "schema": "astrafek-s23b-icl-support-yield-and-block-cap-v1",
        "status": "complete_offline_design_diagnostic",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "image_pixels_read": False,
        "n_parent_lenses": n_parent,
        "n_catalog_only_candidate_upper_bound": n_candidates,
        "n_parent_lenses_after_catalog_upper_bound": n_projected_parent,
        "i2_outer_annulus_kpc": [50, 75],
        "support_projection_scenarios": rows,
        "spatial_power": {
            "observed_block_count_in_discovery": int(power["observed_blocks"]),
            "power_target": target_power,
            "familywise_test_count": comparisons,
            "exploratory_effect_r_scenario": single_endpoint_r,
            "required_blocks_single_predeclared_endpoint_alpha_0p05": single_endpoint_req,
            "required_blocks_bonferroni_30_tests": familywise_req,
            "s23b_candidate_cells_new_vs_discovery": n_new_cells,
            "s23b_candidate_unique_5deg_cells": n_total_cells,
            "discovery_parent_unique_5deg_cells": n_current_cells,
            "combined_target_cell_upper_bound": n_total_cells,
            "independent_s23b_cell_upper_bound_relative_to_discovery": n_new_cells,
            "single_endpoint_shortfall_from_s23b_only_new_cells": max(0, single_endpoint_req - n_new_cells),
            "single_endpoint_shortfall_from_combined_cell_upper_bound": max(0, single_endpoint_req - n_total_cells),
            "familywise_shortfall_from_combined_cell_upper_bound": max(0, familywise_req - n_total_cells),
        },
        "inputs": {name: {"path": str(INPUTS[name]), "sha256": sha256(INPUTS[name])}
                   for name in INPUTS},
        "interpretation": (
            "The 40% support transfer projects more ICL/satellite rows, but its "
            "observed i2 interval includes zero. Even the full S23B candidate "
            "pool contributes only 22 target cells outside the discovery sample "
            "and 100 cells combined. These are catalog and spatial-cell upper "
            "bounds, not available images, independent blocks, or a measurement."
        ),
        "limitations": [
            "Candidate coadd, mask, depth, PSF, member completeness, shape support, and matched controls are unverified.",
            "Support rates are transferred from the selected 800-lens sample and may not apply to S23B candidates.",
            "The 40% support floor changes the satellite estimator and its observed i2 partial-r point interval reaches zero.",
            "Counts of candidate cluster/filter pairs are not counts of independent spatial blocks.",
            "Five-degree target cells only upper-bound connected blocks under the current block-label rule and do not establish statistical independence.",
            "Fisher-z power is a planning approximation using a post-selection effect-size scenario; it is not the block bootstrap or maxT analysis.",
        ],
    }


def render(result: dict) -> str:
    rows = result["support_projection_scenarios"]
    power = result["spatial_power"]
    lines = [
        "# S23B 후보 확장: 지원률과 독립 공간표본 상한 — 2026-11-05",
        "",
        "## 요약",
        "",
        f"현재 {result['n_parent_lenses']}개 표본에 catalog-only S23B 후보 상한 "
        f"{result['n_catalog_only_candidate_upper_bound']}개를 더하면 총 "
        f"{result['n_parent_lenses_after_catalog_upper_bound']}개다. 지원률이 그대로 유지된다고 가정할 때 50–75 kpc `i2` "
        f"pair는 60% 기준 약 {rows['coverage_60pct']['projected_i2_pairs_after_catalog_pool']:.0f}개, "
        f"40% 기준 약 {rows['coverage_40pct']['projected_i2_pairs_after_catalog_pool']:.0f}개로 예상된다. "
        "이는 영상 가용성과 매칭 대조장을 모르는 낙관적 계획값이다.",
        "",
        "| Satellite coverage floor | Current i2 pairs | Current center-eligible | Observed partial r (95% block interval) | Projected pairs after all candidates | Projected center-eligible |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for name, label in (("coverage_60pct", "60%"), ("coverage_40pct", "40%")):
        row = rows[name]
        ci = row["partial_r_ci95"]
        lines.append(
            f"| {label} | {row['observed_i2_pairs']} | {row['observed_center_eligible_pairs']} "
            f"| {row['partial_r']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] "
            f"| {row['projected_i2_pairs_after_catalog_pool']:.1f} "
            f"| {row['projected_center_eligible_pairs_after_catalog_pool']:.1f} |"
        )
    lines.extend([
        "",
        "40% 기준은 현재 표본의 pair 수를 늘리지만 i2 부분상관 구간이 0을 포함하고 control 선택도 달라진다. Baseline으로 채택하거나 60% 결과와 합치지 않는다.",
        "",
        "## 공간표본 상한",
        "",
        f"30개 비교를 보정하는 Fisher-z 계획에서 현재 40% 점추정 r={rows['coverage_40pct']['partial_r']:.3f}를 가정해도 80% 검정력에 필요한 독립 블록은 {rows['coverage_40pct']['required_blocks_bonferroni_30_tests']}개다. 사후 탐색값 r={power['exploratory_effect_r_scenario']:.3f}에 대한 단일 사전 지정 endpoint 요구량은 {power['required_blocks_single_predeclared_endpoint_alpha_0p05']}개, 30개 비교 보정 요구량은 {power['required_blocks_bonferroni_30_tests']}개다.",
        "",
        f"S23B candidate ledger에는 전체 100개 5° 셀이 있지만 기존 discovery와 겹치지 않는 셀은 {power['s23b_candidate_cells_new_vs_discovery']}개뿐이다. 따라서 S23B만으로 독립 재현에 추가할 수 있는 target 셀 상한은 {power['independent_s23b_cell_upper_bound_relative_to_discovery']}개이며, 현재 800개 표본과 합친 셀 상한 {power['combined_target_cell_upper_bound']}개도 가족보정 계획의 {power['required_blocks_bonferroni_30_tests']}개보다 작다. 셀은 실제 독립 공분산 블록이 아니며, 대조장 연결은 상한을 더 줄일 수 있다.",
        "",
        "## 결정",
        "",
        "S23B 이미지를 모두 받아도 이 catalog pool은 confirmatory한 공간 독립성 목표를 충족할 수 없다. 다만 탐색적 stack 확장에는 도움이 될 수 있다. 우선 지원 표본의 PSF/selection 검증을 마치고, 독립 블록이 더 큰 HSC 영역 또는 다른 적격 공개 footprint를 찾은 뒤 새 영상 수신을 계획한다. 본 분석은 카탈로그·지원률 계획만 읽었으며 CAS 조회와 데이터 수신은 없고 `science_measurement=false`다.",
        "",
        "## 재현성",
        "",
        "- Astrafek script: `assess_s23b_support_yield_and_blocks.py`",
        f"- Output: `{DEFAULT_OUTPUT}`",
        "- Machine summary stores every input SHA-256.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    docs = read_inputs(INPUTS)
    result = assess(docs)
    result["script_sha256"] = sha256(Path(__file__).resolve())
    output.mkdir(parents=True, exist_ok=False)
    summary = output / "summary.json"
    summary.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = output / "assessment.md"
    report.write_text(render(result), encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "projected_i2_pairs_60pct": round(result["support_projection_scenarios"]["coverage_60pct"]["projected_i2_pairs_after_catalog_pool"], 1),
        "projected_i2_pairs_40pct": round(result["support_projection_scenarios"]["coverage_40pct"]["projected_i2_pairs_after_catalog_pool"], 1),
        "new_independent_target_cells_upper_bound": result["spatial_power"]["independent_s23b_cell_upper_bound_relative_to_discovery"],
        "combined_target_cell_upper_bound": result["spatial_power"]["combined_target_cell_upper_bound"],
        "required_blocks_single_endpoint": result["spatial_power"]["required_blocks_single_predeclared_endpoint_alpha_0p05"],
        "required_blocks_familywise": result["spatial_power"]["required_blocks_bonferroni_30_tests"],
        "output": str(output),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
