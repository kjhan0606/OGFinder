#!/usr/bin/env python3
"""Compare ICL--satellite block-power needs across effect and test-family sizes.

This is a design diagnostic only. It uses the existing exploratory effect size,
observed connected-block count, and audited CAMIRA target-cell upper bound; it
does not select clusters, query HSC, or make a science measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from plan_icl_spatial_block_power import required_units, two_sided_power


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_POWER = DEFAULT_ROOT / "icl_association_power_plan_2026-10-19_01/summary.json"
DEFAULT_FEASIBILITY = DEFAULT_ROOT / "icl_catalog_expansion_feasibility_2026-10-20_04/feasibility.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inputs(power_path: Path, feasibility_path: Path) -> tuple[dict, dict]:
    power = json.loads(power_path.read_text(encoding="utf-8"))
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (power.get("schema") != "astrafek-icl-spatial-block-power-plan-v1"
            or power.get("status") != "planning_approximation"
            or power.get("science_measurement") is not False):
        raise ValueError("power-plan input is not the expected exploratory approximation")
    if (feasibility.get("schema") != "astrafek-icl-camira-expansion-feasibility-v1"
            or feasibility.get("status") != "complete_design_diagnostic"
            or feasibility.get("science_measurement") is not False
            or feasibility.get("cas_query_submitted") is not False
            or feasibility.get("coordinates_exported") is not False
            or feasibility.get("images_downloaded") is not False):
        raise ValueError("CAMIRA feasibility input is not the expected offline diagnostic")
    if power.get("observed_blocks") != feasibility.get("observed_independent_blocks"):
        raise ValueError("power and CAMIRA diagnostics disagree on observed block count")
    return power, feasibility


def assess(power: dict, feasibility: dict, *, effect_sizes=(0.2, 0.293, 0.4),
           family_sizes=(1, 6, 30), target_power=0.8, n_controls=2) -> dict:
    if (feasibility.get("cas_query_submitted") is not False
            or feasibility.get("coordinates_exported") is not False
            or feasibility.get("images_downloaded") is not False):
        raise ValueError("CAMIRA feasibility input must be an offline diagnostic")
    observed = int(power["observed_blocks"])
    upper_bound = int(feasibility["broadest_pool_target_cell_upper_bound_on_connected_blocks"])
    if observed <= n_controls + 3 or upper_bound < observed:
        raise ValueError("invalid observed block count or expansion upper bound")
    if not 0 < float(target_power) < 1:
        raise ValueError("target_power must be between zero and one")
    if not effect_sizes or any(not 0 < float(r) < 1 for r in effect_sizes):
        raise ValueError("effect sizes must be finite correlations in (0, 1)")
    if not family_sizes or any(int(size) < 1 for size in family_sizes):
        raise ValueError("test-family sizes must be positive integers")

    current_r = float(power["observed_effect_r_scenario"])
    rows = []
    for effect in map(float, effect_sizes):
        for family_size in map(int, family_sizes):
            alpha = 0.05 / family_size
            required = required_units(effect, alpha, target_power, n_controls)
            rows.append({
                "effect_r_scenario": effect,
                "test_family_size_scenario": family_size,
                "two_sided_alpha_after_bonferroni": alpha,
                "target_power": float(target_power),
                "required_independent_blocks": required,
                "observed_blocks": observed,
                "observed_power": two_sided_power(observed, effect, alpha, n_controls),
                "candidate_pool_target_cell_upper_bound": upper_bound,
                "upper_bound_power": two_sided_power(upper_bound, effect, alpha, n_controls),
                "shortfall_even_at_upper_bound": max(0, required - upper_bound),
                "exploratory_effect_size_matches_observed": effect == current_r,
            })
    return {
        "schema": "astrafek-icl-power-family-sensitivity-v1",
        "status": "planning_approximation",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "coordinates_exported": False,
        "images_downloaded": False,
        "method": "two-sided Fisher-z approximation for partial correlation; Bonferroni alpha within each hypothetical test family",
        "observed_connected_blocks": observed,
        "broadest_existing_pool_connected_block_upper_bound": upper_bound,
        "observed_exploratory_effect_r": current_r,
        "n_adjustment_variables": int(power["n_controls"]),
        "target_power": float(target_power),
        "scenarios": rows,
        "limitations": [
            "The observed r=0.293 was selected after exploratory analyses and may overestimate the underlying association.",
            "A one-test family is only a prospective design scenario; it is not permission to choose an endpoint after seeing the data.",
            "The 84 target-cell count is an upper bound on connected blocks, not a count of usable independent blocks.",
            "Fisher-z power does not replace spatial-block maxT inference or account for selection, sky, PSF, mask, or measurement-systematic uncertainty.",
            "Every output remains science_measurement=false and does not establish an ICL--satellite detection.",
        ],
    }


def render(result: dict, power_path: Path, feasibility_path: Path,
           hashes: dict[str, str]) -> str:
    lines = [
        "# ICL–위성 연관성 공간블록 검정력 민감도",
        "",
        "이 설계 진단은 현재 33개 연결블록과 기존 CAMIRA 후보 풀이 제공할 수 있는 target-cell 상한 84개를 사용한다. 관측된 `r=0.293`은 탐색 결과이므로 실제 효과로 간주하지 않는다. 단일 검정은 향후 분석 전에 하나의 주평가변수를 등록하는 가상 시나리오이며, 기존 30개 비교 결과를 사후에 축소하는 근거가 아니다.",
        "",
        "| partial r 시나리오 | 가설적 검정 수 | 관측 33블록 power | 필요 독립블록 (80%) | 84블록 상한에서 power | 상한 부족분 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["scenarios"]:
        lines.append(
            f"| {row['effect_r_scenario']:.3f} | {row['test_family_size_scenario']} "
            f"| {row['observed_power']:.3f} | {row['required_independent_blocks']} "
            f"| {row['upper_bound_power']:.3f} | {row['shortfall_even_at_upper_bound']} |"
        )
    lines.extend([
        "",
        "`r=0.293`에서 Bonferroni 단일 주검정은 약 92블록, 30개 비교는 180블록이 필요하다. 그러나 단일검정의 92블록도 84블록 상한을 넘고, 더 보수적인 `r=0.2` 단일검정은 196블록을 요구한다. 따라서 같은 CAMIRA 발자국의 후보 추가만으로 충분한 검정력을 확보한다고 볼 근거가 없다. 다중비교 수를 줄이려면 새 독립 자료를 보기 전에 과학적으로 정한 한 가지 주평가변수와 분석 규약을 사전등록해야 한다.",
        "",
        "모든 수치는 Fisher-z 정규근사 설계값이며, 공간블록 독립성이나 과학적 검출을 증명하지 않는다. 모든 결과는 `science_measurement=false`다.",
        "",
        f"Power 입력: `{power_path}` (SHA-256 `{hashes[str(power_path)]}`)",
        f"후보풀 입력: `{feasibility_path}` (SHA-256 `{hashes[str(feasibility_path)]}`)",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power", type=Path, default=DEFAULT_POWER)
    parser.add_argument("--feasibility", type=Path, default=DEFAULT_FEASIBILITY)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    power_path, feasibility_path = args.power.resolve(), args.feasibility.resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"refusing to overwrite existing output: {work}")
    power, feasibility = load_inputs(power_path, feasibility_path)
    result = assess(power, feasibility)
    hashes = {str(path): sha256(path) for path in (power_path, feasibility_path)}
    result["input_hashes"] = hashes
    result["script_sha256"] = sha256(Path(__file__).resolve())
    work.mkdir(parents=True, exist_ok=False)
    (work / "sensitivity.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (work / "assessment.md").write_text(
        render(result, power_path, feasibility_path, hashes), encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "observed_blocks": result["observed_connected_blocks"],
        "candidate_upper_bound": result["broadest_existing_pool_connected_block_upper_bound"],
        "scenarios": [{
            "r": row["effect_r_scenario"], "tests": row["test_family_size_scenario"],
            "required_blocks": row["required_independent_blocks"],
            "power_at_upper_bound": round(row["upper_bound_power"], 4),
        } for row in result["scenarios"]],
        "science_measurement": result["science_measurement"],
        "output": str(work),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
