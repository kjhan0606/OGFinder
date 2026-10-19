#!/usr/bin/env python3
"""Compare ICL--satellite block-power needs with the local CAMIRA pool.

This is an offline design diagnostic. It reads the existing exploratory power
summary and the audited CAMIRA expansion-pool summary; it never queries HSC,
selects new clusters, downloads images, or makes a science measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_POWER = DEFAULT_ROOT / "icl_association_power_plan_2026-10-19_01/summary.json"
DEFAULT_EXPANSION = DEFAULT_ROOT / "camira_expansion_pool_2026-10-19_02/summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_summary(path: Path, expected_schema: str) -> dict:
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("schema") != expected_schema:
        raise ValueError(f"unexpected summary schema in {path}")
    if result.get("science_measurement") is not False:
        raise ValueError(f"input must remain explicitly exploratory: {path}")
    return result


def assess(power: dict, expansion: dict) -> dict:
    if power.get("status") != "planning_approximation":
        raise ValueError("power input is not the expected planning approximation")
    if expansion.get("status") != "complete_catalog_pool_diagnostic":
        raise ValueError("CAMIRA expansion input is not a completed pool diagnostic")
    if expansion.get("cas_query_submitted") is not False \
            or expansion.get("coordinates_exported") is not False \
            or expansion.get("images_downloaded") is not False:
        raise ValueError("CAMIRA expansion input must be offline and unselected")

    scenarios = expansion.get("selection_scenarios", {})
    broadest = scenarios.get("wide_z_cluster_redshift")
    current = scenarios.get("current_selection")
    if not broadest or not current:
        raise ValueError("required current and broadest CAMIRA scenarios are missing")
    rows = [r for r in power.get("rows", []) if float(r.get("effect_r_scenario", -1)) == 0.293]
    if len(rows) != 1:
        raise ValueError("expected exactly one r=0.293 power scenario")
    power_row = rows[0]
    needed = int(power_row["required_blocks"]["power_0.8"]["adjusted_alpha"])
    observed = int(power.get("observed_blocks", 0))
    current_cells = int(expansion["current_occupied_5deg_cells"])
    broadest_cells = int(broadest["total_5deg_cells_occupied"])
    broadest_new_cells = int(broadest["additional_5deg_cells_occupied"])
    if broadest_cells < current_cells:
        raise ValueError("broadest candidate pool cannot occupy fewer cells than current sample")

    return {
        "schema": "astrafek-icl-camira-expansion-feasibility-v1",
        "status": "complete_design_diagnostic",
        "science_measurement": False,
        "cas_query_submitted": False,
        "coordinates_exported": False,
        "images_downloaded": False,
        "observed_independent_blocks": observed,
        "exploratory_effect_r_scenario": 0.293,
        "power_target": 0.8,
        "n_comparisons": int(power["n_comparisons"]),
        "required_independent_blocks": needed,
        "block_shortfall_from_observed": max(0, needed - observed),
        "current_sample_occupied_5deg_cells": current_cells,
        "current_selection_additional_candidates": int(current["greedy_additions_after_6arcmin_separation"]),
        "current_selection_new_5deg_cells": int(current["additional_5deg_cells_occupied"]),
        "broadest_pool_additional_candidates": int(broadest["greedy_additions_after_6arcmin_separation"]),
        "broadest_pool_new_target_5deg_cells": broadest_new_cells,
        "broadest_pool_total_5deg_cells": broadest_cells,
        "broadest_pool_target_cell_upper_bound_on_connected_blocks": broadest_cells,
        "minimum_block_gap_even_at_target_cell_upper_bound": max(0, needed - broadest_cells),
        "input_hashes": {},
        "limitations": [
            "The r=0.293 effect is exploratory and may be inflated by post-selection.",
            "Fisher-z power is a design approximation, not the connected-block bootstrap or maxT analysis.",
            "Under pair_block_labels, rows sharing any target cell are joined; therefore distinct target cells upper-bound the number of connected pair blocks, while control cells may merge blocks but cannot create additional ones.",
            "The upper bound assumes the expansion pool and final pairs use the same 5-degree target-cell and pair-connectivity definition; it is not a count of usable or statistically independent blocks.",
            "The CAMIRA expansion scenarios are catalog-pool diagnostics; new targets were not selected, queried, or checked for images, masks, depth, PSF, or member completeness.",
            "The assessment does not authorize any HSC data transfer.",
        ],
    }


def render_report(result: dict, power_path: Path, expansion_path: Path,
                  power_sha: str, expansion_sha: str) -> str:
    return f"""# ICL–위성 연관성: CAMIRA 표본 확장 가능성 — 2026-10-20

## 결과

기존 탐색 효과크기 시나리오 `r={result['exploratory_effect_r_scenario']:.3f}`와 30개 비교의 Bonferroni 계획에서 80% 검정력에 필요한 공간 블록은 **{result['required_independent_blocks']}개**다. 현재 분석은 {result['observed_independent_blocks']}개로, 계획상 부족분은 {result['block_shortfall_from_observed']}개다.

현재 CAMIRA 선택을 그대로 적용하면 최소 6 arcmin 분리 후 후보가 {result['current_selection_additional_candidates']}개 더 있지만 새 5도 셀은 {result['current_selection_new_5deg_cells']}개다. 적색이동을 0.20–0.70으로 넓히고 BCG 분광 일치 조건까지 제거한 가장 큰 오프라인 후보 시나리오는 {result['broadest_pool_additional_candidates']}개 후보를 더하고 새 target 5도 셀 {result['broadest_pool_new_target_5deg_cells']}개, 총 {result['broadest_pool_total_5deg_cells']}개를 차지한다. 현재 pair-block 정의에서는 target 5도 셀을 공유하는 행들을 반드시 같은 연결요소로 합치므로, 각 target 셀이 서로 다른 블록을 만든다고 가정한 상한도 {result['broadest_pool_target_cell_upper_bound_on_connected_blocks']}개다. 이는 180개 계획치보다 최소 {result['minimum_block_gap_even_at_target_cell_upper_bound']}개 적다. matched-control 위치와 영상 지원을 적용하면 블록이 더 합쳐지거나 표본이 제외될 수 있어 실제 유효 블록 수는 상한 이하이며, 이 계산은 독립성 자체를 증명하지 않는다.

따라서 현재 S20A CAMIRA 목록에서 표본 수만 늘리는 것으로는 기존 검정력 계획을 충족할 수 없다. 확장 연구는 더 넓은 독립 HSC 하늘영역의 사전 정의된 군집 표본을 찾아야 한다. PDR3 Wide의 명목 설계 면적은 약 1,200 deg²이고 PDR3 공개자료 페이지는 명목 깊이의 다중대역 자료가 600 deg²를 넘는다고 안내한다. 이 면적 설명은 후보별 영상/필터 지원을 증명하지 않으며, 실제 쿼리나 자료 반출도 하지 않았다.

## 방법과 한계

Astrafek 스크립트가 두 기존 JSON 산출물의 스키마·상태·과학 게이트를 확인한다. pair-block 구현은 같은 target 5도 셀을 공유하는 모든 행을 연결하므로 target 셀 개수가 연결 블록 수의 상한이다. power 계획은 이론적 Fisher-z 근사이며 `r=0.293`은 사후 탐색값이다. 이 상한은 독립성을 뜻하지 않는다. 새 표본의 redshift, HSC 영상 지원, target/control footprint, PSF, 깊이, sky, 위성 완전도와 공분산은 전부 미검증이다. 모든 결과는 `science_measurement=false`다.

PDR3 공식 설명: [PDR3 survey area](https://hsc-release.mtk.nao.ac.jp/doc/index.php/survey__pdr3/), [PDR3 available data](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/).

## 재현성

- 스크립트: `assess_icl_catalog_expansion_feasibility.py`
- power 입력: `{power_path}` (SHA-256 `{power_sha}`)
- CAMIRA 후보 입력: `{expansion_path}` (SHA-256 `{expansion_sha}`)
- 결과 JSON: `feasibility.json`
- 요청/다운로드/좌표 반출: 없음
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power", type=Path, default=DEFAULT_POWER)
    parser.add_argument("--expansion", type=Path, default=DEFAULT_EXPANSION)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"refusing to overwrite existing output: {work}")
    power_path, expansion_path = args.power.resolve(), args.expansion.resolve()
    power = read_summary(power_path, "astrafek-icl-spatial-block-power-plan-v1")
    expansion = read_summary(expansion_path, "camira-offline-expansion-pool-audit-v1")
    result = assess(power, expansion)
    result["input_hashes"] = {str(power_path): sha256(power_path),
                              str(expansion_path): sha256(expansion_path)}
    result["script_sha256"] = sha256(Path(__file__).resolve())
    work.mkdir(parents=True, exist_ok=False)
    (work / "feasibility.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (work / "assessment.md").write_text(
        render_report(result, power_path, expansion_path,
                      result["input_hashes"][str(power_path)],
                      result["input_hashes"][str(expansion_path)]), encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "observed_independent_blocks", "required_independent_blocks",
        "broadest_pool_new_target_5deg_cells", "broadest_pool_target_cell_upper_bound_on_connected_blocks",
        "minimum_block_gap_even_at_target_cell_upper_bound",
        "science_measurement")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
