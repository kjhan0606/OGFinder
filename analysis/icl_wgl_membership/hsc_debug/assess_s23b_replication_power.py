#!/usr/bin/env python3
"""Compare S23B's optimistic new-sky-cell cap with a frozen ICL replication design."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from plan_icl_spatial_block_power import required_units, two_sided_power


DEFAULT_LEDGER = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "icl_s23b_spatial_candidate_ledger_2026-10-22_01/summary.json"
)
DEFAULT_OUTPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "icl_s23b_replication_power_2026-10-22_01"
)
EXPECTED_LEDGER_SHA256 = "f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1"
EFFECTS = (0.20, 0.25, 0.293, 0.30, 0.40)
ALPHA = 0.05
N_CONTROLS = 2  # redshift and richness in the partial correlation
N_EXPLORATORY_COMPARISONS = 30


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def evaluate(ledger: dict) -> dict:
    if (ledger.get("schema") != "astrafek-s23b-spatial-candidate-ledger-v1"
            or ledger.get("status") != "offline_coordinate_and_footprint_ledger_complete"
            or ledger.get("science_measurement") is not False):
        raise ValueError("input must be the completed exploratory S23B spatial ledger")
    stats = ledger.get("statistics", {})
    if (stats.get("candidate_count") != 1333
            or stats.get("association_unique_target_cells", {}).get("i2") != 33
            or stats.get("candidate_cells_new_vs_association", {}).get("i2") != 67
            or stats.get("candidate_cell_count_is_only_a_target_footprint_upper_bound") is not True):
        raise ValueError("S23B target-cell upper-bound inputs differ from the reviewed ledger")
    current_blocks = int(stats["association_unique_target_cells"]["i2"])
    new_cap = int(stats["candidate_cells_new_vs_association"]["i2"])
    rows = []
    for effect in EFFECTS:
        required_one = required_units(effect, ALPHA, 0.80, N_CONTROLS)
        required_family = required_units(
            effect, ALPHA / N_EXPLORATORY_COMPARISONS, 0.80, N_CONTROLS)
        rows.append({
            "assumed_partial_r": effect,
            "required_blocks_80pct_one_frozen_endpoint_alpha_0p05": required_one,
            "nominal_power_at_s23b_new_target_cell_cap": two_sided_power(
                new_cap, effect, ALPHA, N_CONTROLS),
            "additional_blocks_needed_beyond_s23b_cap": max(0, required_one - new_cap),
            "required_blocks_80pct_bonferroni_30_tests": required_family,
            "additional_blocks_needed_beyond_current_plus_s23b_target_cell_cap": max(
                0, required_family - current_blocks - new_cap),
        })
    return {
        "schema": "astrafek-s23b-icl-replication-power-gate-v1",
        "status": "design_capacity_assessment",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "current_i2_discovery_target_cell_count": current_blocks,
        "s23b_new_i2_target_cell_upper_bound": new_cap,
        "combined_target_cell_upper_bound": current_blocks + new_cap,
        "test_design": {
            "primary_endpoint_candidate": "i2; ICL target-minus-control mean at 200-300 kpc versus satellite target-minus-control mean at 50-75 kpc",
            "statistic": "two-sided partial Pearson correlation adjusted for redshift and richness",
            "alpha": ALPHA,
            "exploratory_family_comparisons": N_EXPLORATORY_COMPARISONS,
            "independent_unit": "connected target/control spatial block, not candidate count or target cell",
            "pre_registration_status": "proposed; must be frozen before inspecting S23B outcomes",
        },
        "scenarios": rows,
        "limitations": [
            "The 67 S23B cells are an optimistic target-cell ceiling; image support, masks, depth, PSF, satellite completeness, controls, and connected pair-block merges can only reduce the usable count.",
            "The 33 existing i2 cells belong to the exploratory discovery sample. Combining them with S23B is not an independent replication.",
            "The r=0.293 scenario was selected after looking at the original sample and may be inflated; smaller effects need substantially more blocks.",
            "Fisher-z power is an approximation and does not replace the planned connected-block bootstrap or maxT procedure.",
            "The candidate primary endpoint is proposed for an S23B-only replication and is not yet preregistered or measured.",
        ],
    }


def render(result: dict, ledger_path: Path, ledger_sha: str) -> str:
    lines = [
        "# S23B ICL–위성 독립 재현 표본력 관문 — 2026-10-22",
        "",
        "## 설계 결과",
        "",
        f"현재 탐색 i2 연관표는 target 셀 {result['current_i2_discovery_target_cell_count']}개다. S23B 후보표가 더할 수 있는 새 i2 target 셀은 최대 {result['s23b_new_i2_target_cell_upper_bound']}개여서 합산 상한은 {result['combined_target_cell_upper_bound']}개다. 새 표본은 기존 discovery와 target 셀 기준으로 분리되지만, 실제 matched-control 위치가 같은 공간 블록을 합칠 수 있다.",
        "",
        "향후 S23B-only 재현의 후보 primary endpoint는 i2에서 ICL 200–300 kpc target-minus-control와 위성 50–75 kpc target-minus-control의 부분 Pearson 상관(적색이동·richness 조정), 양측 α=0.05 단일 검정이다. 이 endpoint는 현재 결과를 본 뒤 제안된 것이므로, S23B 결과를 열람하기 전에 고정해야 독립 재현으로 해석할 수 있다.",
        "",
        "| 가정 partial r | 단일 고정 endpoint 80% 검정력 블록 수 | S23B 최대 67블록에서 power | S23B 상한 이후 최소 추가 블록 | 30비교 보정 80% 블록 수 | 현재+S23B target 셀 상한 대비 부족 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["scenarios"]:
        lines.append(
            f"| {row['assumed_partial_r']:.3f} "
            f"| {row['required_blocks_80pct_one_frozen_endpoint_alpha_0p05']} "
            f"| {row['nominal_power_at_s23b_new_target_cell_cap']:.3f} "
            f"| {row['additional_blocks_needed_beyond_s23b_cap']} "
            f"| {row['required_blocks_80pct_bonferroni_30_tests']} "
            f"| {row['additional_blocks_needed_beyond_current_plus_s23b_target_cell_cap']} |"
        )
    lines.extend([
        "",
        "사후 시나리오 r=0.293을 그대로 가정해도 독립 S23B 최대 상한 67블록의 Fisher-z power는 약 0.662로 80%에 못 미친다. 단일 검정 80% 목표는 92 blocks라 최소 25개 더 필요하고, r=0.20이면 196 blocks라 최소 129개가 더 필요하다. 실제 support/대조장/연결요소 손실 전의 상한 계산이므로 새 표본 규모 목표는 이보다 커야 한다. 기존 33 discovery 셀을 더한 100은 독립 재현 표본으로 사용할 수 없다.",
        "",
        "이것은 HSC coverage query나 과학 측정이 아니다. `science_measurement=false`이며 S23B 결과를 보기 전 단일 endpoint와 제외 규칙을 고정하고, 충분한 독립 footprint를 더 확보해야 한다.",
        "",
        f"입력 장부: `{ledger_path}` (SHA-256 `{ledger_sha}`).",
        "",
    ])
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    ledger_path, output = args.ledger.resolve(), args.output_dir.resolve()
    if not ledger_path.is_file() or sha256(ledger_path) != EXPECTED_LEDGER_SHA256:
        parser.error("S23B ledger missing or differs from the pinned reviewed input")
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    result = evaluate(ledger)
    result["input_ledger"] = str(ledger_path)
    result["input_ledger_sha256"] = sha256(ledger_path)
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "assessment.md").write_text(
        render(result, ledger_path, result["input_ledger_sha256"]), encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "s23b_independent_target_cell_cap": result["s23b_new_i2_target_cell_upper_bound"],
        "combined_target_cell_upper_bound": result["combined_target_cell_upper_bound"],
        "r293_single_test_power_at_cap": result["scenarios"][2]["nominal_power_at_s23b_new_target_cell_cap"],
        "r293_blocks_needed_one_test": result["scenarios"][2]["required_blocks_80pct_one_frozen_endpoint_alpha_0p05"],
        "output_dir": str(output),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
