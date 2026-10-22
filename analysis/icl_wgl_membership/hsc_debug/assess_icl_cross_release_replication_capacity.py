#!/usr/bin/env python3
"""Audit how much new i2 sky-cell capacity S20A and S23B candidate pools add."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from plan_icl_spatial_block_power import required_units, two_sided_power
from prestack_core import spatial_labels


DATA = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE = DATA / "sample_800_2026-10-13"
DEFAULT_S20A_PLAN = SAMPLE / "icl_expansion_coverage_plan_2026-10-21_01"
DEFAULT_S23B_LEDGER = SAMPLE / "icl_s23b_spatial_candidate_ledger_2026-10-22_01"
DEFAULT_PARENT = SAMPLE / "sample_manifest.json"
DEFAULT_ASSOCIATION = SAMPLE / "satellite_outer_icl_association_perbin_2026-10-18_03/matched_cluster_values.tsv"
DEFAULT_OUTPUT = SAMPLE / "icl_cross_release_replication_capacity_2026-10-22_01"

EXPECTED = {
    "s20a_plan": "9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006",
    "s20a_centers": "03c6d049a39dd5b271e1cf201129301f1c76e0462010d2f67c4c123773a4344c",
    "s23b_summary": "f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1",
    "s23b_ledger": "ab8df9fce1435d1d9eca15f742f43356cbd430b2159852a14d9b67963d83b775",
    "parent": "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401",
    "association": "9a6befa87a81257f8faf1d34fa7c019a91d64484fbfa64136bd93496dfc7b58c",
}
EFFECTS = (0.20, 0.25, 0.293, 0.30, 0.40)
ALPHA = 0.05
N_CONTROLS = 2
N_COMPARISONS = 30


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def cell(ra: float, dec: float) -> str:
    return str(spatial_labels([ra], [dec], width=5.0)[0])


def summarize_cells(discovery: set[str], s20a: set[str], s23b: set[str]) -> dict:
    union = s20a | s23b
    pools = {"S20A_broad": s20a, "S23B_selection_compatible": s23b, "union_catalog_upper_bound": union}
    rows = []
    for name, cells in pools.items():
        new = cells - discovery
        rows.append({
            "pool": name,
            "n_unique_candidate_cells": len(cells),
            "n_cells_overlapping_i2_discovery": len(cells & discovery),
            "n_new_target_cells_vs_i2_discovery": len(new),
            "n_blocks_80pct_r293_single_endpoint": required_units(0.293, ALPHA, 0.80, N_CONTROLS),
            "power_r293_at_new_target_cell_cap_single_endpoint": two_sided_power(
                len(new), 0.293, ALPHA, N_CONTROLS),
            "shortfall_r293_single_endpoint": max(
                0, required_units(0.293, ALPHA, 0.80, N_CONTROLS) - len(new)),
            "n_blocks_80pct_r293_bonferroni_30": required_units(
                0.293, ALPHA / N_COMPARISONS, 0.80, N_CONTROLS),
        })
    return {
        "discovery_i2_target_cells": len(discovery),
        "s20a_unique_cells_added_beyond_s23b": len(s20a - s23b),
        "s20a_incremental_new_cells_vs_discovery_beyond_s23b": len((s20a - s23b) - discovery),
        "pools": rows,
        "effect_sensitivity": [
            {
                "assumed_partial_r": effect,
                "required_blocks_80pct_single_endpoint": required_units(effect, ALPHA, .8, N_CONTROLS),
                "power_at_union_new_cell_cap": two_sided_power(
                    len(union - discovery), effect, ALPHA, N_CONTROLS),
                "shortfall_at_union_new_cell_cap": max(
                    0, required_units(effect, ALPHA, .8, N_CONTROLS) - len(union - discovery)),
                "required_blocks_80pct_bonferroni_30": required_units(
                    effect, ALPHA / N_COMPARISONS, .8, N_CONTROLS),
            }
            for effect in EFFECTS
        ],
        "union_new_cell_count_is_only_an_upper_bound": True,
    }


def load_inputs(s20a_plan_dir: Path, s23b_dir: Path, parent_path: Path,
                association_path: Path) -> tuple[dict, dict]:
    s20a_plan_path = s20a_plan_dir / "coverage_query_plan.json"
    s20a_centers_path = s20a_plan_dir / "candidate_centers.tsv"
    s23b_summary_path = s23b_dir / "summary.json"
    s23b_ledger_path = s23b_dir / "candidate_ledger.tsv"
    inputs = {
        "s20a_plan": s20a_plan_path, "s20a_centers": s20a_centers_path,
        "s23b_summary": s23b_summary_path, "s23b_ledger": s23b_ledger_path,
        "parent": parent_path, "association": association_path,
    }
    for key, path in inputs.items():
        if not path.is_file() or sha256(path) != EXPECTED[key]:
            raise ValueError(f"missing or changed reviewed input: {key}")
    s20a_plan = json.loads(s20a_plan_path.read_text(encoding="utf-8"))
    s23b_summary = json.loads(s23b_summary_path.read_text(encoding="utf-8"))
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    if (s20a_plan.get("status") != "prepared_not_submitted"
            or s20a_plan.get("cas_query_submitted") is not False
            or s20a_plan.get("n_candidates") != 2596
            or s23b_summary.get("status") != "offline_coordinate_and_footprint_ledger_complete"
            or s23b_summary.get("science_measurement") is not False
            or parent.get("science_measurement") is not False):
        raise ValueError("candidate-pool scope or status changed")
    manifest_cells = {str(row["id"]): cell(float(row["ra_deg"]), float(row["dec_deg"]))
                      for row in parent.get("clusters", [])}
    if len(manifest_cells) != 800:
        raise ValueError("current parent manifest is not the reviewed 800-object sample")
    association = read_tsv(association_path)
    i2_ids = {str(row["cluster_id"]) for row in association if row["filter"] == "i2"}
    if len(i2_ids) != 68 or not i2_ids.issubset(manifest_cells):
        raise ValueError("i2 discovery association does not match the reviewed 68 rows")
    discovery = {manifest_cells[ident] for ident in i2_ids}

    s20a_rows = read_tsv(s20a_centers_path)
    if len(s20a_rows) != 2596:
        raise ValueError("S20A candidate center inventory changed")
    s20a_cells = {cell(float(row["ra_deg"]), float(row["dec_deg"])) for row in s20a_rows}
    s23b_rows = read_tsv(s23b_ledger_path)
    if len(s23b_rows) != 1333:
        raise ValueError("S23B candidate ledger row count changed")
    s23b_cells = set()
    for row in s23b_rows:
        actual = cell(float(row["ra_deg"]), float(row["dec_deg"]))
        if actual != row["spatial_cell_5deg"]:
            raise ValueError("S23B candidate ledger cell label differs from shared helper")
        s23b_cells.add(actual)
    if len(s23b_cells) != 100:
        raise ValueError("S23B candidate sky-cell count changed")
    statistics = summarize_cells(discovery, s20a_cells, s23b_cells)
    return statistics, {key: sha256(path) for key, path in inputs.items()}


def render(statistics: dict, hashes: dict) -> str:
    pools = statistics["pools"]
    lines = [
        "# S20A–S23B 후보 풀의 독립 footprint 상한 — 2026-10-22",
        "",
        "기존 i2 discovery target은 33개 5도 셀에 있다. S23B의 selection-compatible 후보는 67개 새 셀을 제공한다. 기존에 승인되었으나 아직 제출되지 않은 S20A 계획은 더 넓은 `0.20 ≤ z_cl < 0.70` cluster-only pool이며, 그 후보 2,596개는 51개 새 셀을 제공하지만 S23B와 합쳐도 **69개**뿐이다. S20A 후보 셀 중 S23B보다 추가하는 셀은 2개다.",
        "",
        "| 후보 풀 | 고유 후보 셀 | i2 discovery 중첩 | 새 target 셀 상한 | r=0.293 단일검정 power 상한 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in pools:
        lines.append(
            f"| {row['pool']} | {row['n_unique_candidate_cells']} "
            f"| {row['n_cells_overlapping_i2_discovery']} "
            f"| {row['n_new_target_cells_vs_i2_discovery']} "
            f"| {row['power_r293_at_new_target_cell_cap_single_endpoint']:.3f} |"
        )
    union = pools[-1]
    lines.extend([
        "",
        f"union 69-cell 상한은 단일 고정 endpoint에서 r=0.293의 80% 목표 92 blocks에 {union['shortfall_r293_single_endpoint']}개 부족하다. 30개 비교 보정 목표는 180 blocks이므로 현재 discovery 33셀까지 단순 더해도 상한 102에서 78개가 모자란다. r=0.20의 단일 endpoint 목표는 196 blocks이고 69-cell power 상한은 {statistics['effect_sensitivity'][0]['power_at_union_new_cell_cap']:.3f}다.",
        "",
        "S20A cluster-only redshift 선택은 기존의 BCG 분광 일치 selection과 다르므로 이 union은 동일 과학 표본이 아니다. 표의 셀은 image/mask/PSF/satellite/control gate 전의 거친 target-cell ceiling이며 실제 independent connected blocks를 보장하지 않는다. 따라서 이전 승인 S20A metadata query를 수행해도 독립 footprint 부족을 실질적으로 해결하지는 않는다. 더 넓은 별도 sky/catalog footprint가 필요하다.",
        "",
        "조회·다운로드·pixel 측정은 없고 `science_measurement=false`다. 입력 SHA-256:",
    ])
    lines.extend(f"- `{name}`: `{digest}`" for name, digest in sorted(hashes.items()))
    lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s20a-plan-dir", type=Path, default=DEFAULT_S20A_PLAN)
    parser.add_argument("--s23b-dir", type=Path, default=DEFAULT_S23B_LEDGER)
    parser.add_argument("--parent", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    stats, hashes = load_inputs(args.s20a_plan_dir.resolve(), args.s23b_dir.resolve(),
                                args.parent.resolve(), args.association.resolve())
    result = {
        "schema": "astrafek-icl-cross-release-replication-capacity-v1",
        "status": "offline_catalog_cell_upper_bound",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "statistics": stats,
        "input_sha256": hashes,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "S20A and S23B candidate pools use different redshift/spectroscopic selections.",
            "5-degree rectangular target cells are only an upper bound; controls can merge connected blocks and support cuts reduce the usable sample.",
            "The existing S20A coverage plan is prepared and approved but has no returned metadata; no coverage was inferred from its SQL positions.",
            "No HSC request, download, image-pixel measurement, or scientific inference was performed.",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "assessment.md").write_text(render(stats, hashes), encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "s20a_new_cells": stats["pools"][0]["n_new_target_cells_vs_i2_discovery"],
        "s23b_new_cells": stats["pools"][1]["n_new_target_cells_vs_i2_discovery"],
        "union_new_cells": stats["pools"][2]["n_new_target_cells_vs_i2_discovery"],
        "s20a_incremental_cells_beyond_s23b": stats["s20a_incremental_new_cells_vs_discovery_beyond_s23b"],
        "science_measurement": False, "output_dir": str(output),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
