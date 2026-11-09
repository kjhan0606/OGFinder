#!/usr/bin/env python3
"""Compare the 94-pair ICL--satellite association by central-mask eligibility.

This is an offline sample-sensitivity analysis. It uses saved pair summaries
and target/control coordinates only, with separate connected-sky-block
bootstrap intervals for each selected subgroup.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PAIR_TABLE = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03/matched_cluster_values.tsv"
SAMPLE_MANIFEST = ROOT / "sample_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_satellite_center_mask_sensitivity_2026-11-05_02"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_pairs() -> list[dict]:
    with PAIR_TABLE.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 94 or len({(r["cluster_id"], r["filter"]) for r in rows}) != 94:
        raise ValueError("expected 94 unique saved ICL/satellite support pairs")
    return rows


def pair_row(source: dict, targets: dict, controls: dict) -> dict:
    ident = str(source["cluster_id"])
    if ident not in targets or ident not in controls:
        raise ValueError(f"missing target/control manifest row for cluster {ident}")
    target, control = targets[ident], controls[ident]
    return {
        "cluster_id": ident,
        "filter": str(source["filter"]),
        "icl": float(source["icl_delta_200_300"]),
        "satellite": float(source["satellite_delta_50_75"]),
        "redshift": float(source["redshift"]),
        "richness": float(source["richness"]),
        "center_eligible": str(source["center_eligible"]).lower() == "true",
        "target_ra_deg": float(target["ra_deg"]),
        "target_dec_deg": float(target["dec_deg"]),
        "control_ra_deg": float(control["ra_deg"]),
        "control_dec_deg": float(control["dec_deg"]),
    }


def assoc(rows: list[dict]) -> float:
    return partial_pearson(
        [r["icl"] for r in rows], [r["satellite"] for r in rows],
        [[r["redshift"] for r in rows], [r["richness"] for r in rows]],
    )


def subset_bootstrap(rows: list[dict], replicas: int, seed: int) -> dict:
    if len(rows) < 6:
        return {"status": "insufficient_pairs", "n_pairs": len(rows)}
    targets = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    controls = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, nt, nc, nblocks, largest = pair_block_labels(
        targets, controls, np.ones(len(rows), dtype=bool))
    if nblocks < 6:
        return {
            "status": "insufficient_blocks", "n_pairs": len(rows),
            "n_connected_5deg_blocks": int(nblocks),
            "n_target_cells": int(nt), "n_control_cells": int(nc),
            "largest_block_pairs": int(largest), "partial_r": float(assoc(rows)),
        }
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        idx = np.concatenate([by_block[group] for group in chosen])
        try:
            value = assoc([rows[int(i)] for i in idx])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        if np.isfinite(value):
            draws.append(value)
    if len(draws) < 0.9 * replicas:
        raise ValueError(f"too few valid bootstrap draws: {len(draws)}/{replicas}")
    return {
        "status": "complete_conditional_diagnostic",
        "n_pairs": len(rows), "partial_r_z_richness": float(assoc(rows)),
        "block_bootstrap_ci95": np.percentile(draws, [2.5, 97.5]).tolist(),
        "bootstrap_fraction_positive": float(np.mean(np.asarray(draws) > 0)),
        "n_connected_5deg_blocks": int(nblocks),
        "n_target_cells": int(nt), "n_control_cells": int(nc),
        "largest_block_pairs": int(largest),
        "replicas_requested": replicas, "replicas_valid": len(draws), "seed": seed,
        "resampling_scope": "connected 5-degree target/control endpoint blocks, recomputed within this subgroup",
    }


def paired_eligibility_difference(rows: list[dict], replicas: int, seed: int) -> dict:
    eligible = [r for r in rows if r["center_eligible"]]
    ineligible = [r for r in rows if not r["center_eligible"]]
    if len(eligible) < 6 or len(ineligible) < 6:
        return {"status": "insufficient_pairs", "eligible": len(eligible),
                "ineligible": len(ineligible)}
    targets = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    controls = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, _, _, nblocks, _ = pair_block_labels(
        targets, controls, np.ones(len(rows), dtype=bool))
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    by_block = {g: np.flatnonzero(labels == g) for g in groups}
    eligible_flag = np.asarray([r["center_eligible"] for r in rows], dtype=bool)
    rng = np.random.default_rng(seed)
    differences = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        idx = np.concatenate([by_block[g] for g in chosen])
        idx_eligible = idx[eligible_flag[idx]]
        idx_ineligible = idx[~eligible_flag[idx]]
        if len(idx_eligible) < 6 or len(idx_ineligible) < 6:
            continue
        try:
            difference = assoc([rows[int(i)] for i in idx_eligible]) - assoc(
                [rows[int(i)] for i in idx_ineligible])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        if np.isfinite(difference):
            differences.append(float(difference))
    point = float(assoc(eligible) - assoc(ineligible))
    result = {
        "status": "complete_paired_block_sensitivity" if len(differences) >= 0.8 * replicas
                  else "low_valid_fraction",
        "estimand": "partial r(center-mask eligible) minus partial r(ineligible)",
        "point_difference": point,
        "n_connected_5deg_blocks": int(nblocks),
        "replicas_requested": replicas, "replicas_valid": len(differences),
        "valid_fraction": len(differences) / replicas,
        "seed": seed,
        "resampling_scope": "same resampled connected 5-degree blocks for both disjoint central-mask subgroups",
    }
    if len(differences) >= 0.8 * replicas:
        result["paired_difference_ci95"] = np.percentile(differences, [2.5, 97.5]).tolist()
        result["bootstrap_fraction_positive"] = float(np.mean(np.asarray(differences) > 0))
    return result


def run(output: Path, replicas: int, seed: int) -> dict:
    for p in (PAIR_TABLE, SAMPLE_MANIFEST, CONTROL_MANIFEST):
        if not p.is_file():
            raise FileNotFoundError(p)
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    targets_doc = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(r["id"]): r for r in targets_doc["clusters"]}
    controls = {str(r["target_cluster_id"]): r for r in controls_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798:
        raise ValueError("unexpected target/control manifest sizes")

    support = [pair_row(r, targets, controls) for r in read_pairs()]
    output_rows, by_filter = [], {}
    for band in ("i", "i2"):
        rows = [r for r in support if r["filter"] == band]
        expected = 26 if band == "i" else 68
        if len(rows) != expected:
            raise ValueError(f"unexpected {band} support count")
        eligible = [r for r in rows if r["center_eligible"]]
        masked = [r for r in rows if not r["center_eligible"]]
        if len(eligible) != (8 if band == "i" else 15):
            raise ValueError(f"unexpected {band} center-eligible count")
        by_filter[band] = {
            "all_icl_satellite_supported": subset_bootstrap(rows, replicas, seed + (0 if band == "i" else 10)),
            "center_mask_eligible": subset_bootstrap(eligible, replicas, seed + (1 if band == "i" else 11)),
            "center_mask_ineligible": subset_bootstrap(masked, replicas, seed + (2 if band == "i" else 12)),
            "paired_center_eligibility_difference": paired_eligibility_difference(
                rows, replicas, seed + (3 if band == "i" else 13)),
            "selection_counts": {"supported": len(rows), "eligible": len(eligible), "ineligible": len(masked)},
        }
        for label, selected in (("eligible", eligible), ("ineligible", masked)):
            for r in selected:
                output_rows.append({**r, "central_mask_group": label})

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "center_mask_subgroups.tsv"
    fields = list(output_rows[0])
    with pair_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    result = {
        "schema": "astrafek-icl-satellite-center-mask-sensitivity-v1",
        "status": "complete_offline_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "n_pairs": len(support),
        "estimator": "unweighted partial Pearson correlation of raw ICL 200-300 kpc and satellite 50-75 kpc contrasts, controlling linearly for cluster redshift and richness",
        "groups_by_filter": by_filter,
        "block_method": "5-degree connected components of target/control endpoint cells; recomputed within each subgroup and filter",
        "replicas": replicas,
        "seed_base": seed,
        "inputs_sha256": {str(p): sha256(p) for p in (PAIR_TABLE, SAMPLE_MANIFEST, CONTROL_MANIFEST)},
        "subgroup_table": str(pair_path),
        "subgroup_table_sha256": sha256(pair_path),
        "limitations": [
            "The 8 i-band and 15 i2-band center-mask-eligible groups are small and are not representative samples.",
            "The mask flag is a selection diagnostic, not a calibrated completeness or flux correction.",
            "The conditional block intervals omit PSF-wing, sky, extinction, coadd covariance, and selection-model uncertainty.",
            "No causal or full-800-parent inference is supported; science_measurement remains false.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(summary_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicas are required")
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "groups_by_filter": result["groups_by_filter"],
        "output_summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
