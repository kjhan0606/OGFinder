#!/usr/bin/env python3
"""Test whether target-control HSC frame-count offsets track the ICL contrast.

This offline check joins the approved count-only frame audit to the existing
200--300 kpc ICL target-control profiles. It estimates partial correlations
within filter, with connected target/control sky-block bootstrap intervals.
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


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PROFILE_DIR = ROOT / "icl_annulus_support_2026-10-18_01"
PROFILE_PATH = PROFILE_DIR / "annulus_pair_profiles.tsv"
PROFILE_SUMMARY = PROFILE_DIR / "summary.json"
TARGET_MANIFEST = ROOT / "sample_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
COUNT_AUDIT_DIR = ROOT / "icl_frame_metadata_count_audit_2026-11-05_02"
COUNT_AUDIT_PATH = COUNT_AUDIT_DIR / "summary.json"
COUNT_TABLE = COUNT_AUDIT_DIR / "frame_counts_verified.tsv"
DEFAULT_OUTPUT = ROOT / "icl_frame_count_sensitivity_2026-11-05_01"
EXPECTED_COUNT_PLAN_SHA256 = "3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014"
COUNT_FIELDS = ("n_unique_frames", "n_unique_visits", "n_unique_ccds")
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(values.size), covariates])
    beta, _, rank, _ = np.linalg.lstsq(design, values, rcond=None)
    if rank < design.shape[1]:
        raise ValueError("control-variable design is rank deficient")
    return values - design @ beta


def partial_r(rows: list[dict], count_field: str) -> float:
    y = np.asarray([float(r["delta_sb_njy_arcsec2"]) for r in rows])
    x = np.asarray([float(r[f"delta_{count_field}"]) for r in rows])
    controls = np.column_stack([
        np.asarray([float(r["z"]) for r in rows]),
        np.log(np.asarray([float(r["richness"]) for r in rows])),
        np.asarray([float(r["common_coverage"]) for r in rows]),
    ])
    ry, rx = residualize(y, controls), residualize(x, controls)
    if np.std(ry) == 0 or np.std(rx) == 0:
        return float("nan")
    return float(np.corrcoef(ry, rx)[0, 1])


def bootstrap_partial_r(rows: list[dict], labels: np.ndarray, field: str,
                        replicas: int, seed: int) -> dict:
    value = partial_r(rows, field)
    groups = np.unique(labels)
    if len(rows) < 12 or len(groups) < 6:
        return {"status": "insufficient_support", "n_pairs": len(rows),
                "n_connected_blocks": len(groups), "partial_r": value}
    block_rows = {g: np.flatnonzero(labels == g) for g in groups}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        idx = np.concatenate([block_rows[g] for g in chosen])
        sample = [rows[int(i)] for i in idx]
        try:
            r = partial_r(sample, field)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        if np.isfinite(r):
            draws.append(r)
    valid = len(draws)
    result = {
        "status": "complete_conditional_diagnostic" if valid >= 0.8 * replicas else "low_valid_fraction",
        "estimand": f"partial Pearson r(delta ICL SB, delta {field}) controlling z, log richness, and common aperture coverage",
        "n_pairs": len(rows), "n_connected_5deg_blocks": len(groups),
        "partial_r": value, "replicas_requested": replicas,
        "replicas_valid": valid, "valid_fraction": valid / replicas,
        "seed": seed,
        "resampling_scope": "connected 5-degree target/control endpoint blocks, resampled within filter",
    }
    if valid:
        result["block_bootstrap_ci95"] = np.percentile(draws, [2.5, 97.5]).tolist()
        result["bootstrap_fraction_positive"] = float(np.mean(np.asarray(draws) > 0))
    return result


def load_inputs() -> tuple[list[dict], dict]:
    for path in (PROFILE_PATH, PROFILE_SUMMARY, TARGET_MANIFEST, CONTROL_MANIFEST,
                 COUNT_AUDIT_PATH, COUNT_TABLE):
        if not path.is_file():
            raise FileNotFoundError(path)
    profile_summary = json.loads(PROFILE_SUMMARY.read_text(encoding="utf-8"))
    if (profile_summary.get("status") != "complete"
            or profile_summary.get("science_gates", {}).get("science_measurement") is not False
            or profile_summary.get("n_200_300_pairs") != 224
            or profile_summary.get("output_sha256", {}).get(str(PROFILE_PATH)) != sha256(PROFILE_PATH)):
        raise ValueError("saved ICL annulus profile is incomplete or changed")
    count_audit = json.loads(COUNT_AUDIT_PATH.read_text(encoding="utf-8"))
    if (count_audit.get("status") != "complete_exploratory"
            or count_audit.get("query_plan_sha256") != EXPECTED_COUNT_PLAN_SHA256
            or count_audit.get("n_centers") != 1598
            or count_audit.get("centers_without_patch") != 0
            or count_audit.get("science_measurement") is not False
            or count_audit.get("verified_rows_tsv_sha256") != sha256(COUNT_TABLE)):
        raise ValueError("approved count-only frame audit is incomplete or changed")

    targets_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(r["id"]): r for r in targets_doc["clusters"]}
    controls = {str(r["target_cluster_id"]): r for r in controls_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798:
        raise ValueError("unexpected target/control manifests")
    count_rows = {r["center_id"]: r for r in read_tsv(COUNT_TABLE)}
    if len(count_rows) != 1598:
        raise ValueError("expected one audited count row for each of 1,598 centers")

    joined = []
    for source in read_tsv(PROFILE_PATH):
        if source["radial_bin_kpc"] != "200-300":
            continue
        ident, band = str(source["id"]), source["filter"]
        if ident not in targets or ident not in controls or band not in FILTERS:
            raise ValueError(f"profile pair lacks its original matched manifest: {ident}/{band}")
        target_count = count_rows.get(f"target:{ident}")
        control_id = str(source["control_id"])
        control_count = count_rows.get(f"control:{control_id}")
        if target_count is None or control_count is None:
            raise ValueError(f"missing target/control frame count: {ident}/{control_id}")
        for col in ("sample_role",):
            if target_count[col] != "target" or control_count[col] != "control":
                raise ValueError("target/control frame-count role mismatch")
        if (control_id != str(controls[ident]["id"])
                or str(controls[ident]["target_cluster_id"]) != ident):
            raise ValueError(f"control manifest mismatch for cluster {ident}")
        row = {
            "cluster_id": ident, "filter": band, "control_id": control_id,
            "ra_deg": float(source["ra_deg"]), "dec_deg": float(source["dec_deg"]),
            "z": float(source["z"]), "richness": float(source["richness"]),
            "delta_sb_njy_arcsec2": float(source["delta_sb_njy_arcsec2"]),
            "common_coverage": float(source["common_coverage"]),
            "target_ra_deg": float(targets[ident]["ra_deg"]),
            "target_dec_deg": float(targets[ident]["dec_deg"]),
            "control_ra_deg": float(controls[ident]["ra_deg"]),
            "control_dec_deg": float(controls[ident]["dec_deg"]),
        }
        for field in COUNT_FIELDS:
            row[f"target_{field}"] = int(target_count[field])
            row[f"control_{field}"] = int(control_count[field])
            row[f"delta_{field}"] = row[f"target_{field}"] - row[f"control_{field}"]
        joined.append(row)
    if len(joined) != 224 or len({(r["cluster_id"], r["filter"]) for r in joined}) != 224:
        raise ValueError("expected 224 unique 200--300 kpc ICL profile/count pairs")
    return joined, {"profile_summary": profile_summary, "count_audit": count_audit,
                    "targets": targets, "controls": controls}


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    rows, state = load_inputs()
    results, detailed_rows = {}, []
    for filter_index, band in enumerate(FILTERS):
        selected = [r for r in rows if r["filter"] == band]
        expected = 71 if band == "i" else 153
        if len(selected) != expected:
            raise ValueError(f"unexpected {band} supported-pair count: {len(selected)}")
        targets = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in selected]
        controls = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in selected]
        _, labels, nt, nc, nblocks, largest = pair_block_labels(
            targets, controls, np.ones(len(selected), dtype=bool))
        labels = np.asarray(labels, dtype=str)
        outputs = {}
        for field_index, field in enumerate(COUNT_FIELDS):
            outputs[field] = bootstrap_partial_r(
                selected, labels, field, replicas,
                seed + filter_index * 100 + field_index,
            )
        results[band] = {
            "n_pairs": len(selected), "n_connected_5deg_blocks": int(nblocks),
            "n_target_cells": int(nt), "n_control_cells": int(nc),
            "largest_block_pairs": int(largest),
            "frame_count_offset_test": outputs,
            "count_offset_distributions": {
                field: {
                    "median_target_minus_control": float(np.median([r[f"delta_{field}"] for r in selected])),
                    "mean_target_minus_control": float(np.mean([r[f"delta_{field}"] for r in selected])),
                    "fraction_equal": float(np.mean([r[f"delta_{field}"] == 0 for r in selected])),
                } for field in COUNT_FIELDS
            },
        }
        for r, label in zip(selected, labels, strict=True):
            detailed_rows.append({**r, "filter_block_label": label})

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "icl_frame_count_join_200_300.tsv"
    fields = list(detailed_rows[0])
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(detailed_rows)
    inputs = (PROFILE_PATH, PROFILE_SUMMARY, TARGET_MANIFEST, CONTROL_MANIFEST,
              COUNT_AUDIT_PATH, COUNT_TABLE)
    result = {
        "schema": "astrafek-icl-frame-count-sensitivity-v1",
        "status": "complete_offline_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "n_profile_pairs": len(rows),
        "outcome": "existing target-minus-matched-control ICL surface brightness in the 200-300 kpc annulus",
        "predictors": [f"target-minus-control {field}" for field in COUNT_FIELDS],
        "controls": ["cluster redshift", "log richness", "common ICL aperture coverage"],
        "estimator": "unweighted partial Pearson correlation, stratified by image filter",
        "block_method": "connected 5-degree target/control endpoint blocks resampled within filter",
        "replicas": replicas,
        "seed_base": seed,
        "by_filter": results,
        "joined_table": str(table_path),
        "joined_table_sha256": sha256(table_path),
        "inputs_sha256": {str(path): sha256(path) for path in inputs},
        "limitations": [
            "The 224 rows are the ICL-aperture-supported subset, not the full 800-center parent.",
            "The count query measures patch membership, not pixel-level coverage, exposure weighting, or variance.",
            "The 200-300 kpc ICL profile has its own support and source-mask selection; this is not a corrected science estimator.",
            "Block intervals are conditional on the saved 5-degree endpoint-block partition and omit coadd, sky, PSF, and selection uncertainty.",
            "No causal interpretation or science_measurement promotion is supported.",
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
        "by_filter": result["by_filter"],
        "joined_table_sha256": result["joined_table_sha256"],
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
