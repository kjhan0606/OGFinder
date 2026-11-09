#!/usr/bin/env python3
"""Joint sky-residual and coadd-covariance sensitivity for ICL--satellite pairs.

This offline analysis joins the fixed 74-pair ICL--satellite sample to existing
target/control outer-residual and local coadd-correlation diagnostics. It uses
sample-scoped connected 5-degree blocks and never rereads or modifies images.
The covariance statistic is only a nuisance proxy, not estimator covariance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_icl_target_control_outer_sky import partial_pearson, read_json, read_tsv, sha256
from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ASSOCIATION = DATA / "icl_outer_sky_association_reblocked_2026-11-03_01"
OUTER = DATA / "icl_target_control_outer_sky_2026-11-03_01"
COVARIANCE = DATA / "icl_coadd_noise_correlation_full_2026-11-03_01"
DEFAULT_OUTPUT = DATA / "icl_satellite_joint_sky_covariance_2026-11-03_01"


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty table")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def load_joined_rows() -> tuple[dict[str, list[dict]], dict[str, str]]:
    assoc_summary_path = ASSOCIATION / "summary.json"
    assoc_rows_path = ASSOCIATION / "association_rows.tsv"
    assoc_hash_path = ASSOCIATION / "output_hashes.json"
    outer_summary_path = OUTER / "summary.json"
    outer_rows_path = OUTER / "paired_outer_sky_diagnostics.tsv"
    cov_summary_path = COVARIANCE / "summary.json"
    cov_audit_path = COVARIANCE / "independent_audit.json"
    cov_rows_path = COVARIANCE / "frame_noise_correlation.tsv"
    for path in (assoc_summary_path, assoc_rows_path, assoc_hash_path,
                 outer_summary_path, outer_rows_path, cov_summary_path,
                 cov_audit_path, cov_rows_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    assoc_summary = read_json(assoc_summary_path)
    outer_summary = read_json(outer_summary_path)
    cov_summary = read_json(cov_summary_path)
    cov_audit = read_json(cov_audit_path)
    assoc_hashes = read_json(assoc_hash_path)
    if (assoc_summary.get("status") != "complete_exploratory_replay"
            or assoc_summary.get("science_measurement") is not False
            or assoc_summary.get("supersedes_source_summary_association_bootstrap") is not True
            or assoc_summary.get("sky_subtracted") is not False
            or assoc_hashes.get("association_rows_sha256") != sha256(assoc_rows_path)):
        raise ValueError("association replay is incomplete or its row hash is invalid")
    if (outer_summary.get("status") != "complete_diagnostic"
            or outer_summary.get("n_target_control_pairs") != 783
            or outer_summary.get("science_measurement") is not False
            or outer_summary.get("sky_subtracted") is not False):
        raise ValueError("outer-residual diagnostic is incomplete or outside its declared scope")
    if (cov_summary.get("status") != "complete_exploratory"
            or cov_summary.get("science_measurement") is not False
            or cov_audit.get("status") != "verified"
            or cov_audit.get("science_measurement") is not False
            or cov_audit.get("frame_and_lag_arithmetic_reproduced") is not True
            or cov_audit.get("input_acquisition_hashes_matched") is not True
            or cov_audit.get("frame_table_sha256") != sha256(cov_rows_path)):
        raise ValueError("coadd covariance diagnostic is incomplete or failed its independent audit")

    association_rows = read_tsv(assoc_rows_path)
    outer_rows = read_tsv(outer_rows_path)
    covariance_rows = read_tsv(cov_rows_path)
    if len(association_rows) != 74 or len(outer_rows) != 783 or len(covariance_rows) != 1566:
        raise ValueError("unexpected row count in a required input table")

    outer_by_key = {(r["target_id"], r["filter"]): r for r in outer_rows}
    if len(outer_by_key) != len(outer_rows):
        raise ValueError("duplicate target/filter key in outer-residual table")
    covariance_by_key = {(r["cluster_id"], r["role"], r["filter"]): r
                         for r in covariance_rows}
    if len(covariance_by_key) != len(covariance_rows):
        raise ValueError("duplicate cluster/role/filter key in covariance table")

    by_filter = {"i": [], "i2": []}
    seen = set()
    for source in association_rows:
        target_id, control_id, band = (source["target_id"], source["control_id"], source["filter"])
        key = (target_id, band)
        if key in seen or band not in by_filter:
            raise ValueError(f"duplicate or unexpected association key: {key}")
        seen.add(key)
        outer = outer_by_key.get(key)
        target_noise = covariance_by_key.get((target_id, "target", band))
        control_noise = covariance_by_key.get((target_id, "matched_control", band))
        if outer is None or target_noise is None or control_noise is None:
            raise ValueError(f"missing matched diagnostic for {target_id}/{band}")
        if outer["control_id"] != control_id:
            raise ValueError(f"outer-residual control identity mismatch for {target_id}/{band}")
        t_sigma = float(target_noise["top_hat_sigma_inflation_approx_square_lag_window"])
        c_sigma = float(control_noise["top_hat_sigma_inflation_approx_square_lag_window"])
        values = {
            "target_id": target_id,
            "control_id": control_id,
            "filter": band,
            "icl_delta_pixel_sfd_njy_arcsec2": float(source["icl_delta_pixel_sfd_njy_arcsec2"]),
            "satellite_delta_center_sfd_njy_arcsec2": float(source["satellite_delta_center_sfd_njy_arcsec2"]),
            "z": float(source["z"]),
            "richness": float(source["richness"]),
            "outer_plane_delta_njy_arcsec2": float(outer["target_minus_control_plane_constant_njy_arcsec2"]),
            "target_local_coadd_sigma_inflation_proxy": t_sigma,
            "control_local_coadd_sigma_inflation_proxy": c_sigma,
            "target_minus_control_local_coadd_sigma_inflation_proxy": t_sigma - c_sigma,
            "target_ra_deg": float(source["target_ra_deg"]),
            "target_dec_deg": float(source["target_dec_deg"]),
            "control_ra_deg": float(source["control_ra_deg"]),
            "control_dec_deg": float(source["control_dec_deg"]),
            "connected_5deg_block": source["connected_5deg_block"],
        }
        if not all(math.isfinite(v) for v in values.values() if isinstance(v, float)):
            raise ValueError(f"non-finite paired value for {target_id}/{band}")
        by_filter[band].append(values)
    if (len(by_filter["i2"]) != 53 or len(by_filter["i"]) != 21
            or len(seen) != 74):
        raise ValueError("fixed ICL--satellite association membership changed")
    inputs = {str(p.resolve()): sha256(p) for p in (
        assoc_summary_path, assoc_rows_path, assoc_hash_path,
        outer_summary_path, outer_rows_path, cov_summary_path, cov_audit_path, cov_rows_path,
    )}
    return by_filter, inputs


def summarize_filter(rows: list[dict], replicas: int, seed: int) -> dict:
    target = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    control = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        target, control, np.ones(len(rows), dtype=bool))
    saved_labels = np.asarray([r["connected_5deg_block"] for r in rows])
    if len(labels) != len(rows):
        raise ValueError("spatial block helper changed the fixed sample row count")
    # Labels are canonicalized by the helper; only the induced partition must agree.
    for i in range(len(rows)):
        for j in range(len(rows)):
            if (labels[i] == labels[j]) != (saved_labels[i] == saved_labels[j]):
                raise ValueError("recomputed association spatial blocks disagree with saved blocks")

    x = np.asarray([r["icl_delta_pixel_sfd_njy_arcsec2"] for r in rows])
    y = np.asarray([r["satellite_delta_center_sfd_njy_arcsec2"] for r in rows])
    covariates = {
        "z_richness": ["z", "richness"],
        "z_richness_outer_residual": ["z", "richness", "outer_plane_delta_njy_arcsec2"],
        "z_richness_coadd_proxy": ["z", "richness", "target_minus_control_local_coadd_sigma_inflation_proxy"],
        "z_richness_outer_and_coadd_proxy": ["z", "richness", "outer_plane_delta_njy_arcsec2",
                                              "target_minus_control_local_coadd_sigma_inflation_proxy"],
    }
    arrays = {key: np.asarray([r[key] for r in rows], dtype=float)
              for keys in covariates.values() for key in keys}
    points = {name: partial_pearson(x, y, [arrays[k] for k in keys])
              for name, keys in covariates.items()}
    blocks = np.unique(labels)
    by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = {name: np.full(replicas, np.nan) for name in covariates}
    for replica in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([by_block[b] for b in chosen])
        bx, by = x[indices], y[indices]
        for name, keys in covariates.items():
            draws[name][replica] = partial_pearson(
                bx, by, [arrays[k][indices] for k in keys])

    q = lambda values: ([float(v) for v in np.percentile(values, (2.5, 97.5))]
                        if len(values) else [None, None])
    baseline = "z_richness"
    metrics = {}
    base_draw = draws[baseline]
    for name in covariates:
        valid = draws[name][np.isfinite(draws[name])]
        paired = np.isfinite(base_draw) & np.isfinite(draws[name])
        delta = draws[name][paired] - base_draw[paired]
        metrics[name] = {
            "partial_pearson_r": points[name] if math.isfinite(points[name]) else None,
            "block_bootstrap_ci95": q(valid),
            "fraction_positive_valid_replicas": float(np.mean(valid > 0)) if len(valid) else None,
            "valid_replicas": int(len(valid)),
            "paired_delta_from_z_richness_r": 0.0 if name == baseline else points[name] - points[baseline],
            "paired_delta_from_z_richness_block_bootstrap_ci95": q(delta),
            "paired_valid_replicas": int(len(delta)),
        }
    return {
        "n_pairs": len(rows),
        "n_connected_blocks": int(n_blocks),
        "n_target_cells": int(n_target),
        "n_control_cells": int(n_control),
        "largest_block_pairs": int(largest),
        "bootstrap_replicas": replicas,
        "seed": seed,
        "models": metrics,
        "interpretation": "exploratory partial-correlation covariate sensitivity; the coadd proxy is not radial-estimator covariance and outer residual is not pure sky",
    }


def analyze(output: Path, replicas: int, seed: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    rows_by_filter, inputs = load_joined_rows()
    results = {
        band: summarize_filter(rows, replicas, seed + (0 if band == "i2" else 1))
        for band, rows in rows_by_filter.items()
    }
    output.mkdir(parents=True, exist_ok=False)
    saved = [row for band in ("i2", "i") for row in rows_by_filter[band]]
    row_path = output / "association_sky_covariance_rows.tsv"
    write_tsv(row_path, saved)
    summary = {
        "schema": "astrafek-icl-satellite-joint-sky-covariance-sensitivity-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "n_pairs": len(saved),
        "n_pairs_by_filter": {band: len(rows) for band, rows in rows_by_filter.items()},
        "method": "partial Pearson correlation between paired target-control ICL and satellite surface brightness; nuisance covariates are added cumulatively; resample association-sample connected 5-degree target/control blocks within each filter",
        "coadd_covariance_proxy": "per-frame approximate top-hat sigma-inflation factor from variance-normalized local 9x9 residual-lag sum; target minus matched-control value",
        "sky_proxy": "600-900 proper-kpc target-minus-control plane intercept; diagnostic only and not subtracted",
        "filters": results,
        "inputs": inputs,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "output_rows_sha256": sha256(row_path),
        "hsc_network_used": False,
        "images_reread": False,
        "data_requested": False,
        "sky_subtracted": False,
        "covariance_correction_applied": False,
        "limitations": [
            "The local lag-window covariance factor is not the exact covariance of the radial ICL estimator and is used only as a diagnostic covariate.",
            "The outer annulus may contain real cluster/ICL light, PSF wings, cirrus, and missed sources; no outer residual was subtracted.",
            "The ICL--satellite sample is fixed, selected, and exploratory; small filter-specific sample and block counts limit precision.",
            "This regression neither propagates full image-pixel covariance nor validates PSF wings, source completeness, extinction, or upstream HSC sky transfer.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       allow_nan=False, sort_keys=True) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({"status": summary["status"], "pairs": summary["n_pairs_by_filter"],
                      "filters": results, "science_measurement": False,
                      "summary_path": str(summary_path),
                      "summary_sha256": summary["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicas are required")
    analyze(args.output.resolve(), args.replicas, args.seed)


if __name__ == "__main__":
    main()
