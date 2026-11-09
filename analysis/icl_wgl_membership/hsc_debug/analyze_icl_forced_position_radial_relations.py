#!/usr/bin/env python3
"""Test ICL trends in forced, unmatched CAMIRA-position light by radial bin.

This is a sensitivity analysis of forced-position profiles measured from the
existing HSC cutouts. Positions lack unique SEP detections, so the response is
called unmatched-CAMIRA-position light rather than confirmed satellite light.
One shared connected-sky-block Webb wild bootstrap controls the 14
filter-by-radius continuous ICL slopes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from evaluate_prestack_pair_support import pair_block_labels
from familywise_block_wild_association import WEBB_WEIGHTS, holm_adjust

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
FORCED = DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_both_2026-11-06_01"
ICL = DATA / "icl_annulus_support_2026-10-18_01" / "annulus_pair_profiles.tsv"
TARGETS = DATA / "download_01" / "download_manifest.json"
CONTROLS = DATA / "control_download_01" / "download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_forced_position_radial_relations_2026-11-06_01"
RADIAL_BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
FILTERS = ("i", "i2")
REPLICATES = 10_000
SEED = 20261111


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def cluster_t(beta: np.ndarray, residual: np.ndarray, design: np.ndarray,
              bread: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(labels)
    n, k = design.shape
    if n <= k or len(unique) < 2:
        return math.nan
    scores = np.vstack([design[labels == group].T @ residual[labels == group]
                        for group in unique])
    correction = (len(unique) / (len(unique) - 1)) * ((n - 1) / (n - k))
    covariance = correction * bread @ (scores.T @ scores) @ bread
    variance = covariance[-1, -1]
    if not np.isfinite(variance) or variance <= 0:
        return math.nan
    return float(beta[-1] / math.sqrt(variance))


def prepare(rows: list[dict], mode: str = "continuous_icl") -> dict:
    y = np.asarray([float(row["satellite_delta_mean"]) for row in rows])
    icl = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    log_richness = np.log(np.asarray([float(row["richness"]) for row in rows]))
    scale_icl = float(icl.std())
    scale_z = float(z.std())
    scale_richness = float(log_richness.std())
    if min(scale_icl, scale_z, scale_richness) <= 0:
        raise ValueError("ICL or adjustment covariate has zero variance")
    signal = (icl - icl.mean()) / scale_icl
    z = (z - z.mean()) / scale_z
    log_richness = (log_richness - log_richness.mean()) / scale_richness
    nuisance = np.column_stack((np.ones(len(rows)), z, log_richness))
    design = np.column_stack((nuisance, signal))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("ICL radial-relation design is rank deficient")
    bread = np.linalg.inv(design.T @ design)
    beta = bread @ design.T @ y
    residual = y - design @ beta
    nuisance_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_fit = nuisance @ (nuisance_bread @ nuisance.T @ y)
    null_residual = y - null_fit
    leverage = np.sum((nuisance @ nuisance_bread) * nuisance, axis=1)
    hc2_residual = null_residual / np.sqrt(np.maximum(1.0 - leverage, 0.05))
    blocks = np.asarray([str(row["global_connected_5deg_block"]) for row in rows])
    observed_t = cluster_t(beta, residual, design, bread, blocks)
    if not math.isfinite(observed_t):
        raise ValueError("invalid observed cluster-robust statistic")
    return {
        "rows": rows, "design": design, "bread": bread, "beta": beta,
        "blocks": blocks, "null_fit": null_fit, "hc2_residual": hc2_residual,
        "observed_t": observed_t, "estimate": float(beta[-1]),
        "se": float(abs(beta[-1] / observed_t)) if observed_t else math.nan,
        "n_clusters": len(rows), "n_blocks": len(np.unique(blocks)),
        "largest_block": int(max(np.count_nonzero(blocks == b)
                                  for b in np.unique(blocks))),
    }


def run(output: Path, replicates: int = REPLICATES) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicates < 1000:
        raise ValueError("at least 1,000 shared-block bootstrap replicates are required")
    summary_path = FORCED / "summary.json"
    profiles_path = FORCED / "forced_satellite_profiles.tsv"
    required = (summary_path, profiles_path, ICL, TARGETS, CONTROLS)
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    forced_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (forced_summary.get("status") != "complete_exploratory"
            or forced_summary.get("science_measurement") is not False
            or forced_summary.get("method", {}).get("target_coverage_floor") != 0.2
            or forced_summary.get("method", {}).get("control_coverage_floor") != 0.2
            or forced_summary.get("method", {}).get("mask_other_unmatched_camira") is not False):
        raise ValueError("forced profiles do not match the pinned 20/20 unmatched-CAMIRA scope")
    for name, digest in forced_summary.get("output_sha256", {}).items():
        path = FORCED / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"forced-profile output hash mismatch: {name}")

    target_manifest = json.loads(TARGETS.read_text(encoding="utf-8"))
    control_manifest = json.loads(CONTROLS.read_text(encoding="utf-8"))
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    control_by_id = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    icl_rows = [row for row in read_tsv(ICL) if row["radial_bin_kpc"] == "200-300"]
    icl_by_key = {(str(row["id"]), row["filter"]): row for row in icl_rows}
    if len(icl_by_key) != 224:
        raise ValueError(f"expected 224 ICL-supported cluster/filter pairs, got {len(icl_by_key)}")

    block_input = []
    for cluster_id, band in sorted(icl_by_key, key=lambda item: (item[1], int(item[0]))):
        if cluster_id not in target_by_id or cluster_id not in control_by_id:
            raise ValueError(f"missing target/control coordinate for {cluster_id}")
        block_input.append({
            "cluster_id": cluster_id, "filter": band,
            "target": target_by_id[cluster_id], "control": control_by_id[cluster_id],
        })
    targets = [{"ra_deg": float(row["target"]["ra_deg"]), "dec_deg": float(row["target"]["dec_deg"])}
               for row in block_input]
    controls = [{"ra_deg": float(row["control"]["ra_deg"]), "dec_deg": float(row["control"]["dec_deg"])}
                for row in block_input]
    indices, labels, *_ = pair_block_labels(targets, controls, np.ones(len(block_input), dtype=bool))
    if len(indices) != len(block_input):
        raise ValueError("spatial block labeling changed the 224-pair row order")
    block_by_key = {(row["cluster_id"], row["filter"]): str(label)
                    for row, label in zip(block_input, labels)}

    profiles = read_tsv(profiles_path)
    grouped_positions: dict[tuple[str, str, str], list[dict]] = {}
    for row in profiles:
        if row.get("science_measurement", "False").lower() != "false":
            raise ValueError("forced profile table contains science_measurement promotion")
        key = (str(row["cluster_id"]), row["filter"], row["radial_bin_kpc"])
        if key[:2] not in icl_by_key:
            continue
        grouped_positions.setdefault(key, []).append(row)

    by_group: dict[tuple[str, str], list[dict]] = {}
    for (cluster_id, band, radial), member_rows in grouped_positions.items():
        icle = icl_by_key[(cluster_id, band)]
        target_mean = float(np.mean([float(row["delta_sb_njy_arcsec2"]) for row in member_rows]))
        by_group.setdefault((band, radial), []).append({
            "cluster_id": cluster_id, "filter": band, "radial_bin_kpc": radial,
            "redshift": float(icle["z"]), "richness": float(icle["richness"]),
            "icl_delta_200_300": float(icle["delta_sb_njy_arcsec2"]),
            "satellite_delta_mean": target_mean,
            "n_forced_positions": len(member_rows),
            "global_connected_5deg_block": block_by_key[(cluster_id, band)],
        })

    prepared: dict[tuple[str, str], dict] = {}
    support_rows = []
    for band in FILTERS:
        for radial in RADIAL_BINS:
            rows = sorted(by_group.get((band, radial), []),
                          key=lambda row: int(row["cluster_id"]))
            if not rows:
                support_rows.append({"filter": band, "radial_bin_kpc": radial,
                                     "n_clusters": 0, "n_positions": 0,
                                     "n_blocks": 0, "status": "no_support"})
                continue
            n_blocks = len({row["global_connected_5deg_block"] for row in rows})
            support_rows.append({"filter": band, "radial_bin_kpc": radial,
                                 "n_clusters": len(rows),
                                 "n_positions": sum(int(row["n_forced_positions"]) for row in rows),
                                 "n_blocks": n_blocks,
                                 "status": "supported" if n_blocks >= 6 and len(rows) >= 12 else "support_limited"})
            if n_blocks < 6 or len(rows) < 12:
                continue
            try:
                group = prepare(rows)
            except ValueError as exc:
                support_rows[-1]["status"] = f"model_limited: {exc}"
                continue
            prepared[(band, radial)] = group
    if not prepared:
        raise ValueError("no radial bin meets the minimum 12 pairs / 6 spatial blocks")

    all_blocks = sorted(set(block_by_key.values()))
    rng = np.random.default_rng(SEED)
    bootstrap = {key: np.empty(replicates, dtype=float) for key in prepared}
    max_t = np.empty(replicates, dtype=float)
    for draw_index in range(replicates):
        weights = dict(zip(all_blocks, rng.choice(WEBB_WEIGHTS, len(all_blocks), replace=True)))
        maximum = 0.0
        for key, group in prepared.items():
            w = np.asarray([weights[label] for label in group["blocks"]])
            y_star = group["null_fit"] + group["hc2_residual"] * w
            beta_star = group["bread"] @ group["design"].T @ y_star
            residual_star = y_star - group["design"] @ beta_star
            statistic = cluster_t(beta_star, residual_star, group["design"],
                                  group["bread"], group["blocks"])
            if not math.isfinite(statistic):
                raise ValueError(f"invalid wild-bootstrap statistic for {key}")
            bootstrap[key][draw_index] = statistic
            maximum = max(maximum, abs(statistic))
        max_t[draw_index] = maximum

    raw_p = {key: float((1 + np.count_nonzero(np.abs(draw) >= abs(prepared[key]["observed_t"])))
                        / (replicates + 1)) for key, draw in bootstrap.items()}
    holm_values = holm_adjust(np.asarray(list(raw_p.values())))
    holm = {key: float(value) for key, value in zip(prepared, holm_values)}
    records = []
    for key, group in prepared.items():
        band, radial = key
        critical = float(student_t.ppf(0.975, group["n_blocks"] - 1))
        low = group["estimate"] - critical * group["se"]
        high = group["estimate"] + critical * group["se"]
        records.append({
            "filter": band, "radial_bin_kpc": radial,
            "n_clusters": group["n_clusters"], "n_forced_positions": sum(int(r["n_forced_positions"]) for r in group["rows"]),
            "n_connected_5deg_blocks": group["n_blocks"], "largest_block_n": group["largest_block"],
            "icl_slope_per_sd": group["estimate"], "cr1_se": group["se"],
            "cr1_ci95_low": low, "cr1_ci95_high": high,
            "t_observed": group["observed_t"], "p_webb_raw": raw_p[key],
            "p_holm_across_supported": holm[key],
            "p_maxT_across_supported_bins": float((1 + np.count_nonzero(max_t >= abs(group["observed_t"]))) / (replicates + 1)),
        })

    output.mkdir(parents=True, exist_ok=False)
    support_path = output / "radial_support.tsv"
    write_tsv(support_path, support_rows,
              ["filter", "radial_bin_kpc", "n_clusters", "n_positions", "n_blocks", "status"])
    result_path = output / "icl_forced_position_radial_relations.tsv"
    write_tsv(result_path, records, list(records[0]))
    max_t_path = output / "maxT_null.tsv"
    write_tsv(max_t_path, [{"replica": i + 1, "max_abs_t": float(value)}
                           for i, value in enumerate(max_t)], ["replica", "max_abs_t"])
    result = {
        "schema": "astrafek-icl-forced-position-radial-relations-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "hsc_network_used": False, "image_pixels_read": False,
        "n_icl_supported_pairs": len(icl_by_key),
        "n_relation_tests": len(records), "replicates": replicates, "seed": SEED,
        "method": {
            "response": "mean delta surface brightness at unmatched CAMIRA positions within cluster/filter/radial bin",
            "predictor": "continuous 200-300 proper-kpc ICL target-minus-control residual, standardized within filter/radial bin",
            "adjustment": "standardized redshift and log CAMIRA richness",
            "profile_stack": "equal weight per unmatched CAMIRA position within cluster/filter/radial bin; OLS relation gives equal weight per supported cluster/filter row",
            "covariance": "CR1 cluster-robust covariance using connected target/control 5-degree blocks",
            "bootstrap": "restricted HC2 Webb six-point wild-cluster bootstrap-t; shared weights across all eligible filter/radius models and all 224 ICL-supported spatial blocks",
            "familywise": "Holm and single-step max-T across all filter-by-radius continuous slopes meeting the declared support threshold (at most 14)",
            "support_rule": "input forced profiles already require at least three rotated controls and >=20 percent target and control coverage in each annulus",
        },
        "radial_support": support_rows,
        "results": records,
        "input_sha256": {str(path): sha256(path) for path in required},
        "output_sha256": {support_path.name: sha256(support_path),
                          result_path.name: sha256(result_path), max_t_path.name: sha256(max_t_path)},
        "runner_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "forced CAMIRA positions lack unique SEP detections and are not confirmed isolated satellites; the response is unmatched-position light",
            "20/20 coverage is a post-hoc sensitivity; missing-position and cluster support are non-random and no transport correction is applied",
            "inner radial bins have sparse support and are reported in the support ledger; only bins with at least 12 cluster rows and 6 connected blocks enter inference",
            "familywise adjustment covers the supported filter/radius tests in this analysis only, not the broader exploratory sequence or threshold choices",
            "measurement uncertainty, source Poisson, PSF wings, sky transfer, total coadd covariance, and injection response are not propagated",
            "no correction or science promotion; science_measurement remains false",
        ],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    result = run(args.output, args.replicates)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_relation_tests": result["n_relation_tests"],
                      "radial_support": result["radial_support"],
                      "results": result["results"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
