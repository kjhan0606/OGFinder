#!/usr/bin/env python3
"""Check whether coarse HSC frame metadata changes the ICL--satellite relation.

Offline sensitivity analysis of the fixed 94-pair supported sample. It adds
one target/control frame-metadata contrast at a time to the existing
redshift-and-richness adjusted ICL--satellite partial correlation. It does not
convert frame summaries into pixel weights or noise corrections.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_DIR = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
SUPPORT_PATH = SUPPORT_DIR / "supported_pairs.tsv"
SUPPORT_SHA256 = "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"
METADATA_DIR = ROOT / "icl_frame_metadata_row_audit_2026-11-06_04"
METADATA_SUMMARY = METADATA_DIR / "summary.json"
METADATA_PATH = METADATA_DIR / "center_frame_metadata_aggregates.tsv"
METADATA_SHA256 = "2560e8ece327bdd5ae02108f5187b40efa6bf3fa8b602b929287da9b64dffaf5"
PLAN_SHA256 = "49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3"
DEFAULT_OUTPUT = ROOT / "icl_frame_metadata_satellite_relation_2026-11-06_03"
FILTERS = ("i", "i2")
METRICS = (
    "delta_log_frame_count",
    "delta_log_visit_count",
    "delta_log_sky_rms_proxy",
    "delta_log_fluxmag0",
    "delta_log_gain_geomean",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_center_metadata(path: Path) -> dict[str, dict[str, float | str]]:
    output = {}
    for row in read_tsv(path):
        center_id = row["center_id"]
        if center_id in output:
            raise ValueError(f"duplicate center metadata row: {center_id}")
        vals = {key: float(row[key]) for key in (
            "n_unique_frames", "n_unique_visits", "median_sigma_sky",
            "median_fluxmag0", "median_gain1", "median_gain2",
            "median_gain3", "median_gain4")}
        if any(not math.isfinite(value) or value <= 0 for value in vals.values()):
            raise ValueError(f"non-positive or non-finite frame metadata for {center_id}")
        vals["sample_role"] = row["sample_role"]
        vals["sky_rms_proxy"] = vals["median_sigma_sky"] * 3.631e12 / vals["median_fluxmag0"]
        vals["gain_geomean"] = math.exp(float(np.mean(np.log([
            vals[f"median_gain{i}"] for i in range(1, 5)]))))
        output[center_id] = vals
    if len(output) != 1598:
        raise ValueError(f"expected metadata for 1,598 centers; got {len(output)}")
    return output


def join_rows() -> tuple[list[dict], dict[str, str]]:
    for path in (SUPPORT_PATH, METADATA_SUMMARY, METADATA_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(SUPPORT_PATH) != SUPPORT_SHA256 or sha256(METADATA_PATH) != METADATA_SHA256:
        raise ValueError("pinned supported-pair or metadata aggregate input changed")
    support_doc = json.loads((SUPPORT_DIR / "summary.json").read_text(encoding="utf-8"))
    metadata_doc = json.loads(METADATA_SUMMARY.read_text(encoding="utf-8"))
    if (support_doc.get("status") != "complete_diagnostic"
            or support_doc.get("science_measurement") is not False
            or metadata_doc.get("status") != "complete_exploratory"
            or metadata_doc.get("science_measurement") is not False
            or metadata_doc.get("approved_query_submitted") is not True
            or metadata_doc.get("plan_sha256") != PLAN_SHA256
            or metadata_doc.get("n_frame_rows") != 67071
            or metadata_doc.get("n_centers") != 1598):
        raise ValueError("source audit is incomplete or does not match the approved plan")
    if metadata_doc.get("result_files_sha256", {}).get("center_aggregates.tsv") != METADATA_SHA256:
        raise ValueError("metadata audit does not pin the requested aggregate table")

    metadata = read_center_metadata(METADATA_PATH)
    joined = []
    for source in read_tsv(SUPPORT_PATH):
        if source.get("science_measurement") != "False":
            raise ValueError("supported-pair input unexpectedly promotes science_measurement")
        cluster_id, control_id = str(source["cluster_id"]), str(source["control_id"])
        target_key, control_key = f"target:{cluster_id}", f"control:{control_id}"
        target, control = metadata.get(target_key), metadata.get(control_key)
        if target is None or control is None:
            raise ValueError(f"missing metadata for supported pair {cluster_id}/{control_id}")
        if target["sample_role"] != "target" or control["sample_role"] != "control":
            raise ValueError(f"metadata endpoint role mismatch for {cluster_id}")
        row = {key: source[key] for key in (
            "cluster_id", "control_id", "filter", "icl_delta_200_300",
            "satellite_delta_50_75", "redshift", "richness",
            "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg")}
        for field, target_field in (
            ("frame_count", "n_unique_frames"), ("visit_count", "n_unique_visits"),
            ("sky_rms_proxy", "sky_rms_proxy"), ("fluxmag0", "median_fluxmag0"),
            ("gain_geomean", "gain_geomean")):
            a, b = float(target[target_field]), float(control[target_field])
            row[f"target_{field}"] = a
            row[f"control_{field}"] = b
            row[f"delta_log_{field}"] = math.log(a / b)
            row[f"mean_log_{field}"] = 0.5 * (math.log(a) + math.log(b))
        # Keep the metadata table's descriptive label stable in outputs.
        row["delta_log_sky_rms_proxy"] = row.pop("delta_log_sky_rms_proxy")
        row["delta_log_frame_count"] = row.pop("delta_log_frame_count")
        row["delta_log_visit_count"] = row.pop("delta_log_visit_count")
        row["delta_log_fluxmag0"] = row.pop("delta_log_fluxmag0")
        row["delta_log_gain_geomean"] = row.pop("delta_log_gain_geomean")
        joined.append(row)
    if len(joined) != 94 or {(r["filter"], sum(x["filter"] == r["filter"] for x in joined))
                              for r in joined} != {("i", 26), ("i2", 68)}:
        raise ValueError("supported-pair join differs from the fixed 26 i / 68 i2 sample")
    return joined, {
        "supported_pairs_sha256": sha256(SUPPORT_PATH),
        "supported_summary_sha256": sha256(SUPPORT_DIR / "summary.json"),
        "metadata_summary_sha256": sha256(METADATA_SUMMARY),
        "metadata_aggregate_sha256": sha256(METADATA_PATH),
    }


def partial_r(rows: list[dict], metric: str | None, indices: np.ndarray | None = None) -> float:
    chosen = rows if indices is None else [rows[int(i)] for i in indices]
    x = np.asarray([float(row["icl_delta_200_300"]) for row in chosen])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in chosen])
    covariates = [np.asarray([float(row["redshift"]) for row in chosen]),
                  np.asarray([float(row["richness"]) for row in chosen])]
    if metric is not None:
        covariates.append(np.asarray([float(row[metric]) for row in chosen]))
    design = np.column_stack([np.ones(len(chosen)), *covariates])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError(f"rank-deficient controls for {metric or 'baseline'}")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    rows, input_hashes = join_rows()
    estimates: dict[tuple[str, str], float] = {}
    draws: dict[tuple[str, str], list[float]] = {}
    block_info = {}
    rng_by_filter = {}
    block_members_by_filter = {}
    rows_by_filter = {}
    for band_index, band in enumerate(FILTERS):
        selected = [row for row in rows if row["filter"] == band]
        target_coords = [{"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])} for row in selected]
        control_coords = [{"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])} for row in selected]
        _, labels, nt, nc, nblocks, largest = pair_block_labels(
            target_coords, control_coords, np.ones(len(selected), dtype=bool))
        labels = np.asarray(labels, dtype=str)
        groups = np.unique(labels)
        if len(groups) < 6:
            raise ValueError(f"at least six connected spatial blocks required; {band} has {len(groups)}")
        members = {group: np.flatnonzero(labels == group) for group in groups}
        block_info[band] = {"n_pairs": len(selected), "n_blocks": len(groups),
                            "n_target_cells": int(nt), "n_control_cells": int(nc),
                            "largest_block_pairs": int(largest)}
        rows_by_filter[band] = selected
        block_members_by_filter[band] = members
        rng_by_filter[band] = np.random.default_rng(seed + band_index * 1000)
        for metric in ("baseline", *METRICS):
            key_metric = None if metric == "baseline" else metric
            estimates[(band, metric)] = partial_r(selected, key_metric)
            draws[(band, metric)] = []

    for _ in range(replicas):
        for band in FILTERS:
            members = block_members_by_filter[band]
            groups = np.asarray(list(members), dtype=str)
            selected_groups = rng_by_filter[band].choice(groups, size=len(groups), replace=True)
            indices = np.concatenate([members[group] for group in selected_groups])
            selected = rows_by_filter[band]
            values = {metric: partial_r(selected, None if metric == "baseline" else metric, indices)
                      for metric in ("baseline", *METRICS)}
            if all(np.isfinite(value) for value in values.values()):
                for metric, value in values.items():
                    draws[(band, metric)].append(float(value))

    metric_results = {}
    for band in FILTERS:
        for metric in METRICS:
            base_draws = np.asarray(draws[(band, "baseline")])
            adj_draws = np.asarray(draws[(band, metric)])
            valid = min(len(base_draws), len(adj_draws))
            if valid < replicas * 0.9:
                raise ValueError(f"too few paired valid block draws for {band}/{metric}: {valid}/{replicas}")
            delta_draws = adj_draws[:valid] - base_draws[:valid]
            delta = estimates[(band, metric)] - estimates[(band, "baseline")]
            metric_results[(band, metric)] = {
                "baseline_partial_r": estimates[(band, "baseline")],
                "metadata_adjusted_partial_r": estimates[(band, metric)],
                "paired_change_in_r": delta,
                "paired_change_bootstrap_se": float(np.std(delta_draws, ddof=1)),
                "paired_change_ci95_pointwise": np.percentile(delta_draws, [2.5, 97.5]).tolist(),
                "replicas_valid": valid,
            }

    # One max-T family covers five metadata adjustments in each filter.
    family = [(band, metric) for band in FILTERS for metric in METRICS]
    standardized_draws = []
    for replicate in range(replicas):
        values = []
        ok = True
        for band, metric in family:
            b = draws[(band, "baseline")]
            a = draws[(band, metric)]
            if replicate >= min(len(b), len(a)):
                ok = False
                break
            se = metric_results[(band, metric)]["paired_change_bootstrap_se"]
            if not np.isfinite(se) or se <= 0:
                ok = False
                break
            observed_delta = metric_results[(band, metric)]["paired_change_in_r"]
            values.append(abs(((a[replicate] - b[replicate]) - observed_delta) / se))
        if ok:
            standardized_draws.append(max(values))
    if len(standardized_draws) < replicas * 0.9:
        raise ValueError(f"too few complete max-T draws: {len(standardized_draws)}/{replicas}")
    critical = float(np.quantile(standardized_draws, 0.95))
    max_t = np.asarray(standardized_draws)
    for band, metric in family:
        result = metric_results[(band, metric)]
        se = result["paired_change_bootstrap_se"]
        delta = result["paired_change_in_r"]
        result["paired_change_ci95_maxT"] = [delta - critical * se, delta + critical * se]
        result["maxT_adjusted_p"] = float((1 + np.count_nonzero(max_t >= abs(delta / se))) / (len(max_t) + 1))

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "joined_supported_pairs.tsv"
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema": "astrafek-icl-frame-metadata-satellite-relation-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "plan_sha256": PLAN_SHA256,
        "n_pairs": len(rows),
        "filter_counts": {band: sum(row["filter"] == band for row in rows) for band in FILTERS},
        "outcome": "existing ICL 200--300 kpc target-control contrast correlated with satellite 50--75 kpc excess",
        "baseline_controls": ["redshift", "linear CAMIRA richness"],
        "metadata_adjustments": list(METRICS),
        "interpretation": "Each coarse target/control frame-summary offset is added one at a time; adjustment is a sensitivity diagnostic, not a causal or noise correction.",
        "spatial_resampling": {
            "method": "connected 5-degree target/control endpoint blocks, resampled within filter",
            "replicas_requested": replicas,
            "replicas_seed_base": seed,
            "family": "five metadata adjustment changes jointly across i and i2 (10 tests)",
            "maxT_critical_95": critical,
            "replicas_valid_joint": len(standardized_draws),
            "by_filter": block_info,
        },
        "results": {band: {metric: metric_results[(band, metric)] for metric in METRICS}
                    for band in FILTERS},
        "inputs_sha256": input_hashes,
        "runner_sha256": sha256(Path(__file__).resolve()),
        "joined_pairs_path": str(joined_path),
        "joined_pairs_sha256": sha256(joined_path),
        "limitations": [
            "Only 26 i and 68 i2 profile pairs are supported; estimates are low-power and selected-sample diagnostics.",
            "Frame metadata are pooled patch-level summaries and do not resolve pixel-specific input membership or effective weights.",
            "Sky-RMS conversion is a proxy; it is not a coadd-noise or source-Poisson estimate.",
            "Five metadata adjustments are tested in one max-T family; the spatial blocks remain approximate.",
            "No correlation adjustment is applied to the published profiles and science_measurement remains false.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicates are required")
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "results": result["results"], "output": str(args.output.resolve()),
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
