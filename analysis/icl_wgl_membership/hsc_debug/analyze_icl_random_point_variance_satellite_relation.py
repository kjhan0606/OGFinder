#!/usr/bin/env python3
"""Test whether annular HSC random-point conditions change the ICL/satellite relation.

Joins audited random-point medians to the fixed 94-pair ICL/satellite sample.
It fits one local-condition covariate at a time beyond redshift and richness,
using shared connected-sky block resamples and one familywise max-T correction
across all metrics and both filters. This is a selected-sample sensitivity,
not a correction or causal analysis.
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
SUPPORT_SUMMARY = SUPPORT_DIR / "summary.json"
AUDIT_DIR = ROOT / "icl_random_pixel_variance_audit_2026-11-06_01"
AUDIT_SUMMARY = AUDIT_DIR / "summary.json"
CENTER_TABLE = AUDIT_DIR / "center_random_point_summary.tsv"
PLAN_SHA256 = "eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52"
SUPPORT_SHA256 = "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"
AUDIT_SHA256 = "96245d565375635a994c93ec3294fc581d007161c4479232d5091bbd3cf52389"
CENTER_TABLE_SHA256 = "3137ea8653ae96bd0bbe0e4cd4627183b2ea81f799bc785c3854dd4ab932a0b4"
DEFAULT_OUTPUT = ROOT / "icl_random_point_variance_satellite_relation_2026-11-06_01"
FILTERS = ("i", "i2")
REPLICATES = 5000
SEED = 20261106

# Target/control contrasts and pair-average conditions are tested separately.
# All eight features are one-at-a-time additions to the same z/richness model.
FEATURES = (
    "delta_log_pix_variance", "mean_log_pix_variance",
    "delta_log_sky_std", "mean_log_sky_std",
    "delta_log_inputcount", "mean_log_inputcount",
    "delta_log_psf_trace", "mean_log_psf_trace",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> str:
    if not rows:
        raise ValueError(f"refusing empty output table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    return sha256(path)


def load_join() -> tuple[list[dict], dict[str, str]]:
    for path in (SUPPORT_PATH, SUPPORT_SUMMARY, AUDIT_SUMMARY, CENTER_TABLE):
        if not path.is_file():
            raise FileNotFoundError(path)
    observed = {
        "supported_pairs": sha256(SUPPORT_PATH),
        "random_point_audit": sha256(AUDIT_SUMMARY),
        "random_point_centers": sha256(CENTER_TABLE),
    }
    if (observed["supported_pairs"] != SUPPORT_SHA256
            or observed["random_point_audit"] != AUDIT_SHA256
            or observed["random_point_centers"] != CENTER_TABLE_SHA256):
        raise ValueError("pinned support or random-point audit input hash changed")

    support_doc = json.loads(SUPPORT_SUMMARY.read_text(encoding="utf-8"))
    audit_doc = json.loads(AUDIT_SUMMARY.read_text(encoding="utf-8"))
    if (support_doc.get("status") != "complete_diagnostic"
            or support_doc.get("science_measurement") is not False
            or audit_doc.get("schema") != "astrafek-icl-random-point-variance-audit-v1"
            or audit_doc.get("status") != "complete_exploratory_audit"
            or audit_doc.get("science_measurement") is not False
            or audit_doc.get("plan_sha256") != PLAN_SHA256
            or audit_doc.get("n_returned_random_points") != 40048
            or audit_doc.get("n_planned_centers") != 1598
            or audit_doc.get("n_centers_with_points") != 1587
            or audit_doc.get("n_matched_target_control_pairs") != 798
            or audit_doc.get("matched_pairs_with_both_endpoint_statistics", {}).get("i_pix_variance") != 787
            or audit_doc.get("not_returned", {}).get("images") is not True
            or audit_doc.get("not_returned", {}).get("shape_rows") is not True
            or audit_doc.get("not_returned", {}).get("photoz_rows") is not True):
        raise ValueError("source audit is incomplete or outside the approved scope")
    if (audit_doc.get("outputs", {}).get(CENTER_TABLE.name) != observed["random_point_centers"]):
        raise ValueError("random-point center table hash differs from the audit manifest")

    center_rows = read_tsv(CENTER_TABLE)
    centers = {row["center_id"]: row for row in center_rows}
    if len(center_rows) != 1598 or len(centers) != 1598:
        raise ValueError("random-point audit must contain exactly 1,598 centers")
    joined = []
    support_rows = read_tsv(SUPPORT_PATH)
    for source in support_rows:
        if source.get("science_measurement") != "False":
            raise ValueError("supported-pair parent unexpectedly promotes science")
        cluster_id, control_id = source["cluster_id"], source["control_id"]
        target = centers.get(f"target:{cluster_id}")
        control = centers.get(f"control:{control_id}")
        if target is None or control is None:
            raise ValueError(f"missing random-point endpoint for {cluster_id}/{control_id}")
        if (target["sample_role"] != "target" or control["sample_role"] != "control"
                or target["target_cluster_id"] != cluster_id
                or control["target_cluster_id"] != cluster_id):
            raise ValueError(f"random-point target/control identity mismatch for {cluster_id}")

        out = {key: source[key] for key in (
            "cluster_id", "control_id", "filter", "icl_delta_200_300",
            "satellite_delta_50_75", "redshift", "richness", "target_ra_deg",
            "target_dec_deg", "control_ra_deg", "control_dec_deg")}
        if (int(target["n_random_points"]) <= 0 or int(control["n_random_points"]) <= 0):
            raise ValueError(f"supported pair lacks random-point support: {cluster_id}")
        out["target_random_point_count"] = int(target["n_random_points"])
        out["control_random_point_count"] = int(control["n_random_points"])
        metrics = {
            "pix_variance": ("median_i_pix_variance", "median_i_pix_variance"),
            "sky_std": ("median_i_sky_std", "median_i_sky_std"),
            "inputcount": ("median_i_inputcount_value", "median_i_inputcount_value"),
            "psf_trace": ("median_i_sdssshape_psf_shape11", "median_i_sdssshape_psf_shape22"),
        }
        for name, (field_a, field_b) in metrics.items():
            ta, tb = float(target[field_a]), float(target[field_b])
            ca, cb = float(control[field_a]), float(control[field_b])
            tv = ta + tb if name == "psf_trace" else ta
            cv = ca + cb if name == "psf_trace" else ca
            if not (math.isfinite(tv) and math.isfinite(cv) and tv > 0 and cv > 0):
                raise ValueError(f"invalid random-point {name} for {cluster_id}")
            out[f"target_{name}"] = tv
            out[f"control_{name}"] = cv
            out[f"delta_log_{name}"] = math.log(tv / cv)
            out[f"mean_log_{name}"] = 0.5 * (math.log(tv) + math.log(cv))
        for flag in ("interpolatedcenter", "saturatedcenter", "crcenter", "suspectcenter"):
            a = float(target[f"fraction_true_i_pixelflags_{flag}"])
            b = float(control[f"fraction_true_i_pixelflags_{flag}"])
            out[f"target_mask_{flag}_fraction"] = a
            out[f"control_mask_{flag}_fraction"] = b
        joined.append(out)

    counts = {band: sum(row["filter"] == band for row in joined) for band in FILTERS}
    if len(joined) != 94 or counts != {"i": 26, "i2": 68}:
        raise ValueError(f"supported sample changed: {len(joined)} pairs, {counts}")
    return joined, {
        **{str(path): sha256(path) for path in (SUPPORT_PATH, SUPPORT_SUMMARY,
                                                  AUDIT_SUMMARY, CENTER_TABLE)},
        "approved_plan_sha256": PLAN_SHA256,
        "batch_hashes_sha256": audit_doc["input_sha256"]["batch_hashes_sha256"],
    }


def partial_r(rows: list[dict], feature: str | None,
              indices: np.ndarray | None = None) -> float:
    chosen = rows if indices is None else [rows[int(i)] for i in indices]
    x = np.asarray([float(row["icl_delta_200_300"]) for row in chosen])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in chosen])
    columns = [np.ones(len(chosen)),
               np.asarray([float(row["redshift"]) for row in chosen]),
               np.asarray([float(row["richness"]) for row in chosen])]
    if feature is not None:
        columns.append(np.asarray([float(row[feature]) for row in chosen]))
    design = np.column_stack(columns)
    if len(chosen) <= design.shape[1] + 2 or np.linalg.matrix_rank(design) != design.shape[1]:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def block_labels(rows: list[dict]) -> np.ndarray:
    target = [{"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])}
              for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])}
               for row in rows]
    _, labels, *_ = pair_block_labels(target, control, np.ones(len(rows), dtype=bool))
    return np.asarray(labels, dtype=str)


def run(output: Path, replicas: int = REPLICATES, seed: int = SEED) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 spatial block-bootstrap replicates are required")
    joined, input_hashes = load_join()
    by_filter = {band: [row for row in joined if row["filter"] == band] for band in FILTERS}
    labels = {band: block_labels(rows) for band, rows in by_filter.items()}
    block_members = {band: {block: np.flatnonzero(labels[band] == block)
                            for block in np.unique(labels[band])} for band in FILTERS}
    for band in FILTERS:
        if len(block_members[band]) < 6:
            raise ValueError(f"insufficient connected sky blocks in {band}: {len(block_members[band])}")

    baseline = {band: partial_r(rows, None) for band, rows in by_filter.items()}
    observed = {}
    draws = {(band, feature): [] for band in FILTERS for feature in FEATURES}
    for band in FILTERS:
        for feature in FEATURES:
            value = partial_r(by_filter[band], feature)
            if not math.isfinite(value):
                raise ValueError(f"non-finite adjusted partial correlation {band}/{feature}")
            observed[(band, feature)] = value - baseline[band]

    rng = np.random.default_rng(seed)
    valid_replicates = 0
    attempts = 0
    while valid_replicates < replicas:
        attempts += 1
        if attempts > replicas * 3:
            raise RuntimeError("could not obtain enough complete joint block-bootstrap draws")
        estimates = {}
        for band in FILTERS:
            members = block_members[band]
            blocks = np.asarray(list(members), dtype=str)
            chosen = rng.choice(blocks, size=len(blocks), replace=True)
            indices = np.concatenate([members[block] for block in chosen])
            rows = by_filter[band]
            base_value = partial_r(rows, None, indices)
            if not math.isfinite(base_value):
                break
            for feature in FEATURES:
                adjusted = partial_r(rows, feature, indices)
                if not math.isfinite(adjusted):
                    break
                estimates[(band, feature)] = adjusted - base_value
            if len(estimates) != (FILTERS.index(band) + 1) * len(FEATURES):
                break
        if len(estimates) != len(FILTERS) * len(FEATURES):
            continue
        for key, value in estimates.items():
            draws[key].append(value)
        valid_replicates += 1

    ses = {key: float(np.std(values, ddof=1)) for key, values in draws.items()}
    if any(not math.isfinite(se) or se <= 0 for se in ses.values()):
        raise ValueError("degenerate paired block-bootstrap scale")
    max_abs = np.empty(replicas, dtype=float)
    for index in range(replicas):
        max_abs[index] = max(abs((draws[key][index] - observed[key]) / ses[key])
                             for key in observed)
    critical = float(np.quantile(max_abs, 0.95))
    results = {band: {} for band in FILTERS}
    for band in FILTERS:
        for feature in FEATURES:
            key = (band, feature)
            delta = observed[key]
            se = ses[key]
            draws_array = np.asarray(draws[key])
            threshold = abs(delta / se)
            results[band][feature] = {
                "n_pairs": len(by_filter[band]),
                "baseline_partial_r": baseline[band],
                "adjusted_partial_r": baseline[band] + delta,
                "change_from_baseline": delta,
                "paired_block_bootstrap_se": se,
                "paired_change_percentile_ci95": np.percentile(draws_array, [2.5, 97.5]).tolist(),
                "paired_change_maxT_ci95": [delta - critical * se, delta + critical * se],
                "maxT_adjusted_p": float((1 + np.count_nonzero(max_abs >= threshold)) / (replicas + 1)),
            }

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "joined_supported_pairs.tsv"
    joined_hash = write_tsv(joined_path, joined)
    summary = {
        "schema": "astrafek-icl-random-point-variance-satellite-relation-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "plan_sha256": PLAN_SHA256,
        "n_pairs": len(joined),
        "filter_counts": {band: len(by_filter[band]) for band in FILTERS},
        "outcome": "partial correlation between ICL 200-300 kpc contrast and satellite 50-75 kpc excess",
        "baseline_controls": ["redshift", "linear CAMIRA richness"],
        "random_point_features": list(FEATURES),
        "paired_block_resampling": {
            "method": "connected 5-degree target/control endpoint blocks, resampled within filter",
            "replicas_requested": replicas,
            "replicas_valid_joint": valid_replicates,
            "seed": seed,
            "family": "all eight random-point condition adjustments across both filters (16 tests)",
            "maxT_critical_95": critical,
            "blocks_by_filter": {band: len(block_members[band]) for band in FILTERS},
        },
        "results": results,
        "inputs_sha256": input_hashes,
        "joined_pairs_path": str(joined_path),
        "joined_pairs_sha256": joined_hash,
        "limitations": [
            "Only the selected 94-pair support intersection is analyzed (26 i, 68 i2); the sample does not represent all 800 targets.",
            "Random-point medians are sparse local summaries, not pixel variance maps or exact ICL-estimator covariance.",
            "One-at-a-time covariate adjustment is a sensitivity analysis, not confounding control or a causal/mediation model.",
            "Random-point median uncertainty is not propagated; centers have varying random-point counts.",
            "Connected 5-degree cells are approximate spatial blocks rather than validated HSC field regions.",
            "No correction is applied to the saved ICL or satellite stack measurements.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=REPLICATES)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    summary = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "filter_counts": summary["filter_counts"],
                      "minimum_maxT_p": min(result[feature]["maxT_adjusted_p"]
                                             for result in summary["results"].values()
                                             for feature in FEATURES),
                      "output": str(args.output.resolve()),
                      "summary_sha256": summary["summary_sha256"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
