#!/usr/bin/env python3
"""Audit support selection and familywise uncertainty for full-sample ICL covariance.

Reads only the saved full-783 annular covariance tables and target manifest.
It does not read images or contact HSC.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUT = DATA / "icl_annular_coadd_covariance_full783_2026-11-06_02"
OLD_INPUT = DATA / "icl_annular_coadd_covariance_2026-11-06_01"
MANIFEST = DATA / "download_01/download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_annular_covariance_full783_support_audit_2026-11-06_01"
METRICS = ("zero_lag", "sum8", "coverage")
REPLICATES = 10000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def num(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    result = float(value)
    return result if np.isfinite(result) else None


def familywise_block_bootstrap(rows: list[dict[str, str]], band: str,
                               seed: int) -> dict:
    metric_values: dict[str, dict[str, list[float]]] = {}
    observed = {}
    union_blocks: set[str] = set()
    for metric in METRICS:
        column = f"target_minus_control_{metric}"
        by_block: dict[str, list[float]] = {}
        for row in rows:
            value = num(row.get(column))
            if row["filter"] == band and value is not None:
                by_block.setdefault(row["connected_5deg_block"], []).append(value)
        metric_values[metric] = by_block
        union_blocks.update(by_block)
        values = [value for block_values in by_block.values() for value in block_values]
        if values:
            observed[metric] = float(np.median(values))
    blocks = sorted(union_blocks)
    if len(observed) != len(METRICS) or len(blocks) < 3:
        raise ValueError(f"insufficient complete block support for {band}: {len(blocks)}")

    rng = np.random.default_rng(seed)
    boot = np.empty((REPLICATES, len(METRICS)), dtype=np.float64)
    valid = 0
    attempts = 0
    while valid < REPLICATES:
        attempts += 1
        if attempts > REPLICATES * 2:
            raise RuntimeError(f"could not draw complete joint block samples for {band}")
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        estimates = []
        for metric in METRICS:
            values = [value for block in chosen
                      for value in metric_values[metric].get(str(block), [])]
            if not values:
                break
            estimates.append(float(np.median(values)))
        if len(estimates) == len(METRICS):
            boot[valid] = estimates
            valid += 1

    scales = np.std(boot, axis=0, ddof=1)
    if np.any(scales <= 0) or not np.all(np.isfinite(scales)):
        raise ValueError(f"degenerate bootstrap scale in {band}")
    obs = np.asarray([observed[metric] for metric in METRICS])
    centered = (boot - obs) / scales
    max_abs = np.max(np.abs(centered), axis=1)
    critical = float(np.quantile(max_abs, 0.95))
    result = {"n_blocks": len(blocks), "replicates": REPLICATES,
              "seed": seed, "family_size": len(METRICS),
              "familywise_max_t_critical_95": critical,
              "metrics": {}}
    for index, metric in enumerate(METRICS):
        threshold = abs(obs[index] / scales[index])
        result["metrics"][metric] = {
            "n_pairs": sum(len(values) for values in metric_values[metric].values()),
            "median_target_minus_control": float(obs[index]),
            "block_bootstrap_se": float(scales[index]),
            "percentile_ci95": [float(x) for x in np.quantile(boot[:, index], [0.025, 0.975])],
            "max_t_simultaneous_ci95": [
                float(obs[index] - critical * scales[index]),
                float(obs[index] + critical * scales[index]),
            ],
            "familywise_two_sided_max_t_p": float((1 + np.sum(max_abs >= threshold)) / (REPLICATES + 1)),
        }
    return result


def quantile_groups(values: list[float], selected: list[bool], label: str) -> list[dict]:
    ordered = sorted(range(len(values)), key=lambda i: values[i])
    groups = np.array_split(np.asarray(ordered, dtype=int), 4)
    output = []
    for index, group in enumerate(groups, 1):
        vals = [values[int(i)] for i in group]
        chosen = sum(selected[int(i)] for i in group)
        output.append({
            "covariate": label, "quartile": index,
            "lower": float(min(vals)), "upper": float(max(vals)),
            "n_pairs": len(vals), "n_with_complete_covariance_support": chosen,
            "support_fraction": float(chosen / len(vals)) if vals else None,
        })
    return output


def standardized_mean_difference(values: list[float], selected: list[bool]) -> float | None:
    a = np.asarray([v for v, flag in zip(values, selected) if flag], dtype=float)
    b = np.asarray([v for v, flag in zip(values, selected) if not flag], dtype=float)
    if len(a) < 2 or len(b) < 2:
        return None
    pooled = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return float((np.mean(a) - np.mean(b)) / pooled) if pooled > 0 else None


def write_tsv(path: Path, rows: list[dict]) -> str:
    if not rows:
        raise ValueError(f"refusing empty output table: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    return sha256(path)


def run(output: Path, input_dir: Path = INPUT) -> dict:
    output, input_dir = output.resolve(), input_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_path = input_dir / "summary.json"
    per_center_path = input_dir / "per_center_icl_annulus_covariance.tsv"
    pair_path = input_dir / "paired_icl_annulus_covariance.tsv"
    input_summary = read_json(summary_path)
    if (input_summary.get("schema") != "astrafek-icl-full783-annular-coadd-covariance-v1"
            or input_summary.get("status") != "complete_exploratory"
            or input_summary.get("science_measurement") is not False
            or input_summary.get("n_same_filter_pairs") != 783
            or input_summary.get("n_image_records") != 1566
            or input_summary.get("failed_inputs")):
        raise ValueError("full-sample annular covariance input is incomplete or changed scope")
    method = input_summary.get("method", {})
    patch_size = method.get("patch_size_pixels")
    if not isinstance(patch_size, int) or patch_size < 1:
        raise ValueError("input covariance summary has no valid patch-size provenance")
    for path in (per_center_path, pair_path):
        if input_summary.get("outputs_sha256", {}).get(str(path)) != sha256(path):
            raise ValueError(f"full-sample covariance output hash mismatch: {path}")
    center_rows, pair_rows = read_tsv(per_center_path), read_tsv(pair_path)
    if len(center_rows) != 1566 or len(pair_rows) != 783:
        raise ValueError("full-sample covariance table row counts changed")
    manifest = read_json(MANIFEST)
    clusters = {str(row["id"]): row for row in manifest["clusters"]}
    if len(clusters) != 800:
        raise ValueError("target manifest does not contain all 800 clusters")

    center_by_key = {(row["cluster_id"], row["filter"], row["role"]): row
                     for row in center_rows}
    selected_map = {}
    diagnostics = []
    quartiles = []
    smd = {}
    for band in ("i", "i2"):
        band_pairs = [row for row in pair_rows if row["filter"] == band]
        flags = [row["complete_target_control_support"] == "True" for row in band_pairs]
        if len(band_pairs) != input_summary["input_sha256"].get("target_control_pair_counts", {}).get(band):
            raise ValueError(f"unexpected pair count for {band}")
        z = [float(clusters[row["cluster_id"]]["z"]) for row in band_pairs]
        richness = [float(clusters[row["cluster_id"]]["richness"]) for row in band_pairs]
        diagnostics.append({
            "filter": band, "n_same_filter_pairs": len(band_pairs),
            "n_complete_pairs": sum(flags),
            "support_fraction": float(np.mean(flags)),
            "n_complete_pair_blocks": len({row["connected_5deg_block"] for row, flag in zip(band_pairs, flags) if flag}),
            "n_parent_blocks": len({row["connected_5deg_block"] for row in band_pairs}),
            "largest_supported_block_pairs": max((sum(
                1 for row, flag in zip(band_pairs, flags)
                if flag and row["connected_5deg_block"] == block)
                for block in {row["connected_5deg_block"] for row, flag in zip(band_pairs, flags) if flag}),
                default=0),
        })
        quartiles.extend(quantile_groups(z, flags, f"{band}:cluster_redshift"))
        quartiles.extend(quantile_groups(richness, flags, f"{band}:CAMIRA_richness"))
        smd[band] = {
            "redshift_complete_minus_incomplete_smd": standardized_mean_difference(z, flags),
            "richness_complete_minus_incomplete_smd": standardized_mean_difference(richness, flags),
        }
        for pair, flag in zip(band_pairs, flags):
            target = center_by_key[(pair["cluster_id"], band, "target")]
            control = center_by_key[(pair["cluster_id"], band, "matched_control")]
            pair["redshift"] = float(clusters[pair["cluster_id"]]["z"])
            pair["richness"] = float(clusters[pair["cluster_id"]]["richness"])
            pair["target_geometric_patches"] = int(target["n_geometric_whole_patches"])
            pair["control_geometric_patches"] = int(control["n_geometric_whole_patches"])
            pair["target_mask_eligible_patches"] = int(target["n_mask_eligible_whole_patches"])
            pair["control_mask_eligible_patches"] = int(control["n_mask_eligible_whole_patches"])
            pair["target_patches_used"] = int(target["n_patches_used"])
            pair["control_patches_used"] = int(control["n_patches_used"])
            pair["complete_target_control_support"] = bool(flag)
            selected_map[(pair["cluster_id"], band)] = flag
    maxt = {band: familywise_block_bootstrap(pair_rows, band, 20261120 + i)
            for i, band in enumerate(("i", "i2"))}
    output.mkdir(parents=True, exist_ok=False)
    support_path = output / "pair_support_audit.tsv"
    quartile_path = output / "support_by_covariate_quartile.tsv"
    support_hash = write_tsv(support_path, pair_rows)
    quartile_hash = write_tsv(quartile_path, quartiles)
    summary = {
        "schema": "astrafek-icl-full783-annular-covariance-support-audit-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "hsc_network_used": False, "image_pixels_read": False,
        "audit_script_sha256": sha256(Path(__file__).resolve()),
        "input": {
            "summary_sha256": sha256(summary_path),
            "per_center_sha256": sha256(per_center_path),
            "paired_table_sha256": sha256(pair_path),
            "target_manifest_sha256": sha256(MANIFEST),
        },
        "method": {
            "patch_size_pixels": patch_size,
            "max_lag_pixels": method.get("max_lag_pixels"),
            "minimum_patches_per_center_annulus": method.get("minimum_patches_per_center_annulus"),
            "minimum_valid_fraction_per_patch": method.get("minimum_valid_fraction_per_patch"),
        },
        "support_by_filter": diagnostics,
        "support_covariate_smd": smd,
        "support_by_covariate_quartile": quartiles,
        "familywise_max_t_by_filter": maxt,
        "interpretation": [
            f"Results apply only to pairs with both target and control annuli having at least {method.get('minimum_patches_per_center_annulus')} fully contained, at least {method.get('minimum_valid_fraction_per_patch'):.0%}-valid {patch_size}x{patch_size} patches.",
            "The support restriction is substantial and descriptive redshift/richness balance does not establish missing-at-random support.",
            "The lag-sum metric can include astrophysical ICL/cluster residual structure and is not a pure detector-noise estimate.",
            "This is delivered-coadd covariance diagnostics; no variance correction or science promotion is made.",
        ],
        "outputs_sha256": {str(support_path): support_hash, str(quartile_path): quartile_hash},
    }
    result_path = output / "summary.json"
    result_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "support_by_filter": diagnostics,
                      "familywise_max_t": maxt,
                      "summary": str(result_path)}, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--input-dir", type=Path, default=INPUT)
    args = parser.parse_args()
    run(args.output_dir, args.input_dir)


if __name__ == "__main__":
    main()
