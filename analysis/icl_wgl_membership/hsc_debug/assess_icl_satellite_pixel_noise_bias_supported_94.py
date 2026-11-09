#!/usr/bin/env python3
"""Apply the existing finite-lag pixel-noise correction to all 94 supported pairs.

Inputs are the saved 94-pair support table and the already-computed target/
control covariance product. This adapter performs no network access.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from analyze_icl_satellite_supported_94_blockwidth import connected_blocks
from assess_icl_satellite_pixel_noise_bias import (
    corrected_partial_correlation,
    run as run_assessment,
)
from verify_satellite_icl_covariance_support_scope import run as verify_scope


HERE = Path(__file__).resolve().parent
DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_SUPPORT = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
DEFAULT_PREVIOUS = DATA / "icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv"
DEFAULT_ICL_SUPPORT = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
DEFAULT_SATELLITE = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01/satellite_profiles.tsv"
DEFAULT_TARGET_COV = DATA / "satellite_icl_lag_covariance_2026-11-03_02"
DEFAULT_TARGET_CONTROL_COV = DATA / "satellite_icl_target_control_covariance_2026-11-03_02"
DEFAULT_OUTPUT = DATA / "icl_satellite_pixel_noise_bias_supported_94_2026-11-06_04"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def corrected_blockwidth_sensitivity(
    supported: list[dict[str, str]], covariance_dir: Path,
    replicas: int, seed: int,
) -> dict:
    covariance = read_tsv(covariance_dir / "target_control_pixel_covariance.tsv")
    covariance_by_key = {(row["cluster_id"], row["filter"]): row for row in covariance}
    if len(covariance_by_key) != len(covariance):
        raise ValueError("duplicate covariance pair keys")
    output: dict[str, dict] = {}
    for width_index, width in enumerate((2.5, 5.0, 10.0)):
        width_result = {}
        for band_index, band in enumerate(("i", "i2")):
            rows = [row for row in supported if row["filter"] == band]
            if len(rows) != (26 if band == "i" else 68):
                raise ValueError(f"unexpected supported-pair count in {band}")
            labels = connected_blocks(rows, width)
            groups = np.unique(labels)
            by_block = {block: np.flatnonzero(labels == block) for block in groups}
            arrays_by_key = {
                (row["cluster_id"], row["filter"]): covariance_by_key[(row["cluster_id"], row["filter"])]
                for row in rows
            }
            x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
            y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
            z = np.asarray([float(row["redshift"]) for row in rows])
            richness = np.asarray([float(row["richness"]) for row in rows])
            z = (z - z.mean()) / z.std(ddof=1)
            richness = (richness - richness.mean()) / richness.std(ddof=1)
            covariance_vectors = [
                np.asarray([float(arrays_by_key[(row["cluster_id"], row["filter"])][field]) for row in rows])
                for field in (
                    "icl_delta_lag_variance",
                    "target_satellite_lag_variance_njy2_arcsec4",
                    "target_icl_satellite_lag_covariance_njy2_arcsec4",
                )
            ]
            estimate = corrected_partial_correlation(
                x, y, z, richness, *covariance_vectors,
            )
            if estimate is None or estimate["raw_partial_r"] is None or estimate["pixel_noise_corrected_partial_r"] is None:
                raise ValueError(f"invalid corrected point estimate for {band} at {width} degrees")
            rng = np.random.default_rng(seed + width_index * 100 + band_index)
            raw_draws: list[float] = []
            corrected_draws: list[float] = []
            for _ in range(replicas):
                selected = rng.choice(groups, size=len(groups), replace=True)
                weights = np.zeros(len(rows), dtype=float)
                for block in selected:
                    weights[by_block[block]] += 1.0
                draw = corrected_partial_correlation(
                    x, y, z, richness, *covariance_vectors, weights=weights,
                )
                if draw is None:
                    continue
                if draw["raw_partial_r"] is not None:
                    raw_draws.append(float(draw["raw_partial_r"]))
                if draw["pixel_noise_corrected_partial_r"] is not None:
                    corrected_draws.append(float(draw["pixel_noise_corrected_partial_r"]))
            if len(raw_draws) < replicas * 0.9 or len(corrected_draws) < replicas * 0.9:
                raise ValueError(f"too few valid {band} bootstrap replicates at {width} degrees")
            width_result[band] = {
                "n_pairs": len(rows),
                "n_connected_blocks": int(len(groups)),
                "largest_block_pairs": int(max(len(index) for index in by_block.values())),
                "point_estimate": estimate,
                "raw_partial_r_bootstrap_ci95": np.percentile(raw_draws, [2.5, 97.5]).tolist(),
                "pixel_noise_corrected_partial_r_bootstrap_ci95": np.percentile(
                    corrected_draws, [2.5, 97.5],
                ).tolist(),
                "bootstrap_raw_valid": len(raw_draws),
                "bootstrap_noise_corrected_valid": len(corrected_draws),
            }
        output[str(width)] = width_result
    return {
        "block_widths_deg": [2.5, 5.0, 10.0],
        "replicas": replicas,
        "seed": seed,
        "block_definition": "connected components joining pairs that share either a target or control equal-area sky cell",
        "results_by_width_and_filter": output,
    }


def run(args: argparse.Namespace) -> dict:
    for path in (args.support, args.previous_74, args.icl_support, args.satellite_profiles):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    supported = read_tsv(args.support)
    previous = read_tsv(args.previous_74)
    if len(supported) != 94 or len({(r["cluster_id"], r["filter"]) for r in supported}) != 94:
        raise ValueError("expected exactly 94 unique supported cluster/filter pairs")
    if len(previous) != 74:
        raise ValueError("expected the prior fixed 74-pair association sample")
    supported_by_key = {(r["cluster_id"], r["filter"]): r for r in supported}
    for row in previous:
        key = (row["cluster_id"], row["filter"])
        current = supported_by_key.get(key)
        if current is None:
            raise ValueError(f"previous pair missing from the 94-pair support: {key}")
        if (abs(float(row["icl_delta_raw_njy_arcsec2"]) - float(current["icl_delta_200_300"])) > 5e-7
                or abs(float(row["satellite_delta_raw_njy_arcsec2"]) - float(current["satellite_delta_50_75"])) > 1e-12):
            raise ValueError(f"raw pair values changed for the nested 74-pair sample: {key}")

    args.output.mkdir(parents=True, exist_ok=False)
    association_path = args.output / "association_rows.tsv"
    blocks_path = args.output / "reblocked_association_rows.tsv"
    association_rows = [{
        "cluster_id": r["cluster_id"], "target_id": r["cluster_id"], "filter": r["filter"],
        "icl_delta_raw_njy_arcsec2": r["icl_delta_200_300"],
        "satellite_delta_raw_njy_arcsec2": r["satellite_delta_50_75"],
        "z": r["redshift"], "richness": r["richness"],
    } for r in supported]
    block_rows = [{
        "target_id": r["cluster_id"], "filter": r["filter"],
        "connected_5deg_block": r["connected_5deg_block"],
    } for r in supported]
    write_tsv(association_path, list(association_rows[0]), association_rows)
    write_tsv(blocks_path, list(block_rows[0]), block_rows)

    scope = verify_scope(SimpleNamespace(
        icl_support=args.icl_support,
        satellite_profiles=args.satellite_profiles,
        association_rows=association_path,
        target_covariance=args.target_covariance,
        target_control_covariance=args.target_control_covariance,
    ))
    if (scope.get("status") != "verified" or scope.get("n_covariance_output_pairs_verified") != 94
            or scope.get("n_association_pairs_covered") != 94
            or scope.get("missing_association_pairs")):
        raise ValueError("independent covariance support-scope verification failed for 94 pairs")
    scope_path = args.output / "covariance_scope_verification.json"
    scope_path.write_text(json.dumps(scope, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    result = run_assessment(SimpleNamespace(
        association_rows=association_path,
        reblocked_association_rows=blocks_path,
        target_control_covariance=args.target_control_covariance,
        scope_verification=scope_path,
        replicas=args.replicas,
        seed=args.seed,
    ))
    result["schema"] = "astrafek-icl-satellite-pixel-noise-bias-supported-94-v1"
    result["support_scope"] = {
        "n_pairs": 94, "by_filter": {"i": 26, "i2": 68},
        "nested_previous_sample_pairs": 74,
        "support_table_sha256": sha256(args.support),
        "previous_74_pair_value_check": "passed",
        "covariance_scope_verification_sha256": sha256(scope_path),
    }
    result["limitations"].append(
        "the supplied covariance is target/control pixel-noise covariance on fixed supported annuli; it does not account for support selection or PSF wings"
    )
    result["limitations"].append(
        "the connected 5-degree blocks are the saved approximate spatial partition, not independently validated survey fields"
    )
    result["blockwidth_sensitivity"] = corrected_blockwidth_sensitivity(
        supported, args.target_control_covariance, args.replicas, args.seed,
    )
    result["input_sha256"].update({
        str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
        str(args.support.resolve()): sha256(args.support),
        str(args.previous_74.resolve()): sha256(args.previous_74),
        str(args.icl_support.resolve()): sha256(args.icl_support),
        str(args.satellite_profiles.resolve()): sha256(args.satellite_profiles),
        str(association_path.resolve()): sha256(association_path),
        str(blocks_path.resolve()): sha256(blocks_path),
    })
    result["output_sha256"] = {
        str(association_path.resolve()): sha256(association_path),
        str(blocks_path.resolve()): sha256(blocks_path),
        str(scope_path.resolve()): sha256(scope_path),
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=Path, default=DEFAULT_SUPPORT)
    parser.add_argument("--previous-74", type=Path, default=DEFAULT_PREVIOUS)
    parser.add_argument("--icl-support", type=Path, default=DEFAULT_ICL_SUPPORT)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--target-covariance", type=Path, default=DEFAULT_TARGET_COV)
    parser.add_argument("--target-control-covariance", type=Path, default=DEFAULT_TARGET_CONTROL_COV)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 block resamples are required")
    result = run(args)
    print(json.dumps({
        "status": result["status"], "science_measurement": result["science_measurement"],
        "scope": result["support_scope"],
        "filters": {band: values["point_estimate"] for band, values in result["filters"].items()},
        "blockwidth_sensitivity": {
            width: {band: {
                "corrected_ci95": item["pixel_noise_corrected_partial_r_bootstrap_ci95"],
                "blocks": item["n_connected_blocks"],
            } for band, item in by_band.items()}
            for width, by_band in result["blockwidth_sensitivity"]["results_by_width_and_filter"].items()
        },
        "output": str(args.output),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
