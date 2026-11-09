#!/usr/bin/env python3
"""Audit the 200--300 kpc excess against the 300--500 kpc outer profile.

Uses the previously measured 800-target/798-control annular profile archive.
It reconstructs the fixed 50--500 kpc complete-case pair samples, validates
the saved selection and means, then block-bootstraps the area-weighted
200--300 minus 300--500 kpc surface-brightness contrast. This is exploratory;
it does not identify pure ICL or promote science_measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PROFILE = DATA / "prestack_flux_profiles_2026-10-18_02"
DEFAULT_CONTROL = DATA / "prestack_extraction_2026-10-15_01"
DEFAULT_OUTPUT = DATA / "icl_radial_outer_baseline_contrast_2026-11-06_01"
CONFIG = "k2.5_c40_r50-500"
FACTOR_INDEX = 0  # FACTORS[0] = 2.5
EDGES_KPC = [50, 100, 150, 200, 300, 400, 500]
OUTER_CONTRAST_WEIGHTS = np.array([0., 0., 0., 1., -0.4375, -0.5625])
REPLICAS = 20000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_hash(path: Path, expected: str, name: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{name} hash differs from its source summary")
    return actual


def percentile(values: np.ndarray) -> tuple[float, float, float]:
    q = np.quantile(np.asarray(values, dtype=float), [0.025, 0.5, 0.975])
    return float(q[0]), float(q[1]), float(q[2])


def analyze(profile_dir: Path, control_dir: Path, output_dir: Path) -> dict:
    profile_dir, control_dir, output_dir = map(Path.resolve, (profile_dir, control_dir, output_dir))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    summary_path = profile_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete" or summary.get("n_targets_measured") != 800:
        raise ValueError("source full-sample profile run is incomplete")
    expected_outputs = summary.get("output_sha256", {})
    required = {
        "cluster_profiles.npz": profile_dir / "cluster_profiles.npz",
        "stack_covariances.npz": profile_dir / "stack_covariances.npz",
        "stack_profiles.tsv": profile_dir / "stack_profiles.tsv",
        "cluster_profiles.tsv": profile_dir / "cluster_profiles.tsv",
    }
    source_hashes = {str(summary_path): sha256(summary_path)}
    for name, path in required.items():
        expected = expected_outputs.get(str(path))
        if expected is None:
            raise ValueError(f"source summary does not pin {name}")
        source_hashes[str(path)] = verify_hash(path, expected, name)

    sample_path = DATA / "sample_manifest.json"
    extraction_path = control_dir / "extraction_summary.json"
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    if len(sample.get("clusters", [])) != 800 or extraction.get("n_measured") != 798:
        raise ValueError("target/control source manifests are incomplete")
    source_hashes[str(sample_path)] = sha256(sample_path)
    source_hashes[str(extraction_path)] = sha256(extraction_path)
    frames = sample["clusters"]
    target_ids = [str(row["id"]) for row in frames]
    if len(set(target_ids)) != 800:
        raise ValueError("duplicate target ID in sample manifest")
    records = extraction.get("records", [])
    controls = {str(row["target_cluster_id"]): row for row in records}
    if len(controls) != 798 or len(controls) != len(records):
        raise ValueError("duplicate or missing matched-control mapping")
    control_rows = []
    target_control_order = []
    for frame in frames:
        target_id = str(frame["id"])
        record = controls.get(target_id)
        if record is None:
            control_rows.append({"ra_deg": float(frame["ra_deg"]), "dec_deg": float(frame["dec_deg"]),
                                 "filter": ""})
            target_control_order.append(None)
            continue
        baseline = control_dir / record["control_id"] / "baseline.npz"
        source_hashes[str(baseline)] = verify_hash(baseline, record["baseline_sha256"],
                                                   f"control baseline {record['control_id']}")
        control_rows.append({"ra_deg": float(record["ra_deg"]), "dec_deg": float(record["dec_deg"]),
                             "filter": record["header_filter"]})
        target_control_order.append(baseline)

    with np.load(required["cluster_profiles.npz"], allow_pickle=False) as target_archive:
        target = {name: target_archive[name].copy() for name in
                  ("raw", "corrected", "coverage", "nominal_raw", "nominal_corrected")}
    expected_shape = (800, 3, 6)
    if target["raw"].shape != expected_shape or target["corrected"].shape != expected_shape:
        raise ValueError("target profile dimensions do not match the pinned six-annulus stack")
    control = {mode: np.full(expected_shape, np.nan) for mode in ("raw", "corrected")}
    control_coverage = np.full(expected_shape, np.nan)
    for index, baseline_path in enumerate(target_control_order):
        if baseline_path is None:
            continue
        with np.load(baseline_path, allow_pickle=False) as archive:
            for mode in control:
                control[mode][index] = archive[mode]
            control_coverage[index] = archive["coverage"]

    stack_rows = read_tsv(required["stack_profiles.tsv"])
    output_rows = []
    bootstrap_draws = {}
    band_support = {}
    with np.load(required["stack_covariances.npz"], allow_pickle=False) as archive:
        for band, seed in (("i", 2026110601), ("i2", 2026110602)):
            raw_key = f"{CONFIG}_{band}_raw__selected_indices"
            corrected_key = f"{CONFIG}_{band}_corrected__selected_indices"
            selected_indices = archive[raw_key].astype(int)
            if not np.array_equal(selected_indices, archive[corrected_key].astype(int)):
                raise ValueError(f"raw/corrected support differs in {band}")
            if selected_indices.ndim != 1 or len(selected_indices) < 2:
                raise ValueError(f"insufficient saved {band} profile support")
            selected = np.zeros(800, dtype=bool)
            selected[selected_indices] = True
            indices, labels, target_cells, control_cells, n_blocks, largest = pair_block_labels(
                frames, control_rows, selected)
            if not np.array_equal(indices, selected_indices):
                raise ValueError(f"reconstructed {band} complete-case selection changed")
            if n_blocks < 3:
                raise ValueError(f"too few connected sky blocks for {band} bootstrap")
            groups = np.unique(labels)
            group_n = np.array([(labels == group).sum() for group in groups], dtype=float)
            rng = np.random.default_rng(seed)
            draws = rng.multinomial(len(groups), np.full(len(groups), 1 / len(groups)),
                                    size=REPLICAS)
            band_support[band] = {
                "n_pairs": int(len(indices)), "n_connected_blocks": int(n_blocks),
                "n_target_cells": int(target_cells), "n_control_cells": int(control_cells),
                "largest_block_pairs": int(largest),
                "block_sizes": {str(group): int((labels == group).sum()) for group in groups},
            }
            band_contrast_draws = {}
            for mode in ("raw", "corrected"):
                values = target[mode][indices, FACTOR_INDEX, :] - control[mode][indices, FACTOR_INDEX, :]
                if not np.isfinite(values).all():
                    raise ValueError(f"non-finite pair profile in {band}/{mode}")
                mean = values.mean(axis=0)
                contrast = float(OUTER_CONTRAST_WEIGHTS @ mean)
                group_sums = np.array([values[labels == group].sum(axis=0) for group in groups])
                denom = draws @ group_n
                means = (draws @ group_sums) / denom[:, None]
                draws_contrast = means @ OUTER_CONTRAST_WEIGHTS
                band_contrast_draws[mode] = draws_contrast

                saved_mean = archive[f"{CONFIG}_{band}_{mode}__mean"]
                if not np.allclose(mean, saved_mean, rtol=1e-10, atol=1e-10):
                    raise ValueError(f"reconstructed {band}/{mode} mean differs from saved stack")
                saved = next((row for row in stack_rows
                              if row["config"] == CONFIG and row["filter"] == band
                              and row["extraction"] == mode), None)
                if saved is None or int(saved["n_pairs"]) != len(indices):
                    raise ValueError(f"saved stack table lacks the {band}/{mode} sample")
                saved_values = [float(row["mean_pair_difference_njy_arcsec2"])
                                for row in sorted((row for row in stack_rows
                                                   if row["config"] == CONFIG and row["filter"] == band
                                                   and row["extraction"] == mode),
                                                  key=lambda row: float(row["r_min_kpc"]))]
                if not np.allclose(mean, saved_values, rtol=1e-9, atol=1e-9):
                    raise ValueError(f"reconstructed {band}/{mode} annular means differ from TSV")
                lo, med, hi = percentile(draws_contrast)
                output_rows.append({
                    "filter": band, "extraction": mode, "n_pairs": int(len(indices)),
                    "n_connected_blocks": int(n_blocks),
                    "delta_200_300_njy_arcsec2": float(mean[3]),
                    "delta_300_400_njy_arcsec2": float(mean[4]),
                    "delta_400_500_njy_arcsec2": float(mean[5]),
                    "outer_area_weighted_sb_njy_arcsec2": float(.4375 * mean[4] + .5625 * mean[5]),
                    "delta_200_300_minus_outer_njy_arcsec2": contrast,
                    "block_bootstrap_ci95_low": lo,
                    "block_bootstrap_median": med,
                    "block_bootstrap_ci95_high": hi,
                    "interpretation": "exploratory radial contrast; not a pure ICL estimate",
                })
            bootstrap_draws[band + "_raw"] = band_contrast_draws["raw"]
            bootstrap_draws[band + "_corrected"] = band_contrast_draws["corrected"]

    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "outer_baseline_contrast.tsv"
    fields = ["filter", "extraction", "n_pairs", "n_connected_blocks",
              "delta_200_300_njy_arcsec2", "delta_300_400_njy_arcsec2",
              "delta_400_500_njy_arcsec2", "outer_area_weighted_sb_njy_arcsec2",
              "delta_200_300_minus_outer_njy_arcsec2", "block_bootstrap_ci95_low",
              "block_bootstrap_median", "block_bootstrap_ci95_high", "interpretation"]
    write_tsv(table_path, output_rows, fields)
    draws_path = output_dir / "block_bootstrap_outer_contrast_draws.npz"
    np.savez_compressed(draws_path, **bootstrap_draws)
    result = {
        "schema": "astrafek-icl-radial-outer-baseline-contrast-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "image_pixels_read": False,
        "source_profile_directory": str(profile_dir), "control_baseline_directory": str(control_dir),
        "source_hashes": source_hashes,
        "radial_bins_kpc": EDGES_KPC,
        "contrast_definition": "mean target-minus-control SB in 200--300 kpc minus area-weighted mean SB in 300--500 kpc; annular area weights are 0.4375 and 0.5625",
        "coverage_definition": "reuse the pinned 2.5 Kron-mask, 40% common-coverage, six-annulus complete-case samples from the existing exploratory 50--500 kpc stack",
        "bootstrap": {"method": "resample connected paired target/control 5-degree sky blocks with replacement; pair-weight within each draw",
                      "replicas": REPLICAS, "seeds": {"i": 2026110601, "i2": 2026110602}},
        "support": band_support, "results": output_rows,
        "limitations": [
            "20 i and 48 i2 pairs remain after the six-annulus 40%-coverage gate; this is strongly selected support",
            "the 300--500 kpc signal may include real extended cluster light and is not a blank-sky estimate",
            "the 40% coverage floor and post-hoc radial contrast are exploratory, not a preregistered science cut",
            "BCG/diffuse-light separation, source completeness, source-Poisson, exact coadd covariance, sky transfer, and extended PSF wings remain unresolved",
            "raw and outer-plane-corrected profiles share images and are sensitivity estimates, not independent replications",
        ],
        "result_files_sha256": {str(table_path): sha256(table_path), str(draws_path): sha256(draws_path)},
    }
    summary_out = output_dir / "summary.json"
    summary_out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"output": str(output_dir), "support": band_support, "results": output_rows,
                      "summary_sha256": sha256(summary_out)}, ensure_ascii=False, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--controls", type=Path, default=DEFAULT_CONTROL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    analyze(args.profiles, args.controls, args.output)


if __name__ == "__main__":
    main()
