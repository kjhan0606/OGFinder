#!/usr/bin/env python3
"""Relate the six-annulus ICL contrast to satellite-centered outer light.

The analysis joins existing 50--500 proper-kpc ICL profile measurements to
the existing 50--75 kpc satellite-centered residual profile by cluster ID and
filter. It uses only measured intersections and resamples connected paired
5-degree sky blocks. It makes no HSC requests and does not read image pixels.
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
DEFAULT_SATELLITE = DATA / "satellite_outer_profile_perbin_cov040_2026-11-05_01"
DEFAULT_ICL_CONTRAST = DATA / "icl_radial_outer_baseline_contrast_2026-11-06_01"
DEFAULT_OUTPUT = DATA / "icl_satellite_outer_profile_relation_2026-11-06_01"
CONFIG = "k2.5_c40_r50-500"
FACTOR_INDEX = 0
OUTER_WEIGHTS = np.array([0., 0., 0., 1., -0.4375, -0.5625])
BOOTSTRAP_REPLICAS = 20000
BOOTSTRAP_SEEDS = {"i": 2026110611, "i2": 2026110612}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def weighted_partial_r(icl: np.ndarray, satellite: np.ndarray, z: np.ndarray,
                       richness: np.ndarray, weights: np.ndarray) -> float:
    covariates = np.column_stack((np.ones(len(z)), z, np.log(richness)))
    sw = np.sqrt(weights)
    design = covariates * sw[:, None]
    rx = icl - covariates @ np.linalg.lstsq(design, icl * sw, rcond=None)[0]
    ry = satellite - covariates @ np.linalg.lstsq(design, satellite * sw, rcond=None)[0]
    mx = np.average(rx, weights=weights)
    my = np.average(ry, weights=weights)
    dx, dy = rx - mx, ry - my
    vx, vy = np.average(dx * dx, weights=weights), np.average(dy * dy, weights=weights)
    if vx <= 0 or vy <= 0:
        return float("nan")
    return float(np.average(dx * dy, weights=weights) / np.sqrt(vx * vy))


def percentile(values: np.ndarray) -> list[float] | None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return None
    return [float(x) for x in np.quantile(finite, [0.025, 0.5, 0.975])]


def analyze(profile_dir: Path, control_dir: Path, satellite_dir: Path,
            contrast_dir: Path, output_dir: Path) -> dict:
    profile_dir, control_dir, satellite_dir, contrast_dir, output_dir = map(
        Path.resolve, (profile_dir, control_dir, satellite_dir, contrast_dir, output_dir))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    contrast_summary_path = contrast_dir / "summary.json"
    contrast_summary = json.loads(contrast_summary_path.read_text(encoding="utf-8"))
    if (contrast_summary.get("schema") != "astrafek-icl-radial-outer-baseline-contrast-v1"
            or contrast_summary.get("status") != "complete_exploratory"
            or contrast_summary.get("science_measurement") is not False):
        raise ValueError("parent ICL radial-contrast analysis is incomplete")
    for text, digest in contrast_summary.get("source_hashes", {}).items():
        path = Path(text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"parent ICL contrast source hash changed: {path}")

    profile_summary_path = profile_dir / "summary.json"
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    satellite_summary_path = satellite_dir / "summary.json"
    satellite_summary = json.loads(satellite_summary_path.read_text(encoding="utf-8"))
    satellite_path = satellite_dir / "cluster_profiles.tsv"
    if (profile_summary.get("status") != "complete" or profile_summary.get("n_targets_measured") != 800
            or satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("sample_scope") != "full_manifest"
            or satellite_summary.get("n_frames_completed") != 800):
        raise ValueError("parent ICL or satellite profile run is incomplete")
    expected_satellite_hash = satellite_summary.get("output_sha256", {}).get(str(satellite_path))
    if not expected_satellite_hash or sha256(satellite_path) != expected_satellite_hash:
        raise ValueError("satellite cluster profile table hash differs from its summary")

    sample_path = DATA / "sample_manifest.json"
    extraction_path = control_dir / "extraction_summary.json"
    sample = json.loads(sample_path.read_text(encoding="utf-8"))
    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    frames = sample.get("clusters", [])
    if len(frames) != 800 or extraction.get("n_measured") != 798:
        raise ValueError("target/control manifest coverage changed")
    target_rows = read_tsv(profile_dir / "cluster_profiles.tsv")
    targets = {str(row["id"]): row for row in target_rows}
    if len(targets) != 800:
        raise ValueError("target profile ledger is incomplete")
    control_records = extraction.get("records", [])
    controls = {str(row["target_cluster_id"]): row for row in control_records}
    if len(controls) != 798 or len(controls) != len(control_records):
        raise ValueError("matched-control map is duplicated or incomplete")
    controls_by_index = []
    source_hashes = {
        str(sample_path): sha256(sample_path),
        str(extraction_path): sha256(extraction_path),
        str(profile_summary_path): sha256(profile_summary_path),
        str(satellite_summary_path): sha256(satellite_summary_path),
        str(satellite_path): sha256(satellite_path),
        str(contrast_summary_path): sha256(contrast_summary_path),
    }
    for frame in frames:
        record = controls.get(str(frame["id"]))
        if record is None:
            controls_by_index.append(None)
            continue
        baseline = control_dir / record["control_id"] / "baseline.npz"
        if not baseline.is_file() or sha256(baseline) != record["baseline_sha256"]:
            raise ValueError(f"control baseline hash mismatch: {baseline}")
        source_hashes[str(baseline)] = record["baseline_sha256"]
        controls_by_index.append((record, baseline))

    satellite_rows = read_tsv(satellite_path)
    satellite_by_key = {}
    for row in satellite_rows:
        key = (str(row["cluster_id"]), row["filter"], row["radial_bin_kpc"])
        if key in satellite_by_key:
            raise ValueError(f"duplicate satellite cluster profile: {key}")
        satellite_by_key[key] = row

    with np.load(profile_dir / "cluster_profiles.npz", allow_pickle=False) as archive:
        target_raw = archive["raw"].copy()
        target_corrected = archive["corrected"].copy()
    per_pair_rows = []
    result_rows = []
    draws_archive = {}
    support_summary = {}
    with np.load(profile_dir / "stack_covariances.npz", allow_pickle=False) as archive:
        for band in ("i", "i2"):
            selected = archive[f"{CONFIG}_{band}_raw__selected_indices"].astype(int)
            if not np.array_equal(selected, archive[f"{CONFIG}_{band}_corrected__selected_indices"].astype(int)):
                raise ValueError(f"raw/corrected ICL supports differ for {band}")
            icl_values = {mode: [] for mode in ("raw", "corrected")}
            matched_indices, sat_values, zs, richnesses = [], [], [], []
            row_map = {}
            for index in selected:
                frame = frames[int(index)]
                ident = str(frame["id"])
                target = targets[ident]
                control_entry = controls_by_index[int(index)]
                if (target["filter"] != band or control_entry is None
                        or control_entry[0]["header_filter"] != band):
                    raise ValueError(f"selected ICL pair/filter mismatch for {ident}/{band}")
                satellite = satellite_by_key.get((ident, band, "50-75"))
                if satellite is None:
                    continue
                baseline_path = control_entry[1]
                with np.load(baseline_path, allow_pickle=False) as baseline:
                    raw_diff = target_raw[int(index), FACTOR_INDEX] - baseline["raw"][FACTOR_INDEX]
                    corrected_diff = target_corrected[int(index), FACTOR_INDEX] - baseline["corrected"][FACTOR_INDEX]
                raw_contrast = float(OUTER_WEIGHTS @ raw_diff)
                corrected_contrast = float(OUTER_WEIGHTS @ corrected_diff)
                row = {
                    "cluster_id": ident, "filter": band, "z": float(frame["z"]),
                    "richness": float(frame["richness"]),
                    "satellite_outer_delta_sb_50_75_njy_arcsec2": float(satellite["mean_delta_sb_njy_arcsec2"]),
                    "n_satellites_in_50_75_bin": int(satellite["n_satellites"]),
                    "icl_contrast_raw_njy_arcsec2": raw_contrast,
                    "icl_contrast_corrected_njy_arcsec2": corrected_contrast,
                    "icl_delta_200_300_raw_njy_arcsec2": float(raw_diff[3]),
                    "icl_delta_300_400_raw_njy_arcsec2": float(raw_diff[4]),
                    "icl_delta_400_500_raw_njy_arcsec2": float(raw_diff[5]),
                }
                row_map[ident] = row
                per_pair_rows.append(row)
                matched_indices.append(int(index))
                sat_values.append(float(satellite["mean_delta_sb_njy_arcsec2"]))
                zs.append(float(frame["z"]))
                richnesses.append(float(frame["richness"]))
            if len(matched_indices) < 8:
                raise ValueError(f"too few complete ICL/satellite clusters for {band}")
            sat_values = np.asarray(sat_values)
            zs = np.asarray(zs)
            richnesses = np.asarray(richnesses)
            index_mask = np.zeros(800, dtype=bool)
            index_mask[matched_indices] = True
            indices, labels, ntarget, ncontrol, nblocks, largest = pair_block_labels(
                frames, [{"ra_deg": (entry[0]["ra_deg"] if entry else frames[i]["ra_deg"]),
                          "dec_deg": (entry[0]["dec_deg"] if entry else frames[i]["dec_deg"])}
                         for i, entry in enumerate(controls_by_index)], index_mask)
            if list(indices) != matched_indices:
                raise ValueError(f"paired sky-block order changed for {band}")
            groups = np.unique(labels)
            group_n = np.array([(labels == group).sum() for group in groups], dtype=float)
            rng = np.random.default_rng(BOOTSTRAP_SEEDS[band])
            block_draws = rng.multinomial(nblocks, np.full(nblocks, 1 / nblocks),
                                          size=BOOTSTRAP_REPLICAS)
            row_blocks = np.asarray(labels)
            weights = np.ones(len(indices), dtype=float)
            support_summary[band] = {
                "n_icl_complete_pairs": int(len(selected)),
                "n_satellite_50_75_supported_pairs": int(len(matched_indices)),
                "n_excluded_for_missing_satellite_profile": int(len(selected) - len(matched_indices)),
                "n_connected_paired_sky_blocks": int(nblocks),
                "n_target_cells": int(ntarget), "n_control_cells": int(ncontrol),
                "largest_block_pairs": int(largest),
                "satellite_selection": "existing detected-member 50--75 proper-kpc profile with its original 40% annular coverage selection",
            }
            mode_draws = {}
            for mode in ("raw", "corrected"):
                x = np.asarray([row_map[str(frames[i]["id"])][f"icl_contrast_{mode}_njy_arcsec2"]
                                for i in matched_indices])
                observed = weighted_partial_r(x, sat_values, zs, richnesses, weights)
                boot = np.full(BOOTSTRAP_REPLICAS, np.nan)
                for j, multiplicities in enumerate(block_draws):
                    row_weights = np.asarray([multiplicities[np.flatnonzero(groups == label)[0]]
                                              for label in row_blocks], dtype=float)
                    if np.count_nonzero(row_weights) < 8:
                        continue
                    boot[j] = weighted_partial_r(x, sat_values, zs, richnesses, row_weights)
                ci = percentile(boot)
                mode_draws[mode] = boot
                result_rows.append({
                    "filter": band, "icl_extraction": mode,
                    "n_clusters": int(len(matched_indices)), "n_connected_blocks": int(nblocks),
                    "partial_correlation_control_z_logrichness": observed,
                    "block_bootstrap_ci95_low": ci[0] if ci else None,
                    "block_bootstrap_median": ci[1] if ci else None,
                    "block_bootstrap_ci95_high": ci[2] if ci else None,
                    "valid_block_bootstrap_replicates": int(np.isfinite(boot).sum()),
                    "nominal_p_value": None,
                    "interpretation": "exploratory selected-sample association; no measurement-error propagation",
                })
            draws_archive[f"{band}_raw"] = mode_draws["raw"]
            draws_archive[f"{band}_corrected"] = mode_draws["corrected"]

    output_dir.mkdir(parents=True, exist_ok=False)
    pair_path = output_dir / "icl_satellite_common_support_pairs.tsv"
    pair_fields = ["cluster_id", "filter", "z", "richness",
                   "satellite_outer_delta_sb_50_75_njy_arcsec2", "n_satellites_in_50_75_bin",
                   "icl_contrast_raw_njy_arcsec2", "icl_contrast_corrected_njy_arcsec2",
                   "icl_delta_200_300_raw_njy_arcsec2", "icl_delta_300_400_raw_njy_arcsec2",
                   "icl_delta_400_500_raw_njy_arcsec2"]
    write_tsv(pair_path, per_pair_rows, pair_fields)
    result_path = output_dir / "association_results.tsv"
    write_tsv(result_path, result_rows, ["filter", "icl_extraction", "n_clusters", "n_connected_blocks",
               "partial_correlation_control_z_logrichness", "block_bootstrap_ci95_low",
               "block_bootstrap_median", "block_bootstrap_ci95_high", "valid_block_bootstrap_replicates",
               "nominal_p_value", "interpretation"])
    draws_path = output_dir / "block_bootstrap_partial_correlations.npz"
    np.savez_compressed(draws_path, **draws_archive)
    result = {
        "schema": "astrafek-icl-satellite-outer-profile-relation-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "image_pixels_read": False,
        "icl_parent_contrast_summary_sha256": sha256(contrast_summary_path),
        "satellite_profile_summary_sha256": sha256(satellite_summary_path),
        "satellite_profile_table_sha256": sha256(satellite_path),
        "predictor": "measured satellite-centered 50--75 proper-kpc target-minus-azimuthal-control residual profile, averaged by cluster",
        "outcome": "ICL target-minus-matched-control 200--300 kpc contrast minus area-weighted 300--500 kpc contrast",
        "covariates": ["cluster redshift", "log CAMIRA richness"],
        "partial_correlation": "Pearson correlation of covariate-residualized cluster-level outcome and predictor",
        "bootstrap": {"method": "connected paired target/control 5-degree block resampling with replacement",
                      "replicas": BOOTSTRAP_REPLICAS, "seeds": BOOTSTRAP_SEEDS},
        "support": support_summary,
        "results": result_rows,
        "source_hashes": source_hashes,
        "output_sha256": {str(pair_path): sha256(pair_path), str(result_path): sha256(result_path),
                          str(draws_path): sha256(draws_path)},
        "limitations": [
            "the ICL sample requires complete 50--500 kpc support at a post-hoc 40% coverage floor",
            "satellite profiles use detected CAMIRA member matches and a separate 40% annular coverage cut",
            "the satellite azimuthal-control residual does not isolate satellite stellar envelopes from asymmetric ICL or PSF wings",
            "the 300--500 kpc ICL profile can contain real cluster light and is not a blank-sky estimate",
            "the block bootstrap omits measurement uncertainty, PSF/sky transfer uncertainty, and source-Poisson covariance",
            "raw and outer-plane-corrected ICL profiles are shared-data sensitivity estimates, not independent replications",
        ],
    }
    summary_out = output_dir / "summary.json"
    summary_out.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"output": str(output_dir), "support": support_summary,
                      "results": result_rows, "summary_sha256": sha256(summary_out)},
                     ensure_ascii=False, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--controls", type=Path, default=DEFAULT_CONTROL)
    parser.add_argument("--satellite", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--icl-contrast", type=Path, default=DEFAULT_ICL_CONTRAST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    analyze(args.profiles, args.controls, args.satellite, args.icl_contrast, args.output)


if __name__ == "__main__":
    main()
