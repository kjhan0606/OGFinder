#!/usr/bin/env python3
"""Verify an independent-support ICL sensitivity output against its source tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def equal_number(actual: str, expected: float | None, *, atol: float = 1e-8) -> bool:
    if expected is None:
        return actual in ("", "None")
    try:
        return math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=atol)
    except (TypeError, ValueError):
        return False


def verify(work: Path) -> dict:
    summary_path = work / "summary.json"
    pair_path = work / "paired_independent_support_profiles.tsv"
    radial_path = work / "independent_support_radial_profiles.tsv"
    summary = read_json(summary_path)
    require = lambda test, message: (_ for _ in ()).throw(ValueError(message)) if not test else None
    require(summary.get("schema") == "astrafek-hsc-icl-independent-support-v1",
            "unexpected independent-support schema")
    require(summary.get("status") == "complete_exploratory", "independent-support run is incomplete")
    require(summary.get("science_measurement") is False, "science gate must remain false")
    require(summary.get("n_targets") == 800 and summary.get("n_matched_controls") == 798,
            "source sample counts differ from the approved inputs")
    require(summary.get("n_pairs") == 783 and summary.get("n_filter_mismatch_excluded") == 15,
            "pair/filter counts differ from the expected matched sample")
    for source, digest in summary.get("scripts_sha256", {}).items():
        path = Path(source)
        require(path.is_file() and sha256(path) == digest,
                f"run-time source checksum mismatch: {path}")

    common = summary["common_support_reference"]
    common_work = Path(common["work"])
    for filename, key in (("summary.json", "summary_sha256"),
                          ("independent_verification.json", "verification_sha256"),
                          ("paired_radial_profiles.tsv", "pair_profile_sha256")):
        path = common_work / filename
        require(path.is_file() and sha256(path) == common[key],
                f"common-support provenance mismatch: {path}")
    common_summary = read_json(common_work / "summary.json")
    common_verification = read_json(common_work / "independent_verification.json")
    common_pairs_path = common_work / "paired_radial_profiles.tsv"
    require(common_pairs_path.is_file()
            and sha256(common_pairs_path) == common_summary.get("outputs", {}).get("paired_profiles_sha256"),
            "common-support annular pair profiles fail their saved checksum")
    require(common_summary.get("status") == "complete_exploratory"
            and common_summary.get("science_measurement") is False
            and common_verification.get("status") == "verified",
            "common-support baseline is not complete and verified")
    require(float(common_summary.get("mask_factor", 2.5)) == float(summary["mask_factor"]),
            "mask factors differ between estimators")

    pair_rows = read_tsv(pair_path)
    radial_rows = read_tsv(radial_path)
    require(sha256(pair_path) == summary["outputs"]["pair_profiles_sha256"],
            "pair-profile table hash mismatch")
    require(sha256(radial_path) == summary["outputs"]["radial_profiles_sha256"],
            "radial-profile table hash mismatch")
    require(len(pair_rows) == summary["n_pairs"], "pair-profile row count mismatch")
    pair_map = {}
    for row in pair_rows:
        key = (row["target_id"], row["control_id"], row["filter"])
        require(key not in pair_map, f"duplicate target/control pair: {key}")
        profile = json.loads(row["profile_json"])
        require(len(profile) == 20, f"expected 20 annuli for pair {key}")
        for annulus in profile:
            nt, nc, ncommon = (int(annulus[key]) for key in
                               ("n_target_pixels", "n_control_pixels", "n_common_pixels"))
            geom = int(annulus["geometric_pixels"])
            require(0 <= ncommon <= min(nt, nc) <= geom,
                    f"invalid support counts for pair {key}")
            require(math.isclose(float(annulus["target_coverage_fraction"]), nt / geom, abs_tol=1e-12)
                    and math.isclose(float(annulus["control_coverage_fraction"]), nc / geom, abs_tol=1e-12)
                    and math.isclose(float(annulus["common_coverage_fraction"]), ncommon / geom, abs_tol=1e-12),
                    f"coverage fractions disagree with counts for pair {key}")
            target, control, delta = (annulus["target_mean_njy_arcsec2"],
                                      annulus["control_mean_njy_arcsec2"],
                                      annulus["delta_mean_njy_arcsec2"])
            if nt and nc:
                require(target is not None and control is not None and delta is not None
                        and math.isclose(float(target) - float(control), float(delta),
                                         rel_tol=1e-10, abs_tol=1e-8),
                        f"independent subtraction identity failed for pair {key}")
            else:
                require(delta is None, f"unsupported pair has a delta for pair {key}")
        pair_map[key] = profile

    common_rows = read_tsv(common_work / "paired_radial_profiles.tsv")
    common_map = {}
    for row in common_rows:
        key = (row["target_id"], row["control_id"], row["filter"])
        require(key not in common_map, f"duplicate common-support pair: {key}")
        profile = json.loads(row["profile_json"])
        require(len(profile) == 20, f"expected 20 common-support annuli for pair {key}")
        common_map[key] = profile
    require(set(common_map) == set(pair_map),
            "independent and common-support pair sets differ")
    require(len(radial_rows) == 40, "expected 20 radial bins in each filter")
    radial_map = {}
    for row in radial_rows:
        key = (row["filter"], int(row["r_min_kpc"]), int(row["r_max_kpc"]))
        require(key not in radial_map, f"duplicate aggregate radial row: {key}")
        radial_map[key] = row
    for band in ("i", "i2"):
        for bin_index in range(20):
            lo, hi = bin_index * 25, (bin_index + 1) * 25
            key = (band, lo, hi)
            independent_values, shared_independent, shared_common = [], [], []
            n_target = n_control = 0
            target_cov, control_cov, common_cov = [], [], []
            for pair_key, profile in pair_map.items():
                if pair_key[2] != band:
                    continue
                item = profile[bin_index]
                n_target += int(item["n_target_pixels"] > 0)
                n_control += int(item["n_control_pixels"] > 0)
                target_cov.append(float(item["target_coverage_fraction"]))
                control_cov.append(float(item["control_coverage_fraction"]))
                common_cov.append(float(item["common_coverage_fraction"]))
                value = item["delta_mean_njy_arcsec2"]
                if value is None:
                    continue
                independent_values.append(float(value))
                ref = common_map[(pair_key[0], pair_key[1], band)][bin_index].get(
                    "delta_mean_njy_arcsec2")
                if ref is not None:
                    shared_independent.append(float(value))
                    shared_common.append(float(ref))
            row = radial_map.get(key)
            require(row is not None, f"missing aggregate radial bin {key}")
            ncommon = len(shared_common)
            require(int(row["n_pairs_independent_delta"]) == len(independent_values)
                    and int(row["n_pairs_target_support"]) == n_target
                    and int(row["n_pairs_control_support"]) == n_control
                    and int(row["n_pairs_common_support"]) == ncommon
                    and int(row["n_pairs_added_vs_common"]) == len(independent_values) - ncommon,
                    f"pair support counts disagree in {key}")

            def mean(values):
                return float(np.mean(values)) if values else None

            def sem(values):
                return float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else None

            values_by_key = {
                "mean_target_independent_njy_arcsec2": None,
                "mean_control_independent_njy_arcsec2": None,
                "mean_delta_independent_njy_arcsec2": mean(independent_values),
                "pair_sem_delta_independent_njy_arcsec2": sem(independent_values),
                "mean_independent_delta_on_common_pairs_njy_arcsec2": mean(shared_independent),
                "mean_common_delta_on_same_pairs_njy_arcsec2": mean(shared_common),
                "pair_sem_independent_delta_on_common_pairs_njy_arcsec2": sem(shared_independent),
                "pair_sem_common_delta_on_same_pairs_njy_arcsec2": sem(shared_common),
                "median_target_coverage": float(np.median(target_cov)) if target_cov else None,
                "median_control_coverage": float(np.median(control_cov)) if control_cov else None,
                "median_common_coverage": float(np.median(common_cov)) if common_cov else None,
            }
            target_values = [float(profile[bin_index]["target_mean_njy_arcsec2"])
                             for pair_key, profile in pair_map.items()
                             if pair_key[2] == band and profile[bin_index]["delta_mean_njy_arcsec2"] is not None]
            control_values = [float(profile[bin_index]["control_mean_njy_arcsec2"])
                              for pair_key, profile in pair_map.items()
                              if pair_key[2] == band and profile[bin_index]["delta_mean_njy_arcsec2"] is not None]
            values_by_key["mean_target_independent_njy_arcsec2"] = mean(target_values)
            values_by_key["mean_control_independent_njy_arcsec2"] = mean(control_values)
            for column, expected in values_by_key.items():
                require(equal_number(row[column], expected), f"aggregate {column} disagrees in {key}")

    base_inputs = (
        Path(summary["sample"]) / "sample_manifest.json",
        Path(summary["sample"]) / "prestack_cluster_masks_2026-10-16_01/summary.json",
        Path(summary["sample"]) / "prestack_cluster_masks_2026-10-16_01/records.tsv",
        Path(summary["sample"]) / "control_download_01/download_manifest.json",
        Path(summary["sample"]) / "control_download_01/download_status.json",
    )
    for path in base_inputs:
        expected = summary["input_sha256"].get(str(path))
        require(path.is_file() and expected is not None and sha256(path) == expected,
                f"manifest input hash mismatch: {path}")
    base_paths = {str(path) for path in base_inputs}
    image_input_paths = [path for path in summary["input_sha256"] if path not in base_paths]
    require(len(image_input_paths) == 3 * len(pair_rows), "expected one target, mask, and control FITS per pair")
    for path_text in image_input_paths:
        path = Path(path_text)
        require(path.is_file(), f"input FITS is missing: {path}")
    result = {
        "schema": "astrafek-hsc-icl-independent-support-verification-v1",
        "status": "verified", "summary_sha256": sha256(summary_path),
        "n_pairs": len(pair_rows), "n_by_filter": summary["n_by_filter"],
        "n_radial_rows": len(radial_rows), "manifest_paths_rehashed": len(base_inputs),
        "image_paths_exist": len(image_input_paths),
        "common_profile_rows_checked": len(common_rows),
        "science_measurement": False,
        "note": "Verifies input hashes, per-pair support geometry, pairwise target-control subtraction, and independently recomputed radial aggregation. It does not re-read all source FITS pixel arrays.",
    }
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    work = args.work.resolve()
    result = verify(work)
    output = (args.output.resolve() if args.output else work / "independent_verification.json")
    if output.exists():
        parser.error(f"refusing to overwrite verification: {output}")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"verification_sha256={sha256(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
