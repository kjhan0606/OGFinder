#!/usr/bin/env python3
"""Measure an exploratory ICL profile with separate target/control pixel support.

For every existing matched cluster/control pair, this computes an annular mean
from all valid target pixels and independently from all valid control pixels,
then differences those means. It retains a same-pair common-support comparison
to separate area selection from the additional pairs admitted by the estimator.
No network access or data download is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import stack_hsc_icl_target_control_2d as stack


DEFAULT_SAMPLE = stack.DEFAULT_SAMPLE
DEFAULT_COMMON_WORKS = stack.DEFAULT_WORKS
DEFAULT_WORKS = {
    2.5: DEFAULT_SAMPLE / "icl_2d_independent_support_k2p5_2026-10-19_01",
    3.5: DEFAULT_SAMPLE / "icl_2d_independent_support_k3p5_2026-10-19_01",
    5.0: DEFAULT_SAMPLE / "icl_2d_independent_support_k5p0_2026-10-19_01",
}


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def profile_by_pair(common_rows: list[dict[str, str]]) -> dict[tuple, list[dict]]:
    result = {}
    for row in common_rows:
        key = (row["target_id"], row["control_id"], row["filter"])
        if key in result:
            raise ValueError(f"duplicate common-support pair: {key}")
        profile = json.loads(row["profile_json"])
        if len(profile) != len(stack.PROFILE_EDGES_KPC) - 1:
            raise ValueError(f"unexpected common-support annulus count for pair {key}")
        result[key] = profile
    return result


def summarize_pairs(pair_rows: list[dict], common_profiles: dict[tuple, list[dict]]) -> list[dict]:
    output = []
    for band in stack.FILTERS:
        selected_pairs = [row for row in pair_rows if row["filter"] == band]
        for bin_index, (lo, hi) in enumerate(zip(stack.PROFILE_EDGES_KPC[:-1],
                                                  stack.PROFILE_EDGES_KPC[1:])):
            supported = []
            shared = []
            n_target = n_control = 0
            target_cov, control_cov, common_cov = [], [], []
            for pair in selected_pairs:
                profile = pair["independent_profile"][bin_index]
                n_target += int(profile["n_target_pixels"] > 0)
                n_control += int(profile["n_control_pixels"] > 0)
                target_cov.append(float(profile["target_coverage_fraction"]))
                control_cov.append(float(profile["control_coverage_fraction"]))
                common_cov.append(float(profile["common_coverage_fraction"]))
                if profile["delta_mean_njy_arcsec2"] is None:
                    continue
                row = {"pair": pair, "profile": profile}
                supported.append(row)
                common = common_profiles.get((pair["id"], pair["control_id"], band))
                common_delta = (common[bin_index].get("delta_mean_njy_arcsec2")
                                if common is not None else None)
                if common_delta is not None:
                    row["common_delta"] = float(common_delta)
                    shared.append(row)

            def stat(rows, key):
                values = np.asarray([float(row["profile"][key]) for row in rows], dtype=float)
                return {
                    "mean": float(values.mean()) if len(values) else None,
                    "sd": float(values.std(ddof=1)) if len(values) > 1 else None,
                    "sem": float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else None,
                }

            delta = stat(supported, "delta_mean_njy_arcsec2")
            target = stat(supported, "target_mean_njy_arcsec2")
            control = stat(supported, "control_mean_njy_arcsec2")
            shared_independent = stat(shared, "delta_mean_njy_arcsec2")
            common_delta = np.asarray([row["common_delta"] for row in shared], dtype=float)
            output.append({
                "filter": band, "r_min_kpc": int(lo), "r_max_kpc": int(hi),
                "n_pairs_available": len(selected_pairs),
                "n_pairs_target_support": n_target, "n_pairs_control_support": n_control,
                "n_pairs_independent_delta": len(supported),
                "n_pairs_common_support": len(shared),
                "n_pairs_added_vs_common": len(supported) - len(shared),
                "mean_target_independent_njy_arcsec2": target["mean"],
                "mean_control_independent_njy_arcsec2": control["mean"],
                "mean_delta_independent_njy_arcsec2": delta["mean"],
                "pair_sd_delta_independent_njy_arcsec2": delta["sd"],
                "pair_sem_delta_independent_njy_arcsec2": delta["sem"],
                "median_target_coverage": float(np.median(target_cov)) if target_cov else None,
                "median_control_coverage": float(np.median(control_cov)) if control_cov else None,
                "median_common_coverage": float(np.median(common_cov)) if common_cov else None,
                "mean_independent_delta_on_common_pairs_njy_arcsec2": shared_independent["mean"],
                "mean_common_delta_on_same_pairs_njy_arcsec2": float(common_delta.mean()) if len(common_delta) else None,
                "pair_sem_independent_delta_on_common_pairs_njy_arcsec2": shared_independent["sem"],
                "pair_sem_common_delta_on_same_pairs_njy_arcsec2": (
                    float(common_delta.std(ddof=1) / math.sqrt(len(common_delta)))
                    if len(common_delta) > 1 else None),
            })
    return output


def run(sample: Path, common_work: Path, work: Path, mask_factor: float,
        kpc_max: float = 500.0, n_side: int = 201) -> dict:
    if mask_factor not in stack.MASK_FACTORS:
        raise ValueError(f"mask factor must be one of {stack.MASK_FACTORS}")
    if work.exists():
        raise FileExistsError(f"refusing to overwrite existing work directory: {work}")
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("full-sample independent-support processing must run inside Slurm")

    sample = sample.resolve()
    common_work = common_work.resolve()
    sample_manifest_path = sample / "sample_manifest.json"
    mask_dir = sample / "prestack_cluster_masks_2026-10-16_01"
    mask_summary_path = mask_dir / "summary.json"
    mask_records_path = mask_dir / "records.tsv"
    control_manifest_path = sample / "control_download_01/download_manifest.json"
    control_status_path = sample / "control_download_01/download_status.json"
    sample_doc = stack.read_json(sample_manifest_path)
    mask_summary = stack.read_json(mask_summary_path)
    targets = sample_doc.get("clusters", [])
    masks = {str(row["id"]): row for row in stack.read_tsv(mask_records_path)}
    target_by_id = {str(row["id"]): row for row in targets}
    if len(targets) != 800 or len(target_by_id) != 800 or len(masks) != 800:
        raise ValueError("expected 800 unique targets and target-mask records")
    if (mask_summary.get("status") != "complete" or mask_summary.get("n_success") != 800
            or mask_summary.get("science_measurement") is not False
            or stack.sha256(mask_records_path) != mask_summary.get("output_sha256", {}).get(str(mask_records_path))):
        raise ValueError("target-mask inputs are incomplete or fail their saved checksum")

    controls_doc, downloads, control_base, control_hashes = stack.load_controls(sample)
    controls_by_target = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    downloads_by_id = {str(row["id"]): row for row in downloads.values()}
    if len(controls_by_target) != 798 or len(downloads_by_id) != 798:
        raise ValueError("expected 798 unique matched controls and downloaded images")
    if sample_doc.get("science_measurement") is not False or controls_doc.get("science_measurement") is not False:
        raise ValueError("input science gates are not explicitly false")

    common_summary = stack.read_json(common_work / "summary.json")
    common_verification = stack.read_json(common_work / "independent_verification.json")
    common_profile_path = common_work / "paired_radial_profiles.tsv"
    if (common_summary.get("status") != "complete_exploratory"
            or common_summary.get("science_measurement") is not False
            or float(common_summary.get("mask_factor", 2.5)) != mask_factor
            or common_verification.get("status") != "verified"
            or common_verification.get("science_measurement") is not False):
        raise ValueError("common-support reference is incomplete, unverified, or uses another mask")
    common_profiles = profile_by_pair(stack.read_tsv(common_profile_path))
    if len(common_profiles) != int(common_summary["n_pairs_stacked"]):
        raise ValueError("common-support reference has unexpected pair coverage")

    cosmology = stack.Cosmology(H0=float(sample_doc["h0"]), Om=float(sample_doc["omega_m"]))
    control_status_sha = control_hashes[str(control_status_path)]
    exclusions, pair_rows, input_hashes = [], [], {
        str(sample_manifest_path): stack.sha256(sample_manifest_path),
        str(mask_summary_path): stack.sha256(mask_summary_path),
        str(mask_records_path): stack.sha256(mask_records_path),
        str(control_manifest_path): control_hashes[str(control_manifest_path)],
        str(control_status_path): control_status_sha,
    }
    start_time = time.monotonic()
    for index, target in enumerate(targets, 1):
        ident = str(target["id"])
        control = controls_by_target.get(ident)
        if control is None:
            continue
        control_id = str(control["id"])
        download = downloads_by_id.get(control_id)
        if download is None:
            raise ValueError(f"{ident}: matched-control download record missing")
        target_band = masks[ident]["primary_filter"].strip().lower()
        control_band = str(control.get("coadd_header_filter", "")).strip().lower()
        if control_band != str(download.get("review", {}).get("coadd_header_filter", "")).strip().lower():
            raise ValueError(f"{ident}: control filter evidence disagrees")
        if target_band not in stack.FILTERS or control_band not in stack.FILTERS:
            raise ValueError(f"{ident}: filter is not in HSC-I family")
        if target_band != control_band:
            exclusions.append({"target_id": ident, "control_id": control_id,
                               "target_filter": target_band, "control_filter": control_band})
            continue
        pair, _ = stack.process_pair(
            target, masks[ident], mask_summary, control, download, control_base,
            control_status_sha, cosmology, kpc_max, n_side, mask_factor,
            include_independent_profile=True)
        pair_rows.append(pair)
        for key in ("target_image", "target_mask", "control_image"):
            path_key = pair[key]
            hash_key = key.replace("_image", "_image_sha256").replace("_mask", "_mask_sha256")
            input_hashes[path_key] = pair[hash_key]
        if index % 20 == 0 or index == len(targets):
            print(f"Independent-support {index}/{len(targets)} targets; pairs={len(pair_rows)}; "
                  f"elapsed={time.monotonic() - start_time:.1f}s", flush=True)

    expected = len(controls_doc["clusters"]) - len(exclusions)
    if len(pair_rows) != expected or expected != common_summary.get("n_pairs_stacked"):
        raise ValueError(f"target/control pair set mismatch: {len(pair_rows)} vs common {expected}")
    common_pairs = {(row["target_id"], row["control_id"], row["filter"])
                    for row in stack.read_tsv(common_work / "paired_radial_profiles.tsv")}
    independent_pairs = {(row["id"], row["control_id"], row["filter"]) for row in pair_rows}
    if independent_pairs != common_pairs:
        raise ValueError("independent/common pair identities differ")

    work.mkdir(parents=True, exist_ok=False)
    pair_path = work / "paired_independent_support_profiles.tsv"
    pair_records = []
    for pair in pair_rows:
        pair_records.append({
            "target_id": pair["id"], "control_id": pair["control_id"], "filter": pair["filter"],
            "z": pair["z"], "ra_deg": pair["ra_deg"], "dec_deg": pair["dec_deg"],
            "n_common_grid_pixels": pair["n_common_pixels"],
            "common_grid_fraction": pair["common_grid_fraction"],
            "profile_json": json.dumps(pair["independent_profile"], separators=(",", ":"), allow_nan=False),
        })
    write_tsv(pair_path, pair_records)
    radial_rows = summarize_pairs(pair_rows, common_profiles)
    radial_path = work / "independent_support_radial_profiles.tsv"
    write_tsv(radial_path, radial_rows)

    repo_files = [Path(__file__).resolve(), Path(stack.__file__).resolve(),
                  HERE / "download_pilot.py", HERE / "validate_control_images.py"]
    result = {
        "schema": "astrafek-hsc-icl-independent-support-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "sample": str(sample), "mask_factor": float(mask_factor),
        "n_targets": len(targets), "n_matched_controls": len(controls_doc["clusters"]),
        "n_filter_mismatch_excluded": len(exclusions), "n_pairs": len(pair_rows),
        "n_by_filter": {band: sum(row["filter"] == band for row in pair_rows) for band in stack.FILTERS},
        "method": "For each matched pair and annulus, compute target and control means using each frame's own valid warped pixels, then subtract the two means. Compare to strict common-support results on identical pairs.",
        "grid": {"coordinate": "proper_kpc", "kpc_max": kpc_max, "n_side": n_side,
                 "pixel_step_kpc": 2.0 * kpc_max / (n_side - 1)},
        "exclusions": exclusions,
        "common_support_reference": {
            "work": str(common_work), "summary_sha256": stack.sha256(common_work / "summary.json"),
            "verification_sha256": stack.sha256(common_work / "independent_verification.json"),
            "pair_profile_sha256": stack.sha256(common_profile_path),
            "n_pairs": common_summary["n_pairs_stacked"],
        },
        "input_sha256": input_hashes,
        "scripts_sha256": {str(path): stack.sha256(path) for path in repo_files},
        "outputs": {"pair_profiles": str(pair_path), "pair_profiles_sha256": stack.sha256(pair_path),
                    "radial_profiles": str(radial_path), "radial_profiles_sha256": stack.sha256(radial_path)},
        "elapsed_seconds": time.monotonic() - start_time,
        "limitations": [
            "target and control annuli use different valid angular pixels; azimuthal symmetry and matched populations are assumed",
            "masking can depend on source environment, so separate support can change the estimator as well as increase area",
            "the output does not correct residual sky transfer, PSF wings, extinction, or correlated coadd noise",
            "pair SEM is descriptive and omits spatial-block and full image-level covariance",
        ],
    }
    summary_path = work / "summary.json"
    stack.write_json(summary_path, result)
    print(json.dumps({"status": result["status"], "mask_factor": mask_factor,
                      "n_pairs": len(pair_rows), "n_by_filter": result["n_by_filter"],
                      "elapsed_seconds": result["elapsed_seconds"],
                      "science_measurement": False}, indent=2), flush=True)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--mask-factor", type=float, choices=stack.MASK_FACTORS, required=True)
    parser.add_argument("--common-work", type=Path)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--kpc-max", type=float, default=500.0)
    parser.add_argument("--n-side", type=int, default=201)
    args = parser.parse_args(argv)
    common_work = args.common_work or DEFAULT_COMMON_WORKS[args.mask_factor]
    work = args.work or DEFAULT_WORKS[args.mask_factor]
    run(args.sample, common_work, work, args.mask_factor, args.kpc_max, args.n_side)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
