#!/usr/bin/env python3
"""Independently verify a completed exploratory HSC ICL 2-D stack."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


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


def verify(sample: Path, work: Path) -> dict:
    summary_path = work / "summary.json"
    summary = read_json(summary_path)
    require = lambda test, message: (_ for _ in ()).throw(ValueError(message)) if not test else None
    require(summary.get("schema") == "astrafek-hsc-icl-target-control-2d-stack-v1",
            "unexpected stack summary schema")
    require(summary.get("status") == "complete_exploratory", "stack is not complete")
    require(summary.get("science_measurement") is False, "science gate must remain false")

    sample_manifest = sample / "sample_manifest.json"
    mask_summary = sample / "prestack_cluster_masks_2026-10-16_01/summary.json"
    mask_records_path = sample / "prestack_cluster_masks_2026-10-16_01/records.tsv"
    control_manifest = sample / "control_download_01/download_manifest.json"
    control_status = sample / "control_download_01/download_status.json"
    for path in (sample_manifest, mask_summary, mask_records_path, control_manifest, control_status):
        require(path.is_file(), f"missing provenance input: {path}")
        require(summary["input_sha256"].get(str(path)) == sha256(path),
                f"provenance checksum mismatch: {path}")

    targets = read_json(sample_manifest)["clusters"]
    masks = {str(row["id"]): row for row in read_tsv(mask_records_path)}
    controls = read_json(control_manifest)["clusters"]
    downloads = read_json(control_status)
    require(len(targets) == 800 and len(masks) == 800 and len(controls) == 798,
            "input frame counts differ from the approved sample")
    require(downloads.get("status") == "complete" and downloads.get("downloaded") == 798,
            "matched-control downloads are incomplete")

    target_by_id = {str(row["id"]): row for row in targets}
    eligible = {}
    excluded = []
    for control in controls:
        tid = str(control["target_cluster_id"])
        require(tid in target_by_id and tid in masks, f"unknown matched target: {tid}")
        target_filter = masks[tid]["primary_filter"]
        control_filter = control.get("coadd_header_filter")
        if target_filter != control_filter:
            excluded.append((tid, str(control["id"]), target_filter, control_filter))
        else:
            eligible[tid] = (str(control["id"]), control_filter)
    observed_excluded = {
        (str(row["target_id"]), str(row["control_id"]), row["target_filter"], row["control_filter"])
        for row in summary.get("filter_mismatch_exclusions", [])
    }
    require(observed_excluded == set(excluded), "filter-mismatch exclusions differ from input metadata")
    require(summary.get("n_pairs_excluded_filter_mismatch") == len(excluded),
            "filter-mismatch count is wrong")
    require(summary.get("n_pairs_stacked") == len(eligible), "eligible pair count is wrong")

    pair_path = work / "paired_radial_profiles.tsv"
    radial_path = work / "stacked_radial_profiles.tsv"
    require(summary["outputs"]["paired_profiles_sha256"] == sha256(pair_path),
            "paired radial-profile checksum mismatch")
    require(summary["outputs"]["stacked_profiles_sha256"] == sha256(radial_path),
            "stacked radial-profile checksum mismatch")
    pair_rows = read_tsv(pair_path)
    by_target = {row["target_id"]: row for row in pair_rows}
    require(len(by_target) == len(pair_rows) == len(eligible), "paired profile rows are duplicated or incomplete")
    for tid, (control_id, band) in eligible.items():
        row = by_target.get(tid)
        require(row is not None and row["control_id"] == control_id and row["filter"] == band,
                f"pair table mismatch for target {tid}")
        profile = json.loads(row["profile_json"])
        require(len(profile) == 20, f"target {tid} does not have 20 radial bins")
        for annulus in profile:
            if annulus["delta_mean_njy_arcsec2"] is not None:
                require(np.isclose(annulus["target_mean_njy_arcsec2"]
                                   - annulus["control_mean_njy_arcsec2"],
                                   annulus["delta_mean_njy_arcsec2"], rtol=1e-6, atol=1e-6),
                        f"target-control radial subtraction mismatch for {tid}")
    radial_rows = read_tsv(radial_path)
    require(len(radial_rows) == 40, "stacked radial profile must have 20 bins per filter")

    map_checks = {}
    for band in ("i", "i2"):
        product = summary["outputs"]["maps"][band]
        path = Path(product["file"])
        require(path.is_file() and sha256(path) == product["sha256"],
                f"{band} FITS output checksum mismatch")
        with fits.open(path, memmap=False, checksum=True) as hdus:
            required = {"TARGET", "CONTROL", "DELTA", "NPAIR", "TARGET_SEM",
                        "CONTROL_SEM", "DELTA_SEM", "GRID_KPC"}
            require(required.issubset({h.name for h in hdus}), f"{band} FITS extensions are incomplete")
            require(hdus[0].header.get("SCIVAL") is False, f"{band} FITS science gate is not false")
            require(hdus[0].header.get("FILTER") == band, f"{band} FITS filter header mismatch")
            target, control = np.asarray(hdus["TARGET"].data), np.asarray(hdus["CONTROL"].data)
            delta, count = np.asarray(hdus["DELTA"].data), np.asarray(hdus["NPAIR"].data)
            require(target.shape == control.shape == delta.shape == count.shape == (201, 201),
                    f"{band} FITS image shapes are unexpected")
            require(int(count.max()) <= int(summary["n_by_filter"][band]),
                    f"{band} pixel pair count exceeds number of pairs")
            common = np.isfinite(target) & np.isfinite(control) & np.isfinite(delta)
            require(np.array_equal(common, count > 0), f"{band} map coverage and pair-count masks differ")
            residual = delta[common] - (target[common] - control[common])
            max_residual = float(np.max(np.abs(residual))) if residual.size else 0.0
            scale = max(1.0, float(np.max(np.abs(delta[common]))) if residual.size else 1.0)
            require(max_residual <= 2e-5 * scale,
                    f"{band} mean maps do not satisfy target minus control")
            require(int(common.sum()) == product["covered_pixels"], f"{band} covered-pixel count mismatch")
            map_checks[band] = {"n_pairs": int(summary["n_by_filter"][band]),
                                "covered_pixels": int(common.sum()),
                                "max_delta_identity_residual": max_residual}

    for name, expected_hash in summary.get("script_sha256", {}).items():
        path = Path(name)
        current_matches = path.is_file() and sha256(path) == expected_hash
        snapshot_matches = (
            name == str((Path(__file__).resolve().parent / "stack_hsc_icl_target_control_2d.py").resolve())
            and summary.get("executed_script_snapshot_sha256") == expected_hash
            and Path(summary.get("executed_script_snapshot", "")).is_file()
            and sha256(Path(summary["executed_script_snapshot"])) == expected_hash
        )
        require(current_matches or snapshot_matches,
                f"analysis script changed since run and no matching run snapshot exists: {path}")
    if "serialization_correction" in summary:
        correction = summary["serialization_correction"]
        require(summary.get("science_measurement") is False
                and correction.get("data_products_modified") is False
                and correction.get("executed_script_snapshot_sha256")
                == summary.get("executed_script_snapshot_sha256"),
                "summary serialization correction record is invalid")
        repair_script = Path(correction.get("repair_script", ""))
        require(repair_script.is_file()
                and sha256(repair_script) == correction.get("repair_script_sha256"),
                "summary serialization repair provenance is invalid")
    return {
        "schema": "astrafek-hsc-icl-target-control-2d-stack-verification-v1",
        "status": "verified",
        "summary_sha256": sha256(summary_path),
        "n_input_controls": len(controls), "n_filter_mismatch_excluded": len(excluded),
        "n_eligible_pairs": len(eligible), "n_pair_rows": len(pair_rows),
        "n_stacked_radial_rows": len(radial_rows), "map_checks": map_checks,
        "science_measurement": False,
        "input_images_rehashed": False,
        "note": "The processing script rehashed every input frame; this verifier independently rechecks manifests, outputs, maps, pairing and calculations without another full-image read.",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = verify(args.sample.resolve(), args.work.resolve())
    output = args.output.resolve() if args.output else args.work.resolve() / "independent_verification.json"
    if output.exists():
        parser.error(f"refusing to overwrite existing verification: {output}")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"verification_sha256={sha256(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
