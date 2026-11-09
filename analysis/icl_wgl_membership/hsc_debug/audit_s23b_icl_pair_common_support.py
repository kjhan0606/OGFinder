#!/usr/bin/env python3
"""Audit target/control mask overlap without reading science-image pixels.

This is the required pre-stack usable-area review for the approved S23B sample.
Only FITS headers and the prepared mask planes are read. No flux values are
measured or stacked, and the output never promotes the science gate.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import stack_hsc_icl_target_control_2d as stack
import stack_s23b_icl_target_control_2d as s23b


def read_mask(record: dict, summary: dict, factor: float) -> tuple[np.ndarray, Path]:
    path = Path(record["mask_file"])
    expected = summary.get("mask_output_sha256", {}).get(str(path))
    if not expected or s23b.checked_sha256(path) != expected:
        raise ValueError(f"prepared mask checksum mismatch: {path}")
    hdu_name = stack.mask_hdu_name(factor)
    with fits.open(path, memmap=False, checksum=True) as hdul:
        if hdu_name not in hdul:
            raise ValueError(f"prepared mask lacks {hdu_name}: {path}")
        return np.asarray(hdul[hdu_name].data, dtype=bool), path


def radial_support(target_weight: np.ndarray, control_weight: np.ndarray,
                   axis: np.ndarray, edges: np.ndarray) -> list[dict]:
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    target_valid = np.asarray(target_weight) > 0
    control_valid = np.asarray(control_weight) > 0
    common_valid = target_valid & control_valid
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        annulus = (radius >= lo) & (radius < hi)
        geometric = int(annulus.sum())
        nt = int((annulus & target_valid).sum())
        nc = int((annulus & control_valid).sum())
        ncommon = int((annulus & common_valid).sum())
        rows.append({
            "r_min_kpc": float(lo), "r_max_kpc": float(hi),
            "geometric_pixels": geometric,
            "n_target_pixels": nt, "n_control_pixels": nc, "n_common_pixels": ncommon,
            "target_coverage_fraction": nt / geometric if geometric else 0.0,
            "control_coverage_fraction": nc / geometric if geometric else 0.0,
            "common_coverage_fraction": ncommon / geometric if geometric else 0.0,
        })
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(pair_rows: list[dict], edges: np.ndarray) -> list[dict]:
    result = []
    for band in stack.FILTERS:
        pairs = [row for row in pair_rows if row["filter"] == band]
        for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            row = {"filter": band, "r_min_kpc": float(lo), "r_max_kpc": float(hi),
                   "n_pairs_available": len(pairs)}
            for key, stem in (("target_coverage_fraction", "target"),
                              ("control_coverage_fraction", "control"),
                              ("common_coverage_fraction", "common")):
                values = np.asarray([pair["profile"][index][key] for pair in pairs], dtype=float)
                row[f"n_pairs_{stem}_support"] = int(np.sum(values > 0))
                row[f"n_pairs_{stem}_zero_support"] = int(np.sum(values == 0))
                row[f"median_pair_{stem}_coverage"] = float(np.median(values)) if values.size else None
                row[f"p10_pair_{stem}_coverage"] = (
                    float(np.quantile(values, 0.1)) if values.size else None)
            result.append(row)
    return result


def run(args: argparse.Namespace) -> dict:
    if not np.isfinite(args.kpc_max) or args.kpc_max != 1000.0:
        raise ValueError("the approved S23B cutouts require a 1,000-proper-kpc half-width")
    if args.n_side < 101 or args.n_side % 2 != 1:
        raise ValueError("n-side must be odd and at least 101")
    if abs(1000.0 / args.radial_bin_kpc - round(1000.0 / args.radial_bin_kpc)) > 1e-10:
        raise ValueError("radial-bin-kpc must divide the 1,000-kpc range evenly")
    work = args.work.resolve()
    if work.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {work}")

    target_manifest_path = args.target_manifest.resolve()
    control_manifest_path = args.control_manifest.resolve()
    target_manifest = stack.read_json(target_manifest_path)
    target_download, target_status, target_audit, target_download_path, target_status_path, _ = s23b.load_download(
        args.target_download_dir, target_manifest_path, s23b.EXPECTED_TARGET_MANIFEST_SHA256,
        args.target_audit_summary, 1043, "S23B target")
    control_download, control_status, control_audit, control_download_path, control_status_path, control_base = s23b.load_download(
        args.control_download_dir, control_manifest_path, s23b.EXPECTED_CONTROL_MANIFEST_SHA256,
        args.control_audit_summary, 1041, "S23B matched-control")
    target_mask_summary, support_assessment, target_masks, targets_by_id, target_mask_summary_path, target_mask_records_path = s23b.load_target_masks(
        args.target_mask_dir, args.target_support_assessment, target_manifest,
        target_manifest_path, target_download_path, target_status_path)
    controls = control_download["clusters"]
    control_mask_summary, control_masks, control_mask_summary_path, control_mask_records_path = s23b.load_control_masks(
        args.control_mask_dir, control_manifest_path, control_download_path,
        control_status_path, args.control_audit_summary.resolve(), controls,
        control_status["completed"])

    cosmology = stack.Cosmology(H0=float(target_manifest["h0"]), Om=float(target_manifest["omega_m"]))
    controls_by_target = {str(row["target_cluster_id"]): row for row in controls}
    if (len(controls_by_target) != len(controls)
            or not set(controls_by_target).issubset(set(targets_by_id))):
        raise ValueError("S23B controls are not a unique one-to-one mapping onto target IDs")
    filter_exclusions, zero_support_exclusions = [], []
    pair_rows = []
    input_hashes = {str(path): s23b.checked_sha256(path) for path in (
        target_manifest_path, target_download_path, target_status_path,
        args.target_audit_summary.resolve(), target_mask_summary_path,
        target_mask_records_path, args.target_support_assessment.resolve(),
        control_manifest_path, control_download_path, control_status_path,
        args.control_audit_summary.resolve(), control_mask_summary_path,
        control_mask_records_path,
    )}
    input_hashes.update(target_mask_summary.get("source_sha256", {}))
    input_hashes.update(control_mask_summary.get("source_sha256", {}))
    edges = np.arange(0.0, 1000.0 + args.radial_bin_kpc, args.radial_bin_kpc, dtype=float)
    axis = np.linspace(-args.kpc_max, args.kpc_max, args.n_side, dtype=float)
    start = time.monotonic()
    work.mkdir(parents=True, exist_ok=False)

    for index, control in enumerate(controls, 1):
        target_id, control_id = str(control["target_cluster_id"]), str(control["id"])
        target = targets_by_id[target_id]
        target_mask_record = target_masks.get(target_id)
        if target_mask_record is None:
            if target_id not in s23b.EXPECTED_ZERO_SUPPORT_IDS:
                raise ValueError(f"{target_id}: missing target mask is not an audited zero-support target")
            zero_support_exclusions.append({"target_id": target_id, "control_id": control_id,
                                            "reason": "target_frame_zero_unusable_pixels"})
            continue
        control_mask_record = control_masks[control_id]
        target_band = str(target_mask_record["primary_filter"]).strip().lower()
        control_band = str(control.get("coadd_header_filter", "")).strip().lower()
        if target_band != control_band:
            filter_exclusions.append({"target_id": target_id, "control_id": control_id,
                                      "target_filter": target_band,
                                      "control_filter": control_band,
                                      "reason": "target_control_filter_mismatch"})
            continue
        if target_band not in stack.FILTERS:
            raise ValueError(f"{target_id}: unexpected target/control filter {target_band}")

        target_image_path = Path(target_mask_record["image"])
        control_image_path = stack.control_image_path(Path(control_base), control)
        target_info = stack.validate_fits(target_image_path, target)
        control_info = stack.validate_fits(control_image_path, control)
        target_mask, target_mask_path = read_mask(target_mask_record, target_mask_summary, args.mask_factor)
        control_mask, control_mask_path = read_mask(control_mask_record, control_mask_summary, args.mask_factor)
        if target_mask.shape != tuple(target_info["shape"]) or control_mask.shape != tuple(control_info["shape"]):
            raise ValueError(f"{target_id}: prepared mask shape differs from audited image header")
        target_warp = stack.warp_metric(
            np.zeros(target_mask.shape, dtype=np.float32), target_info["center_pixel_xy"],
            target_info["pixel_scale_arcsec"], float(target["z"]), cosmo=cosmology,
            mask=target_mask, kpc_max=args.kpc_max, n_side=args.n_side)
        control_warp = stack.warp_metric(
            np.zeros(control_mask.shape, dtype=np.float32), control_info["center_pixel_xy"],
            control_info["pixel_scale_arcsec"], float(target["z"]), cosmo=cosmology,
            mask=control_mask, kpc_max=args.kpc_max, n_side=args.n_side)
        profile = radial_support(target_warp["weight"], control_warp["weight"], axis, edges)
        pair_rows.append({"target_id": target_id, "control_id": control_id,
                          "filter": target_band, "z": float(target["z"]),
                          "target_mask": str(target_mask_path),
                          "control_mask": str(control_mask_path), "profile": profile})
        input_hashes[str(target_mask_path)] = target_mask_record["mask_sha256"]
        input_hashes[str(control_mask_path)] = control_mask_record["mask_sha256"]
        if index % 20 == 0 or index == len(controls):
            print(f"Support audit {index}/{len(controls)} controls; pairs={len(pair_rows)}; "
                  f"elapsed={time.monotonic() - start:.1f}s", flush=True)

    unmatched = sorted(set(targets_by_id) - set(controls_by_target))
    if len(unmatched) != 2 or len(pair_rows) + len(filter_exclusions) + len(zero_support_exclusions) != len(controls):
        raise ValueError("S23B pair/support accounting does not close")
    zero_common_pairs = [
        {"target_id": pair["target_id"], "control_id": pair["control_id"], "filter": pair["filter"]}
        for pair in pair_rows
        if not any(int(row["n_common_pixels"]) > 0 for row in pair["profile"])
    ]
    pair_path = work / "paired_common_support.tsv"
    flat_rows = []
    for pair in pair_rows:
        for profile in pair["profile"]:
            flat_rows.append({"target_id": pair["target_id"], "control_id": pair["control_id"],
                              "filter": pair["filter"], "z": pair["z"],
                              "target_mask": pair["target_mask"],
                              "control_mask": pair["control_mask"], **profile})
    write_tsv(pair_path, flat_rows)
    summary_rows = summarize(pair_rows, edges)
    summary_path = work / "radial_support_summary.tsv"
    write_tsv(summary_path, summary_rows)
    script_paths = [Path(__file__).resolve(), Path(s23b.__file__).resolve(),
                    Path(stack.__file__).resolve(), HERE / "download_pilot.py",
                    HERE / "prepare_full_cluster_mask_support.py", HERE / "prepare_detection_psf.py",
                    HERE.parents[2] / "icl/wgl.py"]
    result = {
        "schema": "astrafek-s23b-icl-pair-common-support-audit-v2",
        "status": "complete_mask_only_support_audit",
        "n_targets_available": len(targets_by_id), "n_targets_without_control": len(unmatched),
        "targets_without_control": unmatched, "n_controls_available": len(controls),
        "n_pairs_common_filter": len(pair_rows), "n_filter_mismatch_exclusions": len(filter_exclusions),
        "filter_mismatch_exclusions": filter_exclusions,
        "n_pairs_zero_common_support_within_1000_kpc": len(zero_common_pairs),
        "zero_common_support_pairs_within_1000_kpc": zero_common_pairs,
        "n_zero_support_target_exclusions": len(zero_support_exclusions),
        "zero_support_target_exclusions": zero_support_exclusions,
        "n_by_filter": {band: sum(row["filter"] == band for row in pair_rows)
                        for band in stack.FILTERS},
        "mask_factor": float(args.mask_factor),
        "grid": {"coordinate": "proper_kpc", "kpc_max": args.kpc_max,
                 "n_side": args.n_side, "pixel_step_kpc": float(axis[1] - axis[0]),
                 "radial_bin_kpc": args.radial_bin_kpc},
        "products": {"paired_support": str(pair_path),
                     "paired_support_sha256": s23b.checked_sha256(pair_path),
                     "radial_summary": str(summary_path),
                     "radial_summary_sha256": s23b.checked_sha256(summary_path)},
        "input_sha256": input_hashes,
        "script_sha256": {str(path): s23b.checked_sha256(path) for path in script_paths},
        "science_image_pixels_read": False, "mask_pixels_read": True,
        "flux_measured": False, "flux_stacked": False,
        "science_measurement": False,
        "limitations": [
            "This audit reads FITS headers and prepared mask pixels only; it does not inspect surface brightness or residual sky.",
            "Warped support measures interpolation coverage, not completeness or a scientific acceptance threshold.",
            "The four audited zero-support target images are excluded; two other targets have no matched random control.",
            "The `i` and `i2` filter groups are kept separate.",
            "No flux, extinction, PSF, sky-transfer, covariance, or ICL measurement is calculated.",
        ],
        "elapsed_seconds": time.monotonic() - start,
    }
    stack.write_json(work / "summary.json", result)
    print(json.dumps({key: result[key] for key in (
        "status", "n_pairs_common_filter", "n_by_filter", "n_filter_mismatch_exclusions",
        "grid", "science_image_pixels_read", "flux_stacked", "science_measurement")}, indent=2), flush=True)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--target-download-dir", type=Path, required=True)
    parser.add_argument("--target-audit-summary", type=Path, required=True)
    parser.add_argument("--target-mask-dir", type=Path, required=True)
    parser.add_argument("--target-support-assessment", type=Path, required=True)
    parser.add_argument("--control-manifest", type=Path, required=True)
    parser.add_argument("--control-download-dir", type=Path, required=True)
    parser.add_argument("--control-audit-summary", type=Path, required=True)
    parser.add_argument("--control-mask-dir", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--kpc-max", type=float, default=1000.0)
    parser.add_argument("--n-side", type=int, default=401)
    parser.add_argument("--mask-factor", type=float, choices=stack.MASK_FACTORS, default=2.5)
    parser.add_argument("--radial-bin-kpc", type=float, default=25.0)
    args = parser.parse_args(argv)
    if not np.isfinite(args.radial_bin_kpc) or args.radial_bin_kpc <= 0:
        parser.error("radial-bin-kpc must be finite and positive")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
