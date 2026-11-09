#!/usr/bin/env python3
"""Build a provenance-pinned exploratory 2-D ICL stack for the S23B sample.

The original matched-stack script is intentionally fixed to the legacy
800-target/798-control sample. This entry point accepts the separate S23B
target/control products, verifies their exact manifests and audits, excludes
only the four target frames documented by the zero-support audit, and reports
HSC ``i`` and ``i2`` independently. It is an exploratory engineering product:
it does not set ``science_measurement`` true.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import stack_hsc_icl_target_control_2d as stack


EXPECTED_TARGET_MANIFEST_SHA256 = (
    "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
)
EXPECTED_CONTROL_MANIFEST_SHA256 = (
    "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df"
)
EXPECTED_ZERO_SUPPORT_IDS = {
    "S23B_10647", "S23B_11786", "S23B_432", "S23B_8226",
}
PROFILE_EDGES_KPC = np.arange(0.0, 1000.0 + 25.0, 25.0)


def checked_sha256(path: Path) -> str:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    return stack.sha256(path)


def check_audit_inputs(audit: dict, paths: list[Path], label: str) -> None:
    recorded = audit.get("inputs_sha256", {})
    for recorded_path, digest in recorded.items():
        path = Path(recorded_path)
        if not path.is_absolute():
            raise ValueError(f"{label} audit has a non-absolute input path: {recorded_path}")
        if checked_sha256(path) != digest:
            raise ValueError(f"{label} audit recorded input changed: {recorded_path}")
    for path in paths:
        absolute = str(path.resolve())
        expected = recorded.get(absolute)
        if expected is None:
            raise ValueError(f"{label} audit input hash mismatch: {absolute}")


def check_audit_outputs(audit: dict, audit_path: Path, label: str) -> None:
    for output, digest in audit.get("output_sha256", {}).items():
        path = Path(output)
        if not path.is_absolute():
            path = audit_path.resolve().parent / path
        if checked_sha256(path) != digest:
            raise ValueError(f"{label} audit output changed: {path}")


def load_download(download_dir: Path, approved_manifest: Path,
                  expected_manifest_sha256: str, audit_path: Path,
                  expected_n: int, label: str) -> tuple[dict, dict, dict, Path, Path, str]:
    download_dir = download_dir.resolve()
    approved_manifest = approved_manifest.resolve()
    audit_path = audit_path.resolve()
    manifest_path = download_dir / "download_manifest.json"
    status_path = download_dir / "download_status.json"
    for path in (manifest_path, status_path, audit_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    approved_sha = checked_sha256(approved_manifest)
    if approved_sha != expected_manifest_sha256:
        raise ValueError(f"{label} approved-manifest SHA-256 changed")
    manifest = stack.read_json(manifest_path)
    status = stack.read_json(status_path)
    audit = stack.read_json(audit_path)
    frames = manifest.get("clusters", [])
    completed = status.get("completed", [])
    expected = stack.read_json(approved_manifest).get("clusters", [])
    by_id = {str(row.get("id")): row for row in frames}
    completed_by_id = {str(row.get("id")): row for row in completed}
    approved_by_id = {str(row.get("id")): row for row in expected}
    if (manifest.get("science_measurement") is not False
            or status.get("status") != "complete"
            or status.get("manifest_sha256") != approved_sha
            or status.get("n_expected") != expected_n
            or status.get("science_measurement") is not False
            or status.get("credentials_saved") is not False
            or len(frames) != expected_n or len(completed) != expected_n
            or len(by_id) != expected_n or len(completed_by_id) != expected_n
            or len(approved_by_id) != expected_n
            or set(by_id) != set(completed_by_id)
            or set(by_id) != set(approved_by_id)):
        raise ValueError(f"{label} download is incomplete or outside the pinned manifest")
    audit_inputs = audit.get("inputs_sha256", {})
    audited_manifest_sha = audit.get("manifest_sha256") or audit_inputs.get(str(approved_manifest))
    if (audit.get("status") != "complete_header_wcs_file_integrity_audit"
            or audit.get("n_expected_frames") != expected_n
            or audit.get("n_verified_frames") != expected_n
            or audited_manifest_sha != approved_sha
            or audit.get("science_measurement") is not False
            or audit.get("images_ready_for_stacking") is not False):
        raise ValueError(f"{label} integrity/WCS audit is incomplete or unexpected")
    check_audit_inputs(audit, [approved_manifest, manifest_path, status_path], label)
    check_audit_outputs(audit, audit_path, label)

    for ident, frame in by_id.items():
        source = approved_by_id[ident]
        for key in ("ra_deg", "dec_deg", "z"):
            if float(frame[key]) != float(source[key]):
                raise ValueError(f"{label} {ident}: downloaded {key} differs from approved plan")
        record = completed_by_id[ident]
        if (not record.get("file") or Path(str(record["file"])).name != record["file"]
                or int(record.get("bytes", -1)) <= 0
                or not record.get("sha256")):
            raise ValueError(f"{label} {ident}: malformed completed-file record")

    download_base = Path(manifest.get("downloaded_images", "images"))
    if not download_base.is_absolute():
        download_base = download_dir / download_base
    return manifest, status, audit, manifest_path, status_path, str(download_base.resolve())


def load_target_masks(mask_dir: Path, assessment_path: Path,
                      target_manifest: dict, target_manifest_path: Path,
                      target_download_manifest_path: Path,
                      target_status_path: Path) -> tuple[dict, dict, dict, dict, Path, Path]:
    mask_dir = mask_dir.resolve()
    summary_path = mask_dir / "summary.json"
    records_path = mask_dir / "records.tsv"
    assessment_path = assessment_path.resolve()
    for path in (summary_path, records_path, assessment_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    summary = stack.read_json(summary_path)
    assessment = stack.read_json(assessment_path)
    records = {str(row["id"]): row for row in stack.read_tsv(records_path)}
    target_rows = target_manifest.get("clusters", [])
    target_ids = {str(row["id"]) for row in target_rows}
    zero_ids = set(map(str, assessment.get("zero_pixel_exclusion_ids", [])))
    errors = {str(row.get("id")) for row in summary.get("errors", [])}

    if (assessment.get("status") != "support_assessment_complete_science_gate_closed"
            or assessment.get("science_measurement") is not False
            or assessment.get("images_ready_for_stacking") is not False
            or assessment.get("n_expected_targets") != len(target_ids)
            or assessment.get("n_mask_support_success") != len(records)
            or assessment.get("n_zero_pixel_exclusions") != len(zero_ids)
            or zero_ids != EXPECTED_ZERO_SUPPORT_IDS):
        raise ValueError("S23B target-support assessment does not match the pinned four-frame exclusion")
    if (summary.get("n_expected") != len(target_ids)
            or summary.get("n_success") != len(records)
            or summary.get("science_measurement") is not False
            or summary.get("mask_ready_for_science") is not False
            or errors != zero_ids
            or set(records) != target_ids - zero_ids):
        raise ValueError("target mask records/errors do not match the independent zero-support audit")
    if (len(records) + len(errors) != len(target_ids)
            or summary.get("status") != "partial"):
        raise ValueError("target mask preparation has unaccounted records or errors")

    source_hashes = summary.get("source_sha256", {})
    for source_path, digest in source_hashes.items():
        if checked_sha256(Path(source_path)) != digest:
            raise ValueError(f"target mask preparation input changed: {source_path}")
    for path in (target_manifest_path, target_download_manifest_path, target_status_path):
        if str(path.resolve()) not in source_hashes:
            raise ValueError(f"target mask preparation omitted required input: {path}")
    # The summary cannot contain its own digest. The independent support
    # assessment records the target summary hash, while this summary pins the
    # records table and generated mask FITS files.
    for path in (records_path,):
        absolute = str(path.resolve())
        if summary.get("output_sha256", {}).get(absolute) != checked_sha256(path):
            raise ValueError(f"target-mask output checksum mismatch: {absolute}")
    for path, digest in summary.get("mask_output_sha256", {}).items():
        if checked_sha256(Path(path)) != digest:
            raise ValueError(f"target-mask FITS checksum mismatch: {path}")
    assessment_inputs = assessment.get("inputs_sha256", {})
    for source_path, digest in assessment_inputs.items():
        if checked_sha256(Path(source_path)) != digest:
            raise ValueError(f"support-assessment input changed: {source_path}")
    if str(summary_path.resolve()) not in assessment_inputs:
        raise ValueError("support assessment is not tied to this target-mask preparation")
    for ident, row in records.items():
        if (row.get("science_measurement") != "False"
                or row.get("mask_ready_for_science") != "False"
                or row.get("mask_sha256") != summary.get("mask_output_sha256", {}).get(row.get("mask_file"))):
            raise ValueError(f"{ident}: target mask record is not provenance-valid")
    return summary, assessment, records, {str(row["id"]): row for row in target_rows}, summary_path, records_path


def load_control_masks(mask_dir: Path, control_manifest_path: Path,
                       control_download_manifest_path: Path,
                       control_status_path: Path, control_audit_path: Path,
                       controls: list[dict], completed: list[dict]) -> tuple[dict, dict, Path, Path]:
    mask_dir = mask_dir.resolve()
    summary_path, records_path = mask_dir / "summary.json", mask_dir / "records.tsv"
    if not summary_path.is_file() or not records_path.is_file():
        raise FileNotFoundError("completed S23B control-mask summary/records are required")
    summary = stack.read_json(summary_path)
    records = {str(row["id"]): row for row in stack.read_tsv(records_path)}
    control_by_id = {str(row["id"]): row for row in controls}
    completed_by_id = {str(row["id"]): row for row in completed}
    if (summary.get("status") != "complete" or summary.get("n_expected") != len(controls)
            or summary.get("n_success") != len(controls) or summary.get("n_error") != 0
            or summary.get("science_measurement") is not False
            or summary.get("mask_ready_for_science") is not False
            or len(records) != len(controls) or set(records) != set(control_by_id)
            or len(completed_by_id) != len(controls)):
        raise ValueError("S23B control-mask preparation is incomplete or out of scope")
    for source_path, digest in summary.get("source_sha256", {}).items():
        if checked_sha256(Path(source_path)) != digest:
            raise ValueError(f"control-mask preparation input changed: {source_path}")
    for path in (control_manifest_path, control_download_manifest_path,
                 control_status_path, control_audit_path):
        if str(path.resolve()) not in summary.get("source_sha256", {}):
            raise ValueError(f"control-mask preparation omitted required input: {path}")
    # As with the target preparation, the control summary cannot contain its
    # own digest; the pair audit records the summary hash as an input.
    for path in (records_path,):
        if summary.get("output_sha256", {}).get(str(path.resolve())) != checked_sha256(path):
            raise ValueError(f"control-mask output checksum mismatch: {path}")
    output_hashes = summary.get("mask_output_sha256", {})
    if len(output_hashes) != len(controls):
        raise ValueError("control-mask FITS count differs from the approved control sample")
    for mask_path, digest in output_hashes.items():
        if checked_sha256(Path(mask_path)) != digest:
            raise ValueError(f"control-mask FITS checksum mismatch: {mask_path}")
    for ident, record in records.items():
        frame, done = control_by_id[ident], completed_by_id[ident]
        if (record.get("science_measurement") != "False"
                or record.get("mask_ready_for_science") != "False"
                or record.get("image_sha256") != done.get("sha256")
                or record.get("primary_filter") != frame.get("coadd_header_filter")
                or record.get("mask_sha256") != output_hashes.get(record.get("mask_file"))):
            raise ValueError(f"{ident}: control-mask record is not tied to the audited download")
    return summary, records, summary_path, records_path


def read_prepared_mask(record: dict, summary: dict, mask_factor: float) -> tuple[np.ndarray, Path]:
    path = Path(record["mask_file"])
    expected = summary.get("mask_output_sha256", {}).get(str(path))
    if not expected or checked_sha256(path) != expected:
        raise ValueError(f"prepared control mask hash mismatch: {path}")
    hdu_name = stack.mask_hdu_name(mask_factor)
    with fits.open(path, memmap=False, checksum=True) as hdul:
        if hdu_name not in hdul:
            raise ValueError(f"prepared control mask lacks {hdu_name}: {path}")
        mask = np.asarray(hdul[hdu_name].data, dtype=bool)
    return mask, path


def load_pair_support(summary_path: Path, core_inputs: dict[str, str],
                      mask_factor: float, kpc_max: float, n_side: int) -> tuple[dict, dict, Path]:
    summary_path = summary_path.resolve()
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = stack.read_json(summary_path)
    grid = summary.get("grid", {})
    products = summary.get("products", {})
    support_path = Path(str(products.get("paired_support", "")))
    radial_path = Path(str(products.get("radial_summary", "")))
    if (summary.get("status") != "complete_mask_only_support_audit"
            or summary.get("science_image_pixels_read") is not False
            or summary.get("flux_measured") is not False
            or summary.get("flux_stacked") is not False
            or summary.get("science_measurement") is not False
            or float(summary.get("mask_factor", -1)) != mask_factor
            or float(grid.get("kpc_max", -1)) != kpc_max
            or int(grid.get("n_side", -1)) != n_side
            or float(grid.get("radial_bin_kpc", -1)) != 25.0
            or not support_path.is_absolute() or not radial_path.is_absolute()
            or checked_sha256(support_path) != products.get("paired_support_sha256")
            or checked_sha256(radial_path) != products.get("radial_summary_sha256")):
        raise ValueError("pre-stack common-support audit is missing, incomplete, or uses another grid/mask")
    recorded_inputs = summary.get("input_sha256", {})
    for path, digest in core_inputs.items():
        if recorded_inputs.get(path) != digest:
            raise ValueError(f"common-support audit used a different source: {path}")
    pair_rows = stack.read_tsv(support_path)
    grouped: dict[tuple[str, str, str], list[dict]] = {}
    for row in pair_rows:
        key = (str(row["target_id"]), str(row["control_id"]), str(row["filter"]))
        grouped.setdefault(key, []).append(row)
    expected_bins = int(round(kpc_max / float(grid["radial_bin_kpc"])))
    support = {}
    for key, rows in grouped.items():
        if len(rows) != expected_bins:
            raise ValueError(f"common-support audit has incomplete radial bins for {key}")
        if (len({str(row["target_mask"]) for row in rows}) != 1
                or len({str(row["control_mask"]) for row in rows}) != 1):
            raise ValueError(f"common-support audit changes mask inputs across radial bins for {key}")
        support[key] = {
            "n_common_pixels_within_kpc_max": sum(int(row["n_common_pixels"]) for row in rows),
            "target_mask": str(rows[0]["target_mask"]),
            "control_mask": str(rows[0]["control_mask"]),
        }
        for row in rows:
            if (recorded_inputs.get(support[key]["target_mask"]) is None
                    or recorded_inputs.get(support[key]["control_mask"]) is None):
                raise ValueError(f"common-support audit did not pin mask inputs for {key}")
    by_filter = {band: sum(key[2] == band for key in support) for band in stack.FILTERS}
    zero_common = sum(row["n_common_pixels_within_kpc_max"] == 0 for row in support.values())
    if (len(support) != summary.get("n_pairs_common_filter")
            or by_filter != summary.get("n_by_filter")
            or zero_common != summary.get("n_pairs_zero_common_support_within_1000_kpc")):
        raise ValueError("common-support audit pair counts do not close")
    return summary, support, support_path


def profile_on_grid(target: np.ndarray, control: np.ndarray, delta: np.ndarray,
                    target_weight: np.ndarray, control_weight: np.ndarray,
                    axis_kpc: np.ndarray) -> list[dict]:
    yy, xx = np.meshgrid(axis_kpc, axis_kpc, indexing="ij")
    radius = np.hypot(xx, yy)
    target_valid = (np.asarray(target_weight) > 0) & np.isfinite(target)
    control_valid = (np.asarray(control_weight) > 0) & np.isfinite(control)
    common_valid = target_valid & control_valid & np.isfinite(delta)
    rows = []
    for lo, hi in zip(PROFILE_EDGES_KPC[:-1], PROFILE_EDGES_KPC[1:]):
        annulus = (radius >= lo) & (radius < hi)
        target_selected = annulus & target_valid
        control_selected = annulus & control_valid
        selected = annulus & common_valid
        n_target, n_control = int(target_selected.sum()), int(control_selected.sum())
        n_valid, n_geom = int(selected.sum()), int(annulus.sum())
        rows.append({
            "r_min_kpc": float(lo), "r_max_kpc": float(hi),
            "n_target_pixels": n_target, "n_control_pixels": n_control,
            "n_valid_pixels": n_valid, "geometric_pixels": n_geom,
            "coverage_fraction": n_valid / n_geom if n_geom else 0.0,
            "target_coverage_fraction": n_target / n_geom if n_geom else 0.0,
            "control_coverage_fraction": n_control / n_geom if n_geom else 0.0,
            "target_mean_njy_arcsec2": float(np.mean(target[selected])) if n_valid else None,
            "control_mean_njy_arcsec2": float(np.mean(control[selected])) if n_valid else None,
            "delta_mean_njy_arcsec2": float(np.mean(delta[selected])) if n_valid else None,
            "target_independent_mean_njy_arcsec2": (
                float(np.mean(target[target_selected])) if n_target else None),
            "control_independent_mean_njy_arcsec2": (
                float(np.mean(control[control_selected])) if n_control else None),
        })
    return rows


def summarize_profiles(pair_rows: list[dict]) -> list[dict]:
    output = []
    for band in stack.FILTERS:
        selected = [row for row in pair_rows if row["filter"] == band]
        for index, (lo, hi) in enumerate(zip(PROFILE_EDGES_KPC[:-1], PROFILE_EDGES_KPC[1:])):
            record = {"filter": band, "r_min_kpc": float(lo), "r_max_kpc": float(hi)}
            for key, stem in (("target_mean_njy_arcsec2", "target"),
                              ("control_mean_njy_arcsec2", "control"),
                              ("delta_mean_njy_arcsec2", "delta")):
                values = np.asarray([row["profile"][index][key] for row in selected
                                     if row["profile"][index][key] is not None], dtype=float)
                record[f"n_pairs_{stem}"] = int(values.size)
                record[f"mean_{stem}_njy_arcsec2"] = float(values.mean()) if values.size else None
                record[f"pair_sd_{stem}_njy_arcsec2"] = (
                    float(values.std(ddof=1)) if values.size > 1 else None)
                record[f"pair_sem_{stem}_njy_arcsec2"] = (
                    float(values.std(ddof=1) / math.sqrt(values.size)) if values.size > 1 else None)
            for key, stem in (("target_independent_mean_njy_arcsec2", "target_independent"),
                              ("control_independent_mean_njy_arcsec2", "control_independent")):
                values = np.asarray([row["profile"][index][key] for row in selected
                                     if row["profile"][index][key] is not None], dtype=float)
                record[f"n_pairs_{stem}"] = int(values.size)
                record[f"mean_{stem}_njy_arcsec2"] = float(values.mean()) if values.size else None
            for key, stem in (("target_coverage_fraction", "target"),
                              ("control_coverage_fraction", "control"),
                              ("coverage_fraction", "common")):
                coverage = np.asarray([row["profile"][index][key] for row in selected], dtype=float)
                record[f"median_pair_{stem}_coverage"] = float(np.median(coverage)) if coverage.size else None
                record[f"p10_pair_{stem}_coverage"] = (
                    float(np.quantile(coverage, 0.1)) if coverage.size else None)
            output.append(record)
    return output


def write_profile_table(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty radial-profile table")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict:
    if args.mask_factor not in stack.MASK_FACTORS:
        raise ValueError(f"mask factor must be one of {stack.MASK_FACTORS}")
    if args.kpc_max != 1000.0 or args.n_side < 101 or args.n_side % 2 != 1:
        raise ValueError("S23B cutouts require kpc-max=1000 and an odd grid of at least 101 pixels")
    if not np.isfinite(args.kpc_max):
        raise ValueError("grid extent must be finite")
    work = args.work.resolve()
    if work.exists():
        raise FileExistsError(f"refusing to overwrite existing work directory: {work}")

    target_manifest_path = args.target_manifest.resolve()
    control_plan_path = args.control_manifest.resolve()
    target_download_dir = args.target_download_dir.resolve()
    control_download_dir = args.control_download_dir.resolve()
    target_manifest = stack.read_json(target_manifest_path)
    target_sha = checked_sha256(target_manifest_path)
    if (target_sha != EXPECTED_TARGET_MANIFEST_SHA256
            or target_manifest.get("science_measurement") is not False):
        raise ValueError("target manifest is not the pinned S23B image sample")
    control_sha = checked_sha256(control_plan_path)
    if control_sha != EXPECTED_CONTROL_MANIFEST_SHA256:
        raise ValueError("control manifest is not the user-approved v2 scope")

    target_manifest_doc, target_status, target_audit, target_download_manifest_path, target_status_path, target_base = load_download(
        target_download_dir, target_manifest_path, EXPECTED_TARGET_MANIFEST_SHA256,
        args.target_audit_summary, 1043, "S23B target")
    control_manifest_doc, control_status, control_audit, control_download_manifest_path, control_status_path, control_base_string = load_download(
        control_download_dir, control_plan_path, EXPECTED_CONTROL_MANIFEST_SHA256,
        args.control_audit_summary, 1041, "S23B matched-control")
    control_base = Path(control_base_string)
    control_mask_summary, control_mask_records, control_mask_summary_path, control_mask_records_path = load_control_masks(
        args.control_mask_dir, control_plan_path, control_download_manifest_path,
        control_status_path, args.control_audit_summary.resolve(),
        control_manifest_doc["clusters"], control_status["completed"])
    if (float(target_manifest.get("h0")) != float(control_manifest_doc.get("h0"))
            or float(target_manifest.get("omega_m")) != float(control_manifest_doc.get("omega_m"))):
        raise ValueError("target and control cosmologies differ")

    target_mask_summary, support_assessment, mask_records, targets_by_id, mask_summary_path, mask_records_path = load_target_masks(
        args.target_mask_dir, args.target_support_assessment, target_manifest,
        target_manifest_path, target_download_manifest_path, target_status_path)
    if len(targets_by_id) != 1043 or len(mask_records) != 1039:
        raise ValueError("unexpected S23B target/mask count")

    controls = control_manifest_doc.get("clusters", [])
    target_ids = set(targets_by_id)
    controls_by_target = {str(row["target_cluster_id"]): row for row in controls}
    controls_downloads = {str(row["id"]): row for row in control_manifest_doc.get("clusters", [])}
    control_download_records = {str(row["id"]): row for row in control_status.get("completed", [])}
    if (len(controls) != 1041 or len(controls_by_target) != 1041
            or len(controls_downloads) != 1041 or len(control_download_records) != 1041
            or not set(controls_by_target).issubset(target_ids)):
        raise ValueError("S23B control-to-target mapping is incomplete or duplicated")

    cosmology = stack.Cosmology(H0=float(target_manifest["h0"]), Om=float(target_manifest["omega_m"]))
    target_filters = {ident: row for ident, row in mask_records.items()}
    exclusions: list[dict] = []
    pair_rows: list[dict] = []
    accumulators: dict[str, dict] = {}
    counts = {band: 0 for band in stack.FILTERS}
    input_hashes = {
        str(path): checked_sha256(path) for path in (
            target_manifest_path, target_download_manifest_path, target_status_path,
            args.target_audit_summary.resolve(), mask_summary_path, mask_records_path,
            args.target_support_assessment.resolve(), control_plan_path,
            control_download_manifest_path, control_status_path,
            args.control_audit_summary.resolve(),
            control_mask_summary_path, control_mask_records_path,
        )
    }
    input_hashes.update(target_mask_summary.get("source_sha256", {}))
    input_hashes.update(control_mask_summary.get("source_sha256", {}))
    support_audit_path = args.pair_support_audit_summary.resolve()
    pair_support_summary, pair_support, pair_support_table_path = load_pair_support(
        support_audit_path, input_hashes, args.mask_factor, args.kpc_max, args.n_side)
    input_hashes[str(support_audit_path)] = checked_sha256(support_audit_path)
    input_hashes[str(pair_support_table_path)] = checked_sha256(pair_support_table_path)
    control_status_sha = checked_sha256(control_status_path)
    start = time.monotonic()
    work.mkdir(parents=True, exist_ok=False)

    for index, control in enumerate(controls, 1):
        target_id = str(control["target_cluster_id"])
        target = targets_by_id[target_id]
        mask_record = target_filters.get(target_id)
        control_id = str(control["id"])
        if mask_record is None:
            if target_id not in EXPECTED_ZERO_SUPPORT_IDS:
                raise ValueError(f"{target_id}: missing target mask is not in the pinned zero-support audit")
            exclusions.append({"target_id": target_id, "control_id": control_id,
                               "reason": "target_frame_zero_unusable_pixels"})
            continue

        download_record = control_download_records[control_id]
        control_mask_record = control_mask_records[control_id]
        if control.get("image") != controls_downloads[control_id].get("image"):
            raise ValueError(f"{control_id}: control download product identity mismatch")
        target_band = str(mask_record["primary_filter"]).strip().lower()
        control_band = str(control.get("coadd_header_filter", "")).strip().lower()
        if control_band != str(download_record.get("review", {}).get("coadd_header_filter", "")).strip().lower():
            raise ValueError(f"{control_id}: control header-filter evidence differs from its download audit")
        if target_band not in stack.FILTERS or control_band not in stack.FILTERS:
            raise ValueError(f"{target_id}: unexpected HSC-I-family filter evidence")
        if target_band != control_band:
            exclusions.append({"target_id": target_id, "control_id": control_id,
                               "reason": "target_control_filter_mismatch",
                               "target_filter": target_band, "control_filter": control_band})
            continue

        support_key = (target_id, control_id, target_band)
        support_row = pair_support.get(support_key)
        if support_row is None:
            raise ValueError(f"{support_key}: pair is absent from the pre-stack common-support audit")
        target_mask_path_for_support = str(mask_record["mask_file"])
        control_mask_path_for_support = str(control_mask_records[control_id]["mask_file"])
        if (support_row["target_mask"] != target_mask_path_for_support
                or support_row["control_mask"] != control_mask_path_for_support
                or pair_support_summary.get("input_sha256", {}).get(target_mask_path_for_support)
                != mask_record["mask_sha256"]
                or pair_support_summary.get("input_sha256", {}).get(control_mask_path_for_support)
                != control_mask_records[control_id]["mask_sha256"]):
            raise ValueError(f"{support_key}: support audit used different prepared masks")
        input_hashes[target_mask_path_for_support] = mask_record["mask_sha256"]
        input_hashes[control_mask_path_for_support] = control_mask_records[control_id]["mask_sha256"]
        if int(support_row["n_common_pixels_within_kpc_max"]) == 0:
            exclusions.append({"target_id": target_id, "control_id": control_id,
                               "reason": "zero_common_mask_support_within_1000_kpc"})
            continue

        control_mask, control_mask_path = read_prepared_mask(
            control_mask_record, control_mask_summary, args.mask_factor)
        pair, maps = stack.process_pair(
            target, mask_record, target_mask_summary, control, download_record,
            control_base, control_status_sha, cosmology,
            args.kpc_max, args.n_side, args.mask_factor,
            control_mask_override=control_mask)
        pair["control_mask"] = str(control_mask_path)
        pair["control_mask_sha256"] = control_mask_record["mask_sha256"]
        pair["profile"] = profile_on_grid(
            maps["target"], maps["control"], maps["delta"], maps["target_weight"],
            maps["control_weight"], np.asarray(maps["axis_kpc"], dtype=float))
        band = pair["filter"]
        if band not in accumulators:
            shape = maps["delta"].shape
            accumulators[band] = {key: stack.OnlineMap(shape) for key in ("target", "control", "delta")}
        for key in ("target", "control", "delta"):
            accumulators[band][key].add(maps[key], maps["valid"])
        pair_rows.append(pair)
        counts[band] += 1
        input_hashes[pair["target_image"]] = pair["target_image_sha256"]
        input_hashes[pair["target_mask"]] = pair["target_mask_sha256"]
        input_hashes[pair["control_image"]] = pair["control_image_sha256"]
        input_hashes[pair["control_mask"]] = pair["control_mask_sha256"]

        if index % 10 == 0 or index == len(controls):
            stack.write_json(work / "progress.json", {
                "status": "running", "controls_seen": index,
                "pairs_stacked": len(pair_rows), "by_filter": counts,
                "exclusions": len(exclusions), "elapsed_seconds": time.monotonic() - start,
                "science_measurement": False,
            })
            print(f"S23B stack {index}/{len(controls)} controls; pairs={len(pair_rows)}; "
                  f"bands={counts}; excluded={len(exclusions)}", flush=True)

    unmatched_targets = sorted(target_ids - set(controls_by_target))
    if len(unmatched_targets) != 2:
        raise ValueError(f"expected two target clusters without controls, got {len(unmatched_targets)}")
    expected_support_pairs = {key for key, row in pair_support.items()
                              if int(row["n_common_pixels_within_kpc_max"]) > 0}
    stacked_pair_keys = {(row["id"], row["control_id"], row["filter"]) for row in pair_rows}
    if stacked_pair_keys != expected_support_pairs:
        raise ValueError("flux-stack pair identities differ from pre-stack mask-support eligibility")
    if len(pair_rows) + len(exclusions) != len(controls):
        raise ValueError("control-to-pair accounting does not close")
    if not pair_rows or any(counts[band] == 0 for band in stack.FILTERS):
        raise ValueError(f"empty stack or filter group: {counts}")

    kpc_axis = np.linspace(-args.kpc_max, args.kpc_max, args.n_side, dtype=float)
    output_maps = {}
    for band in stack.FILTERS:
        band_pairs = [row for row in pair_rows if row["filter"] == band]
        output_maps[band] = stack.save_band_fits(
            work / f"s23b_icl_target_control_stack_{band}.fits", band, kpc_axis,
            {"n_pairs": len(band_pairs)}, accumulators[band], len(band_pairs),
            {"H0": float(target_manifest["h0"]), "Om0": float(target_manifest["omega_m"])},
            args.kpc_max, args.mask_factor)

    pair_path = work / "paired_radial_profiles.tsv"
    stack.write_pair_table(pair_path, pair_rows)
    profile_path = work / "stacked_radial_profiles.tsv"
    write_profile_table(profile_path, summarize_profiles(pair_rows))
    scripts = [Path(__file__).resolve(), Path(stack.__file__).resolve(),
               HERE / "download_pilot.py", HERE / "validate_control_images.py",
               HERE / "prepare_detection_psf.py", HERE / "prepare_full_cluster_mask_support.py",
               HERE.parents[2] / "icl/wgl.py"]
    result = {
        "schema": "astrafek-s23b-icl-target-control-2d-stack-v1",
        "status": "complete_exploratory",
        "n_targets_available": len(target_ids), "n_matched_controls_available": len(controls),
        "n_control_masks_available": len(control_mask_records),
        "n_targets_without_control": len(unmatched_targets), "targets_without_control": unmatched_targets,
        "n_target_frames_excluded_zero_support": len(EXPECTED_ZERO_SUPPORT_IDS),
        "target_zero_support_exclusions": sorted(EXPECTED_ZERO_SUPPORT_IDS),
        "n_pairs_excluded_filter_mismatch": sum(row["reason"] == "target_control_filter_mismatch" for row in exclusions),
        "n_pairs_excluded_zero_support": sum(row["reason"] == "target_frame_zero_unusable_pixels" for row in exclusions),
        "n_pairs_excluded_zero_common_mask_support": sum(
            row["reason"] == "zero_common_mask_support_within_1000_kpc" for row in exclusions),
        "exclusions": exclusions, "n_pairs_stacked": len(pair_rows), "n_by_filter": counts,
        "mask_factor": float(args.mask_factor),
        "grid": {"coordinate": "proper_kpc", "kpc_max": args.kpc_max,
                 "n_side": args.n_side, "pixel_step_kpc": float(kpc_axis[1] - kpc_axis[0]),
                 "radial_profile_edges_kpc": PROFILE_EDGES_KPC.tolist()},
        "stack_definition": {
            "calibration": "HSC FLUXMAG0 to nJy/pixel, then pixel-area division in icl.wgl.warp_metric",
            "mask": f"target proposed Kron-{args.mask_factor:g} mask; control SEP/Kron-{args.mask_factor:g}; native quality/bright-object pixels retained",
            "control_mask_source": "precomputed control masks tied to the complete v2 download, its FITS/WCS audit, and mask-support audit",
            "support": "strict intersection of valid target and matched-random-control pixels after metric-grid interpolation",
            "weighting": "equal pair weight at each grid pixel; no inverse-variance weighting",
            "radial_profiles": "one common-support annular mean per target/control pair and 25-kpc bin from 0 to 1,000 proper kpc",
            "filter_policy": "HSC i and i2 stacked and reported separately",
        },
        "outputs": {"maps": output_maps, "paired_profiles": str(pair_path),
                    "paired_profiles_sha256": checked_sha256(pair_path),
                    "stacked_profiles": str(profile_path),
                    "stacked_profiles_sha256": checked_sha256(profile_path)},
        "pair_support_audit": {
            "summary": str(support_audit_path),
            "summary_sha256": checked_sha256(support_audit_path),
            "paired_support_table": str(pair_support_table_path),
            "paired_support_table_sha256": checked_sha256(pair_support_table_path),
        },
        "input_sha256": input_hashes,
        "script_sha256": {str(path): checked_sha256(path) for path in scripts},
        "software": {"numpy": np.__version__, "python": sys.version.split()[0]},
        "science_measurement": False,
        "images_ready_for_science": False,
        "limitations": [
            "This is an exploratory S23B stack; a complete FITS/WCS audit is not a pixel-level science-quality approval.",
            "The four target cutouts with no usable pixels are excluded; two targets lack a matched random control.",
            "Pairs with no common mask support anywhere inside 1,000 proper kpc are excluded by the pre-stack mask-only audit.",
            "A strict common-support estimator can change the sampled pixels with radius; inspect the per-pair coverage table.",
            "No extinction correction, residual-sky transfer, coadd covariance correction, or PSF-wing matching is applied.",
            "Automatic center detections are not manually confirmed BCG identities; the masks remain proposed masks.",
            "Between-pair SEM is not the full uncertainty and omits correlated coadd noise, sky, PSF, calibration, and astrophysical covariance.",
            "The control positions are selected random sky positions outside the CAMIRA exclusion; unmodeled foreground/background structure can remain.",
            "This stack is not the weak-lensing source stack and does not validate lensing calibration.",
        ],
        "elapsed_seconds": time.monotonic() - start,
    }
    stack.write_json(work / "summary.json", result)
    (work / "progress.json").unlink(missing_ok=True)
    print(json.dumps({key: result[key] for key in (
        "status", "n_pairs_stacked", "n_by_filter", "n_pairs_excluded_filter_mismatch",
        "grid", "science_measurement", "elapsed_seconds")}, indent=2), flush=True)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--target-download-dir", type=Path, required=True)
    parser.add_argument("--target-audit-summary", type=Path, required=True)
    parser.add_argument("--target-mask-dir", type=Path, required=True,
                        help="completed target mask preparation directory")
    parser.add_argument("--target-support-assessment", type=Path, required=True)
    parser.add_argument("--control-manifest", type=Path, required=True,
                        help="approved v2 control manifest")
    parser.add_argument("--control-download-dir", type=Path, required=True)
    parser.add_argument("--control-audit-summary", type=Path, required=True)
    parser.add_argument("--control-mask-dir", type=Path, required=True,
                        help="complete mask-support preparation directory for the approved controls")
    parser.add_argument("--pair-support-audit-summary", type=Path, required=True,
                        help="completed mask-only paired-area audit reviewed before flux stacking")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--kpc-max", type=float, default=1000.0)
    parser.add_argument("--n-side", type=int, default=401)
    parser.add_argument("--mask-factor", type=float, choices=stack.MASK_FACTORS, default=2.5)
    args = parser.parse_args(argv)
    if not np.isfinite(args.kpc_max) or args.kpc_max != 1000.0 or args.n_side < 101 or args.n_side % 2 != 1:
        parser.error("S23B cutouts require --kpc-max 1000 and an odd --n-side of at least 101")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
