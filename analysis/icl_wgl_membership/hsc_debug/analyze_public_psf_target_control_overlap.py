#!/usr/bin/env python3
"""Compare public extended HSC-I wings with all 1,598 target/control PSF cores.

This is a filter-stratified overlap sensitivity diagnostic only. It never
constructs a hybrid kernel or modifies an observed image.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from audit_public_psf_picker_core_overlap import (
    DEFAULT_RADIAL_EDGES_ARCSEC,
    PIXEL_SCALE_ARCSEC,
    _read_psf,
    fit_profile_scale,
    radial_mean_profile,
)


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_CORE_DIR = ROOT / "target_control_psf_core_diagnostic_2026-10-21_01"
DEFAULT_PROFILE_DIR = ROOT / "prestack_flux_profiles_2026-10-18_02"
DEFAULT_PUBLIC_PSF = ROOT / "extended_psf_public_reference_2026-10-21_01/HSC_PSF_i.fits"
DEFAULT_NEW_CORE_DIR = ROOT / "icl_psf_picker_target_control_centers_download_2026-10-21_01"
DEFAULT_OLD_CORE_DIR = ROOT / "icl_psf_picker_injection_2026-10-19_01"
EXPECTED_CORE_SUMMARY_SHA256 = "cb9244c39e5608932434eec8eb0e7aca0a50ebac585cd769a19f996937843bb0"
EXPECTED_CORE_TABLE_SHA256 = "7c61821b2c73e7a7c3a43bc8c4e58be53e5c7e53ad564cc1d076734a8cb83a27"
EXPECTED_PUBLIC_PSF_SHA256 = "8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85"
DEFAULT_OUTPUT = ROOT / "extended_psf_target_control_overlap_2026-10-21_02"
RADIAL_EDGES = np.asarray(DEFAULT_RADIAL_EDGES_ARCSEC, dtype=float)
FIT_WINDOWS = {
    "inner_1p5_2p4_arcsec": (1.50, 2.40),
    "middle_1p8_2p7_arcsec": (1.80, 2.70),
    "outer_2p4_3p3_arcsec": (2.40, 3.30),
    "full_1p5_3p3_arcsec": (1.50, 3.30),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fit_window_scales(local_profile: np.ndarray, public_profile: np.ndarray,
                      pixel_counts: np.ndarray,
                      radial_edges_arcsec: np.ndarray = RADIAL_EDGES) -> dict:
    local = np.asarray(local_profile, dtype=float)
    public = np.asarray(public_profile, dtype=float)
    counts = np.asarray(pixel_counts, dtype=int)
    edges = np.asarray(radial_edges_arcsec, dtype=float)
    if local.shape != public.shape or local.shape != counts.shape or len(edges) != len(local) + 1:
        raise ValueError("radial profiles, counts, and bin edges are inconsistent")
    result = {}
    for name, (lo, hi) in FIT_WINDOWS.items():
        selected = (edges[:-1] >= lo - 1e-12) & (edges[1:] <= hi + 1e-12)
        if selected.sum() < 3:
            raise ValueError(f"fit window {name} has fewer than three complete bins")
        fit = fit_profile_scale(local[selected], public[selected], counts[selected])
        result[name] = fit["amplitude_scale_local_over_public"]
    return result


def distribution(values) -> dict:
    data = np.asarray(list(values), dtype=float)
    if data.size == 0 or not np.isfinite(data).all():
        raise ValueError("cannot summarize empty or non-finite values")
    return {
        "n": int(data.size),
        "median": float(np.median(data)),
        "q16": float(np.quantile(data, 0.16)),
        "q84": float(np.quantile(data, 0.84)),
        "p05": float(np.quantile(data, 0.05)),
        "p95": float(np.quantile(data, 0.95)),
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def pair_by_actual_filter(targets: dict, controls: dict) -> tuple[list[tuple[str, dict, dict]], list[dict]]:
    """Return exact same-filter pairs and explicit target/control mismatches."""
    matched, mismatched = [], []
    if not controls:
        raise ValueError("matched-control records are empty")
    for ident, control in sorted(controls.items()):
        target = targets.get(ident)
        if target is None:
            raise ValueError(f"control has no corresponding target: {ident}")
        if target["filter"] == control["filter"]:
            matched.append((ident, target, control))
        else:
            mismatched.append({"target_cluster_id": ident,
                               "target_filter": target["filter"],
                               "control_filter": control["filter"]})
    return matched, mismatched


def analyze(core_dir: Path = DEFAULT_CORE_DIR, profile_dir: Path = DEFAULT_PROFILE_DIR,
            public_psf_path: Path = DEFAULT_PUBLIC_PSF,
            new_core_dir: Path = DEFAULT_NEW_CORE_DIR,
            old_core_dir: Path = DEFAULT_OLD_CORE_DIR,
            output_dir: Path = DEFAULT_OUTPUT) -> dict:
    core_dir, profile_dir, public_psf_path, new_core_dir, old_core_dir = map(
        Path, (core_dir, profile_dir, public_psf_path, new_core_dir, old_core_dir)
    )
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_dir}")

    core_summary_path = core_dir / "summary.json"
    core_table_path = core_dir / "center_psf_core_metrics.tsv"
    profile_summary_path = profile_dir / "summary.json"
    profile_table_path = profile_dir / "cluster_profiles.tsv"
    for path in (core_summary_path, core_table_path, profile_summary_path,
                 profile_table_path, public_psf_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(core_summary_path) != EXPECTED_CORE_SUMMARY_SHA256:
        raise ValueError("all-center PSF summary checksum differs from the audited input")
    if sha256(core_table_path) != EXPECTED_CORE_TABLE_SHA256:
        raise ValueError("all-center PSF table checksum differs from the audited input")
    if sha256(public_psf_path) != EXPECTED_PUBLIC_PSF_SHA256:
        raise ValueError("public extended PSF checksum differs from the approved input")

    core_summary = json.loads(core_summary_path.read_text(encoding="utf-8"))
    if (core_summary.get("schema") != "astrafek-hsc-target-control-psf-core-diagnostic-v1"
            or core_summary.get("status") != "complete_descriptive_core_mismatch_diagnostic"
            or core_summary.get("support", {}).get("unique_centers") != 1598
            or core_summary.get("science_measurement") is not False
            or core_summary.get("psf_approved") is not False
            or core_summary.get("tables", {}).get("center_psf_core_metrics", {}).get("sha256")
            != sha256(core_table_path)):
        raise ValueError("all-center PSF input is incomplete or promoted")

    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    if (profile_summary.get("schema") != "hsc-icl-full-sample-profiles-v1"
            or profile_summary.get("status") != "complete"
            or profile_summary.get("science_gates", {}).get("science_measurement") is not False
            or profile_summary.get("output_sha256", {}).get(str(profile_table_path))
            != sha256(profile_table_path)):
        raise ValueError("filter map is not tied to the validated full-sample profile table")
    filter_by_target = {}
    for row in _read_tsv(profile_table_path):
        ident = str(row["id"])
        band = str(row["filter"])
        if band not in {"i", "i2"} or ident in filter_by_target:
            raise ValueError(f"invalid or duplicate target filter mapping: {ident}/{band}")
        filter_by_target[ident] = band

    target_plan_path = Path(core_summary["inputs"]["target_control_plan"])
    if sha256(target_plan_path) != core_summary["inputs"]["target_control_plan_sha256"]:
        raise ValueError("target/control PSF plan checksum differs from the core diagnostic")
    target_plan = json.loads(target_plan_path.read_text(encoding="utf-8"))
    control_manifest_path = Path(target_plan["control_manifest"])
    if sha256(control_manifest_path) != target_plan["control_manifest_sha256"]:
        raise ValueError("matched-control manifest checksum differs from the approved plan")
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_filter_by_id = {}
    control_target_by_id = {}
    for row in control_manifest.get("clusters", []):
        control_id = str(row["id"])
        target_id = str(row["target_cluster_id"])
        band = str(row.get("coadd_header_filter", ""))
        if band not in {"i", "i2"} or control_id in control_filter_by_id:
            raise ValueError(f"invalid or duplicate actual control filter: {control_id}/{band}")
        control_filter_by_id[control_id] = band
        control_target_by_id[control_id] = target_id
    if len(control_filter_by_id) != 798:
        raise ValueError(f"expected 798 actual control filter records, got {len(control_filter_by_id)}")

    public_path_sha = sha256(public_psf_path)
    with fits.open(public_psf_path, memmap=True, checksum=True) as hdus:
        if len(hdus) <= 1 or hdus[1].data is None:
            raise ValueError("public extended PSF lacks image HDU 1")
        public_full = np.asarray(hdus[1].data, dtype=np.float64)
    if public_full.shape != (4001, 4001) or not np.isfinite(public_full).all():
        raise ValueError("public PSF dimensions or finite-pixel check failed")
    public_core = public_full[1979:2022, 1979:2022]
    public_profile, public_counts = radial_mean_profile(public_core, RADIAL_EDGES)

    records = []
    target_by_id = {}
    control_by_target = {}
    core_rows = _read_tsv(core_table_path)
    if len(core_rows) != 1598:
        raise ValueError(f"expected 1,598 unique PSF centers, found {len(core_rows)}")
    for row in core_rows:
        roles = json.loads(row["roles_json"])
        if len(roles) != 1:
            raise ValueError("each deduplicated center must have exactly one sample role")
        role = roles[0]
        if role.get("role") == "target_cluster_center":
            ident = str(role["cluster_id"])
            sample_role = "target"
            if ident in target_by_id:
                raise ValueError(f"duplicate target PSF center: {ident}")
            target_by_id[ident] = row
        elif role.get("role") == "matched_control_center":
            ident = str(role["target_cluster_id"])
            sample_role = "control"
            if ident in control_by_target:
                raise ValueError(f"duplicate matched-control PSF center: {ident}")
            control_id = str(role["control_id"])
            if control_id not in control_filter_by_id or control_target_by_id[control_id] != ident:
                raise ValueError(f"control center is not in the approved control manifest: {control_id}")
            control_by_target[ident] = row
        else:
            raise ValueError(f"unexpected PSF center role: {role}")
        if ident not in filter_by_target:
            raise ValueError(f"PSF center target ID missing from profile table: {ident}")
        band = (filter_by_target[ident] if sample_role == "target"
                else control_filter_by_id[control_id])
        source = str(row["psf_source"])
        if source == "target_control_center_batch":
            local_path = new_core_dir / row["file"]
        elif source == "approved_injection_batch_overlap":
            local_path = old_core_dir / row["file"]
        else:
            raise ValueError(f"unexpected PSF source for {ident}: {source}")
        if not local_path.is_file() or sha256(local_path) != row["sha256"]:
            raise ValueError(f"local PSF core missing or changed: {local_path}")
        local_image = _read_psf(local_path, hdu=0)
        local_profile, local_counts = radial_mean_profile(local_image, RADIAL_EDGES)
        if not np.array_equal(local_counts, public_counts):
            raise ValueError(f"local/public radial pixel supports differ: {local_path.name}")
        full_fit = fit_profile_scale(local_profile, public_profile, local_counts)
        window_scales = fit_window_scales(local_profile, public_profile, local_counts)
        record = {
            "target_cluster_id": ident,
            "filter": band,
            "sample_role": sample_role,
            "ra_deg": row["ra_deg"],
            "dec_deg": row["dec_deg"],
            "psf_source": source,
            "psf_file": row["file"],
            "psf_sha256": row["sha256"],
            "local_shape_yx": row["shape_yx"],
            "full_window_scale_local_over_public": full_fit["amplitude_scale_local_over_public"],
            "shape_residual_rms_fraction": full_fit["weighted_fractional_shape_residual_rms"],
            "max_abs_residual_over_local_rms": full_fit["max_abs_residual_over_weighted_local_rms"],
            **{f"scale_{name}": value for name, value in window_scales.items()},
        }
        records.append(record)

    if len(target_by_id) != 800 or len(control_by_target) != 798:
        raise ValueError(f"expected 800 targets/798 controls, got {len(target_by_id)}/{len(control_by_target)}")
    target_counts = {band: sum(filter_by_target[key] == band for key in target_by_id)
                     for band in ("i", "i2")}
    if target_counts != {"i": 309, "i2": 491}:
        raise ValueError(f"target filter counts differ from validated stack: {target_counts}")

    record_by_role = {(r["target_cluster_id"], r["sample_role"]): r for r in records}
    role_targets = {ident: record_by_role[(ident, "target")] for ident in target_by_id}
    role_controls = {ident: record_by_role[(ident, "control")] for ident in control_by_target}
    matched_pairs, filter_mismatches = pair_by_actual_filter(role_targets, role_controls)
    pairs = []
    for ident, t, c in matched_pairs:
        pair = {"target_cluster_id": ident, "filter": t["filter"],
                "target_scale": t["full_window_scale_local_over_public"],
                "control_scale": c["full_window_scale_local_over_public"],
                "target_minus_control_scale": (
                    t["full_window_scale_local_over_public"]
                    - c["full_window_scale_local_over_public"]),
                "target_minus_control_shape_residual": (
                    t["shape_residual_rms_fraction"] - c["shape_residual_rms_fraction"])}
        pairs.append(pair)
    if len(pairs) != 783 or len(filter_mismatches) != 15:
        raise ValueError(f"expected 783 same-filter pairs and 15 mismatches, got {len(pairs)} and {len(filter_mismatches)}")

    by_filter_role = {}
    for band in ("i", "i2"):
        by_filter_role[band] = {}
        for role in ("target", "control"):
            subset = [r for r in records if r["filter"] == band and r["sample_role"] == role]
            by_filter_role[band][role] = {
                "n": len(subset),
                "full_window_scale": distribution(r["full_window_scale_local_over_public"] for r in subset),
                "shape_residual_rms_fraction": distribution(r["shape_residual_rms_fraction"] for r in subset),
                "window_scale_medians": {
                    name: float(np.median([r[f"scale_{name}"] for r in subset]))
                    for name in FIT_WINDOWS
                },
            }
    paired_summary = {}
    for band in ("i", "i2"):
        subset = [r for r in pairs if r["filter"] == band]
        paired_summary[band] = {
            "n_pairs": len(subset),
            "target_minus_control_scale": distribution(r["target_minus_control_scale"] for r in subset),
            "absolute_target_control_scale_difference_median": float(np.median(
                np.abs([r["target_minus_control_scale"] for r in subset]))),
            "target_minus_control_shape_residual": distribution(
                r["target_minus_control_shape_residual"] for r in subset),
        }

    output_dir.mkdir(parents=True, exist_ok=False)
    positions_path = output_dir / "target_control_public_psf_overlap.tsv"
    pairs_path = output_dir / "paired_target_control_public_psf_overlap.tsv"
    _write_tsv(positions_path, records)
    _write_tsv(pairs_path, pairs)
    result = {
        "schema": "astrafek-public-psf-full-target-control-overlap-v1",
        "status": "complete_filter_stratified_overlap_sensitivity",
        "inputs": {
            "core_summary": str(core_summary_path.resolve()),
            "core_summary_sha256": sha256(core_summary_path),
            "core_table": str(core_table_path.resolve()),
            "core_table_sha256": sha256(core_table_path),
            "profile_summary": str(profile_summary_path.resolve()),
            "profile_summary_sha256": sha256(profile_summary_path),
            "profile_table": str(profile_table_path.resolve()),
            "profile_table_sha256": sha256(profile_table_path),
            "control_manifest": str(control_manifest_path.resolve()),
            "control_manifest_sha256": sha256(control_manifest_path),
            "public_psf": str(public_psf_path.resolve()),
            "public_psf_sha256": public_path_sha,
            "radial_edges_arcsec": RADIAL_EDGES.tolist(),
            "pixel_scale_arcsec": PIXEL_SCALE_ARCSEC,
        },
        "support": {"target_centers": 800, "matched_control_centers": 798,
                    "unique_centers": len(records), "same_filter_pairs": len(pairs),
                    "filter_mismatched_pairs": len(filter_mismatches),
                    "filter_mismatch_records": filter_mismatches,
                    "target_filter_counts": target_counts,
                    "control_filter_counts": {band: sum(r["filter"] == band and r["sample_role"] == "control" for r in records)
                                              for band in ("i", "i2")}},
        "fit_windows_arcsec": {name: list(bounds) for name, bounds in FIT_WINDOWS.items()},
        "by_filter_and_role": by_filter_role,
        "paired_target_control": paired_summary,
        "tables": {
            "positions": {"path": str(positions_path), "sha256": sha256(positions_path), "rows": len(records)},
            "pairs": {"path": str(pairs_path), "sha256": sha256(pairs_path), "rows": len(pairs)},
        },
        "science_measurement": False,
        "psf_approved": False,
        "limitations": [
            "The HSC-I reference is a survey-level PSF model, not a local per-frame wing measurement.",
            "Overlap scales and shape residuals are sensitivity diagnostics; no hybrid kernel is constructed or applied.",
            "The local core stamps end near 3.6 arcsec and do not validate wings at the 200-300 kpc ICL radii.",
            "No image pixels, cluster ICL values, source shapes, or photo-z data are read by this analysis.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, default=DEFAULT_CORE_DIR)
    parser.add_argument("--profile-dir", type=Path, default=DEFAULT_PROFILE_DIR)
    parser.add_argument("--public-psf", type=Path, default=DEFAULT_PUBLIC_PSF)
    parser.add_argument("--new-core-dir", type=Path, default=DEFAULT_NEW_CORE_DIR)
    parser.add_argument("--old-core-dir", type=Path, default=DEFAULT_OLD_CORE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = analyze(args.core_dir, args.profile_dir, args.public_psf,
                     args.new_core_dir, args.old_core_dir, args.output_dir)
    print(json.dumps({"status": result["status"], "support": result["support"],
                      "paired_target_control": result["paired_target_control"],
                      "science_measurement": result["science_measurement"],
                      "psf_approved": result["psf_approved"],
                      "summary_sha256": sha256(Path(args.output_dir) / "summary.json")},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
