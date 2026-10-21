#!/usr/bin/env python3
"""Independently reproduce the all-center public/local PSF overlap audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
EXPECTED_PUBLIC_SHA256 = "8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85"
FIT_WINDOWS = {
    "inner_1p5_2p4_arcsec": (1.50, 2.40),
    "middle_1p8_2p7_arcsec": (1.80, 2.70),
    "outer_2p4_3p3_arcsec": (2.40, 3.30),
    "full_1p5_3p3_arcsec": (1.50, 3.30),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _independent_profile(image: np.ndarray, edges: np.ndarray, pixel_scale: float):
    arr = np.asarray(image, dtype=np.float64)
    yy, xx = np.ogrid[:arr.shape[0], :arr.shape[1]]
    radius = np.sqrt((xx - (arr.shape[1] - 1) / 2) ** 2
                     + (yy - (arr.shape[0] - 1) / 2) ** 2) * pixel_scale
    values, counts = [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = (radius >= lower) & (radius < upper)
        counts.append(int(selected.sum()))
        values.append(float(arr[selected].mean()) if selected.any() else float("nan"))
    return np.asarray(values), np.asarray(counts, dtype=int)


def _independent_fit(local: np.ndarray, public: np.ndarray, counts: np.ndarray,
                     edges: np.ndarray, window: tuple[float, float]):
    bins = (edges[:-1] >= window[0] - 1e-12) & (edges[1:] <= window[1] + 1e-12)
    l, p, n = local[bins], public[bins], counts[bins].astype(float)
    valid = np.isfinite(l) & np.isfinite(p) & (n > 0)
    l, p, n = l[valid], p[valid], n[valid]
    scale = float(np.sum(n * l * p) / np.sum(n * p * p))
    resid = l - scale * p
    frac = float(np.sqrt(np.sum(n * resid * resid) / np.sum(n * l * l)))
    max_scaled = float(np.max(np.abs(resid)) / np.sqrt(np.sum(n * l * l) / np.sum(n)))
    return scale, frac, max_scaled


def _assert_close(actual: float, expected: float, label: str):
    if not np.isclose(float(actual), float(expected), rtol=2e-12, atol=2e-14):
        raise ValueError(f"independent mismatch for {label}: {actual} != {expected}")


def validate(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    summary_path = output_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-public-psf-full-target-control-overlap-v1"
            or summary.get("status") != "complete_filter_stratified_overlap_sensitivity"
            or summary.get("science_measurement") is not False
            or summary.get("psf_approved") is not False):
        raise ValueError("unexpected or promoted analysis summary")

    inputs = summary["inputs"]
    core_summary_path = Path(inputs["core_summary"])
    core_table_path = Path(inputs["core_table"])
    profile_summary_path = Path(inputs["profile_summary"])
    profile_table_path = Path(inputs["profile_table"])
    control_manifest_path = Path(inputs["control_manifest"])
    public_path = Path(inputs["public_psf"])
    for path, key in ((core_summary_path, "core_summary_sha256"),
                      (core_table_path, "core_table_sha256"),
                      (profile_summary_path, "profile_summary_sha256"),
                      (profile_table_path, "profile_table_sha256"),
                      (control_manifest_path, "control_manifest_sha256"),
                      (public_path, "public_psf_sha256")):
        if sha256(path) != inputs[key]:
            raise ValueError(f"input checksum mismatch: {path}")
    if sha256(public_path) != EXPECTED_PUBLIC_SHA256:
        raise ValueError("public PSF is not the approved symmetric HSC-I reference")

    core_summary = json.loads(core_summary_path.read_text(encoding="utf-8"))
    core_table = _read_tsv(core_table_path)
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    profile_rows = _read_tsv(profile_table_path)
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    if (core_summary.get("science_measurement") is not False
            or core_summary.get("psf_approved") is not False
            or profile_summary.get("science_gates", {}).get("science_measurement") is not False):
        raise ValueError("an input science gate was unexpectedly promoted")
    if profile_summary.get("output_sha256", {}).get(str(profile_table_path)) != sha256(profile_table_path):
        raise ValueError("profile table is not tied to its stack summary")

    target_filter = {str(row["id"]): str(row["filter"]) for row in profile_rows}
    control_filter = {}
    control_target = {}
    for row in control_manifest["clusters"]:
        cid = str(row["id"])
        control_filter[cid] = str(row["coadd_header_filter"])
        control_target[cid] = str(row["target_cluster_id"])
    new_audit_dir = Path(core_summary["inputs"]["target_control_audit"]).parent
    old_audit_dir = Path(core_summary["inputs"]["injection_audit"]).parent
    core_table_records = {r["file"]: r for r in core_table}
    positions_path = output_dir / "target_control_public_psf_overlap.tsv"
    pairs_path = output_dir / "paired_target_control_public_psf_overlap.tsv"
    if sha256(positions_path) != summary["tables"]["positions"]["sha256"]:
        raise ValueError("position table checksum mismatch")
    if sha256(pairs_path) != summary["tables"]["pairs"]["sha256"]:
        raise ValueError("paired table checksum mismatch")
    positions = _read_tsv(positions_path)
    pairs = _read_tsv(pairs_path)
    if len(positions) != 1598 or len(pairs) != 783:
        raise ValueError(f"unexpected row counts: positions={len(positions)}, pairs={len(pairs)}")

    edges = np.asarray(inputs["radial_edges_arcsec"], dtype=float)
    public_image = np.asarray(fits.getdata(public_path, ext=1), dtype=np.float64)
    public_crop = public_image[1979:2022, 1979:2022]
    public_profile, public_counts = _independent_profile(
        public_crop, edges, float(inputs["pixel_scale_arcsec"]))
    positions_by_key = {}
    for index, row in enumerate(positions):
        source = row["psf_source"]
        source_dir = new_audit_dir if source == "target_control_center_batch" else old_audit_dir
        local_path = source_dir / row["psf_file"]
        core_record = core_table_records.get(row["psf_file"])
        if (core_record is None or core_record["sha256"] != row["psf_sha256"]
                or sha256(local_path) != row["psf_sha256"]):
            raise ValueError(f"position PSF provenance mismatch: {row['psf_file']}")
        local_image = np.asarray(fits.getdata(local_path, ext=0), dtype=np.float64)
        local_profile, local_counts = _independent_profile(
            local_image, edges, float(inputs["pixel_scale_arcsec"]))
        if not np.array_equal(local_counts, public_counts):
            raise ValueError(f"independent annular support mismatch: {row['psf_file']}")
        for name, window in FIT_WINDOWS.items():
            scale, frac, max_scaled = _independent_fit(
                local_profile, public_profile, local_counts, edges, window)
            _assert_close(row[f"scale_{name}"], scale, f"{index}/{name}")
            if name == "full_1p5_3p3_arcsec":
                _assert_close(row["shape_residual_rms_fraction"], frac, f"{index}/shape")
                _assert_close(row["max_abs_residual_over_local_rms"], max_scaled,
                              f"{index}/max_residual")
        role_data = json.loads(core_record["roles_json"])
        role = role_data[0]
        if role["role"] == "target_cluster_center":
            ident = str(role["cluster_id"])
            expected_filter = target_filter[ident]
            expected_role = "target"
        else:
            ident = str(role["target_cluster_id"])
            cid = str(role["control_id"])
            if control_target.get(cid) != ident:
                raise ValueError(f"control target mapping mismatch: {cid}")
            expected_filter = control_filter[cid]
            expected_role = "control"
        if row["target_cluster_id"] != ident or row["filter"] != expected_filter or row["sample_role"] != expected_role:
            raise ValueError(f"filter/role mapping mismatch: {row['psf_file']}")
        positions_by_key[(ident, expected_role)] = row

    expected_pairs = []
    mismatches = []
    for ident in sorted({key[0] for key in positions_by_key if key[1] == "control"}):
        target = positions_by_key.get((ident, "target"))
        control = positions_by_key[(ident, "control")]
        if target is None:
            raise ValueError(f"matched control missing target: {ident}")
        if target["filter"] != control["filter"]:
            mismatches.append(ident)
            continue
        expected_pairs.append((ident, target, control))
    if len(expected_pairs) != 783 or len(mismatches) != 15:
        raise ValueError("same-filter pair support does not reproduce")
    if len(pairs) != len(expected_pairs):
        raise ValueError("paired result table has the wrong support")
    for result_row, (ident, target, control) in zip(pairs, expected_pairs):
        if str(result_row["target_cluster_id"]) != ident or result_row["filter"] != target["filter"]:
            raise ValueError("paired result ordering or filter mismatch")
        delta = (float(target["full_window_scale_local_over_public"])
                 - float(control["full_window_scale_local_over_public"]))
        _assert_close(result_row["target_minus_control_scale"], delta, f"pair/{ident}")

    return {"status": "independent_validation_passed", "positions_recomputed": len(positions),
            "same_filter_pairs_reproduced": len(expected_pairs), "filter_mismatches_reproduced": len(mismatches),
            "radial_fit_windows_recomputed": len(FIT_WINDOWS),
            "science_measurement": False, "psf_approved": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "extended_psf_target_control_overlap_2026-10-21_02")
    args = parser.parse_args()
    result = validate(args.output_dir)
    validation_path = args.output_dir / "validation.json"
    with validation_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    result["validation_path"] = str(validation_path)
    result["validation_sha256"] = sha256(validation_path)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
