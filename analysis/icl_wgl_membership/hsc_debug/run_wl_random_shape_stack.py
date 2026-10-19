#!/usr/bin/env python3
"""Verify an approved random-shape export and run the local HSC Y3 stack."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import catalogue_reader, sha
from ogfkit.wlcluster import Cosmology
from prepare_wl_random_shape_rows import OUTPUT_COLUMNS
from s19a_shape_stack import (
    decode_hsc_b_mode_mask,
    load_hsc_y3_nz_sacc,
    stack_s19a_tomographic_pairs,
    subtract_random_lens_profile,
)
from submit_wl_random_shape_rows import validate_plan


PAIR_FIELDS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "weight", "rms_e",
    "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12",
)
ID_FIELDS = {"lens_id", "object_id"}


def _arrays_from_csv(reader, masked_value, centers, expected_counts):
    if reader.fieldnames != OUTPUT_COLUMNS:
        raise ValueError(f"unexpected random source-row columns: {reader.fieldnames}")
    center_by_id = {str(row["random_center_id"]): row for row in centers}
    columns = {name: [] for name in (*PAIR_FIELDS, "source_zbin", "b_mode_mask_raw")}
    observed_centers = Counter()
    n_rows = 0
    for row in reader:
        if None in row or not row.get("random_center_id") or not row.get("object_id"):
            raise ValueError("random source-row CSV contains a malformed center or object ID")
        center_id = str(row["random_center_id"])
        center = center_by_id.get(center_id)
        if center is None:
            raise ValueError(f"row contains an unrequested random center: {center_id}")
        if (str(row["matched_lens_id"]) != str(center["matched_lens_id"])
                or abs(float(row["random_ra_deg"]) - float(center["ra_deg"])) > 1e-9
                or abs(float(row["random_dec_deg"]) - float(center["dec_deg"])) > 1e-9
                or abs(float(row["lens_z"]) - float(center["z_lens"])) > 1e-10):
            raise ValueError("source row changed its random-center-to-lens mapping")
        n_rows += 1
        observed_centers[center_id] += 1
        values = {
            "lens_id": center_id,
            "object_id": str(row["object_id"]),
            "lens_ra_deg": float(row["random_ra_deg"]),
            "lens_dec_deg": float(row["random_dec_deg"]),
            "lens_z": float(row["lens_z"]),
        }
        for name in PAIR_FIELDS:
            if name not in values:
                values[name] = float(row[name])
            columns[name].append(values[name])
        columns["source_zbin"].append(float(row["source_zbin"]))
        columns["b_mode_mask_raw"].append(row["b_mode_mask_raw"])
    if any(observed_centers.get(center_id, 0) != int(expected_counts[center_id]["n_shape_sources"])
           for center_id in center_by_id):
        raise ValueError("per-center random shape-row totals differ from the count-only preflight")

    arrays = {}
    for name in PAIR_FIELDS:
        arrays[name] = np.asarray(columns[name], dtype=object if name in ID_FIELDS else float)
    bins = np.asarray(columns["source_zbin"], dtype=float)
    if not np.all(np.isfinite(bins)) or not np.array_equal(bins, bins.astype(int)) or np.any(bins < 0):
        raise ValueError("source_zbin must contain finite non-negative integer values")
    arrays["source_zbin"] = bins.astype(int)
    raw_mask = np.asarray(columns["b_mode_mask_raw"], dtype=float)
    if not np.all(np.isfinite(raw_mask)) or np.any((raw_mask != 0) & (raw_mask != 1)):
        raise ValueError("raw b_mode_mask is not binary; verify its catalog definition")
    arrays["b_mode_mask_raw"] = raw_mask.astype(np.int8)
    arrays["b_mode_pass"] = decode_hsc_b_mode_mask(raw_mask, masked_value=masked_value)
    return arrays, observed_centers, n_rows


def _join_arrays(parts):
    if not parts:
        raise ValueError("random source-row export contains no batches")
    keys = tuple(parts[0])
    if any(tuple(part) != keys for part in parts[1:]):
        raise ValueError("random source-row batches have inconsistent columns")
    return {key: np.concatenate([part[key] for part in parts]) for key in keys}


def verify_and_load_pairs(rows_dir, plan_dir, masked_value):
    rows_dir, plan_dir = Path(rows_dir).resolve(), Path(plan_dir).resolve()
    manifest, centers, batches, counts, count_plan_sha, row_plan_sha = validate_plan(plan_dir)
    state_path, summary_path = rows_dir / "row_status.json", rows_dir / "row_summary.json"
    if not state_path.is_file() or not summary_path.is_file():
        raise ValueError("completed random source-row checkpoint and summary are required")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (state.get("schema") != "hsc-s19a-random-shape-row-run-v1"
            or state.get("status") != "complete"
            or state.get("plan_sha256") != row_plan_sha
            or state.get("count_plan_sha256") != count_plan_sha):
        raise ValueError("random row-export checkpoint is incomplete or belongs to a different plan")
    if (state.get("source_shape_rows_downloaded") is not True
            or state.get("credentials_saved") is not False
            or state.get("photoz_pdf_rows_included") is not False
            or state.get("images_returned") is not False
            or state.get("science_measurement") is not False):
        raise ValueError("row-export checkpoint has unexpected data or gate state")
    if (summary.get("status") != "complete" or summary.get("plan_sha256") != row_plan_sha
            or summary.get("n_source_pairs") != manifest["expected_random_source_pair_rows"]
            or summary.get("source_shape_rows_downloaded") is not True
            or summary.get("photoz_pdf_rows_included") is not False
            or summary.get("images_returned") is not False
            or summary.get("science_measurement") is not False):
        raise ValueError("random row-export summary disagrees with the immutable plan")
    jobs = state.get("jobs", {})
    if len(batches) != manifest.get("n_query_batches") or len(jobs) != len(batches):
        raise ValueError("expected every verified random shape-row batch")

    batch_size = int(manifest["batch_size"])
    parts, row_hashes, observed_total = [], {}, 0
    for index, batch in enumerate(batches):
        name = f"random_shape_rows_{index:03d}"
        item = jobs.get(name, {})
        path = rows_dir / f"{name}.csv"
        if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
            raise ValueError(f"missing, incomplete, or changed random source-row batch: {name}")
        chunk = centers[index * batch_size:(index + 1) * batch_size]
        expected = {center["random_center_id"]: counts[center["random_center_id"]] for center in chunk}
        with path.open(encoding="utf-8-sig", newline="") as stream:
            arrays, observed, n_rows = _arrays_from_csv(
                catalogue_reader(stream), masked_value, chunk, expected,
            )
        if n_rows != int(batch["expected_pair_rows_from_count_only"]) or item.get("rows") != n_rows:
            raise ValueError(f"random source-row total mismatch: {name}")
        parts.append(arrays)
        row_hashes[name] = item["sha256"]
        observed_total += n_rows
    if observed_total != int(manifest["expected_random_source_pair_rows"]):
        raise ValueError("combined random row total differs from the count-only plan")
    return _join_arrays(parts), {
        "random_row_plan_sha256": row_plan_sha,
        "random_row_batch_sha256": row_hashes,
        "count_plan_sha256": count_plan_sha,
        "count_summary_sha256": sha(Path(manifest["source_count_summary"])),
        "n_source_pairs": observed_total,
        "n_random_centers": len(centers),
        "n_matched_lenses": manifest["n_matched_lenses"],
        "b_mode_masked_value": int(masked_value),
    }


def validate_cluster_stack_compatibility(lens_payload, nz_sacc_sha256, masked_value,
                                         radial_edges_mpc, h0, omega_m):
    if (lens_payload.get("schema") != "hsc-s19a-cluster-lens-stack-run-v1"
            or lens_payload.get("random_lens_subtraction_applied") is not False
            or lens_payload.get("spatial_block_covariance_applied") is not False):
        raise ValueError("target stack input is not the verified cluster-only exploratory product")
    lens_input = lens_payload.get("input", {})
    expected_cosmology = {"H0_km_s_mpc": float(h0), "Omega_m": float(omega_m)}
    if (lens_input.get("n_z_sacc_sha256") != nz_sacc_sha256
            or lens_input.get("masked_value") != int(masked_value)
            or lens_input.get("radial_edges_mpc") != [float(x) for x in radial_edges_mpc]
            or lens_input.get("cosmology") != expected_cosmology):
        raise ValueError("target and random stacks do not use identical n(z), mask, bins, and cosmology")


def run(args):
    if args.masked_value not in (0, 1):
        raise ValueError("--b-mode-masked-value must be 0 or 1, verified from the HSC catalog")
    if not args.nz_sacc.is_file():
        raise FileNotFoundError(f"local HSC Y3 n(z) SACC file not found: {args.nz_sacc}")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite random stack output: {output}")
    pairs, provenance = verify_and_load_pairs(args.rows, args.plan_dir, args.masked_value)
    z_nz, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    random_stack = stack_s19a_tomographic_pairs(
        pairs, z_nz, nz_by_bin, args.radial_edges_mpc, cosmology,
    )
    subtracted_profile = None
    target_stack_sha = None
    if args.lens_stack is not None:
        lens_path = args.lens_stack.resolve()
        lens_payload = json.loads(lens_path.read_text(encoding="utf-8"))
        validate_cluster_stack_compatibility(
            lens_payload, sha(args.nz_sacc), args.masked_value,
            args.radial_edges_mpc, args.h0, args.omega_m,
        )
        target_stack = lens_payload.get("stack")
        subtracted = subtract_random_lens_profile(target_stack, random_stack)
        subtracted_profile = subtracted["profile"]
        target_stack_sha = sha(lens_path)

    payload = {
        "schema": "hsc-s19a-random-center-stack-run-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": args.lens_stack is not None,
        "spatial_block_covariance_applied": False,
        "photoz_uncertainty_propagated": False,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_Mpc": args.h0, "Omega_m": args.omega_m},
            "target_stack_path": str(args.lens_stack.resolve()) if args.lens_stack else None,
            "target_stack_sha256": target_stack_sha,
        },
        "random_stack": random_stack,
        "random_subtracted_profile": subtracted_profile,
        "limitations": [
            "Random centers are only as valid as their upstream full-aperture depth, mask, PSF, and sky matching gates.",
            "Random subtraction is a point estimate; field/source-overlap covariance is not included.",
            "Photo-z calibration uncertainty, reduced shear, and spatial-block jackknife are not propagated here.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": "exploratory", "science_measurement": False,
        "random_centers": provenance["n_random_centers"],
        "source_pairs": provenance["n_source_pairs"],
        "used_pairs": random_stack["counts"]["used_pairs"],
        "random_subtraction": args.lens_stack is not None,
        "output": str(output),
    }, ensure_ascii=False))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--b-mode-masked-value", dest="masked_value", type=int, required=True)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lens-stack", type=Path,
                        help="optional verified cluster-only stack JSON for point-estimate subtraction")
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    return parser


def main():
    args = build_parser().parse_args()
    if len(args.radial_edges_mpc) < 2:
        raise SystemExit("provide at least two radial edges")
    run(args)


if __name__ == "__main__":
    main()
