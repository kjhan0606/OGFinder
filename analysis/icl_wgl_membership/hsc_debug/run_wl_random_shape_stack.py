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
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from run_wl_shape_stack import verify_and_load_pairs as verify_and_load_lens_pairs
from s19a_shape_stack import (
    decode_hsc_b_mode_mask,
    load_hsc_y3_nz_sacc,
    paired_tomographic_endpoint_jackknife,
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
        "matched_lens_ids": sorted({str(center["matched_lens_id"]) for center in centers}),
        "random_center_matched_lens_ids": {
            str(center["random_center_id"]): str(center["matched_lens_id"])
            for center in centers
        },
        "random_center_expected_source_counts": {
            str(center["random_center_id"]): int(counts[str(center["random_center_id"])]["n_shape_sources"])
            for center in centers
        },
        "b_mode_masked_value": int(masked_value),
    }


def validate_cluster_stack_compatibility(lens_payload, nz_sacc_sha256, masked_value,
                                         radial_edges_mpc, h0, omega_m):
    if (lens_payload.get("schema") != "hsc-s19a-cluster-lens-stack-run-v1"
            or lens_payload.get("science_measurement") is not False
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


def verify_saved_cluster_profile(saved_stack, recomputed_profile):
    saved = saved_stack.get("profile")
    if not isinstance(saved, list) or len(saved) != len(recomputed_profile):
        raise ValueError("saved cluster profile has an unexpected radial-bin layout")
    for index, (saved_row, recomputed_row) in enumerate(zip(saved, recomputed_profile)):
        recomputed_t = recomputed_row.get(
            "lens_delta_sigma_msun_mpc2",
            recomputed_row.get("delta_sigma_weak_msun_mpc2"),
        )
        recomputed_x = recomputed_row.get(
            "lens_cross_delta_sigma_msun_mpc2",
            recomputed_row.get("cross_delta_sigma_weak_msun_mpc2"),
        )
        if (saved_row.get("r_min_mpc") != recomputed_row["r_min_mpc"]
                or saved_row.get("r_max_mpc") != recomputed_row["r_max_mpc"]
                or saved_row.get("delta_sigma_weak_msun_mpc2") is None
                or saved_row.get("cross_delta_sigma_weak_msun_mpc2") is None
                or not np.isclose(
                    saved_row["delta_sigma_weak_msun_mpc2"],
                    recomputed_t, rtol=1e-10, atol=1e-12,
                )
                or not np.isclose(
                    saved_row["cross_delta_sigma_weak_msun_mpc2"],
                    recomputed_x, rtol=1e-10, atol=1e-12,
                )):
            raise ValueError(f"verified target source rows do not reproduce saved radial bin {index}")


def run(args):
    if args.masked_value not in (0, 1):
        raise ValueError("--b-mode-masked-value must be 0 or 1, verified from the HSC catalog")
    if not args.nz_sacc.is_file():
        raise FileNotFoundError(f"local HSC Y3 n(z) SACC file not found: {args.nz_sacc}")
    lens_row_args = (args.lens_rows, args.lens_preflight, args.lens_plan)
    if args.lens_stack is None and any(value is not None for value in lens_row_args):
        raise ValueError("--lens-rows, --lens-preflight, and --lens-plan require --lens-stack")
    if args.lens_stack is not None and (
            any(value is None for value in lens_row_args)
            or not np.isfinite(args.block_width_deg) or args.block_width_deg <= 0):
        raise ValueError("paired covariance requires verified lens inputs and positive --block-width-deg")
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
    paired_jackknife = None
    matched_random_stack = None
    matched_sample = None
    target_stack_sha = None
    lens_row_provenance = None
    if args.lens_stack is not None:
        if not all((args.lens_rows, args.lens_preflight, args.lens_plan)):
            raise ValueError("paired covariance requires --lens-rows, --lens-preflight, and --lens-plan")
        lens_path = args.lens_stack.resolve()
        lens_payload = json.loads(lens_path.read_text(encoding="utf-8"))
        validate_cluster_stack_compatibility(
            lens_payload, sha(args.nz_sacc), args.masked_value,
            args.radial_edges_mpc, args.h0, args.omega_m,
        )
        target_stack = lens_payload.get("stack")
        lens_pairs, lens_row_provenance = verify_and_load_lens_pairs(
            args.lens_rows, args.lens_preflight, args.lens_plan, args.masked_value,
        )
        lens_input = lens_payload.get("input", {})
        if (lens_row_provenance["row_plan_sha256"] != lens_input.get("row_plan_sha256")
                or lens_row_provenance["masked_value"] != int(args.masked_value)):
            raise ValueError("verified cluster source rows do not match the stored target stack provenance")

        # First reproduce the published 800-lens target stack, then restrict
        # both the cluster and random estimators to the exact lenses that have
        # five eligible random centers.  The random plan currently covers only
        # 620 of the 800 clusters; subtracting its profile from all 800 lenses
        # would mix different lens samples and bias the comparison.
        full_lens_stack = stack_s19a_tomographic_pairs(
            lens_pairs, z_nz, nz_by_bin, args.radial_edges_mpc, cosmology,
        )
        verify_saved_cluster_profile(target_stack, full_lens_stack["profile"])
        matched_lens_ids = set(provenance["matched_lens_ids"])
        available_lens_ids = {str(value) for value in np.unique(lens_pairs["lens_id"])}
        if len(matched_lens_ids) != int(provenance["n_matched_lenses"]):
            raise ValueError("random-center matched-lens IDs do not match the approved random-center plan")
        center_to_lens = provenance["random_center_matched_lens_ids"]
        center_expected_counts = provenance["random_center_expected_source_counts"]
        random_support_lens_ids = {
            center_to_lens[center_id]
            for center_id, n_sources in center_expected_counts.items()
            if int(n_sources) > 0
        }
        common_lens_ids = matched_lens_ids & available_lens_ids & random_support_lens_ids
        no_shape_pair_lens_ids = matched_lens_ids - available_lens_ids
        no_random_source_lens_ids = matched_lens_ids - random_support_lens_ids
        if not common_lens_ids:
            raise ValueError("no random-center matched lens has source rows in the verified cluster export")
        subset_mask = np.isin(lens_pairs["lens_id"].astype(str), list(common_lens_ids))
        matched_lens_pairs = {name: values[subset_mask] for name, values in lens_pairs.items()}
        observed_matched_ids = {str(value) for value in np.unique(matched_lens_pairs["lens_id"])}
        if observed_matched_ids != common_lens_ids:
            raise ValueError("matched cluster source rows do not contain exactly the common-support lens IDs")

        matched_center_ids = {
            center_id for center_id, lens_id in center_to_lens.items()
            if lens_id in common_lens_ids
        }
        matched_center_counts = Counter(
            lens_id for lens_id in center_to_lens.values() if lens_id in common_lens_ids
        )
        if (len(matched_center_counts) != len(common_lens_ids)
                or set(matched_center_counts.values()) != {5}):
            raise ValueError("every common-support lens must have exactly five matched random centers")
        random_subset_mask = np.isin(pairs["lens_id"].astype(str), list(matched_center_ids))
        matched_random_pairs = {
            name: values[random_subset_mask] for name, values in pairs.items()
        }
        observed_center_ids = {str(value) for value in np.unique(matched_random_pairs["lens_id"])}
        expected_center_ids_with_rows = {
            center_id for center_id in matched_center_ids
            if int(center_expected_counts[center_id]) > 0
        }
        if observed_center_ids != expected_center_ids_with_rows:
            raise ValueError("random source rows do not match the count-only positive-support centers")
        matched_random_stack = stack_s19a_tomographic_pairs(
            matched_random_pairs, z_nz, nz_by_bin, args.radial_edges_mpc, cosmology,
        )
        matched_lens_stack = stack_s19a_tomographic_pairs(
            matched_lens_pairs, z_nz, nz_by_bin, args.radial_edges_mpc, cosmology,
        )
        lens_lens_blocks = equal_area_block_ids(
            matched_lens_pairs["lens_ra_deg"], matched_lens_pairs["lens_dec_deg"],
            cell_width_deg=args.block_width_deg,
        )
        lens_source_blocks = equal_area_block_ids(
            matched_lens_pairs["source_ra_deg"], matched_lens_pairs["source_dec_deg"],
            cell_width_deg=args.block_width_deg,
        )
        random_lens_blocks = equal_area_block_ids(
            matched_random_pairs["lens_ra_deg"], matched_random_pairs["lens_dec_deg"],
            cell_width_deg=args.block_width_deg,
        )
        random_source_blocks = equal_area_block_ids(
            matched_random_pairs["source_ra_deg"], matched_random_pairs["source_dec_deg"],
            cell_width_deg=args.block_width_deg,
        )
        paired_jackknife = paired_tomographic_endpoint_jackknife(
            matched_lens_pairs, matched_random_pairs,
            lens_lens_blocks, lens_source_blocks,
            random_lens_blocks, random_source_blocks, z_nz, nz_by_bin,
            args.radial_edges_mpc, cosmology,
            block_method=(f"matched delete-one equal-area sky cells; "
                          f"RA/sin(Dec) width {args.block_width_deg:g} deg"),
        )
        verify_saved_cluster_profile(matched_lens_stack, paired_jackknife["profile"])
        subtracted = subtract_random_lens_profile(matched_lens_stack, matched_random_stack)
        subtracted_profile = subtracted["profile"]
        target_stack_sha = sha(lens_path)
        matched_sample = {
            "n_selected_matched_lenses": len(matched_lens_ids),
            "n_common_shape_support_lenses": len(common_lens_ids),
            "n_selected_lenses_without_cluster_shape_pairs": len(no_shape_pair_lens_ids),
            "n_selected_lenses_without_random_shape_pairs": len(no_random_source_lens_ids),
            "n_selected_lenses_without_common_support": len(matched_lens_ids - common_lens_ids),
            "n_random_centers_per_common_lens": next(iter(matched_center_counts.values())),
            "n_random_centers_in_common_support": len(matched_center_ids),
            "n_random_centers_with_shape_pairs_in_common_support": len(expected_center_ids_with_rows),
            "n_random_centers_without_shape_pairs_in_common_support": len(matched_center_ids - expected_center_ids_with_rows),
            "common_lens_ids_sha256": hashlib.sha256(
                "\n".join(sorted(common_lens_ids)).encode("utf-8")
            ).hexdigest(),
            "no_shape_pair_lens_ids_sha256": hashlib.sha256(
                "\n".join(sorted(no_shape_pair_lens_ids)).encode("utf-8")
            ).hexdigest(),
            "no_random_shape_pair_lens_ids_sha256": hashlib.sha256(
                "\n".join(sorted(no_random_source_lens_ids)).encode("utf-8")
            ).hexdigest(),
            "random_subtraction_sample": "cluster and random profiles restricted to identical matched lenses with at least one verified cluster source pair",
        }

    payload = {
        "schema": "hsc-s19a-random-center-stack-run-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": args.lens_stack is not None,
        "spatial_block_covariance_applied": paired_jackknife is not None,
        "photoz_uncertainty_propagated": False,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_Mpc": args.h0, "Omega_m": args.omega_m},
            "target_stack_path": str(args.lens_stack.resolve()) if args.lens_stack else None,
            "target_stack_sha256": target_stack_sha,
            "lens_row_plan_sha256": lens_row_provenance["row_plan_sha256"] if lens_row_provenance else None,
            "matched_lens_ids_sha256": hashlib.sha256(
                "\n".join(provenance["matched_lens_ids"]).encode("utf-8")
            ).hexdigest() if lens_row_provenance else None,
            "n_matched_lenses_in_cluster_and_random_stacks": (
                matched_sample["n_common_shape_support_lenses"] if lens_row_provenance else None
            ),
            "matched_sample": matched_sample,
            "matched_random_source_pairs": (
                int(matched_random_stack["counts"]["input_pairs"])
                if matched_random_stack is not None else None
            ),
            "full_target_lens_count": target_stack.get("counts", {}).get("n_lenses") if args.lens_stack else None,
            "block_width_deg": args.block_width_deg if paired_jackknife else None,
        },
        "random_stack": random_stack,
        "matched_random_stack": matched_random_stack,
        "random_subtracted_profile": subtracted_profile,
        "paired_spatial_jackknife": paired_jackknife,
        "limitations": [
            "Random centers are only as valid as their upstream full-aperture depth, mask, PSF, and sky matching gates.",
            "Equal-area blocks can cross disconnected HSC regions; block-width sensitivity remains required.",
            "Photo-z calibration uncertainty and reduced shear are not propagated here.",
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
        "common_support_lenses": matched_sample["n_common_shape_support_lenses"] if matched_sample else None,
        "common_support_random_centers": matched_sample["n_random_centers_in_common_support"] if matched_sample else None,
        "common_support_random_pairs": matched_random_stack["counts"]["input_pairs"] if matched_random_stack else None,
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
    parser.add_argument("--lens-rows", type=Path,
                        help="verified cluster source-row export; required with --lens-stack")
    parser.add_argument("--lens-preflight", type=Path,
                        help="verified cluster count/preflight directory; required with --lens-stack")
    parser.add_argument("--lens-plan", type=Path,
                        help="verified cluster row-export plan directory; required with --lens-stack")
    parser.add_argument("--block-width-deg", type=float, default=5.0)
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
