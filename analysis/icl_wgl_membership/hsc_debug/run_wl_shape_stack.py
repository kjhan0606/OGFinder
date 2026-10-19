#!/usr/bin/env python3
"""Run the local HSC Y3 tomographic stack on a verified S19A row export.

This command performs no network access and never downloads data. It requires
an already-complete, checksum-verified row export, a local HSC Y3 SACC n(z)
file, and an explicitly verified stored-value polarity for ``b_mode_mask``.
The cluster-lens-only output remains exploratory until random-lens subtraction,
photo-z uncertainty, and spatial-block covariance are included.
"""
import argparse
from collections import Counter
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
from prepare_wl_shape_rows import OUTPUT_COLUMNS
from s19a_shape_stack import (
    decode_hsc_b_mode_mask,
    load_hsc_y3_nz_sacc,
    stack_s19a_tomographic_pairs,
)
from submit_wl_shape_rows import DEFAULT_PLAN, DEFAULT_PREFLIGHT, validate_export_plan
from submit_wl_shape_counts import load_plan
from ogfkit.wlcluster import Cosmology


PAIR_FIELDS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "weight", "rms_e",
    "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12",
)
ID_FIELDS = {"lens_id", "object_id"}
BIN_FIELDS = {"source_zbin"}


def _arrays_from_csv(reader, masked_value, lens_counts):
    if reader.fieldnames != OUTPUT_COLUMNS:
        raise ValueError(f"unexpected source-row columns: {reader.fieldnames}")
    columns = {name: [] for name in (*PAIR_FIELDS, "source_zbin", "b_mode_mask_raw")}
    observed_lenses = Counter()
    n_rows = 0
    for row in reader:
        if None in row or not row.get("lens_id") or not row.get("object_id"):
            raise ValueError("source-row CSV contains a malformed lens or object ID")
        lens_id = str(row["lens_id"])
        if lens_id not in lens_counts:
            raise ValueError(f"source row has an unrequested lens_id: {lens_id}")
        observed_lenses[lens_id] += 1
        n_rows += 1
        for name in PAIR_FIELDS:
            value = row[name]
            columns[name].append(value if name in ID_FIELDS else float(value))
        columns["source_zbin"].append(float(row["source_zbin"]))
        columns["b_mode_mask_raw"].append(row["b_mode_mask_raw"])

    arrays = {}
    for name in PAIR_FIELDS:
        dtype = object if name in ID_FIELDS else float
        arrays[name] = np.asarray(columns[name], dtype=dtype)
    bins = np.asarray(columns["source_zbin"], dtype=float)
    if not np.all(np.isfinite(bins)) or not np.array_equal(bins, bins.astype(int)):
        raise ValueError("source_zbin must contain finite integer values")
    arrays["source_zbin"] = bins.astype(int)
    raw_mask = np.asarray(columns["b_mode_mask_raw"], dtype=float)
    if not np.all(np.isfinite(raw_mask)) or np.any((raw_mask != 0) & (raw_mask != 1)):
        raise ValueError("raw b_mode_mask is not binary; verify its catalog definition")
    arrays["b_mode_mask_raw"] = raw_mask.astype(np.int8)
    if masked_value is not None:
        arrays["b_mode_pass"] = decode_hsc_b_mode_mask(
            raw_mask, masked_value=masked_value,
        )
    return arrays, observed_lenses, n_rows


def _join_arrays(parts):
    if not parts:
        raise ValueError("source-row export contains no batches")
    keys = tuple(parts[0])
    if any(tuple(part) != keys for part in parts[1:]):
        raise ValueError("source-row batches have inconsistent columns")
    return {key: np.concatenate([part[key] for part in parts]) for key in keys}


def verify_and_load_pairs(rows_dir, preflight, plan_dir, masked_value):
    """Verify the export and load rows; ``None`` validates but does not decode b-mode."""
    rows_dir = Path(rows_dir).resolve()
    preflight = Path(preflight).resolve()
    plan_dir = Path(plan_dir).resolve()
    manifest, lens_counts, count_summary = validate_export_plan(plan_dir, preflight)
    plan_sha = hashlib.sha256((plan_dir / "row_export_plan.json").read_bytes()).hexdigest()
    state_path = rows_dir / "row_status.json"
    summary_path = rows_dir / "row_summary.json"
    if not state_path.is_file() or not summary_path.is_file():
        raise ValueError("a completed row-export checkpoint and summary are required")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if state.get("schema") != "hsc-s19a-shape-row-run-v1" or state.get("status") != "complete":
        raise ValueError("row-export checkpoint is not complete")
    if state.get("plan_sha256") != plan_sha:
        raise ValueError("row-export checkpoint belongs to a different SQL plan")
    if state.get("source_rows_downloaded") is not True or state.get("credentials_saved") is not False:
        raise ValueError("row-export checkpoint has unexpected transfer/credential state")
    expected_total = int(manifest["lens_source_pair_rows_expected"])
    if (summary.get("status") != "complete" or summary.get("n_pair_rows") != expected_total
            or summary.get("science_measurement") is not False
            or summary.get("photoz_pdf_rows_included") is not False):
        raise ValueError("row-export summary does not certify the expected exploratory export")
    if summary.get("count_summary_sum") != count_summary.get("sum_lens_cone_shape_counts"):
        raise ValueError("row-export summary disagrees with the count-only result")
    jobs = state.get("jobs", {})
    batches = manifest.get("batches", [])
    if len(batches) != 32 or len(jobs) != 32:
        raise ValueError("expected all 32 verified shape-row batches")

    query_plan = load_plan(preflight)
    batch_size = int(query_plan["batch_size"])
    lens_ids = [str(item["lens_id"]) for item in query_plan["lenses"]]
    if len(lens_ids) != 800 or len(set(lens_ids)) != 800:
        raise ValueError("count plan does not contain 800 unique lenses")

    parts = []
    row_hashes = {}
    observed_total = 0
    for index, batch in enumerate(batches):
        name = f"shape_rows_{index:03d}"
        item = jobs.get(name, {})
        path = rows_dir / f"{name}.csv"
        if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
            raise ValueError(f"missing, incomplete, or changed source-row batch: {name}")
        start = index * batch_size
        batch_lenses = lens_ids[start:start + batch_size]
        expected = {lens_id: lens_counts[lens_id] for lens_id in batch_lenses}
        with path.open(encoding="utf-8-sig", newline="") as stream:
            arrays, observed, n_rows = _arrays_from_csv(
                catalogue_reader(stream), masked_value, lens_counts,
            )
        if observed != Counter(expected):
            raise ValueError(f"per-lens source-row counts differ from count preflight: {name}")
        if n_rows != int(batch["expected_pair_rows_from_count_only"]) or item.get("rows") != n_rows:
            raise ValueError(f"source-row count mismatch: {name}")
        parts.append(arrays)
        row_hashes[name] = item["sha256"]
        observed_total += n_rows
    if observed_total != expected_total or observed_total != int(state.get("total_pair_rows", -1)):
        raise ValueError("combined source-row total differs from the verified export")
    return _join_arrays(parts), {
        "row_plan_sha256": plan_sha,
        "row_batch_sha256": row_hashes,
        "n_pair_rows": observed_total,
        "n_lenses": 800,
        "count_query_plan_sha256": count_summary["query_plan_sha256"],
        "count_summary_sha256": sha(preflight / "cas_counts" / "count_summary.json"),
        "masked_value": None if masked_value is None else int(masked_value),
    }


def run(args):
    if not 0 <= args.masked_value <= 1:
        raise ValueError("--b-mode-masked-value must be 0 or 1 and verified from the catalog")
    if not args.nz_sacc.is_file():
        raise FileNotFoundError(f"local HSC Y3 n(z) SACC file not found: {args.nz_sacc}")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing stack output: {output}")
    pairs, provenance = verify_and_load_pairs(
        args.rows.resolve(), args.preflight, args.plan_dir, args.masked_value,
    )
    z_nz, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    stack = stack_s19a_tomographic_pairs(
        pairs, z_nz, nz_by_bin, args.radial_edges_mpc, cosmology,
    )
    payload = {
        "schema": "hsc-s19a-cluster-lens-stack-run-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "spatial_block_covariance_applied": False,
        "photoz_uncertainty_propagated": False,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
        },
        "stack": stack,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": "exploratory", "science_measurement": False,
        "lens_source_pairs": provenance["n_pair_rows"],
        "used_pairs": stack["counts"]["used_pairs"],
        "radial_bins": len(stack["profile"]), "output": str(output),
    }, ensure_ascii=False))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True, help="completed shape-row export directory")
    parser.add_argument("--nz-sacc", type=Path, required=True, help="local HSC Y3 SACC file with wl_0..wl_3")
    parser.add_argument("--b-mode-masked-value", dest="masked_value", type=int, required=True,
                        help="verified raw b_mode_mask value to exclude (0 or 1)")
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True, help="new JSON output path; refuses overwrite")
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    return parser


def main():
    args = build_parser().parse_args()
    if len(args.radial_edges_mpc) < 2:
        parser.error("provide at least two radial edges")
    run(args)


if __name__ == "__main__":
    main()
