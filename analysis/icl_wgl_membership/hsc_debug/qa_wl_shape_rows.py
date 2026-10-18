#!/usr/bin/env python3
"""Stream basic quality diagnostics over a verified HSC S19A row export.

This is an offline catalogue audit only. It does not decode B-mode polarity,
apply a lensing estimator, or promote the data to a science measurement.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import platform
import statistics
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_wl_shape_rows import DEFAULT_OUTPUT as DEFAULT_PLAN, OUTPUT_COLUMNS
from submit_wl_shape_counts import load_plan
from submit_wl_shape_rows import DEFAULT_PREFLIGHT, validate_export_plan

DEFAULT_ROWS = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_rows_export_2026-10-18_03"
)
DEFAULT_OUTPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_row_quality_2026-10-19_01/summary.json"
)

NUMERIC_FIELDS = (
    "lens_ra_deg", "lens_dec_deg", "lens_z", "source_ra_deg", "source_dec_deg",
    "e1", "e2", "sigma_e", "rms_e", "weight", "m", "c1", "c2",
    "mag_a10", "resolution", "psf_shape11", "psf_shape22", "psf_shape12",
)


def _empty_stats():
    return {"count": 0, "parse_invalid": 0, "nonfinite": 0,
            "min": None, "max": None, "sum": 0.0, "sum_sq": 0.0}


def _add(stats, raw):
    try:
        value = float(raw)
    except (TypeError, ValueError, OverflowError):
        stats["parse_invalid"] += 1
        return None
    if not math.isfinite(value):
        stats["nonfinite"] += 1
        return None
    stats["count"] += 1
    stats["min"] = value if stats["min"] is None else min(stats["min"], value)
    stats["max"] = value if stats["max"] is None else max(stats["max"], value)
    stats["sum"] += value
    stats["sum_sq"] += value * value
    return value


def _finish(stats):
    count = stats["count"]
    total, total_sq = stats.pop("sum"), stats.pop("sum_sq")
    stats["mean"] = total / count if count else None
    stats["rms_about_zero"] = math.sqrt(total_sq / count) if count else None
    return stats


def audit(rows_dir: Path, preflight: Path, plan_dir: Path):
    rows_dir, preflight, plan_dir = map(Path.resolve, (rows_dir, preflight, plan_dir))
    manifest, lens_counts, count_summary = validate_export_plan(plan_dir, preflight)
    state = json.loads((rows_dir / "row_status.json").read_text(encoding="utf-8"))
    export_summary = json.loads((rows_dir / "row_summary.json").read_text(encoding="utf-8"))
    expected_total = int(manifest["lens_source_pair_rows_expected"])
    if (state.get("status") != "complete" or export_summary.get("status") != "complete"
            or state.get("source_rows_downloaded") is not True
            or state.get("credentials_saved") is not False
            or export_summary.get("n_pair_rows") != expected_total
            or export_summary.get("science_measurement") is not False
            or export_summary.get("photoz_pdf_rows_included") is not False):
        raise ValueError("row export is not a complete, scope-matched exploratory export")

    query_plan = load_plan(preflight)
    batches = manifest["batches"]
    if len(batches) != 32 or len(state.get("jobs", {})) != 32:
        raise ValueError("expected all 32 planned row batches")
    numeric = {name: _empty_stats() for name in NUMERIC_FIELDS}
    lens_hist = Counter()
    zbin_hist = Counter()
    bmode_hist = Counter()
    selected_by_zbin = Counter()
    bad_sky_coords = 0
    nonphysical_lens_z = 0
    ellipticity_norm_gt_one = 0
    ellipticity_norm_ge_two = 0
    sigma_e_outside_published_range = 0
    positive_weight = zero_weight = negative_weight = 0
    calibration_cuts = {"zbin_1_to_4": 0, "mag_a10_le_25p5": 0,
                        "resolution_ge_0p3": 0, "weight_gt_0": 0,
                        "sigma_e_le_0p4": 0, "ellipticity_norm_lt_2": 0,
                        "published_shape_cuts_and_tomographic_bin": 0}
    total_rows = 0
    batch_results = []

    for index, batch in enumerate(batches):
        name = f"shape_rows_{index:03d}"
        job = state["jobs"].get(name, {})
        path = rows_dir / f"{name}.csv"
        if (job.get("status") != "complete" or not path.is_file()
                or path.stat().st_size != job.get("bytes")
                or sha(path) != job.get("sha256")):
            raise ValueError(f"missing, incomplete, or wrong-sized batch: {name}")
        observed = Counter()
        seen_pairs = set()
        batch_rows = 0
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = catalogue_reader(stream)
            if reader.fieldnames != OUTPUT_COLUMNS:
                raise ValueError(f"unexpected columns in {name}: {reader.fieldnames}")
            for row in reader:
                if None in row or not row.get("lens_id") or not row.get("object_id"):
                    raise ValueError(f"malformed row ID in {name}")
                lens_id, object_id = str(row["lens_id"]), str(row["object_id"])
                if lens_id not in lens_counts:
                    raise ValueError(f"unrequested lens in {name}: {lens_id}")
                pair = (lens_id, object_id)
                if pair in seen_pairs:
                    raise ValueError(f"duplicate lens-source pair in {name}: {pair}")
                seen_pairs.add(pair)
                observed[lens_id] += 1
                lens_hist[lens_id] += 1
                batch_rows += 1

                values = {field: _add(numeric[field], row.get(field)) for field in NUMERIC_FIELDS}
                zbin_raw = (row.get("source_zbin") or "").strip()
                mask_raw = (row.get("b_mode_mask_raw") or "").strip()
                zbin_hist[zbin_raw] += 1
                bmode_hist[mask_raw] += 1
                if values["lens_z"] is not None and values["lens_z"] <= 0:
                    nonphysical_lens_z += 1
                if (values["lens_ra_deg"] is not None and values["lens_dec_deg"] is not None
                        and not (0 <= values["lens_ra_deg"] < 360
                                 and -90 <= values["lens_dec_deg"] <= 90)):
                    bad_sky_coords += 1
                if values["source_ra_deg"] is not None and values["source_dec_deg"] is not None:
                    if not (0 <= values["source_ra_deg"] < 360
                            and -90 <= values["source_dec_deg"] <= 90):
                        bad_sky_coords += 1
                if values["e1"] is not None and values["e2"] is not None:
                    e_norm = math.hypot(values["e1"], values["e2"])
                    if e_norm > 1:
                        ellipticity_norm_gt_one += 1
                    if e_norm >= 2:
                        ellipticity_norm_ge_two += 1
                else:
                    e_norm = None
                sigma_e = values["sigma_e"]
                if sigma_e is None or sigma_e < 0 or sigma_e > 0.4:
                    sigma_e_outside_published_range += 1
                weight = values["weight"]
                positive_weight += int(weight is not None and weight > 0)
                zero_weight += int(weight == 0)
                negative_weight += int(weight is not None and weight < 0)

                try:
                    zbin = int(zbin_raw)
                    if str(zbin) != zbin_raw or not 0 <= zbin <= 4:
                        raise ValueError
                except ValueError:
                    raise ValueError(f"invalid HSC source_zbin value in {name}: {zbin_raw!r}") from None
                if 1 <= zbin <= 4:
                    calibration_cuts["zbin_1_to_4"] += 1
                mag, resolution = values["mag_a10"], values["resolution"]
                if mag is not None and mag <= 25.5:
                    calibration_cuts["mag_a10_le_25p5"] += 1
                if resolution is not None and resolution >= 0.3:
                    calibration_cuts["resolution_ge_0p3"] += 1
                if sigma_e is not None and 0 <= sigma_e <= 0.4:
                    calibration_cuts["sigma_e_le_0p4"] += 1
                if e_norm is not None and e_norm < 2:
                    calibration_cuts["ellipticity_norm_lt_2"] += 1
                if weight is not None and weight > 0:
                    calibration_cuts["weight_gt_0"] += 1
                if (1 <= zbin <= 4 and mag is not None and mag <= 25.5
                        and resolution is not None and resolution >= 0.3
                        and sigma_e is not None and 0 <= sigma_e <= 0.4
                        and e_norm is not None and e_norm < 2
                        and weight is not None and weight > 0):
                    selected_by_zbin[str(zbin)] += 1
                    calibration_cuts["published_shape_cuts_and_tomographic_bin"] += 1

        expected_lens_ids = [str(row["lens_id"]) for row in
                             query_plan["lenses"][index * int(query_plan["batch_size"]):
                                                  (index + 1) * int(query_plan["batch_size"])] ]
        expected = {lens_id: lens_counts[lens_id] for lens_id in expected_lens_ids}
        if observed != Counter(expected) or batch_rows != int(batch["expected_pair_rows_from_count_only"]):
            raise ValueError(f"per-lens counts differ from count-only results in {name}")
        if batch_rows != job.get("rows"):
            raise ValueError(f"row_status count differs from CSV in {name}")
        total_rows += batch_rows
        batch_results.append({"name": name, "rows": batch_rows, "bytes": path.stat().st_size,
                              "sha256": job["sha256"]})

    if total_rows != expected_total or total_rows != int(state.get("total_pair_rows", -1)):
        raise ValueError("combined row count differs from the approved plan/checkpoint")
    numeric = {name: _finish(stats) for name, stats in numeric.items()}
    per_lens_values = [lens_hist.get(lens_id, 0) for lens_id in lens_counts]
    result = {
                    "schema": "hsc-s19a-shape-row-quality-v2",
        "status": "complete_descriptive",
        "science_measurement": False,
        "source_rows_downloaded": True,
        "n_lenses": len(lens_counts),
        "n_pair_rows": total_rows,
        "total_csv_bytes": sum(row["bytes"] for row in batch_results),
        "mean_bytes_per_pair": sum(row["bytes"] for row in batch_results) / total_rows,
        "per_lens_pair_count": {
            "min": min(per_lens_values), "median": statistics.median(per_lens_values),
            "max": max(per_lens_values), "zero_candidate_lenses": sum(v == 0 for v in per_lens_values),
        },
        "numeric_fields": numeric,
        "source_zbin_raw_counts": dict(sorted(zbin_hist.items())),
        "b_mode_mask_raw_counts_no_polarity_assumed": dict(sorted(bmode_hist.items())),
        "positive_weight_rows": positive_weight,
        "zero_weight_rows": zero_weight,
        "negative_weight_rows": negative_weight,
        "invalid_sky_coordinate_rows": bad_sky_coords,
        "nonpositive_lens_redshift_rows": nonphysical_lens_z,
        "ellipticity_norm_gt_one_rows": ellipticity_norm_gt_one,
        "ellipticity_norm_ge_two_rows": ellipticity_norm_ge_two,
        "sigma_e_outside_published_range_rows": sigma_e_outside_published_range,
        "published_shape_cut_counts_before_b_mode_and_radial_selection": calibration_cuts,
        "shape_cut_and_tomographic_counts_by_zbin_before_b_mode_and_radial_selection": dict(sorted(selected_by_zbin.items())),
        "credentials_saved": state.get("credentials_saved"),
        "photoz_pdf_rows_included": export_summary.get("photoz_pdf_rows_included"),
        "b_mode_polarity_decoded": False,
        "inputs": {
            "row_status_sha256": sha(rows_dir / "row_status.json"),
            "row_summary_sha256": sha(rows_dir / "row_summary.json"),
            "row_export_plan_sha256": sha(plan_dir / "row_export_plan.json"),
            "count_summary_sha256": sha(preflight / "cas_counts" / "count_summary.json"),
            "analysis_script_sha256": sha(Path(__file__).resolve()),
        },
        "batches": batch_results,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("run the 0.9 GB row audit on lageunha, not the syntax login node")
    plan_dir = args.plan_dir or DEFAULT_PLAN
    result = audit(args.rows, args.preflight, plan_dir)
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_lenses", "n_pair_rows", "total_csv_bytes",
        "source_zbin_raw_counts", "b_mode_mask_raw_counts_no_polarity_assumed",
        "positive_weight_rows", "zero_weight_rows", "negative_weight_rows",
        "invalid_sky_coordinate_rows", "ellipticity_norm_gt_one_rows",
        "ellipticity_norm_ge_two_rows", "sigma_e_outside_published_range_rows",
        "published_shape_cut_counts_before_b_mode_and_radial_selection",
        "science_measurement")}, sort_keys=True))


if __name__ == "__main__":
    main()
