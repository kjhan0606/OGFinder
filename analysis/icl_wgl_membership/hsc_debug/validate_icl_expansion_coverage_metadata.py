#!/usr/bin/env python3
"""Validate and summarize approved HSC PDR3 patch-coverage metadata results.

This produces a coarse 9-point-per-cluster coverage screen only. It cannot
establish usable full-aperture image coverage or scientific image quality.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader, sha
from submit_icl_expansion_coverage_metadata import (
    DEFAULT_PLAN, OUTPUT_COLUMNS, validate_plan, validate_result,
)


def as_float(value):
    if value in (None, "", "NULL", r"\N"):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Non-finite metadata value")
    return result


def aggregate_candidates(candidates, probes, rows_by_probe):
    probes_by_candidate = {}
    for probe in probes:
        probes_by_candidate.setdefault(str(probe["candidate_id"]), []).append(probe)
    result = []
    for candidate in candidates:
        # The original S20A pool stores its key as ``id``; the release-specific
        # S23B plan uses ``candidate_id`` to avoid cross-release ID collisions.
        candidate_id = str(candidate.get("id", candidate.get("candidate_id", "")))
        if not candidate_id:
            raise ValueError("Candidate inventory has no release-qualified ID")
        candidate_probes = probes_by_candidate.get(candidate_id, [])
        if len(candidate_probes) != 9 or {row["probe_id"] for row in candidate_probes} != {
            "center", *(f"r1Mpc_pa{angle:03d}" for angle in range(0, 360, 45))
        }:
            raise ValueError(f"Candidate {candidate_id} does not have the expected nine probes")
        covered = 0
        depths, seeings = [], []
        metadata_rows = 0
        imeas_true = idet_true = 0
        for probe in candidate_probes:
            key = (candidate_id, str(probe["probe_id"]))
            matched = rows_by_probe.get(key, [])
            if matched:
                covered += 1
            metadata_rows += len(matched)
            valid_depths = [value for row in matched
                            if (value := as_float(row["imag_psf_depth"])) is not None]
            valid_seeings = [value for row in matched
                             if (value := as_float(row["iseeing"])) is not None]
            if valid_depths:
                depths.append(max(valid_depths))
            if valid_seeings:
                seeings.append(min(valid_seeings))
            imeas_true += sum(str(row["imeas_exists"]).lower() in {"true", "t", "1"} for row in matched)
            idet_true += sum(str(row["idet_bkgd_exists"]).lower() in {"true", "t", "1"} for row in matched)
        result.append({
            "candidate_id": candidate_id,
            "ra_deg": candidate["ra_deg"], "dec_deg": candidate["dec_deg"],
            "z_cl": candidate["z_cl"], "richness": candidate["richness"],
            "coverage_probes": covered, "total_probes": 9,
            "probe_coverage_fraction": covered / 9.0,
            "all_nine_probes_covered": covered == 9,
            "metadata_rows": metadata_rows,
            "median_best_imag_psf_depth": statistics.median(depths) if depths else None,
            "median_best_iseeing": statistics.median(seeings) if seeings else None,
            "imeas_exists_true_fraction": imeas_true / metadata_rows if metadata_rows else None,
            "idet_bkgd_exists_true_fraction": idet_true / metadata_rows if metadata_rows else None,
            "coverage_is_nine_point_screen_only": True,
            "science_measurement": False,
        })
    return result


def load_tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_and_summarize(plan_dir: Path, output_dir: Path):
    manifest, batches, allowed_probes, plan_sha = validate_plan(plan_dir)
    results_dir = plan_dir / "cas_coverage_metadata"
    checkpoint_path = results_dir / "coverage_status.json"
    summary_path = results_dir / "coverage_summary.json"
    if not checkpoint_path.is_file() or not summary_path.is_file():
        raise ValueError("Coverage metadata run has no terminal checkpoint and summary")
    state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    saved_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if state.get("status") != "complete" or state.get("plan_sha256") != plan_sha:
        raise ValueError("Coverage metadata run is incomplete or belongs to another plan")
    scope_flags = ("images_returned", "object_rows_returned", "shape_rows_returned",
                   "photoz_rows_returned", "science_measurement", "credentials_saved")
    if any(state.get(key) is not False for key in scope_flags):
        raise ValueError("Coverage checkpoint reports data or state outside approved scope")
    if saved_summary.get("status") != "complete" or saved_summary.get("plan_sha256") != plan_sha:
        raise ValueError("Coverage summary is incomplete or mismatched")
    if saved_summary.get("output_columns") != OUTPUT_COLUMNS or any(
        saved_summary.get(key) is not False for key in scope_flags[:-1]
    ):
        raise ValueError("Coverage summary reports fields outside approved scope")
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {output_dir}")

    rows_by_probe = {}
    actual_total = 0
    for batch in batches:
        name = Path(batch["path"]).stem
        item = state.get("jobs", {}).get(name, {})
        csv_path = results_dir / f"{name}.csv"
        if item.get("status") != "complete" or not csv_path.is_file() or sha(csv_path) != item.get("sha256"):
            raise ValueError(f"Coverage metadata batch missing or changed: {name}")
        n_rows = validate_result(csv_path, batch, allowed_probes)
        if n_rows != item.get("rows") or csv_path.stat().st_size != item.get("bytes"):
            raise ValueError(f"Coverage metadata checkpoint disagrees with {name}")
        with csv_path.open(encoding="utf-8-sig", newline="") as stream:
            reader = catalogue_reader(stream)
            if reader.fieldnames != OUTPUT_COLUMNS:
                raise ValueError(f"Unexpected columns in {name}")
            for row in reader:
                rows_by_probe.setdefault((str(row["candidate_id"]), str(row["probe_id"])), []).append(row)
        actual_total += n_rows
    if actual_total != state.get("total_metadata_rows") or actual_total != saved_summary.get("n_metadata_rows"):
        raise ValueError("Coverage metadata row totals disagree")
    candidates = load_tsv(plan_dir / "candidate_centers.tsv")
    probes = load_tsv(plan_dir / "coverage_probes.tsv")
    screen = aggregate_candidates(candidates, probes, rows_by_probe)
    output_dir.mkdir(parents=True, exist_ok=False)
    csv_path = output_dir / "candidate_coverage_screen.tsv"
    fields = list(screen[0]) if screen else []
    with csv_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(screen)
    fully = sum(row["all_nine_probes_covered"] for row in screen)
    summary = {
        "schema": "astrafek-icl-pdr3-coverage-screen-summary-v1",
        "status": "metadata_screen_complete",
        "n_candidates": len(screen), "n_candidates_all_nine_probes": fully,
        "fraction_all_nine_probes": fully / len(screen) if screen else None,
        "n_metadata_rows": actual_total,
        "probe_geometry": "candidate center plus eight directions at proper 1,000 kpc",
        "limitations": [
            "nine-point patch screen only; not proof of full 1-Mpc disk coverage",
            "patch flags/depth/seeing are metadata, not pixel-mask or image-quality validation",
            "does not establish matched-control coverage or scientific ICL usability",
        ],
        "images_returned": False, "object_rows_returned": False,
        "shape_rows_returned": False, "photoz_rows_returned": False,
        "science_measurement": False,
        "input_plan_sha256": plan_sha,
        "candidate_screen_sha256": sha(csv_path),
    }
    (output_dir / "coverage_screen_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    plan_dir = args.plan.resolve()
    output = args.output_dir or (plan_dir / "coverage_screen")
    result = validate_and_summarize(plan_dir, output)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
