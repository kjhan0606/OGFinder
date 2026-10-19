#!/usr/bin/env python3
"""Map the undocumented HSC Y3 b_mode_mask values in the approved WL export.

This diagnostic checks whether either stored value is spatially localized. It
does not infer the flag polarity from the shear signal, and does not perform a
science measurement or any network access.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from collections import defaultdict


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ROWS = DEFAULT_ROOT / "wl_shape_rows_export_2026-10-18_03"
DEFAULT_PREFLIGHT = DEFAULT_ROOT / "wl_shape_preflight_2026-10-18_01"
DEFAULT_PLAN = DEFAULT_ROOT / "wl_shape_rows_plan_2026-10-18_04"
EXPECTED_COLUMNS = [
    "lens_id", "lens_ra_deg", "lens_dec_deg", "lens_z", "object_id",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution",
    "psf_shape11", "psf_shape22", "psf_shape12", "b_mode_mask_raw",
    "source_zbin",
]
QUALITY_COLUMNS = ("e1", "e2", "sigma_e", "resolution", "mag_a10", "source_zbin")
PUBLISHED_BMODE_FOOTPRINT = {"ra_deg": [132.5, 140.0], "dec_deg": [1.6, 5.0]}


def normalize_shape_header(fieldnames: list[str] | None) -> list[str] | None:
    """Remove the comment marker used on the export's first CSV header field."""
    if fieldnames and fieldnames[0].startswith("# "):
        fieldnames = list(fieldnames)
        fieldnames[0] = fieldnames[0][2:]
    return fieldnames


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def grid_cell(ra_deg: float, dec_deg: float, width_deg: float = 1.0) -> tuple[int, int]:
    """Return a deterministic RA/Dec grid cell after validating sky bounds."""
    ra, dec, width = float(ra_deg), float(dec_deg), float(width_deg)
    if not (math.isfinite(ra) and 0 <= ra < 360):
        raise ValueError("RA must be finite and lie in [0, 360)")
    if not (math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError("Dec must be finite and lie in [-90, 90]")
    if not math.isfinite(width) or width <= 0 or width > 180:
        raise ValueError("grid width must lie in (0, 180]")
    return int(math.floor(ra / width)), int(math.floor((dec + 90.0) / width))


def in_published_bmode_footprint(ra_deg: float, dec_deg: float) -> bool:
    """Test the published GAMA09H B-mode rectangle using inclusive bounds."""
    ra, dec = float(ra_deg), float(dec_deg)
    if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError("invalid sky position")
    ra_min, ra_max = PUBLISHED_BMODE_FOOTPRINT["ra_deg"]
    dec_min, dec_max = PUBLISHED_BMODE_FOOTPRINT["dec_deg"]
    return ra_min <= ra <= ra_max and dec_min <= dec <= dec_max


def decode_mask_from_footprint(footprint_pairs: dict[str, list[int]]) -> int | None:
    """Decode the excluded raw bit only for an exact two-way footprint match."""
    inside = footprint_pairs.get("inside", [])
    outside = footprint_pairs.get("outside", [])
    if len(inside) != 2 or len(outside) != 2:
        raise ValueError("expected two raw mask counts inside and outside the published footprint")
    if inside[0] > 0 and inside[1] == 0 and outside[0] == 0 and outside[1] > 0:
        return 0
    if inside[1] > 0 and inside[0] == 0 and outside[1] == 0 and outside[0] > 0:
        return 1
    return None


def passes_published_source_cuts(row: dict[str, str]) -> bool:
    """The catalog-paper source cuts and released Y3 tomographic bins."""
    e1, e2, sigma = (float(row[name]) for name in QUALITY_COLUMNS[:3])
    resolution = float(row["resolution"])
    magnitude = float(row["mag_a10"])
    zbin = int(row["source_zbin"])
    values = (e1, e2, sigma, resolution, magnitude)
    if not all(math.isfinite(value) for value in values):
        return False
    return (math.hypot(e1, e2) < 2.0 and 0.0 <= sigma <= 0.4
            and resolution >= 0.3 and magnitude <= 25.5 and 1 <= zbin <= 4)


def audit(rows_dir: Path, preflight: Path, plan_dir: Path, grid_width_deg: float = 1.0) -> dict:
    rows_dir, preflight, plan_dir = map(Path, (rows_dir, preflight, plan_dir))
    export_plan = json.loads((plan_dir / "row_export_plan.json").read_text())
    checkpoint = json.loads((rows_dir / "row_status.json").read_text())
    summary = json.loads((rows_dir / "row_summary.json").read_text())
    count_summary = json.loads((preflight / "cas_counts" / "count_summary.json").read_text())
    plan_sha = sha256(plan_dir / "row_export_plan.json")
    if checkpoint.get("status") != "complete" or checkpoint.get("plan_sha256") != plan_sha:
        raise ValueError("shape-row export is incomplete or belongs to another plan")
    expected_total = int(count_summary.get("sum_lens_cone_shape_counts", -1))
    if (summary.get("status") != "complete" or summary.get("n_pair_rows") != expected_total
            or summary.get("count_summary_sum") != expected_total):
        raise ValueError("expected the complete 3,069,544-pair exploratory export")
    jobs = checkpoint.get("jobs", {})
    batches = export_plan.get("batches", [])
    if len(batches) != 32 or len(jobs) != 32:
        raise ValueError("expected all 32 checksum-verified shape-row batches")
    for key in ("science_measurement", "photoz_pdf_rows_included"):
        if summary.get(key) is not False:
            raise ValueError(f"unexpected export status for {key}")

    pair_counts = [0, 0]
    unique_counts = [0, 0]
    eligible_pairs = 0
    invalid_or_ineligible = 0
    # Each object should retain one sky position and one catalog mask value
    # across the multiple lens cones in which it can appear.
    seen_sources: dict[str, tuple[tuple[int, int], int, bool]] = {}
    cells = defaultdict(lambda: [0, 0])
    footprint_pairs = {"inside": [0, 0], "outside": [0, 0]}
    footprint_unique_sources = {"inside": [0, 0], "outside": [0, 0]}
    batch_hashes = {}
    observed_total = 0
    for index, batch in enumerate(batches):
        name = f"shape_rows_{index:03d}"
        record = jobs.get(name, {})
        path = rows_dir / f"{name}.csv"
        digest = sha256(path)
        if (record.get("status") != "complete" or record.get("sha256") != digest
                or record.get("rows") != batch.get("expected_pair_rows_from_count_only")):
            raise ValueError(f"shape-row batch hash/status mismatch: {name}")
        batch_hashes[name] = digest
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            reader.fieldnames = normalize_shape_header(reader.fieldnames)
            if reader.fieldnames != EXPECTED_COLUMNS:
                raise ValueError(f"unexpected shape-row columns in {name}")
            for row in reader:
                observed_total += 1
                raw = int(row["b_mode_mask_raw"])
                if raw not in (0, 1):
                    raise ValueError("b_mode_mask is not binary")
                if not passes_published_source_cuts(row):
                    invalid_or_ineligible += 1
                    continue
                ra, dec = float(row["source_ra_deg"]), float(row["source_dec_deg"])
                cell = grid_cell(ra, dec, grid_width_deg)
                eligible_pairs += 1
                pair_counts[raw] += 1
                cells[cell][raw] += 1
                inside_footprint = in_published_bmode_footprint(ra, dec)
                region = "inside" if inside_footprint else "outside"
                footprint_pairs[region][raw] += 1
                ident = str(row["object_id"])
                prior = seen_sources.get(ident)
                if prior is None:
                    seen_sources[ident] = (cell, raw, inside_footprint)
                    unique_counts[raw] += 1
                    footprint_unique_sources[region][raw] += 1
                elif prior != (cell, raw, inside_footprint):
                    raise ValueError(f"object_id {ident} changes sky cell or b_mode_mask")

    total_rows = int(checkpoint.get("total_pair_rows", -1))
    if total_rows != expected_total or observed_total != expected_total or eligible_pairs > total_rows:
        raise ValueError("observed row counts do not match the completed export")
    grid = []
    for (ra_bin, dec_bin), values in sorted(cells.items()):
        grid.append({
            "ra_min_deg": ra_bin * grid_width_deg,
            "ra_max_deg": min(360.0, (ra_bin + 1) * grid_width_deg),
            "dec_min_deg": dec_bin * grid_width_deg - 90.0,
            "dec_max_deg": min(90.0, (dec_bin + 1) * grid_width_deg - 90.0),
            "raw_0_pair_rows": values[0],
            "raw_1_pair_rows": values[1],
            "raw_0_fraction": values[0] / sum(values),
        })
    masked_raw_value = decode_mask_from_footprint(footprint_pairs)
    return {
        "schema": "hsc-s19a-bmode-spatial-distribution-audit-v1",
        "status": "complete_spatial_diagnostic",
        "science_measurement": False,
        "b_mode_polarity_decoded": masked_raw_value is not None,
        "b_mode_masked_raw_value": masked_raw_value,
        "b_mode_polarity_validation": "exact raw-bit match to the published GAMA09H exclusion rectangle; no shear statistic used",
        "raw_value_meaning_inferred_from_shear": False,
        "n_lenses": 800,
        "n_exported_pair_rows": total_rows,
        "n_pair_rows_after_published_shape_and_zbin_cuts": eligible_pairs,
        "n_pair_rows_excluded_by_those_cuts": invalid_or_ineligible,
        "raw_mask_pair_rows_after_cuts": {"0": pair_counts[0], "1": pair_counts[1]},
        "n_unique_sources_after_cuts": len(seen_sources),
        "raw_mask_unique_sources_after_cuts": {"0": unique_counts[0], "1": unique_counts[1]},
        "published_bmode_footprint": PUBLISHED_BMODE_FOOTPRINT,
        "mask_value_counts_by_published_footprint": {
            "pair_rows": {region: {"0": values[0], "1": values[1]}
                          for region, values in footprint_pairs.items()},
            "unique_sources": {region: {"0": values[0], "1": values[1]}
                               for region, values in footprint_unique_sources.items()},
        },
        "grid_width_deg": grid_width_deg,
        "source_sky_grid": grid,
        "provenance": {
            "row_plan_sha256": plan_sha,
            "count_summary_sha256": sha256(preflight / "cas_counts" / "count_summary.json"),
            "row_batch_sha256": batch_hashes,
            "grid_method": "equal-width RA/Dec bins; no area correction at high declination",
        },
        "limitations": [
            "raw-bit polarity is decoded only when each value maps exclusively inside or outside the published rectangular exclusion footprint",
            "the source map is a coordinate-bin diagnostic, not a spherical-area density estimator",
            "lens-source pairs repeat sources across overlapping lens apertures; unique-source counts are also reported",
            "the export covers only the approved 800-lens sample and cone queries",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--grid-width-deg", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    result = audit(args.rows, args.preflight, args.plan_dir, args.grid_width_deg)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "pairs_after_published_cuts": result["n_pair_rows_after_published_shape_and_zbin_cuts"],
        "unique_sources": result["n_unique_sources_after_cuts"],
        "raw_mask_pair_rows": result["raw_mask_pair_rows_after_cuts"],
        "polarity_decoded": False,
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
