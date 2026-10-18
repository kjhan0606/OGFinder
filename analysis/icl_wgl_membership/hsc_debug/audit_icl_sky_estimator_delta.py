#!/usr/bin/env python3
"""Paired spatial-block comparison of ICL–satellite association sky estimators.

Preregistered comparison: retain the 74 existing matched pairs and their
original filter assignments; within each filter, compare the z/richness-
adjusted partial correlation using pixel-SFD ICL against plane-first pixel-SFD
ICL. The matched satellite measurement and every spatial-block bootstrap draw
are held fixed between estimators. The primary result is the i2 difference
``r(plane-first) - r(pixel-SFD)``. This is an exploratory sensitivity audit,
not a confirmatory science measurement.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels


DEFAULT_INPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv"
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def partial_pearson(x, y, z, richness):
    arrays = [np.asarray(values, dtype=float) for values in (x, y, z, richness)]
    if len(arrays[0]) < 5 or any(len(values) != len(arrays[0]) for values in arrays):
        return np.nan
    if any(not np.all(np.isfinite(values)) for values in arrays):
        return np.nan
    design = np.column_stack((np.ones(len(arrays[0])), arrays[2], arrays[3]))
    if np.linalg.matrix_rank(design) < 3:
        return np.nan
    rx = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    ry = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    denom = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / denom) if denom > 0 else np.nan


def percentile(values, q):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.percentile(finite, q)) if finite.size else None


def paired_block_contrast(rows, labels, replicas=5000, seed=20261020):
    """Compare two ICL estimates on identical rows and identical block draws."""
    labels = np.asarray(labels, dtype=str)
    if len(rows) != len(labels):
        raise ValueError("row and block-label counts differ")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("at least six connected spatial blocks are required")
    per_block = {block: np.flatnonzero(labels == block) for block in blocks}
    pixel = np.asarray([r["icl_delta_pixel_sfd_njy_arcsec2"] for r in rows], dtype=float)
    plane = np.asarray([r["icl_delta_pixel_sfd_after_plane_njy_arcsec2"] for r in rows], dtype=float)
    satellite = np.asarray([r["satellite_delta_center_sfd_njy_arcsec2"] for r in rows], dtype=float)
    z = np.asarray([r["z"] for r in rows], dtype=float)
    richness = np.asarray([r["richness"] for r in rows], dtype=float)
    arrays = (pixel, plane, satellite, z, richness)
    if any(not np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("non-finite values in the fixed paired sample")

    r_pixel = partial_pearson(pixel, satellite, z, richness)
    r_plane = partial_pearson(plane, satellite, z, richness)
    draws_pixel = np.full(replicas, np.nan)
    draws_plane = np.full(replicas, np.nan)
    rng = np.random.default_rng(seed)
    for replica in range(replicas):
        selected_blocks = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([per_block[block] for block in selected_blocks])
        draws_pixel[replica] = partial_pearson(
            pixel[indices], satellite[indices], z[indices], richness[indices])
        draws_plane[replica] = partial_pearson(
            plane[indices], satellite[indices], z[indices], richness[indices])
    difference = draws_plane - draws_pixel
    return {
        "n_pairs": len(rows),
        "n_connected_5deg_blocks": int(len(blocks)),
        "largest_block_pairs": int(max(len(index) for index in per_block.values())),
        "pixel_sfd_partial_r": r_pixel,
        "plane_first_partial_r": r_plane,
        "point_delta_plane_minus_pixel": r_plane - r_pixel,
        "paired_block_bootstrap_delta_ci95": [percentile(difference, 2.5),
                                                percentile(difference, 97.5)],
        "pixel_sfd_r_ci95": [percentile(draws_pixel, 2.5), percentile(draws_pixel, 97.5)],
        "plane_first_r_ci95": [percentile(draws_plane, 2.5), percentile(draws_plane, 97.5)],
        "valid_paired_replicas": int(np.isfinite(difference).sum()),
        "bootstrap_replicas": int(replicas),
        "seed": int(seed),
    }


def load_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    expected = {
        "cluster_id", "filter", "icl_delta_pixel_sfd_njy_arcsec2",
        "icl_delta_pixel_sfd_after_plane_njy_arcsec2",
        "satellite_delta_center_sfd_njy_arcsec2", "z", "richness",
        "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg",
    }
    if not rows or not expected.issubset(rows[0]):
        raise ValueError("input table is empty or lacks required columns")
    if len(rows) != 74:
        raise ValueError(f"expected the fixed 74-pair sample, received {len(rows)}")
    keys = [(row["cluster_id"], row["filter"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate cluster/filter rows in fixed sample")
    return rows


def run(args):
    source = args.input.resolve()
    output = args.work.resolve()
    if output.exists():
        raise ValueError(f"refusing to overwrite existing output: {output}")
    if args.replicas < 1000:
        raise ValueError("at least 1000 block-bootstrap replicas are required")
    rows = load_rows(source)
    by_filter = {}
    for band in ("i2", "i"):
        subset = [row for row in rows if row["filter"] == band]
        if not subset:
            raise ValueError(f"expected filter {band} in fixed sample")
        target_coords = [{"ra_deg": float(row["target_ra_deg"]),
                          "dec_deg": float(row["target_dec_deg"])}
                         for row in subset]
        control_coords = [{"ra_deg": float(row["control_ra_deg"]),
                           "dec_deg": float(row["control_dec_deg"])}
                          for row in subset]
        _, labels, _, _, _, _ = pair_block_labels(
            target_coords, control_coords, np.ones(len(subset), dtype=bool))
        by_filter[band] = paired_block_contrast(
            subset, labels, replicas=args.replicas,
            seed=args.seed + (0 if band == "i2" else 1))

    output.mkdir(parents=True, exist_ok=False)
    summary = {
        "schema": "hsc-icl-satellite-sky-estimator-contrast-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "preregistered_primary": "i2: paired spatial-block contrast r(plane-first ICL) - r(pixel-SFD ICL)",
        "method": {
            "fixed_sample": "all 74 matched pairs in the input table; no re-selection",
            "association": "partial Pearson correlation residualized on redshift and richness",
            "satellite_estimator": "the same center-SFD-corrected satellite delta for both ICL estimates",
            "resampling": "resample connected joint target/control 5-degree sky blocks, using identical block draws for both estimators",
            "primary_filter": "i2",
            "secondary_filter": "i, descriptive",
            "replicas": args.replicas,
            "seed": args.seed,
        },
        "n_pairs_total": len(rows),
        "filter_counts": {band: sum(row["filter"] == band for row in rows) for band in ("i", "i2")},
        "contrasts_by_filter": by_filter,
        "limitations": [
            "This quantifies sensitivity of the selected exploratory association to the plane-first sky estimator; it is not a causal test or a confirmatory ICL–satellite detection.",
            "It does not include upstream sky-transfer, PSF-wing, mask-leakage, selection, injection/recovery, or full image-level covariance uncertainty.",
            "The bootstrap intervals condition on the existing 74 matched pairs and chosen connected 5-degree blocks.",
        ],
        "input_sha256": sha256(source),
        "program_sha256": sha256(__file__),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "output": str(output), "contrasts_by_filter": by_filter},
                     ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
