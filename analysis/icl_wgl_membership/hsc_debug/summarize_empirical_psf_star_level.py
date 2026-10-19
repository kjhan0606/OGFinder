#!/usr/bin/env python3
"""Summarize outer empirical-PSF profiles by unique Gaia-supported candidate."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import time
from collections import Counter, defaultdict


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"q16": None, "median": None, "q84": None}
    ordered = sorted(values)

    def percentile(p: float) -> float:
        position = (len(ordered) - 1) * p
        lo = int(position)
        hi = min(lo + 1, len(ordered) - 1)
        return ordered[lo] + (ordered[hi] - ordered[lo]) * (position - lo)

    return {"q16": percentile(0.16), "median": percentile(0.5), "q84": percentile(0.84)}


def run(profiles_path: Path, candidates_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite star-level assessment: {output}")
    candidates = read_tsv(candidates_path)
    if len(candidates) != 480:
        raise ValueError(f"expected 480 Gaia-supported candidates, found {len(candidates)}")
    candidate_by_id = {row["candidate_id"]: row for row in candidates}
    if len(candidate_by_id) != len(candidates):
        raise ValueError("duplicate candidate IDs in input table")
    source_ids = [row["gaia_source_id"] for row in candidates]
    if len(set(source_ids)) != len(source_ids):
        raise ValueError("candidate table does not contain unique Gaia source IDs")

    selected = []
    for row in read_tsv(profiles_path):
        if float(row["r_min_pix"]) != 64.0 or float(row["r_max_pix"]) != 96.0:
            continue
        candidate_id = row["candidate_id"]
        if candidate_id not in candidate_by_id:
            raise ValueError(f"unknown candidate ID in profiles: {candidate_id}")
        selected.append(row)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    keys = set()
    for row in selected:
        key = (row["candidate_id"], row["kron_factor"],
               row["background_r_min_pix"], row["background_r_max_pix"])
        if key in keys:
            raise ValueError(f"duplicate candidate/method outer profile: {key}")
        keys.add(key)
        grouped[row["candidate_id"]].append(row)

    per_star = []
    for candidate_id, candidate in candidate_by_id.items():
        rows = grouped.get(candidate_id, [])
        valid = [row for row in rows if row["radial_coverage_pass"].lower() == "true"]
        if any(row["radial_coverage_pass"].lower() not in ("true", "false") for row in rows):
            raise ValueError("invalid radial_coverage_pass value")
        values = [float(row["mean_per_core_flux"]) for row in valid]
        if any(not math.isfinite(value) for value in values):
            raise ValueError(f"non-finite outer profile for candidate {candidate_id}")
        signs = ["positive" if value > 0 else "nonpositive" for value in values]
        per_star.append({
            "candidate_id": candidate_id,
            "cluster_id": candidate["cluster_id"],
            "gaia_source_id": candidate["gaia_source_id"],
            "gaia_g_mag": candidate["gaia_g_mag"],
            "n_outer_methods_recorded": len(rows),
            "n_outer_methods_coverage_pass": len(valid),
            "n_positive_methods": signs.count("positive"),
            "n_nonpositive_methods": signs.count("nonpositive"),
            "outer_mean_per_core_flux_median": statistics.median(values) if values else "",
            "outer_sign_pattern": (
                "no_coverage" if not signs else
                "all_positive" if all(sign == "positive" for sign in signs) else
                "all_nonpositive" if all(sign == "nonpositive" for sign in signs) else
                "mixed"
            ),
            "all_four_outer_methods_coverage_pass": len(valid) == 4,
        })

    output.mkdir(parents=True)
    fields = list(per_star[0])
    with (output / "per_star_outer_profile.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(per_star)

    with_outer = [row for row in per_star if row["n_outer_methods_coverage_pass"] > 0]
    complete_outer = [row for row in per_star if row["all_four_outer_methods_coverage_pass"]]
    sign_counts = Counter(row["outer_sign_pattern"] for row in with_outer)
    per_source_medians = [float(row["outer_mean_per_core_flux_median"]) for row in with_outer]
    magnitude_bins = {"G_le_19": [], "G_19_to_20": [], "G_gt_20": []}
    for row in with_outer:
        magnitude = float(row["gaia_g_mag"])
        bucket = "G_le_19" if magnitude <= 19 else "G_19_to_20" if magnitude <= 20 else "G_gt_20"
        magnitude_bins[bucket].append(float(row["outer_mean_per_core_flux_median"]))

    scales = [float(row["pixel_scale_arcsec"]) for row in selected]
    summary = {
        "status": "complete",
        "outer_annulus_pixels": [64, 96],
        "outer_annulus_arcsec_range_across_frames": [
            round(64 * min(scales), 3), round(96 * max(scales), 3)
        ] if scales else None,
        "method_dimensions": {"kron_factors": [2.5, 3.5], "background_annuli_pixels": [[120, 150], [160, 190]]},
        "n_input_candidates": len(candidates),
        "n_candidates_with_outer_profile_record": len(grouped),
        "n_candidates_with_outer_coverage": len(with_outer),
        "n_candidates_with_all_four_outer_methods_covered": len(complete_outer),
        "sign_pattern_counts_by_candidate": dict(sorted(sign_counts.items())),
        "unique_candidate_median_profile_per_core_flux": quantiles(per_source_medians),
        "magnitude_bin_candidate_medians": {
            key: {"n_candidates": len(values), **quantiles(values)}
            for key, values in magnitude_bins.items()
        },
        "n_outer_method_profiles_recorded": len(selected),
        "n_outer_method_profiles_coverage_pass": sum(
            row["radial_coverage_pass"].lower() == "true" for row in selected
        ),
        "n_outer_method_profiles_positive": sum(
            row["radial_coverage_pass"].lower() == "true"
            and float(row["mean_per_core_flux"]) > 0 for row in selected
        ),
        "n_outer_method_profiles_nonpositive": sum(
            row["radial_coverage_pass"].lower() == "true"
            and float(row["mean_per_core_flux"]) <= 0 for row in selected
        ),
        "input_sha256": {
            str(profiles_path.resolve()): sha256(profiles_path),
            str(candidates_path.resolve()): sha256(candidates_path),
        },
        "output_sha256": sha256(output / "per_star_outer_profile.tsv"),
        "script_sha256": sha256(Path(__file__)),
        "interpretation": "Per-candidate aggregation avoids treating repeated background/Kron methods as independent stars; coverage is geometric support, not a wing detection.",
        "science_measurement": False,
        "psf_approved": False,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.profiles, args.candidates, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
