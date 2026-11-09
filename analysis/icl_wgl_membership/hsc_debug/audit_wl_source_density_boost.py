#!/usr/bin/env python3
"""Compute a count-only source-density boost diagnostic from saved WL stacks.

This is a point-only audit of the existing matched cluster/random shape stacks.
It does not reread shape rows, images, or query HSC, and it does not apply a
boost correction to DeltaSigma.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUTS = {
    "5deg": ROOT / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json",
    "2p5deg": ROOT / "wl_random_shape_stack_lens_reuse_2026-11-04_02p5/summary.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_stack(path: Path) -> dict:
    result = json.loads(path.read_text(encoding="utf-8"))
    if (result.get("status") != "exploratory" or result.get("science_measurement") is not False
            or result.get("random_lens_subtraction_applied") is not True):
        raise ValueError(f"unexpected or incomplete matched WL stack: {path}")
    sample = result.get("input", {}).get("matched_sample", {})
    n_lens = int(sample.get("n_common_shape_support_lenses", 0))
    n_random = int(sample.get("n_random_centers_in_common_support", 0))
    n_per_lens = int(sample.get("n_random_centers_per_common_lens", 0))
    if not n_lens or n_random != 5 * n_lens or n_per_lens != 5:
        raise ValueError(f"expected five matched random centers per common lens in {path}")
    if result.get("paired_spatial_jackknife", {}).get("science_measurement") is not False:
        raise ValueError(f"missing exploratory paired spatial jackknife: {path}")
    return result


def summarize(results: dict[str, dict]) -> list[dict]:
    reference = results["2p5deg"]
    ref_sample = reference["input"]["matched_sample"]
    n_lens = int(ref_sample["n_common_shape_support_lenses"])
    n_random = int(ref_sample["n_random_centers_in_common_support"])
    ref_profile = reference["paired_spatial_jackknife"]["profile"]
    output = []
    for block, result in results.items():
        sample = result["input"]["matched_sample"]
        if (sample["n_common_shape_support_lenses"] != n_lens
                or sample["n_random_centers_in_common_support"] != n_random):
            raise ValueError("block-scale summaries use different matched lens/random-center counts")
        profile = result["paired_spatial_jackknife"]["profile"]
        if len(profile) != len(ref_profile):
            raise ValueError("block-scale summaries have different radial-bin counts")
        for row, reference_row in zip(profile, ref_profile):
            if any(row[key] != reference_row[key] for key in
                   ("bin", "r_min_mpc", "r_max_mpc", "n_lens_pairs", "n_random_pairs")):
                raise ValueError("block-scale summaries disagree on radial pair counts")
    for row in ref_profile:
        n_lens_pairs = int(row["n_lens_pairs"])
        n_random_pairs = int(row["n_random_pairs"])
        if n_lens_pairs <= 0 or n_random_pairs <= 0:
            raise ValueError("cannot compute a source-density ratio from an empty radial bin")
        per_lens = n_lens_pairs / n_lens
        per_random_center = n_random_pairs / n_random
        boost = per_lens / per_random_center
        output.append({
            "bin": int(row["bin"]), "r_min_mpc": float(row["r_min_mpc"]),
            "r_max_mpc": float(row["r_max_mpc"]),
            "n_common_lenses": n_lens, "n_random_centers": n_random,
            "random_centers_per_lens": n_random / n_lens,
            "lens_source_pairs": n_lens_pairs,
            "random_source_pairs": n_random_pairs,
            "source_pairs_per_lens": per_lens,
            "source_pairs_per_random_center": per_random_center,
            "count_based_boost_candidate": boost,
            "fractional_source_density_excess": boost - 1.0,
            "uncertainty": "not estimated; pair-count totals do not retain radial-by-block counts",
            "science_measurement": False,
        })
    return output


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    results = {name: read_stack(path) for name, path in INPUTS.items()}
    profile = summarize(results)
    output.mkdir(parents=True, exist_ok=False)
    profile_path = output / "source_density_boost_profile.tsv"
    with profile_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(profile[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(profile)
    summary = {
        "schema": "astrafek-wl-source-density-boost-count-audit-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "boost_correction_applied": False,
        "method": "(cluster-lens eligible source-pair count / common lens count) / (matched-random eligible source-pair count / matched random-center count)",
        "source_selection": "same saved S19A shape, HSC Y3 tomographic n(z), B-mode, magnitude, resolution, and positive-efficiency gates used by the matched stacks",
        "normalization": "435 common lenses and all 2,175 planned matched random centers (five per lens); centers with zero exported source rows remain in the random denominator",
        "profile": profile,
        "input_sha256": {name: {"path": str(path.resolve()), "sha256": sha256(path)}
                         for name, path in INPUTS.items()},
        "hsc_network_used": False,
        "shape_rows_reread": False,
        "image_pixels_read": False,
        "limitations": [
            "The count ratio is a source-density boost candidate, not a validated correction to the random-subtracted DeltaSigma profile.",
            "It assumes comparable radial aperture area and usable depth/mask/PSF/sky around cluster and matched random centers; upstream matching is not fully validated.",
            "No radial-by-block source counts are stored in the stack summary, so this product has no uncertainty interval; spatially paired block counts must be recomputed from the approved local row exports.",
            "The count includes only sources passing the existing shape, B-mode, tomographic positive-efficiency, and radial cuts; photo-z dilution and calibration uncertainty remain.",
        ],
        "output_sha256": {str(profile_path.resolve()): sha256(profile_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps({
        "status": result["status"], "n_bins": len(result["profile"]),
        "science_measurement": result["science_measurement"],
        "boost_candidate": [{"radius_mpc": [r["r_min_mpc"], r["r_max_mpc"]],
                             "boost": r["count_based_boost_candidate"]}
                            for r in result["profile"]],
        "output_sha256": result["output_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
