#!/usr/bin/env python3
"""Independently verify output hashes, pair accounting, and science-state gates."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate(output: Path, *, rehash_inputs: bool = False) -> dict:
    output = Path(output).resolve()
    summary_path = output / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "astrafek-empirical-psf-outer-field-wings-v1":
        raise ValueError("unexpected outer-field profile schema")
    if summary.get("status") != "complete_exploratory_diagnostic":
        raise ValueError("outer-field result is not marked as a completed diagnostic")
    if summary.get("science_measurement") is not False or summary.get("psf_approved") is not False:
        raise ValueError("exploratory result crossed a science/PSF approval gate")
    for relative, expected in summary["output_sha256"].items():
        path = output / relative
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"output checksum mismatch: {relative}")
    input_hashes = summary["input_sha256"]
    missing = [name for name in input_hashes if not Path(name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing input files listed in summary: {missing[:3]}")
    verified_inputs = 0
    if rehash_inputs:
        for name, expected in input_hashes.items():
            if _sha(Path(name)) != expected:
                raise ValueError(f"input checksum mismatch: {name}")
            verified_inputs += 1

    candidates = _tsv(output / "selected_candidates.tsv")
    profiles = _tsv(output / "matched_radial_profiles.tsv")
    stack = _tsv(output / "frame_equal_paired_stack.tsv")
    candidate_ids = {row["candidate_id"] for row in candidates}
    if len(candidates) != 18 or len(candidate_ids) != len(candidates):
        raise ValueError("selected candidate table does not contain 18 unique stars")
    if len({row["cluster_id"] for row in candidates}) != len(candidates):
        raise ValueError("selected stars do not represent 18 independent HSC target frames")
    key_fields = ("candidate_id", "role", "sky_r_min_pix", "sky_r_max_pix", "r_min_kpc", "r_max_kpc")
    profile_keys = [tuple(row[name] for name in key_fields) for row in profiles]
    if len(profile_keys) != len(set(profile_keys)):
        raise ValueError("duplicate candidate/role/annulus rows in profile table")
    if len(profiles) != len(candidate_ids) * 2 * 2 * 12:
        raise ValueError("profile rows do not match candidates × roles × sky annuli × radial bins")
    for row in profiles:
        if row["candidate_id"] not in candidate_ids or row["role"] not in ("target", "control"):
            raise ValueError("profile row references an unknown candidate or role")
        if row["science_measurement"].lower() != "false" or row["psf_approved"].lower() != "false":
            raise ValueError("profile row crossed a science/PSF approval gate")
        value = row["profile_per_target_core_per_pixel"]
        if row["status"] == "ok":
            if value == "" or not math.isfinite(float(value)):
                raise ValueError("successful profile row has a missing or non-finite value")
        elif value != "":
            raise ValueError("blocked profile row unexpectedly contains a measurement")

    profiles_by_key: dict[tuple[float, float, float, float, str], dict[str, dict[str, str]]] = {}
    for row in profiles:
        key = (float(row["sky_r_min_pix"]), float(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]), row["candidate_id"])
        profiles_by_key.setdefault(key, {})[row["role"]] = row
    expected_stack: dict[tuple[float, float, float, float], list[float]] = {}
    for key, pair in profiles_by_key.items():
        if pair.get("target", {}).get("status") == "ok" and pair.get("control", {}).get("status") == "ok":
            s_key = key[:4]
            expected_stack.setdefault(s_key, []).append(
                float(pair["target"]["profile_per_target_core_per_pixel"])
                - float(pair["control"]["profile_per_target_core_per_pixel"])
            )
    if len(stack) != 24:
        raise ValueError("paired stack must have 2 sky choices × 12 physical bins")
    for row in stack:
        if row["science_measurement"].lower() != "false" or row["psf_approved"].lower() != "false":
            raise ValueError("stack row crossed a science/PSF approval gate")
        key = (float(row["sky_r_min_pix"]), float(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]))
        values = expected_stack.get(key, [])
        if int(row["n_paired_frames"]) != len(values):
            raise ValueError("paired frame count disagrees with matched profile rows")
        if values:
            observed = float(row["paired_difference_frame_equal_mean"])
            expected = sum(values) / len(values)
            if abs(observed - expected) > max(1e-14, 1e-10 * abs(expected)):
                raise ValueError("frame-equal paired mean does not recompute from the profile table")
        elif row["paired_difference_frame_equal_mean"] != "":
            raise ValueError("empty paired bin unexpectedly has a reported mean")
    max_pairs = max(int(row["n_paired_frames"]) for row in stack)
    result = {
        "schema": "astrafek-empirical-psf-outer-field-wings-validation-v1",
        "status": "validated_exploratory_output",
        "selected_candidates": len(candidates), "independent_target_frames": len({r["cluster_id"] for r in candidates}),
        "profile_rows": len(profiles), "paired_stack_bins": len(stack),
        "max_paired_frames_in_any_radial_bin": max_pairs,
        "output_hashes_verified": len(summary["output_sha256"]),
        "input_files_exist": len(input_hashes), "input_hashes_recomputed": verified_inputs,
        "science_measurement": False, "psf_approved": False,
        "interpretation": "The output and arithmetic are internally consistent; at most four paired frames in a radial bin are insufficient for an ICL science correction.",
    }
    validation_path = output / "validation.json"
    if validation_path.exists():
        if json.loads(validation_path.read_text(encoding="utf-8")) != result:
            raise ValueError(f"existing validation result differs from current checks: {validation_path}")
    else:
        validation_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["validation_sha256"] = _sha(validation_path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rehash-inputs", action="store_true", help="read all source images again to recalculate their hashes")
    args = parser.parse_args()
    print(json.dumps(validate(args.output, rehash_inputs=args.rehash_inputs), indent=2))


if __name__ == "__main__":
    main()
