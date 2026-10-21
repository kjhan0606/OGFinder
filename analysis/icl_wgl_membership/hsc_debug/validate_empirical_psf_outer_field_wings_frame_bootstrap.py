#!/usr/bin/env python3
"""Independently verify the 600-kpc multi-star PSF pilot and outer-ring support."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


EXPECTED_FILES = ("selected_candidates.tsv", "matched_radial_profiles.tsv", "frame_bootstrap_paired_stack.tsv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _number(row, field, *, optional=False):
    raw = row.get(field, "")
    if optional and raw == "":
        return None
    value = float(raw)
    if not np.isfinite(value):
        raise ValueError(f"nonfinite value in {field}")
    return value


def recompute_stack(profile_rows: list[dict[str, str]], *, n_bootstrap: int, seed: int) -> list[dict]:
    star_values = defaultdict(lambda: defaultdict(dict))
    all_keys = set()
    for row in profile_rows:
        key = (str(row["filter"]), float(row["sky_r_min_pix"]), float(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]))
        all_keys.add(key)
        if row["status"] != "ok" or row["profile_per_target_core_per_pixel"] == "":
            continue
        frame_key = (*key, str(row["cluster_id"]))
        star_values[frame_key][str(row["candidate_id"])][row["role"]] = _number(
            row, "profile_per_target_core_per_pixel")
    frame_values = defaultdict(dict)
    for key, candidates in star_values.items():
        paired = [values["target"] - values["control"] for values in candidates.values()
                  if "target" in values and "control" in values]
        if paired:
            frame_values[key[:5]][key[5]] = paired
    rng = np.random.default_rng(seed)
    output = []
    for key in sorted(all_keys):
        grouped = frame_values.get(key, {})
        frames = sorted(grouped)
        means = np.asarray([np.mean(grouped[frame]) for frame in frames], dtype=float)
        n_frames = len(frames)
        n_pairs = sum(len(grouped[frame]) for frame in frames)
        if n_frames and n_bootstrap:
            choices = rng.integers(0, n_frames, size=(n_bootstrap, n_frames))
            draws = np.mean(means[choices], axis=1)
            q16, q84 = map(float, np.quantile(draws, (0.16, 0.84)))
        else:
            q16 = q84 = None
        output.append({
            "filter": key[0], "sky_r_min_pix": key[1], "sky_r_max_pix": key[2],
            "r_min_kpc": key[3], "r_max_kpc": key[4], "r_mid_kpc": (key[3] + key[4]) / 2,
            "n_complete_star_pairs": n_pairs, "n_independent_target_frames": n_frames,
            "n_frames_with_multiple_complete_stars": sum(len(grouped[frame]) > 1 for frame in frames),
            "frame_equal_paired_difference_mean": float(np.mean(means)) if n_frames else None,
            "frame_equal_paired_difference_median": float(np.median(means)) if n_frames else None,
            "bootstrap_q16": q16, "bootstrap_q84": q84,
        })
    return output


def partial_support(profile_rows: list[dict[str, str]], *, area_cut: float, min_sectors: int) -> list[dict]:
    if not 0 < area_cut <= 1 or not 1 <= min_sectors <= 8:
        raise ValueError("invalid partial-support thresholds")
    grouped = defaultdict(dict)
    annuli = set()
    for row in profile_rows:
        key = (str(row["filter"]), float(row["sky_r_min_pix"]), float(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]))
        annuli.add(key)
        identity = (str(row["candidate_id"]), str(row["cluster_id"]))
        grouped[key].setdefault(identity, {})[row["role"]] = row
    output = []
    for key in sorted(annuli):
        if key[3] < 200 or key[4] > 300:
            continue
        complete = []
        for (candidate_id, frame), roles in grouped[key].items():
            if not {"target", "control"} <= roles.keys():
                continue
            accepted = True
            for role in ("target", "control"):
                row = roles[role]
                if (row.get("background", "") == "" or row.get("target_core_flux", "") == ""
                        or row.get("raw_mean", "") == ""
                        or float(row["target_core_flux"]) <= 0
                        or float(row["area_fraction"]) < area_cut
                        or int(row["n_covered_sectors"]) < min_sectors):
                    accepted = False
                    break
            if accepted:
                complete.append((candidate_id, frame))
        output.append({
            "filter": key[0], "sky_r_min_pix": key[1], "sky_r_max_pix": key[2],
            "r_min_kpc": key[3], "r_max_kpc": key[4],
            "area_fraction_min_per_role": area_cut,
            "minimum_covered_sectors_per_role": min_sectors,
            "n_complete_star_pairs": len(complete),
            "n_independent_frames": len({frame for _, frame in complete}),
        })
    return output


def validate(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    summary_path = output_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "astrafek-empirical-psf-outer-field-frame-bootstrap-v1":
        raise ValueError("unexpected frame-bootstrap summary schema")
    if summary.get("science_measurement") is not False or summary.get("psf_approved") is not False:
        raise ValueError("PSF pilot has an unexpected science/approval status")
    hashes = summary.get("output_sha256", {})
    if set(hashes) != set(EXPECTED_FILES):
        raise ValueError("summary output hash manifest is incomplete or unexpected")
    for name in EXPECTED_FILES:
        if sha256(output_dir / name) != hashes[name]:
            raise ValueError(f"output checksum mismatch: {name}")
    selected = read_tsv(output_dir / "selected_candidates.tsv")
    profiles = read_tsv(output_dir / "matched_radial_profiles.tsv")
    reported = read_tsv(output_dir / "frame_bootstrap_paired_stack.tsv")
    selected_ids = [row["candidate_id"] for row in selected]
    if len(selected_ids) != len(set(selected_ids)) or set(selected_ids) != set(summary.get("selected_candidate_ids", selected_ids)):
        raise ValueError("selected candidate IDs are duplicate or disagree with summary")
    if len(selected_ids) != int(summary["counts"]["selected"]):
        raise ValueError("selected candidate row count disagrees with summary")
    unique_frames = len({row["cluster_id"] for row in selected})
    if unique_frames != summary["n_unique_target_frames"]:
        raise ValueError("independent-frame count disagrees with selected candidates")
    if len(profiles) != int(summary["n_profile_rows"]) or len(profiles) != len(selected) * 48:
        raise ValueError("profile row count disagrees with 2 roles × 2 skies × 12 radii")
    if {row["candidate_id"] for row in profiles} != set(selected_ids):
        raise ValueError("profile candidate IDs disagree with selected candidate table")
    if len(reported) != 48:
        raise ValueError("expected 2 filters × 2 sky annuli × 12 physical radial bins in stack table")
    expected = recompute_stack(profiles, n_bootstrap=int(summary["n_bootstrap_draws"]), seed=int(summary["seed"]))
    if len(expected) != len(reported):
        raise ValueError("recomputed stack has a different number of cells")
    float_fields = ("r_mid_kpc", "frame_equal_paired_difference_mean",
                    "frame_equal_paired_difference_median", "bootstrap_q16", "bootstrap_q84")
    integer_fields = ("n_complete_star_pairs", "n_independent_target_frames", "n_frames_with_multiple_complete_stars")
    for expected_row, actual in zip(expected, reported):
        if expected_row["filter"] != actual["filter"]:
            raise ValueError("stack cell filter differs")
        for field in ("sky_r_min_pix", "sky_r_max_pix", "r_min_kpc", "r_max_kpc"):
            if not np.isclose(float(expected_row[field]), float(actual[field]), rtol=0, atol=1e-12):
                raise ValueError(f"stack cell key differs: {field}")
        for field in integer_fields:
            if int(actual[field]) != expected_row[field]:
                raise ValueError(f"recomputed stack differs: {field}")
        for field in float_fields:
            wanted = expected_row[field]
            got = None if actual[field] == "" else float(actual[field])
            if wanted is None and got is None:
                continue
            if wanted is None or got is None or not np.isclose(wanted, got, rtol=1e-12, atol=1e-18):
                raise ValueError(f"recomputed stack differs: {field}")
    partial = []
    for area_cut, sectors in ((0.2, 6), (0.2, 8), (0.3, 6), (0.5, 6), (0.5, 8)):
        partial.extend(partial_support(profiles, area_cut=area_cut, min_sectors=sectors))
    return {
        "schema": "astrafek-empirical-psf-frame-bootstrap-validation-v1",
        "status": "validated_exploratory_diagnostic",
        "science_measurement": False,
        "psf_approved": False,
        "summary_sha256": sha256(summary_path),
        "output_sha256_verified": hashes,
        "selected_star_pairs": len(selected),
        "independent_target_frames": unique_frames,
        "profile_rows_recomputed": len(profiles),
        "stack_cells_recomputed": len(reported),
        "bootstrap_replicates_reproduced": int(summary["n_bootstrap_draws"]),
        "partial_mask_support_sensitivity_200_300_kpc": partial,
        "limitations": [
            "Partial-support thresholds are diagnostics only; they do not correct mask-induced selection bias.",
            "No individual 200-300 kpc profile should be interpreted as a calibrated PSF wing.",
            "science_measurement remains false and psf_approved remains false.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--validation-json", type=Path)
    args = parser.parse_args()
    result = validate(args.output_dir)
    out = args.validation_json.resolve() if args.validation_json else args.output_dir.resolve() / "validation.json"
    if out.exists():
        raise FileExistsError(f"refusing to overwrite validation output: {out}")
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["validation_sha256"] = sha256(out)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
