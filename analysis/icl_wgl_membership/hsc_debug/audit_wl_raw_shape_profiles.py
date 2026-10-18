#!/usr/bin/env python3
"""Descriptive raw HSC tangential/cross ellipticity profiles.

This is a pre-calibration diagnostic, not a shear or DeltaSigma measurement.
It applies the published HSC source cuts, but intentionally does not apply
responsivity, multiplicative/additive or selection-bias corrections, source
redshift efficiency, random-lens subtraction, or covariance. Since the stored
``b_mode_mask`` polarity is undocumented, it reports both possible exclusions
and an all-rows diagnostic without assigning either polarity as correct.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_wl_shape_stack import DEFAULT_PLAN, DEFAULT_PREFLIGHT, verify_and_load_pairs
from s19a_shape_stack import _to_tangential, pair_geometry
from ogfkit.wlcluster import Cosmology

RADIAL_EDGES_MPC = np.array([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0], dtype=float)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def published_source_selection(pairs):
    """Return the HSC Y3 source cuts and the four calibrated tomographic bins."""
    e = np.hypot(np.asarray(pairs["e1"], float), np.asarray(pairs["e2"], float))
    sigma_e = np.asarray(pairs["sigma_e"], float)
    resolution = np.asarray(pairs["resolution"], float)
    magnitude = np.asarray(pairs["mag_a10"], float)
    zbin = np.asarray(pairs["source_zbin"], int)
    quality = ((e < 2.0) & (sigma_e >= 0.0) & (sigma_e <= 0.4)
               & (resolution >= 0.3) & (magnitude <= 25.5))
    return quality, zbin


def radial_bin_index(radius, edges=RADIAL_EDGES_MPC):
    """Assign [lo, hi) bins, including the final upper edge in the last bin."""
    radius = np.asarray(radius, dtype=float)
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)):
        raise ValueError("radial edges must be a finite one-dimensional vector")
    if np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError("radial edges must be positive and strictly increasing")
    index = np.searchsorted(edges, radius, side="right") - 1
    index = np.where(radius == edges[-1], len(edges) - 2, index)
    valid = np.isfinite(radius) & (radius >= edges[0]) & (radius <= edges[-1])
    return np.where(valid, index, -1)


def summarize_profiles(pairs, *, h0=70.0, omega_m=0.3,
                       radial_edges_mpc=RADIAL_EDGES_MPC):
    """Calculate raw shape profiles for both B-mode polarity hypotheses."""
    edges = np.asarray(radial_edges_mpc, float)
    quality, zbins = published_source_selection(pairs)
    raw_mask = np.asarray(pairs["b_mode_mask_raw"], float)
    if not np.all(np.isfinite(raw_mask)) or np.any((raw_mask != 0) & (raw_mask != 1)):
        raise ValueError("raw b_mode_mask must contain only 0 and 1")

    cosmology = Cosmology(H0=h0, Om=omega_m)
    radius, phi = pair_geometry(pairs, cosmology)
    radial_bin = radial_bin_index(radius, edges)
    e_t, e_x = _to_tangential(phi, np.asarray(pairs["e1"], float),
                              np.asarray(pairs["e2"], float))
    weight = np.asarray(pairs["weight"], float)
    lens_ids = np.asarray(pairs["lens_id"], dtype=object)

    selected_z = quality & (zbins >= 1) & (zbins <= 4)
    bmode_options = {
        "no_b_mode_exclusion": np.ones(len(raw_mask), dtype=bool),
        "exclude_raw_zero_hypothesis": raw_mask != 0,
        "exclude_raw_one_hypothesis": raw_mask != 1,
    }
    profiles = []
    groups = [("all_zbins_1_to_4", selected_z)]
    groups.extend((f"source_zbin_{k}", selected_z & (zbins == k)) for k in range(1, 5))
    for mode, mask in bmode_options.items():
        for group_name, group_mask in groups:
            usable = group_mask & mask & (radial_bin >= 0)
            bins = radial_bin[usable]
            w = weight[usable]
            sw = np.bincount(bins, weights=w, minlength=len(edges) - 1)
            sw2 = np.bincount(bins, weights=w * w, minlength=len(edges) - 1)
            ste = np.bincount(bins, weights=w * e_t[usable], minlength=len(edges) - 1)
            sex = np.bincount(bins, weights=w * e_x[usable], minlength=len(edges) - 1)
            count = np.bincount(bins, minlength=len(edges) - 1)
            for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                in_bin = usable & (radial_bin == index)
                n_lenses = int(len(set(lens_ids[in_bin].tolist())))
                denom = float(sw[index])
                profiles.append({
                    "b_mode_scenario": mode,
                    "source_bin": group_name,
                    "radial_bin": index,
                    "r_min_mpc": float(lo),
                    "r_max_mpc": float(hi),
                    "n_lens_source_pairs": int(count[index]),
                    "n_lenses": n_lenses,
                    "sum_shape_weight": denom,
                    "effective_weighted_pairs": float(denom * denom / sw2[index]) if sw2[index] > 0 else None,
                    "mean_raw_e_t": float(ste[index] / denom) if denom > 0 else None,
                    "mean_raw_e_x": float(sex[index] / denom) if denom > 0 else None,
                })

    profile_map = {(row["b_mode_scenario"], row["source_bin"], row["radial_bin"]): row
                   for row in profiles}
    for group_name, _ in groups:
        for index in range(len(edges) - 1):
            all_count = profile_map[("no_b_mode_exclusion", group_name, index)]["n_lens_source_pairs"]
            zero_count = profile_map[("exclude_raw_zero_hypothesis", group_name, index)]["n_lens_source_pairs"]
            one_count = profile_map[("exclude_raw_one_hypothesis", group_name, index)]["n_lens_source_pairs"]
            if zero_count + one_count != all_count:
                raise AssertionError("opposite raw-mask exclusions do not partition the no-mask pair count")

    summary = {
        "schema": "hsc-s19a-raw-shape-profile-diagnostic-v1",
        "status": "complete_descriptive",
        "science_measurement": False,
        "calibrated_shear": False,
        "delta_sigma_calculated": False,
        "b_mode_polarity_decoded": False,
        "random_lens_subtraction_applied": False,
        "covariance_applied": False,
        "source_selection": "|e|<2, 0<=sigma_e<=0.4, resolution>=0.3, mag_A<=25.5, source_zbin 1..4",
        "tangential_component": "e_t=-e1*cos(2phi)-e2*sin(2phi); phi measured east of lens toward north",
        "weighting": "raw HSC shape-catalog derived weight; no responsivity or bias correction",
        "b_mode_scenarios": list(bmode_options),
        "radial_edges_mpc": [float(x) for x in edges],
        "input_pair_count": int(len(radius)),
        "published_source_cuts_and_zbins_count": int(np.count_nonzero(selected_z)),
        "source_bin_counts_before_radial_and_b_mode": {
            str(k): int(np.count_nonzero(selected_z & (zbins == k))) for k in range(1, 5)
        },
        "b_mode_raw_value_counts": {
            str(k): int(np.count_nonzero(raw_mask == k)) for k in (0, 1)
        },
        "b_mode_opposite_exclusions_partition_pair_counts": True,
        "cosmology": {"H0_km_s_mpc": float(h0), "Omega_m": float(omega_m)},
        "limitations": [
            "This is a descriptive raw ellipticity profile, not a calibrated shear estimator.",
            "No HSC responsivity, multiplicative/additive bias, or selection-bias correction is applied.",
            "No source n(z) or lensing-efficiency / foreground-dilution correction is applied.",
            "The two raw b_mode_mask exclusions are sensitivity scenarios, not a polarity decision.",
            "No random-lens subtraction, source-overlap correction, or spatial covariance is applied.",
            "No mass model or significance claim is produced.",
        ],
        "profiles": profiles,
    }
    return summary


def run(args):
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {out}")
    pairs, provenance = verify_and_load_pairs(
        args.rows.resolve(), args.preflight, args.plan_dir, None,
    )
    result = summarize_profiles(pairs, h0=args.h0, omega_m=args.omega_m,
                                radial_edges_mpc=args.radial_edges_mpc)
    provenance.pop("masked_value", None)
    provenance["row_summary_sha256"] = sha256(args.rows / "row_summary.json")
    provenance["row_status_sha256"] = sha256(args.rows / "row_status.json")
    provenance["analysis_code_sha256"] = {
        Path(__file__).name: sha256(__file__),
        "run_wl_shape_stack.py": sha256(HERE / "run_wl_shape_stack.py"),
        "s19a_shape_stack.py": sha256(HERE / "s19a_shape_stack.py"),
    }
    result["input"] = provenance
    result["software"] = {
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    out.mkdir(parents=True, exist_ok=False)
    with (out / "profiles.tsv").open("x", encoding="utf-8") as stream:
        fields = list(result["profiles"][0]) if result["profiles"] else []
        stream.write("\t".join(fields) + "\n")
        for row in result["profiles"]:
            stream.write("\t".join("" if row[k] is None else str(row[k]) for k in fields) + "\n")
    result["profiles_tsv_sha256"] = sha256(out / "profiles.tsv")
    with (out / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "input_pairs": result["input_pair_count"],
        "selected_pairs_before_radial_b_mode": result["published_source_cuts_and_zbins_count"],
        "profiles": len(result["profiles"]),
        "b_mode_polarity_decoded": result["b_mode_polarity_decoded"],
        "output": str(out),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", default=RADIAL_EDGES_MPC.tolist())
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    if len(args.radial_edges_mpc) < 2:
        parser.error("provide at least two radial edges")
    run(args)


if __name__ == "__main__":
    main()
