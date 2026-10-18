#!/usr/bin/env python3
"""Run both stored-bit hypotheses for the HSC Y3 calibrated cluster stack.

This is a sensitivity calculation, not a decoder of the undocumented stored
``b_mode_mask`` polarity. It reads only the previously approved local shape-row
export and SACC n(z) file. The product always remains exploratory.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
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

from fetch_ancillary_catalogues import sha
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import (
    decode_hsc_b_mode_mask,
    load_hsc_y3_nz_sacc,
    stack_s19a_tomographic_pairs,
)
from submit_wl_shape_rows import DEFAULT_PLAN, DEFAULT_PREFLIGHT, DEFAULT_RESULTS
from ogfkit.wlcluster import Cosmology


DEFAULT_NZ_SACC = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
)


def calculate_sensitivity(pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology):
    raw = np.asarray(pairs["b_mode_mask_raw"], dtype=np.int8)
    scenarios = {}
    for masked_value in (0, 1):
        case_pairs = dict(pairs)
        case_pairs["b_mode_pass"] = decode_hsc_b_mode_mask(
            raw, masked_value=masked_value,
        )
        stack = stack_s19a_tomographic_pairs(
            case_pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
        )
        scenarios[f"exclude_raw_{masked_value}_hypothesis"] = {
            "hypothesis_masked_raw_value": masked_value,
            "retained_raw_value": 1 - masked_value,
            "stack": stack,
        }

    zero = scenarios["exclude_raw_0_hypothesis"]["stack"]
    one = scenarios["exclude_raw_1_hypothesis"]["stack"]
    differences = []
    for row0, row1 in zip(zero["profile"], one["profile"]):
        differences.append({
            "bin": row0["bin"],
            "r_min_mpc": row0["r_min_mpc"],
            "r_max_mpc": row0["r_max_mpc"],
            "delta_sigma_exclude_0_minus_exclude_1_msun_mpc2": (
                row0["delta_sigma_weak_msun_mpc2"] - row1["delta_sigma_weak_msun_mpc2"]
                if row0["delta_sigma_weak_msun_mpc2"] is not None
                and row1["delta_sigma_weak_msun_mpc2"] is not None else None
            ),
            "cross_delta_sigma_exclude_0_minus_exclude_1_msun_mpc2": (
                row0["cross_delta_sigma_weak_msun_mpc2"]
                - row1["cross_delta_sigma_weak_msun_mpc2"]
                if row0["cross_delta_sigma_weak_msun_mpc2"] is not None
                and row1["cross_delta_sigma_weak_msun_mpc2"] is not None else None
            ),
        })
    return scenarios, differences


def run(args):
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    for path in (args.nz_sacc, args.rows, args.preflight, args.plan_dir):
        if not path.exists():
            raise FileNotFoundError(path)

    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, masked_value=None,
    )
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    scenarios, differences = calculate_sensitivity(
        pairs, nz_redshift, nz_by_bin, args.radial_edges_mpc, cosmology,
    )
    payload = {
        "schema": "hsc-s19a-bmode-polarity-sensitivity-v1",
        "status": "complete_sensitivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "science_measurement": False,
        "b_mode_polarity_decoded": False,
        "random_lens_subtraction_applied": False,
        "spatial_block_covariance_applied": False,
        "photoz_uncertainty_propagated": False,
        "interpretation": (
            "Both possible stored-bit exclusions are reported. Neither scenario "
            "is selected as the HSC mask convention."
        ),
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(x) for x in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
        },
        "scenarios": scenarios,
        "profile_difference_exclude_0_minus_exclude_1": differences,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    print(json.dumps({
        "status": payload["status"],
        "science_measurement": False,
        "input_pairs": provenance["n_pair_rows"],
        "scenarios": {
            name: {
                "masked_hypothesis": result["hypothesis_masked_raw_value"],
                "used_pairs": result["stack"]["counts"]["used_pairs"],
            }
            for name, result in scenarios.items()
        },
        "output": str(output),
        "sha256": digest,
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ_SACC)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+",
                        default=[0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    if len(args.radial_edges_mpc) < 2:
        parser.error("provide at least two radial edges")
    run(args)


if __name__ == "__main__":
    main()
