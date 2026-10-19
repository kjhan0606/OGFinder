#!/usr/bin/env python3
"""Independently validate the saved five-scenario HSC n(z) shift grid."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


EXPECTED_SHIFTS = {
    "central": (0.043, 0.050),
    "low3_low4": (0.021, 0.038),
    "low3_high4": (0.021, 0.062),
    "high3_low4": (0.063, 0.038),
    "high3_high4": (0.063, 0.062),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def close(left, right, *, rel=1e-10, abs_tol=1e-12):
    return math.isclose(float(left), float(right), rel_tol=rel, abs_tol=abs_tol)


def validate_scenario(name, scenario, top_input):
    if (scenario.get("schema") != "astrafek-hsc-y3-published-nz-shift-sensitivity-v1"
            or scenario.get("status") != "approximate_sensitivity"
            or scenario.get("science_measurement") is not False):
        raise ValueError(f"scenario has invalid schema/status: {name}")
    shifts = scenario.get("applied_mean_redshift_shifts", {})
    expected3, expected4 = EXPECTED_SHIFTS[name]
    for source_bin, expected in (("1", 0.0), ("2", 0.0),
                                 ("3", expected3), ("4", expected4)):
        detail = shifts.get(source_bin)
        if not isinstance(detail, dict) or not close(detail.get("requested_mean_shift"), expected):
            raise ValueError(f"scenario requested shift mismatch: {name}, bin {source_bin}")
        realized = detail.get("realized_mean_shift_after_grid_truncation")
        if not close(realized, expected, rel=0.0, abs_tol=1e-4):
            raise ValueError(f"scenario realized shift outside tolerance: {name}, bin {source_bin}")

    tomography = scenario.get("shifted_tomography", {})
    if (tomography.get("science_measurement") is not False
            or tomography.get("random_lens_subtraction_applied") is not False
            or tomography.get("n_z_uncertainty_propagated") is not False
            or tomography.get("n_occupied_cells") != 52
            or tomography.get("n_input_pairs_in_tomographic_bins") != 2254646):
        raise ValueError(f"tomographic gate/count mismatch: {name}")
    covariance = np.asarray(tomography.get("covariance"), dtype=float)
    if covariance.shape != (56, 56) or not np.isfinite(covariance).all():
        raise ValueError(f"invalid covariance dimensions or values: {name}")
    scale = max(float(np.max(np.abs(covariance))), 1e-300)
    symmetry_error = float(np.max(np.abs(covariance - covariance.T)) / scale)
    if symmetry_error > 1e-12:
        raise ValueError(f"covariance is not symmetric: {name}")
    rank = int(np.linalg.matrix_rank(covariance))
    if rank != 51:
        raise ValueError(f"unexpected full covariance rank: {name} ({rank})")
    if len(tomography.get("profile", [])) != 28:
        raise ValueError(f"profile does not cover 4 bins x 7 radii: {name}")
    tests = scenario.get("baseline_vs_shifted_consistency", [])
    if len(tests) != 7:
        raise ValueError(f"consistency table does not cover seven radial bins: {name}")
    for row in tests:
        for key in ("shifted_chi2", "shifted_p_value", "shifted_holm_p_value_7_radial_tests"):
            value = float(row[key])
            if not math.isfinite(value) or (key.endswith("p_value_7_radial_tests") and not 0 <= value <= 1):
                raise ValueError(f"invalid consistency statistic {key}: {name}")
    input_record = scenario.get("input", {})
    for key in ("row_plan_sha256", "row_batch_sha256", "n_pair_rows", "n_lenses",
                "masked_value", "n_z_sacc_sha256"):
        if input_record.get(key) != top_input.get(key):
            raise ValueError(f"scenario input provenance differs from grid: {name}/{key}")
    return {
        "requested_shift_bin3": expected3,
        "requested_shift_bin4": expected4,
        "covariance_shape": list(covariance.shape),
        "covariance_rank": rank,
        "relative_covariance_symmetry_error": symmetry_error,
        "spatial_cells": 52,
        "tomographic_pair_rows": 2254646,
        "radial_tests": 7,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite validation output: {args.output}")
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-hsc-y3-published-nz-shift-uncertainty-grid-v1"
            or summary.get("status") != "approximate_sensitivity_grid"
            or summary.get("science_measurement") is not False):
        raise ValueError("unexpected five-scenario grid schema/status")
    runner = Path(__file__).with_name("audit_wl_published_nz_shift_sensitivity.py")
    runner_hash = sha256(runner)
    top_input = summary.get("input", {})
    if top_input.get("runner_sha256") != runner_hash:
        raise ValueError("grid runner hash mismatch")
    if sha256(Path(top_input["n_z_sacc_path"])) != top_input.get("n_z_sacc_sha256"):
        raise ValueError("SACC checksum mismatch")
    if sha256(Path(top_input["baseline_summary_path"])) != top_input.get("baseline_summary_sha256"):
        raise ValueError("baseline summary checksum mismatch")
    scenarios = summary.get("scenarios", {})
    if set(scenarios) != set(EXPECTED_SHIFTS):
        raise ValueError("scenario labels differ from the planned central-plus-corners grid")
    scenario_audits = {
        name: validate_scenario(name, scenarios[name], top_input)
        for name in EXPECTED_SHIFTS
    }

    radii = [(float(row["r_min_mpc"]), float(row["r_max_mpc"]))
             for row in scenarios["central"]["baseline_vs_shifted_consistency"]]
    expected_envelope = []
    for lo, hi in radii:
        rows = []
        for scenario in scenarios.values():
            rows.append(next(row for row in scenario["baseline_vs_shifted_consistency"]
                             if (float(row["r_min_mpc"]), float(row["r_max_mpc"])) == (lo, hi)))
        expected_envelope.append({
            "r_min_mpc": lo,
            "r_max_mpc": hi,
            "min_shifted_chi2": min(float(row["shifted_chi2"]) for row in rows),
            "max_shifted_chi2": max(float(row["shifted_chi2"]) for row in rows),
            "min_shifted_holm_p_value": min(float(row["shifted_holm_p_value_7_radial_tests"]) for row in rows),
            "max_shifted_holm_p_value": max(float(row["shifted_holm_p_value_7_radial_tests"]) for row in rows),
        })
    envelope = summary.get("tomographic_consistency_envelope", [])
    if len(envelope) != len(expected_envelope):
        raise ValueError("stored tomography envelope has the wrong number of radial bins")
    for expected, observed in zip(expected_envelope, envelope):
        for key, value in expected.items():
            if key not in observed or not close(observed[key], value):
                raise ValueError(f"stored envelope disagrees with scenarios at {key}")

    audit = {
        "schema": "astrafek-hsc-y3-published-nz-shift-uncertainty-grid-validation-v1",
        "status": "validated_approximate_sensitivity_grid",
        "science_measurement": False,
        "source_summary": str(args.summary.resolve()),
        "source_summary_sha256": sha256(args.summary),
        "runner": str(runner.resolve()),
        "runner_sha256": runner_hash,
        "scenario_audits": scenario_audits,
        "tomographic_consistency_envelope": expected_envelope,
        "n_scenarios": len(scenarios),
        "no_new_hsc_data_requested": True,
        "limitations": [
            "Validation confirms internal schema, hashes, covariance structure, and aggregation only.",
            "The calculation remains a constant-translation approximation, not calibrated DESI spline n(z) or posterior propagation.",
            "It remains cluster-only without random-lens subtraction or connected-footprint covariance.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": audit["status"],
        "science_measurement": audit["science_measurement"],
        "n_scenarios": audit["n_scenarios"],
        "source_summary_sha256": audit["source_summary_sha256"],
        "output": str(args.output),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
