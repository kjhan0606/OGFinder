#!/usr/bin/env python3
"""Independently validate an Astrafek source-reuse covariance summary.

This checks the saved covariance matrices, profile errors, tomographic GLS
consistency tests, Holm adjustment, deletion bookkeeping, and source-reuse
arithmetic without loading or recalculating from the source catalogue. It does
not certify the estimator or turn the cluster-only result into a science
measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2


MODES = ("lens_only", "source_only", "joint_endpoints")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def close(actual: float, expected: float, label: str, *, atol: float = 1e-10,
          rtol: float = 1e-9) -> None:
    if not np.isfinite(actual) or not np.isclose(actual, expected, atol=atol, rtol=rtol):
        raise ValueError(f"{label} differs: {actual!r} != {expected!r}")


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * float(p_values[index]))
        adjusted[index] = min(1.0, running)
    return adjusted


def fit_common(values: np.ndarray, covariance: np.ndarray) -> dict:
    covariance = 0.5 * (covariance + covariance.T)
    inverse = np.linalg.pinv(covariance, hermitian=True)
    ones = np.ones(len(values))
    denominator = float(ones @ inverse @ ones)
    if denominator <= 0:
        raise ValueError("tomographic covariance has no positive common-mode information")
    common = float(ones @ inverse @ values / denominator)
    residual = values - common
    statistic = float(residual @ inverse @ residual)
    rank = int(np.linalg.matrix_rank(covariance))
    dof = rank - 1
    p_value = float(chi2.sf(statistic, dof)) if dof > 0 else None
    return {"common_value": common, "chi2": statistic,
            "degrees_of_freedom": dof, "p_value": p_value,
            "covariance_rank": rank}


def validate_summary(payload: dict) -> dict:
    if payload.get("schema") != "astrafek-hsc-y3-source-reuse-covariance-audit-v1":
        raise ValueError("unexpected source-reuse covariance schema")
    if payload.get("status") != "verified_covariance_sensitivity":
        raise ValueError("producer did not mark covariance sensitivity verified")
    if payload.get("science_measurement") is not False:
        raise ValueError("cluster-only covariance audit must remain science_measurement=false")
    if payload.get("random_lens_subtraction_applied") is not False:
        raise ValueError("unexpected random-lens subtraction flag")

    input_record = payload.get("input", {})
    edges = np.asarray(input_record.get("radial_edges_mpc", []), dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
            or np.any(np.diff(edges) <= 0)):
        raise ValueError("invalid radial-edge provenance")
    n_radial = len(edges) - 1
    if n_radial != 7:
        raise ValueError("source-reuse summary must contain the seven approved radial tests")
    dimension = 4 * 2 * n_radial
    source_bin_indices = np.arange(4) * 2 * n_radial
    modes = payload.get("covariance_modes", {})
    if set(modes) != set(MODES):
        raise ValueError("covariance modes do not match lens/source/joint endpoint audit")

    validated_modes = {}
    for mode_name in MODES:
        mode = modes[mode_name]
        matrix = np.asarray(mode.get("covariance"), dtype=float)
        if matrix.shape != (dimension, dimension) or not np.isfinite(matrix).all():
            raise ValueError(f"{mode_name} covariance has an invalid shape or non-finite values")
        scale = max(float(np.max(np.abs(matrix))), 1e-300)
        symmetry_error = float(np.max(np.abs(matrix - matrix.T)) / scale)
        close(float(mode.get("relative_covariance_symmetry_error")), symmetry_error,
              f"{mode_name} reported symmetry error", atol=1e-12, rtol=1e-6)
        if symmetry_error > 1e-10:
            raise ValueError(f"{mode_name} covariance is not symmetric")
        symmetric = 0.5 * (matrix + matrix.T)
        eigen_min = float(np.linalg.eigvalsh(symmetric)[0])
        psd_tolerance_relative = 1e-12
        if eigen_min < -psd_tolerance_relative * scale:
            raise ValueError(f"{mode_name} covariance is not positive semidefinite")
        rank = int(np.linalg.matrix_rank(symmetric))
        if int(mode.get("covariance_rank", -1)) != rank:
            raise ValueError(f"{mode_name} covariance rank mismatch")
        if int(mode.get("covariance_dimension", -1)) != dimension:
            raise ValueError(f"{mode_name} covariance dimension mismatch")

        profile = mode.get("profile", [])
        if len(profile) != 4 * n_radial:
            raise ValueError(f"{mode_name} profile does not cover all source/radial bins")
        vector = np.empty(dimension, dtype=float)
        for source_index in range(4):
            for radial_index in range(n_radial):
                row_index = source_index * n_radial + radial_index
                row = profile[row_index]
                if int(row.get("source_zbin", -1)) != source_index + 1:
                    raise ValueError(f"{mode_name} profile source-bin ordering mismatch")
                close(float(row.get("r_min_mpc")), float(edges[radial_index]),
                      f"{mode_name} profile lower edge")
                close(float(row.get("r_max_mpc")), float(edges[radial_index + 1]),
                      f"{mode_name} profile upper edge")
                pos = source_index * 2 * n_radial + 2 * radial_index
                vector[pos] = float(row["delta_sigma_weak_msun_mpc2"])
                vector[pos + 1] = float(row["cross_delta_sigma_weak_msun_mpc2"])
                for key, diagonal in (
                    ("delta_sigma_jackknife_error_msun_mpc2", matrix[pos, pos]),
                    ("cross_jackknife_error_msun_mpc2", matrix[pos + 1, pos + 1]),
                ):
                    expected = float(np.sqrt(max(0.0, diagonal)))
                    close(float(row[key]), expected, f"{mode_name} profile {key}")
        if not np.isfinite(vector).all():
            raise ValueError(f"{mode_name} full profile contains non-finite values")

        consistency = mode.get("tomographic_consistency", [])
        if len(consistency) != n_radial:
            raise ValueError(f"{mode_name} has an incomplete tomographic consistency table")
        raw_p = []
        consistency_checks = []
        for radial_index, row in enumerate(consistency):
            if int(row.get("radial_bin", -1)) != radial_index:
                raise ValueError(f"{mode_name} tomographic radial ordering mismatch")
            close(float(row.get("r_min_mpc")), float(edges[radial_index]),
                  f"{mode_name} consistency lower edge")
            close(float(row.get("r_max_mpc")), float(edges[radial_index + 1]),
                  f"{mode_name} consistency upper edge")
            indices = source_bin_indices + 2 * radial_index
            values = vector[indices]
            reported = np.asarray(row.get("source_bin_delta_sigma_msun_mpc2"), dtype=float)
            if reported.shape != (4,):
                raise ValueError(f"{mode_name} consistency row lacks four source-bin values")
            if not np.allclose(reported, values, atol=1e-10, rtol=1e-9):
                raise ValueError(f"{mode_name} consistency values disagree with saved profiles")
            fit = fit_common(values, symmetric[np.ix_(indices, indices)])
            for key in ("common_value", "chi2"):
                close(float(row.get(key)), fit[key], f"{mode_name} consistency {key}")
            if int(row.get("degrees_of_freedom", -1)) != fit["degrees_of_freedom"]:
                raise ValueError(f"{mode_name} consistency degrees of freedom mismatch")
            if int(row.get("covariance_rank", -1)) != fit["covariance_rank"]:
                raise ValueError(f"{mode_name} consistency covariance rank mismatch")
            p_value = row.get("p_value")
            if fit["p_value"] is None:
                if p_value is not None:
                    raise ValueError(f"{mode_name} degenerate consistency p-value should be null")
                raw_p.append(1.0)
            else:
                close(float(p_value), fit["p_value"], f"{mode_name} consistency p-value")
                raw_p.append(float(p_value))
            consistency_checks.append(row)

        adjusted = holm_adjust(np.asarray(raw_p, dtype=float))
        for row, expected in zip(consistency_checks, adjusted):
            close(float(row.get("holm_p_value_7_radial_tests")), float(expected),
                  f"{mode_name} Holm-adjusted p-value")

        deletions = mode.get("block_deletions", [])
        if len(deletions) != int(mode.get("n_deletion_blocks", -1)) or len(deletions) < 3:
            raise ValueError(f"{mode_name} deletion count is inconsistent")
        block_ids = [int(row["block_id"]) for row in deletions]
        if len(set(block_ids)) != len(block_ids):
            raise ValueError(f"{mode_name} deletion block IDs are not unique")
        for row in deletions:
            if int(row.get("deleted_pair_rows", -1)) < 0:
                raise ValueError(f"{mode_name} has a negative deleted-row count")
            remaining = row.get("used_pairs_remaining_by_source_bin", [])
            if len(remaining) != 4 or any(int(value) <= 0 for value in remaining):
                raise ValueError(f"{mode_name} has invalid remaining source-bin counts")

        validated_modes[mode_name] = {
            "n_deletion_blocks": len(deletions),
            "covariance_rank": rank,
            "covariance_dimension": dimension,
            "maximum_absolute_covariance_entry": scale,
            "minimum_covariance_eigenvalue": eigen_min,
            "minimum_covariance_eigenvalue_relative_to_entry_scale": eigen_min / scale,
            "psd_tolerance_relative_to_entry_scale": psd_tolerance_relative,
            "maximum_holm_p_value": float(np.max(adjusted)),
            "profile_rows_checked": len(profile),
        }

    reuse = payload.get("source_reuse", {})
    selected = int(reuse.get("selected_pair_rows", -1))
    unique = int(reuse.get("unique_source_objects", -1))
    repeated = int(reuse.get("source_objects_reused_across_lens_pairs", -1))
    repeated_rows = int(reuse.get("pair_rows_from_reused_sources", -1))
    maximum = int(reuse.get("maximum_lens_pair_multiplicity_per_source", -1))
    if not (selected > 0 and 0 < unique <= selected and 0 <= repeated <= unique
            and repeated_rows <= selected and maximum >= 1):
        raise ValueError("source-reuse counts violate basic bounds")
    if (repeated == 0 and repeated_rows != 0) or (repeated > 0 and repeated_rows < 2 * repeated):
        raise ValueError("reused-source row count is inconsistent with reused object count")
    close(float(reuse.get("fraction_unique_sources_reused")), repeated / unique,
          "fraction of unique sources reused")
    close(float(reuse.get("fraction_pair_rows_from_reused_sources")), repeated_rows / selected,
          "fraction of pair rows from reused sources")
    for key in ("reused_sources_spanning_multiple_lens_endpoint_cells",
                "reused_sources_spanning_multiple_source_endpoint_cells"):
        if not (0 <= int(reuse.get(key, -1)) <= repeated):
            raise ValueError(f"source-reuse cell count out of bounds: {key}")
    if reuse.get("source_endpoint_cell_assignment_is_unique_per_object_id") is not True:
        raise ValueError("source IDs do not have unique source-endpoint cell assignments")

    return {
        "schema": "astrafek-hsc-y3-source-reuse-covariance-validation-v1",
        "status": "verified",
        "science_measurement": False,
        "source_summary_schema": payload["schema"],
        "selected_pair_rows": selected,
        "unique_source_objects": unique,
        "reused_source_objects": repeated,
        "covariance_modes": validated_modes,
        "limitations": [
            "Validates internal numerical consistency of the saved output; it does not rerun the estimator from source rows.",
            "Does not supply random-lens subtraction or certify connected HSC survey-region covariance.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary_path = args.summary.resolve()
    output_path = args.output.resolve()
    if not summary_path.is_file():
        parser.error(f"summary does not exist: {summary_path}")
    if output_path.exists():
        parser.error(f"refusing to overwrite validation output: {output_path}")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    result = validate_summary(payload)
    result.update({"source_summary": str(summary_path),
                   "source_summary_sha256": sha256(summary_path),
                   "validator_sha256": sha256(Path(__file__).resolve())})
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"],
                      "covariance_modes": result["covariance_modes"],
                      "science_measurement": result["science_measurement"],
                      "output": str(output_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
