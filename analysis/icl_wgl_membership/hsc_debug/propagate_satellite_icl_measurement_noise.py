#!/usr/bin/env python3
"""Propagate measured HSC aperture-noise covariance into ICL--satellite r.

This conditional Monte Carlo perturbs the 57 association rows with verified
target/control pixel-noise covariance. It reports diagonal and finite-lag
noise models separately from spatial-block uncertainty; it is not a science
measurement, a total confidence interval, or a correction to the full 94-row
association sample.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_COVARIANCE = ROOT / "satellite_icl_target_control_covariance_2026-10-19_03"
DEFAULT_OUTPUT = ROOT / "satellite_icl_measurement_noise_mc_2026-10-19_01"
FILTERS = ("i2", "i")
VARIANTS = ("diagonal", "finite_lag")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_pearson(x, y, redshift, richness) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    redshift = np.asarray(redshift, dtype=float)
    richness = np.asarray(richness, dtype=float)
    if (x.ndim != 1 or y.shape != x.shape or redshift.shape != x.shape
            or richness.shape != x.shape or len(x) < 5
            or not all(np.isfinite(a).all() for a in (x, y, redshift, richness))):
        raise ValueError("partial correlation requires >=5 aligned finite rows")
    design = np.column_stack((np.ones(len(x)), redshift, richness))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("redshift/richness nuisance design is rank deficient")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not math.isfinite(denom) or denom <= 0:
        raise ValueError("partial correlation is undefined")
    return float(rx @ ry / denom)


def partial_pearson_batch(x, y, redshift, richness) -> np.ndarray:
    """Evaluate one partial correlation per row of aligned x/y matrices."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    z, richness = np.asarray(redshift, float), np.asarray(richness, float)
    if x.ndim != 2 or y.shape != x.shape or x.shape[1] != len(z) or len(z) != len(richness):
        raise ValueError("batch arrays have inconsistent dimensions")
    design = np.column_stack((np.ones(len(z)), z, richness))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("redshift/richness nuisance design is rank deficient")
    pinv = np.linalg.pinv(design)
    xr = x - (x @ pinv.T) @ design.T
    yr = y - (y @ pinv.T) @ design.T
    denom = np.linalg.norm(xr, axis=1) * np.linalg.norm(yr, axis=1)
    if np.any(~np.isfinite(denom)) or np.any(denom <= 0):
        raise ValueError("one or more partial correlations are undefined")
    return np.einsum("ij,ij->i", xr, yr) / denom


def covariance_cholesky(rows: list[dict], variant: str) -> np.ndarray:
    if variant not in VARIANTS:
        raise ValueError(f"unknown covariance variant: {variant}")
    suffix = "diagonal" if variant == "diagonal" else "lag"
    matrices = []
    for row in rows:
        vx = float(row[f"icl_delta_{suffix}_variance"])
        vy = float(row[f"target_satellite_{suffix}_variance_njy2_arcsec4"])
        cov = float(row[f"target_icl_satellite_{suffix}_covariance_njy2_arcsec4"])
        matrix = np.array([[vx, cov], [cov, vy]], dtype=float)
        if not np.isfinite(matrix).all() or vx <= 0 or vy <= 0:
            raise ValueError("noise covariance requires finite, positive variances")
        eig = np.linalg.eigvalsh(matrix)
        tolerance = 1e-10 * max(vx, vy, 1.0)
        if eig[0] < -tolerance:
            raise ValueError(f"noise covariance is not positive semidefinite: {eig[0]}")
        # Clamp only tiny roundoff-level negative eigenvalues, preserving the
        # supplied covariance entries to numerical precision.
        if eig[0] <= 0:
            values, vectors = np.linalg.eigh(matrix)
            matrix = (vectors * np.maximum(values, 0.0)) @ vectors.T
        matrices.append(np.linalg.cholesky(matrix + np.eye(2) * np.finfo(float).eps))
    return np.asarray(matrices)


def simulate_noise_r(partial_rows: list[dict], covariance_rows: list[dict],
                     variant: str, replicas: int, seed: int) -> np.ndarray:
    if len(partial_rows) != len(covariance_rows) or replicas < 100:
        raise ValueError("need aligned rows and at least 100 Monte Carlo draws")
    x = np.asarray([float(row["icl_delta_200_300"]) for row in partial_rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in partial_rows])
    z = np.asarray([float(row["redshift"]) for row in partial_rows])
    richness = np.asarray([float(row["richness"]) for row in partial_rows])
    chol = covariance_cholesky(covariance_rows, variant)
    rng = np.random.default_rng(seed)
    standard = rng.standard_normal((replicas, len(x), 2))
    error = np.einsum("nij,rnj->rni", chol, standard, optimize=True)
    return partial_pearson_batch(x[None, :] + error[:, :, 0],
                                 y[None, :] + error[:, :, 1], z, richness)


def load_inputs(association_dir: Path, covariance_dir: Path):
    association_dir, covariance_dir = association_dir.resolve(), covariance_dir.resolve()
    association_summary_path = association_dir / "summary.json"
    association_path = association_dir / "matched_cluster_values.tsv"
    covariance_summary_path = covariance_dir / "summary.json"
    covariance_path = covariance_dir / "target_control_pixel_covariance.tsv"
    association_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    covariance_summary = json.loads(covariance_summary_path.read_text(encoding="utf-8"))
    if (association_summary.get("schema") != "hsc-satellite-outer-icl-correlation-v1"
            or association_summary.get("status") != "complete"
            or association_summary.get("science_measurement") is not False
            or covariance_summary.get("schema") != "hsc-icl-satellite-target-control-pixel-covariance-v1"
            or covariance_summary.get("status") != "complete_exploratory"
            or covariance_summary.get("science_measurement") is not False):
        raise ValueError("input summaries are missing, changed, or not exploratory")
    assoc_sha = sha256(association_path)
    if association_summary.get("outputs_sha256", {}).get(association_path.name) != assoc_sha:
        raise ValueError("association input checksum mismatch")
    cov_sha = sha256(covariance_path)
    if covariance_summary.get("output_sha256", {}).get(str(covariance_path)) != cov_sha:
        raise ValueError("pixel-covariance input checksum mismatch")
    assoc, cov = read_tsv(association_path), read_tsv(covariance_path)
    if len(assoc) != 94 or len(cov) != 57:
        raise ValueError("expected 94 association rows and 57 covariance rows")
    assoc_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in assoc}
    cov_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in cov}
    if len(assoc_by_key) != len(assoc) or len(cov_by_key) != len(cov):
        raise ValueError("duplicate cluster/filter key in input tables")
    if not set(cov_by_key).issubset(assoc_by_key):
        raise ValueError("covariance input includes rows outside the association sample")
    if len(set(cov_by_key) & set(assoc_by_key)) != 57:
        raise ValueError("expected 57 exactly matched noise-covariance pairs")
    joined = {}
    for band in FILTERS:
        keys = sorted(key for key in cov_by_key if key[1] == band)
        expected = 44 if band == "i2" else 13
        if len(keys) != expected:
            raise ValueError(f"unexpected covariance support for {band}: {len(keys)}")
        if any(cov_by_key[key].get("assumed_target_control_noise_independent") != "True"
               for key in keys):
            raise ValueError("target/control noise-independence assumption is not recorded")
        joined[band] = ([assoc_by_key[key] for key in keys],
                        [cov_by_key[key] for key in keys])
    provenance = {
        "association_summary": str(association_summary_path),
        "association_summary_sha256": sha256(association_summary_path),
        "association_table": str(association_path),
        "association_table_sha256": assoc_sha,
        "covariance_summary": str(covariance_summary_path),
        "covariance_summary_sha256": sha256(covariance_summary_path),
        "covariance_table": str(covariance_path),
        "covariance_table_sha256": cov_sha,
        "n_association_rows_full_sample": len(assoc),
        "n_covariance_rows_joined": len(cov),
        "n_covariance_rows_by_filter": {band: len(joined[band][0]) for band in FILTERS},
    }
    return joined, provenance


def analyze(association_dir: Path, covariance_dir: Path, output: Path,
            replicas: int = 10000, seed: int = 20261019) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if replicas < 1000:
        raise ValueError("replicas must be >=1000")
    joined, provenance = load_inputs(association_dir, covariance_dir)
    by_filter = {}
    output_rows = []
    for filter_index, band in enumerate(FILTERS):
        assoc_rows, cov_rows = joined[band]
        point = partial_pearson(
            [float(row["icl_delta_200_300"]) for row in assoc_rows],
            [float(row["satellite_delta_50_75"]) for row in assoc_rows],
            [float(row["redshift"]) for row in assoc_rows],
            [float(row["richness"]) for row in assoc_rows],
        )
        results = {"n_pairs": len(assoc_rows), "partial_r_measured": point}
        for variant_index, variant in enumerate(VARIANTS):
            draw_seed = seed + 100 * filter_index + variant_index
            draws = simulate_noise_r(assoc_rows, cov_rows, variant, replicas, draw_seed)
            result = {
                "model": "fixed-pair Gaussian pixel-noise perturbation; no spatial resampling",
                "replicas": replicas,
                "seed": draw_seed,
                "partial_r_quantiles_2p5_16_50_84_97p5": np.percentile(
                    draws, [2.5, 16, 50, 84, 97.5]).tolist(),
                "fraction_partial_r_positive": float(np.mean(draws > 0)),
                "fraction_partial_r_abs_ge_measured": float(np.mean(np.abs(draws) >= abs(point))),
            }
            results[variant] = result
            output_rows.append({"filter": band, "covariance_model": variant,
                                "n_pairs": len(assoc_rows), "partial_r_measured": point,
                                **result})
        by_filter[band] = results
    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "measurement_noise_partial_r.tsv"
    fields = list(output_rows[0])
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        for row in output_rows:
            row = {key: json.dumps(value, separators=(",", ":")) if isinstance(value, (list, dict)) else value
                   for key, value in row.items()}
            writer.writerow(row)
    summary = {
        "schema": "astrafek-icl-satellite-measurement-noise-mc-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "population_scope": "57 exact association rows with audited target/control lag covariance; full association has 94 rows",
            "covariance_models": {
                "diagonal": "target plus matched-control ICL variance and target-image ICL/satellite covariance from delivered variance planes",
                "finite_lag": "same terms contracted through the measured symmetric 9x9 finite-lag residual covariance kernels",
            },
            "simulation": "independent per-pair bivariate Gaussian perturbations of measured ICL and satellite values; redshift and richness held fixed",
            "spatial_resampling": False,
            "replicas": replicas, "seed": seed,
            "reported_intervals": "conditional Monte Carlo quantiles under the stated noise model; not confidence intervals or a total uncertainty",
        },
        "provenance": provenance,
        "by_filter": by_filter,
        "results_table": str(table_path),
        "results_table_sha256": sha256(table_path),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "only the 57 covariance-available pairs are used; i has 13 and i2 has 44, not the full 94-row relation",
            "measurement perturbations assume Gaussian independent pair noise; spatial correlations between fields are not modeled",
            "the empirical lag kernel is stationary and conditional on fixed masks and detections",
            "source Poisson noise, upstream sky transfer, PSF-wing, extinction, selection, and astrophysical covariance are excluded",
            "this noise-only sensitivity is separate from spatial-block intervals to avoid double-counting field scatter",
            "the test does not select an estimator, mask, or scientific conclusion from observed association strength",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       allow_nan=False) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_pairs": provenance["n_covariance_rows_joined"],
                      "by_filter": by_filter,
                      "results_table_sha256": summary["results_table_sha256"],
                      "summary_sha256": summary["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-dir", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--covariance-dir", type=Path, default=DEFAULT_COVARIANCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    analyze(args.association_dir, args.covariance_dir, args.output,
            args.replicas, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
