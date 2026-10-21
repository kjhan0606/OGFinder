#!/usr/bin/env python3
"""Compare author-published DESI/HSC n(z) curves to the HSC Y3 SACC baseline.

The audit reruns the same cluster-only tomographic spatial jackknife for each
published scale cut and each published mean/magnification-corrected curve. The
pointwise lower/upper curves are treated only as marginal envelope sensitivities;
the NPZ does not contain joint posterior samples. No catalog rows are fetched.
"""
from __future__ import annotations

import argparse
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

from audit_wl_tomographic_consistency import tomographic_spatial_jackknife
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


SCALE_CUTS = {
    "0.3_3_v_1p1": (0.3, 3.0),
    "1_5_v_1p1": (1.0, 5.0),
}
CURVE_KINDS = ("mean", "mean_mu", "lower", "upper", "lower_mu", "upper_mu")
PINNED_COMMIT = "abc70af29910e49fce17a63ef3083160d515fb6c"


def load_desi_curves(directory: Path) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray]], dict]:
    """Validate the four approved files and load finite, normalized four-bin curves."""
    directory = Path(directory).resolve()
    curves: dict[str, dict[str, np.ndarray]] = {}
    hashes = {}
    z_common = None
    expected_variants = set(SCALE_CUTS)
    if not directory.is_dir():
        raise FileNotFoundError(directory)
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("commit") != PINNED_COMMIT:
        raise ValueError("DESI n(z) manifest is not pinned to the approved commit")
    manifest_rows = manifest.get("files", [])
    expected_paths = {
        f"results/distributions/{variant}/{filename}"
        for variant in expected_variants
        for filename in (f"distributions_{variant}.npz", f"nz_res_metadata_{variant}.json")
    }
    if {row.get("path") for row in manifest_rows} != expected_paths or len(manifest_rows) != 4:
        raise ValueError("DESI manifest does not contain exactly the four approved files")
    manifest_by_path = {row["path"]: row for row in manifest_rows}
    for variant in sorted(expected_variants):
        folder = directory / "results" / "distributions" / variant
        npz_path = folder / f"distributions_{variant}.npz"
        metadata_path = folder / f"nz_res_metadata_{variant}.json"
        if not npz_path.is_file() or not metadata_path.is_file():
            raise FileNotFoundError(f"missing approved DESI n(z) artifact for {variant}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if tuple(float(v) for v in metadata.get("scale_cuts", ())) != SCALE_CUTS[variant]:
            raise ValueError(f"metadata scale cut does not match {variant}")
        with np.load(npz_path, allow_pickle=False) as archive:
            z = np.asarray(archive["z"], dtype=float)
            if (z.ndim != 1 or len(z) < 2 or not np.all(np.isfinite(z))
                    or np.any(z < 0) or np.any(np.diff(z) <= 0)):
                raise ValueError(f"invalid redshift grid in {npz_path}")
            if z_common is None:
                z_common = z.copy()
            elif not np.array_equal(z_common, z):
                raise ValueError("DESI scale-cut files do not share the same redshift grid")
            by_kind = {}
            for kind in CURVE_KINDS:
                rows = []
                for source_bin in range(1, 5):
                    key = f"{source_bin}/nz_{kind}"
                    if key not in archive.files:
                        raise ValueError(f"missing required DESI n(z) array {key}")
                    curve = np.asarray(archive[key], dtype=float)
                    if curve.shape != z.shape or not np.all(np.isfinite(curve)) or np.any(curve < 0):
                        raise ValueError(f"invalid DESI n(z) array {key}")
                    integral = float(np.trapezoid(curve, z))
                    if not np.isfinite(integral) or integral <= 0:
                        raise ValueError(f"non-positive integral for DESI n(z) array {key}")
                    rows.append(curve / integral)
                by_kind[kind] = np.vstack(rows)
        curves[variant] = by_kind
        for file_path, rel in (
            (npz_path, f"results/distributions/{variant}/{npz_path.name}"),
            (metadata_path, f"results/distributions/{variant}/{metadata_path.name}"),
        ):
            row = manifest_by_path[rel]
            expected_url = f"https://raw.githubusercontent.com/JeanCHDJdev/desi-y3-hsc/{PINNED_COMMIT}/{rel}"
            observed_hash = sha(file_path)
            if (row.get("url") != expected_url or row.get("sha256") != observed_hash
                    or row.get("bytes") != file_path.stat().st_size):
                raise ValueError(f"DESI artifact does not match its pinned manifest row: {rel}")
            hashes[str(file_path)] = observed_hash
    hashes[str(manifest_path)] = sha(manifest_path)
    if z_common is None:
        raise ValueError("no DESI n(z) curves loaded")
    return z_common, curves, hashes


def profile_comparison(baseline: dict, candidate: dict) -> list[dict]:
    def keyed(summary):
        return {
            (int(row["source_zbin"]), float(row["r_min_mpc"]), float(row["r_max_mpc"])): row
            for row in summary["profile"]
        }
    base, current = keyed(baseline), keyed(candidate)
    if base.keys() != current.keys():
        raise ValueError("DESI and SACC profiles have different source-bin/radius keys")
    rows = []
    for key in sorted(base):
        a, b = base[key], current[key]
        av = float(a["delta_sigma_weak_msun_mpc2"])
        bv = float(b["delta_sigma_weak_msun_mpc2"])
        rows.append({
            "source_zbin": key[0], "r_min_mpc": key[1], "r_max_mpc": key[2],
            "sacc_delta_sigma_msun_mpc2": av,
            "desi_delta_sigma_msun_mpc2": bv,
            "change_msun_mpc2": bv - av,
            "fractional_change": (bv / av - 1.0) if av != 0 else None,
            "sacc_spatial_jackknife_error_msun_mpc2": float(a["delta_sigma_jackknife_error_msun_mpc2"]),
            "desi_spatial_jackknife_error_msun_mpc2": (
                float(b["delta_sigma_jackknife_error_msun_mpc2"])
                if b.get("delta_sigma_jackknife_error_msun_mpc2") is not None else None
            ),
        })
    return rows


def load_verified_baseline(path: Path, provenance: dict, nz_sacc: Path,
                           radial_edges: np.ndarray, cell_width_deg: float,
                           h0: float, omega_m: float) -> tuple[dict, str]:
    """Reuse a persisted SACC jackknife only when every analysis input matches."""
    path = Path(path).resolve()
    baseline = json.loads(path.read_text(encoding="utf-8"))
    if (baseline.get("schema") != "astrafek-hsc-y3-tomographic-consistency-jackknife-v1"
            or baseline.get("science_measurement") is not False
            or baseline.get("random_lens_subtraction_applied") is not False):
        raise ValueError("baseline summary has an unexpected schema or science state")
    record = baseline.get("input", {})
    comparisons = {
        "row_plan_sha256": provenance["row_plan_sha256"],
        "row_batch_sha256": provenance["row_batch_sha256"],
        "n_pair_rows": provenance["n_pair_rows"],
        "n_lenses": provenance["n_lenses"],
        "masked_value": provenance["masked_value"],
        "n_z_sacc_sha256": sha(nz_sacc),
        "radial_edges_mpc": [float(value) for value in radial_edges],
        "cell_width_equatorial_deg": float(cell_width_deg),
        "cosmology": {"H0_km_s_mpc": float(h0), "Omega_m": float(omega_m)},
    }
    for key, expected in comparisons.items():
        observed = baseline.get("cell_width_equatorial_deg") if key == "cell_width_equatorial_deg" else record.get(key)
        if observed != expected:
            raise ValueError(f"persisted SACC baseline input mismatch at {key}")
    if not baseline.get("profile") or baseline.get("n_occupied_cells") != 52:
        raise ValueError("persisted SACC baseline is missing the expected 52-cell profiles")
    return baseline, sha(path)


def run(args) -> dict:
    if args.masked_value not in (0, 1):
        raise ValueError("masked_value must be the previously verified catalog value 0 or 1")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    load_hsc_y3_nz_sacc(args.nz_sacc)  # validate that the approved SACC file remains readable
    desi_z, desi_curves, desi_hashes = load_desi_curves(args.desi_directory)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    edges = np.asarray(args.radial_edges_mpc, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges))
            or np.any(edges <= 0) or np.any(np.diff(edges) <= 0)):
        raise ValueError("radial edges must be positive, finite, and strictly increasing")

    baseline, baseline_summary_sha256 = load_verified_baseline(
        args.baseline_summary, provenance, args.nz_sacc, edges,
        args.cell_width_deg, args.h0, args.omega_m,
    )
    scenarios = {}
    for variant, scale_cut in SCALE_CUTS.items():
        for kind in CURVE_KINDS:
            label = f"{variant}__{kind}"
            if kind in ("mean", "mean_mu"):
                result = tomographic_spatial_jackknife(
                    pairs, desi_z, desi_curves[variant][kind], edges, cosmology,
                    cell_width_deg=args.cell_width_deg,
                )
                uncertainty_method = "same spatial jackknife as SACC baseline"
            else:
                result = {
                    "schema": "hsc-s19a-tomographic-nz-envelope-stack-v1",
                    "status": "pointwise_marginal_envelope_sensitivity",
                    "science_measurement": False,
                    "profile": [],
                }
                assigned = np.asarray(pairs["source_zbin"], dtype=int)
                for source_bin in range(1, 5):
                    mask = assigned == source_bin
                    if not np.any(mask):
                        raise ValueError(f"no verified shape rows in source bin {source_bin}")
                    subset = {key: np.asarray(value)[mask] for key, value in pairs.items()}
                    stack = stack_s19a_tomographic_pairs(
                        subset, desi_z, desi_curves[variant][kind], edges, cosmology,
                    )
                    result["profile"].extend(
                        {"source_zbin": source_bin, **row} for row in stack["profile"]
                    )
                uncertainty_method = "none; pointwise published lower/upper curves only, no paired covariance"
            scenarios[label] = {
                "scale_cut_mpc_over_h": list(scale_cut),
                "curve_kind": kind,
                "is_central_estimate": kind in ("mean", "mean_mu"),
                "is_magnification_corrected": kind.endswith("_mu"),
                "is_pointwise_marginal_envelope": kind.startswith(("lower", "upper")),
                "uncertainty_method": uncertainty_method,
                "spatial_jackknife": result,
                "vs_sacc_baseline": profile_comparison(baseline, result),
            }

    output.mkdir(parents=True, exist_ok=False)
    payload = {
        "schema": "astrafek-hsc-pdr3-desi-nz-spatial-jackknife-sensitivity-v1",
        "status": "complete_cluster_only_nz_sensitivity",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "photoz_posterior_samples_available": False,
        "n_z_uncertainty_treatment": (
            "Recomputes each stack with the authors' pointwise lower and upper n(z) curves; "
            "these are marginal-envelope sensitivities, not joint posterior draws or a covariance."
        ),
        "method": "same verified 800-lens shape-row export; reused its exactly provenance-matched SACC spatial-jackknife baseline; recomputed DESI mean-curve spatial jackknives, and computed pointwise envelope profiles without assigning them covariance",
        "input": {
            **provenance,
            "shape_rows_path": str(args.rows.resolve()),
            "preflight_path": str(args.preflight.resolve()),
            "plan_path": str(args.plan_dir.resolve()),
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "baseline_summary_path": str(args.baseline_summary.resolve()),
            "baseline_summary_sha256": baseline_summary_sha256,
            "desi_calibration_directory": str(args.desi_directory.resolve()),
            "desi_artifact_sha256": desi_hashes,
            "radial_edges_mpc": [float(v) for v in edges],
            "cell_width_deg": float(args.cell_width_deg),
            "cosmology": {"H0_km_s_Mpc": args.h0, "Omega_m": args.omega_m},
        },
        "baseline_sacc": baseline,
        "scenarios": scenarios,
        "limitations": [
            "No matched random-lens subtraction is available; all profiles are cluster-only diagnostics.",
            "The published lower/upper curves do not encode correlated posterior draws across redshift or tomographic bins.",
            "The authors' calibrated n(z) curves are sample-level distributions, not individual source P(z) PDFs.",
            "The equal-area cells are an exploratory spatial partition and can join disconnected HSC fields.",
            "The public repository reports no SPDX license; the associated paper is CC BY 4.0. This data import retains source citation and is not a license determination.",
        ],
        "sources": [
            "Choppin de Janvry et al. 2026, JCAP 05 073, arXiv:2511.18133v3",
            "https://github.com/JeanCHDJdev/desi-y3-hsc",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"], "science_measurement": False,
        "n_scenarios": len(scenarios),
        "n_pair_rows": provenance["n_pair_rows"],
        "n_occupied_cells": baseline["n_occupied_cells"],
        "output": str(output), "summary_sha256": sha(summary_path),
    }, ensure_ascii=False))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--desi-directory", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", default=[0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
