#!/usr/bin/env python3
"""Explore the existing i2 ICL–satellite–weak-lensing triple overlap.

The input is restricted to already saved, audited HSC-derived tables. This
script performs no network access or HSC query and always labels the results
exploratory (science_measurement=false).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SATELLITE = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
WL = DATA / "wl_continuous_icl_relation_i2_2026-11-06_02/per_lens_profiles.tsv"
OVERLAP = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03/summary.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_wl_triple_overlap_2026-11-06_01"
EXPECTED_HASHES = {
    SATELLITE: "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da",
    WL: "d3e07ea875511027739a1a519a7beac765bd00a4ccd3c6b40e23dcc3beb789d8",
    OVERLAP: "98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12",
}
N_BOOT = 5000
SEED = 20261106


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_corr(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]) -> float:
    arrays = [np.asarray(x, float), np.asarray(y, float), *[np.asarray(c, float) for c in covariates]]
    if len(arrays[0]) < 6 or any(a.ndim != 1 or len(a) != len(arrays[0]) for a in arrays):
        return float("nan")
    if any(not np.isfinite(a).all() for a in arrays):
        return float("nan")
    design = np.column_stack((np.ones(len(arrays[0])), *arrays[2:]))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return float("nan")
    rx = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    ry = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    norm = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / norm) if norm > 0 else float("nan")


def block_bootstrap(rows: list[dict[str, str]], xkey: str, ykey: str,
                    covarkeys: list[str], replicas: int, seed: int) -> dict:
    coords_t = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])} for r in rows]
    coords_c = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])} for r in rows]
    _, labels, _, _, n_blocks, largest = pair_block_labels(
        coords_t, coords_c, np.ones(len(rows), dtype=bool)
    )
    blocks = np.unique(labels)
    if n_blocks < 6:
        raise ValueError(f"Only {n_blocks} spatial blocks for {xkey} vs {ykey}")
    by_block = {b: np.flatnonzero(labels == b) for b in blocks}

    def coefficient(indices: np.ndarray) -> float:
        return partial_corr(
            np.asarray([float(rows[i][xkey]) for i in indices]),
            np.asarray([float(rows[i][ykey]) for i in indices]),
            [np.asarray([float(rows[i][k]) for i in indices]) for k in covarkeys],
        )

    point = coefficient(np.arange(len(rows)))
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        sampled = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([by_block[b] for b in sampled])
        value = coefficient(indices)
        if math.isfinite(value):
            draws.append(value)
    if not draws:
        lo = hi = float("nan")
    else:
        lo, hi = np.quantile(draws, [0.025, 0.975]).tolist()
    return {
        "partial_r": point,
        "block_bootstrap_95_low": float(lo),
        "block_bootstrap_95_high": float(hi),
        "bootstrap_replicates": replicas,
        "bootstrap_valid": len(draws),
        "bootstrap_valid_fraction": len(draws) / replicas,
        "n_blocks": n_blocks,
        "largest_block_n": largest,
    }


def finite_row(row: dict[str, str], keys: list[str]) -> bool:
    try:
        return all(math.isfinite(float(row[k])) for k in keys)
    except (KeyError, TypeError, ValueError):
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite", type=Path, default=SATELLITE)
    parser.add_argument("--wl", type=Path, default=WL)
    parser.add_argument("--overlap", type=Path, default=OVERLAP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap", type=int, default=N_BOOT)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"Refusing to overwrite existing output: {args.output}")
    if args.bootstrap < 100:
        parser.error("Use at least 100 block-bootstrap replicates")
    paths = {"satellite": args.satellite, "wl": args.wl, "overlap": args.overlap}
    hashes = {name: sha256(path) for name, path in paths.items()}
    for name, path in paths.items():
        if path.resolve() == {"satellite": SATELLITE, "wl": WL, "overlap": OVERLAP}[name].resolve():
            if hashes[name] != EXPECTED_HASHES[path.resolve()]:
                raise ValueError(f"Pinned input hash mismatch: {path}")

    overlap = json.loads(args.overlap.read_text())
    if overlap.get("science_measurement") is not False:
        raise ValueError("Overlap audit no longer has science_measurement=false")
    expected_ids = set(overlap["joint_icl_satellite_support_by_filter"]["i2"]["cluster_ids_in_common_support"])
    sat_rows = [r for r in read_tsv(args.satellite) if r["filter"] == "i2"]
    if any(r.get("science_measurement", "False").lower() != "false" for r in sat_rows):
        raise ValueError("Satellite support table contains a science_measurement promotion")
    sat_by_id = {r["cluster_id"]: r for r in sat_rows}
    if len(sat_by_id) != len(sat_rows):
        raise ValueError("Duplicate i2 cluster IDs in satellite support table")
    wl_rows = read_tsv(args.wl)
    wl_by_id = {r["lens_id"]: r for r in wl_rows}
    if len(wl_by_id) != len(wl_rows):
        raise ValueError("Duplicate lens IDs in continuous WL table")
    if not expected_ids <= set(sat_by_id) or not expected_ids <= set(wl_by_id):
        raise ValueError("Input tables do not cover the audited 27-lens i2 intersection")
    if len(expected_ids) != 27:
        raise ValueError(f"Expected 27 audited i2 lenses, got {len(expected_ids)}")

    rows = []
    for cluster_id in sorted(expected_ids, key=int):
        sat, wl = sat_by_id[cluster_id], wl_by_id[cluster_id]
        for key, skey, wkey in (
            ("redshift", "redshift", "redshift"),
            ("richness", "richness", "richness"),
            ("icl", "icl_delta_200_300", "icl_delta_sb_200_300_njy_arcsec2"),
        ):
            a, b = float(sat[skey]), float(wl[wkey])
            if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-8):
                raise ValueError(f"{key} mismatch for lens {cluster_id}: {a} vs {b}")
        merged = dict(sat)
        merged.update({f"wl_{k}": v for k, v in wl.items() if k != "lens_id"})
        merged["log_richness"] = f"{math.log(float(sat['richness'])):.16g}"
        rows.append(merged)

    metrics = []
    idx = 0
    # Direct ICL–satellite association in the triple-overlap sample.
    pair = block_bootstrap(rows, "icl_delta_200_300", "satellite_delta_50_75",
                           ["redshift", "log_richness"], args.bootstrap, SEED)
    metrics.append({"relationship": "icl_vs_satellite", "component": "na", "radial_bin": "na",
                    "n_lenses": len(rows), "n_target_pairs_median": "na",
                    "n_random_pairs_median": "na", **pair})

    # Conditional WL relations; require usable target and random support per lens/bin.
    for radial_bin in range(7):
        for component, ykey in (
            ("tangential", f"wl_tangential_bin{radial_bin}_msun_mpc2"),
            ("cross", f"wl_cross_bin{radial_bin}_msun_mpc2"),
        ):
            target_key = f"wl_target_pairs_bin{radial_bin}"
            random_key = f"wl_random_pairs_bin{radial_bin}"
            usable = [r for r in rows if finite_row(r, [ykey, "icl_delta_200_300",
                        "satellite_delta_50_75", "redshift", "log_richness", target_key, random_key])
                      and int(float(r[target_key])) > 0 and int(float(r[random_key])) > 0]
            support = {
                "n_lenses": len(usable),
                "n_target_pairs_median": float(np.median([float(r[target_key]) for r in usable])) if usable else None,
                "n_random_pairs_median": float(np.median([float(r[random_key]) for r in usable])) if usable else None,
            }
            for predictor, xkey, covars in (
                ("icl", "icl_delta_200_300", ["satellite_delta_50_75", "redshift", "log_richness"]),
                ("satellite", "satellite_delta_50_75", ["icl_delta_200_300", "redshift", "log_richness"]),
            ):
                if len(usable) < 6:
                    stats = {"partial_r": float("nan"), "block_bootstrap_95_low": float("nan"),
                             "block_bootstrap_95_high": float("nan"), "bootstrap_replicates": args.bootstrap,
                             "bootstrap_valid": 0, "bootstrap_valid_fraction": 0.0,
                             "n_blocks": 0, "largest_block_n": 0}
                else:
                    idx += 1
                    stats = block_bootstrap(usable, xkey, ykey, covars, args.bootstrap, SEED + idx)
                metrics.append({"relationship": f"wl_{predictor}", "component": component,
                                "radial_bin": radial_bin, **support, **stats})

    args.output.mkdir(parents=True)
    joined_fields = ["cluster_id", "filter", "control_id", "icl_delta_200_300",
                     "satellite_delta_50_75", "redshift", "richness", "log_richness",
                     "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg"]
    joined_fields += [f"wl_{key}" for key in wl_rows[0] if key != "lens_id"]
    with (args.output / "triple_overlap_lenses.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=joined_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "metrics.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(metrics)
    summary = {
        "schema": "astrafek-icl-satellite-wl-triple-overlap-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "scope": "existing i2 lenses jointly supported by ICL, satellite profiles, and random-subtracted WL profiles",
        "n_triple_overlap_lenses": len(rows),
        "n_i2_satellite_supported": len(sat_rows),
        "n_i2_continuous_wl_lenses": len(wl_rows),
        "radial_edges_proper_mpc": [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0],
        "n_reported_relationships": len(metrics),
        "expected_overlap_ids_from_audit": sorted(expected_ids, key=int),
        "connected_blocks_5deg_recomputed_per_metric": True,
        "bootstrap": {"method": "resample connected target/control 5-degree blocks with replacement",
                      "replicates": args.bootstrap, "base_seed": SEED},
        "wl_support_rule": "per-bin lens included only if target and random source-pair counts are both positive and all modeled values finite",
        "wl_profile_units": "Msun Mpc^-2; profiles are precomputed random-subtracted target-minus-random DeltaSigma",
        "covariates": {"icl_vs_satellite": ["redshift", "ln(richness)"],
                       "wl_vs_icl": ["satellite ICL contrast", "redshift", "ln(richness)"],
                       "wl_vs_satellite": ["ICL contrast", "redshift", "ln(richness)"]},
        "limitations": ["selected, small common-support sample; inference is exploratory",
                        "5-degree equal-area connected cells are approximate sky blocks, not validated HSC regions",
                        "WL per-lens profiles inherit the saved estimator and its unresolved photo-z, shear, magnification, masking, PSF, and covariance systematics",
                        "bootstrap intervals are descriptive and do not establish causal or physical detection",
                        "no multiplicity-adjusted significance is claimed across radial bins or relationships"],
        "input_sha256": {str(path): hashes[name] for name, path in paths.items()},
        "output_sha256": {name: sha256(args.output / name) for name in
                          ("triple_overlap_lenses.tsv", "metrics.tsv")},
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("status", "science_measurement", "n_triple_overlap_lenses",
                                              "output_sha256")}, indent=2))


if __name__ == "__main__":
    main()
