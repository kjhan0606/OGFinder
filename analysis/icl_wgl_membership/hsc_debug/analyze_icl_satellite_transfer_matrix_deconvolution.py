#!/usr/bin/env python3
"""Explore pairwise inversion of the measured ICL/satellite transfer matrix.

The response coefficients come from target-only injections. Holding the real
matched random controls fixed gives the incremental response of the
target-minus-control ICL contrast. This diagnostic inverts each pair's measured
2x2 response matrix, then compares the inferred component association with the
observed association. The bootstrap is conditional on the measured response
matrices and does not propagate injection-response uncertainty, selection, or
unmodeled source Poisson, PSF-wing, and sky-transfer errors. The inferred values
are a sensitivity diagnostic, not corrected science measurements.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUT_DIR = DATA / "icl_satellite_injection_cross_talk_2026-11-06_02"
INPUT_SUMMARY = INPUT_DIR / "summary.json"
INPUT_ROWS = INPUT_DIR / "association_injection_cross_talk.tsv"
DEFAULT_OUTPUT = DATA / "icl_satellite_transfer_matrix_deconvolution_2026-11-06_01"
EXPECTED_INPUT_SCHEMA = "astrafek-icl-satellite-injection-cross-talk-v2"
BOOTSTRAP_REPLICAS = 10000
BOOTSTRAP_SEED = 20261107


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray,
              richness: np.ndarray) -> float:
    if len(x) < 6 or not all(np.isfinite(a).all() for a in (x, y, z, richness)):
        return float("nan")
    design = np.column_stack((np.ones(len(x)), z, richness))
    if np.linalg.matrix_rank(design) != 3:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    return float(rx @ ry / denom) if denom > 0 else float("nan")


def pct_interval(values: list[float]) -> list[float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        raise ValueError("no finite bootstrap values")
    return np.quantile(array, [0.025, 0.975]).tolist()


def summarize_filter(rows: list[dict[str, str]], band: str,
                     replicas: int, seed: int) -> tuple[dict, list[dict]]:
    selected = [row for row in rows if row["filter"] == band]
    if len(selected) < 6:
        raise ValueError(f"too few transfer-matrix pairs in {band}: {len(selected)}")
    blocks = sorted({row["connected_5deg_block"] for row in selected})
    if len(blocks) < 6:
        raise ValueError(f"too few connected blocks in {band}: {len(blocks)}")

    matrices, observed = [], []
    for row in selected:
        matrices.append([
            [float(row["icl_from_diffuse_response"]),
             float(row["icl_from_satellite_response"])],
            [float(row["satellite_from_diffuse_response"]),
             float(row["satellite_from_satellite_response"])],
        ])
        observed.append([float(row["icl_delta_njy_arcsec2"]),
                         float(row["satellite_delta_njy_arcsec2"])])
    matrices = np.asarray(matrices, dtype=float)
    observed = np.asarray(observed, dtype=float)
    conditions = np.linalg.cond(matrices)
    determinants = np.linalg.det(matrices)
    latent = np.empty_like(observed)
    for index, (matrix, values) in enumerate(zip(matrices, observed)):
        latent[index] = np.linalg.solve(matrix, values)

    redshift = np.asarray([float(row["redshift"]) for row in selected])
    log_richness = np.asarray([float(row["log_richness"]) for row in selected])
    raw_r = partial_r(observed[:, 0], observed[:, 1], redshift, log_richness)
    inverted_r = partial_r(latent[:, 0], latent[:, 1], redshift, log_richness)

    by_block = {block: np.asarray([i for i, row in enumerate(selected)
                                   if row["connected_5deg_block"] == block], dtype=int)
                for block in blocks}
    rng = np.random.default_rng(seed)
    raw_draws, inverted_draws, change_draws = [], [], []
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        index = np.concatenate([by_block[block] for block in chosen])
        rr = partial_r(observed[index, 0], observed[index, 1],
                       redshift[index], log_richness[index])
        ir = partial_r(latent[index, 0], latent[index, 1],
                       redshift[index], log_richness[index])
        if np.isfinite(rr) and np.isfinite(ir):
            raw_draws.append(rr)
            inverted_draws.append(ir)
            change_draws.append(ir - rr)
    if len(raw_draws) < 0.9 * replicas:
        raise ValueError(f"too few valid block-bootstrap replicates for {band}: {len(raw_draws)}")

    derived = []
    for row, matrix, det, cond, values in zip(selected, matrices, determinants,
                                               conditions, latent):
        derived.append({
            "cluster_id": row["cluster_id"], "filter": band,
            "connected_5deg_block": row["connected_5deg_block"],
            "icl_observed_njy_arcsec2": row["icl_delta_njy_arcsec2"],
            "satellite_observed_njy_arcsec2": row["satellite_delta_njy_arcsec2"],
            "response_m11_icl_from_diffuse": float(matrix[0, 0]),
            "response_m12_icl_from_satellite": float(matrix[0, 1]),
            "response_m21_satellite_from_diffuse": float(matrix[1, 0]),
            "response_m22_satellite_from_satellite": float(matrix[1, 1]),
            "response_determinant": float(det),
            "response_condition_number": float(cond),
            "inverted_diffuse_amplitude_njy_arcsec2": float(values[0]),
            "inverted_satellite_amplitude_njy_arcsec2": float(values[1]),
        })

    def amp_stats(column: int) -> dict:
        amp = latent[:, column]
        return {
            "q16_median_q84": np.quantile(amp, [0.16, 0.5, 0.84]).tolist(),
            "minimum": float(np.min(amp)), "maximum": float(np.max(amp)),
            "maximum_absolute": float(np.max(np.abs(amp))),
            "n_negative": int(np.count_nonzero(amp < 0)),
        }

    result = {
        "n_pairs": len(selected), "n_connected_5deg_blocks": len(blocks),
        "largest_block_pairs": int(max(len(index) for index in by_block.values())),
        "observed_partial_r_z_logrichness": raw_r,
        "pairwise_matrix_inverted_partial_r_z_logrichness": inverted_r,
        "block_bootstrap_observed_r_ci95": pct_interval(raw_draws),
        "block_bootstrap_inverted_r_ci95_conditional_on_response": pct_interval(inverted_draws),
        "block_bootstrap_inverted_minus_observed_r_ci95_conditional_on_response": pct_interval(change_draws),
        "bootstrap_replicas_requested": replicas, "bootstrap_replicas_valid": len(raw_draws),
        "bootstrap_seed": seed,
        "matrix_condition_q16_median_q84_max": [
            *np.quantile(conditions, [0.16, 0.5, 0.84]).tolist(), float(np.max(conditions))],
        "matrix_determinant_min_median_max": [
            float(np.min(determinants)), float(np.median(determinants)), float(np.max(determinants))],
        "inverted_diffuse_amplitude": amp_stats(0),
        "inverted_satellite_amplitude": amp_stats(1),
        "interpretation": (
            "The pairwise matrix inversion is an exploratory sensitivity only. "
            "Intervals condition on the measured injection matrices and omit "
            "response uncertainty, support-selection uncertainty, source Poisson, "
            "extended PSF wings, and upstream sky transfer. Large or negative "
            "amplitudes diagnose unstable/nonphysical inversions and are not corrected measurements."
        ),
    }
    return result, derived


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary = json.loads(INPUT_SUMMARY.read_text(encoding="utf-8"))
    if (summary.get("schema") != EXPECTED_INPUT_SCHEMA
            or summary.get("status") != "complete_exploratory_diagnostic"
            or summary.get("science_measurement") is not False
            or summary.get("n_pairs_with_valid_injection_transfer") != 58):
        raise ValueError("input transfer summary is incomplete or has unexpected scope")
    if summary.get("output_sha256", {}).get("association_injection_cross_talk.tsv") != sha256(INPUT_ROWS):
        raise ValueError("input transfer table hash differs from its summary")
    rows = read_rows(INPUT_ROWS)
    if len(rows) != 58 or len({(row["cluster_id"], row["filter"]) for row in rows}) != 58:
        raise ValueError("expected 58 unique audited injection/association pairs")
    results, derived_rows = {}, []
    for band_index, band in enumerate(("i", "i2")):
        results[band], derived = summarize_filter(rows, band, replicas, seed + band_index * 1000)
        derived_rows.extend(derived)

    output.mkdir(parents=True, exist_ok=False)
    details_path = output / "pairwise_deconvolution.tsv"
    with details_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(derived_rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(sorted(derived_rows, key=lambda row: (row["filter"], int(row["cluster_id"]))))
    result = {
        "schema": "astrafek-icl-satellite-transfer-matrix-deconvolution-v1",
        "status": "complete_exploratory_diagnostic",
        "science_measurement": False,
        "data_requested": False, "hsc_network_used": False,
        "matrix_definition": "[[ICL contrast response to diffuse, ICL contrast response to satellite], [satellite profile response to diffuse, satellite profile response to satellite]]",
        "target_only_injection_interpretation": "The target-side incremental response applies to the target-minus-control ICL contrast when the matched random control is held fixed.",
        "filters": results,
        "limitations": [
            "Only the selected 58-pair injection overlap is analyzed, not the complete 94-pair association sample.",
            "Bootstrap intervals condition on fixed measured transfer matrices and omit injection-response uncertainty.",
            "The inversion assumes a linear two-component model with zero unmodeled intercept and fixed observed masks/support.",
            "Source Poisson noise, coadd pixel weights, extended PSF wings, and upstream HSC sky transfer are omitted.",
            "Inverted amplitudes and correlations are stress-test diagnostics, not corrected physical measurements.",
        ],
        "input_sha256": {str(INPUT_SUMMARY): sha256(INPUT_SUMMARY), str(INPUT_ROWS): sha256(INPUT_ROWS)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {"pairwise_deconvolution.tsv": sha256(details_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    (output / "output_hashes.json").write_text(json.dumps({
        "summary_sha256": sha256(summary_path),
        "pairwise_deconvolution_sha256": sha256(details_path),
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "output": str(output), "filters": {
                          band: {"n": value["n_pairs"],
                                 "observed_r": value["observed_partial_r_z_logrichness"],
                                 "inverted_r": value["pairwise_matrix_inverted_partial_r_z_logrichness"],
                                 "inverted_r_ci95": value["block_bootstrap_inverted_r_ci95_conditional_on_response"],
                                 "matrix_condition": value["matrix_condition_q16_median_q84_max"]}
                          for band, value in results.items()}}, ensure_ascii=False, indent=2), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=BOOTSTRAP_REPLICAS)
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicates are required")
    run(args.output, args.replicas, args.seed)


if __name__ == "__main__":
    main()
