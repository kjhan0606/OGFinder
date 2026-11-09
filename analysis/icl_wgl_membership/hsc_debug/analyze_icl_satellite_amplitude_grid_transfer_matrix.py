#!/usr/bin/env python3
"""Propagate the 0.5/1/2 satellite-injection grid through transfer-matrix inversion.

This offline diagnostic uses the same 137 PSF-supported target images at all
three satellite amplitudes. It joins the exact 58-pair overlap with the
observed ICL/satellite association sample, verifies that amplitude-1 responses
match the audited cross-talk table, then repeats pairwise matrix inversion at
each injection amplitude. Spatial-block bootstrap intervals condition on the
measured transfer matrices. The analysis does not propagate response
uncertainty and remains a sensitivity diagnostic, not a science correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RUNS = {
    0.5: DATA / "joint_icl_satellite_psf_redetection_icl_supported224_amp0p5_2026-11-06_01",
    1.0: DATA / "joint_icl_satellite_psf_redetection_icl_supported224_2026-11-06_01",
    2.0: DATA / "joint_icl_satellite_psf_redetection_icl_supported224_amp2p0_2026-11-06_01",
}
ASSOCIATION_DIR = DATA / "icl_satellite_injection_cross_talk_2026-11-06_02"
ASSOCIATION_SUMMARY = ASSOCIATION_DIR / "summary.json"
ASSOCIATION_ROWS = ASSOCIATION_DIR / "association_injection_cross_talk.tsv"
DEFAULT_OUTPUT = DATA / "icl_satellite_amplitude_grid_transfer_matrix_2026-11-06_01"
BOOTSTRAP_REPLICAS = 10000
BOOTSTRAP_SEED = 20261108
RESPONSE_COLUMNS = {
    "a_icl_from_diffuse": "icl_response_to_diffuse_unit_injection",
    "b_icl_from_satellite": "icl_response_to_all_satellite_envelopes",
    "c_satellite_from_diffuse": "satellite_response_to_smooth_diffuse",
    "d_satellite_from_satellite": "satellite_response_to_satellite_unit_injection",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray,
              richness: np.ndarray) -> float:
    if len(x) < 6 or not all(np.isfinite(v).all() for v in (x, y, z, richness)):
        return float("nan")
    design = np.column_stack((np.ones(len(x)), z, richness))
    if np.linalg.matrix_rank(design) != 3:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    return float(rx @ ry / denom) if denom > 0 else float("nan")


def bootstrap_median_difference(rows: list[dict], field: str,
                                band: str, amplitude: float,
                                replicas: int, seed: int) -> dict:
    selected = [row for row in rows if row["filter"] == band]
    blocks = sorted({row["connected_5deg_block"] for row in selected})
    by_block = {block: np.asarray([i for i, row in enumerate(selected)
                                   if row["connected_5deg_block"] == block], dtype=int)
                for block in blocks}
    values = np.asarray([float(row[field]) for row in selected])
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        index = np.concatenate([by_block[block] for block in chosen])
        draws.append(float(np.median(values[index])))
    return {
        "n_pairs": len(selected), "n_connected_5deg_blocks": len(blocks),
        "median_difference": float(np.median(values)),
        "block_bootstrap_ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "replicas": replicas, "seed": seed,
    }


def validate_batches() -> tuple[dict, dict[float, dict[tuple[str, str], dict]]]:
    loaded = {}
    batch_summaries = {}
    reference_keys = None
    reference_image_hashes = None
    reference_inputs = None
    reference_baselines = None
    for amp, directory in RUNS.items():
        summary_path = directory / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if (summary.get("status") != "complete"
                or summary.get("science_measurement") is not False
                or summary.get("n_completed") != 137
                or summary.get("n_errors") != 0
                or summary.get("injection", {}).get("detection_rerun") is not True
                or summary.get("injection", {}).get("psf_convolution_applied") is not True
                or summary.get("injection", {}).get("satellite_amplitude_at_50_per_galaxy_njy_arcsec2") != amp
                or summary.get("injection", {}).get("diffuse_amplitude_at_100_njy_arcsec2") != 1.0):
            raise ValueError(f"injection batch has wrong status/scope/amplitude: {directory}")
        records = {}
        image_hashes = {}
        for record in summary.get("records", []):
            key = (str(record["cluster_id"]), str(record["filter"]))
            path = Path(record["path"])
            if not path.is_file() or sha256(path) != record.get("sha256"):
                raise ValueError(f"injection record missing or changed: {key} at {amp}")
            detail = json.loads(path.read_text(encoding="utf-8"))
            if (detail.get("science_measurement") is not False
                    or detail.get("injection", {}).get("satellite_amplitude_at_50_per_galaxy_njy_arcsec2") != amp
                    or detail.get("cluster_id") != key[0]
                    or detail.get("filter") != key[1]):
                raise ValueError(f"per-frame injection provenance mismatch: {key} at {amp}")
            image_paths = [name for name in detail.get("input_sha256", {})
                           if name.endswith(f"/images/cluster_{key[0]}.fits")]
            if len(image_paths) != 1:
                raise ValueError(f"missing unique image hash for {key} at {amp}")
            image_hashes[key] = detail["input_sha256"][image_paths[0]]
            records[key] = {"record": record, "detail": detail}
        if len(records) != 137:
            raise ValueError(f"expected 137 unique pair records at amplitude {amp}")
        if reference_keys is None:
            reference_keys = set(records)
            reference_image_hashes = image_hashes
            reference_inputs = {key: records[key]["detail"]["input_sha256"] for key in records}
            reference_baselines = {key: records[key]["detail"]["redetection_by_scenario"]["baseline"]
                                  for key in records}
        else:
            if set(records) != reference_keys or image_hashes != reference_image_hashes:
                raise ValueError("injection amplitudes do not use the exact same pair/image set")
            for key in records:
                if records[key]["detail"]["input_sha256"] != reference_inputs[key]:
                    raise ValueError(f"input hashes differ across amplitudes for {key}")
                if (records[key]["detail"]["redetection_by_scenario"]["baseline"]
                        != reference_baselines[key]):
                    raise ValueError(f"baseline redetection differs across amplitudes for {key}")
        loaded[amp] = records
        batch_summaries[amp] = {"path": str(summary_path), "sha256": sha256(summary_path),
                                "n_records": len(records), "injection": summary["injection"]}
    return batch_summaries, loaded


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    batch_summaries, batches = validate_batches()
    transfer_summary = json.loads(ASSOCIATION_SUMMARY.read_text(encoding="utf-8"))
    if (transfer_summary.get("schema") != "astrafek-icl-satellite-injection-cross-talk-v2"
            or transfer_summary.get("status") != "complete_exploratory_diagnostic"
            or transfer_summary.get("science_measurement") is not False
            or transfer_summary.get("output_sha256", {}).get("association_injection_cross_talk.tsv")
            != sha256(ASSOCIATION_ROWS)):
        raise ValueError("corrected association transfer table failed provenance checks")
    association = read_tsv(ASSOCIATION_ROWS)
    if len(association) != 58:
        raise ValueError(f"expected 58 common injection/association rows; got {len(association)}")
    association_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in association}
    if len(association_by_key) != 58:
        raise ValueError("duplicate association transfer key")
    for key in association_by_key:
        if not all(key in batch for batch in batches.values()):
            raise ValueError(f"a 58-pair association row is missing from one amplitude batch: {key}")

    # At unit amplitude, the current 137-frame sweep must reproduce the audited
    # cross-talk response coefficients on every shared pair.
    comparison = {
        "icl_from_diffuse_response": "icl_response_to_diffuse_unit_injection",
        "icl_from_satellite_response": "icl_response_to_all_satellite_envelopes",
        "satellite_from_diffuse_response": "satellite_response_to_smooth_diffuse",
        "satellite_from_satellite_response": "satellite_response_to_satellite_unit_injection",
    }
    maximum_unit_amplitude_difference = {column: 0.0 for column in comparison}
    responses_by_amp = {amp: {} for amp in batches}
    excluded_by_amp = {amp: {} for amp in batches}
    for key in association_by_key:
        for amp, batch in batches.items():
            detail = batch[key]["detail"]
            response = detail.get("transfer_response", {})
            values = [response.get(RESPONSE_COLUMNS[name]) for name in RESPONSE_COLUMNS]
            reason = None
            if detail.get("satellite_support_gaps"):
                reason = "satellite_support_gap"
            elif any(value is None or not np.isfinite(float(value)) for value in values):
                reason = "nonfinite_transfer_response"
            if reason is not None:
                excluded_by_amp[amp][key] = reason
                continue
            responses_by_amp[amp][key] = [float(value) for value in values]
    for key, observed in association_by_key.items():
        if key not in responses_by_amp[1.0]:
            raise ValueError(f"unit-amplitude response missing for audited association pair {key}")
        unit_response = responses_by_amp[1.0][key]
        for assoc_column, index in zip(comparison, range(4)):
            difference = abs(float(observed[assoc_column]) - unit_response[index])
            maximum_unit_amplitude_difference[assoc_column] = max(
                maximum_unit_amplitude_difference[assoc_column], difference)
            if difference > 1e-12:
                raise ValueError(f"unit-amplitude response does not reproduce cross-talk table: {key}, {assoc_column}")
    common_keys = set.intersection(*(set(responses_by_amp[amp]) for amp in batches))
    if not common_keys:
        raise ValueError("no common response support across all three satellite amplitudes")
    support_counts = {
        f"{amp:g}": {
            "n_valid_association_pairs": len(responses_by_amp[amp]),
            "n_valid_by_filter": {band: sum(key[1] == band for key in responses_by_amp[amp])
                                  for band in ("i", "i2")},
            "excluded_reasons": {reason: sum(value == reason for value in excluded_by_amp[amp].values())
                                 for reason in sorted(set(excluded_by_amp[amp].values()))},
        } for amp in batches
    }
    common_by_filter = {band: sum(key[1] == band for key in common_keys) for band in ("i", "i2")}
    records = []
    for key in sorted(common_keys, key=lambda item: (item[1], int(item[0]))):
        observed = association_by_key[key]
        vector = np.asarray([float(observed["icl_delta_njy_arcsec2"]),
                             float(observed["satellite_delta_njy_arcsec2"])])
        for amp, response_by_key in responses_by_amp.items():
            a, b, c, d = response_by_key[key]
            matrix = np.asarray([[a, b], [c, d]], dtype=float)
            latent = np.linalg.solve(matrix, vector)
            records.append({
                "cluster_id": key[0], "filter": key[1], "amplitude": amp,
                "connected_5deg_block": observed["connected_5deg_block"],
                "redshift": float(observed["redshift"]),
                "log_richness": float(observed["log_richness"]),
                "icl_observed_njy_arcsec2": vector[0],
                "satellite_observed_njy_arcsec2": vector[1],
                "a_icl_from_diffuse": a, "b_icl_from_satellite": b,
                "c_satellite_from_diffuse": c, "d_satellite_from_satellite": d,
                "matrix_determinant": float(np.linalg.det(matrix)),
                "matrix_condition_number": float(np.linalg.cond(matrix)),
                "inverted_diffuse_amplitude_njy_arcsec2": float(latent[0]),
                "inverted_satellite_amplitude_njy_arcsec2": float(latent[1]),
            })

    filter_results = {}
    for band_index, band in enumerate(("i", "i2")):
        observed_rows = [association_by_key[key] for key in common_keys if key[1] == band]
        observed_rows.sort(key=lambda row: int(row["cluster_id"]))
        if len(observed_rows) < 6:
            raise ValueError(f"too few association rows in {band}")
        z = np.asarray([float(row["redshift"]) for row in observed_rows])
        richness = np.asarray([float(row["log_richness"]) for row in observed_rows])
        blocks = sorted({row["connected_5deg_block"] for row in observed_rows})
        by_block = {block: np.asarray([i for i, row in enumerate(observed_rows)
                                       if row["connected_5deg_block"] == block], dtype=int)
                    for block in blocks}
        if len(blocks) < 6:
            raise ValueError(f"too few connected sky blocks in {band}")
        observed_vector = np.asarray([[float(row["icl_delta_njy_arcsec2"]),
                                       float(row["satellite_delta_njy_arcsec2"])]
                                      for row in observed_rows])
        result = {"n_pairs": len(observed_rows), "n_connected_5deg_blocks": len(blocks),
                  "observed_partial_r_z_logrichness": partial_r(
                      np.asarray([float(row["icl_delta_njy_arcsec2"]) for row in observed_rows]),
                      np.asarray([float(row["satellite_delta_njy_arcsec2"]) for row in observed_rows]), z, richness),
                  "maximum_unit_amplitude_response_difference": maximum_unit_amplitude_difference}
        for amp_index, amp in enumerate((0.5, 1.0, 2.0)):
            amp_rows = [row for row in records if row["filter"] == band and row["amplitude"] == amp]
            amp_rows.sort(key=lambda row: (int(row["cluster_id"])))
            selected = observed_rows
            if [row["cluster_id"] for row in amp_rows] != [row["cluster_id"] for row in selected]:
                raise ValueError(f"amplitude row order/key mismatch for {band} {amp}")
            latent = np.asarray([[row["inverted_diffuse_amplitude_njy_arcsec2"],
                                 row["inverted_satellite_amplitude_njy_arcsec2"]]
                                for row in amp_rows])
            matrices_cond = np.asarray([row["matrix_condition_number"] for row in amp_rows])
            raw_r = partial_r(observed_vector[:, 0], observed_vector[:, 1], z, richness)
            deconv_r = partial_r(latent[:, 0], latent[:, 1], z, richness)
            rng = np.random.default_rng(seed + 1000 * band_index + int(amp * 100))
            raw_draws, inverse_draws, delta_draws = [], [], []
            for _ in range(replicas):
                chosen = rng.choice(blocks, size=len(blocks), replace=True)
                index = np.concatenate([by_block[block] for block in chosen])
                rr = partial_r(observed_vector[index, 0], observed_vector[index, 1],
                               z[index], richness[index])
                ir = partial_r(latent[index, 0], latent[index, 1], z[index], richness[index])
                if np.isfinite(rr) and np.isfinite(ir):
                    raw_draws.append(rr); inverse_draws.append(ir); delta_draws.append(ir-rr)
            if len(inverse_draws) < 0.9 * replicas:
                raise ValueError(f"too few finite bootstrap draws for {band} amp={amp}")
            result[str(amp)] = {
                "observed_partial_r": raw_r,
                "inverted_partial_r": deconv_r,
                "conditional_inverted_r_block_ci95": np.quantile(inverse_draws, [0.025, 0.975]).tolist(),
                "conditional_inverted_minus_observed_ci95": np.quantile(delta_draws, [0.025, 0.975]).tolist(),
                "median_condition_number": float(np.median(matrices_cond)),
                "maximum_condition_number": float(np.max(matrices_cond)),
                "negative_diffuse_amplitudes": int(np.count_nonzero(latent[:, 0] < 0)),
                "negative_satellite_amplitudes": int(np.count_nonzero(latent[:, 1] < 0)),
                "diffuse_amplitude_q16_median_q84": np.quantile(latent[:, 0], [0.16, 0.5, 0.84]).tolist(),
                "satellite_amplitude_q16_median_q84": np.quantile(latent[:, 1], [0.16, 0.5, 0.84]).tolist(),
                "diffuse_maximum_absolute": float(np.max(np.abs(latent[:, 0]))),
                "satellite_maximum_absolute": float(np.max(np.abs(latent[:, 1]))),
                "bootstrap_requested": replicas, "bootstrap_valid": len(inverse_draws),
            }
        result["pairwise_coefficient_changes_vs_amplitude_1"] = {}
        for coefficient in ("b_icl_from_satellite", "d_satellite_from_satellite"):
            result["pairwise_coefficient_changes_vs_amplitude_1"][coefficient] = {}
            for amp_index, amp in enumerate((0.5, 2.0)):
                by_key = {(row["cluster_id"], row["filter"], row["amplitude"]): row for row in records}
                differences = []
                for assoc in observed_rows:
                    key = (assoc["cluster_id"], assoc["filter"])
                    differences.append({
                        "cluster_id": key[0], "filter": band,
                        "connected_5deg_block": assoc["connected_5deg_block"],
                        "delta": (by_key[(key[0], key[1], amp)][coefficient]
                                  - by_key[(key[0], key[1], 1.0)][coefficient]),
                    })
                result["pairwise_coefficient_changes_vs_amplitude_1"][coefficient][str(amp)] = \
                    bootstrap_median_difference(differences, "delta", band, amp, replicas,
                                                seed + 5000 + 1000 * band_index + int(amp * 100))
        filter_results[band] = result

    output.mkdir(parents=True, exist_ok=False)
    rows_path = output / "amplitude_pairwise_deconvolution.tsv"
    with rows_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda row: (row["filter"], float(row["amplitude"]), int(row["cluster_id"]))))
    result = {
        "schema": "astrafek-icl-satellite-amplitude-transfer-matrix-v1",
        "status": "complete_exploratory_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False, "data_requested": False,
        "amplitudes_njy_arcsec2": [0.5, 1.0, 2.0],
        "n_source_frames_per_amplitude": 137,
        "n_observed_pairs_parent": len(association_by_key),
        "n_observed_pairs_with_complete_amplitude_matrix": len(common_keys),
        "common_pairs_by_filter": common_by_filter,
        "valid_support_by_amplitude": support_counts,
        "filter_results": filter_results,
        "matrix_definition": "[[ICL target-minus-control contrast response to diffuse, ICL contrast response to satellite], [satellite profile response to diffuse, satellite response to satellite]]",
        "limitations": [
            "All amplitudes use the same selected target images, positions, and PSF files; the 58-pair overlap remains selected from the 94-pair association sample.",
            "The 1.0-amplitude responses reproduce the corrected cross-talk table exactly; response changes at 0.5 and 2.0 arise under injection-dependent redetection/masking.",
            "Spatial-block bootstrap intervals condition on measured response matrices and omit injection-response uncertainty.",
            "The inversion assumes a zero-intercept linear two-component model; negative and extreme amplitudes indicate instability, not physical measurements.",
            "Injected source Poisson noise, coadd pixel weights, extended PSF wings, and upstream HSC sky transfer remain unmodeled.",
        ],
        "input_sha256": {
            **{f"batch_{amp:g}_summary": summary["sha256"] for amp, summary in batch_summaries.items()},
            "association_summary": sha256(ASSOCIATION_SUMMARY),
            "association_rows": sha256(ASSOCIATION_ROWS),
        },
        "batch_summaries": {f"{amp:g}": summary for amp, summary in batch_summaries.items()},
        "maximum_unit_amplitude_response_difference": maximum_unit_amplitude_difference,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {"amplitude_pairwise_deconvolution.tsv": sha256(rows_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    (output / "output_hashes.json").write_text(json.dumps({
        "summary_sha256": sha256(summary_path),
        "amplitude_pairwise_deconvolution_sha256": sha256(rows_path),
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "output": str(output),
                      "filter_results": {band: {
                          amp: {"inverted_r": values[amp]["inverted_partial_r"],
                                "conditional_ci95": values[amp]["conditional_inverted_r_block_ci95"],
                                "median_condition": values[amp]["median_condition_number"],
                                "negative_components": [values[amp]["negative_diffuse_amplitudes"],
                                                        values[amp]["negative_satellite_amplitudes"]]}
                          for amp in ("0.5", "1.0", "2.0")}
                          for band, values in filter_results.items()}}, ensure_ascii=False, indent=2), flush=True)
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
