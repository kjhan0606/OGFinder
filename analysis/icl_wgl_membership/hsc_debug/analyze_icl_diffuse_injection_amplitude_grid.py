#!/usr/bin/env python3
"""Paired diffuse-ICL injection linearity and mask-response audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SAMPLE_MANIFEST = ROOT / "sample_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
PSF_PLAN_SHA256 = "1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_run(path: Path, expected_amplitude: float) -> tuple[dict, dict]:
    summary_path = path / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete" or summary.get("n_completed") != 137
            or summary.get("n_errors") != 0 or summary.get("science_measurement") is not False
            or summary.get("observational_validation_passed") is not False
            or summary.get("flux_covariance_available") is not False):
        raise ValueError(f"unexpected run status/scope: {summary_path}")
    injection = summary["injection"]
    if (not np.isclose(injection["diffuse_amplitude_at_100_njy_arcsec2"], expected_amplitude)
            or not np.isclose(injection["satellite_amplitude_at_50_per_galaxy_njy_arcsec2"], 1.0)
            or injection["diffuse_scale_kpc"] != 100.0
            or injection["satellite_scale_kpc"] != 15.0
            or injection["detection_rerun"] is not True
            or injection["psf_convolution_applied"] is not True
            or summary.get("supported_satellite_psf_request_approval", {}).get("approved_plan_sha256") != PSF_PLAN_SHA256):
        raise ValueError(f"injection configuration or PSF approval mismatch: {summary_path}")
    records = {}
    for item in summary["records"]:
        key = (str(item["cluster_id"]), item["filter"])
        if key in records:
            raise ValueError(f"duplicate frame in {summary_path}: {key}")
        result_path = Path(item["path"])
        if not result_path.is_file() or sha(result_path) != item["sha256"]:
            raise ValueError(f"missing or changed frame result: {result_path}")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if (result.get("status") != "complete" or result.get("science_measurement") is not False
                or str(result.get("cluster_id")) != key[0] or result.get("filter") != key[1]):
            raise ValueError(f"unexpected frame result scope or identity: {result_path}")
        records[key] = (result_path, result)
    if len(records) != 137:
        raise ValueError(f"expected 137 frame results in {summary_path}; found {len(records)}")
    return summary, records


def response(result: dict, name: str) -> float | None:
    value = result["transfer_response"].get(name)
    return None if value is None else float(value)


def finite(value) -> bool:
    return value is not None and np.isfinite(float(value))


def block_bootstrap_delta(rows: list[dict], low_tag: str, ref_tag: str,
                          replicas: int, seed: int) -> dict:
    by_block: dict[str, list[dict]] = {}
    for row in rows:
        if finite(row[f"response_{low_tag}"]) and finite(row[f"response_{ref_tag}"]):
            by_block.setdefault(row["spatial_block_5deg"], []).append(row)
    blocks = sorted(by_block)
    paired = [row for block in blocks for row in by_block[block]]
    deltas = np.asarray([float(row[f"response_{low_tag}"]) - float(row[f"response_{ref_tag}"])
                         for row in paired], dtype=float)
    if len(blocks) < 6 or not len(deltas):
        raise ValueError(f"insufficient matched spatial blocks for bootstrap: {len(blocks)}")
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    for i in range(replicas):
        picked = rng.choice(blocks, size=len(blocks), replace=True)
        values = np.concatenate([
            np.asarray([float(row[f"response_{low_tag}"]) - float(row[f"response_{ref_tag}"])
                        for row in by_block[block]], dtype=float)
            for block in picked
        ])
        draws[i] = np.median(values)
    return {
        "paired_frames": int(len(deltas)), "spatial_blocks": int(len(blocks)),
        "largest_block_frames": int(max(len(value) for value in by_block.values())),
        "median_delta": float(np.median(deltas)),
        "spatial_block_bootstrap_95pct_ci": np.quantile(draws, [0.025, 0.975]).tolist(),
        "replicas": replicas, "seed": seed,
    }


def quantiles(values: list[float]) -> dict:
    x = np.asarray(values, dtype=float)
    return {"n": int(len(x)), "median": float(np.median(x)),
            "q16": float(np.quantile(x, 0.16)), "q84": float(np.quantile(x, 0.84))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amp-0p5", type=Path, required=True)
    parser.add_argument("--amp-1p0", type=Path, required=True)
    parser.add_argument("--amp-2p0", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10000)
    args = parser.parse_args()
    if args.bootstrap < 1:
        raise ValueError("bootstrap count must be positive")
    paths = {0.5: args.amp_0p5.resolve(), 1.0: args.amp_1p0.resolve(),
             2.0: args.amp_2p0.resolve()}
    runs = {amp: load_run(path, amp) for amp, path in paths.items()}
    keys = set(runs[1.0][1])
    if any(set(runs[amp][1]) != keys for amp in runs):
        raise ValueError("the diffuse-amplitude runs do not contain identical frame sets")

    sample = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(x["id"]): x for x in sample["clusters"]}
    controls = {str(x["target_cluster_id"]): x for x in controls_doc["clusters"]}
    ordered = sorted(keys, key=lambda key: int(key[0]))
    target_coords, control_coords = [], []
    for ident, _ in ordered:
        if ident not in targets or ident not in controls:
            raise ValueError(f"missing target/control coordinate for {ident}")
        target_coords.append({"ra_deg": targets[ident]["ra_deg"],
                              "dec_deg": targets[ident]["dec_deg"]})
        control_coords.append({"ra_deg": controls[ident]["ra_deg"],
                               "dec_deg": controls[ident]["dec_deg"]})
    _, labels, _, _, n_blocks, largest_block = pair_block_labels(
        target_coords, control_coords, np.ones(len(ordered), dtype=bool))
    block_by_key = dict(zip(ordered, map(str, labels)))

    rows = []
    for key in ordered:
        samples = {amp: runs[amp][1][key][1] for amp in runs}
        reference = samples[1.0]
        for amp, result in samples.items():
            if (result["z"] != reference["z"] or result["filter"] != reference["filter"]
                    or result["icl_profile_bin_kpc"] != reference["icl_profile_bin_kpc"]
                    or result["satellite_profile_bin_kpc"] != reference["satellite_profile_bin_kpc"]
                    or result["cluster_radius_selection_kpc"] != reference["cluster_radius_selection_kpc"]):
                raise ValueError(f"frame configuration mismatch for {key}")
            inj, ref_inj = result["injection"], reference["injection"]
            for field in ("satellite_amplitude_at_50_per_galaxy_njy_arcsec2",
                          "diffuse_scale_kpc", "satellite_scale_kpc",
                          "n_injected_satellite_centers", "psf_file_count", "psf_model",
                          "unresolved_bcg_policy", "source_detection_and_masks_rerun"):
                if inj.get(field) != ref_inj.get(field):
                    raise ValueError(f"injection configuration mismatch in {field} for {key}")
            def image_hash(doc):
                return next((v for name, v in doc["input_sha256"].items()
                             if name.endswith(f"/images/cluster_{key[0]}.fits")), None)
            def psf_hashes(doc):
                return sorted(v for name, v in doc["input_sha256"].items()
                              if "/psf_" in name and name.endswith("_HSC-I.fits"))
            if image_hash(result) != image_hash(reference) or psf_hashes(result) != psf_hashes(reference):
                raise ValueError(f"image or position-dependent PSF content changed for {key}")
            if result["redetection_by_scenario"]["baseline"] != reference["redetection_by_scenario"]["baseline"]:
                raise ValueError(f"no-injection baseline changed for {key}")

        row = {"cluster_id": key[0], "filter": key[1], "z": reference["z"],
               "spatial_block_5deg": block_by_key[key]}
        for amp, result in samples.items():
            tag = str(amp).replace(".", "p")
            row[f"response_{tag}"] = response(result, "icl_response_to_diffuse_unit_injection")
            row[f"plane_corrected_response_{tag}"] = response(
                result, "icl_plane_corrected_response_to_diffuse_unit_injection")
            for scenario in ("icl_only", "both"):
                measures = result["redetection_by_scenario"][scenario]
                baseline = result["redetection_by_scenario"]["baseline"]
                row[f"detections_{scenario}_minus_baseline_{tag}"] = int(measures["n_detections"] - baseline["n_detections"])
                row[f"mask_pixels_{scenario}_minus_baseline_{tag}"] = int(measures["icl_mask_pixels"] - baseline["icl_mask_pixels"])
        rows.append(row)

    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)
    pair_path = args.output_dir / "paired_amplitude_frames.tsv"
    with pair_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "schema": "astrafek-icl-diffuse-injection-amplitude-grid-v1",
        "status": "complete_diagnostic", "science_measurement": False,
        "n_frames": len(rows), "frames_by_filter": {
            band: sum(row["filter"] == band for row in rows)
            for band in sorted({row["filter"] for row in rows})},
        "diffuse_amplitudes_njy_arcsec2_at_100kpc": [0.5, 1.0, 2.0],
        "fixed_satellite_amplitude_njy_arcsec2_per_galaxy_at_50kpc": 1.0,
        "image_hashes_identical": True, "position_psf_hashes_identical": True,
        "baseline_detection_masks_identical": True,
        "spatial_block_definition": "connected pairs sharing either endpoint's equal-area 5-degree sky cell",
        "spatial_blocks": int(n_blocks), "largest_spatial_block_frames": int(largest_block),
        "metrics_by_filter": {}, "bootstrap_replicas": args.bootstrap,
        "input_summaries_sha256": {str(amp): sha(path / "summary.json") for amp, path in paths.items()},
        "input_manifests_sha256": {str(SAMPLE_MANIFEST): sha(SAMPLE_MANIFEST),
                                   str(CONTROL_MANIFEST): sha(CONTROL_MANIFEST)},
        "paired_table_sha256": None,
        "limitations": [
            "This is a selected-frame injection/recovery transfer diagnostic, not a science stack or completeness correction.",
            "Injected-source Poisson noise, coadd pixel weights, upstream HSC sky transfer, and extended PSF wings are not modeled.",
            "The 5-degree connected cells are an approximate spatial partition, not validated independent HSC fields.",
            "The response is conditional on the fixed observed baseline sample; no population-level inference is made.",
        ],
    }
    for band in [*sorted({row["filter"] for row in rows}), "all"]:
        selected = [row for row in rows if band == "all" or row["filter"] == band]
        metrics = {}
        for metric in ("response", "plane_corrected_response"):
            metrics[metric] = {
                str(amp): quantiles([float(row[f"{metric}_{str(amp).replace('.', 'p')}"])
                                     for row in selected
                                     if finite(row[f"{metric}_{str(amp).replace('.', 'p')}"])])
                for amp in (0.5, 1.0, 2.0)}
            metrics[metric]["paired_delta_vs_1"] = {
                str(amp): block_bootstrap_delta(selected, str(amp).replace(".", "p"), "1p0",
                                                args.bootstrap, 20261116 + int(amp * 10))
                for amp in (0.5, 2.0)}
        metrics["mask_and_detection_changes"] = {}
        for amp in (0.5, 1.0, 2.0):
            tag = str(amp).replace(".", "p")
            metrics["mask_and_detection_changes"][str(amp)] = {}
            for scenario in ("icl_only", "both"):
                for kind in ("detections", "mask_pixels"):
                    name = f"{kind}_{scenario}_minus_baseline_{tag}"
                    values = [row[name] for row in selected]
                    metrics["mask_and_detection_changes"][str(amp)][name] = {
                        **quantiles(values), "n_frames_changed": int(sum(value != 0 for value in values))}
        summary["metrics_by_filter"][band] = metrics
    summary["paired_table_sha256"] = sha(pair_path)
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = ["# Diffuse ICL injection amplitude response audit", "",
             f"The same {len(rows)} PSF-convolved HSC frames were reprocessed at diffuse amplitudes 0.5, 1.0, and 2.0 nJy arcsec⁻² at 100 kpc, with a fixed satellite injection and regenerated SEP masks. All outputs retain `science_measurement=false`.", "",
             "| Filter | Input amplitude | ICL annulus response per unit input (median, q16–q84) | Paired Δ vs 1.0 (spatial-block 95% interval) |", "|---|---:|---:|---:|"]
    for band in sorted({row["filter"] for row in rows}):
        metrics = summary["metrics_by_filter"][band]
        for amp in (0.5, 1.0, 2.0):
            item = metrics["response"][str(amp)]
            delta = metrics["response"]["paired_delta_vs_1"].get(str(amp))
            delta_text = "—" if delta is None else f"{delta['median_delta']:.3f} [{delta['spatial_block_bootstrap_95pct_ci'][0]:.3f}, {delta['spatial_block_bootstrap_95pct_ci'][1]:.3f}]"
            lines.append(f"| {band} | {amp:.1f} | {item['median']:.3f} ({item['q16']:.3f}–{item['q84']:.3f}) | {delta_text} |")
    lines += ["", "The paired response is approximately stable across this two-fold amplitude range if the spatial-block intervals include zero. This does not establish absolute recovery: median response is around one quarter of injected light, and the simulation omits source Poisson noise, coadd weights, upstream HSC sky processing, and extended PSF wings. It is a pipeline-transfer sensitivity, not an astrophysical ICL measurement or a correction factor.", ""]
    (args.output_dir / "assessment.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir), "frames": len(rows),
                      "summary_sha256": sha(summary_path), "paired_table_sha256": sha(pair_path),
                      "assessment_sha256": sha(args.output_dir / "assessment.md")}, indent=2))


if __name__ == "__main__":
    main()
