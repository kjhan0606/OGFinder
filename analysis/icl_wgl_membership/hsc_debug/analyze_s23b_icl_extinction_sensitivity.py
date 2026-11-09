#!/usr/bin/env python3
"""Measure the effect of nominal HSC SFD extinction correction on the S23B stack.

This is a paired-profile sensitivity analysis. It applies the hscPipe i/i2
extinction coefficients separately to each target and its matched random
control, then estimates the mean profile and inner-minus-outer contrasts with
connected target/control sky-block bootstrap intervals. It reads saved radial
tables and metadata only; no image pixels or network data are used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_s23b_icl_spatial_blocks import connected_pair_blocks, largest_gap_origin


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
CONTROL_MANIFEST = DATA / "s23b_icl_control_image_plan_2026-11-07_02/control_manifest.json"
STACK_SUMMARY = DATA / "s23b_icl_target_control_stack_2026-11-07_01/summary.json"
PAIRED_PROFILES = DATA / "s23b_icl_target_control_stack_2026-11-07_01/paired_radial_profiles.tsv"
DEFAULT_OUTPUT = DATA / "s23b_icl_extinction_sensitivity_2026-11-07_04"
SPATIAL_SUMMARY = DATA / "s23b_icl_spatial_blocks_2026-11-07_02/summary.json"

EXPECTED_HASHES = {
    TARGET_MANIFEST: "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44",
    CONTROL_MANIFEST: "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df",
    STACK_SUMMARY: "c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b",
    PAIRED_PROFILES: "285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43",
    SPATIAL_SUMMARY: "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc",
}

# hscPipe extinctionCoeffs.py values, also used by the existing S23B HSC
# ancillary-preparation path. They are a nominal filter correction here.
COEFFICIENTS = {"i": 1.641, "i2": 1.617}
COEFFICIENT_SOURCE = "https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7_e/colorterms.html"
EXPECTED_PAIRS = {"i": 241, "i2": 781}
PROFILE_WINDOWS = ((0, 50), (50, 100), (100, 300), (50, 300),
                   (200, 300), (300, 800), (800, 1000))
PLOT_WINDOWS = ((0, 50), (50, 100), (100, 300), (300, 800), (800, 1000))
RELATIVE_WINDOWS = ((50, 300), (100, 300), (200, 300))
OUTER_WINDOW = (300, 800)
BLOCK_WIDTHS_DEG = (2.5, 5.0, 10.0)
BLOCK_PHASES = (0.0, 0.5)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def annulus_mean(profile: list[dict], low: int, high: int,
                 target_factor: float, control_factor: float) -> tuple[float | None, float | None, int]:
    selected = [item for item in profile
                if float(item["r_min_kpc"]) >= low and float(item["r_max_kpc"]) <= high]
    if not selected:
        raise ValueError(f"profile does not contain the {low}-{high} kpc window")
    raw_sum = corrected_sum = 0.0
    n_pixels = 0
    for item in selected:
        count = int(item["n_valid_pixels"])
        if count <= 0:
            continue
        target = float(item["target_mean_njy_arcsec2"])
        control = float(item["control_mean_njy_arcsec2"])
        raw = float(item["delta_mean_njy_arcsec2"])
        if not all(math.isfinite(x) for x in (target, control, raw)):
            raise ValueError("nonfinite pixel-supported annular profile")
        if not math.isclose(raw, target - control, rel_tol=0.0, abs_tol=2e-7):
            raise ValueError("saved raw target-minus-control profile does not reproduce")
        raw_sum += raw * count
        corrected_sum += (target * target_factor - control * control_factor) * count
        n_pixels += count
    if not n_pixels:
        return None, None, 0
    return raw_sum / n_pixels, corrected_sum / n_pixels, n_pixels


def prepare_pairs(target_doc: dict, control_doc: dict, profile_rows: list[dict]) -> dict[str, list[dict]]:
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    if len(targets) != 1043 or len(controls) != 1041:
        raise ValueError("pinned S23B target/control manifest row counts changed")
    if len(controls) != len(control_doc["clusters"]):
        raise ValueError("duplicate target-to-control mapping")

    result = {band: [] for band in ("i", "i2")}
    seen = set()
    for row in profile_rows:
        target_id, control_id = str(row["target_id"]), str(row["control_id"])
        band = str(row["filter"])
        key = (target_id, band)
        if band not in result or key in seen:
            raise ValueError(f"invalid or duplicate target/filter row: {key}")
        seen.add(key)
        target, control = targets.get(target_id), controls.get(target_id)
        if target is None or control is None:
            raise ValueError(f"missing target or matched control for {target_id}")
        if (str(control["id"]) != control_id
                or float(row["ra_deg"]) != float(target["ra_deg"])
                or float(row["dec_deg"]) != float(target["dec_deg"])
                or str(control["target_cluster_id"]) != target_id):
            raise ValueError(f"target/control identity or coordinates differ for {target_id}")
        ebv_target = float(control["target_ebv_sfd"])
        ebv_control = float(control["candidate_ebv_sfd"])
        delta_ebv = float(control["delta_ebv_sfd"])
        if (not all(math.isfinite(x) and x >= 0 for x in (ebv_target, ebv_control))
                or not math.isclose(delta_ebv, ebv_control - ebv_target,
                                    rel_tol=0.0, abs_tol=1e-12)):
            raise ValueError(f"invalid SFD reddening values for {target_id}")
        coefficient = COEFFICIENTS[band]
        target_factor = 10.0 ** (0.4 * coefficient * ebv_target)
        control_factor = 10.0 ** (0.4 * coefficient * ebv_control)
        profile = json.loads(row["profile_json"])
        record = {
            "target_id": target_id,
            "control_id": control_id,
            "filter": band,
            "ra_deg": float(row["ra_deg"]),
            "dec_deg": float(row["dec_deg"]),
            "z": float(row["z"]),
            "target_ebv_sfd": ebv_target,
            "control_ebv_sfd": ebv_control,
            "delta_ebv_sfd_control_minus_target": delta_ebv,
            "extinction_coefficient_mag_per_ebv": coefficient,
            "target_flux_factor": target_factor,
            "control_flux_factor": control_factor,
            "windows": {},
        }
        for low, high in PROFILE_WINDOWS:
            raw, corrected, pixels = annulus_mean(
                profile, low, high, target_factor, control_factor)
            record["windows"][(low, high)] = {
                "raw": raw, "corrected": corrected, "n_common_pixels": pixels,
            }
        for low, high in RELATIVE_WINDOWS:
            inner = record["windows"][(low, high)]
            outer = record["windows"][OUTER_WINDOW]
            if (inner["raw"] is None or outer["raw"] is None
                    or inner["corrected"] is None or outer["corrected"] is None):
                raw_contrast = corrected_contrast = None
            else:
                raw_contrast = inner["raw"] - outer["raw"]
                corrected_contrast = inner["corrected"] - outer["corrected"]
            record["windows"][(low, high, "relative")] = {
                "raw": raw_contrast,
                "corrected": corrected_contrast,
                "n_common_pixels": min(inner["n_common_pixels"], outer["n_common_pixels"]),
            }
        result[band].append(record)

    if {band: len(rows) for band, rows in result.items()} != EXPECTED_PAIRS:
        raise ValueError("filter-specific pair counts differ from the pinned exploratory stack")
    return result


def pair_metric_rows(pair_data: dict[str, list[dict]]) -> list[dict]:
    rows = []
    for band, pairs in pair_data.items():
        for pair in pairs:
            for low, high in PROFILE_WINDOWS:
                value = pair["windows"][(low, high)]
                rows.append({
                    "target_id": pair["target_id"], "control_id": pair["control_id"],
                    "filter": band, "metric": "profile", "r_min_kpc": low,
                    "r_max_kpc": high, "ra_deg": pair["ra_deg"], "dec_deg": pair["dec_deg"],
                    "z": pair["z"], "target_ebv_sfd": pair["target_ebv_sfd"],
                    "control_ebv_sfd": pair["control_ebv_sfd"],
                    "delta_ebv_sfd_control_minus_target": pair["delta_ebv_sfd_control_minus_target"],
                    "extinction_coefficient_mag_per_ebv": pair["extinction_coefficient_mag_per_ebv"],
                    "n_common_pixels": value["n_common_pixels"],
                    "raw_njy_arcsec2": value["raw"],
                    "extinction_corrected_njy_arcsec2": value["corrected"],
                    "correction_njy_arcsec2": (None if value["raw"] is None
                                                else value["corrected"] - value["raw"]),
                    "science_measurement": False,
                })
            for low, high in RELATIVE_WINDOWS:
                value = pair["windows"][(low, high, "relative")]
                rows.append({
                    "target_id": pair["target_id"], "control_id": pair["control_id"],
                    "filter": band, "metric": "inner_minus_300_800", "r_min_kpc": low,
                    "r_max_kpc": high, "ra_deg": pair["ra_deg"], "dec_deg": pair["dec_deg"],
                    "z": pair["z"], "target_ebv_sfd": pair["target_ebv_sfd"],
                    "control_ebv_sfd": pair["control_ebv_sfd"],
                    "delta_ebv_sfd_control_minus_target": pair["delta_ebv_sfd_control_minus_target"],
                    "extinction_coefficient_mag_per_ebv": pair["extinction_coefficient_mag_per_ebv"],
                    "n_common_pixels": value["n_common_pixels"],
                    "raw_njy_arcsec2": value["raw"],
                    "extinction_corrected_njy_arcsec2": value["corrected"],
                    "correction_njy_arcsec2": (None if value["raw"] is None
                                                else value["corrected"] - value["raw"]),
                    "science_measurement": False,
                })
    return rows


def summarize_layout(rows: list[dict], control_by_target: dict, width: float,
                     phase: float, ra_origin: float, replicas: int,
                     seed: int) -> list[dict]:
    labels, blocks = connected_pair_blocks(rows, control_by_target, width, phase, ra_origin)
    groups = np.unique(labels)
    group_index = {group: index for index, group in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)), size=replicas)
    metrics = [("profile", lo, hi) for lo, hi in PROFILE_WINDOWS]
    metrics += [("inner_minus_300_800", lo, hi) for lo, hi in RELATIVE_WINDOWS]
    out = []
    for kind, low, high in metrics:
        key = (low, high) if kind == "profile" else (low, high, "relative")
        values = [pair["windows"][key] for pair in rows]
        raw_all = np.asarray([np.nan if value["raw"] is None else value["raw"]
                              for value in values], dtype=float)
        corr_all = np.asarray([np.nan if value["corrected"] is None else value["corrected"]
                               for value in values], dtype=float)
        valid = np.isfinite(raw_all) & np.isfinite(corr_all)
        raw, corr = raw_all[valid], corr_all[valid]
        valid_groups = row_group[valid]
        n_blocks = int(np.unique(valid_groups).size)
        if not raw.size:
            continue
        group_counts = np.bincount(valid_groups, minlength=len(groups))
        raw_sums = np.bincount(valid_groups, weights=raw, minlength=len(groups))
        corr_sums = np.bincount(valid_groups, weights=corr, minlength=len(groups))
        denominators = draws @ group_counts
        raw_boot = (draws @ raw_sums)[denominators > 0] / denominators[denominators > 0]
        corr_boot = (draws @ corr_sums)[denominators > 0] / denominators[denominators > 0]
        change_boot = corr_boot - raw_boot
        ci = lambda x: np.percentile(x, [2.5, 97.5]).tolist() if len(x) else [None, None]
        raw_ci, corr_ci, change_ci = ci(raw_boot), ci(corr_boot), ci(change_boot)
        out.append({
            "filter": rows[0]["filter"], "metric": kind,
            "r_min_kpc": low, "r_max_kpc": high,
            "block_width_deg": width, "grid_phase_fraction": phase,
            "n_supported_pairs": int(raw.size), "n_supported_blocks": n_blocks,
            "n_all_spatial_blocks": len(groups), "largest_block_pairs": max(b["n_pairs"] for b in blocks),
            "raw_mean_njy_arcsec2": float(raw.mean()),
            "raw_block_ci95_low": raw_ci[0], "raw_block_ci95_high": raw_ci[1],
            "corrected_mean_njy_arcsec2": float(corr.mean()),
            "corrected_block_ci95_low": corr_ci[0], "corrected_block_ci95_high": corr_ci[1],
            "correction_mean_njy_arcsec2": float((corr - raw).mean()),
            "correction_block_ci95_low": change_ci[0], "correction_block_ci95_high": change_ci[1],
            "bootstrap_replicas_requested": replicas,
            "bootstrap_replicas_valid": int(len(raw_boot)),
            "science_measurement": False,
        })
    return out


def render_figure(path: Path, summaries: list[dict]) -> None:
    centers = np.arange(len(PLOT_WINDOWS), dtype=float)
    labels = ["0–50", "50–100", "100–300", "300–800", "800–1,000"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharex=True, sharey=True)
    for ax, band in zip(axes, ("i", "i2")):
        rows = [r for r in summaries if r["filter"] == band and r["metric"] == "profile"
                and r["block_width_deg"] == 5.0 and r["grid_phase_fraction"] == 0.0
                and (r["r_min_kpc"], r["r_max_kpc"]) in PLOT_WINDOWS]
        order = {window: index for index, window in enumerate(PLOT_WINDOWS)}
        rows.sort(key=lambda r: order[(r["r_min_kpc"], r["r_max_kpc"])])
        raw = np.asarray([r["raw_mean_njy_arcsec2"] for r in rows])
        corr = np.asarray([r["corrected_mean_njy_arcsec2"] for r in rows])
        raw_low = np.asarray([r["raw_block_ci95_low"] for r in rows])
        raw_high = np.asarray([r["raw_block_ci95_high"] for r in rows])
        corr_low = np.asarray([r["corrected_block_ci95_low"] for r in rows])
        corr_high = np.asarray([r["corrected_block_ci95_high"] for r in rows])
        ax.errorbar(centers - 0.06, raw, yerr=np.vstack([raw - raw_low, raw_high - raw]),
                    color="#1f5a85", marker="o", linewidth=1.2, capsize=3,
                    label="Raw; 5° block 95% CI")
        ax.errorbar(centers + 0.06, corr, yerr=np.vstack([corr - corr_low, corr_high - corr]),
                    color="#bc5c22", marker="s", linewidth=1.2, capsize=3,
                    label="SFD-corrected; 5° block 95% CI")
        ax.axhline(0.0, color="0.35", linewidth=0.9)
        ax.set_xticks(centers, labels)
        ax.set_xlim(-0.45, len(PLOT_WINDOWS) - 0.55)
        ax.set_yscale("symlog", linthresh=5, linscale=0.8)
        ax.grid(axis="y", alpha=0.22)
        ax.text(0.04, 0.95, f"HSC-{band}", transform=ax.transAxes,
                ha="left", va="top", fontsize=10)
        ax.set_xlabel("Radial annulus [proper kpc]")
    axes[0].set_ylabel("Target − control surface brightness [nJy arcsec⁻²]")
    axes[1].legend(loc="upper right", fontsize=8.5, frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 connected-block bootstrap replicates are required")
    input_paths = tuple(EXPECTED_HASHES)
    for path, expected in EXPECTED_HASHES.items():
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"pinned input hash mismatch: {path}")
    target_doc = read_json(TARGET_MANIFEST)
    control_doc = read_json(CONTROL_MANIFEST)
    stack = read_json(STACK_SUMMARY)
    spatial = read_json(SPATIAL_SUMMARY)
    if (stack.get("status") != "complete_exploratory"
            or stack.get("science_measurement") is not False
            or stack.get("n_pairs_stacked") != 1022
            or stack.get("n_by_filter") != EXPECTED_PAIRS
            or spatial.get("status") != "complete_spatial_block_sensitivity"
            or spatial.get("science_measurement") is not False):
        raise ValueError("pinned exploratory S23B stack summary changed")
    spatial_runner = spatial.get("runner", {})
    spatial_runner_path = Path(spatial_runner.get("path", ""))
    if (not spatial_runner_path.is_file()
            or spatial_runner.get("sha256") != sha256(spatial_runner_path)):
        raise ValueError("pinned spatial-block helper source changed")
    pair_rows = read_tsv(PAIRED_PROFILES)
    pair_data = prepare_pairs(target_doc, control_doc, pair_rows)
    control_by_target = {str(row["target_cluster_id"]): row
                         for row in control_doc["clusters"]}
    summary_rows = []
    endpoints = []
    for rows in pair_data.values():
        for row in rows:
            control = control_by_target[row["target_id"]]
            endpoints.extend((float(row["ra_deg"]), float(control["ra_deg"])))
    ra_origin = largest_gap_origin(np.asarray(endpoints, dtype=float))
    if not math.isclose(ra_origin, float(spatial["ra_grid_origin_deg"]),
                        rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("largest-gap RA origin differs from the pinned spatial-block analysis")

    layout_index = 0
    for width in BLOCK_WIDTHS_DEG:
        for phase in BLOCK_PHASES:
            for band_index, band in enumerate(("i", "i2")):
                rows = pair_data[band]
                labels, blocks = connected_pair_blocks(
                    rows, control_by_target, width, phase, ra_origin)
                expected_layout = next(item for item in spatial["layout_summary"]
                                       if item["filter"] == band
                                       and item["block_width_deg"] == width
                                       and item["grid_phase_fraction"] == phase)
                if (len(blocks) != expected_layout["n_connected_blocks"]
                        or max(block["n_pairs"] for block in blocks)
                        != expected_layout["largest_block_pairs"]):
                    raise ValueError("connected-block labels do not reproduce the pinned audit")
                summary_rows.extend(summarize_layout(
                    rows, control_by_target, width, phase, ra_origin, replicas,
                    seed + layout_index * 101 + band_index))
                layout_index += 1

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "pair_extinction_profiles.tsv"
    summary_path = output / "block_bootstrap_summary.tsv"
    plot_path = output / "extinction_corrected_profiles.png"
    write_tsv(pair_path, pair_metric_rows(pair_data))
    write_tsv(summary_path, summary_rows)
    render_figure(plot_path, summary_rows)
    hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    result = {
        "schema": "astrafek-s23b-icl-extinction-sensitivity-v1",
        "status": "complete_exploratory",
        "n_pairs_by_filter": {band: len(rows) for band, rows in pair_data.items()},
        "extinction_coefficients_mag_per_ebv": COEFFICIENTS,
        "coefficient_source": COEFFICIENT_SOURCE,
        "block_widths_deg": list(BLOCK_WIDTHS_DEG),
        "block_phases_fraction": list(BLOCK_PHASES),
        "ra_grid_origin_deg": ra_origin,
        "bootstrap_replicates_per_filter_layout_metric": replicas,
        "bootstrap_seed": seed,
        "results_by_filter": {
            band: {
                "median_abs_delta_ebv": float(np.median(np.abs([
                    row["delta_ebv_sfd_control_minus_target"] for row in rows]))),
                "p95_abs_delta_ebv": float(np.percentile(np.abs([
                    row["delta_ebv_sfd_control_minus_target"] for row in rows]), 95)),
                "max_abs_delta_ebv": float(np.max(np.abs([
                    row["delta_ebv_sfd_control_minus_target"] for row in rows]))),
            } for band, rows in pair_data.items()
        },
        "input_sha256": hashes,
        "output_sha256": {
            pair_path.name: sha256(pair_path),
            summary_path.name: sha256(summary_path),
            plot_path.name: sha256(plot_path),
        },
        "script_sha256": sha256(Path(__file__).resolve()),
        "science_measurement": False,
        "science_image_pixels_read": False,
        "network_accessed": False,
        "limitations": [
            "This is an extinction sensitivity on saved common-support target/control profiles, not a full ICL measurement.",
            "The hscPipe i/i2 coefficients are nominal here; residual sky terms are not physical extragalactic flux and are also rescaled.",
            "No PSF-wing, sky-transfer, source-Poisson, coadd-covariance, mask, or member-completeness correction is applied.",
            "The 300-800 kpc reference annulus can contain cluster light and is not assumed to be blank sky.",
            "All outputs retain science_measurement=false.",
        ],
    }
    result_path = output / "summary.json"
    write_json(result_path, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261108)
    args = parser.parse_args()
    summary = run(args.output, args.replicas, args.seed)
    print(json.dumps({key: summary[key] for key in (
        "status", "n_pairs_by_filter", "extinction_coefficients_mag_per_ebv",
        "results_by_filter", "science_measurement", "science_image_pixels_read")},
        indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
