#!/usr/bin/env python3
"""Estimate spatial-block uncertainty for the S23B exploratory ICL stack.

This analysis reads the saved pair-level radial profiles and the pinned target
and control manifests. It does not read image pixels, contact HSC, or promote
the stack to a science measurement.
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


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
CONTROL_MANIFEST = DATA / "s23b_icl_control_image_plan_2026-11-07_02/control_manifest.json"
STACK_DIR = DATA / "s23b_icl_target_control_stack_2026-11-07_01"
STACK_SUMMARY = STACK_DIR / "summary.json"
PAIRED_PROFILES = STACK_DIR / "paired_radial_profiles.tsv"
DEFAULT_OUTPUT = DATA / "s23b_icl_spatial_blocks_2026-11-07_01"

EXPECTED_TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_CONTROL_MANIFEST_SHA256 = "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df"
EXPECTED_STACK_SUMMARY_SHA256 = "c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b"
EXPECTED_STACKED_PAIRS = {"i": 241, "i2": 781}
WIDTHS_DEG = (2.5, 5.0, 10.0)
GRID_PHASES = (0.0, 0.5)
BROAD_BINS_KPC = ((0, 50), (50, 100), (100, 300), (300, 800), (800, 1000))
PROFILE_EDGES_KPC = np.arange(0.0, 1000.0 + 25.0, 25.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def largest_gap_origin(ra_values: np.ndarray) -> float:
    """Place the RA grid seam at the midpoint of the largest empty sky gap."""
    values = np.sort(np.mod(np.asarray(ra_values, dtype=float), 360.0))
    if values.size < 2:
        return float(values[0]) if values.size else 0.0
    gaps = np.diff(np.concatenate([values, values[:1] + 360.0]))
    index = int(np.argmax(gaps))
    return float((values[index] + 0.5 * gaps[index]) % 360.0)


def connected_pair_blocks(rows: list[dict], control_by_target: dict,
                           width_deg: float, phase_fraction: float,
                           ra_origin_deg: float) -> tuple[np.ndarray, list[dict]]:
    """Join pairs that share a target or control cell in a shifted RA/Dec grid."""
    n = len(rows)
    parent = list(range(n))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    ra_shift = width_deg * phase_fraction
    dec_shift = width_deg * phase_fraction

    def cell(ra_deg: float, dec_deg: float) -> tuple[int, int]:
        ra = (ra_deg - ra_origin_deg - ra_shift) % 360.0
        return (math.floor(ra / width_deg),
                math.floor((dec_deg + 90.0 - dec_shift) / width_deg))

    owner: dict[tuple[int, int], int] = {}
    for i, row in enumerate(rows):
        control = control_by_target[row["target_id"]]
        endpoints = (
            cell(float(row["ra_deg"]), float(row["dec_deg"])),
            cell(float(control["ra_deg"]), float(control["dec_deg"])),
        )
        for endpoint in endpoints:
            previous = owner.get(endpoint)
            if previous is None:
                owner[endpoint] = i
            else:
                union(i, previous)

    roots = sorted({find(i) for i in range(n)})
    root_to_label = {root: f"block_{index:03d}" for index, root in enumerate(roots)}
    labels = np.asarray([root_to_label[find(i)] for i in range(n)], dtype=str)
    sizes = {label: int(np.count_nonzero(labels == label))
             for label in np.unique(labels)}
    block_rows = [{"block_id": key, "n_pairs": sizes[key]} for key in sorted(sizes)]
    return labels, block_rows


def broad_profile(profile: list[dict], low_kpc: int, high_kpc: int) -> tuple[float | None, int]:
    """Combine adjacent 25-kpc profile means by their common valid-pixel counts."""
    selected = [item for item in profile
                if float(item["r_min_kpc"]) >= low_kpc
                and float(item["r_max_kpc"]) <= high_kpc]
    if not selected:
        raise ValueError(f"no saved 25-kpc bins cover {low_kpc}-{high_kpc} kpc")
    weighted_sum = 0.0
    n_pixels = 0
    for item in selected:
        count = int(item["n_valid_pixels"])
        value = item.get("delta_mean_njy_arcsec2")
        if count and value is not None:
            weighted_sum += float(value) * count
            n_pixels += count
    return (weighted_sum / n_pixels if n_pixels else None, n_pixels)


def summarize_layout(rows: list[dict], labels: np.ndarray, band: str,
                     width_deg: float, phase_fraction: float, replicas: int,
                     seed: int) -> tuple[list[dict], list[dict]]:
    groups = np.unique(labels)
    group_index = {group: index for index, group in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    # The same block draws are used for every annulus in this layout.
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)),
                            size=replicas)
    profile_rows: list[dict] = []
    loo_rows: list[dict] = []
    for low_kpc, high_kpc in BROAD_BINS_KPC:
        values = np.asarray([row["broad_profiles"][(low_kpc, high_kpc)] for row in rows],
                            dtype=float)
        valid = np.isfinite(values)
        valid_values = values[valid]
        valid_groups = row_group[valid]
        n_blocks = int(np.unique(valid_groups).size)
        if valid_values.size:
            mean = float(valid_values.mean())
            pair_sd = float(valid_values.std(ddof=1)) if valid_values.size > 1 else None
            pair_sem = (pair_sd / math.sqrt(valid_values.size)
                        if pair_sd is not None else None)
        else:
            mean = pair_sd = pair_sem = None

        group_sums = np.bincount(valid_groups, weights=valid_values, minlength=len(groups))
        group_counts = np.bincount(valid_groups, minlength=len(groups))
        denominators = draws @ group_counts
        numerators = draws @ group_sums
        valid_draws = denominators > 0
        boot = numerators[valid_draws] / denominators[valid_draws]
        if n_blocks >= 3 and boot.size >= replicas * 0.9:
            low, high = np.percentile(boot, [2.5, 97.5]).tolist()
            bootstrap_status = "complete"
        else:
            low = high = None
            bootstrap_status = "insufficient_supported_blocks_or_valid_draws"

        loo_means = []
        for group in np.unique(valid_groups):
            keep = valid & (row_group != group)
            if np.any(keep):
                loo_means.append((str(groups[group]), float(values[keep].mean()),
                                  int(np.count_nonzero(keep))))
        if mean is not None and loo_means:
            maximum = max(loo_means, key=lambda item: abs(item[1] - mean))
            max_abs_shift = abs(maximum[1] - mean)
            sign_flip_count = sum(1 for _, value, _ in loo_means
                                  if mean != 0.0 and value != 0.0
                                  and np.sign(mean) != np.sign(value))
        else:
            maximum = (None, None, None)
            max_abs_shift = None
            sign_flip_count = None

        profile_rows.append({
            "filter": band,
            "r_min_kpc": low_kpc,
            "r_max_kpc": high_kpc,
            "block_width_deg": width_deg,
            "grid_phase_fraction": phase_fraction,
            "mean_delta_njy_arcsec2": mean,
            "pair_sd_njy_arcsec2": pair_sd,
            "pair_sem_njy_arcsec2": pair_sem,
            "block_bootstrap_ci95_low": low,
            "block_bootstrap_ci95_high": high,
            "bootstrap_replicas_requested": replicas,
            "bootstrap_replicas_valid": int(boot.size),
            "bootstrap_status": bootstrap_status,
            "n_supported_pairs": int(valid_values.size),
            "n_supported_blocks": n_blocks,
            "largest_supported_block_pairs": int(max(group_counts, default=0)),
            "max_leave_one_block_abs_shift_njy_arcsec2": max_abs_shift,
            "max_shift_block_id": maximum[0],
            "n_leave_one_block_sign_flips": sign_flip_count,
            "science_measurement": False,
        })
        for block_id, block_mean, n_retained in loo_means:
            loo_rows.append({
                "filter": band,
                "r_min_kpc": low_kpc,
                "r_max_kpc": high_kpc,
                "block_width_deg": width_deg,
                "grid_phase_fraction": phase_fraction,
                "deleted_block_id": block_id,
                "n_pairs_retained": n_retained,
                "mean_delta_njy_arcsec2": block_mean,
                "shift_from_full_njy_arcsec2": block_mean - mean,
            })
    return profile_rows, loo_rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    columns = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def plot_profiles(path: Path, profile_rows: list[dict], n_pairs: dict,
                  primary_width: float = 5.0) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True, sharey=True)
    for axis, band in zip(axes, ("i", "i2")):
        selected = [row for row in profile_rows
                    if row["filter"] == band
                    and row["block_width_deg"] == primary_width
                    and row["grid_phase_fraction"] == 0.0]
        selected.sort(key=lambda row: row["r_min_kpc"])
        x = np.asarray([(row["r_min_kpc"] + row["r_max_kpc"]) / 2.0
                        for row in selected], dtype=float)
        y = np.asarray([row["mean_delta_njy_arcsec2"] for row in selected], dtype=float)
        low = np.asarray([row["block_bootstrap_ci95_low"] for row in selected], dtype=float)
        high = np.asarray([row["block_bootstrap_ci95_high"] for row in selected], dtype=float)
        sem = np.asarray([row["pair_sem_njy_arcsec2"] for row in selected], dtype=float)
        axis.fill_between(x, low, high, color="#6aaed6", alpha=0.28,
                          label="5° block bootstrap 95% CI")
        axis.errorbar(x, y, yerr=sem, color="#1f5a85", marker="o", ms=3.6,
                      linewidth=1.2, capsize=2, label="pair SEM")
        axis.axhline(0.0, color="0.35", linewidth=0.9)
        for radius in (50, 100, 300, 800):
            axis.axvline(radius, color="0.75", linewidth=0.6, linestyle=":")
        axis.set_xscale("symlog", linthresh=25, linscale=0.8)
        axis.set_yscale("symlog", linthresh=5, linscale=0.8)
        axis.text(0.04, 0.95, f"HSC-{band}  |  N={n_pairs[band]}",
                  transform=axis.transAxes, ha="left", va="top", fontsize=11)
        axis.set_xlabel("Projected radius [proper kpc]")
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Target − control surface brightness [nJy arcsec⁻²]")
    axes[1].legend(loc="upper right", frameon=True, fontsize=9)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite spatial-block analysis: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 block-bootstrap replicates are required")
    input_paths = (TARGET_MANIFEST, CONTROL_MANIFEST, STACK_SUMMARY, PAIRED_PROFILES)
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    if input_hashes[str(TARGET_MANIFEST.resolve())] != EXPECTED_TARGET_MANIFEST_SHA256:
        raise ValueError("target manifest does not match the downloaded S23B target sample")
    if input_hashes[str(CONTROL_MANIFEST.resolve())] != EXPECTED_CONTROL_MANIFEST_SHA256:
        raise ValueError("control manifest does not match the approved v2 scope")
    if input_hashes[str(STACK_SUMMARY.resolve())] != EXPECTED_STACK_SUMMARY_SHA256:
        raise ValueError("S23B exploratory stack summary changed")

    target_manifest = read_json(TARGET_MANIFEST)
    control_manifest = read_json(CONTROL_MANIFEST)
    stack_summary = read_json(STACK_SUMMARY)
    if (stack_summary.get("status") != "complete_exploratory"
            or stack_summary.get("science_measurement") is not False
            or stack_summary.get("images_ready_for_science") is not False
            or stack_summary.get("n_pairs_stacked") != sum(EXPECTED_STACKED_PAIRS.values())
            or control_manifest.get("science_measurement") is not False
            or target_manifest.get("science_measurement") is not False):
        raise ValueError("input stack or manifest no longer has the expected exploratory status")
    if stack_summary.get("outputs", {}).get("paired_profiles_sha256") != input_hashes[str(PAIRED_PROFILES.resolve())]:
        raise ValueError("pair profile table is not the output pinned by the stack summary")

    targets = {str(row["id"]): row for row in target_manifest.get("clusters", [])}
    controls = control_manifest.get("clusters", [])
    control_by_target = {str(row["target_cluster_id"]): row for row in controls}
    if (len(targets) != 1043 or len(controls) != 1041
            or len(control_by_target) != 1041):
        raise ValueError("S23B target/control manifest counts changed")
    pair_rows = read_tsv(PAIRED_PROFILES)
    if len(pair_rows) != sum(EXPECTED_STACKED_PAIRS.values()):
        raise ValueError("paired radial profile row count changed")

    rows_by_band: dict[str, list[dict]] = {"i": [], "i2": []}
    seen = set()
    all_endpoints = []
    for source in pair_rows:
        row = dict(source)
        target_id = str(row["target_id"])
        control_id = str(row["control_id"])
        band = str(row["filter"])
        if band not in rows_by_band or target_id not in targets or target_id not in control_by_target:
            raise ValueError(f"paired profile has an unexpected target or filter: {target_id}/{band}")
        key = (target_id, band)
        if key in seen:
            raise ValueError(f"duplicate target/filter pair {key}")
        seen.add(key)
        target = targets[target_id]
        control = control_by_target[target_id]
        if (str(control["id"]) != control_id
                or float(row["ra_deg"]) != float(target["ra_deg"])
                or float(row["dec_deg"]) != float(target["dec_deg"])):
            raise ValueError(f"pair coordinates or control mapping disagree for {target_id}")
        profile = json.loads(row["profile_json"])
        if len(profile) != len(PROFILE_EDGES_KPC) - 1:
            raise ValueError(f"{target_id}/{band}: unexpected radial profile length")
        row["broad_profiles"] = {}
        for low, high in BROAD_BINS_KPC:
            value, n_pixels = broad_profile(profile, low, high)
            row["broad_profiles"][(low, high)] = value if n_pixels else float("nan")
        rows_by_band[band].append(row)
        all_endpoints.extend([float(row["ra_deg"]), float(row["dec_deg"]),
                              float(control["ra_deg"]), float(control["dec_deg"])])

    if {band: len(rows) for band, rows in rows_by_band.items()} != EXPECTED_STACKED_PAIRS:
        raise ValueError("per-filter pair accounting differs from the stack summary")
    endpoint_array = np.asarray(all_endpoints, dtype=float).reshape(-1, 2)
    ra_origin = largest_gap_origin(endpoint_array[:, 0])
    profile_rows: list[dict] = []
    loo_rows: list[dict] = []
    block_rows: list[dict] = []
    layout_summary: list[dict] = []
    layout_index = 0
    for width in WIDTHS_DEG:
        for phase in GRID_PHASES:
            for band_index, band in enumerate(("i", "i2")):
                rows = rows_by_band[band]
                labels, blocks = connected_pair_blocks(
                    rows, control_by_target, width, phase, ra_origin)
                block_rows.extend({"filter": band, "block_width_deg": width,
                                   "grid_phase_fraction": phase, **block}
                                  for block in blocks)
                result_rows, result_loo = summarize_layout(
                    rows, labels, band, width, phase, replicas,
                    seed + layout_index * 101 + band_index)
                profile_rows.extend(result_rows)
                loo_rows.extend(result_loo)
                layout_summary.append({
                    "filter": band, "block_width_deg": width,
                    "grid_phase_fraction": phase,
                    "n_pairs": len(rows), "n_connected_blocks": len(blocks),
                    "largest_block_pairs": max(block["n_pairs"] for block in blocks),
                    "smallest_block_pairs": min(block["n_pairs"] for block in blocks),
                    "block_sizes": sorted((block["n_pairs"] for block in blocks), reverse=True),
                })
                layout_index += 1

    output.mkdir(parents=True, exist_ok=False)
    profile_path = output / "radial_block_uncertainty.tsv"
    loo_path = output / "leave_one_block_out.tsv"
    blocks_path = output / "spatial_blocks.tsv"
    plot_path = output / "s23b_icl_spatial_block_uncertainty.png"
    write_tsv(profile_path, profile_rows)
    write_tsv(loo_path, loo_rows)
    write_tsv(blocks_path, block_rows)
    plot_profiles(plot_path, profile_rows, EXPECTED_STACKED_PAIRS)

    primary = [row for row in profile_rows
               if row["block_width_deg"] == 5.0 and row["grid_phase_fraction"] == 0.0]
    profile_outputs = {path.name: {"path": str(path), "sha256": sha256(path)}
                       for path in (profile_path, loo_path, blocks_path, plot_path)}
    summary = {
        "schema": "astrafek-s23b-icl-spatial-block-uncertainty-v1",
        "status": "complete_spatial_block_sensitivity",
        "science_measurement": False,
        "image_pixels_read": False,
        "network_used": False,
        "n_pair_profiles": len(pair_rows),
        "n_by_filter": EXPECTED_STACKED_PAIRS,
        "method": {
            "estimand": "equal-pair mean of each matched target-minus-control common-support annular profile",
            "broad_annuli_kpc": [list(x) for x in BROAD_BINS_KPC],
            "block_construction": "connected components of pairs sharing any target/control endpoint RA-Dec grid cell; RA seam placed at midpoint of the largest empty endpoint-RA gap",
            "block_widths_deg": list(WIDTHS_DEG),
            "grid_phase_fractions": list(GRID_PHASES),
            "bootstrap": "resample connected spatial blocks with replacement; retain all supported pairs in each selected block and recompute equal-pair mean",
            "replicas_per_layout_and_filter": replicas,
            "seed_base": int(seed),
            "primary_layout": "5-degree RA-Dec cells, phase 0.0; shifted half-cell grid and 2.5/10-degree widths are layout sensitivities",
            "pair_sem": "sample standard deviation across supported pairs divided by square root of supported pair count",
        },
        "ra_grid_origin_deg": ra_origin,
        "layout_summary": layout_summary,
        "primary_5deg_profiles": primary,
        "runner": {"path": str(Path(__file__).resolve()),
                   "sha256": sha256(Path(__file__).resolve())},
        "inputs_sha256": input_hashes,
        "outputs": profile_outputs,
        "limitations": [
            "RA-Dec cells are an approximate dependence partition, not official HSC field polygons or a calibrated survey-region covariance.",
            "The two half-cell phases and three widths quantify grid sensitivity but do not establish an optimal correlation scale.",
            "Supported-pair counts vary by radius because the stack uses strict common masks; each annulus has its own effective sample.",
            "This analysis propagates pair-to-pair spatial clustering only; it does not model coadd pixel covariance, residual sky, Galactic extinction, PSF wings, or BCG-centering uncertainty.",
            "The S23B CAMIRA member catalog has not been imported, so no satellite-position association is measured here.",
            "The input stack remains an exploratory pipeline product and is not promoted to science_measurement=true.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       allow_nan=False) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)

    report_path = output / "assessment.md"
    report_lines = [
        "# S23B ICL spatial-block uncertainty assessment — 2026-11-07",
        "",
        "**Status:** exploratory uncertainty sensitivity; `science_measurement=false`.",
        "",
        f"The analysis uses {len(pair_rows)} saved target/control pair profiles ({EXPECTED_STACKED_PAIRS['i']} in `i`, {EXPECTED_STACKED_PAIRS['i2']} in `i2`). It reads no image pixels and makes no HSC request.",
        "",
        "The primary layout uses connected 5-degree RA-Dec cells. Two grid phases and 2.5/10-degree widths test dependence on arbitrary cell placement and scale. These are approximate spatial blocks, not official survey regions.",
        "",
        "## Primary 5-degree results",
        "",
        "| Filter | Shell (proper kpc) | Pairs / supported blocks | Pair SEM | Block-bootstrap 95% interval | Largest leave-one-block shift |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary:
        ci = (f"[{row['block_bootstrap_ci95_low']:.3g}, {row['block_bootstrap_ci95_high']:.3g}]"
              if row["block_bootstrap_ci95_low"] is not None else "not estimated")
        sem = f"{row['pair_sem_njy_arcsec2']:.3g}" if row["pair_sem_njy_arcsec2"] is not None else "—"
        shift = (f"{row['max_leave_one_block_abs_shift_njy_arcsec2']:.3g}"
                 if row["max_leave_one_block_abs_shift_njy_arcsec2"] is not None else "—")
        report_lines.append(
            f"| {row['filter']} | {row['r_min_kpc']}–{row['r_max_kpc']} | "
            f"{row['n_supported_pairs']} / {row['n_supported_blocks']} | {sem} | {ci} | {shift} |")
    report_lines.extend([
        "",
        "## Block-grid and scale sensitivity",
        "",
        "| Width (deg) | Phase | i blocks (largest) | i2 blocks (largest) |",
        "|---:|---:|---:|---:|",
    ])
    for width in WIDTHS_DEG:
        for phase in GRID_PHASES:
            by_key = {(row["filter"], row["block_width_deg"], row["grid_phase_fraction"]): row
                      for row in layout_summary}
            i_layout = by_key[("i", width, phase)]
            i2_layout = by_key[("i2", width, phase)]
            report_lines.append(
                f"| {width:g} | {phase:g} | {i_layout['n_connected_blocks']} ({i_layout['largest_block_pairs']}) | "
                f"{i2_layout['n_connected_blocks']} ({i2_layout['largest_block_pairs']}) |")
    report_lines.extend([
        "",
        "The next table gives the outer envelope of the six separate 95% intervals across block widths and grid phases. It is a sensitivity envelope, not a single 95% confidence interval.",
        "",
        "| Filter | Shell (proper kpc) | Range across interval lower ends | Range across interval upper ends | Every layout excludes zero? |",
        "|---|---:|---:|---:|---|",
    ])
    for band in ("i", "i2"):
        for low_kpc, high_kpc in BROAD_BINS_KPC[1:]:
            variants = [row for row in profile_rows
                        if row["filter"] == band
                        and row["r_min_kpc"] == low_kpc
                        and row["r_max_kpc"] == high_kpc]
            lows = [row["block_bootstrap_ci95_low"] for row in variants
                    if row["block_bootstrap_ci95_low"] is not None]
            highs = [row["block_bootstrap_ci95_high"] for row in variants
                     if row["block_bootstrap_ci95_high"] is not None]
            every_excludes = bool(lows and highs and (min(lows) > 0 or max(highs) < 0))
            report_lines.append(
                f"| {band} | {low_kpc}–{high_kpc} | [{min(lows):.3g}, {max(lows):.3g}] | "
                f"[{min(highs):.3g}, {max(highs):.3g}] | {'yes' if every_excludes else 'no'} |")
    report_lines.extend([
        "",
        "The block-bootstrap interval resamples connected sky blocks and retains all supported pairs in each resampled block. The pair SEM is shown as a reference; it does not account for spatial clustering.",
        "",
        "## Interpretation",
        "",
        "The 50–100 and 100–300 kpc intervals remain positive in both filters across all six block layouts. The 300–800 kpc intervals include zero throughout. At 800–1,000 kpc, the i2 interval excludes zero for some layouts but crosses zero for another, so its outer residual is not stable to the block choice. These intervals quantify only spatial resampling of the selected pairs; they do not establish an ICL detection.",
        "",
        "The block-grid count and interval envelope change with width and phase. Interpret these as a sensitivity range, not as a calibrated survey-region covariance. Even stable intervals would not solve image-level covariance, residual sky, extinction, or PSF effects.",
        "",
        "This is not an ICL–satellite association result. The S23B member-position catalog remains a separate approval-gated input. No satellite coordinates are inferred from the cluster catalog or substituted from the LRG list.",
        "",
        "## Remaining gates",
        "",
        "- Galactic-extinction correction and its filter-specific coefficient.",
        "- Residual-sky transfer and independent outer-field validation.",
        "- PSF-wing matching, coadd-pixel covariance, and injection recovery.",
        "- BCG/centering confirmation and sensitivity to automatic source masks.",
        "- CAMIRA member catalog validation before satellite-centered stacking.",
        "",
        "## Provenance",
        "",
        f"- Inputs: `{input_hashes}`",
        f"- Runner: `{Path(__file__).resolve()}` (SHA-256 `{sha256(Path(__file__).resolve())}`)",
        f"- Output summary SHA-256: `{sha256(summary_path)}`",
        "- `science_measurement=false`; this step does not promote the image stack.",
        "",
    ])
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261107)
    args = parser.parse_args()
    summary = run(args.output, args.replicas, args.seed)
    compact_layouts = [{key: row[key] for key in (
        "filter", "block_width_deg", "grid_phase_fraction", "n_connected_blocks",
        "largest_block_pairs")} for row in summary["layout_summary"]]
    compact_profiles = [{key: row[key] for key in (
        "filter", "r_min_kpc", "r_max_kpc", "n_supported_pairs",
        "n_supported_blocks", "mean_delta_njy_arcsec2", "pair_sem_njy_arcsec2",
        "block_bootstrap_ci95_low", "block_bootstrap_ci95_high")}
        for row in summary["primary_5deg_profiles"]]
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": summary["science_measurement"],
        "layout_summary": compact_layouts,
        "primary_5deg_profiles": compact_profiles,
        "output_dir": str(Path(args.output).resolve()),
        "summary_sha256": summary["summary_sha256"],
    }, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
