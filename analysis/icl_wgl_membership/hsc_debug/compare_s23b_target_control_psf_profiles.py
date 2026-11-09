#!/usr/bin/env python3
"""Compare diagnostic Gaia-candidate profiles in matched S23B target/control frames.

Profiles are first averaged within each frame, then target and control frames
are paired by the previously prepared seeing/filter match. This is a local
diagnostic only: it does not calibrate an effective coadd PSF or an ICL wing
correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PILOT = DATA / "s23b_psf_control_candidates_2026-11-07_01"
EXPANSION = DATA / "s23b_psf_expansion_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_psf_target_control_profile_comparison_2026-11-07_03"
BOOTSTRAP_DRAWS = 5000
SEED = 20261108


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def profile_key(row: dict[str, str]) -> tuple:
    return (
        float(row["kron_factor"]), row["mask_name"],
        int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"]),
        float(row["r_min_pix"]), float(row["r_max_pix"]),
    )


def frame_profiles(rows: list[dict[str, str]], *, band: str, side: str,
                   allowed_frames: set[str]) -> tuple[dict, dict, dict]:
    values: dict[tuple, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    row_counts: dict[tuple, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    statuses: dict[str, int] = defaultdict(int)
    for row in rows:
        if row["frame_id"] not in allowed_frames:
            raise ValueError(f"profile row frame is outside paired selection: {row['frame_id']}")
        if row["science_measurement"] != "False" or row["psf_approved"] != "False":
            raise ValueError("profile input carries an unexpected science or PSF approval state")
        statuses[row["status"]] += 1
        key = profile_key(row)
        values[key]
        row_counts[key][row["frame_id"]] += 1
        value = row["profile_per_core_per_pixel"]
        if row["status"] == "ok" and value != "":
            values[key][row["frame_id"]].append(float(value))
    return values, row_counts, dict(statuses)


def paired_rows(target_rows: list[dict[str, str]], control_rows: list[dict[str, str]],
                *, band: str, subset: str, pairs: list[dict[str, str]],
                rng: np.random.Generator) -> tuple[list[dict], dict]:
    target_frames = {row["target_cluster_id"] for row in pairs if row["target_primary_filter"] == band}
    control_to_target = {
        row["control_id"]: row["target_cluster_id"] for row in pairs
        if row["target_primary_filter"] == band
    }
    control_frames = set(control_to_target)
    target_values, _, target_statuses = frame_profiles(
        target_rows, band=band, side="target", allowed_frames=target_frames
    )
    control_values, _, control_statuses = frame_profiles(
        control_rows, band=band, side="control", allowed_frames=control_frames
    )
    target_keys = set(target_values)
    control_keys = set(control_values)
    if target_keys != control_keys:
        raise ValueError(f"target/control profile grids differ for {band}/{subset}")

    output_rows = []
    n_paired_by_key = []
    radii_by_key: dict[tuple, list[float]] = defaultdict(list)
    for row in [*target_rows, *control_rows]:
        radii_by_key[profile_key(row)].append(float(row["r_mid_arcsec"]))
    for key in sorted(target_keys):
        factor, mask_name, sky_lo, sky_hi, rlo, rhi = key
        radius = float(np.mean(radii_by_key[key]))
        target_frame_mean = {
            frame: float(np.mean(values)) for frame, values in target_values[key].items() if values
        }
        control_frame_mean = {
            control_to_target[frame]: float(np.mean(values))
            for frame, values in control_values[key].items() if values
        }
        paired = sorted(set(target_frame_mean) & set(control_frame_mean))
        differences = np.asarray([
            control_frame_mean[pair_id] - target_frame_mean[pair_id] for pair_id in paired
        ], dtype=float)
        if differences.size:
            estimate = float(np.mean(differences))
            median = float(np.median(differences))
        else:
            estimate = median = None
        if differences.size >= 2:
            draw_idx = rng.integers(0, differences.size, size=(BOOTSTRAP_DRAWS, differences.size))
            boot_means = np.mean(differences[draw_idx], axis=1)
            boot_q16, boot_q84, boot_q025, boot_q975 = (
                float(x) for x in np.quantile(boot_means, [0.16, 0.84, 0.025, 0.975])
            )
        else:
            boot_q16 = boot_q84 = boot_q025 = boot_q975 = None
        target_mean = float(np.mean([target_frame_mean[p] for p in paired])) if paired else None
        control_mean = float(np.mean([control_frame_mean[p] for p in paired])) if paired else None
        output_rows.append({
            "filter": band,
            "subset": subset,
            "kron_factor": factor,
            "mask_name": mask_name,
            "sky_r_min_pix": sky_lo,
            "sky_r_max_pix": sky_hi,
            "r_min_pix": rlo,
            "r_max_pix": rhi,
            "r_mid_arcsec": radius,
            "n_target_frames_with_valid_profile": len(target_frame_mean),
            "n_control_frames_with_valid_profile": len(control_frame_mean),
            "n_matched_frame_pairs": int(differences.size),
            "target_frame_equal_mean_on_paired_support": target_mean,
            "control_frame_equal_mean_on_paired_support": control_mean,
            "control_minus_target_mean": estimate,
            "control_minus_target_median": median,
            "bootstrap_mean_q16": boot_q16,
            "bootstrap_mean_q84": boot_q84,
            "bootstrap_mean_q025": boot_q025,
            "bootstrap_mean_q975": boot_q975,
            "science_measurement": False,
            "psf_approved": False,
        })
        n_paired_by_key.append(int(differences.size))

    summary = {
        "filter": band,
        "subset": subset,
        "n_target_frames_in_selection": len(target_frames),
        "n_control_frames_in_selection": len(control_frames),
        "target_profile_row_status_counts": target_statuses,
        "control_profile_row_status_counts": control_statuses,
        "n_profile_bins": len(output_rows),
        "paired_frame_count_by_profile_bin": {
            "min": min(n_paired_by_key) if n_paired_by_key else 0,
            "median": float(np.median(n_paired_by_key)) if n_paired_by_key else 0.0,
            "max": max(n_paired_by_key) if n_paired_by_key else 0,
            "bins_with_at_least_two_pairs": sum(n >= 2 for n in n_paired_by_key),
        },
    }
    return output_rows, summary


def plot_support(rows: list[dict], path: Path) -> None:
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    colors = {(120, 150): "#1f77b4", (160, 190): "#d95f02"}
    max_pairs = max((int(row["n_matched_frame_pairs"]) for row in rows), default=0)
    ymax = max(1.0, max_pairs + 0.4)
    for band, row_idx in (("i", 0), ("i2", 1)):
        for subset, col_idx in (("all", 0), ("outer500", 1)):
            ax = axes[row_idx, col_idx]
            for sky in ((120, 150), (160, 190)):
                selected = [
                    row for row in rows
                    if row["filter"] == band and row["subset"] == subset
                    and float(row["kron_factor"]) == 2.5
                    and float(row["r_min_pix"]) >= 24.0
                    and (int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"])) == sky
                ]
                selected.sort(key=lambda row: float(row["r_mid_arcsec"]))
                if not selected:
                    continue
                x = np.asarray([float(row["r_mid_arcsec"]) for row in selected])
                y = np.asarray([int(row["n_matched_frame_pairs"]) for row in selected])
                ax.step(x, y, where="mid", marker="o", ms=3, lw=1.1, color=colors[sky],
                        label=f"sky {sky[0]}–{sky[1]} px")
            ax.set_xlim(0, 16.2)
            ax.set_ylim(0, ymax)
            ax.set_yticks(range(max_pairs + 1))
            title = f"{band} · {'all radii' if subset == 'all' else 'center distance ≥500 kpc'}"
            ax.text(0.5, 0.96, title, transform=ax.transAxes, ha="center", va="top",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
            if row_idx == 0 and col_idx == 0 and ax.lines:
                ax.legend(loc="center right", framealpha=0.9, fontsize=9)
            ax.grid(alpha=0.18)
            if row_idx == 1:
                ax.set_xlabel("Radius [arcsec]")
            if col_idx == 0:
                ax.set_ylabel("Valid matched pairs [N]")
    for ax in axes.flat:
        ax.label_outer()
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.10, top=0.99, wspace=0, hspace=0)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def plot_flux_differences(rows: list[dict], path: Path) -> None:
    """Plot paired control-minus-target profiles for the primary Kron mask."""
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    sky_specs = (((120, 150), "#1f77b4"), ((160, 190), "#d95f02"))
    for band, row_idx in (("i", 0), ("i2", 1)):
        for subset, col_idx in (("all", 0), ("outer500", 1)):
            ax = axes[row_idx, col_idx]
            plotted = False
            for sky, color in sky_specs:
                selected = [
                    row for row in rows
                    if row["filter"] == band and row["subset"] == subset
                    and float(row["kron_factor"]) == 2.5
                    and (int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"])) == sky
                    and int(row["n_matched_frame_pairs"]) >= 2
                    and row["bootstrap_mean_q16"] not in (None, "")
                ]
                selected.sort(key=lambda row: float(row["r_mid_arcsec"]))
                if not selected:
                    continue
                x = np.asarray([float(row["r_mid_arcsec"]) for row in selected])
                y = 1e6 * np.asarray([float(row["control_minus_target_mean"]) for row in selected])
                low = 1e6 * np.asarray([
                    float(row["control_minus_target_mean"]) - float(row["bootstrap_mean_q16"])
                    for row in selected
                ])
                high = 1e6 * np.asarray([
                    float(row["bootstrap_mean_q84"]) - float(row["control_minus_target_mean"])
                    for row in selected
                ])
                ax.errorbar(x, y, yerr=np.vstack((low, high)), color=color,
                            marker="o", ms=3.5, lw=1.0, capsize=2,
                            label=f"sky {sky[0]}–{sky[1]} px")
                plotted = True
            ax.axhline(0, color="0.35", lw=0.8)
            ax.set_xlim(3.8, 16.2)
            ax.set_ylim(-50, 50)
            ax.text(0.97, 0.96,
                    f"{band} · {'all radii' if subset == 'all' else 'center distance ≥500 kpc'}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=10,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 2})
            if plotted and band == "i2" and subset == "all":
                ax.legend(loc="upper left", framealpha=0.9, fontsize=8)
            if band == "i" or not plotted:
                ax.text(0.5, 0.56, "No matched bin with N ≥ 2", transform=ax.transAxes,
                        ha="center", va="center", fontsize=9)
            ax.grid(alpha=0.18)
            if row_idx == 1:
                ax.set_xlabel("Radius [arcsec]")
            if col_idx == 0:
                ax.set_ylabel("Control − target [10⁻⁶ core flux pixel⁻¹]")
    for ax in axes.flat:
        ax.label_outer()
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.10, top=0.99,
                        wspace=0.0, hspace=0.0)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def analyze(output: Path, *, cohort: str = "pilot") -> dict:
    if cohort not in ("pilot", "expansion"):
        raise ValueError("cohort must be 'pilot' or 'expansion'")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite profile comparison: {output}")
    if cohort == "expansion":
        selection_path = EXPANSION / "plan/selection.tsv"
        pilot_summary_path = EXPANSION / "plan/summary.json"
    else:
        selection_path = PILOT / "selection.tsv"
        pilot_summary_path = PILOT / "summary.json"
    selection = read_tsv(selection_path)
    pilot_summary = read_json(pilot_summary_path)
    if cohort == "expansion":
        if (len(selection) != 200
                or len({row["target_cluster_id"] for row in selection}) != 200
                or len({row["control_id"] for row in selection}) != 200
                or sha256(selection_path) != pilot_summary.get("selection_table_sha256")
                or pilot_summary.get("n_pairs_selected") != 200):
            raise ValueError("matched expansion selection is incomplete or differs from its receipt")
        pairs = [{
            "target_cluster_id": row["target_cluster_id"],
            "control_id": row["control_id"],
            "target_primary_filter": row["filter"],
        } for row in selection]
    else:
        if (len(selection) != 49 or len({row["target_cluster_id"] for row in selection}) != 49
                or sha256(selection_path) != pilot_summary.get("selection_sha256")):
            raise ValueError("matched image selection is incomplete, duplicated, or differs from its receipt")
        pairs = selection

    output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    all_rows = []
    summaries = []
    input_paths = [selection_path, pilot_summary_path]
    for band in ("i", "i2"):
        for subset in ("all", "outer500"):
            if cohort == "expansion":
                target_dir = EXPANSION / f"target_wing_stack_{band}_{subset}_2026-11-07_01"
                control_dir = EXPANSION / f"control_wing_stack_{band}_{subset}_2026-11-07_01"
            else:
                target_dir = DATA / f"s23b_psf_wing_stack_{band}_{subset}_2026-11-07_02"
                control_dir = DATA / f"s23b_psf_control_wing_stack_{band}_{subset}_2026-11-07_01"
            target_profile_path = target_dir / "star_radial_profiles.tsv"
            control_profile_path = control_dir / "star_radial_profiles.tsv"
            target_summary_path = target_dir / "summary.json"
            control_summary_path = control_dir / "summary.json"
            for path in (target_profile_path, control_profile_path, target_summary_path, control_summary_path):
                if not path.is_file():
                    raise FileNotFoundError(path)
            target_summary = read_json(target_summary_path)
            control_summary = read_json(control_summary_path)
            if cohort == "expansion":
                target_stamp_summary_path = EXPANSION / "target_wing_stamps_2026-11-07_01/summary.json"
                control_stamp_summary_path = EXPANSION / "control_wing_stamps_2026-11-07_01/summary.json"
                target_stamp_summary = read_json(target_stamp_summary_path)
                control_stamp_summary = read_json(control_stamp_summary_path)
                target_ledger_path = (EXPANSION / "target_wing_stamps_2026-11-07_01/stack_ledgers"
                                      / f"{band}_{subset}/star_ledger.tsv")
                control_ledger_path = (EXPANSION / "control_wing_stamps_2026-11-07_01/stack_ledgers"
                                       / f"{band}_{subset}/star_ledger.tsv")
                target_ledger_hash = sha256(target_ledger_path)
                control_ledger_hash = sha256(control_ledger_path)
                if (target_stamp_summary.get("cohort") != "expansion"
                        or target_stamp_summary.get("population") != "target"
                        or control_stamp_summary.get("cohort") != "expansion"
                        or control_stamp_summary.get("population") != "control"
                        or target_stamp_summary.get("science_measurement") is not False
                        or control_stamp_summary.get("science_measurement") is not False
                        or target_stamp_summary.get("psf_approved") is not False
                        or control_stamp_summary.get("psf_approved") is not False
                        or target_stamp_summary.get("subset_ledgers", {}).get(f"{band}_{subset}", {}).get("ledger_sha256") != target_ledger_hash
                        or control_stamp_summary.get("subset_ledgers", {}).get(f"{band}_{subset}", {}).get("ledger_sha256") != control_ledger_hash
                        or target_summary.get("input_sha256", {}).get(str(target_ledger_path.resolve())) != target_ledger_hash
                        or control_summary.get("input_sha256", {}).get(str(control_ledger_path.resolve())) != control_ledger_hash):
                    raise ValueError(f"stack ledgers do not match their expansion stamp receipts for {band}/{subset}")
                input_paths.extend((target_stamp_summary_path, control_stamp_summary_path,
                                    target_ledger_path, control_ledger_path))
            if (sha256(target_profile_path) != target_summary.get("profile_tsv_sha256")
                    or sha256(control_profile_path) != control_summary.get("profile_tsv_sha256")
                    or target_summary.get("science_measurement") is not False
                    or control_summary.get("science_measurement") is not False
                    or target_summary.get("psf_approved") is not False
                    or control_summary.get("psf_approved") is not False):
                raise ValueError(f"target/control profile receipt mismatch for {band}/{subset}")
            target_profiles = read_tsv(target_profile_path)
            control_profiles = read_tsv(control_profile_path)
            radial_target = {profile_key(row) for row in target_profiles}
            radial_control = {profile_key(row) for row in control_profiles}
            if radial_target != radial_control:
                raise ValueError(f"target/control radial grid differs for {band}/{subset}")
            band_rows, band_summary = paired_rows(
                target_profiles, control_profiles, band=band, subset=subset,
                pairs=pairs, rng=rng,
            )
            all_rows.extend(band_rows)
            summaries.append(band_summary)
            input_paths.extend((target_profile_path, control_profile_path,
                                target_summary_path, control_summary_path))

    profile_path = output / "paired_profile_differences.tsv"
    write_tsv(profile_path, all_rows, list(all_rows[0]))
    figure_path = output / "paired_profile_support.png"
    plot_support(all_rows, figure_path)
    flux_figure_path = output / "paired_profile_flux_difference.png"
    plot_flux_differences(all_rows, flux_figure_path)
    support_text = ", ".join(
        f"{row['filter']} {row['subset']} {row['paired_frame_count_by_profile_bin']['bins_with_at_least_two_pairs']}/{row['n_profile_bins']}"
        for row in summaries
    )
    has_estimable_bins = any(
        row["paired_frame_count_by_profile_bin"]["bins_with_at_least_two_pairs"] > 0
        for row in summaries
    )
    if has_estimable_bins:
        uncertainty_scope = (
            "Frame-pair bootstrap intervals are computed only in bins with at least two matched pairs. "
            f"Supported bins by filter and subset are {support_text}. "
            "The intervals omit sky-annulus transfer, coadd covariance, and selection uncertainty."
        )
    else:
        uncertainty_scope = (
            "No radial bin has at least two matched pairs, so no paired bootstrap interval is estimable."
        )
    summary = {
        "schema": "astrafek-s23b-target-control-psf-profile-comparison-v2",
        "status": "complete_paired_diagnostic_only",
        "cohort": cohort,
        "pairing_unit": "matched target/control HSC image pair; stars averaged within each frame",
        "control_minus_target_definition": "frame-mean normalized profile in control minus that in its matched target frame",
        "bootstrap_unit": "matched frame pair",
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "random_seed": SEED,
        "filters_and_subsets": summaries,
        "plot_profile": "number of matched frame pairs with valid profiles per radial bin; Kron 2.5; two local sky annuli",
        "paired_flux_difference_estimable": any(
            band_summary["paired_frame_count_by_profile_bin"]["bins_with_at_least_two_pairs"] > 0
            for band_summary in summaries
        ),
        "uncertainty_scope": uncertainty_scope,
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs_sha256": {str(path.resolve()): sha256(path) for path in input_paths},
        "output_sha256": {
            profile_path.name: sha256(profile_path),
            figure_path.name: sha256(figure_path),
            flux_figure_path.name: sha256(flux_figure_path),
        },
        "science_measurement": False,
        "psf_approved": False,
        "external_queries_submitted": False,
        "limitations": [
            "Gaia-to-HSC epoch propagation remains unresolved; candidate coordinates were not motion-corrected.",
            "The radial profile is background-subtracted by local annuli that may contain PSF wings, cluster light, or residual sky.",
            "Control center distance is measured from its random sky position; outer500 is not a control-cluster-centric radius.",
            "Only frame pairs with valid profiles on both sides enter each difference bin; support changes with radius and sky annulus.",
            "This diagnostic does not establish an effective coadd PSF or an ICL correction.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=("pilot", "expansion"), default="pilot")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.cohort == "expansion" and args.output == DEFAULT_OUTPUT:
        output = EXPANSION / "target_control_profile_comparison_2026-11-07_01"
    result = analyze(output, cohort=args.cohort)
    print(json.dumps({
        "status": result["status"],
        "bootstrap_draws": result["bootstrap_draws"],
        "filters_and_subsets": result["filters_and_subsets"],
        "science_measurement": result["science_measurement"],
        "psf_approved": result["psf_approved"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
