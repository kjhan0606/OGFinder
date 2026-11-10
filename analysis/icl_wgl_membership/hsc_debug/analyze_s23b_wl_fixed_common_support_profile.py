#!/usr/bin/env python3
"""Build a fixed-sample S23B lens-weighted profile from audited per-set stacks.

This reweights existing calibrated lens-set profile rows to the subset that
has lens support and at least four positive-support random centers in every
one of seven radial annuli. It reads only prior Astrafek TSV/JSON outputs; it
does not reload source CSVs, query CAS, download data, access images, estimate
covariance, or fit a mass.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_DIR = DATA / "wl_s23b_matched_radial_support_2026-11-11_02"
SENSITIVITY_DIR = DATA / "wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03"
QA_DIR = DATA / "wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01"
DEFAULT_OUTPUT = DATA / "wl_s23b_fixed_common_support_profile_2026-11-11_02"
EXPECTED = {
    "support_summary.json": "c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234",
    "per_lens_radial_support.tsv": "6744ece545781faffebef65fc7a0143284da3a234425b6bed540bd1e9d19625c",
    "sensitivity_summary.json": "9c66b281e0920b79a5c7854e40081b873621867cfaca2a9991d5e93e031f1aa5",
    "equal_lens_control_profiles.tsv": "9618569b5e4da15dafbff8d1be84e86ab29113bde45a14d6d79d1c314595d97e",
    "per_center.tsv": "e630442cbaeb054763897a46d68d3b4baea7c0884e80ac22f46ece906583f48e",
}
EXPECTED_FIXED_LENS_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"
EXPECTED_RADIAL_COUNTS = [60, 67, 69, 73, 74, 74, 74]
EXPECTED_FIXED_N = 52
SCALE = 1.0e13


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path: Path, payload: dict) -> None:
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def run(support_dir: Path, sensitivity_dir: Path, qa_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    paths = {
        "support_summary.json": support_dir / "summary.json",
        "per_lens_radial_support.tsv": support_dir / "per_lens_radial_support.tsv",
        "sensitivity_summary.json": sensitivity_dir / "summary.json",
        "equal_lens_control_profiles.tsv": sensitivity_dir / "equal_lens_control_profiles.tsv",
        "per_center.tsv": qa_dir / "per_center.tsv",
    }
    for label, path in paths.items():
        observed = sha256(path)
        if observed != EXPECTED[label]:
            raise ValueError(f"{label} input hash mismatch: {observed}")

    support = read_json(paths["support_summary.json"])
    sensitivity = read_json(paths["sensitivity_summary.json"])
    if support.get("science_measurement") is not False or sensitivity.get("science_measurement") is not False:
        raise ValueError("input record unexpectedly promotes the analysis to a science measurement")
    if support.get("status") != "complete_local_radial_matched_support_audit":
        raise ValueError("radial support audit is incomplete")

    fixed_design = support.get("fixed_all_radius_support", {}).get("4", {})
    if (fixed_design.get("n_lens_sets_supported_in_all_seven_annuli") != EXPECTED_FIXED_N
            or fixed_design.get("lens_ids_sha256") != EXPECTED_FIXED_LENS_SHA256):
        raise ValueError("fixed ≥4-control all-radius sample membership changed")

    support_rows = read_tsv(paths["per_lens_radial_support.tsv"])
    grouped_support = {}
    for row in support_rows:
        grouped_support.setdefault(row["matched_lens_id"], []).append(row)
    fixed_lens_ids = {
        lens_id for lens_id, rows in grouped_support.items()
        if len(rows) == 7 and all(
            int(row["n_positive_controls_z234"]) >= 4
            and int(row["lens_pair_rows_z234"]) > 0 for row in rows
        )
    }
    fixed_hash = hashlib.sha256("\n".join(sorted(fixed_lens_ids)).encode()).hexdigest()
    if len(fixed_lens_ids) != EXPECTED_FIXED_N or fixed_hash != EXPECTED_FIXED_LENS_SHA256:
        raise ValueError("fixed sample derived from per-lens table does not match its summary")

    center_rows = read_tsv(paths["per_center.tsv"])
    if len(center_rows) != 480:
        raise ValueError("per-center QA table no longer has 480 matched apertures")
    lens_metadata = {row["center_id"]: row for row in center_rows if row["role"] == "lens"}
    positive_controls: dict[str, int] = {}
    for row in center_rows:
        if row["role"] == "control" and int(row["n_calibration_valid"]) > 0:
            positive_controls[row["matched_lens_id"]] = positive_controls.get(row["matched_lens_id"], 0) + 1
    global_lens_ids = {
        lens_id for lens_id, row in lens_metadata.items()
        if int(row["n_calibration_valid"]) > 0 and positive_controls.get(lens_id, 0) >= 4
    }
    global_hash = hashlib.sha256("\n".join(sorted(global_lens_ids)).encode()).hexdigest()
    if (len(global_lens_ids) != 74
            or global_hash != support["global_sample"]["lens_id_sha256"]
            or not fixed_lens_ids.issubset(global_lens_ids)):
        raise ValueError("fixed and global lens samples do not match their pinned membership")
    composition = {}
    groups = {
        "global_74": global_lens_ids,
        "fixed_52": fixed_lens_ids,
        "outside_fixed_22": global_lens_ids - fixed_lens_ids,
    }
    for label, ids in groups.items():
        selected = [lens_metadata[lens_id] for lens_id in ids]
        composition[label] = {
            "n_lens_sets": len(selected),
            "z_lens_q16_median_q84": {
                key: quantile([float(row["z_lens"]) for row in selected], probability)
                for key, probability in (("q16", 0.16), ("median", 0.5), ("q84", 0.84))
            },
            "retained_catalog_pair_rows_q16_median_q84": {
                key: quantile([float(row["n_calibration_valid"]) for row in selected], probability)
                for key, probability in (("q16", 0.16), ("median", 0.5), ("q84", 0.84))
            },
            "catalog_shape_weight_sum_q16_median_q84": {
                key: quantile([float(row["sum_shape_weight"]) for row in selected], probability)
                for key, probability in (("q16", 0.16), ("median", 0.5), ("q84", 0.84))
            },
        }

    group_rows = read_tsv(paths["equal_lens_control_profiles.tsv"])
    if len(group_rows) != 518:
        raise ValueError("equal-lens group table no longer has its pinned 518 rows")
    by_radial: dict[int, list[dict[str, str]]] = {index: [] for index in range(7)}
    for row in group_rows:
        radial_bin = int(row["radial_bin"])
        if not 0 <= radial_bin < 7:
            raise ValueError("unexpected radial-bin index")
        if row["matched_lens_id"] in fixed_lens_ids:
            by_radial[radial_bin].append(row)

    radial_summary = []
    profile_rows = []
    for radial_bin in range(7):
        fixed_rows = by_radial[radial_bin]
        if len(fixed_rows) != EXPECTED_FIXED_N:
            raise ValueError(f"fixed profile has {len(fixed_rows)} rows in radial bin {radial_bin}")
        paired_t = []
        paired_x = []
        varying_t = []
        varying_x = []
        pooled = []
        for row in fixed_rows:
            required = ("lens_delta_sigma_msun_mpc2", "random_delta_sigma_msun_mpc2",
                        "lens_cross_msun_mpc2", "random_cross_msun_mpc2")
            if any(row[name] == "" for name in required):
                raise ValueError("fixed all-radius sample contains a missing per-set profile value")
            paired_t.append(float(row["lens_delta_sigma_msun_mpc2"])
                            - float(row["random_delta_sigma_msun_mpc2"]))
            paired_x.append(float(row["lens_cross_msun_mpc2"])
                            - float(row["random_cross_msun_mpc2"]))

        support_sets = [row for row in group_rows if int(row["radial_bin"]) == radial_bin
                        and row["lens_delta_sigma_msun_mpc2"] != ""
                        and row["random_delta_sigma_msun_mpc2"] != ""]
        if len(support_sets) != EXPECTED_RADIAL_COUNTS[radial_bin]:
            raise ValueError("radial-varying comparison support differs from the pinned audit")
        varying_t = [float(row["lens_delta_sigma_msun_mpc2"])
                     - float(row["random_delta_sigma_msun_mpc2"]) for row in support_sets]
        varying_x = [float(row["lens_cross_msun_mpc2"])
                     - float(row["random_cross_msun_mpc2"]) for row in support_sets]

        sensitivity_profile = sensitivity["equal_lens_weight"]["profile"][radial_bin]
        pooled_value = float(sensitivity_profile["pooled_delta_sigma_msun_mpc2"])
        sensitivity_mean = float(sensitivity_profile["equal_paired_lens_set_delta_sigma_msun_mpc2"])
        varying_mean = statistics.mean(varying_t)
        if not (abs(varying_mean - sensitivity_mean) <= max(1e-5, abs(sensitivity_mean) * 1e-10)):
            raise ValueError("recomputed radial-varying mean disagrees with sensitivity summary")

        fixed_mean = statistics.mean(paired_t)
        fixed_cross_mean = statistics.mean(paired_x)
        radial = fixed_rows[0]
        profile_rows.append({
            "radial_bin": radial_bin,
            "r_min_mpc": float(radial["r_min_mpc"]),
            "r_max_mpc": float(radial["r_max_mpc"]),
            "n_fixed_lens_sets": len(fixed_rows),
            "delta_sigma_fixed_equal_lens_1e13_msun_mpc2": fixed_mean / SCALE,
            "cross_fixed_equal_lens_1e13_msun_mpc2": fixed_cross_mean / SCALE,
            "delta_sigma_per_set_sd_descriptive_1e13_msun_mpc2": statistics.stdev(paired_t) / SCALE,
            "delta_sigma_per_set_q16_1e13_msun_mpc2": quantile(paired_t, 0.16) / SCALE,
            "delta_sigma_per_set_median_1e13_msun_mpc2": quantile(paired_t, 0.50) / SCALE,
            "delta_sigma_per_set_q84_1e13_msun_mpc2": quantile(paired_t, 0.84) / SCALE,
            "n_radial_varying_lens_sets": len(support_sets),
            "delta_sigma_radial_varying_equal_lens_1e13_msun_mpc2": varying_mean / SCALE,
            "delta_sigma_fixed_vs_radial_varying_fractional_change": fixed_mean / varying_mean - 1.0
            if varying_mean else None,
            "delta_sigma_common_support_pooled_1e13_msun_mpc2": pooled_value / SCALE,
        })
        radial_summary.append({
            "radial_bin": radial_bin,
            "r_min_mpc": float(radial["r_min_mpc"]),
            "r_max_mpc": float(radial["r_max_mpc"]),
            "n_lens_sets": len(fixed_rows),
            "delta_sigma_mean_msun_mpc2": fixed_mean,
            "cross_delta_sigma_mean_msun_mpc2": fixed_cross_mean,
            "per_set_delta_sigma_sd_msun_mpc2_descriptive": statistics.stdev(paired_t),
            "radial_varying_sample_mean_msun_mpc2": varying_mean,
            "radial_varying_sample_mean_cross_msun_mpc2": statistics.mean(varying_x),
            "radial_varying_sample_n_lens_sets": len(support_sets),
            "common_support_pooled_delta_sigma_msun_mpc2": pooled_value,
        })

    runner = Path(__file__).resolve()
    summary = {
        "schema": "astrafek-s23b-wl-fixed-common-support-profile-v1",
        "status": "complete_exploratory_fixed_sample_profile_reweighting",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "source_rows_reloaded": False,
        "covariance_estimated": False,
        "halo_mass_fit_performed": False,
        "scope": "equal lens-set mean of previously calibrated per-lens and matched-random profile rows for one fixed all-radius support sample",
        "selection": {
            "global_sample_size": 74,
            "fixed_sample_size": len(fixed_lens_ids),
            "fixed_sample_rule": "positive lens source rows and at least four positive-support matched controls in each of seven proper-radius annuli, source bins 2-4",
            "fixed_sample_lens_id_sha256": fixed_hash,
            "radial_varying_sample_size_by_bin": EXPECTED_RADIAL_COUNTS,
        },
        "sample_composition": composition,
        "profile": profile_rows,
        "interpretation": [
            "values are equal-weight means of previously calibrated per-center lens-minus-random profile estimates",
            "per-set scatter is descriptive and is not a profile uncertainty or covariance",
            "fixed sample is a support-selected sensitivity subset and may not represent the full 74-lens sample",
            "the profile is not promoted to a science measurement",
        ],
        "provenance": {
            "runner_path": str(runner),
            "runner_sha256": sha256(runner),
            "support_summary_sha256": sha256(paths["support_summary.json"]),
            "per_lens_radial_support_sha256": sha256(paths["per_lens_radial_support.tsv"]),
            "sensitivity_summary_sha256": sha256(paths["sensitivity_summary.json"]),
            "equal_lens_control_profiles_sha256": sha256(paths["equal_lens_control_profiles.tsv"]),
            "per_center_qa_sha256": sha256(paths["per_center.tsv"]),
            "fixed_sample_lens_id_sha256": fixed_hash,
        },
        "output_directory": str(output),
        "radial_profile_msun_mpc2": radial_summary,
    }
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "profile.tsv", profile_rows, list(profile_rows[0]))
    assessment = [
        "# S23B fixed common-support weak-lensing profile (2026-11-11)",
        "",
        "This diagnostic reweights previously calculated per-lens/per-random calibrated profile rows. Within each lens set, the positive random apertures remain pooled by source-pair lensing weights; only the outer mean across lens sets is equal-weighted. It does not reopen source CSVs or images, query CAS, download data, estimate covariance, or fit mass. `science_measurement=false`.",
        "",
        "The fixed sample contains 52 lens sets with positive target rows and at least four positive-support controls in every one of seven annuli (0.1–3 proper Mpc), using HSC Y3 source bins 2–4. The broader 4/5 globally selected sample has 74 sets, but its per-radius common support changes from 60 to 74 sets.",
        "",
        "| R proper (Mpc) | Fixed N | Fixed equal-lens ΔΣt | Radial-varying ΔΣt | Fixed vs varying | Cross ΔΣ | Per-set SD (descriptive) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in profile_rows:
        lines = (f"{row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f}", f"{row['n_fixed_lens_sets']}",
                 f"{row['delta_sigma_fixed_equal_lens_1e13_msun_mpc2']:.4f}",
                 f"{row['delta_sigma_radial_varying_equal_lens_1e13_msun_mpc2']:.4f}",
                 f"{100.0*row['delta_sigma_fixed_vs_radial_varying_fractional_change']:+.2f}%",
                 f"{row['cross_fixed_equal_lens_1e13_msun_mpc2']:.4f}",
                 f"{row['delta_sigma_per_set_sd_descriptive_1e13_msun_mpc2']:.4f}")
        assessment.append("| " + " | ".join(lines) + " |")
    assessment.extend([
        "",
        "All profile and cross-profile values are in `10^13 Msun Mpc^-2`. The per-set standard deviation describes cluster-to-cluster spread only; it is not an error bar. No spatial covariance was estimated, so neither a significance nor a halo mass is reported.",
        "",
        "The fixed sample changes the radial trend's target population compared with annulus-specific support. This table is a selection sensitivity only. The next inferential step needs a paired spatial covariance for the same fixed sample, with repeated sources deleted together and block-scale/footprint sensitivity assessed.",
        "",
        "| Lens-set sample | N | z_lens median (q16–q84) | Retained catalog rows median (q16–q84) | Shape-weight sum median (q16–q84) |",
        "|---|---:|---:|---:|---:|",
    ])
    for label, display in (("global_74", "Global 74"), ("fixed_52", "Fixed 52"),
                           ("outside_fixed_22", "Outside fixed 22")):
        sample = composition[label]
        z = sample["z_lens_q16_median_q84"]
        rows = sample["retained_catalog_pair_rows_q16_median_q84"]
        weights = sample["catalog_shape_weight_sum_q16_median_q84"]
        assessment.append(
            f"| {display} | {sample['n_lens_sets']} | {z['median']:.4f} ({z['q16']:.4f}–{z['q84']:.4f}) | "
            f"{rows['median']:.0f} ({rows['q16']:.0f}–{rows['q84']:.0f}) | "
            f"{weights['median']:.1f} ({weights['q16']:.1f}–{weights['q84']:.1f}) |"
        )
    assessment.extend([
        "",
        "The fixed sample's median lens redshift is close to the global 74-set sample, but its median retained source-pair count and catalog shape-weight sum are higher. This confirms that the all-radius requirement selects better-supported apertures; redshift similarity alone does not establish representativeness or matched image conditions.",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Fixed-support audit summary SHA-256: `{summary['provenance']['support_summary_sha256']}`",
        f"- Equal-lens sensitivity summary SHA-256: `{summary['provenance']['sensitivity_summary_sha256']}`",
        f"- Per-center QA SHA-256: `{summary['provenance']['per_center_qa_sha256']}`",
        f"- Output directory: `{output}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    summary["output_sha256"] = {
        name: sha256(output / name) for name in ("profile.tsv", "assessment.md")
    }
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "fixed_lens_sets": len(fixed_lens_ids),
        "fixed_lens_id_sha256": fixed_hash,
        "profile": profile_rows,
        "output_directory": str(output),
        "summary_sha256": sha256(output / "summary.json"),
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-dir", type=Path, default=SUPPORT_DIR)
    parser.add_argument("--sensitivity-dir", type=Path, default=SENSITIVITY_DIR)
    parser.add_argument("--qa-dir", type=Path, default=QA_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.support_dir, args.sensitivity_dir, args.qa_dir, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
