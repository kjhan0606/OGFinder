#!/usr/bin/env python3
"""Audit matched lens/control support separately in each WL radial bin.

This catalog-row audit uses only the already-approved local QA tables. It
checks, for HSC Y3 source bins 2--4, how many matched controls remain supported
in each annulus after the global 4/5-control sample selection. It does not
reload source CSVs, query CAS, access images, estimate shear, or fit masses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
QA_DIR = DATA / "wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01"
SENSITIVITY_DIR = DATA / "wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03"
DEFAULT_OUTPUT = DATA / "wl_s23b_matched_radial_support_2026-11-11_02"

EXPECTED_SHA256 = {
    "qa_summary.json": "eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445",
    "per_center.tsv": "e630442cbaeb054763897a46d68d3b4baea7c0884e80ac22f46ece906583f48e",
    "per_center_radial_tomographic_support.tsv": "f8d604e91c43705cb24152c7e98553efd8838933034a8657bcc603c2f8c11bed",
    "equal_lens_control_profiles.tsv": "9618569b5e4da15dafbff8d1be84e86ab29113bde45a14d6d79d1c314595d97e",
    "sensitivity_summary.json": "9c66b281e0920b79a5c7854e40081b873621867cfaca2a9991d5e93e031f1aa5",
}
EXPECTED_RADIAL_SUPPORT = [60, 67, 69, 73, 74, 74, 74]
EXPECTED_RADIAL_EDGES = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0]
SOURCE_BINS = {2, 3, 4}


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


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def load_inputs(qa_dir: Path, sensitivity_dir: Path) -> tuple[dict, list[dict], list[dict], list[dict]]:
    paths = {
        "qa_summary.json": qa_dir / "summary.json",
        "per_center.tsv": qa_dir / "per_center.tsv",
        "per_center_radial_tomographic_support.tsv": qa_dir / "per_center_radial_tomographic_support.tsv",
        "equal_lens_control_profiles.tsv": sensitivity_dir / "equal_lens_control_profiles.tsv",
        "sensitivity_summary.json": sensitivity_dir / "summary.json",
    }
    for label, path in paths.items():
        observed = sha256(path)
        if observed != EXPECTED_SHA256[label]:
            raise ValueError(f"{label} hash differs from the audited input: {observed}")

    qa_summary = read_json(paths["qa_summary.json"])
    if (qa_summary.get("science_measurement") is not False
            or qa_summary.get("network_used") is not False
            or qa_summary.get("cas_query_submitted") is not False
            or qa_summary.get("additional_data_downloaded") is not False
            or qa_summary.get("images_included") is not False):
        raise ValueError("QA input scope or science flag changed")
    if qa_summary.get("row_plan_sha256") != "7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc":
        raise ValueError("QA input does not match the approved shape-row plan")

    centers = read_tsv(paths["per_center.tsv"])
    radial = read_tsv(paths["per_center_radial_tomographic_support.tsv"])
    sensitivity_rows = read_tsv(paths["equal_lens_control_profiles.tsv"])
    sensitivity_summary = read_json(paths["sensitivity_summary.json"])
    if sensitivity_summary.get("science_measurement") is not False:
        raise ValueError("weight sensitivity input has an unexpected science flag")
    if len(centers) != 480 or len(radial) != 13_440:
        raise ValueError("input QA tables do not have 480 centers and 13,440 radial/z-bin cells")
    if len(sensitivity_rows) != 518:
        raise ValueError("weight sensitivity table is not the pinned 518-row artifact")
    return qa_summary, centers, radial, sensitivity_rows


def global_sample(centers: list[dict[str, str]]) -> tuple[set[str], dict[str, list[str]], set[str]]:
    lens_support: dict[str, bool] = {}
    controls_by_lens: dict[str, list[tuple[str, bool]]] = {}
    for row in centers:
        center_id = row["center_id"]
        matched_lens_id = row["matched_lens_id"]
        has_support = int(row["n_calibration_valid"]) > 0
        if row["role"] == "lens":
            if center_id != matched_lens_id or center_id in lens_support:
                raise ValueError("invalid or duplicated lens center in support table")
            lens_support[center_id] = has_support
            controls_by_lens[center_id] = []
        elif row["role"] != "control":
            raise ValueError(f"unexpected center role: {row['role']}")
    for row in centers:
        if row["role"] == "control":
            lens_id = row["matched_lens_id"]
            if lens_id not in controls_by_lens:
                raise ValueError(f"control has unknown matched lens: {lens_id}")
            controls_by_lens[lens_id].append((row["center_id"], int(row["n_calibration_valid"]) > 0))
    if len(lens_support) != 80 or sum(len(rows) for rows in controls_by_lens.values()) != 400:
        raise ValueError("support table must contain 80 lenses and five controls per lens")

    eligible_lenses = {
        lens_id for lens_id, supported in lens_support.items()
        if supported and sum(control_supported for _, control_supported in controls_by_lens[lens_id]) >= 4
    }
    selected_controls = {
        control_id
        for lens_id in eligible_lenses
        for control_id, supported in controls_by_lens[lens_id]
        if supported
    }
    lens_hash = hashlib.sha256("\n".join(sorted(eligible_lenses)).encode()).hexdigest()
    control_hash = hashlib.sha256("\n".join(sorted(selected_controls)).encode()).hexdigest()
    if (len(eligible_lenses) != 74 or len(selected_controls) != 367
            or lens_hash != "78e35dcab6bc000160d391ae3b16516bc342ee9f74a550a6106bb256afef7a08"
            or control_hash != "583fd154ebb799338dd354eb3e043583f9247cef87200323c5574f21bc6c5894"):
        raise ValueError("global 4/5 sample membership differs from the pinned S23B profile")
    return eligible_lenses, {
        lens_id: [control_id for control_id, supported in rows if supported]
        for lens_id, rows in controls_by_lens.items()
    }, selected_controls


def radial_support_rows(radial_cells: list[dict[str, str]], eligible_lenses: set[str],
                        sensitivity_rows: list[dict[str, str]]) -> tuple[list[dict], list[dict]]:
    centers: dict[tuple[str, str, int], dict] = {}
    expected_cells: set[tuple[str, str, int, int]] = set()
    for row in radial_cells:
        center_id = row["center_id"]
        matched_lens_id = row["matched_lens_id"]
        role = row["role"]
        radial_bin = int(row["radial_bin"])
        zbin = int(row["source_zbin"])
        key4 = (role, center_id, radial_bin, zbin)
        if key4 in expected_cells:
            raise ValueError(f"duplicate center/radius/z-bin cell: {key4}")
        expected_cells.add(key4)
        if not 0 <= radial_bin < 7 or not 1 <= zbin <= 4:
            raise ValueError("unexpected radial or tomographic bin index")
        if role == "lens" and center_id != matched_lens_id:
            raise ValueError("lens center ID does not match its set ID")
        key = (role, center_id, radial_bin)
        cell = centers.setdefault(key, {
            "matched_lens_id": matched_lens_id,
            "r_min_mpc": float(row["r_min_mpc"]),
            "r_max_mpc": float(row["r_max_mpc"]),
            "n_rows_z234": 0,
            "sum_shape_weight_z234": 0.0,
            "sum_geometry_proxy_z234": 0.0,
        })
        if (cell["matched_lens_id"] != matched_lens_id
                or cell["r_min_mpc"] != float(row["r_min_mpc"])
                or cell["r_max_mpc"] != float(row["r_max_mpc"])):
            raise ValueError(f"inconsistent center metadata in {key}")
        if zbin in SOURCE_BINS:
            n_rows = int(row["n_rows"])
            shape_weight = float(row["sum_shape_weight"])
            geometry_proxy = float(row["sum_sigma_crit_weight_proxy"])
            if n_rows < 0 or shape_weight < 0 or geometry_proxy < 0:
                raise ValueError("negative row or support weight in QA input")
            cell["n_rows_z234"] += n_rows
            cell["sum_shape_weight_z234"] += shape_weight
            cell["sum_geometry_proxy_z234"] += geometry_proxy
    if len(expected_cells) != 480 * 7 * 4 or len(centers) != 480 * 7:
        raise ValueError("radial support table does not cover every center × radius × z-bin")

    table_rows: list[dict] = []
    summary_rows: list[dict] = []
    all_lenses = {row["center_id"] for row in radial_cells if row["role"] == "lens"}
    control_ids_by_lens: dict[str, list[str]] = {}
    for row in radial_cells:
        if row["role"] == "control":
            control_ids_by_lens.setdefault(row["matched_lens_id"], [])
            if row["center_id"] not in control_ids_by_lens[row["matched_lens_id"]]:
                control_ids_by_lens[row["matched_lens_id"]].append(row["center_id"])
    for lens_id, control_ids in control_ids_by_lens.items():
        if len(control_ids) != 5:
            raise ValueError(f"matched lens {lens_id} does not have five controls")

    for radial_bin in range(7):
        lens_ids_with_any_support = []
        n_control_center_total = 0
        n_control_center_selected = 0
        set_rows = []
        for lens_id in sorted(all_lenses):
            lens = centers[("lens", lens_id, radial_bin)]
            controls = [centers[("control", control_id, radial_bin)]
                        for control_id in control_ids_by_lens[lens_id]]
            positive = [row for row in controls if row["n_rows_z234"] > 0]
            n_control_center_total += len(positive)
            if lens_id in eligible_lenses:
                n_control_center_selected += len(positive)
                if lens["n_rows_z234"] > 0:
                    lens_ids_with_any_support.append(lens_id)
            if lens_id not in eligible_lenses:
                continue
            lens_count = int(lens["n_rows_z234"])
            control_count = len(positive)
            mean_control_rows = (statistics.mean(row["n_rows_z234"] for row in positive)
                                  if positive else None)
            mean_control_shape = (statistics.mean(row["sum_shape_weight_z234"] for row in positive)
                                  if positive else None)
            mean_control_proxy = (statistics.mean(row["sum_geometry_proxy_z234"] for row in positive)
                                  if positive else None)
            row = {
                "matched_lens_id": lens_id,
                "radial_bin": radial_bin,
                "r_min_mpc": lens["r_min_mpc"],
                "r_max_mpc": lens["r_max_mpc"],
                "lens_pair_rows_z234": lens_count,
                "n_positive_controls_z234": control_count,
                "control_pair_rows_z234_total": sum(item["n_rows_z234"] for item in controls),
                "mean_positive_control_pair_rows_z234": mean_control_rows,
                "lens_shape_weight_sum_z234": lens["sum_shape_weight_z234"],
                "mean_positive_control_shape_weight_sum_z234": mean_control_shape,
                "lens_geometry_proxy_sum_z234": lens["sum_geometry_proxy_z234"],
                "mean_positive_control_geometry_proxy_sum_z234": mean_control_proxy,
                "control_to_lens_pair_count_ratio": mean_control_rows / lens_count
                if lens_count > 0 and mean_control_rows is not None else None,
                "control_to_lens_shape_weight_ratio": mean_control_shape / lens["sum_shape_weight_z234"]
                if lens["sum_shape_weight_z234"] > 0 and mean_control_shape is not None else None,
                "control_to_lens_geometry_proxy_ratio": mean_control_proxy / lens["sum_geometry_proxy_z234"]
                if lens["sum_geometry_proxy_z234"] > 0 and mean_control_proxy is not None else None,
            }
            set_rows.append(row)
            table_rows.append(row)

        selected_rows = [row for row in set_rows if row["lens_pair_rows_z234"] > 0]
        thresholds = {}
        for minimum in (1, 3, 4, 5):
            eligible = [row for row in selected_rows if row["n_positive_controls_z234"] >= minimum]
            thresholds[str(minimum)] = {
                "n_lens_sets": len(eligible),
                "n_positive_control_centers": sum(row["n_positive_controls_z234"] for row in eligible),
                "lens_pair_rows": sum(row["lens_pair_rows_z234"] for row in eligible),
                "median_control_to_lens_pair_count_ratio": quantile(
                    [row["control_to_lens_pair_count_ratio"] for row in eligible
                     if row["control_to_lens_pair_count_ratio"] is not None], 0.5),
            }
        any_common = [row for row in selected_rows if row["n_positive_controls_z234"] > 0]
        summary_rows.append({
            "radial_bin": radial_bin,
            "r_min_mpc": centers[("lens", sorted(all_lenses)[0], radial_bin)]["r_min_mpc"],
            "r_max_mpc": centers[("lens", sorted(all_lenses)[0], radial_bin)]["r_max_mpc"],
            "global_selected_lens_sets": len(eligible_lenses),
            "selected_lens_sets_with_z234_support": len(lens_ids_with_any_support),
            "selected_control_centers_with_z234_support": n_control_center_selected,
            "all_80_lens_sets_with_z234_support": sum(
                centers[("lens", lens_id, radial_bin)]["n_rows_z234"] > 0 for lens_id in all_lenses),
            "all_400_control_centers_with_z234_support": n_control_center_total,
            "matched_sets_with_both_roles_supported": len(any_common),
            "radial_control_thresholds_within_global_74": thresholds,
        })

    observed_counts = [row["matched_sets_with_both_roles_supported"] for row in summary_rows]
    reference_counts = []
    sensitivity_by_bin: dict[int, set[str]] = {}
    for row in sensitivity_rows:
        radial_bin = int(row["radial_bin"])
        if row["lens_delta_sigma_msun_mpc2"] and row["random_delta_sigma_msun_mpc2"]:
            sensitivity_by_bin.setdefault(radial_bin, set()).add(row["matched_lens_id"])
    reference_counts = [len(sensitivity_by_bin.get(index, set())) for index in range(7)]
    if observed_counts != EXPECTED_RADIAL_SUPPORT or reference_counts != EXPECTED_RADIAL_SUPPORT:
        raise ValueError(
            "radial common-support counts do not reconcile: "
            f"QA={observed_counts}, sensitivity={reference_counts}, "
            f"expected={EXPECTED_RADIAL_SUPPORT}"
        )
    for row in summary_rows:
        radial_bin = int(row["radial_bin"])
        row["common_support_crosscheck_n"] = reference_counts[radial_bin]
    return table_rows, summary_rows


def fixed_all_radius_support(per_lens_rows: list[dict]) -> dict:
    """Count lens sets that retain the requested control support in every bin."""
    by_lens: dict[str, list[dict]] = {}
    for row in per_lens_rows:
        by_lens.setdefault(str(row["matched_lens_id"]), []).append(row)
    result = {}
    for minimum in (1, 3, 4, 5):
        eligible = sorted(
            lens_id for lens_id, rows in by_lens.items()
            if len(rows) == 7
            and all(int(row["lens_pair_rows_z234"]) > 0
                    and int(row["n_positive_controls_z234"]) >= minimum
                    for row in rows)
        )
        result[str(minimum)] = {
            "minimum_positive_controls_in_each_annulus": minimum,
            "requires_positive_lens_rows_in_each_annulus": True,
            "n_lens_sets_supported_in_all_seven_annuli": len(eligible),
            "lens_ids_sha256": hashlib.sha256("\n".join(eligible).encode()).hexdigest(),
        }
    return result


def write_assessment(path: Path, summary: dict) -> None:
    lines = [
        "# S23B matched-control radial support audit (2026-11-11)",
        "",
        "This offline audit reaggregates only the previously verified per-center QA table for HSC Y3 source bins 2–4. It reads no raw source CSVs or images and makes no network/CAS request. It does not calculate shear or mass; `science_measurement=false`.",
        "",
        "The parent analysis selects 74 lens sets and 367 positive-support controls using integrated 0.1–3 proper-Mpc support and the 4/5-control rule. This audit asks whether each annulus itself retains the same matched-control support.",
        "",
        "| R proper (Mpc) | 74-set lenses with z2–4 rows | 74-set controls with z2–4 rows | Both roles supported | ≥3 controls | ≥4 controls | ≥5 controls |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["radial_summary"]:
        thresholds = row["radial_control_thresholds_within_global_74"]
        lines.append(
            f"| {row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f} | "
            f"{row['selected_lens_sets_with_z234_support']}/74 | "
            f"{row['selected_control_centers_with_z234_support']}/367 | "
            f"{row['matched_sets_with_both_roles_supported']} | "
            f"{thresholds['3']['n_lens_sets']} | {thresholds['4']['n_lens_sets']} | "
            f"{thresholds['5']['n_lens_sets']} |"
        )
    lines.extend([
        "",
        "Counts are center/source pair rows after the existing shared cuts; they are not unique source-object counts and are not calibrated shear weights. `sum_shape_weight` and the `D_l beta` geometry proxy are retained only as support diagnostics.",
        "",
        "The per-radius `≥4` and `≥5` columns are stricter annulus-specific diagnostics inside the globally selected 74-set sample. They do not silently replace the sample rule used in the existing profile. The primary equal-lens estimator must state whether it uses radial common support or requires a fixed support set across all bins; that choice affects the estimand and covariance.",
        "",
        "The baseline `both roles supported` counts are cross-checked against the equal-lens sensitivity table and must match [60, 67, 69, 73, 74, 74, 74]. This audit does not choose among the annulus-specific thresholds or promote the profile. Keep `science_measurement=false`.",
        "",
        "A fixed lens-set profile across all seven annuli requires positive lens rows and the stated number of positive-control centers in every annulus. The 4/5-control fixed-set intersection retains "
        f"{summary['fixed_all_radius_support']['4']['n_lens_sets_supported_in_all_seven_annuli']} of the 74 globally selected sets. This is a separate support sensitivity; the radial-varying sample and fixed sample estimate different population profiles.",
        "",
        "## Provenance",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- QA summary SHA-256: `{summary['provenance']['qa_summary_sha256']}`",
        f"- Per-center support SHA-256: `{summary['provenance']['per_center_sha256']}`",
        f"- Per-center radial/z-bin table SHA-256: `{summary['provenance']['radial_support_sha256']}`",
        f"- Equal-lens sensitivity table SHA-256: `{summary['provenance']['equal_lens_table_sha256']}`",
        f"- Output directory: `{summary['output_directory']}`",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def run(qa_dir: Path, sensitivity_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    qa_summary, centers, radial_cells, sensitivity_rows = load_inputs(qa_dir, sensitivity_dir)
    eligible_lenses, controls_by_lens, selected_controls = global_sample(centers)
    per_lens, radial_summary = radial_support_rows(radial_cells, eligible_lenses, sensitivity_rows)

    runner = Path(__file__).resolve()
    summary = {
        "schema": "astrafek-s23b-matched-radial-control-support-v1",
        "status": "complete_local_radial_matched_support_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "shear_profile_calculated": False,
        "halo_mass_fit_performed": False,
        "scope": "per-annulus matched-control support for the existing S23B 80-lens/400-control shape-row QA, HSC Y3 source bins 2-4",
        "global_sample": {
            "input_lens_sets": 80,
            "input_control_centers": 400,
            "selection_rule": "positive integrated lens support over 0.1-3 proper Mpc and at least four of five globally positive matched controls",
            "selected_lens_sets": len(eligible_lenses),
            "selected_positive_control_centers": len(selected_controls),
            "lens_id_sha256": hashlib.sha256("\n".join(sorted(eligible_lenses)).encode()).hexdigest(),
            "control_id_sha256": hashlib.sha256("\n".join(sorted(selected_controls)).encode()).hexdigest(),
        },
        "source_bins": sorted(SOURCE_BINS),
        "radial_edges_proper_mpc": EXPECTED_RADIAL_EDGES,
        "radial_summary": radial_summary,
        "fixed_all_radius_support": fixed_all_radius_support(per_lens),
        "counts_are_center_source_pairs": True,
        "limitations": [
            "pair rows can repeat the same source across lens/control apertures",
            "this is support QA, not a calibrated shear estimate or covariance",
            "annulus-specific minimum-control rules are diagnostics; they are not selected policies",
            "shape-weight and geometry-proxy sums omit the final estimator response and are not estimator weights",
        ],
        "provenance": {
            "runner_path": str(runner),
            "runner_sha256": sha256(runner),
            "qa_summary_path": str(qa_dir / "summary.json"),
            "qa_summary_sha256": sha256(qa_dir / "summary.json"),
            "per_center_path": str(qa_dir / "per_center.tsv"),
            "per_center_sha256": sha256(qa_dir / "per_center.tsv"),
            "radial_support_path": str(qa_dir / "per_center_radial_tomographic_support.tsv"),
            "radial_support_sha256": sha256(qa_dir / "per_center_radial_tomographic_support.tsv"),
            "equal_lens_sensitivity_path": str(sensitivity_dir / "equal_lens_control_profiles.tsv"),
            "equal_lens_table_sha256": sha256(sensitivity_dir / "equal_lens_control_profiles.tsv"),
            "sensitivity_summary_sha256": sha256(sensitivity_dir / "summary.json"),
        },
        "output_directory": str(output),
    }
    output.mkdir(parents=True, exist_ok=False)
    fields = list(per_lens[0])
    write_tsv(output / "per_lens_radial_support.tsv", per_lens, fields)
    write_assessment(output / "assessment.md", summary)
    summary["output_sha256"] = {
        name: sha256(output / name)
        for name in ("assessment.md", "per_lens_radial_support.tsv")
    }
    write_json(output / "summary.json", summary)
    summary_sha = sha256(output / "summary.json")
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "selected_lens_sets": len(eligible_lenses),
        "selected_control_centers": len(selected_controls),
        "radial_support": [
            {"r_mpc": [row["r_min_mpc"], row["r_max_mpc"]],
             "both_roles": row["matched_sets_with_both_roles_supported"],
             "n_ge3": row["radial_control_thresholds_within_global_74"]["3"]["n_lens_sets"],
             "n_ge4": row["radial_control_thresholds_within_global_74"]["4"]["n_lens_sets"],
             "n_ge5": row["radial_control_thresholds_within_global_74"]["5"]["n_lens_sets"]}
            for row in radial_summary
        ],
        "fixed_all_radius_support": summary["fixed_all_radius_support"],
        "output_directory": str(output),
        "summary_sha256": summary_sha,
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa-dir", type=Path, default=QA_DIR)
    parser.add_argument("--sensitivity-dir", type=Path, default=SENSITIVITY_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.qa_dir, args.sensitivity_dir, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
