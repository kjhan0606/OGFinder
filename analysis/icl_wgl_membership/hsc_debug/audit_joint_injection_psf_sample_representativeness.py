#!/usr/bin/env python3
"""Compare the PSF-supported joint-injection sample with its 224-pair parent."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INJECTION = DEFAULT_DATA / "icl_satellite_injection_amplitude_grid_2026-11-05_02/paired_amplitude_frames.tsv"
DEFAULT_ICL = DEFAULT_DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
DEFAULT_SATELLITE = DEFAULT_DATA / "satellite_outer_profile_perbin_cov040_2026-11-05_01/satellite_profiles.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def pair_key(row: dict[str, str], id_field: str) -> tuple[str, str]:
    return str(row[id_field]), str(row["filter"])


def metrics(rows: list[dict[str, str]]) -> dict:
    if not rows:
        raise ValueError("cannot summarize an empty pair sample")
    result = {"n_pairs": len(rows), "n_unique_clusters": len({row["id"] for row in rows})}
    for column in ("z", "richness", "coverage_200_300"):
        values = np.asarray([float(row[column]) for row in rows], dtype=float)
        result[column] = {
            "mean": float(values.mean()),
            "median": float(np.median(values)),
            "std_sample": float(values.std(ddof=1)) if len(values) > 1 else None,
        }
    return result


def support_count(rows: list[dict[str, str]], supported: set[tuple[str, str]]) -> dict:
    n = len(rows)
    selected = sum(pair_key(row, "id") in supported for row in rows)
    return {"n_supported": selected, "n_total": n, "fraction": selected / n if n else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--injection-frames", type=Path, default=DEFAULT_INJECTION)
    parser.add_argument("--icl-support", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    injection_rows = read_tsv(args.injection_frames)
    icl_rows_all = read_tsv(args.icl_support)
    satellite_rows = read_tsv(args.satellite_profiles)
    if len(injection_rows) != 57 or len({pair_key(row, "cluster_id") for row in injection_rows}) != 57:
        raise ValueError("expected 57 unique PSF-supported joint-injection cluster/filter pairs")
    if not {"id", "filter", "z", "richness", "coverage_200_300"}.issubset(icl_rows_all[0]):
        raise ValueError("ICL support table is missing required fields")
    parent = [row for row in icl_rows_all if float(row["coverage_200_300"]) >= 0.40]
    parent_by_key = {pair_key(row, "id"): row for row in parent}
    if len(parent) != 224 or len(parent_by_key) != 224:
        raise ValueError(f"expected 224 unique ICL-supported pairs, found {len(parent_by_key)}")
    psf_keys = {pair_key(row, "cluster_id") for row in injection_rows}
    if not psf_keys.issubset(parent_by_key):
        raise ValueError("PSF-supported injection sample contains pairs outside the parent support")

    sat60: set[tuple[str, str]] = set()
    sat40: set[tuple[str, str]] = set()
    for row in satellite_rows:
        if row.get("radial_bin_kpc") != "50-75":
            continue
        key = pair_key(row, "cluster_id")
        target, control = float(row["target_coverage"]), float(row["control_coverage"])
        if target >= 0.60 and control >= 0.60:
            sat60.add(key)
        if target >= 0.40 and control >= 0.40:
            sat40.add(key)

    grouped = {"all": parent}
    for filt in sorted({row["filter"] for row in parent}):
        grouped[filt] = [row for row in parent if row["filter"] == filt]
    report = {
        "schema": "astrafek-joint-injection-psf-sample-representativeness-v1",
        "status": "complete_offline_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "parent_definition": "224 unique cluster/filter pairs with ICL 200-300 kpc common coverage >= 0.40",
        "selected_definition": "intersection of parent with the 57 existing PSF-convolved joint-injection pairs",
        "n_parent_pairs": len(parent),
        "n_psf_injection_pairs": len(psf_keys),
        "n_other_parent_pairs": len(parent_by_key) - len(psf_keys),
        "selected_is_subset_of_parent": True,
        "by_filter": {},
        "input_sha256": {
            str(args.injection_frames.resolve()): sha256(args.injection_frames),
            str(args.icl_support.resolve()): sha256(args.icl_support),
            str(args.satellite_profiles.resolve()): sha256(args.satellite_profiles),
        },
        "code_sha256": {str(Path(__file__).resolve()): sha256(Path(__file__).resolve())},
        "interpretation_limit": "Descriptive covariate and support comparison; it does not correct selection or establish population representativeness.",
    }
    table_rows = []
    for filt, band_rows in grouped.items():
        band_keys = {pair_key(row, "id") for row in band_rows}
        selected = [parent_by_key[key] for key in sorted(psf_keys & band_keys)]
        remainder = [row for row in band_rows if pair_key(row, "id") not in psf_keys]
        if not selected or not remainder:
            raise ValueError(f"empty selected/remainder group for filter {filt}")
        selected_metrics, remainder_metrics = metrics(selected), metrics(remainder)
        comparisons = {}
        for covariate in ("z", "richness", "coverage_200_300"):
            a = np.asarray([float(row[covariate]) for row in selected])
            b = np.asarray([float(row[covariate]) for row in remainder])
            pooled = float(np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2))
            comparisons[covariate] = {
                "selected_minus_remainder_mean": float(a.mean() - b.mean()),
                "standardized_mean_difference": float((a.mean() - b.mean()) / pooled) if pooled > 0 else None,
            }
        report["by_filter"][filt] = {
            "selected": selected_metrics,
            "remainder": remainder_metrics,
            "selected_minus_remainder": comparisons,
            "satellite_50_75_support_at_60pct": {
                "selected": support_count(selected, sat60),
                "parent": support_count(band_rows, sat60),
            },
            "satellite_50_75_support_at_40pct": {
                "selected": support_count(selected, sat40),
                "parent": support_count(band_rows, sat40),
            },
        }
        for group_name, group, group_metrics in (
            ("psf_selected", selected, selected_metrics),
            ("other_parent", remainder, remainder_metrics),
        ):
            sat_60 = support_count(group, sat60)
            sat_40 = support_count(group, sat40)
            for row in group:
                key = pair_key(row, "id")
                table_rows.append({
                    "cluster_id": key[0], "filter": key[1], "group": group_name,
                    "z": row["z"], "richness": row["richness"],
                    "icl_coverage_200_300": row["coverage_200_300"],
                    "satellite_50_75_supported_60pct": int(key in sat60),
                    "satellite_50_75_supported_40pct": int(key in sat40),
                })

    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    (output / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "pair_selection.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, list(table_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(table_rows)

    lines = [
        "# Representativeness of the PSF-supported joint-injection subset", "",
        "This offline Astrafek audit compares the 57 pairs with existing position-dependent PSF injection results against their 224-pair ICL-supported parent. It reads saved summary/profile tables only, makes no HSC request, reads no image pixels, and retains `science_measurement=false`.", "",
        "| Filter | PSF subset / parent | Δ mean redshift (SMD) | Δ mean richness (SMD) | Δ mean ICL coverage (SMD) | Satellite support 60%, subset / parent |", "|---|---:|---:|---:|---:|---:|",
    ]
    for filt in ("i", "i2", "all"):
        item = report["by_filter"][filt]
        comps = item["selected_minus_remainder"]
        n = item["selected"]["n_pairs"]
        m = item["selected"]["n_pairs"] + item["remainder"]["n_pairs"]
        sat = item["satellite_50_75_support_at_60pct"]
        lines.append(
            f"| {filt} | {n}/{m} | {comps['z']['selected_minus_remainder_mean']:+.4f} ({comps['z']['standardized_mean_difference']:+.3f}) | "
            f"{comps['richness']['selected_minus_remainder_mean']:+.3f} ({comps['richness']['standardized_mean_difference']:+.3f}) | "
            f"{comps['coverage_200_300']['selected_minus_remainder_mean']:+.4f} ({comps['coverage_200_300']['standardized_mean_difference']:+.3f}) | "
            f"{sat['selected']['n_supported']}/{sat['selected']['n_total']} vs {sat['parent']['n_supported']}/{sat['parent']['n_total']} |"
        )
    lines += [
        "", "The PSF subset is selected toward higher ICL-annulus coverage and has higher outer-satellite support than the parent. Its injection transfer should not be generalized to all 224 supported pairs. This comparison is descriptive and does not provide an inverse-probability correction, injection completeness, or population inference.", "",
        f"Summary SHA-256: `{sha256(output / 'summary.json')}`.",
    ]
    (output / "assessment.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "output_dir": str(output),
        "n_parent_pairs": report["n_parent_pairs"],
        "n_psf_injection_pairs": report["n_psf_injection_pairs"],
        "summary_sha256": sha256(output / "summary.json"),
        "pair_table_sha256": sha256(output / "pair_selection.tsv"),
        "assessment_sha256": sha256(output / "assessment.md"),
    }, indent=2))


if __name__ == "__main__":
    main()
