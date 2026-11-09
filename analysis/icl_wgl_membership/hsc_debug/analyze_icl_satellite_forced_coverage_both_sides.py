#!/usr/bin/env python3
"""Paired-block comparison of detected-only and 20/20%-coverage satellite stacks.

Uses only already-downloaded HSC coadds, masks, CAMIRA member positions, and
saved profiles. It decomposes the change in the ICL--satellite partial
correlation into (i) adding newly supported cluster/filter pairs and (ii)
adding forced CAMIRA positions within pairs that already had detected-satellite
profiles. The coverage floor is post-hoc; every result remains exploratory.
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
FORCED = DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_both_2026-11-06_01"
ICL = DATA / "icl_annulus_support_2026-10-18_01" / "annulus_pair_profiles.tsv"
TARGETS = DATA / "download_01" / "download_manifest.json"
CONTROLS = DATA / "control_download_01" / "download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_forced_coverage_both_sides_2026-11-06_03"
REPLICATES = 10_000
SEED = 20261106


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def partial_r(rows: list[dict], y_key: str) -> float:
    if len(rows) < 8:
        return math.nan
    x = np.asarray([float(r["icl_delta_200_300"]) for r in rows])
    y = np.asarray([float(r[y_key]) for r in rows])
    z = np.asarray([float(r["redshift"]) for r in rows])
    richness = np.log(np.asarray([float(r["richness"]) for r in rows]))
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if np.linalg.matrix_rank(design) < 3:
        return math.nan
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    norm = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(np.dot(rx, ry) / norm) if norm > 0 else math.nan


def bootstrap_comparison(rows: list[dict], band: str, block_labels: np.ndarray,
                         replicas: int, seed: int) -> dict:
    selected = [r for r in rows if r["filter"] == band]
    if not selected:
        return {"status": "no_support", "filter": band}
    group_by_key = {(r["cluster_id"], r["filter"]): str(block_labels[i])
                    for i, r in enumerate(rows) if r["filter"] == band}
    groups = np.asarray(sorted(set(group_by_key.values())))
    by_group: dict[int, list[dict]] = {}
    for row in selected:
        by_group.setdefault(group_by_key[(row["cluster_id"], band)], []).append(row)
    rng = np.random.default_rng(seed)
    changes_expansion, changes_measurement, changes_total = [], [], []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        boot = [row for group in chosen for row in by_group[str(group)]]
        full = [r for r in boot if r["satellite_delta_combined_50_75"] is not None]
        common = [r for r in full if r["satellite_delta_detected_50_75"] is not None]
        detected = [dict(r, satellite_delta_combined_50_75=r["satellite_delta_detected_50_75"])
                    for r in common]
        r_full = partial_r(full, "satellite_delta_combined_50_75")
        r_common = partial_r(common, "satellite_delta_combined_50_75")
        r_detected = partial_r(detected, "satellite_delta_combined_50_75")
        if np.isfinite([r_full, r_common, r_detected]).all():
            changes_expansion.append(r_full - r_common)
            changes_measurement.append(r_common - r_detected)
            changes_total.append(r_full - r_detected)
    if len(changes_total) < 0.9 * replicas:
        raise ValueError(f"{band}: only {len(changes_total)}/{replicas} valid paired draws")

    def summarize(values: list[float]) -> dict:
        a = np.asarray(values)
        return {"median": float(np.median(a)),
                "ci95": np.quantile(a, [0.025, 0.975]).tolist(),
                "fraction_positive": float(np.mean(a > 0))}

    full = [r for r in selected if r["satellite_delta_combined_50_75"] is not None]
    common = [r for r in full if r["satellite_delta_detected_50_75"] is not None]
    detected = [dict(r, satellite_delta_combined_50_75=r["satellite_delta_detected_50_75"])
                for r in common]
    return {
        "status": "complete_exploratory", "filter": band,
        "n_20_20_combined_pairs": len(full), "n_detected_common_pairs": len(common),
        "n_connected_5deg_blocks_full": len(groups),
        "partial_r_20_20_full": partial_r(full, "satellite_delta_combined_50_75"),
        "partial_r_20_20_common": partial_r(common, "satellite_delta_combined_50_75"),
        "partial_r_detected_only_common": partial_r(detected, "satellite_delta_combined_50_75"),
        "paired_change_sample_expansion_full_minus_common": summarize(changes_expansion),
        "paired_change_added_forced_measurements_common_minus_detected": summarize(changes_measurement),
        "paired_change_total_full_minus_detected": summarize(changes_total),
        "replicas_requested": replicas, "replicas_valid": len(changes_total), "seed": seed,
        "bootstrap_note": "shared connected 5-degree target/control blocks defined on the 20/20-percent sample; exploratory paired uncertainty",
    }


def run(output: Path, replicas: int = REPLICATES) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    forced_summary_path = FORCED / "summary.json"
    cluster_path = FORCED / "cluster_profiles.tsv"
    association_path = FORCED / "icl_satellite_association.tsv"
    for path in (forced_summary_path, cluster_path, association_path, ICL, TARGETS, CONTROLS):
        if not path.is_file():
            raise FileNotFoundError(path)
    forced_summary = read_json(forced_summary_path)
    if (forced_summary.get("status") != "complete_exploratory"
            or forced_summary.get("science_measurement") is not False
            or forced_summary.get("method", {}).get("target_coverage_floor") != 0.2
            or forced_summary.get("method", {}).get("control_coverage_floor") != 0.2
            or forced_summary.get("method", {}).get("mask_other_unmatched_camira") is not False):
        raise ValueError("forced-profile summary does not match the pinned 20/20-percent scope")
    for name, expected in forced_summary.get("output_sha256", {}).items():
        path = FORCED / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"forced-profile output hash mismatch: {name}")

    detected = {(str(r["cluster_id"]), r["filter"]): r
                for r in read_tsv(FORCED / "icl_satellite_association.tsv")}
    icl = {(str(r["id"]), r["filter"]): r for r in read_tsv(ICL)
           if r["radial_bin_kpc"] == "200-300"}
    rows = []
    for key, ref in sorted(icl.items(), key=lambda item: (item[0][1], int(item[0][0]))):
        r = detected.get(key, {})
        detected_value = None if not r.get("satellite_delta_detected_50_75") else float(r["satellite_delta_detected_50_75"])
        combined_value = None if not r.get("satellite_delta_combined_50_75") else float(r["satellite_delta_combined_50_75"])
        forced_value = None if not r.get("satellite_delta_forced_50_75") else float(r["satellite_delta_forced_50_75"])
        rows.append({
            "cluster_id": key[0], "filter": key[1], "redshift": float(ref["z"]),
            "richness": float(ref["richness"]),
            "icl_delta_200_300": float(ref["delta_sb_njy_arcsec2"]),
            "satellite_delta_detected_50_75": detected_value,
            "satellite_delta_forced_50_75": forced_value,
            "satellite_delta_combined_50_75": combined_value,
            "n_detected_positions": int(r.get("n_detected_positions", 0)),
            "n_forced_positions": int(r.get("n_forced_positions", 0)),
            "selected_20_20": combined_value is not None,
        })
    if len({(r["cluster_id"], r["filter"]) for r in rows}) != len(rows):
        raise ValueError("duplicate ICL/satellite cluster-filter pair")

    target_manifest, control_manifest = read_json(TARGETS), read_json(CONTROLS)
    target_by_id = {str(r["id"]): r for r in target_manifest["clusters"]}
    control_by_id = {str(r["target_cluster_id"]): r for r in control_manifest["clusters"]}
    labeled_rows = []
    for row in sorted(rows, key=lambda r: (r["filter"], int(r["cluster_id"]))):
        target, control = target_by_id.get(row["cluster_id"]), control_by_id.get(row["cluster_id"])
        if target is None or control is None:
            raise ValueError(f"missing paired target/control coordinates for {row['cluster_id']}")
        labeled_rows.append(row)
    labels = np.full(len(labeled_rows), "", dtype=object)
    for band in ("i", "i2"):
        idx = np.asarray([i for i, r in enumerate(labeled_rows) if r["filter"] == band], dtype=int)
        coords_t = [{"ra_deg": float(target_by_id[labeled_rows[i]["cluster_id"]]["ra_deg"]),
                     "dec_deg": float(target_by_id[labeled_rows[i]["cluster_id"]]["dec_deg"])} for i in idx]
        coords_c = [{"ra_deg": float(control_by_id[labeled_rows[i]["cluster_id"]]["ra_deg"]),
                     "dec_deg": float(control_by_id[labeled_rows[i]["cluster_id"]]["dec_deg"])} for i in idx]
        _, band_labels, *_ = pair_block_labels(coords_t, coords_c,
                                                np.ones(len(idx), dtype=bool))
        labels[idx] = [str(value) for value in band_labels]
    comparisons = {}
    for i, band in enumerate(("i", "i2")):
        selected_idx = [j for j, row in enumerate(labeled_rows)
                        if row["filter"] == band and row["selected_20_20"]]
        selected_rows = [labeled_rows[j] for j in selected_idx]
        selected_labels = labels[selected_idx]
        comparisons[band] = bootstrap_comparison(selected_rows, band, selected_labels,
                                                  replicas, SEED + i)

    output.mkdir(parents=True, exist_ok=False)
    fields = ["cluster_id", "filter", "redshift", "richness", "icl_delta_200_300",
              "satellite_delta_detected_50_75", "satellite_delta_forced_50_75",
              "satellite_delta_combined_50_75", "n_detected_positions",
              "n_forced_positions", "selected_20_20"]
    write_tsv(output / "paired_cluster_comparison.tsv", labeled_rows, fields)
    selection = {}
    for band in ("i", "i2"):
        group = [r for r in labeled_rows if r["filter"] == band]
        selected = [r for r in group if r["selected_20_20"]]
        selection[band] = {
            "icl_supported_pairs": len(group), "20_20_selected_pairs": len(selected),
            "excluded_pairs": len(group) - len(selected),
            "selected_clusters": len({r["cluster_id"] for r in selected}),
            "excluded_clusters": len({r["cluster_id"] for r in group if not r["selected_20_20"]}),
            "selected_target_positions": sum(r["n_detected_positions"] + r["n_forced_positions"] for r in selected),
            "selected_median_redshift": float(np.median([r["redshift"] for r in selected])) if selected else None,
            "excluded_median_redshift": float(np.median([r["redshift"] for r in group if not r["selected_20_20"]])) if len(selected) < len(group) else None,
            "selected_median_richness": float(np.median([r["richness"] for r in selected])) if selected else None,
            "excluded_median_richness": float(np.median([r["richness"] for r in group if not r["selected_20_20"]])) if len(selected) < len(group) else None,
            "selected_median_icl_delta_200_300": float(np.median([r["icl_delta_200_300"] for r in selected])) if selected else None,
            "excluded_median_icl_delta_200_300": float(np.median([r["icl_delta_200_300"] for r in group if not r["selected_20_20"]])) if len(selected) < len(group) else None,
        }
    result = {
        "schema": "astrafek-icl-satellite-forced-coverage-paired-comparison-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "scope": "existing HSC cutouts only; target and rotated-control annular coverage floors both 20 percent",
        "method": "detected-only versus combined detected+forced cluster profiles; decompose sample expansion and forced-position measurement changes on a common sample using paired connected-5-degree blocks",
        "selection_support": selection,
        "association_by_filter": comparisons,
        "limitations": [
            "The 20/20-percent rule is post-hoc and incomplete-sample support remains selected; these results do not replace the 60-percent primary sensitivity.",
            "Forced CAMIRA centers without SEP matches are not confirmed isolated galaxies; mean light may include blends or projected contaminants.",
            "The statistic is a selected-position stack, not a calibrated satellite stellar-envelope measurement or ICL detection.",
            "PSF wings, sky transfer, source-Poisson variance, full coadd covariance, and validated survey-region uncertainty remain unresolved.",
            "No multiplicity adjustment across coverage choices or model specifications is included.",
        ],
        "inputs_sha256": {str(p): sha256(p) for p in (
            forced_summary_path, cluster_path, association_path, ICL, TARGETS, CONTROLS)},
        "replicas_requested": replicas,
    }
    result["output_sha256"] = {
        "paired_cluster_comparison.tsv": sha256(output / "paired_cluster_comparison.tsv")
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                                       allow_nan=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    if args.replicates < 1000:
        parser.error("at least 1,000 paired block-bootstrap replicates are required")
    result = run(args.output, args.replicates)
    print(json.dumps({"status": result["status"], "selection_support": result["selection_support"],
                      "association_by_filter": result["association_by_filter"],
                      "science_measurement": result["science_measurement"]}, indent=2))


if __name__ == "__main__":
    main()
