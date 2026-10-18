#!/usr/bin/env python3
"""Re-evaluate the per-annulus ICL--satellite relation with spatial blocks.

This joins already exported cluster/filter rows to the exact ICL target/control
coordinates, then compares ordinary pair bootstrap intervals with a connected
5-degree target/control block bootstrap. It is exploratory and never promotes
the science gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_ICL = ROOT / "icl_annulus_support_2026-10-18_01"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty result table")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def partial_pearson(rows: list[dict], y_key: str) -> float:
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[y_key]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    if len(rows) < 5 or not all(np.all(np.isfinite(a)) for a in (x, y, z, richness)):
        raise ValueError("partial correlation requires >=5 finite cluster/filter rows")
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("redshift/richness nuisance design is rank deficient")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if denominator <= 0:
        raise ValueError("partial correlation is undefined for a constant residual")
    return float(rx @ ry / denominator)


def spearman_rho(rows: list[dict], y_key: str) -> float:
    # Spearman rho is Pearson correlation of average ranks; avoid an optional
    # SciPy dependency in this small, reproducible analysis utility.
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[y_key]) for row in rows])
    if len(rows) < 3 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Spearman correlation requires >=3 finite rows")

    def ranks(values):
        order = np.argsort(values, kind="mergesort")
        sorted_values = values[order]
        output = np.empty(len(values), dtype=float)
        start = 0
        while start < len(values):
            stop = start + 1
            while stop < len(values) and sorted_values[stop] == sorted_values[start]:
                stop += 1
            output[order[start:stop]] = (start + stop - 1) / 2.0
            start = stop
        return output

    rx, ry = ranks(x), ranks(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if denominator <= 0:
        raise ValueError("Spearman correlation is undefined for a constant vector")
    return float(rx @ ry / denominator)


def coordinate_labels(rows: list[dict]) -> np.ndarray:
    targets = [{"ra_deg": float(row["target_ra_deg"]),
                "dec_deg": float(row["target_dec_deg"])} for row in rows]
    controls = [{"ra_deg": float(row["control_ra_deg"]),
                 "dec_deg": float(row["control_dec_deg"])} for row in rows]
    selected, labels, *_ = pair_block_labels(targets, controls,
                                               np.ones(len(rows), dtype=bool))
    if len(selected) != len(rows) or len(labels) != len(rows):
        raise ValueError("spatial-block helper changed the selected row count")
    return labels


def bootstrap(rows: list[dict], y_key: str, labels: np.ndarray | None,
              replicas: int, seed: int) -> dict:
    n = len(rows)
    if n < 5:
        return {"status": "insufficient_pairs", "n_pairs": n,
                "replicas_requested": replicas, "replicas_valid": 0, "ci95": None}
    if labels is None:
        groups = np.arange(n)
        method = "ordinary_pair_bootstrap"
    else:
        labels = np.asarray(labels, dtype=str)
        if len(labels) != n:
            raise ValueError("one spatial block label is required per input row")
        groups = np.unique(labels)
        method = "connected_target_control_5deg_block_bootstrap"
        if len(groups) < 6:
            return {"status": "insufficient_spatial_blocks", "n_pairs": n,
                    "n_blocks": int(len(groups)), "replicas_requested": replicas,
                    "replicas_valid": 0, "ci95": None}
    point_partial = partial_pearson(rows, y_key)
    point_spearman = spearman_rho(rows, y_key)
    indices_by_group = ({group: np.flatnonzero(labels == group) for group in groups}
                        if labels is not None else None)
    rng = np.random.default_rng(seed)
    partial_draws, spearman_draws = [], []
    invalid = 0
    for _ in range(replicas):
        selected = rng.choice(groups, size=len(groups), replace=True)
        if indices_by_group is None:
            indices = selected.astype(int)
        else:
            indices = np.concatenate([indices_by_group[group] for group in selected])
        sample = [rows[index] for index in indices]
        try:
            partial_draws.append(partial_pearson(sample, y_key))
            spearman_draws.append(spearman_rho(sample, y_key))
        except ValueError:
            invalid += 1
    if not partial_draws:
        return {"status": "no_valid_replicates", "method": method, "n_pairs": n,
                "n_blocks": int(len(groups)) if labels is not None else None,
                "replicas_requested": replicas, "replicas_valid": 0,
                "replicas_invalid": invalid, "ci95_partial_r": None,
                "ci95_spearman_rho": None}
    quantile = [2.5, 97.5]
    block_sizes = [len(indices) for indices in indices_by_group.values()] if labels is not None else []
    return {
        "status": "complete_exploratory",
        "method": method,
        "n_pairs": n,
        "n_blocks": int(len(groups)) if labels is not None else None,
        "largest_block_pairs": max(block_sizes) if block_sizes else None,
        "partial_r_point": point_partial,
        "ci95_partial_r": np.percentile(partial_draws, quantile).tolist(),
        "spearman_rho_point": point_spearman,
        "ci95_spearman_rho": np.percentile(spearman_draws, quantile).tolist(),
        "replicas_requested": replicas,
        "replicas_valid": len(partial_draws),
        "replicas_invalid": invalid,
        "seed": seed,
    }


def load_verified_rows(association_dir: Path, icl_dir: Path,
                       extraction_summary_path: Path) -> tuple[list[dict], dict]:
    association_dir = Path(association_dir).resolve()
    icl_dir = Path(icl_dir).resolve()
    extraction_summary_path = Path(extraction_summary_path).resolve()
    association_summary_path = association_dir / "summary.json"
    association_table_path = association_dir / "matched_cluster_values.tsv"
    icl_summary_path = icl_dir / "summary.json"
    annulus_table_path = icl_dir / "annulus_pair_profiles.tsv"
    for path in (association_summary_path, association_table_path, icl_summary_path,
                 annulus_table_path, extraction_summary_path):
        if not path.is_file():
            raise ValueError(f"missing required input: {path}")
    assoc_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    icl_summary = json.loads(icl_summary_path.read_text(encoding="utf-8"))
    extraction = json.loads(extraction_summary_path.read_text(encoding="utf-8"))
    if (assoc_summary.get("schema") != "hsc-satellite-outer-icl-correlation-v1"
            or assoc_summary.get("status") != "complete"
            or assoc_summary.get("science_measurement") is not False):
        raise ValueError("unexpected association input or science gate")
    if assoc_summary.get("outputs_sha256", {}).get(association_table_path.name) != sha256(association_table_path):
        raise ValueError("association table hash differs from its summary")
    if (icl_summary.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_summary.get("status") != "complete"
            or icl_summary.get("science_gates", {}).get("science_measurement") is not False):
        raise ValueError("unexpected ICL support input or science gate")
    recorded_icl_hash = icl_summary.get("output_sha256", {}).get(str(annulus_table_path))
    if recorded_icl_hash != sha256(annulus_table_path):
        raise ValueError("ICL annulus coordinate table hash differs from its summary")
    if (extraction.get("status") != "complete"
            or str(extraction.get("science_measurement", "")).lower() != "false"
            or len(extraction.get("records", [])) != int(extraction.get("n_measured", -1))):
        raise ValueError("control-coordinate extraction summary is incomplete or science-promoted")

    control_records = {}
    for record in extraction["records"]:
        control_id = str(record["control_id"])
        if control_id in control_records:
            raise ValueError(f"duplicate control ID in extraction summary: {control_id}")
        if record.get("frame_role") != "matched_control":
            raise ValueError(f"non-control frame in matched-control extraction: {control_id}")
        ra, dec = float(record["ra_deg"]), float(record["dec_deg"])
        if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
            raise ValueError(f"invalid control coordinates for {control_id}")
        control_records[control_id] = record

    annulus_rows = read_tsv(annulus_table_path)
    icl_by_key = {}
    for row in annulus_rows:
        if row["radial_bin_kpc"] != "200-300":
            continue
        key = (str(row["id"]), row["filter"])
        if key in icl_by_key:
            raise ValueError(f"duplicate ICL 200-300 kpc pair: {key}")
        icl_by_key[key] = row
    association_rows = read_tsv(association_table_path)
    joined, seen = [], set()
    for row in association_rows:
        key = (str(row["cluster_id"]), row["filter"])
        if key in seen:
            raise ValueError(f"duplicate association row: {key}")
        seen.add(key)
        icl = icl_by_key.get(key)
        if icl is None:
            raise ValueError(f"association row lacks ICL 200-300 kpc data: {key}")
        control_id = str(icl["control_id"])
        control = control_records.get(control_id)
        if control is None:
            raise ValueError(f"control coordinate missing for {control_id}")
        if control.get("header_filter") != row["filter"]:
            raise ValueError(f"control filter mismatch for {key}")
        if not math.isclose(float(row["icl_delta_200_300"]),
                            float(icl["delta_sb_njy_arcsec2"]), rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"association ICL value does not reproduce from annulus table: {key}")
        if not math.isclose(float(row["redshift"]), float(icl["z"]), rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"redshift mismatch between association and ICL table: {key}")
        if not math.isclose(float(row["richness"]), float(icl["richness"]), rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"richness mismatch between association and ICL table: {key}")
        target_ra, target_dec = float(icl["ra_deg"]), float(icl["dec_deg"])
        if not (math.isfinite(target_ra) and 0 <= target_ra < 360
                and math.isfinite(target_dec) and -90 <= target_dec <= 90):
            raise ValueError(f"invalid target coordinates for {key}")
        joined.append({
            **row,
            "target_ra_deg": target_ra,
            "target_dec_deg": target_dec,
            "control_id": control_id,
            "control_ra_deg": float(control["ra_deg"]),
            "control_dec_deg": float(control["dec_deg"]),
        })
    if len(joined) != 94:
        raise ValueError(f"fixed per-annulus association sample changed: expected 94, found {len(joined)}")
    return joined, {
        "association_summary_sha256": sha256(association_summary_path),
        "association_table_sha256": sha256(association_table_path),
        "icl_support_summary_sha256": sha256(icl_summary_path),
        "icl_annulus_table_sha256": sha256(annulus_table_path),
        "control_extraction_summary_sha256": sha256(extraction_summary_path),
        "association_rows": len(association_rows),
        "joined_rows": len(joined),
        "coordinate_join": "association cluster/filter -> ICL annulus target and control_id -> validated matched-control extraction RA/Dec",
    }


def analyze(association_dir: Path, icl_dir: Path, extraction_summary: Path,
            output_dir: Path, replicas: int = 5000, seed: int = 20261020) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 bootstrap replicas are required")
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_dir}")
    rows, provenance = load_verified_rows(association_dir, icl_dir, extraction_summary)
    results, flat_results = [], []
    for band_index, band in enumerate(("i2", "i")):
        band_rows = [row for row in rows if row["filter"] == band]
        for selection, selected_rows in (
            ("all_supported_pairs", band_rows),
            ("center_eligible", [row for row in band_rows
                                 if str(row["center_eligible"]).lower() == "true"]),
        ):
            if len(selected_rows) < 5:
                results.append({"filter": band, "selection": selection,
                                "status": "insufficient_pairs", "n_pairs": len(selected_rows)})
                continue
            labels = coordinate_labels(selected_rows)
            spatial = bootstrap(selected_rows, "satellite_delta_50_75", labels,
                                replicas, seed + 2 * band_index + int(selection == "center_eligible"))
            ordinary = bootstrap(selected_rows, "satellite_delta_50_75", None,
                                 replicas, seed + 2 * band_index + int(selection == "center_eligible") + 101)
            unique, counts = np.unique(labels, return_counts=True)
            result = {
                "filter": band,
                "selection": selection,
                "status": "complete_exploratory",
                "n_pairs": len(selected_rows),
                "n_connected_5deg_blocks": int(len(unique)),
                "largest_connected_block_pairs": int(counts.max()),
                "partial_pearson_r_control_redshift_richness": spatial.get("partial_r_point"),
                "partial_spatial_block_bootstrap": spatial,
                "partial_ordinary_pair_bootstrap": ordinary,
                "spearman_rho": spatial.get("spearman_rho_point"),
                "spearman_spatial_block_ci95": spatial.get("ci95_spearman_rho"),
                "interpretation": "spatial block interval is primary; ordinary pair bootstrap is shown only to expose the iid sensitivity",
            }
            results.append(result)
            flat_results.append({
                "filter": band,
                "selection": selection,
                "n_pairs": len(selected_rows),
                "n_connected_5deg_blocks": len(unique),
                "largest_block_pairs": int(counts.max()),
                "partial_r": spatial.get("partial_r_point"),
                "spatial_ci95_low": (spatial.get("ci95_partial_r") or [None, None])[0],
                "spatial_ci95_high": (spatial.get("ci95_partial_r") or [None, None])[1],
                "ordinary_ci95_low": (ordinary.get("ci95_partial_r") or [None, None])[0],
                "ordinary_ci95_high": (ordinary.get("ci95_partial_r") or [None, None])[1],
                "spearman_rho": spatial.get("spearman_rho_point"),
                "spearman_spatial_ci95_low": (spatial.get("ci95_spearman_rho") or [None, None])[0],
                "spearman_spatial_ci95_high": (spatial.get("ci95_spearman_rho") or [None, None])[1],
                "valid_spatial_replicas": spatial.get("replicas_valid", 0),
                "science_measurement": False,
            })
    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "spatial_bootstrap_results.tsv"
    write_tsv(table_path, flat_results)
    summary = {
        "schema": "hsc-satellite-icl-spatial-block-bootstrap-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "method": {
            "fixed_input_sample": "all 94 per-annulus matched cluster/filter rows; no re-selection",
            "primary_filter": "i2",
            "primary_uncertainty": "connected target/control 5-degree spatial-block bootstrap",
            "comparison": "ordinary row bootstrap shown as a sensitivity comparison",
            "association": "partial Pearson residualized on redshift and richness; Spearman is an unadjusted descriptive companion",
            "replicas": replicas,
            "seed": seed,
        },
        "provenance": provenance,
        "results": results,
        "limitations": [
            "All results remain conditional on the selected per-annulus support sample, not the full 800-cluster population.",
            "This does not correct mask, PSF, sky transfer, extinction, member completeness, shared pixels, or image-level covariance.",
            "The association was exploratory and post hoc; intervals are not confirmatory significance tests.",
        ],
        "output_table": str(table_path),
        "output_table_sha256": sha256(table_path),
        "program_sha256": sha256(Path(__file__).resolve()),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    try:
        result = analyze(args.association, args.icl, args.extraction_summary,
                         args.output, args.replicas, args.seed)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        parser.error(str(error))
    print(json.dumps({"status": result["status"],
                      "science_measurement": result["science_measurement"],
                      "output_table": result["output_table"],
                      "results": [{"filter": row["filter"], "selection": row["selection"],
                                   "n_pairs": row["n_pairs"],
                                   "n_blocks": row.get("n_connected_5deg_blocks"),
                                   "partial_r": row.get("partial_pearson_r_control_redshift_richness"),
                                   "spatial_ci95": row.get("partial_spatial_block_bootstrap", {}).get("ci95_partial_r"),
                                   "ordinary_ci95": row.get("partial_ordinary_pair_bootstrap", {}).get("ci95_partial_r")}
                                  for row in result["results"]]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
