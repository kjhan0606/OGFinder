#!/usr/bin/env python3
"""Audit fixed cluster/filter support for every ICL x satellite radial pair.

This is a measurement-support audit, not an association test. It uses only
already downloaded HSC profile tables and reports the exact common rows needed
before any radial cross-correlation analysis can be interpreted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ICL = DATA / "icl_annulus_support_2026-10-18_01"
DEFAULT_SATELLITE = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
ICL_BINS = ("50-100", "100-150", "150-200", "200-300", "300-400", "400-500")
SATELLITE_BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_summary(summary_path: Path, expected_schema: str,
                     table_path: Path, *, require_exploratory: bool) -> dict:
    summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    if summary.get("schema") != expected_schema or summary.get("status") != "complete":
        raise ValueError(f"unexpected or incomplete source summary: {summary_path}")
    if require_exploratory and summary.get("science_measurement") is not False:
        raise ValueError("satellite source unexpectedly has a science-measurement flag")
    if not require_exploratory and summary.get("science_gates", {}).get("science_measurement") is not False:
        raise ValueError("ICL source unexpectedly passes the science-measurement gate")
    outputs = summary.get("output_sha256", {})
    expected = outputs.get(str(table_path.resolve()))
    if expected is None or sha256(table_path) != expected:
        raise ValueError(f"source table hash is missing or mismatched: {table_path}")
    return summary


def support_matrix(icl_rows: list[dict[str, str]], satellite_rows: list[dict[str, str]]) -> list[dict]:
    """Return intersection counts for each filter and pair of radial bins."""
    for name, rows, id_field, radial_field in (
        ("ICL", icl_rows, "id", "radial_bin_kpc"),
        ("satellite", satellite_rows, "cluster_id", "radial_bin_kpc"),
    ):
        keys = [(str(row[id_field]), row["filter"], row[radial_field]) for row in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f"duplicate {name} cluster/filter/radial row")
    icl_sets = {
        (band, radial): {
            str(row["id"]) for row in icl_rows
            if row["filter"] == band and row["radial_bin_kpc"] == radial
        }
        for band in FILTERS for radial in ICL_BINS
    }
    satellite_sets = {
        (band, radial): {
            str(row["cluster_id"]) for row in satellite_rows
            if row.get("cohort") == "full_sample" and row["filter"] == band
            and row["radial_bin_kpc"] == radial
        }
        for band in FILTERS for radial in SATELLITE_BINS
    }
    output = []
    for band in FILTERS:
        for icl_bin in ICL_BINS:
            icl_ids = icl_sets[(band, icl_bin)]
            for satellite_bin in SATELLITE_BINS:
                satellite_ids = satellite_sets[(band, satellite_bin)]
                common = icl_ids & satellite_ids
                output.append({
                    "filter": band,
                    "icl_radial_bin_kpc": icl_bin,
                    "satellite_radial_bin_kpc": satellite_bin,
                    "n_icl_supported": len(icl_ids),
                    "n_satellite_supported": len(satellite_ids),
                    "n_common": len(common),
                    "common_fraction_of_icl": len(common) / len(icl_ids) if icl_ids else None,
                    "common_fraction_of_satellite": len(common) / len(satellite_ids) if satellite_ids else None,
                    "cluster_ids_sha256": hashlib.sha256(
                        "\n".join(sorted(common)).encode("utf-8")
                    ).hexdigest(),
                })
    return output


def analyze(icl_dir: Path, satellite_dir: Path, output_dir: Path) -> dict:
    icl_dir, satellite_dir, output_dir = map(Path, (icl_dir, satellite_dir, output_dir))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output: {output_dir}")
    icl_summary_path, satellite_summary_path = icl_dir / "summary.json", satellite_dir / "summary.json"
    icl_table, satellite_table = icl_dir / "annulus_pair_profiles.tsv", satellite_dir / "cluster_profiles.tsv"
    icl_summary = validate_summary(icl_summary_path, "hsc-icl-per-annulus-support-v1", icl_table,
                                   require_exploratory=False)
    satellite_summary = validate_summary(satellite_summary_path,
                                         "hsc-satellite-outer-profile-per-annulus-v1",
                                         satellite_table, require_exploratory=True)
    if satellite_summary.get("n_frames_completed") != 800:
        raise ValueError("satellite source is not the fixed 800-frame full sample")
    icl_rows = [row for row in read_tsv(icl_table) if row.get("radial_bin_kpc") in ICL_BINS]
    satellite_rows = read_tsv(satellite_table)
    matrix = support_matrix(icl_rows, satellite_rows)
    if len(matrix) != len(FILTERS) * len(ICL_BINS) * len(SATELLITE_BINS):
        raise ValueError("support matrix has an unexpected number of cells")
    if any(row["n_common"] == 0 for row in matrix):
        raise ValueError("one or more expected radial combinations have no common support")

    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "radial_support_matrix.tsv"
    fields = list(matrix[0])
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matrix)
    table_sha = sha256(table_path)
    maxima = {}
    for band in FILTERS:
        cells = [row for row in matrix if row["filter"] == band]
        maximum = max(cells, key=lambda row: row["n_common"])
        outer = next(row for row in cells if row["icl_radial_bin_kpc"] == "200-300"
                     and row["satellite_radial_bin_kpc"] == "50-75")
        maxima[band] = {"max_common_support": maximum,
                        "outer_candidate_200_300_x_50_75": outer}
    result = {
        "schema": "astrafek-icl-satellite-radial-support-matrix-v1",
        "status": "complete_support_audit",
        "science_measurement": False,
        "inference_performed": False,
        "method": {
            "unit": "unique cluster/filter pair",
            "join_key": ["cluster ID", "filter"],
            "scope": "six ICL annuli crossed with seven satellite annuli, i and i2 filters",
            "support_rule": "exact intersection of the existing independently selected ICL and satellite profile tables; no reselection or weighting",
        },
        "source": {
            "icl_summary": {"path": str(icl_summary_path.resolve()), "sha256": sha256(icl_summary_path)},
            "icl_table": {"path": str(icl_table.resolve()), "sha256": sha256(icl_table)},
            "satellite_summary": {"path": str(satellite_summary_path.resolve()), "sha256": sha256(satellite_summary_path)},
            "satellite_table": {"path": str(satellite_table.resolve()), "sha256": sha256(satellite_table)},
        },
        "n_cells": len(matrix),
        "by_filter": maxima,
        "limitations": [
            "support counts do not establish representativeness or missing-at-random selection",
            "no association coefficients, significance tests, or selection correction are calculated",
            "all source profiles remain exploratory; image-level PSF wings, sky transfer, and full covariance remain unresolved",
        ],
        "outputs_sha256": {str(table_path.resolve()): table_sha},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl-dir", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.icl_dir, args.satellite_dir, args.output_dir)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_cells": result["n_cells"], "by_filter": result["by_filter"],
                      "output": str(args.output_dir.resolve())}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
