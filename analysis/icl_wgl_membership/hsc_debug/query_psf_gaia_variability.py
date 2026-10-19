#!/usr/bin/env python3
"""Query Gaia DR3 variability summaries for already matched PSF candidates."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

ENDPOINT = "https://gea.esac.esa.int/tap-server/tap/sync"
SUMMARY_FIELDS = [
    "source_id", "num_selected_g_fov", "time_duration_g_fov",
    "mean_mag_g_fov", "median_mag_g_fov", "std_dev_mag_g_fov",
    "mad_mag_g_fov", "range_mag_g_fov", "trimmed_range_mag_g_fov",
]
CLASSIFIER_FIELDS = ["source_id", "classifier_name", "best_class_name", "best_class_score"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tap(query: str) -> bytes:
    form = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "MAXREC": "5000", "QUERY": query,
    }).encode()
    request = urllib.request.Request(ENDPOINT, form)
    with urllib.request.urlopen(request, timeout=180) as response:
        raw = response.read(5 * 1024 * 1024 + 1)
    if len(raw) > 5 * 1024 * 1024:
        raise ValueError("Gaia response exceeded 5 MiB safety limit")
    return raw


def parse(raw: bytes, required: list[str], label: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not set(required).issubset(reader.fieldnames or []):
        raise ValueError(f"unexpected {label} response columns: {reader.fieldnames}")
    rows = list(reader)
    if len(rows) >= 5000:
        raise ValueError(f"possible MAXREC truncation in {label} response")
    return rows


def run(candidates_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia variability output: {output}")
    with candidates_path.open(encoding="utf-8", newline="") as stream:
        candidates = list(csv.DictReader(stream, delimiter="\t"))
    if not candidates or len(candidates) != 480:
        raise ValueError(f"expected the approved 480-candidate review list; got {len(candidates)}")
    by_source: dict[str, dict[str, str]] = {}
    for row in candidates:
        source_id = row["gaia_source_id"]
        if not source_id or source_id in by_source:
            raise ValueError("candidate list has missing or duplicate Gaia source IDs")
        by_source[source_id] = row
    ids = sorted(by_source, key=int)
    id_list = ",".join(ids)
    summary_query = (
        "SELECT " + ", ".join(SUMMARY_FIELDS) +
        " FROM gaiadr3.vari_summary WHERE source_id IN (" + id_list + ")"
    )
    classifier_query = (
        "SELECT " + ", ".join(CLASSIFIER_FIELDS) +
        " FROM gaiadr3.vari_classifier_result WHERE source_id IN (" + id_list + ")"
    )
    output.mkdir(parents=True)
    (output / "query_vari_summary.adql").write_text(summary_query + "\n", encoding="utf-8")
    (output / "query_vari_classifier_result.adql").write_text(classifier_query + "\n", encoding="utf-8")
    summary_raw = tap(summary_query)
    classifier_raw = tap(classifier_query)
    (output / "vari_summary.csv").write_bytes(summary_raw)
    (output / "vari_classifier_result.csv").write_bytes(classifier_raw)
    summary_rows = parse(summary_raw, SUMMARY_FIELDS, "vari_summary")
    classifier_rows = parse(classifier_raw, CLASSIFIER_FIELDS, "vari_classifier_result")
    summary_by_id = {row["source_id"]: row for row in summary_rows}
    if len(summary_by_id) != len(summary_rows) or not set(summary_by_id).issubset(by_source):
        raise ValueError("unexpected or duplicate source IDs in variability summary")
    classifier_by_id: dict[str, list[dict[str, str]]] = {}
    for row in classifier_rows:
        if row["source_id"] not in by_source:
            raise ValueError("unexpected source ID in variability classifier response")
        classifier_by_id.setdefault(row["source_id"], []).append(row)
    result_fields = ["candidate_id", "cluster_id", "gaia_source_id", "gaia_g_mag",
                     "gaia_phot_variable_flag", *SUMMARY_FIELDS[1:],
                     "classifier_name", "best_class_name", "best_class_score"]
    merged = []
    for candidate in candidates:
        source_id = candidate["gaia_source_id"]
        summaries = summary_by_id.get(source_id)
        classifiers = classifier_by_id.get(source_id, [])
        if len(classifiers) > 1:
            raise ValueError(f"multiple classifier rows for source {source_id}")
        classifier = classifiers[0] if classifiers else {}
        merged.append({
            "candidate_id": candidate["candidate_id"],
            "cluster_id": candidate["cluster_id"],
            "gaia_source_id": source_id,
            "gaia_g_mag": candidate["gaia_g_mag"],
            "gaia_phot_variable_flag": candidate["gaia_phot_variable_flag"],
            **{field: summaries.get(field, "") if summaries else "" for field in SUMMARY_FIELDS[1:]},
            **{field: classifier.get(field, "") for field in CLASSIFIER_FIELDS[1:]},
        })
    with (output / "candidate_variability.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=result_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(merged)
    payload = {
        "status": "complete", "catalogue": "Gaia DR3",
        "tables": ["gaiadr3.vari_summary", "gaiadr3.vari_classifier_result"],
        "n_input_candidates": len(candidates), "n_unique_source_ids": len(ids),
        "n_vari_summary_rows": len(summary_rows),
        "n_classifier_rows": len(classifier_rows),
        "n_candidates_with_vari_summary": len(summary_by_id),
        "n_candidates_with_classifier": len(classifier_by_id),
        "input_sha256": sha256(candidates_path),
        "summary_query_sha256": sha256(output / "query_vari_summary.adql"),
        "classifier_query_sha256": sha256(output / "query_vari_classifier_result.adql"),
        "summary_response_sha256": sha256(output / "vari_summary.csv"),
        "classifier_response_sha256": sha256(output / "vari_classifier_result.csv"),
        "result_sha256": sha256(output / "candidate_variability.tsv"),
        "images_requested": False, "hsc_catalogue_rows_requested": False,
        "science_measurement": False, "psf_approved": False,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (output / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                                          encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True,
                        help="the existing 480-row Gaia-supported candidate TSV")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.candidates, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
