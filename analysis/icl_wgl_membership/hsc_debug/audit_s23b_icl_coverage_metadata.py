#!/usr/bin/env python3
"""Audit the approved S23B HSC-I footprint metadata export.

The audit reads only the pinned plan and returned coverage metadata. It makes
no network request and reads no image pixels; it cannot establish full-aperture
WCS/mask support or matched-control availability.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


PLAN_SHA256 = "62c28446f553bf45a940ea0126af1521d4ee0c1356c12cc73f42712df482c0c7"
COLUMNS = [
    "candidate_id", "probe_id", "skymap_id", "tract", "patch", "band",
    "naxis1", "naxis2", "seeing", "zeropt", "imag_psf_depth", "iseeing",
    "imeas_exists", "idet_bkgd_exists",
]
PROBE_IDS = ["center"] + [f"r1Mpc_pa{bearing:03d}" for bearing in range(0, 360, 45)]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = [name.strip().removeprefix("#").strip() for name in (reader.fieldnames or [])]
        if len(fields) != len(set(fields)):
            raise ValueError(f"duplicate returned columns: {path.name}")
        reader.fieldnames = fields
        return fields, list(reader)


def as_bool(value: str, name: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "t", "1"}:
        return True
    if normalized in {"false", "f", "0"}:
        return False
    raise ValueError(f"invalid boolean {name}={value!r}")


def run(plan_dir: Path, output_dir: Path) -> dict:
    plan_dir = plan_dir.resolve()
    result_dir = plan_dir / "cas_coverage_metadata"
    plan_path = plan_dir / "coverage_query_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_sha = sha(plan_path)
    if (plan_sha != PLAN_SHA256 or plan.get("schema") != "astrafek-s23b-pdr3-coverage-query-plan-v1"
            or plan.get("n_candidates") != 1333 or plan.get("n_probes") != 11997
            or plan.get("n_sql_batches") != 27
            or plan.get("coverage_test", {}).get("output_columns") != COLUMNS):
        raise ValueError("S23B plan differs from the approved metadata scope")
    status_path = result_dir / "coverage_status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    summary_path = result_dir / "coverage_summary.json"
    returned_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (status.get("status") != "complete" or status.get("plan_sha256") != plan_sha
            or len(status.get("jobs", {})) != 27
            or any(item.get("status") != "complete" for item in status["jobs"].values())
            or returned_summary.get("plan_sha256") != plan_sha
            or returned_summary.get("output_columns") != COLUMNS
            or returned_summary.get("n_candidates") != 1333
            or returned_summary.get("n_probes") != 11997
            or returned_summary.get("n_metadata_rows") != sum(item.get("rows", 0) for item in status.get("jobs", {}).values())
            or any(returned_summary.get(key) is not False for key in (
                "images_returned", "object_rows_returned", "shape_rows_returned",
                "photoz_rows_returned", "science_measurement"))):
        raise ValueError("S23B metadata export is incomplete or exceeded its approved scope")

    centers_path = plan_dir / "candidate_centers.tsv"
    probes_path = plan_dir / "coverage_probes.tsv"
    centers = read_tsv(centers_path)
    probes = read_tsv(probes_path)
    candidate_ids = {row["candidate_id"] for row in centers}
    expected_probes = {(row["candidate_id"], row["probe_id"]) for row in probes}
    if (len(centers) != 1333 or len(candidate_ids) != 1333 or len(probes) != 11997
            or len(expected_probes) != 11997
            or any(not candidate.startswith("S23B:") for candidate in candidate_ids)):
        raise ValueError("candidate or probe inventory is incomplete, duplicated, or unprefixed")
    if {(candidate, probe) for candidate in candidate_ids for probe in PROBE_IDS} != expected_probes:
        raise ValueError("every candidate must have exactly the center and eight boundary probes")

    rows: list[dict[str, str]] = []
    for index, batch in enumerate(plan["batches"]):
        batch_name = f"coverage_mosaic_{index:03d}"
        item = status["jobs"].get(batch_name)
        path = result_dir / f"{batch_name}.csv"
        if (item is None or not path.is_file() or sha(path) != item.get("sha256")
                or path.stat().st_size != item.get("bytes")):
            raise ValueError(f"batch checksum/size mismatch: {batch_name}")
        fields, batch_rows = read_csv(path)
        if fields != COLUMNS or len(batch_rows) != item.get("rows"):
            raise ValueError(f"batch schema/row count mismatch: {batch_name}")
        first = sum(int(b["n_candidates"]) for b in plan["batches"][:index])
        stop = first + int(batch["n_candidates"])
        batch_candidates = {row["candidate_id"] for row in centers[first:stop]}
        allowed = {(candidate, probe) for candidate, probe in expected_probes
                   if candidate in batch_candidates}
        seen = set()
        for row in batch_rows:
            key = (row["candidate_id"], row["probe_id"])
            unique_key = (*key, row["skymap_id"])
            if (key not in allowed or unique_key in seen or row["band"] != "HSC-I"
                    or not row["candidate_id"].startswith("S23B:")
                    or not row["tract"].isdigit() or not row["patch"]):
                raise ValueError(f"unexpected/duplicate/non-HSC-I row in {batch_name}")
            seen.add(unique_key)
            for field in ("naxis1", "naxis2"):
                if int(row[field]) <= 0:
                    raise ValueError(f"invalid {field} in {batch_name}")
            for field in ("seeing", "zeropt"):
                if float(row[field]) <= 0:
                    raise ValueError(f"invalid {field} in {batch_name}")
            for field in ("imag_psf_depth", "iseeing"):
                if row[field].strip().lower() not in {"", "null", r"\n", r"\N"} and float(row[field]) <= 0:
                    raise ValueError(f"invalid {field} in {batch_name}")
            as_bool(row["imeas_exists"], "imeas_exists")
            as_bool(row["idet_bkgd_exists"], "idet_bkgd_exists")
            rows.append(row)

    probe_rows: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        probe_rows.setdefault((row["candidate_id"], row["probe_id"]), []).append(row)
    by_candidate = {candidate: [] for candidate in candidate_ids}
    for (candidate, probe), matches in probe_rows.items():
        by_candidate[candidate].append((probe, matches))
    per_candidate = []
    for center in centers:
        candidate = center["candidate_id"]
        matches = {probe: table for probe, table in by_candidate[candidate]}
        covered = set(matches)
        image_points = {
            probe for probe, table in matches.items()
            if any(as_bool(row["imeas_exists"], "imeas_exists") for row in table)
        }
        background_points = {
            probe for probe, table in matches.items()
            if any(as_bool(row["idet_bkgd_exists"], "idet_bkgd_exists") for row in table)
        }
        values = [float(row["imag_psf_depth"]) for table in matches.values()
                  for row in table if row["imag_psf_depth"].strip().lower() not in {"", "null", r"\N"}]
        seeing = [float(row["seeing"]) for table in matches.values() for row in table]
        per_candidate.append({
            "candidate_id": candidate,
            "covered_probe_count": len(covered),
            "covered_boundary_probe_count": sum(probe in covered for probe in PROBE_IDS[1:]),
            "center_covered": "center" in covered,
            "all_nine_probes_covered": len(covered) == 9,
            "image_product_probe_count": len(image_points),
            "all_nine_have_image_product_flag": len(image_points) == 9,
            "background_product_probe_count": len(background_points),
            "all_nine_have_background_product_flag": len(background_points) == 9,
            "n_returned_patch_rows": sum(map(len, matches.values())),
            "min_imag_psf_depth": min(values) if values else None,
            "median_imag_psf_depth": sorted(values)[len(values)//2] if values else None,
            "median_seeing_arcsec": sorted(seeing)[len(seeing)//2] if seeing else None,
            "science_measurement": False,
        })

    # Summarize overlap-aware probe counts; overlapping adjacent patches may
    # return multiple rows for a single sky position.
    fully_covered = [row for row in per_candidate if row["all_nine_probes_covered"]]
    all_image = [row for row in per_candidate if row["all_nine_have_image_product_flag"]]
    all_background = [row for row in per_candidate if row["all_nine_have_background_product_flag"]]
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    per_candidate_path = output_dir / "per_candidate_coverage.tsv"
    fields = list(per_candidate[0])
    with per_candidate_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(per_candidate)
    summary = {
        "schema": "astrafek-s23b-icl-coverage-metadata-audit-v1",
        "status": "complete_exploratory_audit",
        "science_measurement": False,
        "runner_sha256": sha(Path(__file__).resolve()),
        "hsc_network_used_by_audit": False,
        "image_pixels_read": False,
        "plan_sha256": plan_sha,
        "n_candidates": len(centers),
        "n_probes": len(probes),
        "n_returned_patch_metadata_rows": len(rows),
        "n_unique_probes_with_match": len(probe_rows),
        "n_center_positions_with_match": sum((candidate, "center") in probe_rows for candidate in candidate_ids),
        "n_candidates_all_9_probes_covered": len(fully_covered),
        "n_candidates_all_9_probes_with_imeas_exists": len(all_image),
        "n_candidates_all_9_probes_with_idet_bkgd_exists": len(all_background),
        "returned_patch_row_count_per_probe": {
            "min": min(map(len, probe_rows.values())),
            "median": sorted(map(len, probe_rows.values()))[len(probe_rows)//2],
            "max": max(map(len, probe_rows.values())),
        },
        "limitations": [
            "nine points are a footprint screen, not proof of full-disk or pixel-level WCS/mask coverage",
            "patch_qa flags and depth/seeing are metadata, not image validation",
            "no matched controls were constructed or checked for these new candidates",
            "the S23B catalog candidates are not merged with S20A by numeric ID",
            "image downloads require their separate explicit user approval",
        ],
        "input_sha256": {
            str(plan_path): plan_sha,
            str(status_path): sha(status_path),
            str(summary_path): sha(summary_path),
            str(centers_path): sha(centers_path),
            str(probes_path): sha(probes_path),
        },
        "batch_sha256": {
            str(result_dir / f"coverage_mosaic_{index:03d}.csv"): status["jobs"][f"coverage_mosaic_{index:03d}"]["sha256"]
            for index in range(27)
        },
        "output_sha256": {per_candidate_path.name: sha(per_candidate_path)},
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in (
        "status", "n_candidates", "n_probes", "n_returned_patch_metadata_rows",
        "n_unique_probes_with_match", "n_center_positions_with_match",
        "n_candidates_all_9_probes_covered", "n_candidates_all_9_probes_with_imeas_exists",
        "n_candidates_all_9_probes_with_idet_bkgd_exists")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.plan_dir, args.output_dir)


if __name__ == "__main__":
    main()
