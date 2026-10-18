#!/usr/bin/env python3
"""Describe how ICL annulus support cuts select the 800-cluster parent sample.

This is a selection audit, not an ICL measurement or a correction for selection.
It compares the fixed all-annulus K2.5/40% sample with the independent
200--300 kpc annulus sample, within each observed filter, and records all input
hashes so later analyses can trace the exact selection evidence.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


STRICT_FIELDS = ("delta_sb_50_100", "delta_sb_100_150",
                 "delta_sb_150_200", "delta_sb_200_300")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def parse_bool(value):
    value = str(value).strip().lower()
    if value in {"true", "1"}:
        return True
    if value in {"false", "0"}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def describe(rows, reference):
    """Summarize a selected group against same-filter parent values."""
    if not rows:
        raise ValueError("Cannot summarize an empty selection")
    metrics = {}
    for field in ("z", "richness"):
        selected = [float(row[field]) for row in rows]
        parent = [float(row[field]) for row in reference]
        parent_sd = statistics.stdev(parent) if len(parent) > 1 else 0.
        metrics[field] = {
            "mean": statistics.fmean(selected),
            "median": statistics.median(selected),
            "sd": statistics.stdev(selected) if len(selected) > 1 else 0.,
            "parent_mean": statistics.fmean(parent),
            "parent_sd": parent_sd,
            "standardized_mean_difference": ((statistics.fmean(selected) - statistics.fmean(parent)) / parent_sd
                                             if parent_sd > 0 else None),
        }
    center_ok = [parse_bool(row["auto_unique_center_candidate"]) for row in rows]
    center_masked = [parse_bool(row["center_masked_k2p5"]) for row in rows]
    cells = [f"{int(float(row['ra_deg']) // 5)}_{int((float(row['dec_deg']) + 90) // 5)}"
             for row in rows]
    counts = {}
    for cell in cells:
        counts[cell] = counts.get(cell, 0) + 1
    metrics.update({
        "n_clusters": len(rows),
        "unique_5deg_sky_blocks": len(counts),
        "largest_block_n": max(counts.values()),
        "largest_block_fraction": max(counts.values()) / len(rows),
        "auto_unique_center_fraction": sum(center_ok) / len(rows),
        "center_masked_k2p5_fraction": sum(center_masked) / len(rows),
    })
    return metrics


def validate_inputs(manifest, profiles, strict, annulus):
    frames = manifest.get("clusters", [])
    if len(frames) != 800 or len({str(row["id"]) for row in frames}) != 800:
        raise ValueError("Parent sample must contain 800 unique clusters")
    profile_by_id = {str(row["id"]): row for row in profiles}
    parent_by_id = {str(row["id"]): row for row in frames}
    if len(profile_by_id) != 800 or set(profile_by_id) != set(parent_by_id):
        raise ValueError("Cluster-profile table does not match the 800-cluster manifest")

    selected_strict = []
    for row in strict:
        if all(row.get(field, "").strip() for field in STRICT_FIELDS):
            if row["id"] in parent_by_id:
                selected_strict.append(row)
    if len(selected_strict) != 137 or len({row["id"] for row in selected_strict}) != 137:
        raise ValueError(f"Expected 137 unique strict all-annulus pairs; got {len(selected_strict)}")

    selected_outer = [row for row in annulus
                      if float(row["kron_factor"]) == 2.5
                      and float(row["coverage_floor"]) == .40
                      and row["radial_bin_kpc"] == "200-300"
                      and float(row["common_coverage"]) >= .40]
    if len(selected_outer) != 224 or len({row["id"] for row in selected_outer}) != 224:
        raise ValueError(f"Expected 224 unique 200--300 kpc annulus pairs; got {len(selected_outer)}")

    def enrich(row, source):
        ident = str(row["id"])
        if ident not in profile_by_id:
            raise ValueError(f"Selected cluster missing from parent profile table: {ident}")
        profile = profile_by_id[ident]
        if row["filter"] != profile["filter"]:
            raise ValueError(f"Filter mismatch for {ident}")
        if not math.isclose(float(row["z"]), float(profile["z"]), abs_tol=1e-7):
            raise ValueError(f"Redshift mismatch for {ident}")
        if not math.isclose(float(row["richness"]), float(profile["richness"]), abs_tol=1e-6):
            raise ValueError(f"Richness mismatch for {ident}")
        enriched = dict(row)
        enriched["ra_deg"] = profile["ra_deg"]
        enriched["dec_deg"] = profile["dec_deg"]
        if "center_masked_k2p5" not in enriched:
            # The annulus-sensitivity table stores the already-resolved flag,
            # while the strict table uses an explicit K2.5 column.
            enriched["center_masked_k2p5"] = enriched["center_masked_by_kron"]
        return enriched

    selected_strict = [enrich(row, "strict") for row in selected_strict]
    selected_outer = [enrich(row, "outer") for row in selected_outer]
    parent_rows = []
    for ident, row in profile_by_id.items():
        flags = json.loads(row["center_masked_by_kron"])
        parent_rows.append({
            "id": ident, "filter": row["filter"], "z": row["z"], "richness": row["richness"],
            "ra_deg": row["ra_deg"], "dec_deg": row["dec_deg"],
            "auto_unique_center_candidate": row["auto_unique_center_candidate"],
            "center_masked_k2p5": flags["2.5"],
        })
    return parent_rows, selected_strict, selected_outer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--strict-pairs", type=Path, required=True)
    parser.add_argument("--annulus-pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Refusing to overwrite existing output: {output}")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    profiles = read_tsv(args.profiles)
    strict = read_tsv(args.strict_pairs)
    annulus = read_tsv(args.annulus_pairs)
    parent, strict_rows, outer_rows = validate_inputs(manifest, profiles, strict, annulus)

    groups = {"parent": parent, "strict_50_300_k2p5_c40": strict_rows,
              "annulus_200_300_k2p5_c40": outer_rows}
    stats = []
    result = {"schema": "hsc-icl-sample-selection-audit-v1", "status": "complete_descriptive",
              "science_measurement": False, "selection_correction_applied": False,
              "analysis_script_sha256": sha256(Path(__file__)),
              "n_parent": len(parent), "n_strict": len(strict_rows), "n_outer_annulus": len(outer_rows),
              "filters": {}, "input_sha256": {
                  str(args.manifest.resolve()): sha256(args.manifest),
                  str(args.profiles.resolve()): sha256(args.profiles),
                  str(args.strict_pairs.resolve()): sha256(args.strict_pairs),
                  str(args.annulus_pairs.resolve()): sha256(args.annulus_pairs),
              }, "limitations": [
                  "Selection descriptors do not correct the stack or identify the physical ICL signal.",
                  "Strict and annulus-only samples use different radial support criteria and are not interchangeable.",
                  "Five-degree cells summarize footprint concentration; they are not an uncertainty model.",
                  "Central masking and automatic center flags remain imperfect image-quality indicators.",
              ]}
    for band in ("i", "i2"):
        parent_band = [row for row in parent if row["filter"] == band]
        result["filters"][band] = {"parent_n": len(parent_band)}
        for name, rows in groups.items():
            selected = [row for row in rows if row["filter"] == band]
            if not selected:
                continue
            summary = describe(selected, parent_band)
            result["filters"][band][name] = summary
            for metric in ("z", "richness"):
                stats.append({"filter": band, "selection": name, "metric": metric,
                              **summary[metric]})

    output.mkdir(parents=True, exist_ok=False)
    stats_path = output / "selection_stats.tsv"
    fields = ["filter", "selection", "metric", "mean", "median", "sd", "parent_mean",
              "parent_sd", "standardized_mean_difference"]
    with stats_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader(); writer.writerows(stats)
    result["output_sha256"] = {stats_path.name: sha256(stats_path)}
    (output / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2,
        allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "parent": len(parent),
        "strict": len(strict_rows), "outer_annulus": len(outer_rows), "output": str(output)},
        ensure_ascii=False))


if __name__ == "__main__":
    main()
