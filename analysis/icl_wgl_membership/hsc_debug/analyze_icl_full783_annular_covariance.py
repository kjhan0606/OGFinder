#!/usr/bin/env python3
"""Measure ICL-annulus coadd covariance for all 783 same-filter pairs.

This offline Astrafek analysis extends the 94-pair target/control annular
covariance diagnostic to the complete existing same-filter image sample. It
reads only locally saved FITS images, masks, and manifests; it makes no HSC
requests and does not promote science_measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import analyze_icl_annular_coadd_covariance as annular
import audit_icl_full_sample_coadd_covariance as full


DEFAULT_OUTPUT = full.DEFAULT_SAMPLE / "icl_annular_coadd_covariance_full783_2026-11-06_02"
EXPECTED_FILTER_COUNTS = {"i": 302, "i2": 481}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> str:
    if not rows:
        raise ValueError(f"refusing empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        fieldnames = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    return sha256(path)


def read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def append_jsonl(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("a", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()


def write_progress(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def build_pairs(sample: Path) -> tuple[list[dict], list[dict], dict]:
    records, setup = full.build_records(sample)
    manifest_path = sample / "download_01/download_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in manifest["clusters"]}
    pairs = []
    seen = set()
    for record in records:
        if record["role"] != "target":
            continue
        ident = str(record["cluster_id"])
        if ident in seen or ident not in targets:
            raise ValueError(f"duplicate or missing target manifest row: {ident}")
        seen.add(ident)
        target = targets[ident]
        pair = next(row for row in setup["pair_rows"] if str(row["cluster_id"]) == ident)
        pair["redshift"] = float(target["z"])
        pairs.append(pair)
    counts = {band: sum(row["filter"] == band for row in pairs) for band in ("i", "i2")}
    if len(pairs) != 783 or counts != EXPECTED_FILTER_COUNTS:
        raise ValueError(f"same-filter pair scope changed: {len(pairs)} rows, {counts}")
    paired = {(row["cluster_id"], row["filter"]) for row in pairs}
    selected = [record for record in records
                if (str(record["cluster_id"]), record["filter"]) in paired]
    if len(selected) != 1566:
        raise ValueError(f"expected 1,566 target/control image records, got {len(selected)}")
    return pairs, selected, {
        **setup["provenance"],
        str(manifest_path): sha256(manifest_path),
        "target_control_pair_counts": counts,
        "pair_count": len(pairs),
        "connected_5deg_blocks": len({row["connected_5deg_block"] for row in pairs}),
    }


def run(sample: Path, output: Path, *, resume: bool = False,
        patch_size: int = 64, max_patches: int = 32,
        min_patches: int = 4, min_valid_fraction: float = 0.80,
        max_lag: int = 8) -> dict:
    sample, output = sample.resolve(), output.resolve()
    if patch_size < 2 * max_lag + 1:
        raise ValueError("patch size must exceed the full lag-kernel width")
    annular.PATCH_SIZE = int(patch_size)
    annular.MAX_PATCHES = int(max_patches)
    annular.MIN_PATCHES = int(min_patches)
    annular.MIN_VALID_FRACTION = float(min_valid_fraction)
    annular.MAX_LAG = int(max_lag)
    pairs, records, provenance = build_pairs(sample)
    provenance["patch_configuration"] = {
        "patch_size_pixels": annular.PATCH_SIZE,
        "max_patches_per_center_annulus": annular.MAX_PATCHES,
        "minimum_patches_per_center_annulus": annular.MIN_PATCHES,
        "minimum_valid_fraction_per_patch": annular.MIN_VALID_FRACTION,
        "max_lag_pixels": annular.MAX_LAG,
    }
    progress_path = output / "progress.json"
    center_jsonl = output / "center_results.jsonl"
    lag_jsonl = output / "lag_results.jsonl"
    failure_jsonl = output / "failures.jsonl"
    processed_jsonl = output / "processed_records.jsonl"
    if output.exists():
        if not resume or not progress_path.is_file():
            raise FileExistsError(f"refusing existing output without --resume checkpoint: {output}")
        previous = json.loads(progress_path.read_text(encoding="utf-8"))
        if (previous.get("schema") != "astrafek-icl-full783-annular-progress-v1"
                or previous.get("input_sha256") != provenance
                or previous.get("status") not in {"running", "interrupted"}):
            raise ValueError("existing full-sample annular checkpoint does not match this exact input")
    else:
        if resume:
            raise FileNotFoundError(f"cannot resume missing output directory: {output}")
        output.mkdir(parents=True, exist_ok=False)
        write_progress(progress_path, {
            "schema": "astrafek-icl-full783-annular-progress-v1",
            "status": "running", "input_sha256": provenance,
            "expected_records": len(records), "processed_records": 0,
        })
    rows = read_jsonl(center_jsonl)
    lag_rows = read_jsonl(lag_jsonl)
    failures = read_jsonl(failure_jsonl)
    processed = {str(row["record_key"]) for row in read_jsonl(processed_jsonl)}
    processed.update(f"{row['cluster_id']}|{row['filter']}|{row['role']}" for row in rows)
    processed.update(f"{row['cluster_id']}|{row['filter']}|{row['role']}" for row in failures)
    pair_lookup = {(str(row["cluster_id"]), row["filter"]): row for row in pairs}
    for index, record in enumerate(records, 1):
        ident = str(record["cluster_id"])
        pair = pair_lookup[(ident, record["filter"])]
        record_key = f"{ident}|{record['filter']}|{record['role']}"
        if record_key in processed:
            continue
        print(f"[{index}/{len(records)}] {record['filter']} {ident} {record['role']}", flush=True)
        try:
            prepared = annular.prepare_center(record, pair)
            summary, lags = annular.center_annulus(
                prepared, "icl_200_300", *annular.ANNULI["icl_200_300"])
            append_jsonl(center_jsonl, [summary])
            append_jsonl(lag_jsonl, lags)
            rows.append(summary)
            lag_rows.extend(lags)
        except Exception as exc:
            failure = {"cluster_id": ident, "filter": record["filter"],
                       "role": record["role"],
                       "error": f"{type(exc).__name__}: {exc}"}
            append_jsonl(failure_jsonl, [failure])
            failures.append(failure)
        append_jsonl(processed_jsonl, [{"record_key": record_key}])
        processed.add(record_key)
        if index % 10 == 0 or index == len(records):
            write_progress(progress_path, {
                "schema": "astrafek-icl-full783-annular-progress-v1",
                "status": "running", "input_sha256": provenance,
                "expected_records": len(records), "processed_records": len(processed),
                "last_record": record_key,
            })
    lookup = {(str(row["cluster_id"]), row["filter"], row["role"]): row for row in rows}
    comparisons = []
    for pair in pairs:
        ident, band = str(pair["cluster_id"]), pair["filter"]
        target = lookup.get((ident, band, "target"))
        control = lookup.get((ident, band, "matched_control"))
        usable = (target is not None and control is not None
                  and target["status"] == "ok" and control["status"] == "ok")
        row = {"cluster_id": ident, "filter": band,
               "connected_5deg_block": pair["connected_5deg_block"],
               "complete_target_control_support": usable}
        for key, item in (("target", target), ("control", control)):
            for metric, column in (("zero_lag", "zero_lag_variance_ratio"),
                                   ("sum8", "square_covariance_sum_8px"),
                                   ("coverage", "valid_annulus_coverage"),
                                   ("n_patches", "n_patches_used")):
                row[f"{key}_{metric}"] = item.get(column) if item else None
        row["target_minus_control_zero_lag"] = (
            target["zero_lag_variance_ratio"] - control["zero_lag_variance_ratio"]
            if usable else None)
        row["target_minus_control_sum8"] = (
            target["square_covariance_sum_8px"] - control["square_covariance_sum_8px"]
            if usable else None)
        row["target_minus_control_coverage"] = (
            target["valid_annulus_coverage"] - control["valid_annulus_coverage"]
            if usable else None)
        comparisons.append(row)
    bootstrap = {}
    for index, band in enumerate(("i", "i2")):
        selected = [row for row in comparisons
                    if row["filter"] == band and row["complete_target_control_support"]]
        bootstrap[band] = {}
        for metric_index, metric in enumerate(("zero_lag", "sum8", "coverage")):
            vals = [(row["connected_5deg_block"], row[f"target_minus_control_{metric}"])
                    for row in selected if row[f"target_minus_control_{metric}"] is not None]
            bootstrap[band][metric] = annular.block_bootstrap(
                vals, annular.SEED + 70000 + index * 100 + metric_index)
        bootstrap[band]["support"] = {
            "n_same_filter_pairs": EXPECTED_FILTER_COUNTS[band],
            "n_complete_target_control_pairs": len(selected),
            "n_connected_blocks_with_complete_pairs": len({r["connected_5deg_block"] for r in selected}),
            "n_target_insufficient_patch_support": sum(
                1 for r in rows if r["filter"] == band and r["role"] == "target"
                and r["status"] == "insufficient_patch_support"),
            "n_control_insufficient_patch_support": sum(
                1 for r in rows if r["filter"] == band and r["role"] == "matched_control"
                and r["status"] == "insufficient_patch_support"),
        }
    per_center_path = output / "per_center_icl_annulus_covariance.tsv"
    lag_path = output / "icl_annulus_lag_covariance.tsv"
    pair_path = output / "paired_icl_annulus_covariance.tsv"
    output_hashes = {
        str(per_center_path): write_tsv(per_center_path, rows),
        str(lag_path): write_tsv(lag_path, lag_rows),
        str(pair_path): write_tsv(pair_path, comparisons),
    }
    summary = {
        "schema": "astrafek-icl-full783-annular-coadd-covariance-v1",
        "status": "complete_exploratory" if not failures else "complete_with_data_gaps",
        "science_measurement": False, "hsc_network_used": False,
        "image_pixels_read": True,
        "n_same_filter_pairs": len(pairs), "n_image_records": len(records),
        "n_analyzed_center_annuli": len(rows),
        "n_successful_center_annuli": sum(row["status"] == "ok" for row in rows),
        "failed_inputs": failures,
        "annulus_proper_kpc": list(annular.ANNULI["icl_200_300"]),
        "connected_5deg_blocks": provenance["connected_5deg_blocks"],
        "method": {
            "base_method": "analyze_icl_annular_coadd_covariance.py",
            "patch_size_pixels": annular.PATCH_SIZE,
            "max_lag_pixels": annular.MAX_LAG,
            "max_patches_per_center_annulus": annular.MAX_PATCHES,
            "minimum_valid_fraction_per_patch": annular.MIN_VALID_FRACTION,
            "minimum_patches_per_center_annulus": annular.MIN_PATCHES,
            "patches": f"non-overlapping {annular.PATCH_SIZE}x{annular.PATCH_SIZE} fully contained in the 200-300 proper-kpc annulus",
            "estimator": "local affine-plane removal, clipping, and variance-normalized lag covariance",
            "bootstrap": "paired target-control differences resampled by connected 5-degree blocks; 5,000 replicates",
        },
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "paired_block_bootstrap": bootstrap,
        "limitations": [
            "variance-normalized covariance is a delivered-coadd diagnostic, not exact aperture covariance",
            "the 783 same-filter pairs share only 38 approximate connected sky blocks, with uneven block sizes",
            "does not separate resampling covariance, residual sky/source structure, variance calibration, or astrophysical covariance",
            "does not model source-Poisson noise, PSF wings, sky transfer, or total ICL measurement covariance",
            "no empirical covariance correction is applied to ICL or weak-lensing profiles",
        ],
        "input_sha256": provenance,
        "outputs_sha256": output_hashes,
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    write_progress(progress_path, {
        "schema": "astrafek-icl-full783-annular-progress-v1",
        "status": "complete", "input_sha256": provenance,
        "expected_records": len(records), "processed_records": len(processed),
        "summary_sha256": sha256(summary_path),
    })
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "pairs": len(pairs), "successful_center_annuli": summary["n_successful_center_annuli"],
                      "complete_pair_support": {band: data["support"] for band, data in bootstrap.items()},
                      "summary": str(summary_path)}, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=full.DEFAULT_SAMPLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patch-size", type=int, default=64)
    parser.add_argument("--max-patches", type=int, default=32)
    parser.add_argument("--min-patches", type=int, default=4)
    parser.add_argument("--min-valid-fraction", type=float, default=0.80)
    parser.add_argument("--max-lag", type=int, default=8)
    parser.add_argument("--resume", action="store_true",
                        help="resume only an exact matching interrupted run")
    args = parser.parse_args()
    run(args.sample, args.output_dir, resume=args.resume,
        patch_size=args.patch_size, max_patches=args.max_patches,
        min_patches=args.min_patches, min_valid_fraction=args.min_valid_fraction,
        max_lag=args.max_lag)


if __name__ == "__main__":
    main()
