#!/usr/bin/env python3
"""Pair full-sample mask coverage with existing controls; do not measure flux."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from prestack_core import EDGES, FACTORS, spatial_labels


THRESHOLDS = (0.25, 0.35, 0.40, 0.50)
SCOPES = {
    "50-150": [0, 1],
    "50-200": [0, 1, 2],
    "50-300": [0, 1, 2, 3],
    "50-500": [0, 1, 2, 3, 4, 5],
    **{f"{lo}-{hi}": [i] for i, (lo, hi) in enumerate(
        zip((50, 100, 150, 200, 300, 400), (100, 150, 200, 300, 400, 500))
    )},
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path):
    return json.loads(Path(path).read_text())


def cell(ra, dec):
    return str(spatial_labels([ra], [dec], 5.0)[0])


def joint_block_count(target_rows, control_rows):
    """Connected groups link pairs sharing a 5-degree cell on either side."""
    parent = list(range(len(target_rows)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    owner = {}
    for i, (target, control) in enumerate(zip(target_rows, control_rows)):
        for key in (("target", cell(target["ra_deg"], target["dec_deg"])),
                    ("control", cell(control["ra_deg"], control["dec_deg"]))):
            if key in owner:
                union(i, owner[key])
            else:
                owner[key] = i
    return len({find(i) for i in range(len(parent))})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--cluster-masks", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output directory: {args.work}")

    target_manifest_path = args.sample / "download_01/download_manifest.json"
    control_manifest_path = args.sample / "control_patches_2026-10-14_02/control_download_manifest.json"
    status_path = args.sample / "download_01/download_status.json"
    extraction_summary_path = args.extraction / "extraction_summary.json"
    mask_summary_path = args.cluster_masks / "summary.json"
    coverage_path = args.cluster_masks / "coverage.tsv"
    records_path = args.cluster_masks / "records.tsv"
    target_rows = read_json(target_manifest_path)["clusters"]
    control_rows = read_json(control_manifest_path)["clusters"]
    download_status = read_json(status_path)
    extraction_summary = read_json(extraction_summary_path)
    mask_summary = read_json(mask_summary_path)
    if len(target_rows) != 800 or mask_summary.get("status") != "complete" or mask_summary.get("n_success") != 800:
        raise ValueError("Expected complete 800-target mask support")
    if sha256(coverage_path) != mask_summary["output_sha256"].get(str(coverage_path)):
        raise ValueError("Cluster coverage table hash changed")
    if sha256(records_path) != mask_summary["output_sha256"].get(str(records_path)):
        raise ValueError("Cluster mask record table hash changed")

    coverage_rows = read_tsv(coverage_path)
    by_target_coverage = {}
    for row in coverage_rows:
        ident = row["id"]
        fi = FACTORS.index(float(row["kron_factor"]))
        edge_pair = (float(row["r_min_kpc"]), float(row["r_max_kpc"]))
        ring_pairs = list(zip(EDGES[:-1], EDGES[1:]))
        if edge_pair not in ring_pairs:
            raise ValueError(f"Unexpected coverage ring: {edge_pair}")
        ri = ring_pairs.index(edge_pair)
        by_target_coverage.setdefault(ident, np.full((3, 6), np.nan))[fi, ri] = float(row["coverage"])
    if set(by_target_coverage) != {str(row["id"]) for row in target_rows}:
        raise ValueError("Cluster coverage IDs differ from the 800-target manifest")

    cluster_records = {}
    for row in read_tsv(records_path):
        cluster_records[row["id"]] = dict(
            unique=bool(row["auto_unique_center_candidate"] == "True"),
            base_masked=bool(row["center_base_masked"] == "True"),
            center_masked=ast.literal_eval(row["center_masked_by_kron"]),
        )
    if set(cluster_records) != set(by_target_coverage):
        raise ValueError("Cluster mask record IDs differ from coverage IDs")

    controls_by_target = {}
    for row in control_rows:
        key = str(row["target_cluster_id"])
        if key in controls_by_target:
            raise ValueError(f"Duplicate matched control for target {key}")
        controls_by_target[key] = row
    validation_targets = {str(row["id"]) for row in target_rows}
    if not set(controls_by_target) <= validation_targets or len(controls_by_target) != 798:
        raise ValueError("Expected 798 controls matched to unique sample targets")

    extraction_records = {row["control_id"]: row for row in extraction_summary["records"]}
    control_coverage = {}
    control_header_filter = {}
    input_hashes = {str(path): sha256(path) for path in
                    (target_manifest_path, control_manifest_path, status_path,
                     extraction_summary_path, mask_summary_path, coverage_path, records_path)}
    for target_id, control in controls_by_target.items():
        control_id = control["id"]
        baseline = args.extraction / control_id / "baseline.npz"
        metadata_path = args.extraction / control_id / "metadata.json"
        digest = sha256(baseline)
        record = extraction_records.get(control_id)
        metadata = read_json(metadata_path)
        if record is None or record.get("baseline_sha256") != digest or metadata.get("baseline_sha256") != digest:
            raise ValueError(f"Control baseline hash mismatch for {control_id}")
        if str(metadata.get("target_cluster_id")) != target_id:
            raise ValueError(f"Control-to-target link mismatch for {control_id}")
        input_hashes[str(baseline)] = digest
        input_hashes[str(metadata_path)] = sha256(metadata_path)
        control_header_filter[target_id] = metadata["header_filter"]
        with np.load(baseline) as archive:
            control_coverage[target_id] = archive["coverage"]

    by_target = {str(row["id"]): row for row in target_rows}
    support_rows = []
    stratified_rows = []
    selected_lists = {}
    for fi, factor in enumerate(FACTORS):
        for threshold in THRESHOLDS:
            for scope, bins in SCOPES.items():
                eligible_ids = []
                center_ids = []
                filter_ids = []
                for target_id, control in controls_by_target.items():
                    c_cov = by_target_coverage[target_id][fi, bins]
                    r_cov = control_coverage[target_id][fi, bins]
                    if not (np.isfinite(c_cov).all() and np.isfinite(r_cov).all()):
                        continue
                    if not (np.all(c_cov >= threshold) and np.all(r_cov >= threshold)):
                        continue
                    eligible_ids.append(target_id)
                    meta = cluster_records[target_id]
                    if meta["unique"] and not meta["base_masked"] and not meta["center_masked"][str(factor)]:
                        center_ids.append(target_id)
                    if by_target[target_id]["coadd_header_filter"] == control_header_filter[target_id]:
                        filter_ids.append(target_id)
                selected_targets = [by_target[ident] for ident in eligible_ids]
                selected_controls = [controls_by_target[ident] for ident in eligible_ids]
                p = len(bins)
                blocks = joint_block_count(selected_targets, selected_controls) if eligible_ids else 0
                key = f"k{factor}_c{int(threshold*100)}_r{scope.replace('-', '_')}"
                selected_lists[key] = eligible_ids
                support_rows.append(dict(
                    name=key, kron_factor=factor, coverage_floor=threshold, radial_scope_kpc=scope,
                    n_bins=p, n_matched_pairs=len(eligible_ids),
                    n_center_unmasked_auto_candidates=len(center_ids), n_matching_filter=len(filter_ids),
                    n_joint_5deg_blocks=blocks,
                    area_sample_size_screen=(len(eligible_ids) >= 5*p and blocks >= 2*p+1),
                    flux_covariance_available=False, science_measurement=False,
                ))
                for field, getter in (
                    ("redshift_bin", lambda r: r.get("redshift_bin", "unknown")),
                    ("richness_bin", lambda r: r.get("richness_bin", "unknown")),
                    ("target_filter", lambda r: r.get("coadd_header_filter", "unknown")),
                ):
                    levels = sorted({getter(by_target[ident]) for ident in eligible_ids})
                    for level in levels:
                        stratified_rows.append(dict(
                            name=key, factor=factor, coverage_floor=threshold, radial_scope_kpc=scope,
                            stratum=field, value=level,
                            n_pairs=sum(getter(by_target[ident]) == level for ident in eligible_ids),
                            science_measurement=False,
                        ))

    args.work.mkdir(parents=True, exist_ok=False)
    support_out = args.work / "support_counts.tsv"
    strata_out = args.work / "strata_counts.tsv"
    selection_out = args.work / "selected_pair_ids.json"
    write_tsv(support_out, support_rows)
    write_tsv(strata_out, stratified_rows)
    selection_out.write_text(json.dumps(selected_lists, indent=2, sort_keys=True) + "\n")
    summary = dict(
        status="complete", n_targets=len(target_rows), n_matched_controls=len(controls_by_target),
        n_unmatched_targets=len(target_rows) - len(controls_by_target),
        thresholds=list(THRESHOLDS), radial_scopes_kpc=SCOPES,
        n_support_configurations=len(support_rows),
        purpose="Paired usable-area support only; no cluster profiles or flux covariance computed.",
        limitations=[
            "Automatic center detections are positional candidates only and have no full-sample manual association review.",
            "Coverage floors and Kron factors are sensitivity cases; this report does not select a science cut.",
            "A pair count and sky-block count do not validate flux covariance or source-mask completeness.",
            "The existing 798 controls have assigned target angular scales, not physical control redshifts.",
        ],
        science_measurement=False, common_pair_covariance_approved=False,
        source_sha256=input_hashes,
        output_sha256={str(path): sha256(path) for path in (support_out, strata_out, selection_out)},
        script_sha256={str(Path(__file__).resolve()): sha256(Path(__file__).resolve())},
    )
    (args.work / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Full paired area support complete: {len(controls_by_target)} pairs; "
          f"{len(support_rows)} area configurations; flux covariance not computed", flush=True)


if __name__ == "__main__":
    main()
