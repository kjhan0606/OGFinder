#!/usr/bin/env python3
"""Describe sky, redshift, and richness selection in the S16A count-only audit.

This is a local-only analysis. It reads the approved aggregate count table,
the pinned S23B cluster manifest, and hash-pinned IDs from the existing fixed
52-lens ICL/WL sample. It makes no CAS/network request and reads no source-row,
image-pixel, shape-value, or photo-z data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
DATA_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN_DIR = DATA_ROOT / "wl_s23b_s16a_shape_count_plan_2026-11-11_01"
DEFAULT_OUTPUT = DATA_ROOT / "s23b_s16a_shape_coverage_selection_audit_2026-11-11_04"
S19A_COUNTS = DATA_ROOT / "wl_s23b_shape_count_audit_2026-11-07_01" / "per_center_counts.tsv"
S19A_SUMMARY = DATA_ROOT / "wl_s23b_shape_count_audit_2026-11-07_01" / "summary.json"
S19A_COUNT_SUMMARY = DATA_ROOT / "wl_s23b_shape_count_plan_2026-11-07_01" / "results" / "count_summary.json"
FIXED52_WL = DATA_ROOT / "wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_01" / "per_lens_set_profile.tsv"
ICL_PAIRS = DATA_ROOT / "s23b_icl_satellite_relation_2026-11-07_03" / "paired_icl_satellite_profiles.tsv"

EXPECTED_PLAN_SHA256 = "c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0"
EXPECTED_TARGET_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_COUNT_SHA256 = "0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12"
EXPECTED_COUNT_SUMMARY_SHA256 = "0cf90929c2a69eda0e62f3a9350359fd430b37bf869b2ead53246129ac51bb8b"
EXPECTED_S19A_COUNTS_SHA256 = "d650564b4870d4af93ce04b9ebd82475a99366aa68c2dfe236e7f070826c0f7f"
EXPECTED_S19A_PLAN_SHA256 = "b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984"
EXPECTED_S19A_AUDIT_SUMMARY_SHA256 = "e8b25a5d6b939f0b9602fd0fdd622862f5f5b3d898becf526cf384f0cc84db84"
EXPECTED_S19A_COUNT_SUMMARY_SHA256 = "748fea398f144165a25210071a5f6ab1866792df3be231ccfdb75282829935b8"
EXPECTED_FIXED52_SHA256 = "5bc599d97123264058011313b0ad35c3def49b08bba98a016cc5edf6ff03906c"
EXPECTED_ICL_PAIRS_SHA256 = "98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def quantiles(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return {"n": 0, "min": None, "q16": None, "median": None, "q84": None, "max": None}
    return {
        "n": int(array.size),
        "min": float(np.min(array)),
        "q16": float(np.quantile(array, 0.16)),
        "median": float(np.median(array)),
        "q84": float(np.quantile(array, 0.84)),
        "max": float(np.max(array)),
    }


def read_counts(path: Path) -> dict[str, int]:
    rows = read_tsv(path)
    if not rows or list(rows[0]) != ["lens_id", "n_shape_sources"]:
        raise ValueError("count table must contain exactly lens_id and n_shape_sources")
    counts: dict[str, int] = {}
    for row in rows:
        lens_id = row["lens_id"]
        try:
            count = int(row["n_shape_sources"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"non-integer count for {lens_id}") from error
        if not lens_id or count < 0 or lens_id in counts:
            raise ValueError(f"invalid or duplicate count row for {lens_id}")
        counts[lens_id] = count
    return counts


def read_batch_csv(path: Path) -> dict[str, int]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        header[0] = header[0].lstrip("# ").strip()
        if header != ["lens_id", "n_shape_sources"]:
            raise ValueError(f"unexpected count columns in {path.name}: {header}")
        counts: dict[str, int] = {}
        for values in reader:
            if len(values) != 2:
                raise ValueError(f"malformed count row in {path.name}")
            lens_id, raw_count = values
            try:
                count = int(raw_count)
            except ValueError as error:
                raise ValueError(f"non-integer count for {lens_id} in {path.name}") from error
            if not lens_id or count < 0 or lens_id in counts:
                raise ValueError(f"invalid or duplicate count row for {lens_id} in {path.name}")
            counts[lens_id] = count
        return counts


def read_verified_s19a_counts(targets: dict[str, dict]) -> tuple[dict[str, int], dict]:
    if sha256(S19A_COUNTS) != EXPECTED_S19A_COUNTS_SHA256:
        raise ValueError("approved S19A count table hash changed")
    if sha256(S19A_SUMMARY) != EXPECTED_S19A_AUDIT_SUMMARY_SHA256:
        raise ValueError("S19A local audit summary hash changed")
    if sha256(S19A_COUNT_SUMMARY) != EXPECTED_S19A_COUNT_SUMMARY_SHA256:
        raise ValueError("approved S19A query summary hash changed")
    summary = read_json(S19A_SUMMARY)
    if (summary.get("status") != "complete_exploratory_count_support_audit"
            or summary.get("count_only") is not True
            or summary.get("science_measurement") is not False
            or summary.get("source_rows_downloaded") is not False
            or summary.get("images_returned") is not False
            or summary.get("photoz_or_shape_values_returned") is not False
            or summary.get("inputs", {}).get("count_summary_sha256") != EXPECTED_S19A_COUNT_SUMMARY_SHA256
            or summary.get("inputs", {}).get("plan_sha256") != EXPECTED_S19A_PLAN_SHA256
            or summary.get("inputs", {}).get("source_target_manifest_sha256") != EXPECTED_TARGET_SHA256
            or summary.get("output_table_sha256") != EXPECTED_S19A_COUNTS_SHA256):
        raise ValueError("S19A count audit scope or provenance changed")
    records = read_tsv(S19A_COUNTS)
    required = ["lens_id", "catalog_id", "z_lens", "richness", "n_shape_sources"]
    if not records or list(records[0]) != required:
        raise ValueError("S19A count table columns changed")
    counts = {}
    for row in records:
        lens_id = str(row["lens_id"])
        try:
            count = int(row["n_shape_sources"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid S19A count for {lens_id}") from error
        if lens_id not in targets or lens_id in counts or count < 0:
            raise ValueError(f"unexpected or duplicate S19A count ID: {lens_id}")
        if (not np.isclose(float(row["z_lens"]), float(targets[lens_id]["z_cl"]), atol=5e-7, rtol=0)
                or not np.isclose(float(row["richness"]), float(targets[lens_id]["richness"]), atol=5e-4, rtol=0)):
            raise ValueError(f"S19A and S23B target metadata differ for {lens_id}")
        counts[lens_id] = count
    if (len(counts) != summary.get("n_centers")
            or len(counts) != 427
            or sum(counts.values()) != summary.get("sum_lens_cone_shape_counts")):
        raise ValueError("S19A count rows disagree with the approved audit summary")
    return counts, summary


def count_group(rows: list[dict], variable: str) -> list[dict]:
    ordered = sorted(rows, key=lambda row: (float(row[variable]), str(row["lens_id"])))
    groups = np.array_split(np.asarray(ordered, dtype=object), 5)
    output = []
    for index, group in enumerate(groups, start=1):
        selected = list(group)
        counts = [int(row["n_shape_sources"]) for row in selected]
        values = [float(row[variable]) for row in selected]
        output.append({
            "variable": variable,
            "quintile": index,
            "n_centers": len(selected),
            "value_min": min(values),
            "value_max": max(values),
            "positive_centers": sum(value > 0 for value in counts),
            "positive_fraction": sum(value > 0 for value in counts) / len(counts),
            "zero_centers": sum(value == 0 for value in counts),
            "median_aperture_count": float(np.median(counts)),
        })
    return output


def summarize_group(rows: list[dict]) -> dict:
    counts = [int(row["n_shape_sources"]) for row in rows]
    positive = [value for value in counts if value > 0]
    return {
        "n_centers": len(rows),
        "positive_centers": len(positive),
        "zero_centers": len(counts) - len(positive),
        "positive_fraction": len(positive) / len(counts) if counts else None,
        "aperture_count": quantiles(counts),
        "positive_count": quantiles(positive),
        "redshift_z_cl": quantiles([float(row["z_cl"]) for row in rows]),
        "camira_richness": quantiles([float(row["richness"]) for row in rows]),
    }


def verify_inputs(plan_dir: Path) -> tuple[dict, dict[str, int], list[dict], dict, dict]:
    # Reuse the guarded submitter's exact-plan and target-manifest validator.
    import sys
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    from submit_s23b_s16a_shape_counts import validate_plan

    plan, plan_sha = validate_plan(plan_dir)
    if plan_sha != EXPECTED_PLAN_SHA256:
        raise ValueError("count plan hash differs from the approved plan")
    manifest_path = Path(plan["source_target_manifest"])
    if sha256(manifest_path) != EXPECTED_TARGET_SHA256:
        raise ValueError("pinned target manifest hash changed")
    manifest = read_json(manifest_path)

    results = plan_dir / "results"
    count_path = results / "count_results.tsv"
    summary_path = results / "count_summary.json"
    status_path = results / "count_status.json"
    if sha256(count_path) != EXPECTED_COUNT_SHA256:
        raise ValueError("approved aggregate count result hash changed")
    if sha256(summary_path) != EXPECTED_COUNT_SUMMARY_SHA256:
        raise ValueError("approved count summary hash changed")
    summary = read_json(summary_path)
    status = read_json(status_path)
    if (summary.get("status") != "complete"
            or summary.get("query_plan_sha256") != EXPECTED_PLAN_SHA256
            or status.get("status") != "complete"
            or status.get("plan_sha256") != EXPECTED_PLAN_SHA256
            or status.get("credentials_saved") is not False
            or any(status.get(key) is not False for key in (
                "source_rows_downloaded", "images_returned", "photoz_returned", "science_measurement"))):
        raise ValueError("count run scope/status is inconsistent")
    counts = read_counts(count_path)
    expected_ids = {str(center["lens_id"]) for center in plan["centers"]}
    if len(counts) != 1043 or set(counts) != expected_ids:
        raise ValueError("aggregate count IDs do not match all approved centers")
    if sum(counts.values()) != int(summary["sum_lens_cone_shape_counts"]):
        raise ValueError("combined count sum differs from the saved summary")
    jobs = status.get("jobs", {})
    if len(jobs) != len(plan["batches"]) or any(item.get("status") != "complete" for item in jobs.values()):
        raise ValueError("not all approved count-query batches are complete")
    batch_counts: dict[str, int] = {}
    batch_size = int(plan["batch_size"])
    for index in range(len(plan["batches"])):
        name = f"s16a_shape_counts_{index:03d}"
        item = jobs[name]
        batch_path = results / f"{name}.csv"
        if sha256(batch_path) != item.get("sha256"):
            raise ValueError(f"batch result hash changed: {name}")
        result = read_batch_csv(batch_path)
        start = index * batch_size
        expected_batch_ids = {
            str(center["lens_id"]) for center in plan["centers"][start:start + batch_size]
        }
        if (len(result) != int(item.get("rows", -1))
                or set(result) != expected_batch_ids
                or sum(result.values()) != int(item.get("total_sources", -1))
                or sum(value == 0 for value in result.values()) != int(item.get("zero_source_centers", -1))):
            raise ValueError(f"batch row IDs or aggregate values changed: {name}")
        batch_counts.update(result)
    if batch_counts != counts:
        raise ValueError("combined table differs from the 42 verified batch results")

    target_map = {str(row["id"]): row for row in manifest["clusters"]}
    s19_counts, s19_summary = read_verified_s19a_counts(target_map)

    if sha256(FIXED52_WL) != EXPECTED_FIXED52_SHA256:
        raise ValueError("fixed-52 WL profile hash changed")
    fixed_rows = read_tsv(FIXED52_WL)
    fixed_ids = {str(row["matched_lens_id"]) for row in fixed_rows}
    if len(fixed_ids) != 52:
        raise ValueError("fixed WL profile does not contain exactly 52 lenses")

    if sha256(ICL_PAIRS) != EXPECTED_ICL_PAIRS_SHA256:
        raise ValueError("ICL/satellite pair-table hash changed")
    icl_rows = read_tsv(ICL_PAIRS)
    band_ids = {"i": set(), "i2": set()}
    for row in icl_rows:
        lens_id, band = str(row.get("target_id", "")), str(row.get("filter", ""))
        if lens_id in fixed_ids and band in band_ids:
            band_ids[band].add(lens_id)
    if (len(band_ids["i"]) != 25 or len(band_ids["i2"]) != 27
            or band_ids["i"] & band_ids["i2"]
            or band_ids["i"] | band_ids["i2"] != fixed_ids):
        raise ValueError("fixed-52 ICL filter membership differs from the pinned audit")
    return plan, counts, manifest["clusters"], summary, {
        "fixed_ids": fixed_ids,
        "band_ids": band_ids,
        "status_path": status_path,
        "count_path": count_path,
        "summary_path": summary_path,
        "manifest_path": manifest_path,
        "s19_counts": s19_counts,
        "s19_summary": s19_summary,
        "icl_rows": icl_rows,
    }


def write_table(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_selection(rows: list[dict], quintiles: list[dict], path: Path) -> None:
    zero = [row for row in rows if int(row["n_shape_sources"]) == 0]
    positive = [row for row in rows if int(row["n_shape_sources"]) > 0]
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.0), gridspec_kw={"width_ratios": [1.25, 1.0, 0.9]})
    fig.subplots_adjust(left=0.065, right=0.985, top=0.88, bottom=0.17, wspace=0.34)

    sky = axes[0]
    sky.scatter([r["ra_deg"] for r in zero], [r["dec_deg"] for r in zero],
                s=12, color="#c9cdd3", alpha=0.70, linewidths=0, label="0 catalog objects")
    sky.scatter([r["ra_deg"] for r in positive], [r["dec_deg"] for r in positive],
                s=19, color="#1769aa", alpha=0.9, linewidths=0, label=">0 catalog objects")
    sky.set_xlim(360, 0)
    dec = [float(row["dec_deg"]) for row in rows]
    sky.set_ylim(min(dec) - 2.0, max(dec) + 10.0)
    sky.set_xlabel("Right ascension (deg)")
    sky.set_ylabel("Declination (deg)")
    sky.set_title("S16A aperture-count coverage")
    sky.grid(alpha=0.18)
    sky.legend(loc="upper center", fontsize=8, framealpha=0.92, ncol=1)

    rz = axes[1]
    rz.scatter([r["z_cl"] for r in zero], [r["richness"] for r in zero],
               s=12, color="#c9cdd3", alpha=0.65, linewidths=0)
    rz.scatter([r["z_cl"] for r in positive], [r["richness"] for r in positive],
               s=19, color="#1769aa", alpha=0.85, linewidths=0)
    rz.set_xlabel("Cluster redshift")
    rz.set_ylabel("CAMIRA richness")
    rz.set_title("Coverage vs. cluster properties")
    rz.grid(alpha=0.18)

    frac = axes[2]
    for variable, label, color, marker in (
        ("richness", "Richness quintile", "#bc4b51", "o"),
        ("z_cl", "Redshift quintile", "#2a7f62", "s"),
    ):
        groups = [row for row in quintiles if row["variable"] == variable]
        frac.plot([row["quintile"] for row in groups],
                  [100 * row["positive_fraction"] for row in groups],
                  color=color, marker=marker, linewidth=1.6, markersize=4.5, label=label)
    frac.set_xticks(range(1, 6), ["Q1", "Q2", "Q3", "Q4", "Q5"])
    frac.set_ylim(0, 40)
    frac.set_xlabel("Within-sample quintile")
    frac.set_ylabel("Centers with >0 count (%)")
    frac.set_title("Positive-count fraction")
    frac.grid(alpha=0.18)
    frac.legend(loc="upper right", fontsize=8, framealpha=0.92)

    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)


def analyze(plan_dir: Path, output: Path) -> dict:
    plan_dir, output = plan_dir.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output}")
    plan, counts, manifest_rows, count_summary, refs = verify_inputs(plan_dir)
    if sha256(refs["manifest_path"]) != EXPECTED_TARGET_SHA256:
        raise ValueError("target manifest changed after input validation")

    meta = {str(row["id"]): row for row in manifest_rows}
    if set(meta) != set(counts) or any(row.get("sample_role") != "target" for row in manifest_rows):
        raise ValueError("target metadata does not match the 1,043 count IDs")
    rows = []
    for lens_id in sorted(counts):
        source = meta[lens_id]
        bands = [band for band, ids in refs["band_ids"].items() if lens_id in ids]
        s19_count = refs["s19_counts"].get(lens_id)
        rows.append({
            "lens_id": lens_id,
            "ra_deg": float(source["ra_deg"]),
            "dec_deg": float(source["dec_deg"]),
            "z_cl": float(source["z_cl"]),
            "richness": float(source["richness"]),
            "n_shape_sources": counts[lens_id],
            "s16a_count_positive": counts[lens_id] > 0,
            "s19a_count_status": "queried" if s19_count is not None else "outside_427_center_audit",
            "n_shape_sources_s19a": "" if s19_count is None else s19_count,
            "fixed52_s19a_lens": lens_id in refs["fixed_ids"],
            "icl_filter": bands[0] if bands else "",
        })

    release_overlap = []
    overlap_counts = {
        "both_positive": 0,
        "s16a_only_positive": 0,
        "s19a_only_positive": 0,
        "both_zero": 0,
        "both_at_least_1000": 0,
        "s16a_at_least_1000_only": 0,
        "s19a_at_least_1000_only": 0,
    }
    for lens_id in sorted(refs["s19_counts"]):
        s16_count, s19_count = counts[lens_id], refs["s19_counts"][lens_id]
        s16_positive, s19_positive = s16_count > 0, s19_count > 0
        s16_1000, s19_1000 = s16_count >= 1000, s19_count >= 1000
        if s16_positive and s19_positive:
            overlap_counts["both_positive"] += 1
        elif s16_positive:
            overlap_counts["s16a_only_positive"] += 1
        elif s19_positive:
            overlap_counts["s19a_only_positive"] += 1
        else:
            overlap_counts["both_zero"] += 1
        if s16_1000 and s19_1000:
            overlap_counts["both_at_least_1000"] += 1
        elif s16_1000:
            overlap_counts["s16a_at_least_1000_only"] += 1
        elif s19_1000:
            overlap_counts["s19a_at_least_1000_only"] += 1
        meta_row = meta[lens_id]
        release_overlap.append({
            "lens_id": lens_id,
            "z_cl": float(meta_row["z_cl"]),
            "richness": float(meta_row["richness"]),
            "n_shape_sources_s16a": s16_count,
            "s16a_positive": s16_positive,
            "n_shape_sources_s19a": s19_count,
            "s19a_positive": s19_positive,
            "s16a_ge_1000": s16_1000,
            "s19a_ge_1000": s19_1000,
        })
    release_crosscheck = {
        "interpretation": "same-target aggregate count comparison only; rows are not combined",
        "n_s19a_centers_audited": len(refs["s19_counts"]),
        "n_common_target_centers": len(release_overlap),
        "s19a_positive_centers": sum(value > 0 for value in refs["s19_counts"].values()),
        "s19a_centers_ge_1000": sum(value >= 1000 for value in refs["s19_counts"].values()),
        "s16a_positive_in_common_427": overlap_counts["both_positive"] + overlap_counts["s16a_only_positive"],
        "s16a_centers_ge_1000_in_common_427": overlap_counts["both_at_least_1000"] + overlap_counts["s16a_at_least_1000_only"],
        "s16a_positive_outside_s19a_427_audit": sum(
            value > 0 for lens_id, value in counts.items() if lens_id not in refs["s19_counts"]
        ),
        "s16a_ge_1000_outside_s19a_427_audit": sum(
            value >= 1000 for lens_id, value in counts.items() if lens_id not in refs["s19_counts"]
        ),
        "contingency_427": overlap_counts,
    }

    expansion_candidates = {
        row["lens_id"] for row in rows
        if row["s19a_count_status"] == "outside_427_center_audit"
        and int(row["n_shape_sources"]) >= 1000
    }
    control_rows_by_target = {}
    for pair in refs["icl_rows"]:
        lens_id = str(pair.get("target_id", ""))
        if lens_id not in expansion_candidates:
            continue
        if lens_id in control_rows_by_target:
            raise ValueError(f"multiple saved ICL control pairs for expansion candidate {lens_id}")
        source = meta[lens_id]
        if (not np.isclose(float(pair["z_cl"]), float(source["z_cl"]), atol=5e-7, rtol=0)
                or not np.isclose(float(pair["richness"]), float(source["richness"]), atol=5e-4, rtol=0)
                or not np.isclose(float(pair["ra_deg"]), float(source["ra_deg"]), atol=5e-7, rtol=0)
                or not np.isclose(float(pair["dec_deg"]), float(source["dec_deg"]), atol=5e-7, rtol=0)):
            raise ValueError(f"saved ICL pair target does not match manifest for {lens_id}")
        control_rows_by_target[lens_id] = pair
    control_ids = [str(pair.get("control_id", "")) for pair in control_rows_by_target.values()]
    if any(not control_id for control_id in control_ids) or len(set(control_ids)) != len(control_ids):
        raise ValueError("ICL control IDs are missing or reused among expansion candidates")
    expansion_controls = []
    for lens_id in sorted(control_rows_by_target):
        pair = control_rows_by_target[lens_id]
        source = meta[lens_id]
        expansion_controls.append({
            "target_id": lens_id,
            "control_id": str(pair["control_id"]),
            "icl_filter": str(pair["filter"]),
            "target_ra_deg": float(source["ra_deg"]),
            "target_dec_deg": float(source["dec_deg"]),
            "control_ra_deg": float(pair["control_ra_deg"]),
            "control_dec_deg": float(pair["control_dec_deg"]),
            "z_cl": float(source["z_cl"]),
            "richness": float(source["richness"]),
            "s16a_target_count_3pmpc": int(counts[lens_id]),
            "delta_ebv_sfd_target_minus_control": float(pair["delta_ebv_sfd_target_minus_control"]),
            "control_s16a_count_status": "not_queried",
        })
    filter_counts = {}
    for row in expansion_controls:
        filter_counts[row["icl_filter"]] = filter_counts.get(row["icl_filter"], 0) + 1
    expansion_control_summary = {
        "selection": "outside the 427-center S19A count audit and S16A raw count >= 1000",
        "raw_count_threshold_is_final_wl_eligibility": False,
        "n_candidate_targets": len(expansion_candidates),
        "n_targets_with_saved_icl_control_pair": len(expansion_controls),
        "n_targets_without_saved_icl_control_pair": len(expansion_candidates) - len(expansion_controls),
        "n_unique_control_centers": len(set(control_ids)),
        "saved_icl_filter_rows": filter_counts,
        "abs_delta_ebv_median": float(np.median([
            abs(row["delta_ebv_sfd_target_minus_control"]) for row in expansion_controls
        ])) if expansion_controls else None,
        "abs_delta_ebv_max": max((
            abs(row["delta_ebv_sfd_target_minus_control"]) for row in expansion_controls
        ), default=None),
        "control_counts_queried": False,
        "cas_query_submitted": False,
    }

    positive = [row for row in rows if row["n_shape_sources"] > 0]
    zero = [row for row in rows if row["n_shape_sources"] == 0]
    categories = {
        "all_targets": rows,
        "positive_count": positive,
        "zero_count": zero,
        "fixed52_positive_count": [row for row in rows if row["fixed52_s19a_lens"] and row["n_shape_sources"] > 0],
        "fixed52_zero_count": [row for row in rows if row["fixed52_s19a_lens"] and row["n_shape_sources"] == 0],
    }
    quintiles = count_group(rows, "richness") + count_group(rows, "z_cl")
    band_summary = []
    for band in ("i", "i2"):
        selected = [row for row in rows if row["icl_filter"] == band]
        n_positive = sum(int(row["n_shape_sources"]) > 0 for row in selected)
        positive_counts = [int(row["n_shape_sources"]) for row in selected if int(row["n_shape_sources"]) > 0]
        band_summary.append({
            "icl_filter": band,
            "n_fixed_lenses": len(selected),
            "positive_count_lenses": n_positive,
            "zero_count_lenses": len(selected) - n_positive,
            "positive_fraction": n_positive / len(selected) if selected else None,
            "median_positive_count": float(np.median(positive_counts)) if positive_counts else None,
        })

    count_values = [int(row["n_shape_sources"]) for row in rows]
    rho_z = float(spearmanr(count_values, [row["z_cl"] for row in rows]).statistic)
    rho_richness = float(spearmanr(count_values, [row["richness"] for row in rows]).statistic)
    output.mkdir(parents=True, exist_ok=False)

    per_cluster_path = output / "per_cluster_coverage.tsv"
    per_fields = ["lens_id", "ra_deg", "dec_deg", "z_cl", "richness", "n_shape_sources",
                  "s16a_count_positive", "s19a_count_status", "n_shape_sources_s19a",
                  "fixed52_s19a_lens", "icl_filter"]
    write_table(per_cluster_path, rows, per_fields)
    release_overlap_path = output / "s16a_s19a_count_overlap_427.tsv"
    write_table(release_overlap_path, release_overlap, [
        "lens_id", "z_cl", "richness", "n_shape_sources_s16a", "s16a_positive",
        "n_shape_sources_s19a", "s19a_positive", "s16a_ge_1000", "s19a_ge_1000",
    ])
    expansion_controls_path = output / "s16a_ge1000_outside_s19a_icl_control_centers.tsv"
    write_table(expansion_controls_path, expansion_controls, [
        "target_id", "control_id", "icl_filter", "target_ra_deg", "target_dec_deg",
        "control_ra_deg", "control_dec_deg", "z_cl", "richness",
        "s16a_target_count_3pmpc", "delta_ebv_sfd_target_minus_control",
        "control_s16a_count_status",
    ])
    quintile_path = output / "positive_fraction_quintiles.tsv"
    write_table(quintile_path, quintiles, ["variable", "quintile", "n_centers", "value_min", "value_max",
                                          "positive_centers", "positive_fraction", "zero_centers", "median_aperture_count"])
    band_path = output / "fixed52_filter_coverage.tsv"
    write_table(band_path, band_summary, ["icl_filter", "n_fixed_lenses", "positive_count_lenses",
                                         "zero_count_lenses", "positive_fraction", "median_positive_count"])
    figure_path = output / "s16a_count_coverage_selection.png"
    plot_selection(rows, quintiles, figure_path)

    files = [per_cluster_path, release_overlap_path, expansion_controls_path,
             quintile_path, band_path, figure_path]
    summary = {
        "schema": "astrafek-s23b-s16a-count-selection-audit-v1",
        "status": "complete_local_descriptive_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "source_rows_downloaded": False,
        "images_read": False,
        "photoz_or_shape_values_read": False,
        "n_centers": len(rows),
        "positive_centers": len(positive),
        "zero_centers": len(zero),
        "positive_fraction": len(positive) / len(rows),
        "sum_lens_cone_shape_counts": sum(count_values),
        "positive_count_distribution": quantiles([int(row["n_shape_sources"]) for row in positive]),
        "groups": {name: summarize_group(group) for name, group in categories.items()},
        "fixed52_by_icl_filter": band_summary,
        "cross_release_count_comparison": release_crosscheck,
        "expansion_candidate_icl_control_pairs": expansion_control_summary,
        "spearman_descriptive_only": {
            "rho_count_vs_z_cl": rho_z,
            "rho_count_vs_richness": rho_richness,
            "p_values_used": False,
        },
        "quintile_method": "rank-sort all 1,043 centers by the named variable and split into five near-equal groups; ties ordered by lens ID",
        "quintiles": quintiles,
        "interpretation": [
            "Counts are distinct S16A shape-catalog object counts within each 3-proper-Mpc lens cone before redshift, quality, radial, mask, weight, and estimator cuts.",
            "The summed count is a sum of lens-cone counts; an object may appear in multiple cones, so this is not a unique-source total.",
            "The positive/zero split is spatially selected and cannot be treated as a random S23B subsample; the displayed property trends are descriptive and do not separate sky footprint from cluster properties.",
            "S16A/S19A cross-release rows compare only aggregate aperture counts for the same targets; no source rows are concatenated or treated as a single catalog.",
            "No source-row transfer, lensing stack, S16A/S19A merge, or science measurement is performed by this audit.",
        ],
        "inputs": {
            "analysis_script_sha256": sha256(Path(__file__).resolve()),
            "query_plan_sha256": EXPECTED_PLAN_SHA256,
            "count_results_sha256": EXPECTED_COUNT_SHA256,
            "count_summary_sha256": EXPECTED_COUNT_SUMMARY_SHA256,
            "s19a_count_results_sha256": EXPECTED_S19A_COUNTS_SHA256,
            "s19a_count_plan_sha256": EXPECTED_S19A_PLAN_SHA256,
            "s19a_audit_summary_sha256": EXPECTED_S19A_AUDIT_SUMMARY_SHA256,
            "s19a_count_summary_sha256": EXPECTED_S19A_COUNT_SUMMARY_SHA256,
            "target_manifest_sha256": EXPECTED_TARGET_SHA256,
            "fixed52_wl_profile_sha256": EXPECTED_FIXED52_SHA256,
            "icl_satellite_pair_table_sha256": EXPECTED_ICL_PAIRS_SHA256,
        },
        "output_directory": str(output),
        "outputs_sha256": {path.name: sha256(path) for path in files},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = analyze(args.plan_dir, args.output)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "n_centers": summary["n_centers"],
        "positive_centers": summary["positive_centers"],
        "zero_centers": summary["zero_centers"],
        "fixed52_positive": summary["groups"]["fixed52_positive_count"]["n_centers"],
        "s16a_only_positive_in_s19a_427": summary["cross_release_count_comparison"]["contingency_427"]["s16a_only_positive"],
        "s16a_positive_outside_s19a_427": summary["cross_release_count_comparison"]["s16a_positive_outside_s19a_427_audit"],
        "rho_count_vs_richness": summary["spearman_descriptive_only"]["rho_count_vs_richness"],
        "network_used": False,
        "cas_query_submitted": False,
        "output": summary["output_directory"],
    }, indent=2))


if __name__ == "__main__":
    main()
