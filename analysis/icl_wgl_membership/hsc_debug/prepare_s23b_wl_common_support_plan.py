#!/usr/bin/env python3
"""Audit local S23B weak-lensing support and prepare, but never submit, plans.

This offline planner joins the approved S19A count-only results to the saved
S23B target and random-point candidate manifests. It selects provisional
common-support lens/control candidates and writes exact count-only and shape-row
SQL plans. It makes no network request and saves no credentials.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

from astropy.cosmology import FlatLambdaCDM

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prepare_wl_random_shape_count_plan import normalize_centers, shape_count_query
from prepare_wl_random_shape_rows import shape_rows_query
from prepare_wl_shape_rows import OUTPUT_COLUMNS as LENS_ROW_COLUMNS
from prepare_wl_shape_rows import lens_query
from submit_wl_random_shape_counts import validate_plan as validate_random_count_plan
from submit_s23b_wl_shape_counts import read_count_csv, validate_plan

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
COUNT_PLAN = DATA / "wl_s23b_shape_count_plan_2026-11-07_01"
COUNT_AUDIT = DATA / "wl_s23b_shape_count_audit_2026-11-07_01"
TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
CONTROL_AUDIT = DATA / "icl_s23b_control_candidate_audit_2026-11-07_01"
SOURCE_UNION_AUDIT = DATA / "wl_s23b_existing_source_union_2026-11-07_01"
OLD_LENS_EXPORT = DATA / "wl_shape_rows_export_2026-10-18_03"
OLD_RANDOM_EXPORT = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01" / "results"
OLD_RANDOM_SUPPORT = DATA / "wl_random_shape_count_support_audit_2026-11-03_01"
CAMIRA_CATALOG = Path(
    "/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/"
    "camira_s23b_wide_sm.tbl"
)
DEFAULT_OUTPUT = DATA / "wl_s23b_wl_common_support_plan_2026-11-10_01"

EXPECTED_COUNT_PLAN_SHA256 = "b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984"
EXPECTED_TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_CONTROL_CANDIDATE_AUDIT_PLAN_SHA256 = "143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714"
EXPECTED_CAMIRA_CATALOG_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_RANDOM_SUPPORT_PLAN_SHA256 = "18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625"

MIN_LENS_COUNT = 1000
CONTROLS_PER_LENS = 5
MIN_TARGET_CONTROL_SEPARATION_MPC = 6.0
BATCH_SIZE = 25
H0 = 70.0
OMEGA_M = 0.3
ROW_BYTES_LOW = 256
ROW_BYTES_HIGH = 480

LENS_COLUMNS = list(LENS_ROW_COLUMNS)
RANDOM_COLUMNS = [
    "random_center_id", "matched_lens_id", "random_ra_deg", "random_dec_deg", "lens_z",
    "object_id", "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_mask_raw", "source_zbin",
]
COUNT_COLUMNS = ["random_center_id", "matched_lens_id", "n_shape_sources"]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha(document: dict) -> str:
    payload = dict(document)
    payload.pop("query_plan_sha256", None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def write_json(path: Path, document: dict) -> None:
    path.write_text(json.dumps(document, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_verified_counts() -> tuple[dict, dict, dict[str, int]]:
    plan, plan_sha = validate_plan(COUNT_PLAN)
    if plan_sha != EXPECTED_COUNT_PLAN_SHA256:
        raise ValueError("S23B count plan differs from its approved hash")
    results = COUNT_PLAN / "results"
    checkpoint = read_json(results / "count_status.json")
    summary = read_json(results / "count_summary.json")
    if (checkpoint.get("status") != "complete" or summary.get("status") != "complete"
            or checkpoint.get("plan_sha256") != plan_sha
            or summary.get("plan_sha256") not in (None, plan_sha)):
        raise ValueError("approved 427-center count results are incomplete or mismatched")
    for key in ("credentials_saved", "source_rows_downloaded", "images_returned",
                "photoz_returned", "science_measurement"):
        if checkpoint.get(key) is not False:
            raise ValueError(f"count checkpoint has an unexpected scope flag: {key}")
    if (any(summary.get(key) is not False for key in (
            "source_rows_downloaded", "images_returned", "science_measurement"))
            or any(key in summary and summary.get(key) is not False for key in (
                "photoz_returned", "photoz_or_shape_values_returned"))):
        raise ValueError("count summary contains data outside its declared scope")

    centers = plan["centers"]
    batch_size = int(plan["batch_size"])
    jobs = checkpoint.get("jobs", {})
    counts: dict[str, int] = {}
    for index, start in enumerate(range(0, len(centers), batch_size)):
        name = f"shape_counts_{index:03d}"
        path = results / f"{name}.csv"
        item = jobs.get(name, {})
        batch_centers = centers[start:start + batch_size]
        ids = [str(row["lens_id"]) for row in batch_centers]
        if (item.get("status") != "complete" or not path.is_file()
                or sha(path) != item.get("sha256")):
            raise ValueError(f"missing or changed approved count batch: {name}")
        rows = read_count_csv(path, ids)
        if item.get("rows") != len(rows):
            raise ValueError(f"count receipt row total changed: {name}")
        counts.update({row["lens_id"]: row["n_shape_sources"] for row in rows})
    if (len(counts) != 427 or sum(counts.values())
            != int(summary["sum_lens_cone_shape_counts"])):
        raise ValueError("verified count rows disagree with the saved summary")
    return plan, summary, counts


def select_unique_controls(groups: dict[str, list[dict]], quota: int) -> dict[tuple[str, int], dict]:
    """Assign distinct sky positions to control slots with deterministic matching."""
    eligible = {lens_id: options for lens_id, options in groups.items() if len(options) >= quota}
    if not eligible:
        raise ValueError("no target has enough matched control candidates")
    slots = [(lens_id, index) for lens_id in eligible for index in range(quota)]
    slots.sort(key=lambda slot: (len(eligible[slot[0]]), slot[0], slot[1]))
    owner: dict[tuple[float, float], tuple[str, int]] = {}
    assignment: dict[tuple[str, int], dict] = {}

    def visit(slot: tuple[str, int], seen: set[tuple[float, float]]) -> bool:
        lens_id = slot[0]
        for candidate in eligible[lens_id]:
            key = (round(float(candidate["ra_deg"]), 12),
                   round(float(candidate["dec_deg"]), 12))
            if key in seen:
                continue
            seen.add(key)
            previous = owner.get(key)
            if previous is None or visit(previous, seen):
                owner[key] = slot
                assignment[slot] = candidate
                return True
        return False

    for slot in slots:
        if not visit(slot, set()):
            raise ValueError(f"could not assign {quota} globally unique controls per lens: {slot[0]}")
    if len(assignment) != len(slots) or len(owner) != len(slots):
        raise ValueError("control assignment is not globally one-to-one")
    return assignment


def load_camira_catalog(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            values = line.strip().split()
            if not values or line.lstrip().startswith("#"):
                continue
            if len(values) < 6:
                raise ValueError("CAMIRA cluster catalog row has fewer than six columns")
            rows.append({
                "catalog_id": values[0], "ra_deg": float(values[2]),
                "dec_deg": float(values[3]), "z": float(values[4]),
            })
    return rows


def angular_separation_arcmin(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    hav = (math.sin((dec2 - dec1) / 2) ** 2
           + math.cos(dec1) * math.cos(dec2) * math.sin((ra2 - ra1) / 2) ** 2)
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, max(0.0, hav))))) * 60


def count_existing_control_contamination(control_manifest: dict, targets: dict,
                                         count_ids: set[str], catalog: list[dict]) -> dict:
    cosmo = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Tcmb0=0)
    by_target = {str(row["target_cluster_id"]).replace(":", "_"): row
                 for row in control_manifest.get("clusters", [])}
    distances = []
    target_nearby_counts = {radius: 0 for radius in (3.0, 6.0)}
    allz_nearby_counts = {radius: 0 for radius in (3.0, 6.0)}
    dz_nearby_counts = {radius: 0 for radius in (3.0, 6.0)}
    allz_pair_counts = {radius: 0 for radius in (3.0, 6.0)}
    dz_pair_counts = {radius: 0 for radius in (3.0, 6.0)}
    for lens_id in sorted(count_ids):
        control = by_target.get(lens_id)
        if control is None:
            continue
        target = targets[lens_id]
        redshift = float(target["z"])
        da = float(cosmo.angular_diameter_distance(redshift).value)
        sep_arcmin = angular_separation_arcmin(
            float(target["ra_deg"]), float(target["dec_deg"]),
            float(control["ra_deg"]), float(control["dec_deg"]),
        )
        distance_mpc = sep_arcmin * da * math.pi / 10800.0
        distances.append(distance_mpc)
        for radius in target_nearby_counts:
            target_nearby_counts[radius] += distance_mpc < radius

        center_ra, center_dec = float(control["ra_deg"]), float(control["dec_deg"])
        center_dec_rad = math.radians(center_dec)
        control_allz = {radius: False for radius in allz_nearby_counts}
        control_dz = {radius: False for radius in dz_nearby_counts}
        for cluster in catalog:
            if str(cluster["catalog_id"]) == str(target["catalog_id"]):
                continue
            dec_delta = float(cluster["dec_deg"]) - center_dec
            if abs(dec_delta) > 7.0 / da * 180.0 / math.pi:
                continue
            ra_delta = (float(cluster["ra_deg"]) - center_ra + 180.0) % 360.0 - 180.0
            if abs(ra_delta * math.cos(center_dec_rad)) > 7.0 / da * 180.0 / math.pi:
                continue
            separation = angular_separation_arcmin(
                center_ra, center_dec, float(cluster["ra_deg"]), float(cluster["dec_deg"]),
            )
            projected_mpc = separation * da * math.pi / 10800.0
            same_slice = abs(float(cluster["z"]) - redshift) <= 0.05
            for radius in allz_nearby_counts:
                if projected_mpc < radius:
                    allz_pair_counts[radius] += 1
                    control_allz[radius] = True
                    if same_slice:
                        dz_pair_counts[radius] += 1
                        control_dz[radius] = True
        for radius in allz_nearby_counts:
            allz_nearby_counts[radius] += control_allz[radius]
            dz_nearby_counts[radius] += control_dz[radius]

    n_controls = len(distances)
    return {
        "n_controls_with_matched_lens": n_controls,
        "matched_target_control_separation_mpc": {
            "p10": sorted(distances)[int(0.1 * (len(distances) - 1))] if distances else None,
            "median": statistics.median(distances) if distances else None,
            "p90": sorted(distances)[int(0.9 * (len(distances) - 1))] if distances else None,
            "below_3_mpc": target_nearby_counts[3.0],
            "below_6_mpc": target_nearby_counts[6.0],
        },
        "catalog_clusters": len(catalog),
        "other_camira_centers_near_existing_controls": {
            "projected_radius_mpc_at_matched_lens_redshift": {
                str(radius): {
                    "controls_with_at_least_one_all_redshifts": allz_nearby_counts[radius],
                    "control_cluster_pairs_all_redshifts": allz_pair_counts[radius],
                    "controls_with_at_least_one_abs_delta_z_le_0_05": dz_nearby_counts[radius],
                    "control_cluster_pairs_abs_delta_z_le_0_05": dz_pair_counts[radius],
                }
                for radius in (3.0, 6.0)
            },
        },
        "catalog_used_for_control_contamination_counts": True,
        "note": "Projected separations use each matched target redshift. This is a CAMIRA-center proximity diagnostic, not a full projected-mass or shear contamination model.",
    }


def write_assessment_markdown(path: Path, assessment: dict, row_plan: dict,
                              count_plan: dict, count_plan_sha256: str,
                              row_plan_file_sha256: str) -> None:
    support = assessment["common_support_candidate"]
    parent = assessment["parent_lens_support"]
    estimates = row_plan["estimated_uncompressed_csv_bytes"]
    pair_proxy = row_plan["combined_source_pair_rows_planning_proxies"]
    random_proxy = row_plan["random_source_pair_rows_planning_proxies"]
    old_controls = assessment["existing_selected_icl_control_diagnostic"]
    old_sep = old_controls["matched_target_control_separation_mpc"]
    old_camira = old_controls["other_camira_centers_near_existing_controls"][
        "projected_radius_mpc_at_matched_lens_redshift"]
    fmt_gib = lambda value: f"{int(value) / 1024**3:.2f} GiB"
    lines = [
        "# S23B weak-lensing common-support plan",
        "",
        "Date: 2026-11-10. This is an offline candidate-support assessment and plan, not a lensing result.",
        "",
        "## Sample support",
        "",
        f"The approved 3-proper-Mpc count-only audit covers {parent['n_lenses']} S23B targets. "
        f"{parent['zero_raw_counts']} have zero released S19A shape-table entries; "
        f"{parent['at_least_1000']} have at least 1,000 before this project's source cuts. "
        "This count is only a raw availability screen.",
        "",
        f"The provisional sample contains {support['n_lenses']} lenses and "
        f"{support['n_controls']} globally distinct HSC random-table control positions "
        f"(five per lens). Their exact lens-side candidate-pair count is "
        f"{support['exact_lens_pair_rows']:,}; per-lens min/median/max are "
        f"{support['lens_count_min']:,}/{support['lens_count_median']:,.0f}/{support['lens_count_max']:,}. "
        f"The 80 lenses are distributed across {len(support['lens_counts_by_ra_30deg'])} occupied 30-degree RA bins.",
        "",
        f"The other {parent['at_least_1000'] - support['n_lenses']} lenses above the 1,000-count screen "
        "were excluded because the saved candidate pool did not provide five qualifying controls each. "
        "Every selected control has center flags and depth/seeing/SFD-EBV metadata matching, a prior "
        "CAMIRA 2-Mpc overlap screen of zero, and is more than 6 proper Mpc from its matched lens. "
        "The latter keeps paired 3-Mpc query cones from geometrically overlapping. The CAMIRA screen "
        "is incomplete and does not remove all environmental or large-scale shear contributions.",
        "",
        "The previously selected ICL controls are unsuitable as a weak-lensing control set: "
        f"{old_controls['n_controls_with_matched_lens']} mapped positions include "
        f"{old_sep['below_3_mpc']} within 3 Mpc and {old_sep['below_6_mpc']} within 6 Mpc of their matched lens. "
        f"Other S23B CAMIRA centers occur within 3 Mpc of {old_camira['3.0']['controls_with_at_least_one_all_redshifts']} "
        f"of those old controls at any redshift and "
        f"{old_camira['3.0']['controls_with_at_least_one_abs_delta_z_le_0_05']} at |dz| <= 0.05. "
        "These are center-proximity diagnostics, not a mass or shear contamination model.",
        "",
        f"The existing S19A lens/random source-row union has rows in only "
        f"{assessment['existing_source_union_selected_lens_rows_with_any_sources']} of the selected 80 lenses, "
        f"and none reaches 95% estimated query-cone area coverage; it cannot replace a new export.",
        "",
        "## Prepared actions (not submitted)",
        "",
        f"- Count-only control preflight: {count_plan['n_random_centers']} centers, "
        f"{count_plan['n_query_batches']} batches of {count_plan['batch_size']}, 3 output columns "
        "(`random_center_id`, `matched_lens_id`, `n_shape_sources`); no shape values, source coordinates, "
        "photo-z, or images. "
        f"SQL payload: {sum(batch['sql_bytes'] for batch in count_plan['batches']):,} bytes. "
        f"Exact plan-file SHA-256: `{count_plan_sha256}`.",
        f"- If later approved, source-row plan: {row_plan['n_lenses']} lens centers in "
        f"{row_plan['n_lens_row_sql_batches']} batches and {row_plan['n_random_centers']} controls in "
        f"{row_plan['n_random_row_sql_batches']} batches. Expected lens rows from the verified counts are "
        f"{row_plan['expected_lens_source_pair_rows_from_approved_counts']:,}; confirm per-center row counts "
        "when returned. Random rows await the count-only result.",
        f"- Random-pair planning proxies are {random_proxy['q16_proxy']:,}/"
        f"{random_proxy['q50_proxy']:,}/{random_proxy['q84_proxy']:,}; combined lens+random proxies are "
        f"{pair_proxy['q16_proxy']:,}/{pair_proxy['q50_proxy']:,}/{pair_proxy['q84_proxy']:,} "
        "(q16/median/q84). These use a different completed 620-lens sample and are not observed S23B counts.",
        f"- Estimated combined uncompressed CSV size: {fmt_gib(estimates['low'])} low, "
        f"{fmt_gib(estimates['point'])} point, {fmt_gib(estimates['high'])} high. "
        "The point estimate uses measured historical bytes per pair; low/high use 256/480 bytes per pair.",
        "",
        "Lens rows contain 22 columns: `lens_id`, lens RA/Dec/redshift, `object_id`, source RA/Dec, "
        "`e1`, `e2`, `sigma_e`, `rms_e`, `weight`, `m`, `c1`, `c2`, `mag_a10`, `resolution`, "
        "three PSF-shape moments, raw B-mode mask, and HSC Y3 tomographic bin. Control rows contain "
        "the corresponding 23 columns, adding `random_center_id` and `matched_lens_id` and substituting "
        "the random center coordinates/redshift. Full per-object photo-z PDFs and images are excluded.",
        "",
        "## Analysis gates",
        "",
        "The 400-center count-only query is a separate action and needs approval for the exact plan hash above. "
        "After its counts, revise the control balance and exact transfer estimate. Any shape-row query and "
        "transfer needs a second, separate approval tied to the final row-plan hash. Apply identical shape, "
        "source-redshift, B-mode, and radial cuts to lenses and controls; verify the raw B-mode polarity and "
        "supply calibrated source n(z) and lensing efficiency. Then evaluate footprint/depth/PSF support, "
        "control overlap, and spatial covariance before stacking. `science_measurement=false` remains in force.",
        "",
        f"Output directory: `{path.parent}`  ",
        f"Row-plan file SHA-256: `{row_plan_file_sha256}`  ",
        "Count plan and row plan are `prepared_not_submitted`; this planner used no network and returned no data.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_plan(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing plan directory: {output}")

    count_plan, count_summary, counts = load_verified_counts()
    plan_sha = count_plan["query_plan_sha256"]
    target_manifest_path = TARGET_PLAN / "target_manifest.json"
    if sha(target_manifest_path) != EXPECTED_TARGET_MANIFEST_SHA256:
        raise ValueError("S23B target manifest differs from its pinned checksum")
    target_doc = read_json(target_manifest_path)
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    centers = {str(row["lens_id"]): row for row in count_plan["centers"]}
    if set(counts) != set(centers) or not set(centers).issubset(targets):
        raise ValueError("approved lens IDs do not map one-to-one to target-manifest clusters")
    for lens_id, center in centers.items():
        target = targets[lens_id]
        for key1, key2 in (("ra_deg", "ra_deg"), ("dec_deg", "dec_deg"), ("z_lens", "z")):
            if not math.isclose(float(center[key1]), float(target[key2]), rel_tol=0, abs_tol=1e-9):
                raise ValueError(f"count center and target manifest disagree: {lens_id}/{key1}")

    control_summary_path = CONTROL_AUDIT / "summary.json"
    candidate_path = CONTROL_AUDIT / "candidate_audit.csv"
    control_summary = read_json(control_summary_path)
    if (control_summary.get("schema") != "astrafek-s23b-icl-random-control-candidate-audit-v1"
            or control_summary.get("status") != "complete_provisional_candidate_audit"
            or control_summary.get("plan_sha256") != EXPECTED_CONTROL_CANDIDATE_AUDIT_PLAN_SHA256
            or control_summary.get("science_measurement") is not False):
        raise ValueError("S23B random-control candidate audit is incomplete or changed")
    if sha(candidate_path) != control_summary["outputs_sha256"]["candidate_audit.csv"]:
        raise ValueError("S23B candidate table checksum mismatch")

    source_union_summary_path = SOURCE_UNION_AUDIT / "summary.json"
    source_union_table_path = SOURCE_UNION_AUDIT / "per_target_source_union.tsv"
    source_union_summary = read_json(source_union_summary_path)
    if (source_union_summary.get("status") != "complete_local_only"
            or sha(source_union_table_path) != source_union_summary.get("per_target_tsv_sha256")
            or source_union_summary.get("s23b_target_count") != 427
            or source_union_summary.get("science_measurement") is not False):
        raise ValueError("previous source-union audit does not validate")
    union_rows = {row["target_id"]: row for row in read_tsv(source_union_table_path)}
    if set(union_rows) != set(counts):
        raise ValueError("previous source union does not cover the exact 427-center set")

    random_support_path = OLD_RANDOM_SUPPORT / "lens_paired_support_audit_ratio_v25.json"
    random_support = read_json(random_support_path)
    ratios = random_support["lens_support_comparison"][
        "positive_lens_random_mean_to_lens_count_ratio_q16_q50_q84"
    ]
    if (random_support.get("status") != "complete"
            or random_support.get("plan_sha256") != EXPECTED_RANDOM_SUPPORT_PLAN_SHA256
            or random_support.get("science_measurement") is not False):
        raise ValueError("previous matched-random support audit is not valid planning evidence")

    cosmo = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Tcmb0=0)
    quota_by_separation: dict[float, dict[str, list[dict]]] = {
        value: {} for value in (6.0, 8.0, 10.0)
    }
    support_lenses = {lens_id for lens_id, value in counts.items() if value >= MIN_LENS_COUNT}
    with candidate_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            lens_id = str(row["target_cluster_id"]).replace(":", "_")
            if lens_id not in support_lenses:
                continue
            if row.get("center_flags_pass") != "True" or row.get("metadata_match") != "True":
                continue
            if int(row["cluster_overlap_2mpc"]) != 0:
                continue
            z = float(row["target_z"])
            if not math.isclose(z, float(centers[lens_id]["z_lens"]), rel_tol=0, abs_tol=1e-9):
                raise ValueError(f"control candidate redshift does not match its lens: {lens_id}")
            da = float(cosmo.angular_diameter_distance(z).value)
            offset_mpc = float(row["offset_from_target_arcmin"]) * da * math.pi / 10800.0
            candidate = {
                "random_object_id": str(row["random_object_id"]),
                "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
                "target_z": z, "offset_from_target_mpc": offset_mpc,
                "offset_from_target_arcmin": float(row["offset_from_target_arcmin"]),
                "match_score": float(row["match_score"]), "tract": str(row["tract"]),
                "patch": str(row["patch"]), "skymap_id": str(row["skymap_id"]),
                "delta_depth_mag": float(row["delta_depth_mag"]),
                "seeing_ratio": float(row["seeing_ratio"]),
                "delta_ebv_sfd": float(row["delta_ebv_sfd"]),
                "cluster_overlap_2mpc": int(row["cluster_overlap_2mpc"]),
            }
            for minimum in quota_by_separation:
                if offset_mpc > minimum:
                    group = quota_by_separation[minimum].setdefault(lens_id, [])
                    group.append(candidate)

    for groups in quota_by_separation.values():
        for lens_id, options in groups.items():
            unique = {}
            for option in options:
                position = (round(option["ra_deg"], 12), round(option["dec_deg"], 12))
                if position not in unique or (option["match_score"], option["random_object_id"]) < (
                        unique[position]["match_score"], unique[position]["random_object_id"]):
                    unique[position] = option
            groups[lens_id] = sorted(unique.values(),
                                     key=lambda item: (item["match_score"], item["random_object_id"]))

    groups6 = quota_by_separation[MIN_TARGET_CONTROL_SEPARATION_MPC]
    eligible_lenses = {lens_id for lens_id in support_lenses
                       if len(groups6.get(lens_id, [])) >= CONTROLS_PER_LENS}
    assignment = select_unique_controls(groups6, CONTROLS_PER_LENS)
    selected_lens_ids = set(eligible_lenses)
    if {slot[0] for slot in assignment} != selected_lens_ids:
        raise ValueError("global control assignment changed the common-support lens set")

    ordered_lenses = [str(row["lens_id"]) for row in count_plan["centers"]
                      if str(row["lens_id"]) in selected_lens_ids]
    selected_centers = [dict(centers[lens_id]) for lens_id in ordered_lenses]
    lens_rows_expected = sum(counts[lens_id] for lens_id in ordered_lenses)
    if len(selected_centers) != 80 or len(assignment) != 400:
        raise ValueError("reproducible sample size differs from the audited 80 x 5 design")

    random_centers = []
    used_random_ids = set()
    for lens_id in ordered_lenses:
        for slot_index in range(CONTROLS_PER_LENS):
            candidate = assignment[(lens_id, slot_index)]
            random_id = candidate["random_object_id"]
            if random_id in used_random_ids:
                raise ValueError("random table object ID was reused across control lenses")
            used_random_ids.add(random_id)
            center_id = "RC-" + hashlib.sha256(
                f"{lens_id}\0{random_id}\0{candidate['ra_deg']:.12f}\0{candidate['dec_deg']:.12f}".encode()
            ).hexdigest()[:20]
            lens = centers[lens_id]
            random = {
                "random_center_id": center_id,
                "matched_lens_id": lens_id,
                "ra_deg": candidate["ra_deg"], "dec_deg": candidate["dec_deg"],
                "z_lens": float(lens["z_lens"]),
                "radius_mpc_proper": 3.0,
                "radius_arcsec": float(lens["radius_arcsec"]),
                "random_object_id": random_id,
                "control_slot": slot_index + 1,
                **{key: candidate[key] for key in (
                    "offset_from_target_mpc", "offset_from_target_arcmin", "match_score",
                    "tract", "patch", "skymap_id", "delta_depth_mag", "seeing_ratio",
                    "delta_ebv_sfd", "cluster_overlap_2mpc",
                )},
            }
            random_centers.append(random)

    if len({(round(row["ra_deg"], 12), round(row["dec_deg"], 12))
            for row in random_centers}) != len(random_centers):
        raise ValueError("selected random positions are not globally unique")

    output.mkdir(parents=True, exist_ok=False)
    lens_center_path = output / "selected_lenses.tsv"
    lens_fields = ["lens_id", "catalog_id", "ra_deg", "dec_deg", "z_lens", "richness",
                   "n_raw_shape_sources", "eligible_control_candidates_6mpc",
                   "eligible_control_candidates_8mpc", "eligible_control_candidates_10mpc",
                   "old_union_sources", "old_union_area_fraction_4096"]
    lens_output = []
    for lens_id in ordered_lenses:
        target = targets[lens_id]
        union = union_rows[lens_id]
        lens_output.append({
            "lens_id": lens_id, "catalog_id": target["catalog_id"],
            "ra_deg": target["ra_deg"], "dec_deg": target["dec_deg"],
            "z_lens": target["z"], "richness": target["richness"],
            "n_raw_shape_sources": counts[lens_id],
            "eligible_control_candidates_6mpc": len(groups6.get(lens_id, [])),
            "eligible_control_candidates_8mpc": len(quota_by_separation[8.0].get(lens_id, [])),
            "eligible_control_candidates_10mpc": len(quota_by_separation[10.0].get(lens_id, [])),
            "old_union_sources": union["n_union_unique_sources"],
            "old_union_area_fraction_4096": union["all_query_union_area_fraction_4096"],
        })
    write_tsv(lens_center_path, lens_output, lens_fields)

    control_path = output / "selected_random_controls.tsv"
    control_fields = ["random_center_id", "matched_lens_id", "control_slot", "random_object_id",
                      "ra_deg", "dec_deg", "z_lens", "radius_arcsec", "offset_from_target_mpc",
                      "match_score", "tract", "patch", "skymap_id", "delta_depth_mag",
                      "seeing_ratio", "delta_ebv_sfd", "cluster_overlap_2mpc"]
    control_output = [{key: row[key] for key in control_fields} for row in random_centers]
    write_tsv(control_path, control_output, control_fields)

    random_manifest_path = output / "random_centers.json"
    random_manifest = {
        "schema": "astrafek-wl-random-centers-v1",
        "status": "selected_for_shape_count_preflight",
        "h0_km_s_mpc": H0, "omega_m": OMEGA_M,
        "science_measurement": False,
        "selection": {
            "method": "global one-to-one deterministic candidate matching; candidates ranked by metadata match_score and random table object ID",
        "seed": "S23B-WL-common-support-2026-11-10",
            "requested_per_lens": CONTROLS_PER_LENS,
            "minimum_candidates_per_lens": CONTROLS_PER_LENS,
            "reuse_scope": "global_sky_position",
            "lens_source_count_threshold": MIN_LENS_COUNT,
            "minimum_target_control_center_separation_mpc": MIN_TARGET_CONTROL_SEPARATION_MPC,
            "metadata_thresholds": control_summary["match_thresholds"],
            "known_cluster_screen": "candidate_audit.cluster_overlap_2mpc == 0; prior screen uses S20A+S23B CAMIRA centers and the declared 1-Mpc image-footprint cap",
            "n_matched_lenses": len(ordered_lenses),
            "n_random_centers": len(random_centers),
            "per_lens_counts": {lens_id: {"available": len(groups6[lens_id]), "selected": CONTROLS_PER_LENS}
                                for lens_id in ordered_lenses},
        },
        "full_aperture_footprint_verified": False,
        "depth_psf_sky_matched": False,
        "random_center_overlap_audited": True,
        "spatial_covariance_validated": False,
        "execution_requires_explicit_user_approval": True,
        "centers": [{key: row[key] for key in (
            "random_center_id", "matched_lens_id", "ra_deg", "dec_deg", "z_lens",
            "radius_mpc_proper", "radius_arcsec", "random_object_id", "offset_from_target_mpc",
            "match_score", "tract", "patch", "skymap_id", "delta_depth_mag", "seeing_ratio",
            "delta_ebv_sfd", "cluster_overlap_2mpc", "control_slot",
        )} for row in random_centers],
    }
    write_json(random_manifest_path, random_manifest)

    normalized_random_centers, _, _, n_matched_lenses = normalize_centers(
        random_manifest, radius_mpc=3.0, batch_size=BATCH_SIZE,
    )
    if len(normalized_random_centers) != len(random_centers) or n_matched_lenses != len(ordered_lenses):
        raise ValueError("random-center normalization changed the selected control design")
    count_batches = []
    for index, start in enumerate(range(0, len(normalized_random_centers), BATCH_SIZE)):
        chunk = normalized_random_centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(shape_count_query(row) for row in chunk) + "\n"
        name = f"random_shape_counts_{index:03d}.sql"
        path = output / name
        path.write_text(sql, encoding="utf-8")
        count_batches.append({
            "path": str(path), "n_random_centers": len(chunk),
            "first_random_center_id": chunk[0]["random_center_id"],
            "last_random_center_id": chunk[-1]["random_center_id"],
            "sha256": sha(path), "sql_bytes": path.stat().st_size,
        })
    count_plan = {
        "schema": "hsc-s19a-random-shape-count-plan-v1",
        "status": "prepared_not_submitted",
        "source_centers": str(random_manifest_path),
        "source_centers_sha256": sha(random_manifest_path),
        "release": "pdr3", "query_tables": [
            "s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas",
        ],
        "n_random_centers": len(random_centers), "n_matched_lenses": len(ordered_lenses),
        "n_query_batches": len(count_batches), "batch_size": BATCH_SIZE,
        "radius_mpc_proper": 3.0, "cosmology": {"h0_km_s_mpc": H0, "omega_m": OMEGA_M},
        "output_columns": COUNT_COLUMNS,
        "query_semantics": "raw distinct S19A shape-object count inside each random-center proper-radius cone; before shape, redshift, B-mode, and radial cuts",
        "source_shape_rows_returned": False, "photoz_rows_returned": False,
        "images_returned": False, "network_used": False, "query_submitted": False,
        "credentials_saved": False, "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "batches": count_batches,
        "later_gates": [
            "This is a new control-center count-only action, separate from the approved 427 lens-center count plan.",
            "The counts are required before an exact random source-row total can be stated.",
            "Passing a raw count screen does not establish matched source density after redshift, shape-quality, B-mode, or radial selections.",
        ],
    }
    count_plan_path = output / "random_shape_count_plan.json"
    write_json(count_plan_path, count_plan)
    count_plan_sha256 = sha(count_plan_path)
    validated_manifest, validated_centers, validated_batches, validated_sha = validate_random_count_plan(output)
    if (validated_sha != count_plan_sha256 or len(validated_centers) != 400
            or len(validated_batches) != 16
            or validated_manifest.get("query_submitted") is not False):
        raise ValueError("existing random count-plan validator rejected the offline plan")
    (output / "RANDOM_COUNT_PLAN_SHA256").write_text(count_plan_sha256 + "\n", encoding="ascii")

    lens_batches = []
    for index, start in enumerate(range(0, len(selected_centers), BATCH_SIZE)):
        chunk = selected_centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(lens_query(row) for row in chunk) + "\n"
        name = f"lens_shape_rows_{index:03d}.sql"
        path = output / name
        path.write_text(sql, encoding="utf-8")
        lens_batches.append({"name": name, "n_lenses": len(chunk),
                             "expected_rows_from_approved_counts": sum(counts[row["lens_id"]] for row in chunk),
                             "sha256": sha(path), "sql_bytes": path.stat().st_size})

    random_row_batches = []
    for index, start in enumerate(range(0, len(random_centers), BATCH_SIZE)):
        chunk = random_centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(shape_rows_query(row) for row in chunk) + "\n"
        name = f"random_shape_rows_{index:03d}.sql"
        path = output / name
        path.write_text(sql, encoding="utf-8")
        random_row_batches.append({"name": name, "n_random_centers": len(chunk),
                                  "expected_rows_from_count_preflight": None,
                                  "sha256": sha(path), "sql_bytes": path.stat().st_size})

    target_pairs = lens_rows_expected
    q16, q50, q84 = (float(ratios[key]) for key in ("q16", "q50", "q84"))
    control_pair_estimates = {
        key: int(round(target_pairs * CONTROLS_PER_LENS * ratio))
        for key, ratio in (("q16_proxy", q16), ("q50_proxy", q50), ("q84_proxy", q84))
    }
    total_pair_estimates = {key: target_pairs + value
                            for key, value in control_pair_estimates.items()}
    old_lens_summary = read_json(OLD_LENS_EXPORT / "row_summary.json")
    old_random_summary_path = OLD_RANDOM_EXPORT / "row_summary.json"
    old_random_summary = read_json(old_random_summary_path)
    old_lens_n_rows = int(old_lens_summary["n_pair_rows"])
    old_lens_bytes = sum(
        (OLD_LENS_EXPORT / f"shape_rows_{index:03d}.csv").stat().st_size
        for index in range(32)
    )
    if (old_lens_summary.get("status") != "complete"
            or old_lens_summary.get("source_rows_downloaded") is not True
            or old_lens_summary.get("science_measurement") is not False
            or old_lens_n_rows != 3_069_544 or old_lens_bytes != 909_833_598
            or int(old_lens_summary["count_summary_sum"]) != old_lens_n_rows):
        raise ValueError("completed lens export row or byte total differs from its audited values")
    old_lens_bytes_per_row = old_lens_bytes / old_lens_n_rows
    random_n_rows = int(old_random_summary["n_source_pairs"])
    random_n_bytes = int(old_random_summary["n_csv_bytes"])
    if (random_n_rows != 9_010_324 or random_n_bytes != 2_962_522_918
            or old_random_summary.get("status") != "complete"
            or old_random_summary.get("source_shape_rows_downloaded") is not True
            or old_random_summary.get("science_measurement") is not False):
        raise ValueError("completed random export row or byte total differs from its audited values")
    old_random_bytes_per_row = random_n_bytes / random_n_rows
    estimated_bytes = {
        "low": int(round(total_pair_estimates["q16_proxy"] * ROW_BYTES_LOW)),
        "point": int(round(target_pairs * old_lens_bytes_per_row
                            + control_pair_estimates["q50_proxy"] * old_random_bytes_per_row)),
        "high": int(round(total_pair_estimates["q84_proxy"] * ROW_BYTES_HIGH)),
        "planning_bytes_per_row_range": [ROW_BYTES_LOW, ROW_BYTES_HIGH],
        "point_basis_bytes_per_row": {
            "lens": old_lens_bytes_per_row, "random_control": old_random_bytes_per_row,
        },
        "unit": "uncompressed CSV bytes",
        "basis": "q16/q50/q84 control-to-lens mean-count ratios from a different, completed 620-lens/3,100-random-center sample; planning proxy only, not S23B control counts",
    }
    row_plan = {
        "schema": "astrafek-s23b-wl-common-support-shape-row-plan-v1",
        "status": "prepared_not_submitted",
        "count_only_preflight_plan_sha256": count_plan_sha256,
        "lens_count_plan_sha256": plan_sha,
        "n_lenses": len(selected_centers), "n_random_centers": len(random_centers),
        "random_centers_per_lens": CONTROLS_PER_LENS,
        "n_lens_row_sql_batches": len(lens_batches),
        "n_random_row_sql_batches": len(random_row_batches),
        "expected_lens_source_pair_rows_from_approved_counts": target_pairs,
        "expected_random_source_pair_rows": None,
        "random_source_pair_rows_planning_proxies": control_pair_estimates,
        "combined_source_pair_rows_planning_proxies": total_pair_estimates,
        "estimated_uncompressed_csv_bytes": estimated_bytes,
        "lens_output_columns": LENS_COLUMNS,
        "random_output_columns": RANDOM_COLUMNS,
        "source_release": "HSC PDR3 S19A Wide",
        "source_tables": ["s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas"],
        "row_semantics": {
            "lens": "one row per (lens_id, object_id); shared source IDs across different lenses are retained",
            "random": "one row per (random_center_id, object_id); shared source IDs across random centers are retained",
        },
        "selection_rules": [
            "S23B CAMIRA target lens has at least 1,000 released S19A shape-catalog objects in the pre-project 3-proper-Mpc count; this is a provisional availability screen, not a post-cut source count.",
            "Retain only targets with five globally unique matched HSC random-table control candidates from the saved candidate audit.",
            "Each control candidate passes the audited center flags and metadata match: |depth difference| <= 0.30 mag, seeing ratio 1/1.20 to 1.20, and |SFD E(B-V) difference| <= 0.02.",
            "Each random-control center is more than 6 proper Mpc from its matched target center, so the two 3-proper-Mpc query apertures do not geometrically overlap.",
            "Each selected control has candidate_audit.cluster_overlap_2mpc == 0; this previous S20A+S23B CAMIRA screen includes a declared image-footprint margin, but it is an incomplete cluster-contamination tracer.",
            "After an approved row export, apply exactly the same shape-quality, source-redshift-bin, B-mode, and radial cuts to lenses and controls; verify B-mode raw-value polarity and use calibrated HSC Y3 n(z) before any lensing measurement.",
        ],
        "row_selection_fields_preserved": ["sigma_e", "rms_e", "weight", "m", "c1", "c2",
                                           "mag_a10", "resolution", "psf_shape11", "psf_shape22",
                                           "psf_shape12", "b_mode_mask_raw", "source_zbin"],
        "photoz_pdf_rows_included": False, "images_included": False,
        "count_preflight_required_before_exact_control_row_total": True,
        "network_used": False, "query_submitted": False,
        "source_rows_downloaded": False, "credentials_saved": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "lens_batches": lens_batches,
        "random_batches": random_row_batches,
        "notes": [
            "No source-row or count-only query is submitted by this planner.",
            "The 400-control count-only plan is a separate data request and needs plan-hash-specific approval.",
            "A passing raw count preflight still will not establish matched source support after project source cuts.",
        ],
    }
    row_plan["query_plan_sha256"] = canonical_sha(row_plan)
    row_plan_path = output / "row_export_plan.json"
    write_json(row_plan_path, row_plan)
    row_plan_file_sha256 = sha(row_plan_path)
    (output / "ROW_EXPORT_PLAN_SHA256").write_text(row_plan_file_sha256 + "\n",
                                                   encoding="ascii")

    camira_catalog = load_camira_catalog(CAMIRA_CATALOG)
    if sha(CAMIRA_CATALOG) != EXPECTED_CAMIRA_CATALOG_SHA256:
        raise ValueError("pinned CAMIRA catalog changed before contamination diagnostics")
    control_manifest = read_json(DATA / "s23b_icl_control_image_plan_2026-11-07_01" /
                                 "control_manifest.json")
    control_manifest_path = DATA / "s23b_icl_control_image_plan_2026-11-07_01" / "control_manifest.json"
    old_control = count_existing_control_contamination(
        control_manifest, targets, set(counts), camira_catalog,
    )
    selected_offsets = [row["offset_from_target_mpc"] for row in random_centers]
    selected_ra_bins = {}
    for lower in range(0, 360, 30):
        ids = [lens_id for lens_id in ordered_lenses
               if lower <= float(targets[lens_id]["ra_deg"]) < lower + 30]
        if ids:
            selected_ra_bins[f"{lower:03d}-{lower + 30:03d}"] = {
                "n_lenses": len(ids), "lens_source_pairs": sum(counts[lens_id] for lens_id in ids),
            }
    parent_stats = {
        "n_lenses": len(counts),
        "zero_raw_counts": sum(value == 0 for value in counts.values()),
        "positive_raw_counts": sum(value > 0 for value in counts.values()),
        "at_least_1000": sum(value >= MIN_LENS_COUNT for value in counts.values()),
        "pair_rows": sum(counts.values()),
    }
    selected_props = [targets[lens_id] for lens_id in ordered_lenses]
    sample_stats = {
        "n_lenses": len(selected_props), "n_controls": len(random_centers),
        "exact_lens_pair_rows": target_pairs,
        "median_z": statistics.median(float(row["z"]) for row in selected_props),
        "median_richness": statistics.median(float(row["richness"]) for row in selected_props),
        "lens_count_min": min(counts[lens_id] for lens_id in ordered_lenses),
        "lens_count_median": statistics.median(counts[lens_id] for lens_id in ordered_lenses),
        "lens_count_max": max(counts[lens_id] for lens_id in ordered_lenses),
        "median_control_match_score": statistics.median(row["match_score"] for row in random_centers),
        "median_control_offset_mpc": statistics.median(selected_offsets),
        "p90_control_offset_mpc": sorted(selected_offsets)[int(0.9 * (len(selected_offsets) - 1))],
        "lens_counts_by_ra_30deg": selected_ra_bins,
    }

    assessment = {
        "schema": "astrafek-s23b-wl-common-support-assessment-v1",
        "status": "complete_local_candidate_support_assessment",
        "science_measurement": False,
        "network_used": False,
        "count_only_query_submitted": False,
        "source_row_query_submitted": False,
        "source_rows_downloaded": False,
        "credentials_saved": False,
        "approved_lens_count_plan_sha256": plan_sha,
        "count_audit_summary_sha256": sha(COUNT_AUDIT / "summary.json"),
        "target_manifest_sha256": sha(target_manifest_path),
        "control_candidate_plan_sha256": EXPECTED_CONTROL_CANDIDATE_AUDIT_PLAN_SHA256,
        "control_candidate_table_sha256": sha(candidate_path),
        "existing_source_union_summary_sha256": sha(source_union_summary_path),
        "existing_source_union_selected_lens_rows_with_any_sources": sum(
            int(union_rows[lens_id]["n_union_unique_sources"]) > 0 for lens_id in ordered_lenses
        ),
        "existing_source_union_selected_lens_rows_geometric_coverage_95pct": sum(
            float(union_rows[lens_id]["all_query_union_area_fraction_4096"]) >= 0.95
            for lens_id in ordered_lenses
        ),
        "parent_lens_support": parent_stats,
        "candidate_lens_support_threshold": MIN_LENS_COUNT,
        "candidate_control_quota": CONTROLS_PER_LENS,
        "candidate_control_separation_screen_mpc": {
            "6": sum(len(quota_by_separation[6.0].get(lens_id, [])) >= 5 for lens_id in support_lenses),
            "8": sum(len(quota_by_separation[8.0].get(lens_id, [])) >= 5 for lens_id in support_lenses),
            "10": sum(len(quota_by_separation[10.0].get(lens_id, [])) >= 5 for lens_id in support_lenses),
        },
        "common_support_candidate": sample_stats,
        "common_support_exclusions": {
            "n_lenses_with_raw_count_at_least_1000": len(support_lenses),
            "n_lenses_with_five_valid_control_candidates": len(selected_lens_ids),
            "excluded_for_fewer_than_five_controls": sorted(support_lenses - selected_lens_ids),
        },
        "existing_selected_icl_control_diagnostic": old_control,
        "previous_random_control_count_ratio_proxy": {
            "prior_plan_sha256": random_support["plan_sha256"],
            "prior_matched_lenses": random_support["n_complete_lenses"],
            "prior_random_centers": random_support["n_verified_centers"],
            "mean_random_count_to_lens_count_ratio_q16_q50_q84": ratios,
            "interpretation": "planning proxy only; this is a different 620-lens sample and does not replace S23B control counts",
        },
        "row_export_plan_sha256": row_plan_file_sha256,
        "row_export_plan_canonical_sha256": row_plan["query_plan_sha256"],
        "random_control_count_plan_sha256": count_plan_sha256,
        "limitations": [
            "The 1,000-object target threshold is a provisional raw availability screen before project source cuts; it does not establish final signal-to-noise.",
            "The selected control centers have not been queried in the S19A shape table. They are geometrically and metadata matched candidates, not confirmed common source support.",
            "The 2-Mpc CAMIRA candidate screen is based on incomplete cluster tracers; environment and large-scale lensing can persist beyond the geometric aperture.",
            "Control-count expectations use a different completed sample only as a planning proxy; exact row and volume totals require the separately approved 400-center count-only query.",
            "B-mode value polarity, project source cuts, calibrated n(z), lensing efficiency, and spatial covariance are not validated here.",
        ],
        "inputs_sha256": {
            str(COUNT_PLAN / "shape_count_plan.json"): sha(COUNT_PLAN / "shape_count_plan.json"),
            str(COUNT_PLAN / "results" / "count_status.json"): sha(COUNT_PLAN / "results" / "count_status.json"),
            str(COUNT_PLAN / "results" / "count_summary.json"): sha(COUNT_PLAN / "results" / "count_summary.json"),
            str(target_manifest_path): sha(target_manifest_path),
            str(control_summary_path): sha(control_summary_path),
            str(candidate_path): sha(candidate_path),
            str(source_union_summary_path): sha(source_union_summary_path),
            str(source_union_table_path): sha(source_union_table_path),
            str(CAMIRA_CATALOG): sha(CAMIRA_CATALOG),
            str(control_manifest_path): sha(control_manifest_path),
            str(random_support_path): sha(random_support_path),
        },
    }
    write_json(output / "support_assessment.json", assessment)
    write_json(output / "random_count_plan_summary.json", {
        "status": count_plan["status"],
        "query_plan_sha256": count_plan_sha256,
        "n_random_centers": count_plan["n_random_centers"],
        "n_matched_lenses": count_plan["n_matched_lenses"],
        "n_query_batches": count_plan["n_query_batches"],
        "output_columns": count_plan["output_columns"],
        "query_submitted": False, "science_measurement": False,
    })
    write_json(output / "approval_summary.json", {
        "status": "prepared_not_submitted",
        "random_control_count_plan_sha256": count_plan_sha256,
        "row_export_plan_sha256": row_plan_file_sha256,
        "row_export_plan_canonical_sha256": row_plan["query_plan_sha256"],
        "n_count_only_centers": count_plan["n_random_centers"],
        "n_count_only_batches": count_plan["n_query_batches"],
        "n_provisional_lenses": len(selected_centers),
        "n_provisional_random_centers": len(random_centers),
        "expected_lens_pair_rows": target_pairs,
        "expected_control_pair_rows": None,
        "control_row_estimates": control_pair_estimates,
        "network_used": False, "query_submitted": False,
        "source_rows_downloaded": False, "science_measurement": False,
        "requires_explicit_user_approval": True,
    })
    write_assessment_markdown(
        output / "assessment.md", assessment, row_plan, count_plan,
        count_plan_sha256, row_plan_file_sha256,
    )
    return assessment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_plan(args.output)
    print(json.dumps({
        "status": result["status"],
        "lenses": result["common_support_candidate"]["n_lenses"],
        "random_centers": result["common_support_candidate"]["n_controls"],
        "lens_pair_rows_exact": result["common_support_candidate"]["exact_lens_pair_rows"],
        "count_plan_sha256": result["random_control_count_plan_sha256"],
        "row_plan_sha256": result["row_export_plan_sha256"],
        "network_used": result["network_used"],
        "query_submitted": result["source_row_query_submitted"],
        "output": str(Path(args.output).resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
