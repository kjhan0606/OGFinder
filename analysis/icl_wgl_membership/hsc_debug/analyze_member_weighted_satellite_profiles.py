#!/usr/bin/env python3
"""Attach CAMIRA member stellar-mass/weight fields to exploratory HSC profiles.

This is a descriptive sensitivity analysis. The headerless CAMIRA member file's
first six columns are cross-identified against the cluster manifest and saved
member positions. The acquisition parser and its hash-matched summary identify
column 8 as ``weight``; column 7 remains undocumented, though its numeric range
is consistent with log10(M*/Msun). Weight results use the identified CAMIRA
weight, while mass-bin results remain provisional. It does not download data
or make a science measurement.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

EXPECTED_MEMBER_SHA = "bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e"
MASS_BINS = ((8.0, 10.5, "<10.5"), (10.5, 11.0, "10.5-11.0"),
             (11.0, 11.5, "11.0-11.5"), (11.5, 13.5, ">=11.5"))
MEMBER_COLUMNS = ["ra_cl_deg", "dec_cl_deg", "N_mem", "z_cl",
                  "ra_deg", "dec_deg", "field7", "weight"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def validate_member_metadata(summary_path, member_path):
    """Require the acquisition summary that documents this exact member file."""
    summary_path = Path(summary_path)
    if not summary_path.is_file():
        raise ValueError(f"Missing CAMIRA member acquisition summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("source_file") != Path(member_path).name:
        raise ValueError("CAMIRA acquisition summary names a different source file")
    if summary.get("sha256") != sha256(member_path):
        raise ValueError("CAMIRA acquisition summary does not match the member-file SHA-256")
    if summary.get("columns") != MEMBER_COLUMNS:
        raise ValueError("CAMIRA acquisition summary has an unexpected column mapping")
    return summary


def rounded_coord(ra, dec):
    return round(float(ra), 6), round(float(dec), 6)


def read_catalog(member_path, clusters):
    """Read only the 800 selected cluster groups and verify manifest identity."""
    by_cluster = {
        (round(float(c["ra_deg"]), 5), round(float(c["dec_deg"]), 5)): c
        for c in clusters
    }
    if len(by_cluster) != len(clusters):
        raise ValueError("Manifest has duplicate rounded cluster coordinates")
    grouped = {key: [] for key in by_cluster}
    with gzip.open(member_path, "rt", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            values = [float(value) for value in line.split()]
            if len(values) != 8 or not all(math.isfinite(v) for v in values):
                raise ValueError(f"Malformed/non-finite CAMIRA row {line_number}")
            key = (round(values[0], 5), round(values[1], 5))
            if key in grouped:
                grouped[key].append(values)

    indexed = {}
    total = 0
    for key, cluster in by_cluster.items():
        rows = grouped[key]
        if not rows:
            raise ValueError(f"No CAMIRA member rows for cluster {cluster['id']}")
        if any(abs(r[2] - float(cluster["richness"])) > 0.001 for r in rows):
            raise ValueError(f"Richness mismatch for cluster {cluster['id']}")
        # The member file stores CAMIRA's z_cl (photometric), while manifest
        # ``z`` is the spectroscopic redshift used for physical radii.
        cluster_z = float(cluster.get("z_cl", cluster["z"]))
        if any(abs(r[3] - cluster_z) > 0.001 for r in rows):
            raise ValueError(f"Redshift mismatch for cluster {cluster['id']}")
        catalog_coords = [rounded_coord(r[4], r[5]) for r in rows]
        manifest_coords = [rounded_coord(*xy) for xy in cluster["members_radec_deg"]]
        if len(set(catalog_coords)) != len(catalog_coords):
            raise ValueError(f"Duplicate member coordinates for cluster {cluster['id']}")
        if set(catalog_coords) != set(manifest_coords):
            raise ValueError(f"Member-coordinate set mismatch for cluster {cluster['id']}")
        for row in rows:
            mass, weight = row[6], row[7]
            if not 8.0 < mass < 13.5:
                raise ValueError(f"Field 7 outside plausible log-mass domain: {mass}")
            if not 0.0 <= weight <= 1.0:
                raise ValueError(f"Field 8 outside CAMIRA weight domain: {weight}")
            coord = rounded_coord(row[4], row[5])
            indexed[(str(cluster["id"]), coord)] = {
                "candidate_log10_mstar": mass,
                "camira_weight": weight,
            }
        total += len(rows)
    return indexed, total


def read_preparation_crosswalk(prep_dir, cluster_ids, catalog_index):
    """Map persisted member indices to source-catalog rows by coordinates."""
    crosswalk = {}
    for ident in sorted(cluster_ids, key=str):
        member_path = Path(prep_dir) / f"cluster_{ident}" / "members.tsv"
        rows = read_tsv(member_path)
        local_indices = set()
        local_coords = set()
        for row in rows:
            index = int(row["INDEX"])
            coord = rounded_coord(row["RA"], row["DEC"])
            if index in local_indices or coord in local_coords:
                raise ValueError(f"Duplicate prepared member index/coordinate for {ident}")
            local_indices.add(index)
            local_coords.add(coord)
            source = catalog_index.get((str(ident), coord))
            if source is None:
                raise ValueError(f"Prepared member has no catalog match: {ident}/{index}")
            crosswalk[(str(ident), index)] = source
        if len(rows) != len(catalog_index_for_cluster(catalog_index, ident)):
            raise ValueError(f"Prepared member row count differs for cluster {ident}")
    return crosswalk


def catalog_index_for_cluster(index, ident):
    return [key for key in index if key[0] == str(ident)]


def aggregate(profile_rows, crosswalk):
    unique_keys = set()
    records = []
    for row in profile_rows:
        key = (str(row["cluster_id"]), row["filter"], int(row["member_index"]),
               row["radial_bin_kpc"])
        if key in unique_keys:
            raise ValueError(f"Duplicate profile row {key}")
        unique_keys.add(key)
        if float(row["member_detection_separation_arcsec"]) > 1.0:
            raise ValueError(f"Member/detection match exceeds 1 arcsec: {key}")
        props = crosswalk.get((key[0], key[2]))
        if props is None:
            raise ValueError(f"No member properties for profile row {key}")
        records.append({
            "cluster_id": key[0], "filter": key[1], "radial_bin_kpc": key[3],
            "member_index": key[2],
            "delta": float(row["delta_sb_njy_arcsec2"]),
            "candidate_log10_mstar": float(props["candidate_log10_mstar"]),
            "camira_weight": float(props["camira_weight"]),
        })

    out = []
    group_keys = sorted({(r["filter"], r["radial_bin_kpc"]) for r in records})
    selections = [("all", None, None)] + [
        (label, lo, hi) for lo, hi, label in MASS_BINS
    ]
    for band, radial_bin in group_keys:
        for selection, lo, hi in selections:
            selected = [r for r in records if r["filter"] == band
                        and r["radial_bin_kpc"] == radial_bin
                        and (lo is None or lo <= r["candidate_log10_mstar"] < hi)]
            by_cluster = defaultdict(list)
            for row in selected:
                by_cluster[row["cluster_id"]].append(row)
            if not selected:
                continue
            equal_cluster = []
            weighted_cluster = []
            for rows in by_cluster.values():
                values = np.asarray([r["delta"] for r in rows], float)
                weights = np.asarray([r["camira_weight"] for r in rows], float)
                equal_cluster.append(float(values.mean()))
                weighted_cluster.append(float(np.average(values, weights=weights))
                                          if weights.sum() > 0 else math.nan)
            valid_weighted = np.asarray(weighted_cluster, float)
            valid_weighted = valid_weighted[np.isfinite(valid_weighted)]
            all_weights = np.asarray([r["camira_weight"] for r in selected], float)
            neff = (float(all_weights.sum() ** 2 / np.sum(all_weights ** 2))
                    if np.sum(all_weights ** 2) > 0 else 0.0)
            out.append({
                "filter": band, "radial_bin_kpc": radial_bin,
                "mass_selection": selection,
                "n_clusters": len(by_cluster), "n_satellites": len(selected),
                "n_eff_weight": f"{neff:.6g}",
                "mean_log10_mstar": f"{np.mean([r['candidate_log10_mstar'] for r in selected]):.8g}",
                "mean_camira_weight": f"{all_weights.mean():.8g}",
                "equal_cluster_mean_delta_njy_arcsec2": f"{np.mean(equal_cluster):.10g}",
                "camira_weighted_equal_cluster_mean_delta_njy_arcsec2": (
                    f"{np.mean(valid_weighted):.10g}" if len(valid_weighted) else ""),
            })
    return records, out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--satellite", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
                                     "satellite_outer_profile_perbin_full_2026-10-18_01"))
    parser.add_argument("--catalog", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/catalog/"
                                     "camira_s20a_wide_sm_member.dat.gz"))
    parser.add_argument("--member-summary", type=Path,
                        help="SHA-matched member_summary.json produced by 06_fetch_members.py")
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    data, satellite, catalog, work = (p.resolve() for p in
                                       (args.data, args.satellite, args.catalog, args.work))
    member_summary = (args.member_summary.resolve() if args.member_summary else
                      catalog.parent.parent / "logs" / "member_summary.json")
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    summary_path = satellite / "summary.json"
    profile_path = satellite / "satellite_profiles.tsv"
    manifest_path = data / "sample_manifest.json"
    for path in (summary_path, profile_path, manifest_path, catalog):
        if not path.is_file():
            parser.error(f"Missing input: {path}")
    if sha256(catalog) != EXPECTED_MEMBER_SHA:
        parser.error("CAMIRA member catalog hash differs from the reviewed source")
    try:
        acquisition_metadata = validate_member_metadata(member_summary, catalog)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    satellite_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("n_frames_completed") != 800):
        parser.error("Expected the complete 800-frame exploratory profile input")
    if satellite_summary.get("output_sha256", {}).get(str(profile_path)) != sha256(profile_path):
        parser.error("Satellite profile table hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    clusters = manifest.get("clusters", [])
    if len(clusters) != 800 or len({str(c["id"]) for c in clusters}) != 800:
        parser.error("Expected 800 unique manifest clusters")

    catalog_index, n_catalog_rows = read_catalog(catalog, clusters)
    profiles = read_tsv(profile_path)
    cluster_ids = {str(row["cluster_id"]) for row in profiles}
    crosswalk = read_preparation_crosswalk(
        satellite_summary["preparation_directory"], cluster_ids, catalog_index)
    records, summaries = aggregate(profiles, crosswalk)

    work.mkdir(parents=True, exist_ok=False)
    detail_fields = ["cluster_id", "filter", "member_index", "radial_bin_kpc", "delta",
                     "candidate_log10_mstar", "camira_weight"]
    detail_rows = [{"cluster_id": r["cluster_id"], "filter": r["filter"],
                    "member_index": r["member_index"],
                    "radial_bin_kpc": r["radial_bin_kpc"], "delta": r["delta"],
                    "candidate_log10_mstar": r["candidate_log10_mstar"],
                    "camira_weight": r["camira_weight"]} for r in records]
    summary_fields = list(summaries[0]) if summaries else []
    write_tsv(work / "member_weighted_profiles.tsv", summaries, summary_fields)
    write_tsv(work / "profile_member_properties.tsv", detail_rows, detail_fields)
    result = {
        "schema": "hsc-camira-weighted-satellite-profile-audit-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "n_manifest_clusters": len(clusters), "n_catalog_clusters_matched": len(clusters),
        "n_catalog_member_rows_selected": n_catalog_rows,
        "n_member_profile_rows": len(records),
        "n_profile_clusters_crosswalked": len(cluster_ids),
        "acquisition_schema": {
            "source_summary": str(member_summary),
            "source_summary_sha256": sha256(member_summary),
            "columns": acquisition_metadata["columns"],
            "weight_field": "column 8, explicitly labeled weight by the original fetch/parser",
            "field7_note": acquisition_metadata.get("field7_note"),
        },
        "column_mapping": {
            "columns_1_to_6": "verified against manifest cluster RA/Dec/richness/z and member coordinate sets",
            "column_7": "undocumented field; analyzed only as a candidate log10(M_star/M_sun) because of its value domain and CAMIRA documentation",
            "column_8": "weight, identified by the original acquisition parser and its SHA-matched member_summary.json; interpreted as CAMIRA member weight w",
            "certainty": "column 8 is identified as weight; column 7 mass interpretation remains provisional",
        },
        "method": {
            "mass_bins_log10_msun": [list(x) for x in MASS_BINS],
            "weighting": "within each cluster and mass/radius/filter cell, weight satellite profiles by CAMIRA w; then equal weight clusters",
            "uncertainty": "none reported; these are descriptive sensitivity summaries, not inferential estimates",
            "limitations": [
                "CAMIRA w includes the algorithm's stellar-mass and spatial filters, so it is not an independent, unbiased satellite membership probability for this science question.",
                "The catalog's M_star derives from cmodel flux and an SPS mass-to-light calibration; it is not independent HSC satellite photometry.",
                "The input satellite envelopes are selected by image detection, masks, local controls, and annulus support; completeness and profile covariance remain uncalibrated.",
                "No extinction, passband, PSF-wing, sky-transfer, or coadd-correlated-noise correction is applied.",
            ],
        },
        "input_sha256": {str(path): sha256(path) for path in
                          (manifest_path, catalog, member_summary, summary_path, profile_path)},
        "script_sha256": sha256(Path(__file__)),
        "outputs_sha256": {},
    }
    result_path = work / "summary.json"
    for path in sorted(work.iterdir()):
        result["outputs_sha256"][path.name] = sha256(path)
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("status", "science_measurement", "n_catalog_member_rows_selected",
                       "n_member_profile_rows", "n_profile_clusters_crosswalked")},
                     ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
