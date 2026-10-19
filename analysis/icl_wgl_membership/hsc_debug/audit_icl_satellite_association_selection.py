#!/usr/bin/env python3
"""Describe how the 94-row ICL--satellite sample differs from ICL support rows."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_region_fixed_effects import (
    DEFAULT_ASSOC, DEFAULT_EXTRACTION, DEFAULT_ICL,
)
from audit_satellite_icl_spatial_blocks import load_verified_rows, read_tsv
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels, tract_centers_from_patch_qa, tract_region_labels,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"
DEFAULT_SATELLITE = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
NUMERIC_FIELDS = (
    "z", "richness", "delta_sb_njy_arcsec2", "target_coverage",
    "control_coverage", "common_coverage", "n_member_positions", "member_density_mpc2",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_value(value) -> bool:
    text = str(value).strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    raise ValueError(f"invalid saved boolean: {value!r}")


def finite_values(rows, field):
    values = []
    for row in rows:
        value = row.get(field)
        if value in (None, ""):
            continue
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"non-finite {field}")
        values.append(number)
    return np.asarray(values, dtype=float)


def describe(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("summary values must be finite and one-dimensional")
    if len(values) == 0:
        return {"n": 0, "mean": None, "sd": None, "median": None, "q25": None, "q75": None}
    return {
        "n": int(len(values)),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)) if len(values) > 1 else None,
        "median": float(np.median(values)),
        "q25": float(np.percentile(values, 25)),
        "q75": float(np.percentile(values, 75)),
    }


def standardized_mean_difference(selected, unselected):
    if len(selected) < 2 or len(unselected) < 2:
        return None
    variance = ((len(selected) - 1) * np.var(selected, ddof=1)
                + (len(unselected) - 1) * np.var(unselected, ddof=1))
    pooled_sd = math.sqrt(float(variance / (len(selected) + len(unselected) - 2)))
    return float((np.mean(selected) - np.mean(unselected)) / pooled_sd) if pooled_sd else None


def summarize_groups(parent, selected):
    selected_keys = {(str(row["id"]), str(row["filter"])) for row in selected}
    parent_keys = {(str(row["id"]), str(row["filter"])) for row in parent}
    if not selected_keys or len(selected_keys) != len(selected):
        raise ValueError("selected association keys are empty or duplicated")
    if not selected_keys.issubset(parent_keys) or len(parent_keys) != len(parent):
        raise ValueError("association rows do not uniquely match the ICL support parent")
    unselected = [row for row in parent
                  if (str(row["id"]), str(row["filter"])) not in selected_keys]
    output = {
        "n_icl_support_rows": len(parent),
        "n_association_rows": len(selected),
        "n_not_in_association": len(unselected),
        "association_fraction": len(selected) / len(parent),
        "numeric": {},
        "center_eligible_fraction": {},
    }
    for field in NUMERIC_FIELDS:
        all_selected = finite_values(selected, field)
        all_unselected = finite_values(unselected, field)
        output["numeric"][field] = {
            "association": describe(all_selected),
            "not_in_association": describe(all_unselected),
            "standardized_mean_difference": standardized_mean_difference(
                all_selected, all_unselected),
        }
    for name, group in (("association", selected), ("not_in_association", unselected)):
        output["center_eligible_fraction"][name] = (
            sum(bool_value(row["auto_unique_center_candidate"])
                and not bool_value(row["center_masked_k2p5"]) for row in group) / len(group)
            if group else None)
    return output


def summarize_support_overlap(icl_keys, satellite_keys, association_keys):
    icl_keys, satellite_keys, association_keys = map(set,
        (icl_keys, satellite_keys, association_keys))
    overlap = icl_keys & satellite_keys
    if association_keys != overlap:
        raise ValueError("association rows are not the exact ICL/satellite support intersection")
    result = {}
    for band in ("i", "i2"):
        i_keys = {key for key in icl_keys if key[1] == band}
        s_keys = {key for key in satellite_keys if key[1] == band}
        common = i_keys & s_keys
        result[band] = {
            "icl_supported": len(i_keys),
            "satellite_profile_supported": len(s_keys),
            "both_supported": len(common),
            "icl_supported_without_satellite_profile": len(i_keys - s_keys),
            "satellite_profile_without_icl_support": len(s_keys - i_keys),
            "association_rows": sum(key[1] == band for key in association_keys),
            "association_fraction_of_icl_support": len(common) / len(i_keys) if i_keys else None,
            "association_fraction_of_satellite_support": len(common) / len(s_keys) if s_keys else None,
        }
    return result


def analyze(association_dir: Path, icl_dir: Path, extraction_summary: Path,
            satellite_dir: Path, patch_qa: Path, output_path: Path) -> dict:
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    selected, provenance = load_verified_rows(association_dir, icl_dir, extraction_summary)
    icl_dir = Path(icl_dir).resolve()
    satellite_dir = Path(satellite_dir).resolve()
    annulus_path = icl_dir / "annulus_pair_profiles.tsv"
    parent = [row for row in read_tsv(annulus_path)
              if row["radial_bin_kpc"] == "200-300"]
    satellite_summary_path = satellite_dir / "summary.json"
    satellite_table_path = satellite_dir / "cluster_profiles.tsv"
    satellite_summary = json.loads(satellite_summary_path.read_text(encoding="utf-8"))
    if (satellite_summary.get("schema") != "hsc-satellite-outer-profile-per-annulus-v1"
            or satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("n_frames_completed") != 800
            or satellite_summary.get("output_sha256", {}).get(str(satellite_table_path)) != sha256(satellite_table_path)):
        raise ValueError("unexpected or hash-mismatched full-sample satellite profile input")
    satellite_rows = [row for row in read_tsv(satellite_table_path)
                      if row["cohort"] == "full_sample" and row["radial_bin_kpc"] == "50-75"]
    satellite_keys = {(str(row["cluster_id"]), str(row["filter"]))
                      for row in satellite_rows}
    if len(satellite_keys) != len(satellite_rows):
        raise ValueError("duplicate 50-75 kpc satellite cluster/filter profile")
    if len(parent) != 224 or len(selected) != 94:
        raise ValueError(f"support sample changed: expected 224 parent and 94 association rows; got {len(parent)}, {len(selected)}")
    if len({(str(row["id"]), row["filter"]) for row in parent}) != len(parent):
        raise ValueError("duplicate 200-300 kpc parent key")
    selected_keys = {(str(row["cluster_id"]), str(row["filter"])) for row in selected}
    if any(band not in ("i", "i2") for _, band in selected_keys):
        raise ValueError("unexpected selected filter")
    selected_parent = [row for row in parent
                       if (str(row["id"]), str(row["filter"])) in selected_keys]
    if len(selected_parent) != len(selected):
        raise ValueError("association rows do not map one-to-one onto ICL support parent")
    icl_keys = {(str(row["id"]), str(row["filter"])) for row in parent}
    support_overlap = summarize_support_overlap(icl_keys, satellite_keys, selected_keys)
    expected_overlap_counts = {
        "i": {"icl_supported": 71, "satellite_profile_supported": 42,
              "both_supported": 26, "icl_supported_without_satellite_profile": 45,
              "satellite_profile_without_icl_support": 16, "association_rows": 26},
        "i2": {"icl_supported": 153, "satellite_profile_supported": 140,
               "both_supported": 68, "icl_supported_without_satellite_profile": 85,
               "satellite_profile_without_icl_support": 72, "association_rows": 68},
    }
    for band, expected in expected_overlap_counts.items():
        if any(support_overlap[band][key] != value for key, value in expected.items()):
            raise ValueError(f"fixed {band} support intersection changed")
    by_filter = {}
    expected_by_filter = {"i": (71, 26), "i2": (153, 68)}
    for band in ("i", "i2"):
        parent_band = [row for row in parent if row["filter"] == band]
        selected_band = [row for row in selected_parent if row["filter"] == band]
        if (len(parent_band), len(selected_band)) != expected_by_filter[band]:
            raise ValueError(f"fixed {band} parent/association counts changed")
        by_filter[band] = summarize_groups(parent_band, selected_band)

    keys, centers = tract_centers_from_patch_qa(Path(patch_qa))
    region_labels = tract_region_labels(centers, 2.0)
    ra = [float(row["ra_deg"]) for row in parent]
    dec = [float(row["dec_deg"]) for row in parent]
    assigned = nearest_region_labels(ra, dec, centers, region_labels)
    region_id = {(str(row["id"]), str(row["filter"])): int(label)
                 for row, label in zip(parent, assigned)}
    region_table = {}
    for band in ("i", "i2"):
        region_table[band] = {}
        for label in sorted(set(int(v) for v in region_labels)):
            in_region = [row for row in parent if row["filter"] == band
                         and region_id[(str(row["id"]), str(row["filter"]))] == label]
            selected_n = sum((str(row["id"]), str(row["filter"])) in selected_keys
                             for row in in_region)
            region_table[band][str(label)] = {
                "n_icl_support_rows": len(in_region),
                "n_association_rows": selected_n,
                "association_fraction": selected_n / len(in_region) if in_region else None,
            }
    output = {
        "schema": "astrafek-icl-satellite-association-selection-audit-v1",
        "status": "complete_exploratory_selection_audit",
        "science_measurement": False,
        "method": {
            "parent": "hash-verified 200-300 kpc target-control ICL support rows",
            "selected": "94 verified cluster/filter rows with satellite outer profiles",
            "numeric_balance": "descriptive standardized mean differences; no inferential p-values",
            "region_labels": "nearest patch-QA tract-center connected component at 2 degrees; approximate regional composition diagnostic",
        },
        "n_patch_qa_tracts": len(keys),
        "n_patch_qa_components": int(len(np.unique(region_labels))),
        "by_filter": by_filter,
        "support_overlap": support_overlap,
        "approximate_region_counts": region_table,
        "provenance": provenance,
        "inputs": {
            "icl_annulus_table": str(annulus_path),
            "icl_annulus_table_sha256": sha256(annulus_path),
            "icl_summary_sha256": sha256(icl_dir / "summary.json"),
            "satellite_summary": str(satellite_summary_path),
            "satellite_summary_sha256": sha256(satellite_summary_path),
            "satellite_profile_table": str(satellite_table_path),
            "satellite_profile_table_sha256": sha256(satellite_table_path),
            "patch_qa": str(Path(patch_qa).resolve()),
            "patch_qa_sha256": sha256(Path(patch_qa)),
        },
        "limitations": [
            "only 200-300 kpc ICL-supported rows define this comparison parent; it is not all 800 clusters",
            "standardized differences describe selection and do not correct it",
            "patch-center regions are approximate and are not survey polygons",
            "no images were read and no new HSC data were queried or downloaded",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(output, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--satellite", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.association, args.icl, args.extraction_summary,
                     args.satellite, args.patch_qa, args.output)
    print(json.dumps({"status": result["status"],
        "science_measurement": result["science_measurement"],
        "by_filter": {band: {key: value for key, value in summary.items()
                            if key != "numeric"}
                      for band, summary in result["by_filter"].items()},
        "approximate_region_counts": result["approximate_region_counts"],
        "support_overlap": result["support_overlap"]},
        ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
