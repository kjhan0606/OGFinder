#!/usr/bin/env python3
"""Paired sensitivity comparison of PSF-supported 40% vs 60% ICL-satellite profiles.

This compares already-completed injection/re-detection runs on identical frames.
It is a selection/support diagnostic, not an astrophysical detection or correction.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

METRICS = {
    "n_satellites": ("measurements_by_scenario", "baseline", "satellite_support", "n_satellites"),
    "n_control_profiles": ("measurements_by_scenario", "baseline", "satellite_support", "n_control_profiles"),
    "satellite_unit_response": ("transfer_response", "satellite_response_to_satellite_unit_injection"),
    "satellite_envelopes_to_icl": ("transfer_response", "icl_response_to_all_satellite_envelopes"),
    "diffuse_unit_to_icl": ("transfer_response", "icl_response_to_diffuse_unit_injection"),
    "diffuse_to_satellite": ("transfer_response", "satellite_response_to_smooth_diffuse"),
    "icl_plane_corrected_diffuse_to_icl": ("transfer_response", "icl_plane_corrected_response_to_diffuse_unit_injection"),
    "icl_plane_corrected_satellites_to_icl": ("transfer_response", "icl_plane_corrected_response_to_all_satellite_envelopes"),
    "satellite_delta_baseline": ("measurements_by_scenario", "baseline", "satellite_delta_50_75_njy_arcsec2"),
    "satellite_delta_satellite_only": ("measurements_by_scenario", "satellite_only", "satellite_delta_50_75_njy_arcsec2"),
    "satellite_delta_both": ("measurements_by_scenario", "both", "satellite_delta_50_75_njy_arcsec2"),
}


def get(d, path):
    for k in path:
        d = d[k]
    return d


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_run(folder, expected_floor):
    rows = {}
    for path in sorted(folder.glob("cluster_*.json")):
        d = json.loads(path.read_text())
        if d.get("status") != "complete" or d.get("science_measurement") is not False:
            raise ValueError(f"Unexpected status/science_measurement in {path}")
        floor = d.get("satellite_coverage_floor")
        if expected_floor is not None and floor is not None and not np.isclose(floor, expected_floor):
            raise ValueError(f"Wrong coverage floor in {path}: {floor}")
        key = (str(d["cluster_id"]), d["filter"])
        if key in rows:
            raise ValueError(f"Duplicate frame key: {key}")
        rows[key] = (path, d)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--low-support", type=Path, required=True, help="40% run directory")
    p.add_argument("--high-support", type=Path, required=True, help="60% run directory")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--high-floor", type=float, default=0.6,
                   help="documented support floor for the legacy run (its JSON predates this field)")
    p.add_argument("--bootstrap", type=int, default=10000)
    a = p.parse_args()
    low = load_run(a.low_support, 0.4)
    high = load_run(a.high_support, 0.6)
    keys = sorted(set(low) & set(high))
    if set(low) != set(high):
        raise ValueError(f"Frame sets differ: 40% only {len(set(low)-set(high))}; 60% only {len(set(high)-set(low))}")
    if not keys:
        raise ValueError("No paired frames")

    records = []
    for key in keys:
        lp, l = low[key]
        hp, h = high[key]
        # The 100 kpc diffuse profile, satellite profile, PSF convolution, and
        # re-detection settings must agree. Paths may differ across run dates;
        # compare content hashes for image and PSF inputs instead.
        for field in ("filter", "z", "icl_profile_bin_kpc", "satellite_profile_bin_kpc", "cluster_radius_selection_kpc"):
            if l.get(field) != h.get(field):
                raise ValueError(f"Incompatible {field} at {key}: {l.get(field)} != {h.get(field)}")
        for field in ("diffuse_scale_kpc", "satellite_scale_kpc", "satellite_amplitude_at_50_per_galaxy_njy_arcsec2", "diffuse_amplitude_at_100_njy_arcsec2"):
            if l["injection"].get(field) != h["injection"].get(field):
                raise ValueError(f"Incompatible injection {field} at {key}")
        li, hi = l["input_sha256"], h["input_sha256"]
        l_image = next((v for k, v in li.items() if k.endswith(f"/images/cluster_{key[0]}.fits")), None)
        h_image = next((v for k, v in hi.items() if k.endswith(f"/images/cluster_{key[0]}.fits")), None)
        l_psf = sorted(v for k, v in li.items() if "/psf_" in k and k.endswith("_HSC-I.fits"))
        h_psf = sorted(v for k, v in hi.items() if "/psf_" in k and k.endswith("_HSC-I.fits"))
        if l_image != h_image or l_psf != h_psf:
            raise ValueError(f"Image/PSF content differs at {key}")
        rec = {"cluster_id": key[0], "filter": key[1], "z": l["z"],
               "image_sha256": l_image, "n_position_psfs": len(l_psf),
               "coverage_floor_40": l.get("satellite_coverage_floor"),
               "coverage_floor_60": a.high_floor,
               "n_recovered_members_40": l["support"]["n_matched_members"],
               "n_recovered_members_60": h["support"]["n_matched_members"],
               "n_targets_40": l["support"]["n_injected_targets_in_50_500_kpc"],
               "n_targets_60": h["support"]["n_injected_targets_in_50_500_kpc"]}
        for name, path in METRICS.items():
            rec[f"{name}_40"] = get(l, path)
            rec[f"{name}_60"] = get(h, path)
            rec[f"{name}_delta_40_minus_60"] = rec[f"{name}_40"] - rec[f"{name}_60"]
        # Pixel-support and recovered-source metrics should not change merely
        # because the outer-profile support floor changed.
        for scenario in ("baseline", "icl_only", "satellite_only", "both"):
            lr = l["redetection_by_scenario"][scenario]
            hr = h["redetection_by_scenario"][scenario]
            for field in ("n_matched_members", "n_targets_in_50_500_kpc", "n_detections", "icl_mask_pixels"):
                if lr.get(field) != hr.get(field):
                    raise ValueError(f"Redetection mismatch {scenario}.{field} at {key}")
        records.append(rec)

    a.output_dir.mkdir(parents=True, exist_ok=False)
    fields = list(records[0])
    pair_path = a.output_dir / "paired_frames.tsv"
    with pair_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        w.writeheader(); w.writerows(records)

    rng = np.random.default_rng(20261105)
    cluster_ids = sorted({r["cluster_id"] for r in records})
    nboot = max(0, a.bootstrap)
    summary = {"schema": "icl-satellite-support-threshold-comparison-v1",
               "science_measurement": False,
               "interpretation": "paired selection/support sensitivity on same image and PSF contents; not a completeness correction or astrophysical detection",
               "low_support_run": str(a.low_support), "high_support_run": str(a.high_support),
               "paired_frames": len(records), "paired_cluster_ids": len(cluster_ids),
               "filters": {f: sum(r["filter"] == f for r in records) for f in sorted({r["filter"] for r in records})},
               "image_and_psf_hashes_verified_equal": True,
               "re_detection_outputs_verified_equal": True,
               "metrics": {}}
    for metric in METRICS:
        summary["metrics"][metric] = {}
        for filt in sorted({r["filter"] for r in records} | {"all"}):
            rr = [r for r in records if filt == "all" or r["filter"] == filt]
            x40 = np.asarray([r[f"{metric}_40"] for r in rr], dtype=float)
            x60 = np.asarray([r[f"{metric}_60"] for r in rr], dtype=float)
            dd = x40 - x60
            unique = sorted({r["cluster_id"] for r in rr})
            boots = []
            if nboot:
                by_cluster = {cid: np.asarray([r[f"{metric}_delta_40_minus_60"] for r in rr if r["cluster_id"] == cid], dtype=float) for cid in unique}
                for _ in range(nboot):
                    sampled = rng.choice(unique, len(unique), replace=True)
                    vals = np.concatenate([by_cluster[c] for c in sampled])
                    boots.append(float(np.median(vals)))
            summary["metrics"][metric][filt] = {
                "n_frames": len(rr), "median_40": float(np.median(x40)), "q16_40": float(np.quantile(x40, .16)), "q84_40": float(np.quantile(x40, .84)),
                "median_60": float(np.median(x60)), "q16_60": float(np.quantile(x60, .16)), "q84_60": float(np.quantile(x60, .84)),
                "median_paired_delta_40_minus_60": float(np.median(dd)),
                "paired_cluster_bootstrap_95pct_ci": [float(np.quantile(boots, .025)), float(np.quantile(boots, .975))] if boots else None,
                "mean_paired_delta_40_minus_60": float(np.mean(dd))}
    summary["inputs"] = {
        "low_run_json_sha256": {f"cluster_{cid}_{filt}.json": sha256(low[(cid, filt)][0]) for cid, filt in keys},
        "high_run_json_sha256": {f"cluster_{cid}_{filt}.json": sha256(high[(cid, filt)][0]) for cid, filt in keys}}
    (a.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (a.output_dir / "assessment.md").write_text(
        "# ICL–satellite profile support threshold comparison\n\n"
        f"Paired {len(records)} frame results across {len(cluster_ids)} cluster IDs compare the 40% run (`{a.low_support}`) with the 60% run (`{a.high_support}`). Image and PSF file content hashes match within every pair. Re-detection counts and masks are identical across runs.\n\n"
        "The support threshold changes which satellite profiles/control azimuths enter the estimator. Reported differences therefore quantify selection/estimator sensitivity for this PSF-supported subset; they are not a population-completeness correction, confidence interval for an astrophysical effect, or detection claim. The injection records state `science_measurement=false`. Frame quantiles describe sample spread. Bootstrap intervals resample cluster IDs and quantify only paired sample uncertainty for this diagnostic.\n\n"
        "See `summary.json` for filter-specific paired medians, frame spread, and cluster-bootstrap intervals; `paired_frames.tsv` contains the auditable per-frame differences.\n")
    print(json.dumps({"output_dir": str(a.output_dir), "paired_frames": len(records), "paired_cluster_ids": len(cluster_ids),
                      "output_sha256": {p.name: sha256(p) for p in a.output_dir.iterdir() if p.is_file()}}, indent=2))

if __name__ == "__main__":
    main()
