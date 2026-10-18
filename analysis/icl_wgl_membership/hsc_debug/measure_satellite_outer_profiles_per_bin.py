#!/usr/bin/env python3
"""Measure exploratory satellite-centered profiles with annulus-specific support.

Unlike the complete-profile estimator, this retains a target independently in
each radial annulus when the target and at least the configured number of
azimuthal controls pass the same coverage floor in that annulus. This avoids
discarding all valid radii because one annulus is masked, while recording the
changing sample support explicitly. It does not download data or certify a
science measurement.
"""
import argparse
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import measure_satellite_outer_profiles as core
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT


PROFILE_FIELDS = [
    "cluster_id", "filter", "redshift", "richness", "cohort", "member_index",
    "detection_number", "member_detection_separation_arcsec", "cluster_radius_kpc",
    "radial_bin_kpc", "n_controls", "target_sb_njy_arcsec2",
    "control_sb_njy_arcsec2", "delta_sb_njy_arcsec2",
    "target_err_nominal_njy_arcsec2", "target_coverage", "control_coverage",
    "target_n_pixels",
]
DIAGNOSTIC_FIELDS = [
    "cluster_id", "filter", "cohort", "member_index", "detection_number",
    "cluster_radius_kpc", "member_detection_separation_arcsec",
    "target_min_coverage", "target_coverage_by_bin", "candidate_control_positions",
    "valid_controls_by_bin", "n_annuli_retained", "status",
]


def summarize_per_bin(rows, work, replicas=2000, seed=20261018):
    """Summarize satellite rows within clusters, then equally weight clusters."""
    bin_names = [f"{lo:g}-{hi:g}" for lo, hi in
                 zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])]
    grouped = {}
    for row in rows:
        key = (str(row["cluster_id"]), row["filter"], row["cohort"],
               row["radial_bin_kpc"])
        grouped.setdefault(key, []).append(row)

    cluster_rows = []
    for (cluster_id, band, cohort, radial_bin), group in sorted(grouped.items()):
        cluster_rows.append({
            "cluster_id": cluster_id, "filter": band, "cohort": cohort,
            "radial_bin_kpc": radial_bin, "n_satellites": len(group),
            "mean_target_sb_njy_arcsec2": float(np.mean(
                [float(row["target_sb_njy_arcsec2"]) for row in group])),
            "mean_control_sb_njy_arcsec2": float(np.mean(
                [float(row["control_sb_njy_arcsec2"]) for row in group])),
            "mean_delta_sb_njy_arcsec2": float(np.mean(
                [float(row["delta_sb_njy_arcsec2"]) for row in group])),
        })
    cluster_fields = ["cluster_id", "filter", "cohort", "radial_bin_kpc",
                      "n_satellites", "mean_target_sb_njy_arcsec2",
                      "mean_control_sb_njy_arcsec2", "mean_delta_sb_njy_arcsec2"]
    core.write_tsv(work / "cluster_profiles.tsv", cluster_rows, cluster_fields)

    filter_summary = {}
    filter_rows = []
    cohorts = sorted({row["cohort"] for row in cluster_rows})
    for cohort_index, cohort in enumerate(cohorts):
        filter_summary[cohort] = {}
        for band_index, band in enumerate(("i", "i2")):
            band_clusters = [row for row in cluster_rows
                             if row["cohort"] == cohort and row["filter"] == band]
            if not band_clusters:
                continue
            filter_summary[cohort][band] = {"bins": {}}
            for bin_index, radial_bin in enumerate(bin_names):
                selected = [row for row in band_clusters
                            if row["radial_bin_kpc"] == radial_bin]
                if not selected:
                    continue
                cluster_values = np.asarray(
                    [row["mean_delta_sb_njy_arcsec2"] for row in selected], dtype=float)
                rng = np.random.default_rng(
                    seed + cohort_index * 100000 + band_index * 10000 + bin_index)
                boot = np.empty(replicas, dtype=float)
                for replica in range(replicas):
                    draw = rng.integers(0, len(cluster_values), size=len(cluster_values))
                    boot[replica] = np.mean(cluster_values[draw])
                low, high = np.percentile(boot, (2.5, 97.5))
                n_satellites = sum(int(row["n_satellites"]) for row in selected)
                record = {
                    "n_clusters": len(selected), "n_satellites": n_satellites,
                    "mean_delta": float(np.mean(cluster_values)),
                    "bootstrap_low": float(low), "bootstrap_high": float(high),
                }
                filter_summary[cohort][band]["bins"][radial_bin] = record
                filter_rows.append({
                    "cohort": cohort, "filter": band, "radial_bin_kpc": radial_bin,
                    "n_clusters": record["n_clusters"],
                    "n_satellites": n_satellites,
                    "equal_cluster_mean_delta_sb_njy_arcsec2": record["mean_delta"],
                    "cluster_bootstrap_95_low": record["bootstrap_low"],
                    "cluster_bootstrap_95_high": record["bootstrap_high"],
                })
    filter_fields = ["cohort", "filter", "radial_bin_kpc", "n_clusters",
                     "n_satellites", "equal_cluster_mean_delta_sb_njy_arcsec2",
                     "cluster_bootstrap_95_low", "cluster_bootstrap_95_high"]
    core.write_tsv(work / "filter_profiles.tsv", filter_rows, filter_fields)
    return filter_summary


def plot_profiles(filter_summary, work):
    from matplotlib import pyplot as plt

    bin_names = [f"{lo:g}-{hi:g}" for lo, hi in
                 zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])]
    edges = core.R_EDGES_KPC
    centers = np.sqrt(edges[:-1] * edges[1:])
    xerr = np.vstack([centers - edges[:-1], edges[1:] - centers])
    cohorts = sorted(filter_summary)
    fig, axes = plt.subplots(max(1, len(cohorts)), 2,
                             figsize=(11, max(4, 3.8 * len(cohorts))),
                             sharex=True, squeeze=False)
    for row_index, cohort in enumerate(cohorts):
        for col_index, band in enumerate(("i", "i2")):
            ax = axes[row_index, col_index]
            record = filter_summary.get(cohort, {}).get(band, {}).get("bins", {})
            points = [(index, record[name]) for index, name in enumerate(bin_names)
                      if name in record]
            ax.axhline(0., color="0.4", lw=1, ls="--")
            ax.set_xscale("log")
            ax.set_yscale("symlog", linthresh=1.)
            ax.grid(True, which="both", alpha=.2)
            ax.set_title(f"{cohort}; {band}")
            if points:
                indices = np.asarray([index for index, _ in points], dtype=int)
                means = np.asarray([item["mean_delta"] for _, item in points])
                lows = np.asarray([item["bootstrap_low"] for _, item in points])
                highs = np.asarray([item["bootstrap_high"] for _, item in points])
                yerr = np.vstack([means - lows, highs - means])
                ax.errorbar(centers[indices], means, xerr=xerr[:, indices], yerr=yerr,
                            fmt="o", capsize=2, color="#245a87")
                for index, item in points:
                    ax.annotate(f"{item['n_clusters']} cl",
                                (centers[index], item["mean_delta"]),
                                xytext=(3, 5), textcoords="offset points", fontsize=7)
            else:
                ax.text(.5, .5, "no supported annuli", ha="center", va="center",
                        transform=ax.transAxes)
            ax.set_xlabel("Satellite-centric radius (kpc)")
            ax.set_ylabel(r"Target − control (nJy arcsec$^{-2}$)")
    fig.suptitle("HSC satellite outer-light annuli with independent support\n"
                 "Exploratory only; sample changes by radius")
    fig.tight_layout()
    fig.savefig(work / "satellite_profile_review.png", dpi=160)
    plt.close(fig)


def process_frame(frame, prep_dir, cosmology, coverage_floor, min_controls):
    ident = str(frame["id"])
    prep_path = prep_dir / f"cluster_{ident}"
    image_path = Path(frame["image"])
    detections_path, member_path = prep_path / "detections.tsv", prep_path / "members.tsv"
    mask_path = prep_path / "masks.fits.gz"
    prep_json_path = prep_path / "preparation.json"
    required = [image_path, detections_path, member_path, mask_path, prep_json_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input(s): " + ", ".join(missing))
    prep_record = core.read_json(prep_json_path)
    if prep_record.get("science_measurement") is not False:
        raise ValueError(f"Unexpected preparation gate for {ident}")
    info = core.validate_fits(image_path, frame)
    if info["shape"] != tuple(prep_record.get("shape", info["shape"])):
        raise ValueError(f"Image shape/provenance mismatch for {ident}")
    if info["coadd_header_filter"] not in ("i", "i2"):
        raise ValueError(f"Unexpected filter for {ident}: {info['coadd_header_filter']}")
    detections = core.parse_detections(core.read_tsv(detections_path))
    member_rows = core.read_tsv(member_path)
    matched, match_stats = core.match_members(member_rows, detections,
                                               info["pixel_scale_arcsec"])
    if not detections:
        raise ValueError(f"No detections available for {ident}")

    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
    with fits.open(mask_path, memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if not (image.shape == variance.shape == bits.shape == segmentation.shape == base.shape):
        raise ValueError(f"Image/variance/mask shapes differ for {ident}")
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    bad_bits = sum(1 << bit for bit in BAD_BIT_POSITIONS) | BRIGHT_OBJECT_BIT
    if not np.array_equal(base, ((bits & bad_bits) != 0) | invalid):
        raise ValueError(f"Saved base mask no longer matches original quality bits for {ident}")

    in_frame_rows = [row for row in member_rows if row["IN_FRAME"] == "True"]
    all_member_xy = np.asarray([[float(row["X_IMAGE"]) - 1.,
                                 float(row["Y_IMAGE"]) - 1.] for row in in_frame_rows], dtype=float)
    matched_indices = {item["member_index"] for item in matched}
    fallback_xy = np.asarray([[float(row["X_IMAGE"]) - 1.,
                               float(row["Y_IMAGE"]) - 1.] for row in in_frame_rows
                              if int(row["INDEX"]) not in matched_indices], dtype=float)
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(all_member_xy) if len(all_member_xy) else None
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    pixel_scale = float(info["pixel_scale_arcsec"])
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(
        float(frame["z"])).value * pixel_scale / 60.)
    conversion = float(core.nJy_factor(info))
    dets = core.detection_arrays(detections, pixel_scale)

    selected = []
    for item in matched:
        det = detections[item["detection_index"]]
        radius = math.hypot(det["x"] - center_xy[0], det["y"] - center_xy[1]) * kpc_per_pixel
        if core.MIN_CLUSTER_RADIUS_KPC <= radius < core.MAX_CLUSTER_RADIUS_KPC:
            selected.append((item, det, radius))

    rows, diagnostics = [], []
    targets_with_bin = controls_attempted = 0
    for item, target_det, cluster_radius in selected:
        diagnostic = {
            "cluster_id": ident, "filter": info["coadd_header_filter"],
            "cohort": core.cohort_name(frame), "member_index": item["member_index"],
            "detection_number": target_det["NUMBER"],
            "cluster_radius_kpc": cluster_radius,
            "member_detection_separation_arcsec": item["match_arcsec"],
            "target_min_coverage": None, "target_coverage_by_bin": "",
            "candidate_control_positions": 0, "valid_controls_by_bin": "",
            "n_annuli_retained": 0, "status": None,
        }
        detection_index = item["detection_index"]
        if not dets["shape_valid"][detection_index]:
            diagnostic["status"] = "invalid_target_shape"
            diagnostics.append(diagnostic)
            continue
        target_xy = np.asarray([target_det["x"], target_det["y"]])
        target = core.profile_at(
            image, variance, base, segmentation, dets, fallback_tree, target_xy,
            pixel_scale, kpc_per_pixel, conversion, target_det["NUMBER"], coverage_floor)
        if target is None:
            diagnostic["status"] = "target_stamp_unavailable"
            diagnostics.append(diagnostic)
            continue
        coverage = target["coverage"]
        target_good = coverage >= coverage_floor
        diagnostic["target_min_coverage"] = float(np.min(coverage))
        diagnostic["target_coverage_by_bin"] = ",".join(f"{value:.3f}" for value in coverage)
        if not np.any(target_good):
            diagnostic["status"] = "no_usable_annulus"
            diagnostics.append(diagnostic)
            continue
        targets_with_bin += 1

        points = core.control_centers(
            target_xy, center_xy, cluster_radius / kpc_per_pixel, dets, member_tree,
            pixel_scale, kpc_per_pixel, image.shape)
        diagnostic["candidate_control_positions"] = len(points)
        controls = []
        for point in points:
            control = core.profile_at(
                image, variance, base, segmentation, dets, fallback_tree, point,
                pixel_scale, kpc_per_pixel, conversion, None, coverage_floor)
            if control is not None:
                controls.append(control)
        controls_attempted += len(controls)

        valid_by_bin = [sum(control["coverage"][j] >= coverage_floor for control in controls)
                        for j in range(len(core.R_EDGES_KPC) - 1)]
        diagnostic["valid_controls_by_bin"] = ",".join(str(int(value)) for value in valid_by_bin)
        for j, (lo, hi) in enumerate(zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])):
            if not target_good[j]:
                continue
            selected_controls = [control for control in controls
                                 if control["coverage"][j] >= coverage_floor]
            if len(selected_controls) < min_controls:
                continue
            control_sb = float(np.mean([control["sb"][j] for control in selected_controls]))
            control_coverage = float(np.mean([control["coverage"][j]
                                              for control in selected_controls]))
            target_sb = float(target["sb"][j])
            rows.append({
                "cluster_id": ident, "filter": info["coadd_header_filter"],
                "redshift": float(frame["z"]), "richness": float(frame["richness"]),
                "cohort": core.cohort_name(frame), "member_index": item["member_index"],
                "detection_number": target_det["NUMBER"],
                "member_detection_separation_arcsec": item["match_arcsec"],
                "cluster_radius_kpc": cluster_radius,
                "radial_bin_kpc": f"{lo:g}-{hi:g}",
                "n_controls": len(selected_controls),
                "target_sb_njy_arcsec2": target_sb,
                "control_sb_njy_arcsec2": control_sb,
                "delta_sb_njy_arcsec2": target_sb - control_sb,
                "target_err_nominal_njy_arcsec2": float(target["err"][j]),
                "target_coverage": float(coverage[j]),
                "control_coverage": control_coverage,
                "target_n_pixels": int(target["n_pixels"][j]),
            })
        diagnostic["n_annuli_retained"] = sum(
            row["cluster_id"] == ident and row["member_index"] == item["member_index"]
            for row in rows)
        diagnostic["status"] = ("retained" if diagnostic["n_annuli_retained"]
                                 else "insufficient_bin_controls")
        diagnostics.append(diagnostic)

    frame_record = {
        "cluster_id": ident, "filter": info["coadd_header_filter"],
        "redshift": float(frame["z"]), "cohort": core.cohort_name(frame),
        "n_member_positions": len(member_rows),
        "n_member_positions_in_frame": match_stats["in_frame"],
        "n_unique_positional_matches_le_1arcsec": len(matched),
        "n_ambiguous_matches_rejected": match_stats["ambiguous"],
        "n_duplicate_detection_assignments_rejected": match_stats["duplicate_assignments"],
        "n_matches_50_500kpc": len(selected),
        "n_targets_with_any_usable_annulus": targets_with_bin,
        "n_satellite_annulus_profiles": len(rows),
        "n_control_positions_attempted": controls_attempted,
        "image_bytes": image_path.stat().st_size, "image_sha256": core.sha256(image_path),
        "mask_bytes": mask_path.stat().st_size, "mask_sha256": core.sha256(mask_path),
        "detection_table_sha256": core.sha256(detections_path),
        "member_table_sha256": core.sha256(member_path),
    }
    return rows, frame_record, diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--prep-dir", type=Path)
    parser.add_argument("--artifact-audit", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--limit", type=int,
                        help="Optional leading subset for an engineering run")
    parser.add_argument("--coverage-floor", type=float, default=.60)
    parser.add_argument("--min-controls", type=int, default=3)
    parser.add_argument("--bootstrap-replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run CPU image processing on lageunha")
    if (not 0 < args.coverage_floor <= 1 or args.min_controls < 1
            or args.bootstrap_replicas < 100 or (args.limit is not None and args.limit < 1)):
        parser.error("Invalid coverage/control/bootstrap configuration")

    data, work = args.data.resolve(), args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    manifest_path = (args.manifest or data / "sample_manifest.json").resolve()
    prep_dir = (args.prep_dir or data / "satellite_detection_psf_2026-10-18_01").resolve()
    prep_summary_path = prep_dir / "summary.json"
    if not manifest_path.is_file() or not prep_summary_path.is_file():
        parser.error("Input manifest or detection/mask preparation summary is missing")
    manifest, prep_summary = core.read_json(manifest_path), core.read_json(prep_summary_path)
    if (manifest.get("science_measurement") is not False
            or prep_summary.get("science_measurement") is not False):
        parser.error("Expected explicitly exploratory, non-science inputs")
    all_frames = manifest.get("clusters", [])
    all_ids = [str(frame.get("id")) for frame in all_frames]
    prep_ids = [str(record.get("id")) for record in prep_summary.get("records", [])]
    if (not all_frames or len(set(all_ids)) != len(all_ids)
            or len(set(prep_ids)) != len(prep_ids)
            or prep_summary.get("status") != "complete"
            or int(prep_summary.get("n_prepared", -1)) != len(prep_ids)
            or not set(all_ids).issubset(set(prep_ids))):
        parser.error("Manifest and prepared cluster IDs do not form a complete, unique sample")
    frames = all_frames[:args.limit] if args.limit else all_frames
    if not frames:
        parser.error("No frames selected")
    frames = core.resolve_image_paths(frames, manifest_path)
    audit_path = (args.artifact_audit.resolve() if args.artifact_audit
                  else prep_dir / "artifact_audit.json")
    audit = None
    if len(all_frames) > 50:
        if not audit_path.is_file():
            parser.error(f"Full-sample measurement requires an artifact audit: {audit_path}")
        audit = core.read_json(audit_path)
        if (audit.get("passed") is not True or audit.get("science_measurement") is not False
                or audit.get("n_frames") != len(all_frames)
                or audit.get("input_manifest_sha256") != core.sha256(manifest_path)
                or audit.get("prepare_script_sha256") != prep_summary.get("preparation_script_sha256")):
            parser.error("Full-sample artifact audit does not match the selected manifest/preparation")

    started = time.time()
    cosmology = FlatLambdaCDM(H0=float(manifest.get("h0", 70.)),
                              Om0=float(manifest.get("omega_m", .3)), Tcmb0=0 * u.K)
    all_rows, frame_records, all_diagnostics, errors = [], [], [], []
    for ordinal, frame in enumerate(frames, 1):
        print(f"[{ordinal}/{len(frames)}] measuring cluster {frame['id']}", flush=True)
        try:
            rows, record, diagnostics = process_frame(
                frame, prep_dir, cosmology, args.coverage_floor, args.min_controls)
            all_rows.extend(rows)
            frame_records.append(record)
            all_diagnostics.extend(diagnostics)
            print(f"  {record['filter']}: targets with usable annulus="
                  f"{record['n_targets_with_any_usable_annulus']}, "
                  f"satellite-annulus rows={record['n_satellite_annulus_profiles']}", flush=True)
        except Exception as exc:
            errors.append({"cluster_id": str(frame.get("id")),
                           "error_type": type(exc).__name__, "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
    if not frame_records:
        parser.error("No frame completed; no output directory created")

    work.mkdir(parents=True, exist_ok=False)
    core.write_tsv(work / "satellite_profiles.tsv", all_rows, PROFILE_FIELDS)
    core.write_tsv(work / "candidate_diagnostics.tsv", all_diagnostics, DIAGNOSTIC_FIELDS)
    filter_summary = summarize_per_bin(all_rows, work, args.bootstrap_replicas, args.seed)
    plot_profiles(filter_summary, work)

    n_unique_satellites = len({(row["cluster_id"], row["filter"], row["member_index"])
                               for row in all_rows})
    summary = {
        "schema": "hsc-satellite-outer-profile-per-annulus-v1",
        "status": "complete" if not errors else "partial",
        "science_measurement": False,
        "sample_scope": "full_manifest" if len(frames) == len(all_frames) else "manifest_subset",
        "n_frames_requested": len(frames), "n_frames_completed": len(frame_records),
        "n_frame_errors": len(errors), "n_satellites_with_any_bin": n_unique_satellites,
        "n_satellite_annulus_profiles": len(all_rows),
        "n_candidate_diagnostics": len(all_diagnostics), "filters": filter_summary,
        "input_manifest": str(manifest_path), "preparation_directory": str(prep_dir),
        "artifact_audit": str(audit_path) if audit else None,
        "artifact_audit_sha256": core.sha256(audit_path) if audit else None,
        "method": {
            "support": "each annulus is retained independently; radial bins may have different satellite and cluster samples",
            "member_selection": "CAMIRA positions matched to a unique saved SEP detection within 1 arcsec; no mass or CAMIRA weight applied",
            "cluster_radius_kpc": [core.MIN_CLUSTER_RADIUS_KPC, core.MAX_CLUSTER_RADIUS_KPC],
            "satellite_radius_bins_kpc": core.R_EDGES_KPC.tolist(),
            "target_masks": "preserve matched SEP label and Kron ellipse; mask all other detected sources, fallback unmatched CAMIRA positions, and original HSC quality/bright-object pixels",
            "control_positions": f"same cluster-centric projected radius; up to {core.N_CONTROLS_MAX} of {core.N_CONTROL_ANGLES} azimuthal rotations; reject positions within {core.CONTROL_SOURCE_CLEARANCE_ARCSEC} arcsec of any detection or {core.CONTROL_MEMBER_CLEARANCE_KPC} kpc of any CAMIRA member",
            "coverage_floor": args.coverage_floor, "minimum_controls_per_annulus": args.min_controls,
            "stack_weighting": "mean satellites within cluster for each annulus, then equal-weight clusters; filters kept separate",
            "flux": "FLUXMAG0 converted to nJy/arcsec^2; no extinction or passband correction",
            "uncertainty": "independent cluster bootstrap for each annulus; nominal pixel variance excludes coadd correlated noise",
        },
        "limitations": [
            "The satellite and cluster samples can change with radius; profile shape can therefore reflect selection changes.",
            "Automatic source masks, BCG centers, PSF wings, sky transfer and coadd covariance are not certified for science use.",
            "CAMIRA member identities are positional only; exact membership weights, stellar masses, colors and photometric redshifts are not used.",
            "The azimuthal controls estimate a local baseline but do not remove asymmetric substructure or establish a satellite/ICL decomposition.",
            "Extinction, passband response, member contamination, completeness, random-field systematics and spatial covariance remain unvalidated.",
        ],
        "manifest_sha256": core.sha256(manifest_path),
        "preparation_summary_sha256": core.sha256(prep_summary_path),
        "script_sha256": core.sha256(Path(__file__).resolve()),
        "frame_records": frame_records, "errors": errors,
        "elapsed_seconds": time.time() - started,
        "bootstrap": {"replicas": args.bootstrap_replicas, "seed": args.seed},
        "output_sha256": {},
    }
    summary_path = work / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    for path in sorted(work.iterdir()):
        if path.is_file() and path.name != "summary.json":
            summary["output_sha256"][str(path)] = core.sha256(path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "frames": len(frame_records),
        "satellites_with_any_bin": n_unique_satellites,
        "satellite_annulus_profiles": len(all_rows),
        "filters": {cohort: {band: {name: {
            "n_clusters": row["n_clusters"], "n_satellites": row["n_satellites"]}
            for name, row in value["bins"].items()} for band, value in bands.items()}
            for cohort, bands in filter_summary.items()},
        "output": str(work),
    }, ensure_ascii=False), flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
