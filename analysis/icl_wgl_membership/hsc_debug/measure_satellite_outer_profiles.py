#!/usr/bin/env python3
"""Exploratory satellite-centered outer-light pilot using existing HSC cutouts.

This reuses the 50-frame validation sample and saved SEP segmentation maps. It
does not download data, infer CAMIRA membership weights, or certify a science
measurement. Each matched member is measured against source-free positions at
the same projected cluster-centric radius to estimate the local diffuse-light
and sky baseline.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import sep
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy import units as u
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from download_pilot import (BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, file_sha256,
                            validate_fits)
from validate_control_images import nJy_factor

R_EDGES_KPC = np.array([5., 10., 15., 20., 30., 40., 50., 75.])
MATCH_ARCSEC = 1.0
KRON_FACTOR = 2.5
MIN_CLUSTER_RADIUS_KPC = 50.
MAX_CLUSTER_RADIUS_KPC = 500.
N_CONTROL_ANGLES = 32
N_CONTROLS_MAX = 8
CONTROL_SOURCE_CLEARANCE_ARCSEC = 3.0
CONTROL_MEMBER_CLEARANCE_KPC = 15.


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path):
    return file_sha256(path)


def parse_detections(rows):
    """Return zero-origin SEP coordinates and the source fields needed for masks."""
    names = ("NUMBER", "X_IMAGE", "Y_IMAGE", "A_IMAGE", "B_IMAGE", "THETA_RAD",
             "KRON_RADIUS", "KRON_FLAGS")
    parsed = []
    for row in rows:
        vals = {name: float(row[name]) for name in names if name != "NUMBER"}
        vals["NUMBER"] = int(row["NUMBER"])
        vals["x"] = vals["X_IMAGE"] - 1.
        vals["y"] = vals["Y_IMAGE"] - 1.
        parsed.append(vals)
    return parsed


def detection_arrays(detections):
    arrays = {name: np.asarray([row[key] for row in detections], dtype=dtype)
              for name, key, dtype in (
                  ("number", "NUMBER", np.int32), ("x", "x", np.float64),
                  ("y", "y", np.float64), ("a", "A_IMAGE", np.float64),
                  ("b", "B_IMAGE", np.float64), ("theta", "THETA_RAD", np.float64),
                  ("kron", "KRON_RADIUS", np.float64), ("kron_flags", "KRON_FLAGS", np.int32))}
    valid_kron = (arrays["kron_flags"] == 0) & np.isfinite(arrays["kron"]) & (arrays["kron"] > 0)
    arrays["kron_scale"] = KRON_FACTOR * np.where(valid_kron, arrays["kron"], 1.2)
    arrays["extent"] = arrays["kron_scale"] * np.maximum(np.maximum(arrays["a"], .5), np.maximum(arrays["b"], .5))
    arrays["tree"] = cKDTree(np.column_stack([arrays["x"], arrays["y"]])) if len(detections) else None
    return arrays


def match_members(member_rows, detections, pixel_scale):
    """Use one-to-one positional matches; reject blends and reused detections."""
    if not detections:
        return [], {"in_frame": 0, "within_match_radius": 0, "ambiguous": 0,
                    "duplicate_assignments": 0}
    positions = np.array([[item["x"], item["y"]] for item in detections])
    tree = cKDTree(positions)
    in_frame = [row for row in member_rows if row["IN_FRAME"] == "True"]
    member_xy = np.array([[float(row["X_IMAGE"]) - 1.,
                           float(row["Y_IMAGE"]) - 1.] for row in in_frame])
    if not len(member_xy):
        return [], {"in_frame": 0, "within_match_radius": 0, "ambiguous": 0,
                    "duplicate_assignments": 0}
    k = min(2, len(detections))
    distance_pix, indices = tree.query(member_xy, k=k)
    if k == 1:
        distance_pix = distance_pix[:, None]
        indices = indices[:, None]
    distance_arcsec = distance_pix * pixel_scale
    candidate = distance_arcsec[:, 0] <= MATCH_ARCSEC
    ambiguous = candidate & (distance_arcsec[:, 1] <= MATCH_ARCSEC) if k > 1 else np.zeros(len(candidate), bool)
    selected = []
    for i, row in enumerate(in_frame):
        if not candidate[i] or ambiguous[i]:
            continue
        det_index = int(indices[i, 0])
        selected.append({"member_index": int(row["INDEX"]),
                         "member_x": float(member_xy[i, 0]),
                         "member_y": float(member_xy[i, 1]),
                         "detection_index": det_index,
                         "detection_number": detections[det_index]["NUMBER"],
                         "match_arcsec": float(distance_arcsec[i, 0])})
    uses = {}
    for item in selected:
        uses.setdefault(item["detection_index"], []).append(item)
    duplicate_indices = {index for index, rows in uses.items() if len(rows) > 1}
    duplicate_assignments = sum(len(uses[index]) - 1 for index in duplicate_indices)
    selected = [item for item in selected if item["detection_index"] not in duplicate_indices]
    return selected, {"in_frame": len(in_frame),
                      "within_match_radius": int(candidate.sum()),
                      "ambiguous": int(ambiguous.sum()),
                      "duplicate_assignments": int(duplicate_assignments)}


def add_circle(mask, x, y, radius):
    ny, nx = mask.shape
    x0 = max(0, int(math.floor(x - radius)))
    x1 = min(nx, int(math.ceil(x + radius)) + 1)
    y0 = max(0, int(math.floor(y - radius)))
    y1 = min(ny, int(math.ceil(y + radius)) + 1)
    if x0 >= x1 or y0 >= y1:
        return
    yy, xx = np.ogrid[y0:y1, x0:x1]
    mask[y0:y1, x0:x1] |= (xx - x) ** 2 + (yy - y) ** 2 <= radius ** 2


def make_stamp_mask(base, segmentation, detections, fallback_tree, center_xy,
                    x0, y0, target_number, pixel_scale, half_size, kron_factor):
    """Mask bad pixels and every neighboring source while preserving one target."""
    cx, cy = map(float, center_xy)
    y1 = min(segmentation.shape[0], int(math.ceil(cy + half_size)) + 1)
    x1 = min(segmentation.shape[1], int(math.ceil(cx + half_size)) + 1)
    mask = np.array(base[y0:y1, x0:x1], dtype=bool, copy=True)
    labels = segmentation[y0:y1, x0:x1]
    if target_number is None:
        mask |= labels > 0
    else:
        mask |= (labels > 0) & (labels != target_number)

    xs, ys = detections["x"], detections["y"]
    near = ((xs + detections["extent"] >= x0) & (xs - detections["extent"] < x1)
            & (ys + detections["extent"] >= y0) & (ys - detections["extent"] < y1))
    if target_number is not None:
        near &= detections["number"] != target_number
    selected = np.flatnonzero(near)
    if len(selected):
        aa = np.maximum(detections["a"][selected], .5)
        bb = np.maximum(detections["b"][selected], .5)
        sep.mask_ellipse(mask, xs[selected] - x0, ys[selected] - y0, aa, bb,
                         detections["theta"][selected],
                         r=detections["kron_scale"][selected])

    fallback_radius = max(2. / pixel_scale, 1.5)
    if fallback_tree is not None:
        near = fallback_tree.query_ball_point(center_xy, half_size + fallback_radius)
        for index in near:
            point = fallback_tree.data[index]
            add_circle(mask, float(point[0] - x0), float(point[1] - y0), fallback_radius)
    return mask


def profile_at(image, variance, base, segmentation, detections, fallback_tree,
               center_xy, pixel_scale, kpc_per_pixel, conversion, target_number,
               coverage_floor, kron_factor=KRON_FACTOR):
    half_size = int(math.ceil(R_EDGES_KPC[-1] / kpc_per_pixel)) + 2
    cx, cy = map(float, center_xy)
    x0, x1 = max(0, int(math.floor(cx - half_size))), min(image.shape[1], int(math.ceil(cx + half_size)) + 1)
    y0, y1 = max(0, int(math.floor(cy - half_size))), min(image.shape[0], int(math.ceil(cy + half_size)) + 1)
    if x0 >= x1 or y0 >= y1:
        return None
    yy, xx = np.indices((y1 - y0, x1 - x0), dtype=np.float32)
    radius_kpc = np.hypot(xx + x0 - cx, yy + y0 - cy) * kpc_per_pixel
    mask = make_stamp_mask(base, segmentation, detections, fallback_tree, (cx, cy),
                           x0, y0, target_number, pixel_scale, half_size, kron_factor)
    data = image[y0:y1, x0:x1]
    var = variance[y0:y1, x0:x1]
    valid_pixels = (~mask & np.isfinite(data) & np.isfinite(var) & (var > 0))
    sb, err, coverage, n_pixels = [], [], [], []
    for lo, hi in zip(R_EDGES_KPC[:-1], R_EDGES_KPC[1:]):
        annulus = (radius_kpc >= lo) & (radius_kpc < hi)
        valid = annulus & valid_pixels
        expected = math.pi * (hi ** 2 - lo ** 2) / kpc_per_pixel ** 2
        count = int(valid.sum())
        coverage.append(min(1., count / expected) if expected > 0 else 0.)
        n_pixels.append(count)
        if count:
            sb.append(float(np.mean(data[valid]) * conversion))
            # Nominal independent-pixel uncertainty only; coadd covariance is not included.
            err.append(float(math.sqrt(float(np.sum(var[valid]))) * conversion / count))
        else:
            sb.append(float("nan")); err.append(float("nan"))
    good = bool(np.all(np.asarray(coverage) >= coverage_floor))
    return {"sb": np.asarray(sb), "err": np.asarray(err),
            "coverage": np.asarray(coverage), "n_pixels": np.asarray(n_pixels),
            "passes_coverage": good}


def control_centers(target_xy, cluster_xy, cluster_radius_pix, detections,
                    member_tree, pixel_scale, kpc_per_pixel, shape):
    dx, dy = target_xy[0] - cluster_xy[0], target_xy[1] - cluster_xy[1]
    angle = math.atan2(dy, dx)
    accepted = []
    for k in range(1, N_CONTROL_ANGLES + 1):
        phi = angle + 2 * math.pi * k / N_CONTROL_ANGLES
        point = np.array([cluster_xy[0] + cluster_radius_pix * math.cos(phi),
                          cluster_xy[1] + cluster_radius_pix * math.sin(phi)])
        if (point[0] < 0 or point[1] < 0 or point[0] >= shape[1] or point[1] >= shape[0]):
            continue
        if detections["tree"] is not None:
            source_distance = float(detections["tree"].query(point)[0]) * pixel_scale
            if source_distance < CONTROL_SOURCE_CLEARANCE_ARCSEC:
                continue
        if member_tree is not None:
            member_distance = float(member_tree.query(point)[0]) * kpc_per_pixel
            if member_distance < CONTROL_MEMBER_CLEARANCE_KPC:
                continue
        accepted.append(point)
        if len(accepted) >= N_CONTROLS_MAX:
            break
    return accepted


def process_frame(frame, prep_dir, cosmo, coverage_floor, min_controls):
    ident = str(frame["id"])
    prep_path = prep_dir / f"cluster_{ident}"
    image_path = Path(frame["image"])
    det_path, member_path = prep_path / "detections.tsv", prep_path / "members.tsv"
    mask_path, prep_json_path = prep_path / "masks.fits.gz", prep_path / "preparation.json"
    required = [image_path, det_path, member_path, mask_path, prep_json_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input(s): " + ", ".join(missing))
    prep_record = read_json(prep_json_path)
    if prep_record.get("science_measurement") is not False:
        raise ValueError(f"Unexpected preparation gate for {ident}")
    info = validate_fits(image_path, frame)
    if info["shape"] != tuple(prep_record.get("shape", info["shape"])):
        raise ValueError(f"Image shape/provenance mismatch for {ident}")
    if info["coadd_header_filter"] not in ("i", "i2"):
        raise ValueError(f"Unexpected filter for {ident}: {info['coadd_header_filter']}")
    detections = parse_detections(read_tsv(det_path))
    member_rows = read_tsv(member_path)
    matched, match_stats = match_members(member_rows, detections, info["pixel_scale_arcsec"])
    if not detections:
        raise ValueError(f"No detections available for {ident}")

    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
    with fits.open(mask_path, memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if image.shape != variance.shape or image.shape != bits.shape or image.shape != segmentation.shape or image.shape != base.shape:
        raise ValueError(f"Image/variance/mask shapes differ for {ident}")
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    bad_bits = sum(1 << bit for bit in BAD_BIT_POSITIONS) | BRIGHT_OBJECT_BIT
    expected_base = ((bits & bad_bits) != 0) | invalid
    if not np.array_equal(base, expected_base):
        raise ValueError(f"Saved base mask no longer matches original quality bits for {ident}")

    all_member_xy = np.array([[float(row["X_IMAGE"]) - 1., float(row["Y_IMAGE"]) - 1.]
                              for row in member_rows if row["IN_FRAME"] == "True"], dtype=float)
    matched_member_indices = {item["member_index"] for item in matched}
    fallback_xy = np.array([[float(row["X_IMAGE"]) - 1., float(row["Y_IMAGE"]) - 1.]
                            for row in member_rows if row["IN_FRAME"] == "True"
                            and int(row["INDEX"]) not in matched_member_indices], dtype=float)
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(all_member_xy) if len(all_member_xy) else None
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    pixel_scale = float(info["pixel_scale_arcsec"])
    kpc_per_pixel = float(cosmo.kpc_proper_per_arcmin(float(frame["z"])).value * pixel_scale / 60.)
    conversion = float(nJy_factor(info))
    dets = detection_arrays(detections)
    selected = []
    for item in matched:
        det = detections[item["detection_index"]]
        x, y = det["x"], det["y"]
        r_cluster_kpc = math.hypot(x - center_xy[0], y - center_xy[1]) * kpc_per_pixel
        if not MIN_CLUSTER_RADIUS_KPC <= r_cluster_kpc < MAX_CLUSTER_RADIUS_KPC:
            continue
        selected.append((item, det, r_cluster_kpc))

    rows, diagnostics = [], []
    n_target_coverage_pass = n_control_coverage_fail = n_controls_selected = 0
    n_min_control_pass = 0
    for item, target_det, r_cluster_kpc in selected:
        target_xy = np.array([target_det["x"], target_det["y"]])
        target = profile_at(image, variance, base, segmentation, dets, fallback_tree,
                            target_xy, pixel_scale, kpc_per_pixel, conversion,
                            target_det["NUMBER"], coverage_floor)
        diag = {"cluster_id": ident, "filter": info["coadd_header_filter"],
                "cohort": str(frame.get("preparation_cohort", "unspecified")),
                "member_index": item["member_index"],
                "detection_number": target_det["NUMBER"],
                "cluster_radius_kpc": r_cluster_kpc,
                "member_detection_separation_arcsec": item["match_arcsec"],
                "target_min_coverage": None, "target_coverage_by_bin": "",
                "candidate_control_positions": 0, "valid_control_profiles": 0,
                "valid_controls_by_bin": "", "status": None}
        if target is None:
            diag["status"] = "target_stamp_unavailable"
            diagnostics.append(diag)
            continue
        diag["target_min_coverage"] = float(np.min(target["coverage"]))
        diag["target_coverage_by_bin"] = ",".join(f"{value:.3f}" for value in target["coverage"])
        if not target["passes_coverage"]:
            diag["status"] = "target_coverage_below_floor"
            diagnostics.append(diag)
            continue
        n_target_coverage_pass += 1
        points = control_centers(target_xy, center_xy,
                                 r_cluster_kpc / kpc_per_pixel, dets,
                                 member_tree, pixel_scale, kpc_per_pixel, image.shape)
        diag["candidate_control_positions"] = len(points)
        controls, attempted_controls = [], []
        for point in points:
            control = profile_at(image, variance, base, segmentation, dets,
                                 fallback_tree, point, pixel_scale, kpc_per_pixel,
                                 conversion, None, coverage_floor)
            if control is not None:
                attempted_controls.append(control)
            if control is not None and control["passes_coverage"]:
                controls.append(control)
            else:
                n_control_coverage_fail += 1
        diag["valid_control_profiles"] = len(controls)
        diag["valid_controls_by_bin"] = ",".join(
            str(int(sum(control["coverage"][j] >= coverage_floor for control in attempted_controls)))
            for j in range(len(R_EDGES_KPC) - 1))
        n_controls_selected += len(controls)
        if len(controls) < min_controls:
            diag["status"] = "insufficient_valid_controls"
            diagnostics.append(diag)
            continue
        n_min_control_pass += 1
        control_sb = np.mean([control["sb"] for control in controls], axis=0)
        control_cov = np.mean([control["coverage"] for control in controls], axis=0)
        delta = target["sb"] - control_sb
        member_index = item["member_index"]
        rows.append({
            "cluster_id": ident, "filter": info["coadd_header_filter"],
            "redshift": float(frame["z"]), "richness": float(frame["richness"]),
            "cohort": str(frame.get("preparation_cohort", "unspecified")),
            "member_index": member_index, "detection_number": target_det["NUMBER"],
            "member_detection_separation_arcsec": item["match_arcsec"],
            "cluster_radius_kpc": r_cluster_kpc,
            "n_controls": len(controls), "kpc_per_pixel": kpc_per_pixel,
            "target_sb": target["sb"].tolist(), "control_sb": control_sb.tolist(),
            "delta_sb": delta.tolist(), "target_err_nominal": target["err"].tolist(),
            "target_coverage": target["coverage"].tolist(),
            "control_coverage": control_cov.tolist(),
            "target_n_pixels": target["n_pixels"].tolist(),
        })
        diag["status"] = "retained"
        diagnostics.append(diag)

    frame_record = {
        "cluster_id": ident, "filter": info["coadd_header_filter"],
        "redshift": float(frame["z"]), "cohort": str(frame.get("preparation_cohort", "unspecified")),
        "n_member_positions": len(member_rows), "n_member_positions_in_frame": match_stats["in_frame"],
        "n_unique_positional_matches_le_1arcsec": len(matched),
        "n_ambiguous_matches_rejected": match_stats["ambiguous"],
        "n_duplicate_detection_assignments_rejected": match_stats["duplicate_assignments"],
        "n_matches_50_500kpc": len(selected), "n_profiles_retained": len(rows),
        "n_targets_passing_coverage": n_target_coverage_pass,
        "n_targets_passing_min_controls": n_min_control_pass,
        "n_valid_control_profiles": n_controls_selected,
        "n_control_profiles_below_coverage": n_control_coverage_fail,
        "image_bytes": image_path.stat().st_size, "image_sha256": sha256(image_path),
        "mask_bytes": mask_path.stat().st_size, "mask_sha256": sha256(mask_path),
        "detection_table_sha256": sha256(det_path), "member_table_sha256": sha256(member_path),
    }
    return rows, frame_record, diagnostics


def summarize(rows, work, replicas, seed):
    bin_names = [f"{lo:g}-{hi:g}" for lo, hi in zip(R_EDGES_KPC[:-1], R_EDGES_KPC[1:])]
    cluster_rows = []
    by_cluster = {}
    for row in rows:
        by_cluster.setdefault((row["cluster_id"], row["filter"], row["cohort"]), []).append(row)
    for (cluster_id, band, cohort), group in sorted(by_cluster.items()):
        arrays = np.asarray([row["delta_sb"] for row in group], dtype=float)
        tg = np.asarray([row["target_sb"] for row in group], dtype=float)
        cg = np.asarray([row["control_sb"] for row in group], dtype=float)
        for j, name in enumerate(bin_names):
            cluster_rows.append({"cluster_id": cluster_id, "filter": band, "cohort": cohort,
                                 "radial_bin_kpc": name, "n_satellites": len(group),
                                 "mean_target_sb_njy_arcsec2": float(np.mean(tg[:, j])),
                                 "mean_control_sb_njy_arcsec2": float(np.mean(cg[:, j])),
                                 "mean_delta_sb_njy_arcsec2": float(np.mean(arrays[:, j]))})
    write_tsv(work / "cluster_profiles.tsv", cluster_rows,
              ["cluster_id", "filter", "cohort", "radial_bin_kpc", "n_satellites",
               "mean_target_sb_njy_arcsec2", "mean_control_sb_njy_arcsec2",
               "mean_delta_sb_njy_arcsec2"])

    filtered = {}
    filter_rows = []
    cohorts = sorted({row["cohort"] for row in cluster_rows})
    for cohort_index, cohort in enumerate(cohorts):
        filtered[cohort] = {}
        for band_index, band in enumerate(("i", "i2")):
            clusters = sorted({row["cluster_id"] for row in cluster_rows
                               if row["filter"] == band and row["cohort"] == cohort})
            if not clusters:
                continue
            per_cluster = {}
            for row in cluster_rows:
                if row["filter"] == band and row["cohort"] == cohort:
                    per_cluster.setdefault(row["cluster_id"], {})[row["radial_bin_kpc"]] = row["mean_delta_sb_njy_arcsec2"]
            matrix = np.array([[per_cluster[c][name] for name in bin_names] for c in clusters], dtype=float)
            rng = np.random.default_rng(seed + cohort_index * 100000 + band_index * 10000)
            boot = np.empty((replicas, len(bin_names)), dtype=float)
            for replica in range(replicas):
                draw = rng.integers(0, len(clusters), size=len(clusters))
                boot[replica] = np.mean(matrix[draw], axis=0)
            mean = np.mean(matrix, axis=0)
            lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
            n_satellites = sum(1 for row in rows if row["filter"] == band and row["cohort"] == cohort)
            filtered[cohort][band] = {"n_clusters": len(clusters), "n_satellites": n_satellites,
                                      "mean_delta": mean.tolist(), "bootstrap_low": lo.tolist(),
                                      "bootstrap_high": hi.tolist()}
            for j, name in enumerate(bin_names):
                filter_rows.append({"cohort": cohort, "filter": band, "radial_bin_kpc": name,
                                    "n_clusters": len(clusters), "n_satellites": n_satellites,
                                    "equal_cluster_mean_delta_sb_njy_arcsec2": float(mean[j]),
                                    "cluster_bootstrap_95_low": float(lo[j]),
                                    "cluster_bootstrap_95_high": float(hi[j])})
    write_tsv(work / "filter_profiles.tsv", filter_rows,
              ["cohort", "filter", "radial_bin_kpc", "n_clusters", "n_satellites",
               "equal_cluster_mean_delta_sb_njy_arcsec2", "cluster_bootstrap_95_low",
               "cluster_bootstrap_95_high"])
    return filtered


def plot_profiles(filter_summary, work):
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    edges = R_EDGES_KPC
    centers = np.sqrt(edges[:-1] * edges[1:])
    xerr = np.vstack([centers - edges[:-1], edges[1:] - centers])
    for ax, (cohort, band) in zip(axes.flat,
                                  [(c, b) for c in ("area_qualified", "mask_stress")
                                   for b in ("i", "i2")]):
        record = filter_summary.get(cohort, {}).get(band)
        ax.axhline(0., color="0.4", lw=1, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1.)
        ax.grid(True, which="both", alpha=.2)
        ax.set_title(f"{cohort}; {band}")
        if record:
            mean = np.asarray(record["mean_delta"], dtype=float)
            low = np.asarray(record["bootstrap_low"], dtype=float)
            high = np.asarray(record["bootstrap_high"], dtype=float)
            yerr = np.vstack([mean - low, high - mean])
            ax.errorbar(centers, mean, xerr=xerr, yerr=yerr, fmt="o-",
                        capsize=2, color="#245a87")
            ax.text(.03, .04,
                    f"clusters={record['n_clusters']}, satellites={record['n_satellites']}\n"
                    "95% cluster bootstrap; exploratory",
                    transform=ax.transAxes, fontsize=8, va="bottom")
        else:
            ax.text(.5, .5, "no retained profiles", ha="center", va="center",
                    transform=ax.transAxes)
        ax.set_xlabel("Satellite-centric radius (kpc)")
        ax.set_ylabel(r"Target − control (nJy arcsec$^{-2}$)")
    fig.suptitle("HSC satellite-centered outer-light pilot\n"
                 "Validation subset only; not a science measurement")
    fig.tight_layout()
    fig.savefig(work / "satellite_profile_review.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--limit", type=int, help="Engineering subset; omit to use all 50 validation frames")
    parser.add_argument("--coverage-floor", type=float, default=.60)
    parser.add_argument("--min-controls", type=int, default=3)
    parser.add_argument("--bootstrap-replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run CPU image processing on lageunha")
    if not 0 < args.coverage_floor <= 1 or args.min_controls < 1 or args.bootstrap_replicas < 100:
        parser.error("Invalid coverage/control/bootstrap configuration")
    data, work = args.data.resolve(), args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    manifest_path = data / "validation_sample_2026-10-13_01/validation_manifest.json"
    prep_dir = data / "detection_psf_2026-10-13_01"
    prep_summary_path = prep_dir / "summary.json"
    if not manifest_path.is_file() or not prep_summary_path.is_file():
        parser.error("Validation manifest or existing detection/mask summary is missing")
    manifest = read_json(manifest_path)
    prep_summary = read_json(prep_summary_path)
    if manifest.get("science_measurement") is not False or prep_summary.get("science_measurement") is not False:
        parser.error("Expected explicitly exploratory, non-science inputs")
    if prep_summary.get("status") != "complete" or prep_summary.get("n_prepared") != 50:
        parser.error("Expected the complete existing 50-frame preparation set")
    frames = manifest["clusters"][:args.limit] if args.limit else manifest["clusters"]
    if not frames:
        parser.error("No frames selected")
    started = time.time()
    cosmology = FlatLambdaCDM(H0=float(manifest.get("h0", 70.)),
                              Om0=float(manifest.get("omega_m", .3)), Tcmb0=0 * u.K)
    all_rows, frame_records, all_diagnostics, errors = [], [], [], []
    for ordinal, frame in enumerate(frames, 1):
        print(f"[{ordinal}/{len(frames)}] measuring cluster {frame['id']}", flush=True)
        try:
            rows, record, diagnostics = process_frame(frame, prep_dir, cosmology,
                                                       args.coverage_floor, args.min_controls)
            all_rows.extend(rows); frame_records.append(record); all_diagnostics.extend(diagnostics)
            print(f"  {record['filter']}: matched={record['n_unique_positional_matches_le_1arcsec']}, "
                  f"50-500 kpc={record['n_matches_50_500kpc']}, "
                  f"target-coverage-pass={record['n_targets_passing_coverage']}, "
                  f"profiles={len(rows)}", flush=True)
        except Exception as exc:
            errors.append({"cluster_id": str(frame.get("id")), "error_type": type(exc).__name__,
                           "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
    if not frame_records:
        parser.error("No frame completed; no output directory created")
    work.mkdir(parents=True, exist_ok=False)
    radial_names = [f"{lo:g}-{hi:g}" for lo, hi in zip(R_EDGES_KPC[:-1], R_EDGES_KPC[1:])]
    profile_rows = []
    for row in all_rows:
        for j, name in enumerate(radial_names):
            profile_rows.append({
                "cluster_id": row["cluster_id"], "filter": row["filter"],
                "redshift": row["redshift"], "richness": row["richness"],
                "cohort": row["cohort"], "member_index": row["member_index"],
                "detection_number": row["detection_number"],
                "member_detection_separation_arcsec": row["member_detection_separation_arcsec"],
                "cluster_radius_kpc": row["cluster_radius_kpc"], "n_controls": row["n_controls"],
                "radial_bin_kpc": name,
                "target_sb_njy_arcsec2": row["target_sb"][j],
                "control_sb_njy_arcsec2": row["control_sb"][j],
                "delta_sb_njy_arcsec2": row["delta_sb"][j],
                "target_err_nominal_njy_arcsec2": row["target_err_nominal"][j],
                "target_coverage": row["target_coverage"][j],
                "control_coverage": row["control_coverage"][j],
                "target_n_pixels": row["target_n_pixels"][j],
            })
    write_tsv(work / "satellite_profiles.tsv", profile_rows,
              ["cluster_id", "filter", "redshift", "richness", "cohort", "member_index",
               "detection_number", "member_detection_separation_arcsec", "cluster_radius_kpc",
               "n_controls", "radial_bin_kpc", "target_sb_njy_arcsec2", "control_sb_njy_arcsec2",
               "delta_sb_njy_arcsec2", "target_err_nominal_njy_arcsec2", "target_coverage",
               "control_coverage", "target_n_pixels"])
    write_tsv(work / "candidate_diagnostics.tsv", all_diagnostics,
              ["cluster_id", "filter", "cohort", "member_index", "detection_number", "cluster_radius_kpc",
               "member_detection_separation_arcsec", "target_min_coverage", "target_coverage_by_bin",
               "candidate_control_positions", "valid_control_profiles", "valid_controls_by_bin", "status"])
    filter_summary = summarize(all_rows, work, args.bootstrap_replicas, args.seed) if all_rows else {}
    plot_profiles(filter_summary, work)
    script_path = Path(__file__).resolve()
    summary = {
        "schema": "hsc-satellite-outer-profile-pilot-v1",
        "status": "complete" if not errors else "partial",
        "science_measurement": False,
        "n_frames_requested": len(frames), "n_frames_completed": len(frame_records),
        "n_frame_errors": len(errors), "n_satellite_profiles": len(all_rows),
        "filters": filter_summary,
        "method": {
            "member_selection": "CAMIRA positions matched to a unique saved SEP detection within 1 arcsec; no mass or CAMIRA weight applied",
            "cluster_radius_kpc": [MIN_CLUSTER_RADIUS_KPC, MAX_CLUSTER_RADIUS_KPC],
            "satellite_radius_bins_kpc": R_EDGES_KPC.tolist(),
            "target_masks": "preserve matched SEP label and Kron ellipse; mask all other detected sources, fallback unmatched CAMIRA positions, and original HSC quality/bright-object pixels",
            "control_positions": f"same cluster-centric projected radius; up to {N_CONTROLS_MAX} of {N_CONTROL_ANGLES} azimuthal rotations; reject positions within {CONTROL_SOURCE_CLEARANCE_ARCSEC} arcsec of any detection or {CONTROL_MEMBER_CLEARANCE_KPC} kpc of any CAMIRA member",
            "coverage_floor": args.coverage_floor, "minimum_controls": args.min_controls,
            "stack_weighting": "average satellite profiles within each cluster, then equal-weight clusters; filters kept separate",
            "flux": "FLUXMAG0 converted to nJy/arcsec^2; no extinction or passband correction",
            "uncertainty": "cluster bootstrap interval across the selected validation clusters; nominal pixel variance excludes coadd correlated noise",
        },
        "limitations": [
            "The 50 validation frames are an area-selected engineering subset, not a representative or science-complete sample.",
            "CAMIRA member identities are positional only; exact member-file column semantics, membership weights, stellar mass, color and photometric redshifts are not used.",
            "The same-radius azimuthal controls estimate a local diffuse-light/sky baseline but do not remove asymmetric substructure or prove a satellite/ICL decomposition.",
            "Source masks, HSC coadd correlated noise, PSF wings, sky transfer, extinction, passband response, completeness and member contamination remain unvalidated.",
            "The exploratory profiles cannot be reported as a satellite outer-light detection or a physical correlation with ICL.",
        ],
        "manifest_sha256": sha256(manifest_path),
        "preparation_summary_sha256": sha256(prep_summary_path),
        "script_sha256": sha256(script_path),
        "frame_records": frame_records, "errors": errors,
        "elapsed_seconds": time.time() - started,
        "bootstrap": {"replicas": args.bootstrap_replicas, "seed": args.seed},
        "output_sha256": {},
    }
    summary_path = work / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                            encoding="utf-8")
    for path in sorted(work.iterdir()):
        if path.is_file() and path.name != "summary.json":
            summary["output_sha256"][str(path)] = sha256(path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "frames": len(frame_records),
                      "profiles": len(all_rows), "filters": filter_summary,
                      "output": str(work)}, ensure_ascii=False), flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
