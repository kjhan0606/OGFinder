#!/usr/bin/env python3
"""Prepare automatic source-mask candidates and radial area for the 800 clusters.

This creates masks and support diagnostics only. It does not measure or stack
surface brightness, and automatic BCG matches are not science approvals.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

import prepare_detection_psf as prep
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, validate_fits
from prestack_core import EDGES, FACTORS, geometry


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def write_tsv(path, rows, fields):
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def prepare_one(frame, image_path, file_record, output_dir, cosmology):
    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] != frame.get("coadd_header_filter"):
        raise ValueError("stored primary FILTER differs from the download manifest")
    with fits.open(image_path, memmap=False) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
        header = hdul[info["science_hdu"]].header.copy()
        header.extend(hdul[0].header, unique=True)

    wcs = WCS(header).celestial
    scale = info["pixel_scale_arcsec"]
    center = np.asarray(info["center_pixel_xy"], dtype=float)
    member_sky = np.asarray(frame["members_radec_deg"], dtype=float).reshape(-1, 2)
    members = (np.column_stack(wcs.all_world2pix(member_sky[:, 0], member_sky[:, 1], 0))
               if len(member_sky) else np.empty((0, 2)))
    inside = (np.isfinite(members).all(axis=1) & (members[:, 0] >= 0) &
              (members[:, 0] < image.shape[1]) & (members[:, 1] >= 0) &
              (members[:, 1] < image.shape[0]))
    members_in = members[inside]
    bad = (bits & sum(1 << bit for bit in BAD_BIT_POSITIONS)) != 0
    bad |= ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    bright = (bits & BRIGHT_OBJECT_BIT) != 0
    base = bad | bright

    objects, segmentation, kron, kron_flags, bcg, _, detection = prep.detect(
        image, variance, bad, center, members_in, scale)
    masks, fallback = prep.source_masks(
        objects, segmentation, kron, kron_flags, base, bcg, members_in, scale)
    if not all(np.all(mask[base]) for mask in masks):
        raise AssertionError("a proposed source mask removed a quality exclusion")
    if not (np.all(~masks[0] | masks[1]) and np.all(~masks[1] | masks[2])):
        raise AssertionError("Kron masks are not nested")

    kpc_per_pixel = (cosmology.kpc_proper_per_arcmin(frame["z"]).value * scale / 60.0)
    geom = geometry(image.shape, center, kpc_per_pixel)
    radial_bins = geom[3]
    annulus_area = np.pi * np.diff(EDGES ** 2) / kpc_per_pixel ** 2
    ix, iy = int(round(center[0])), int(round(center[1]))
    center_in = 0 <= ix < image.shape[1] and 0 <= iy < image.shape[0]
    output_dir.mkdir(parents=True, exist_ok=False)
    wh = wcs.to_header()
    wh["SCIVAL"] = False
    wh["AUTOMASK"] = True
    wh["BCGCAND"] = bool(bcg >= 0)
    hdus = [fits.PrimaryHDU(), fits.ImageHDU(base.astype(np.uint8), header=wh, name="BASE_MASK")]
    for factor, mask in zip(FACTORS, masks):
        hdus.append(fits.ImageHDU(mask.astype(np.uint8), header=wh,
                                  name="MASK_K" + str(factor).replace(".", "P")))
    mask_path = output_dir / "proposed_masks.fits.gz"
    fits.HDUList(hdus).writeto(mask_path, overwrite=False, checksum=True)

    coverage_rows = []
    coverages = []
    counts = []
    for factor, mask in zip(FACTORS, masks):
        use = (~mask) & (radial_bins >= 0)
        count = np.bincount(radial_bins[use], minlength=len(EDGES) - 1)
        coverage = np.minimum(1.0, count / annulus_area)
        counts.append(count)
        coverages.append(coverage)
        for ring, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
            coverage_rows.append(dict(id=str(frame["id"]), kron_factor=factor,
                                      r_min_kpc=lo, r_max_kpc=hi, n_valid=int(count[ring]),
                                      coverage=float(coverage[ring])))

    center_masked = ([bool(mask[iy, ix]) for mask in masks] if center_in else [True] * len(masks))
    record = dict(
        id=str(frame["id"]), image=str(image_path), image_bytes=int(file_record["bytes"]),
        image_sha256=file_record["sha256"], primary_filter=info["coadd_header_filter"],
        n_detections=int(len(objects)), bcg_match_count=int(detection["bcg_matches"]),
        auto_unique_center_candidate=bool(bcg >= 0),
        nearest_detection_arcsec=detection["nearest_bcg_detection_arcsec"],
        center_base_masked=bool(base[iy, ix]) if center_in else True,
        center_masked_by_kron=dict(zip((str(f) for f in FACTORS), center_masked)),
        n_members=int(len(member_sky)), n_members_in_frame=int(inside.sum()),
        member_fallback_count=int(fallback),
        min_coverage_by_kron={str(f): float(np.min(c)) for f, c in zip(FACTORS, coverages)},
        coverage_by_kron={str(f): [float(x) for x in c] for f, c in zip(FACTORS, coverages)},
        counts_by_kron={str(f): [int(x) for x in c] for f, c in zip(FACTORS, counts)},
        mask_file=str(mask_path), mask_sha256=sha256(mask_path),
        science_measurement=False, mask_ready_for_science=False,
    )
    return record, coverage_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output directory: {args.work}")

    manifest_path = args.sample / "download_01/download_manifest.json"
    status_path = args.sample / "download_01/download_status.json"
    manifest = json.loads(manifest_path.read_text())
    download_status = json.loads(status_path.read_text())
    frames = manifest.get("clusters", [])
    completed = {str(row["id"]): row for row in download_status.get("completed", [])}
    frame_ids = [str(frame["id"]) for frame in frames]
    if download_status.get("status") != "complete" or len(frames) != 800:
        raise ValueError("Expected the complete 800-image download manifest")
    if len(set(frame_ids)) != 800 or set(frame_ids) != set(completed):
        raise ValueError("Download manifest and completed ID sets differ")

    cosmology = FlatLambdaCDM(H0=float(manifest["h0"]), Om0=float(manifest["omega_m"]), Tcmb0=0)
    args.work.mkdir(parents=True, exist_ok=False)
    (args.work / "clusters").mkdir()
    write_json(args.work / "protocol.json", dict(
        method="prepare_detection_psf.detect + source_masks; native HSC masks retained",
        kron_factors=list(FACTORS), radial_edges_kpc=EDGES.tolist(),
        source_manifest=str(manifest_path), source_status=str(status_path),
        n_expected=800, no_surface_brightness_measurement=True,
        bcg_status="automatic positional candidate only; no new manual association review",
        area_status="engineering support only; no threshold or mask approved for science",
        science_measurement=False, mask_ready_for_science=False,
    ))
    summary_path = args.work / "summary.json"
    records_path = args.work / "records.tsv"
    coverage_path = args.work / "coverage.tsv"
    record_fields = ["id", "image", "image_bytes", "image_sha256", "primary_filter",
                     "n_detections", "bcg_match_count", "auto_unique_center_candidate",
                     "nearest_detection_arcsec", "center_base_masked", "center_masked_by_kron",
                     "n_members", "n_members_in_frame", "member_fallback_count",
                     "min_coverage_by_kron", "coverage_by_kron", "counts_by_kron",
                     "mask_file", "mask_sha256", "science_measurement", "mask_ready_for_science"]
    coverage_fields = ["id", "kron_factor", "r_min_kpc", "r_max_kpc", "n_valid", "coverage"]
    records = []
    coverage_rows = []
    errors = []
    start = time.monotonic()
    for index, frame in enumerate(frames, 1):
        ident = str(frame["id"])
        file_record = completed[ident]
        image_path = args.sample / "download_01/images" / file_record["file"]
        try:
            if not image_path.is_file() or image_path.stat().st_size != int(file_record["bytes"]):
                raise ValueError("downloaded FITS is missing or its byte size changed")
            record, field_coverage = prepare_one(
                frame, image_path, file_record, args.work / "clusters" / ("cluster_" + ident), cosmology)
            records.append(record)
            coverage_rows.extend(field_coverage)
            print(f"[{index}/800] cluster {ident}: detections={record['n_detections']} "
                  f"unique_center={record['auto_unique_center_candidate']}", flush=True)
        except Exception as exc:
            errors.append(dict(id=ident, error=f"{type(exc).__name__}: {exc}"))
            print(f"[{index}/800] cluster {ident}: ERROR {type(exc).__name__}: {exc}", flush=True)
        if index % 10 == 0 or index == len(frames):
            write_tsv(records_path, records, record_fields) if records else None
            write_tsv(coverage_path, coverage_rows, coverage_fields) if coverage_rows else None
            write_json(args.work / "progress.json", dict(
                completed_index=index, n_success=len(records), n_error=len(errors),
                elapsed_seconds=time.monotonic() - start,
            ))

    hashes = {str(manifest_path): sha256(manifest_path), str(status_path): sha256(status_path),
              str(args.work / "protocol.json"): sha256(args.work / "protocol.json")}
    summary = dict(
        status="complete" if len(records) == 800 and not errors else "partial",
        n_expected=800, n_success=len(records), n_error=len(errors), errors=errors,
        n_auto_unique_center_candidates=sum(row["auto_unique_center_candidate"] for row in records),
        n_center_base_masked=sum(row["center_base_masked"] for row in records),
        n_center_masked_by_kron={str(f): sum(row["center_masked_by_kron"][str(f)] for row in records)
                                 for f in FACTORS},
        n_by_filter={f: sum(row["primary_filter"] == f for row in records) for f in ("i", "i2")},
        elapsed_seconds=time.monotonic() - start,
        source_sha256=hashes,
        output_sha256={str(path): sha256(path) for path in (records_path, coverage_path)
                       if path.exists()},
        mask_output_sha256={row["mask_file"]: row["mask_sha256"] for row in records},
        script_sha256={str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
                       str(Path(prep.__file__).resolve()): sha256(Path(prep.__file__).resolve())},
        science_measurement=False, mask_ready_for_science=False,
    )
    write_json(summary_path, summary)
    print(f"Mask support preparation {summary['status']}: {len(records)}/800; "
          f"{len(errors)} errors; elapsed={summary['elapsed_seconds']:.1f}s", flush=True)
    if errors:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
