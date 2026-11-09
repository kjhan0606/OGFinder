#!/usr/bin/env python3
"""Prepare automatic source-mask candidates and radial area for HSC cutouts.

This creates masks and support diagnostics only. It does not measure or stack
surface brightness, and automatic BCG matches are not science approvals. The
legacy ``--sample`` input remains fixed to the original 800-target sample;
``--download-dir`` accepts a completed arbitrary-size download when paired with
its pinned approved manifest and completed header/WCS audit.
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
from download_pilot import BAD_MASK_PLANE_NAMES, mask_plane_bitmask, validate_fits
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
    bad = (bits & mask_plane_bitmask(info["mask_planes"], BAD_MASK_PLANE_NAMES)) != 0
    bad |= ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    bright = (bits & mask_plane_bitmask(info["mask_planes"], ("MP_BRIGHT_OBJECT",))) != 0
    base = bad | bright

    objects, segmentation, kron, kron_flags, bcg, _, detection = prep.detect(
        image, variance, bad, center, members_in, scale)
    role = str(frame.get("frame_role", frame.get("sample_role", "target")))
    is_control = role in {"matched_control", "control"}
    # A random control center is not a BCG. Do not protect the central source
    # from the control mask even if SEP finds a compact object there.
    masks, fallback = prep.source_masks(
        objects, segmentation, kron, kron_flags, base,
        -1 if is_control else bcg, members_in, scale)
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
    wh["BCGCAND"] = bool(bcg >= 0) and not is_control
    wh["CTRLPT"] = is_control
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
        id=str(frame["id"]), frame_role=role, image=str(image_path), image_bytes=int(file_record["bytes"]),
        image_sha256=file_record["sha256"], primary_filter=info["coadd_header_filter"],
        n_detections=int(len(objects)), center_source_match_count=int(detection["bcg_matches"]),
        bcg_match_count=0 if is_control else int(detection["bcg_matches"]),
        auto_unique_center_candidate=bool(bcg >= 0) and not is_control,
        center_source_protected=bool(bcg >= 0) and not is_control,
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
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--sample", type=Path,
                        help="legacy 800-target sample directory containing download_01/")
    inputs.add_argument("--download-dir", type=Path,
                        help="completed download work directory containing status, manifest, and images/")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--expected-n", type=int,
                        help="required frame count with --download-dir; defaults to 800 for --sample")
    parser.add_argument("--approved-manifest", type=Path,
                        help="pinned target/control manifest used to authorize this download")
    parser.add_argument("--expected-manifest-sha256",
                        help="expected SHA-256 of --approved-manifest")
    parser.add_argument("--audit-summary", type=Path,
                        help="completed independent FITS/header/WCS audit summary")
    args = parser.parse_args()
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output directory: {args.work}")

    if args.download_dir is not None and (args.expected_n is None or args.expected_n <= 0):
        parser.error("--download-dir requires a positive --expected-n")
    if args.sample is not None and args.expected_n not in (None, 800):
        parser.error("--sample is the fixed legacy 800-target input")
    if (args.approved_manifest is None) != (args.expected_manifest_sha256 is None):
        parser.error("--approved-manifest and --expected-manifest-sha256 must be supplied together")
    if args.download_dir is not None and (args.approved_manifest is None or args.audit_summary is None):
        parser.error("--download-dir requires its pinned approved manifest and --audit-summary")
    download_dir = (args.download_dir.resolve() if args.download_dir is not None
                    else (args.sample / "download_01").resolve())
    expected_n = args.expected_n if args.expected_n is not None else 800
    manifest_path = download_dir / "download_manifest.json"
    status_path = download_dir / "download_status.json"
    manifest = json.loads(manifest_path.read_text())
    download_status = json.loads(status_path.read_text())
    frames = manifest.get("clusters", [])
    completed = {str(row["id"]): row for row in download_status.get("completed", [])}
    frame_ids = [str(frame["id"]) for frame in frames]
    if (download_status.get("status") != "complete"
            or download_status.get("downloaded", len(download_status.get("completed", []))) != expected_n
            or download_status.get("science_measurement") is not False
            or download_status.get("credentials_saved") is not False
            or len(frames) != expected_n):
        raise ValueError(f"Expected the complete {expected_n}-image download manifest")
    if (len(set(frame_ids)) != expected_n
            or len(download_status.get("completed", [])) != expected_n
            or set(frame_ids) != set(completed)):
        raise ValueError("Download manifest and completed ID sets differ")
    approved_manifest_sha = None
    if args.approved_manifest is not None:
        approved_manifest = args.approved_manifest.resolve()
        approved_manifest_sha = sha256(approved_manifest)
        if approved_manifest_sha != args.expected_manifest_sha256:
            raise ValueError("approved manifest SHA-256 differs from the pinned value")
        if download_status.get("manifest_sha256") != approved_manifest_sha:
            raise ValueError("download checkpoint is not tied to the approved manifest")
        approved_doc = json.loads(approved_manifest.read_text(encoding="utf-8"))
        approved_rows = {str(row["id"]): row for row in approved_doc.get("clusters", [])}
        if len(approved_rows) != expected_n or set(approved_rows) != set(frame_ids):
            raise ValueError("downloaded frame identities differ from the approved manifest")
        for frame in frames:
            source = approved_rows[str(frame["id"])]
            for key in ("ra_deg", "dec_deg", "z"):
                if float(frame[key]) != float(source[key]):
                    raise ValueError(f"downloaded frame {frame['id']} changed approved {key}")
            approved_role = str(source.get("frame_role", source.get("sample_role", "")))
            downloaded_role = str(frame.get("frame_role", frame.get("sample_role", "")))
            if approved_role != downloaded_role:
                raise ValueError(f"downloaded frame {frame['id']} changed its approved role")
    audit_summary = None
    if args.audit_summary is not None:
        audit_path = args.audit_summary.resolve()
        audit_summary = json.loads(audit_path.read_text(encoding="utf-8"))
        audit_manifest_sha = audit_summary.get("manifest_sha256")
        audit_inputs = audit_summary.get("inputs_sha256", {})
        if audit_manifest_sha is None and approved_manifest_sha is not None:
            audit_manifest_sha = audit_inputs.get(str(args.approved_manifest.resolve()))
        if (audit_summary.get("status") != "complete_header_wcs_file_integrity_audit"
                or audit_summary.get("n_expected_frames") != expected_n
                or audit_summary.get("n_verified_frames") != expected_n
                or audit_summary.get("science_measurement") is not False
                or audit_summary.get("images_ready_for_stacking") is not False
                or audit_inputs.get(str(status_path)) != sha256(status_path)
                or audit_inputs.get(str(manifest_path)) != sha256(manifest_path)
                or (approved_manifest_sha is not None
                    and audit_manifest_sha != approved_manifest_sha)):
            raise ValueError("independent image audit is incomplete or belongs to different inputs")

    cosmology = FlatLambdaCDM(H0=float(manifest["h0"]), Om0=float(manifest["omega_m"]), Tcmb0=0)
    args.work.mkdir(parents=True, exist_ok=False)
    (args.work / "clusters").mkdir()
    write_json(args.work / "protocol.json", dict(
        method="prepare_detection_psf.detect + source_masks; native HSC masks retained",
        kron_factors=list(FACTORS), radial_edges_kpc=EDGES.tolist(),
        source_manifest=str(manifest_path), source_status=str(status_path),
        n_expected=expected_n, no_surface_brightness_measurement=True,
        bcg_status="automatic positional candidate only; no new manual association review",
        area_status="engineering support only; no threshold or mask approved for science",
        approved_manifest_sha256=approved_manifest_sha,
        audit_summary=str(args.audit_summary.resolve()) if args.audit_summary else None,
        science_measurement=False, mask_ready_for_science=False,
    ))
    summary_path = args.work / "summary.json"
    records_path = args.work / "records.tsv"
    coverage_path = args.work / "coverage.tsv"
    record_fields = ["id", "frame_role", "image", "image_bytes", "image_sha256", "primary_filter",
                     "n_detections", "center_source_match_count", "bcg_match_count",
                     "auto_unique_center_candidate", "center_source_protected",
                     "nearest_detection_arcsec", "center_base_masked", "center_masked_by_kron",
                     "n_members", "n_members_in_frame", "member_fallback_count",
                     "min_coverage_by_kron", "coverage_by_kron", "counts_by_kron",
                     "mask_file", "mask_sha256", "science_measurement", "mask_ready_for_science"]
    coverage_fields = ["id", "kron_factor", "r_min_kpc", "r_max_kpc", "n_valid", "coverage"]
    records = []
    coverage_rows = []
    errors = []
    image_dir = download_dir / "images"
    start = time.monotonic()
    for index, frame in enumerate(frames, 1):
        ident = str(frame["id"])
        file_record = completed[ident]
        image_path = image_dir / file_record["file"]
        try:
            if not image_path.is_file() or image_path.stat().st_size != int(file_record["bytes"]):
                raise ValueError("downloaded FITS is missing or its byte size changed")
            record, field_coverage = prepare_one(
                frame, image_path, file_record, args.work / "clusters" / ("cluster_" + ident), cosmology)
            records.append(record)
            coverage_rows.extend(field_coverage)
            print(f"[{index}/{expected_n}] cluster {ident}: detections={record['n_detections']} "
                  f"unique_center={record['auto_unique_center_candidate']}", flush=True)
        except Exception as exc:
            errors.append(dict(id=ident, error=f"{type(exc).__name__}: {exc}"))
            print(f"[{index}/{expected_n}] cluster {ident}: ERROR {type(exc).__name__}: {exc}", flush=True)
        if index % 10 == 0 or index == len(frames):
            write_tsv(records_path, records, record_fields) if records else None
            write_tsv(coverage_path, coverage_rows, coverage_fields) if coverage_rows else None
            write_json(args.work / "progress.json", dict(
                completed_index=index, n_success=len(records), n_error=len(errors),
                elapsed_seconds=time.monotonic() - start,
            ))

    hashes = {str(manifest_path): sha256(manifest_path), str(status_path): sha256(status_path),
              str(args.work / "protocol.json"): sha256(args.work / "protocol.json")}
    if args.approved_manifest is not None:
        hashes[str(args.approved_manifest.resolve())] = approved_manifest_sha
    if args.audit_summary is not None:
        hashes[str(args.audit_summary.resolve())] = sha256(args.audit_summary.resolve())
    summary = dict(
        status="complete" if len(records) == expected_n and not errors else "partial",
        n_expected=expected_n, n_success=len(records), n_error=len(errors), errors=errors,
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
    print(f"Mask support preparation {summary['status']}: {len(records)}/{expected_n}; "
          f"{len(errors)} errors; elapsed={summary['elapsed_seconds']:.1f}s", flush=True)
    if errors:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
