#!/usr/bin/env python3
"""Independently audit the PSF-convolved joint ICL/satellite injection batch."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def coord_key(ra: float, dec: float) -> tuple[float, float]:
    return round(float(ra), 9), round(float(dec), 9)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--test-center-psf", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path,
                        help="same support-selected redetection run without PSF convolution")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    plan = json.loads((args.plan_dir / "psf_picker_plan.json").read_text())
    mapping_path = args.plan_dir / "position_mapping.json"
    audit_path = args.download_dir / "psf_picker_audit.json"
    mapping = json.loads(mapping_path.read_text())
    psf_audit = json.loads(audit_path.read_text())
    summary_path = args.output_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    if plan.get("user_approved_requests") != 744 or plan.get("filter") != "HSC-I":
        raise ValueError("unexpected approved PSF plan scope")
    if sha256(mapping_path) != plan.get("position_mapping_sha256"):
        raise ValueError("approved mapping checksum mismatch")
    if psf_audit.get("status") != "verified" or psf_audit.get("n_returned") != 744:
        raise ValueError("PSF download is not the audited 744-file batch")
    if summary.get("science_measurement") is not False:
        raise ValueError("injection batch incorrectly claims a science measurement")
    if not summary.get("injection", {}).get("psf_convolution_applied"):
        raise ValueError("batch summary does not confirm PSF convolution")
    code_hashes = summary.get("code_sha256", {})
    if not code_hashes or any(sha256(Path(path)) != digest
                              for path, digest in code_hashes.items()):
        raise ValueError("run-time code hashes do not match the current repository files")

    expected_frames: set[tuple[str, str]] = set()
    expected_indices: dict[tuple[str, str], set[int]] = {}
    role_coordinates: dict[tuple[str, str, str, int | None], tuple[float, float]] = {}
    for position in mapping:
        coord = coord_key(position["ra_deg"], position["dec_deg"])
        for role in position["roles"]:
            cluster, filt = str(role["cluster_id"]), str(role["frame_filter"]).lower()
            role_name = str(role["role"])
            index = None if role.get("target_index") is None else int(role["target_index"])
            role_key = (cluster, filt, role_name, index)
            if role_key in role_coordinates:
                raise ValueError(f"duplicate planned PSF role {role_key}")
            role_coordinates[role_key] = coord
            if role_name == "cluster_center":
                expected_frames.add((cluster, filt))
            elif role_name == "satellite_injection_center":
                expected_frames.add((cluster, filt))
                expected_indices.setdefault((cluster, filt), set()).add(index)
            else:
                raise ValueError(f"unknown PSF plan role: {role_name}")
    test_frame = ("10", "i2")
    expected_frames.add(test_frame)
    test_key = ("10", "i2", "cluster_center", None)
    if test_key in role_coordinates:
        raise ValueError("cluster-10 test center should be supplied separately")
    test_hash = "1276d46113d389116ded6f63bae4c9865bb42d24244b22247a84d2c6929cb202"
    if sha256(args.test_center_psf) != test_hash:
        raise ValueError("test-center PSF hash differs from the approved test record")

    approved_files = {}
    for row in psf_audit["files"]:
        key = coord_key(row["ra_deg_filename"], row["dec_deg_filename"])
        if key in approved_files:
            raise ValueError("duplicate coordinate in PSF audit manifest")
        approved_files[key] = (str((args.download_dir / row["file"]).resolve()), row["sha256"])
    if len(approved_files) != 744:
        raise ValueError("expected 744 distinct returned PSF positions")

    expected_pairs = sorted(expected_frames, key=lambda pair: (int(pair[0]), pair[1]))
    listed_pairs = {(str(row["cluster_id"]), str(row["filter"]).lower())
                    for row in summary["selection"]["selected_pairs"]}
    if listed_pairs != expected_frames:
        raise ValueError("batch selected frames do not exactly match the PSF-approved frame set")

    reference = None
    if args.reference_dir:
        reference = json.loads((args.reference_dir / "summary.json").read_text())
        if reference.get("science_measurement") is not False:
            raise ValueError("reference run must remain exploratory")
        reference_by_pair = {(str(row["cluster_id"]), str(row["filter"]).lower()): row
                             for row in reference.get("records", [])}
    else:
        reference_by_pair = {}

    checked_records = []
    checked_pairs = set()
    response_deltas = {key: [] for key in (
        "icl_response_to_diffuse_unit_injection",
        "icl_plane_corrected_response_to_diffuse_unit_injection",
        "icl_response_to_all_satellite_envelopes",
        "icl_plane_corrected_response_to_all_satellite_envelopes",
        "satellite_response_to_smooth_diffuse",
        "satellite_response_to_satellite_unit_injection",
    )}
    response_deltas_by_filter = {
        filt: {key: [] for key in response_deltas} for filt in ("i", "i2")
    }
    output_root = args.output_dir.resolve()
    for row in summary.get("records", []):
        pair = (str(row["cluster_id"]), str(row["filter"]).lower())
        if pair in checked_pairs:
            raise ValueError(f"duplicate frame result {pair}")
        checked_pairs.add(pair)
        result_path = Path(row["path"]).resolve()
        if result_path.parent != output_root or sha256(result_path) != row["sha256"]:
            raise ValueError(f"frame result path/hash mismatch for {pair}")
        result = json.loads(result_path.read_text())
        frame_code_hashes = result.get("code_sha256", {})
        if (not frame_code_hashes
                or any(code_hashes.get(path) != digest
                       for path, digest in frame_code_hashes.items())):
            raise ValueError(f"frame {pair} used a different core code snapshot than the batch")
        if result.get("science_measurement") is not False:
            raise ValueError(f"frame {pair} incorrectly claims a science measurement")
        injection = result["injection"]
        n_sat = len(expected_indices.get(pair, set()))
        if (injection.get("psf_convolution_applied") is not True
                or injection.get("source_detection_and_masks_rerun") is not True
                or injection.get("unresolved_bcg_policy") != "mask_all"
                or injection.get("variance_and_poisson_noise_changed") is not False
                or injection.get("psf_file_count") != n_sat + 1):
            raise ValueError(f"PSF/redetection/noise flags mismatch for {pair}")
        indices = expected_indices.get(pair, set())
        if indices != set(range(n_sat)):
            raise ValueError(f"satellite PSF target indices are not contiguous for {pair}")

        input_hashes = result.get("input_sha256", {})
        used_psf_paths = {Path(path).resolve(): digest for path, digest in input_hashes.items()
                          if ((Path(path).suffix.lower() == ".fits"
                               and Path(path).name.startswith("psf_"))
                              or Path(path).name == args.test_center_psf.name)}
        expected_paths = {}
        if pair == test_frame:
            expected_paths[str(args.test_center_psf.resolve())] = test_hash
        center_role = (pair[0], pair[1], "cluster_center", None)
        if pair != test_frame and center_role not in role_coordinates:
            raise ValueError(f"missing approved center PSF role for {pair}")
        roles = [center_role] + [(pair[0], pair[1], "satellite_injection_center", i)
                                 for i in sorted(indices)]
        for role in roles:
            if pair == test_frame and role == center_role:
                continue
            coord = role_coordinates[role]
            path, digest = approved_files[coord]
            expected_paths[path] = digest
        if used_psf_paths != {Path(path): digest for path, digest in expected_paths.items()}:
            raise ValueError(f"used PSF files differ from approved position mapping for {pair}")

        deltas = {}
        ref = reference_by_pair.get(pair)
        if ref:
            for key in response_deltas:
                value = float(row["transfer_response"][key]) - float(ref["transfer_response"][key])
                if not np.isfinite(value):
                    raise ValueError(f"non-finite paired PSF response difference for {pair}/{key}")
                response_deltas[key].append(value)
                if pair[1] in response_deltas_by_filter:
                    response_deltas_by_filter[pair[1]][key].append(value)
                deltas[key] = value
        checked_records.append({"cluster_id": pair[0], "filter": pair[1],
                                "n_satellite_psfs": n_sat,
                                "n_psf_files_verified": len(used_psf_paths),
                                "paired_no_psf_response_deltas": deltas})

    error_pairs = {(str(row["cluster_id"]), str(row["filter"]).lower())
                   for row in summary.get("errors", [])}
    if (checked_pairs & error_pairs or checked_pairs | error_pairs != expected_frames
            or len(error_pairs) != int(summary.get("n_errors", -1))):
        raise ValueError("frame results and errors do not partition the approved PSF sample")
    paired_summary = {}
    def summarize_deltas(values):
        return ({
            "n_paired_frames": len(values),
            "median_psf_minus_no_psf": float(np.median(values)),
            "q16": float(np.quantile(values, 0.16)),
            "q84": float(np.quantile(values, 0.84)),
        } if values else {"n_paired_frames": 0, "median_psf_minus_no_psf": None,
                          "q16": None, "q84": None})

    for key, values in response_deltas.items():
        paired_summary[key] = summarize_deltas(values)
    paired_by_filter = {
        filt: {key: summarize_deltas(values) for key, values in metrics.items()}
        for filt, metrics in response_deltas_by_filter.items()
    }

    result = {
        "schema": "astrafek-joint-icl-satellite-psf-redetection-audit-v1",
        "status": "verified" if not summary.get("n_errors") and len(checked_records) == len(expected_frames)
                  else "verified_partial",
        "science_measurement": False,
        "n_approved_frames": len(expected_frames),
        "n_frame_results_checked": len(checked_records),
        "n_errors_reported_by_batch": summary.get("n_errors"),
        "approved_frame_pairs": [{"cluster_id": c, "filter": f} for c, f in expected_pairs],
        "batch_summary_sha256": sha256(summary_path),
        "psf_plan_sha256": sha256(args.plan_dir / "psf_picker_plan.json"),
        "psf_mapping_sha256": sha256(mapping_path),
        "psf_download_audit_sha256": sha256(audit_path),
        "audit_script_sha256": sha256(Path(__file__)),
        "code_sha256": code_hashes,
        "records": checked_records,
        "paired_no_psf_comparison": paired_summary,
        "paired_no_psf_comparison_by_filter": paired_by_filter,
        "limitations": [
            "The test uses compact model-PSF cores only; extended wings are not included.",
            "Injected source Poisson noise, HSC upstream sky transfer, and correlated coadd noise are not modeled.",
            "Results are support-selected injection/mask sensitivities, not observational science measurements.",
        ],
    }
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in (
        "status", "n_approved_frames", "n_frame_results_checked",
        "n_errors_reported_by_batch", "paired_no_psf_comparison",
        "science_measurement",
    )}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
