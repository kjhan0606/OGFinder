#!/usr/bin/env python3
"""Audit an ICL-supported redetection batch against its approved PSF bank."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from audit_joint_icl_satellite_psf_redetection import validate_support_accounting
from position_psf import sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", type=Path, required=True)
    parser.add_argument("--bank-audit", type=Path, required=True)
    parser.add_argument("--approved-plan-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    batch_dir = args.batch_dir.resolve()
    summary_path = batch_dir / "summary.json"
    bank_path = args.bank_audit.resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    bank = json.loads(bank_path.read_text(encoding="utf-8"))
    approved_hash = str(args.approved_plan_sha256).lower()
    if bank.get("status") != "verified" or bank.get("request_approval_verified") is not True:
        raise ValueError("combined position-PSF bank is not verified")
    if bank.get("approved_satellite_plan_sha256") != approved_hash:
        raise ValueError("bank audit does not match the explicitly approved satellite plan")
    if summary.get("status") != "complete" or summary.get("science_measurement") is not False:
        raise ValueError("batch is incomplete or incorrectly claims a science measurement")
    if summary.get("n_errors") != 0 or summary.get("n_completed") != len(summary.get("records", [])):
        raise ValueError("batch error/completion totals are inconsistent")
    if summary.get("supported_satellite_psf_request_approval") != {
            "verified": True, "approved_plan_sha256": approved_hash}:
        raise ValueError("batch does not record the exact approved PSF request")
    injection = summary.get("injection", {})
    if (injection.get("psf_convolution_applied") is not True
            or injection.get("detection_rerun") is not True
            or injection.get("original_fits_modified") is not False
            or injection.get("poisson_noise_added") is not False):
        raise ValueError("batch injection flags do not match the requested sensitivity")
    if summary.get("selection", {}).get("frame_selection_mode") != "icl_supported_parent":
        raise ValueError("batch selection is not the ICL-supported parent sample")

    bank_frames = {
        (str(row["cluster_id"]), str(row["filter"]).lower()): row
        for row in bank.get("frames", [])
    }
    if len(bank_frames) != int(bank.get("n_supported_frames", -1)):
        raise ValueError("combined bank frame manifest has duplicate/missing frames")
    psf_hash_by_name = {Path(path).name: digest
                        for path, digest in bank.get("psf_file_sha256", {}).items()}
    if len(psf_hash_by_name) != len(bank.get("psf_file_sha256", {})):
        raise ValueError("combined bank contains ambiguous PSF basenames")

    selected = {(str(row["cluster_id"]), str(row["filter"]).lower())
                for row in summary.get("selection", {}).get("selected_pairs", [])}
    if len(selected) != summary["selection"].get("n_selected_cluster_filter_pairs"):
        raise ValueError("selected pair count differs from the selection manifest")
    if not selected.issubset(bank_frames):
        raise ValueError("batch selected a pair outside the verified PSF bank")

    checked, checked_pairs = [], set()
    for record in summary["records"]:
        pair = (str(record["cluster_id"]), str(record["filter"]).lower())
        if pair in checked_pairs or pair not in selected:
            raise ValueError(f"duplicate or unselected batch record: {pair}")
        checked_pairs.add(pair)
        result_path = Path(record["path"]).resolve()
        if result_path.parent != batch_dir or not result_path.is_file():
            raise ValueError(f"frame result is outside the batch directory: {pair}")
        if sha256(result_path) != record.get("sha256"):
            raise ValueError(f"frame result checksum mismatch: {pair}")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("science_measurement") is not False:
            raise ValueError(f"frame incorrectly claims a science measurement: {pair}")
        gaps = validate_support_accounting(result)
        if record.get("satellite_support_gaps") != gaps:
            raise ValueError(f"support-gap manifest differs from frame result: {pair}")
        if result.get("injection", {}).get("psf_file_count") != len(bank_frames[pair]["psf_names"]):
            raise ValueError(f"PSF file count differs from bank roles: {pair}")
        if (result["injection"].get("source_detection_and_masks_rerun") is not True
                or result["injection"].get("unresolved_bcg_policy") != "mask_all"
                or result["injection"].get("variance_and_poisson_noise_changed") is not False
                or result["injection"].get("psf_negative_pixels_preserved") is not True):
            raise ValueError(f"frame injection flags are inconsistent: {pair}")

        used_psfs = {
            Path(path).name: digest
            for path, digest in result.get("input_sha256", {}).items()
            if (Path(path).suffix.lower() == ".fits"
                and (Path(path).name.startswith("psf_")
                     or Path(path).name == "cluster10_i2_center_psf.fits"))
        }
        expected_names = set(bank_frames[pair]["psf_names"])
        expected_psfs = {name: psf_hash_by_name[name] for name in expected_names}
        if used_psfs != expected_psfs:
            raise ValueError(f"frame PSF files/hashes differ from the approved position roles: {pair}")

        frame_code = result.get("code_sha256", {})
        if not frame_code or any(summary.get("code_sha256", {}).get(path) != digest
                                 for path, digest in frame_code.items()):
            raise ValueError(f"frame code hashes differ from batch provenance: {pair}")
        checked.append({"cluster_id": pair[0], "filter": pair[1],
                        "n_psf_files_verified": len(used_psfs),
                        "satellite_support_gaps": gaps})

    if checked_pairs != selected:
        raise ValueError("frame records do not cover the selected pair manifest")
    gap_counts = {
        scenario: sum(scenario in row["satellite_support_gaps"] for row in checked)
        for scenario in ("baseline", "icl_only", "satellite_only", "both")
    }
    frames_with_gaps = sum(bool(row["satellite_support_gaps"]) for row in checked)
    if (summary.get("n_frames_with_satellite_support_gaps") != frames_with_gaps
            or summary.get("satellite_support_gaps_by_scenario") != gap_counts):
        raise ValueError("batch-level satellite support-gap totals differ from frame records")

    result = {
        "schema": "astrafek-icl-supported-pairs-psf-redetection-audit-v1",
        "status": "verified",
        "science_measurement": False,
        "approved_satellite_plan_sha256": approved_hash,
        "bank_audit_sha256": sha256(bank_path),
        "batch_summary_sha256": sha256(summary_path),
        "n_bank_frames_available": len(bank_frames),
        "n_selected_frames": len(selected),
        "n_frame_results_checked": len(checked),
        "n_errors_reported_by_batch": summary["n_errors"],
        "n_frames_with_satellite_support_gaps": frames_with_gaps,
        "satellite_support_gaps_by_scenario": gap_counts,
        "support_status": "partial_satellite_support" if frames_with_gaps else "complete",
        "records": checked,
        "limitations": [
            "This validates injected-source transfer and redetection bookkeeping; it is not an observational science measurement.",
            "Only compact HSC model-PSF cores are included; extended PSF wings are not validated.",
            "Injected-source Poisson noise and HSC upstream sky transfer are not modeled.",
        ],
    }
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_bank_frames_available", "n_selected_frames",
        "n_frame_results_checked", "n_errors_reported_by_batch",
        "n_frames_with_satellite_support_gaps", "support_status", "science_measurement")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
