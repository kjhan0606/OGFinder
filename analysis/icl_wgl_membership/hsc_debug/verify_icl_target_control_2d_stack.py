#!/usr/bin/env python3
"""Independently verify a completed exploratory HSC ICL 2-D stack."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from astropy.io import fits


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def resolve_mask_factor(summary: dict) -> float:
    """Resolve the target/control mask factor, including repaired legacy JSON.

    One historical 2.5 stack summary omitted the numeric ``mask_factor``.
    Its immutable stack definition still records both the target mask HDU and
    the matched-control Kron factor. Accept that fallback only when both
    declarations agree exactly.
    """
    value = summary.get("mask_factor")
    if value is not None:
        factor = float(value)
    else:
        definition = str(summary.get("stack_definition", {}).get("mask", ""))
        target = re.search(r"\bMASK_K(\d+)P(\d+)\b", definition)
        control = re.search(r"\bKron\s+(\d+(?:\.\d+)?)\b", definition)
        if target is None or control is None:
            raise ValueError("summary has no numeric mask factor or recoverable mask definition")
        target_factor = float(f"{target.group(1)}.{target.group(2)}")
        control_factor = float(control.group(1))
        if target_factor != control_factor:
            raise ValueError("target and matched-control mask definitions disagree")
        factor = target_factor
    if factor not in (2.5, 3.5, 5.0):
        raise ValueError(f"unsupported mask factor: {factor}")
    return factor


def candidate_audit_pairs(audit: dict, controls_by_target: dict,
                         control_manifest_sha: str, control_status_sha: str) -> dict[str, tuple[str, str, list[str]]]:
    """Independently reconstruct the audited candidate/control pair exclusions."""
    require = lambda test, message: (_ for _ in ()).throw(ValueError(message)) if not test else None
    require(audit.get("schema") == "astrafek-s23b-candidate-existing-cutout-support-v1"
            and audit.get("status") == "header_wcs_audit_complete"
            and audit.get("science_measurement") is False
            and audit.get("hsc_query_submitted") is False
            and audit.get("data_downloaded") is False
            and audit.get("image_pixels_read") is False,
            "candidate-control audit schema or scope is invalid")
    require(audit.get("n_candidates") == 1333
            and audit.get("candidate_table_sha256")
            == "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c",
            "candidate-control audit does not use the pinned S23B table")
    target = audit.get("target", {})
    require(target.get("n_images") == 800
            and target.get("n_images_with_candidate_centers") == 0
            and target.get("n_unique_candidates_inside") == 0
            and target.get("matches") == [],
            "target candidate overlaps are outside this control-only sensitivity audit")
    control_result = audit.get("matched_controls", {})
    hashes = control_result.get("inputs_sha256", {})
    require(hashes.get("manifest") == control_manifest_sha
            and hashes.get("download_status") == control_status_sha,
            "candidate-control audit input checksums differ from the stack inputs")
    matches = control_result.get("matches")
    require(isinstance(matches, list) and control_result.get("n_images") == 798
            and control_result.get("n_images_with_candidate_centers") == len(matches),
            "candidate-control audit match counts are inconsistent")
    excluded = {}
    all_candidates = set()
    for match in matches:
        target_id = str(match.get("target_cluster_id", ""))
        control_id = str(match.get("frame_id", ""))
        control = controls_by_target.get(target_id)
        candidates = match.get("candidates")
        candidate_ids = [str(row.get("candidate_id", "")) for row in candidates or []]
        require(target_id not in excluded and bool(target_id) and bool(control_id)
                and control is not None and str(control.get("id")) == control_id
                and str(control.get("coadd_header_filter", "")).lower() == match.get("filter"),
                "candidate-control audit references an unknown or duplicate matched pair")
        require(bool(candidate_ids) and all(candidate_ids)
                and len(set(candidate_ids)) == len(candidate_ids)
                and match.get("n_candidate_centers_inside") == len(candidate_ids),
                "candidate-control audit candidate rows are malformed")
        excluded[target_id] = (control_id, str(match["filter"]), candidate_ids)
        all_candidates.update(candidate_ids)
    require(control_result.get("n_unique_candidates_inside") == len(all_candidates),
            "candidate-control audit unique candidate count is inconsistent")
    return excluded


def git_version_matches(path: Path, expected_sha256: str, repo_root: Path) -> bool:
    """Check whether an exact historical source version is retained in Git."""
    try:
        relative = path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return False
    commits = subprocess.run(
        ["git", "-C", str(repo_root), "log", "--all", "--format=%H", "--", relative],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    for commit in commits:
        blob = subprocess.run(
            ["git", "-C", str(repo_root), "show", f"{commit}:{relative}"],
            capture_output=True,
        )
        if blob.returncode == 0 and hashlib.sha256(blob.stdout).hexdigest() == expected_sha256:
            return True
    return False


def verify(sample: Path, work: Path) -> dict:
    summary_path = work / "summary.json"
    summary = read_json(summary_path)
    require = lambda test, message: (_ for _ in ()).throw(ValueError(message)) if not test else None
    require(summary.get("schema") == "astrafek-hsc-icl-target-control-2d-stack-v1",
            "unexpected stack summary schema")
    require(summary.get("status") == "complete_exploratory", "stack is not complete")
    require(summary.get("science_measurement") is False, "science gate must remain false")
    mask_factor = resolve_mask_factor(summary)

    sample_manifest = sample / "sample_manifest.json"
    mask_summary = sample / "prestack_cluster_masks_2026-10-16_01/summary.json"
    mask_records_path = sample / "prestack_cluster_masks_2026-10-16_01/records.tsv"
    control_manifest = sample / "control_download_01/download_manifest.json"
    control_status = sample / "control_download_01/download_status.json"
    for path in (sample_manifest, mask_summary, mask_records_path, control_manifest, control_status):
        require(path.is_file(), f"missing provenance input: {path}")
        require(summary["input_sha256"].get(str(path)) == sha256(path),
                f"provenance checksum mismatch: {path}")

    targets = read_json(sample_manifest)["clusters"]
    masks = {str(row["id"]): row for row in read_tsv(mask_records_path)}
    controls = read_json(control_manifest)["clusters"]
    downloads = read_json(control_status)
    require(len(targets) == 800 and len(masks) == 800 and len(controls) == 798,
            "input frame counts differ from the approved sample")
    require(downloads.get("status") == "complete" and downloads.get("downloaded") == 798,
            "matched-control downloads are incomplete")

    target_by_id = {str(row["id"]): row for row in targets}
    controls_by_target = {str(row["target_cluster_id"]): row for row in controls}
    require(len(controls_by_target) == len(controls), "matched-control target mapping is duplicated")
    eligible = {}
    excluded = []
    for control in controls:
        tid = str(control["target_cluster_id"])
        require(tid in target_by_id and tid in masks, f"unknown matched target: {tid}")
        target_filter = masks[tid]["primary_filter"]
        control_filter = control.get("coadd_header_filter")
        if target_filter != control_filter:
            excluded.append((tid, str(control["id"]), target_filter, control_filter))
        else:
            eligible[tid] = (str(control["id"]), control_filter)
    observed_excluded = {
        (str(row["target_id"]), str(row["control_id"]), row["target_filter"], row["control_filter"])
        for row in summary.get("filter_mismatch_exclusions", [])
    }
    require(observed_excluded == set(excluded), "filter-mismatch exclusions differ from input metadata")
    require(summary.get("n_pairs_excluded_filter_mismatch") == len(excluded),
            "filter-mismatch count is wrong")

    audit_ref = summary.get("candidate_control_audit")
    if audit_ref is None:
        candidate_audit_excluded = {}
    else:
        audit_path = Path(audit_ref["path"])
        require(audit_ref.get("schema") == "astrafek-s23b-candidate-existing-cutout-support-v1"
                and audit_path.is_file()
                and sha256(audit_path) == audit_ref.get("sha256")
                and summary["input_sha256"].get(str(audit_path)) == audit_ref.get("sha256"),
                "candidate-control audit checksum or provenance is invalid")
        audit = read_json(audit_path)
        candidate_audit_excluded = candidate_audit_pairs(
            audit, controls_by_target, sha256(control_manifest), sha256(control_status))
    filter_excluded_by_target = {row[0]: row for row in excluded}
    candidate_excluded = {tid: pair for tid, pair in candidate_audit_excluded.items()
                          if tid not in filter_excluded_by_target}
    candidate_already_filter_excluded = {
        tid: pair for tid, pair in candidate_audit_excluded.items()
        if tid in filter_excluded_by_target
    }
    for tid, (control_id, band, _) in candidate_already_filter_excluded.items():
        observed_filter = filter_excluded_by_target[tid]
        require(observed_filter[1] == control_id and observed_filter[2] != band,
                "candidate-control audit overlaps an inconsistent filter exclusion")
    observed_candidate = {}
    for row in summary.get("s23b_candidate_control_exclusions", []):
        tid = str(row.get("target_id", ""))
        require(tid not in observed_candidate, "candidate-control exclusion target is duplicated")
        observed_candidate[tid] = (str(row.get("control_id", "")), str(row.get("filter", "")),
                                   [str(value) for value in row.get("candidate_ids", [])])
    require(observed_candidate == candidate_excluded,
            "stack candidate-control exclusions differ from independent audit reconstruction")
    require(summary.get("n_pairs_excluded_s23b_candidate_control", 0) == len(candidate_excluded),
            "candidate-control exclusion count is wrong")
    observed_already_filter_excluded = {}
    for row in summary.get("candidate_controls_already_filter_excluded", []):
        tid = str(row.get("target_id", ""))
        require(tid not in observed_already_filter_excluded,
                "already-filter-excluded candidate control is duplicated")
        observed_already_filter_excluded[tid] = (
            str(row.get("control_id", "")), str(row.get("filter", "")),
            [str(value) for value in row.get("candidate_ids", [])])
    expected_already_filter_excluded = {
        tid: (control_id, band, candidate_ids)
        for tid, (control_id, band, candidate_ids) in candidate_already_filter_excluded.items()
    }
    require(observed_already_filter_excluded == expected_already_filter_excluded
            and summary.get("n_candidate_controls_already_filter_excluded", 0)
            == len(expected_already_filter_excluded)
            and summary.get("n_candidate_control_audit_matches", 0)
            == len(candidate_audit_excluded),
            "already-filter-excluded candidate controls differ from the audit")
    for tid, (control_id, band, candidate_ids) in candidate_excluded.items():
        require(tid in eligible and eligible[tid] == (control_id, band),
                "candidate-control exclusion does not map to an otherwise eligible pair")
        del eligible[tid]
    require(summary.get("n_pairs_stacked") == len(eligible), "eligible pair count is wrong")

    pair_path = work / "paired_radial_profiles.tsv"
    radial_path = work / "stacked_radial_profiles.tsv"
    require(summary["outputs"]["paired_profiles_sha256"] == sha256(pair_path),
            "paired radial-profile checksum mismatch")
    require(summary["outputs"]["stacked_profiles_sha256"] == sha256(radial_path),
            "stacked radial-profile checksum mismatch")
    pair_rows = read_tsv(pair_path)
    by_target = {row["target_id"]: row for row in pair_rows}
    require(len(by_target) == len(pair_rows) == len(eligible), "paired profile rows are duplicated or incomplete")
    for tid, (control_id, band) in eligible.items():
        row = by_target.get(tid)
        require(row is not None and row["control_id"] == control_id and row["filter"] == band,
                f"pair table mismatch for target {tid}")
        profile = json.loads(row["profile_json"])
        require(len(profile) == 20, f"target {tid} does not have 20 radial bins")
        for annulus in profile:
            if annulus["delta_mean_njy_arcsec2"] is not None:
                require(np.isclose(annulus["target_mean_njy_arcsec2"]
                                   - annulus["control_mean_njy_arcsec2"],
                                   annulus["delta_mean_njy_arcsec2"], rtol=1e-6, atol=1e-6),
                        f"target-control radial subtraction mismatch for {tid}")
    radial_rows = read_tsv(radial_path)
    require(len(radial_rows) == 40, "stacked radial profile must have 20 bins per filter")

    map_checks = {}
    for band in ("i", "i2"):
        product = summary["outputs"]["maps"][band]
        path = Path(product["file"])
        require(path.is_file() and sha256(path) == product["sha256"],
                f"{band} FITS output checksum mismatch")
        with fits.open(path, memmap=False, checksum=True) as hdus:
            required = {"TARGET", "CONTROL", "DELTA", "NPAIR", "TARGET_SEM",
                        "CONTROL_SEM", "DELTA_SEM", "GRID_KPC"}
            require(required.issubset({h.name for h in hdus}), f"{band} FITS extensions are incomplete")
            require(hdus[0].header.get("SCIVAL") is False, f"{band} FITS science gate is not false")
            require(hdus[0].header.get("FILTER") == band, f"{band} FITS filter header mismatch")
            maskfac_header = hdus[0].header.get("MASKFAC")
            if maskfac_header is not None:
                require(float(maskfac_header) == mask_factor,
                        f"{band} FITS mask-factor header mismatch")
            else:
                # The first 2.5 run predates MASKFAC; retain verification only
                # when both its human-readable FITS card and summary definition
                # identify the same factor.
                mask_header = str(hdus[0].header.get("MASK", ""))
                require(f"Kron-{mask_factor:g}" in mask_header,
                        f"{band} FITS header does not identify the mask factor")
            target, control = np.asarray(hdus["TARGET"].data), np.asarray(hdus["CONTROL"].data)
            delta, count = np.asarray(hdus["DELTA"].data), np.asarray(hdus["NPAIR"].data)
            require(target.shape == control.shape == delta.shape == count.shape == (201, 201),
                    f"{band} FITS image shapes are unexpected")
            require(int(count.max()) <= int(summary["n_by_filter"][band]),
                    f"{band} pixel pair count exceeds number of pairs")
            common = np.isfinite(target) & np.isfinite(control) & np.isfinite(delta)
            require(np.array_equal(common, count > 0), f"{band} map coverage and pair-count masks differ")
            residual = delta[common] - (target[common] - control[common])
            max_residual = float(np.max(np.abs(residual))) if residual.size else 0.0
            scale = max(1.0, float(np.max(np.abs(delta[common]))) if residual.size else 1.0)
            require(max_residual <= 2e-5 * scale,
                    f"{band} mean maps do not satisfy target minus control")
            require(int(common.sum()) == product["covered_pixels"], f"{band} covered-pixel count mismatch")
            map_checks[band] = {"n_pairs": int(summary["n_by_filter"][band]),
                                "covered_pixels": int(common.sum()),
                                "max_delta_identity_residual": max_residual}

    code_provenance = {}
    repo_root = Path(__file__).resolve().parents[3]
    for name, expected_hash in summary.get("script_sha256", {}).items():
        path = Path(name)
        current_matches = path.is_file() and sha256(path) == expected_hash
        snapshot_matches = (
            name == str((Path(__file__).resolve().parent / "stack_hsc_icl_target_control_2d.py").resolve())
            and summary.get("executed_script_snapshot_sha256") == expected_hash
            and Path(summary.get("executed_script_snapshot", "")).is_file()
            and sha256(Path(summary["executed_script_snapshot"])) == expected_hash
        )
        history_matches = not (current_matches or snapshot_matches) and git_version_matches(
            path, expected_hash, repo_root)
        require(current_matches or snapshot_matches or history_matches,
                f"analysis script changed since run and no matching run snapshot exists: {path}")
        code_provenance[name] = (
            "current" if current_matches else "run_snapshot" if snapshot_matches else "git_history"
        )
    if "serialization_correction" in summary:
        correction = summary["serialization_correction"]
        require(summary.get("science_measurement") is False
                and correction.get("data_products_modified") is False
                and correction.get("executed_script_snapshot_sha256")
                == summary.get("executed_script_snapshot_sha256"),
                "summary serialization correction record is invalid")
        repair_script = Path(correction.get("repair_script", ""))
        require(repair_script.is_file()
                and sha256(repair_script) == correction.get("repair_script_sha256"),
                "summary serialization repair provenance is invalid")
    return {
        "schema": "astrafek-hsc-icl-target-control-2d-stack-verification-v1",
        "status": "verified",
        "summary_sha256": sha256(summary_path),
        "n_input_controls": len(controls), "n_filter_mismatch_excluded": len(excluded),
        "n_s23b_candidate_control_excluded": len(candidate_excluded),
        "n_eligible_pairs": len(eligible), "n_pair_rows": len(pair_rows),
        "n_stacked_radial_rows": len(radial_rows), "map_checks": map_checks,
        "mask_factor": mask_factor, "code_provenance": code_provenance,
        "science_measurement": False,
        "input_images_rehashed": False,
        "note": "The processing script rehashed every input frame; this verifier independently rechecks manifests, outputs, maps, pairing and calculations without another full-image read.",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = verify(args.sample.resolve(), args.work.resolve())
    output = args.output.resolve() if args.output else args.work.resolve() / "independent_verification.json"
    if output.exists():
        parser.error(f"refusing to overwrite existing verification: {output}")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"verification_sha256={sha256(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
