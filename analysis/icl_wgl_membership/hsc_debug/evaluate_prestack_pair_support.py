#!/usr/bin/env python3
"""Screen common-area support for existing matched cluster/control profiles.

This is a pilot support analysis only. It never changes the established 50%
coverage cut or marks the science measurement as approved.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from prestack_core import FACTORS, covariance_diagnostics, spatial_labels


THRESHOLDS = (0.25, 0.35, 0.40, 0.50)
SCOPES = {
    "50-150": [0, 1],
    "50-200": [0, 1, 2],
    "50-300": [0, 1, 2, 3],
    "50-500": [0, 1, 2, 3, 4, 5],
    **{f"{lo}-{hi}": [i] for i, (lo, hi) in enumerate(
        zip((50, 100, 150, 200, 300, 400), (100, 150, 200, 300, 400, 500))
    )},
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def spatial_cell(ra, dec):
    return str(spatial_labels([ra], [dec], 5.0)[0])


def pair_block_labels(target_rows, control_rows, selected):
    """Join pairs sharing either a target or control 5-degree sky cell."""
    selected = np.flatnonzero(selected)
    parent = {int(i): int(i) for i in selected}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    owner = {}
    for i in selected:
        target_cell = spatial_cell(target_rows[i]["ra_deg"], target_rows[i]["dec_deg"])
        control_cell = spatial_cell(control_rows[i]["ra_deg"], control_rows[i]["dec_deg"])
        for key in (("target", target_cell), ("control", control_cell)):
            if key in owner:
                union(int(i), owner[key])
            else:
                owner[key] = int(i)
    roots = {root: label for label, root in enumerate(sorted({find(int(i)) for i in selected}))}
    labels = np.array([f"pair_block_{roots[find(int(i))]}" for i in selected])
    target_cells = {spatial_cell(target_rows[i]["ra_deg"], target_rows[i]["dec_deg"]) for i in selected}
    control_cells = {spatial_cell(control_rows[i]["ra_deg"], control_rows[i]["dec_deg"]) for i in selected}
    block_sizes = {}
    for label in labels:
        block_sizes[label] = block_sizes.get(label, 0) + 1
    return selected, labels, len(target_cells), len(control_cells), len(block_sizes), max(block_sizes.values(), default=0)


def recovery_sensitivity(path):
    by_group = {}
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["cohort"] != "area_qualified" or row["extraction"] != "raw":
                continue
            key = (float(row["kron_factor"]), int(row["model_index"]), row["control_id"])
            by_group.setdefault(key, {})[int(float(row["r_min_kpc"]))] = row
    result = []
    bins = (50, 100, 150, 200)
    for factor in FACTORS:
        for model in range(6):
            for threshold in THRESHOLDS:
                eligible = []
                for (kron, model_index, _), rows in by_group.items():
                    if kron != factor or model_index != model or not all(b in rows for b in bins):
                        continue
                    if all(float(rows[b]["before_coverage"]) >= threshold and
                           float(rows[b]["after_coverage"]) >= threshold for b in bins):
                        eligible.append(rows)
                bias = []
                mask_aware = []
                for b in bins:
                    values = [float(rows[b]["recovery_fraction"]) for rows in eligible]
                    aware = [float(rows[b]["mask_aware_recovery_fraction"]) for rows in eligible]
                    bias.append(float(np.median(values) - 1) if values else None)
                    mask_aware.append(float(np.median(aware) - 1) if aware else None)
                result.append(dict(kron_factor=factor, model_index=model, coverage_floor=threshold,
                                   n_complete_injection_fields=len(eligible),
                                   median_fractional_bias_50_100_to_200_300=bias,
                                   mask_aware_median_fractional_bias=mask_aware,
                                   status="sensitivity_only"))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output directory: {args.work}")

    validation_path = args.sample / "validation_sample_2026-10-13_01/validation_manifest.json"
    control_path = args.sample / "control_patches_2026-10-14_02/control_download_manifest.json"
    extraction_summary_path = args.extraction / "extraction_summary.json"
    cluster_summary_path = args.clusters / "summary.json"
    cluster_profiles_path = args.clusters / "profiles.npz"
    recovery_path = args.analysis / "recovery.tsv"
    validation_rows = read_json(validation_path)["clusters"]
    control_rows = read_json(control_path)["clusters"]
    extraction_summary = read_json(extraction_summary_path)
    cluster_summary = read_json(cluster_summary_path)
    if len(validation_rows) != 50 or len(cluster_summary["records"]) != 50:
        raise ValueError("This pilot expects the existing 50 validation pairs")
    if cluster_summary.get("status") != "complete" or sha256(cluster_profiles_path) != cluster_summary["output_sha256"]:
        raise ValueError("Cluster profile archive is incomplete or changed")

    target_to_control = {}
    for row in control_rows:
        key = str(row["target_cluster_id"])
        if key in target_to_control:
            raise ValueError(f"Multiple controls for validation target {key}")
        target_to_control[key] = row
    if [str(row["id"]) for row in validation_rows] != [str(row["id"]) for row in cluster_summary["records"]]:
        raise ValueError("Cluster profile order differs from the validation manifest")

    controls = []
    control_ids = []
    input_hashes = {
        str(path): sha256(path) for path in
        (validation_path, control_path, extraction_summary_path, cluster_summary_path,
         cluster_profiles_path, recovery_path)
    }
    extraction_records = {row["control_id"]: row for row in extraction_summary.get("records", [])}
    for target in validation_rows:
        control = target_to_control.get(str(target["id"]))
        if control is None:
            raise ValueError(f"No matched control for target {target['id']}")
        control_id = control["id"]
        baseline = args.extraction / control_id / "baseline.npz"
        metadata_path = args.extraction / control_id / "metadata.json"
        metadata = read_json(metadata_path)
        digest = sha256(baseline)
        extraction_record = extraction_records.get(control_id)
        if metadata["target_cluster_id"] != str(target["id"]):
            raise ValueError(f"Control target mismatch for {control_id}")
        if extraction_record is None or extraction_record.get("baseline_sha256") != digest or metadata["baseline_sha256"] != digest:
            raise ValueError(f"Control baseline hash mismatch for {control_id}")
        input_hashes[str(baseline)] = digest
        input_hashes[str(metadata_path)] = sha256(metadata_path)
        with np.load(baseline) as archive:
            controls.append({key: archive[key] for key in ("coverage", "raw", "nominal_raw")})
        control_ids.append(control_id)

    cluster_records = cluster_summary["records"]
    cluster_ra = np.array([row["ra_deg"] for row in cluster_records])
    cluster_dec = np.array([row["dec_deg"] for row in cluster_records])
    with np.load(cluster_profiles_path) as archive:
        cluster_coverage = archive["coverage"]
        cluster_raw = archive["raw"]
        cluster_nominal = archive["nominal_raw"]
    control_coverage = np.array([row["coverage"] for row in controls])
    control_raw = np.array([row["raw"] for row in controls])
    control_nominal = np.array([row["nominal_raw"] for row in controls])
    if control_coverage.shape != cluster_coverage.shape or control_coverage.shape != (50, 3, 6):
        raise ValueError("Unexpected profile shapes")

    full_rows = []
    covariance_archive = {}
    for factor_index, factor in enumerate(FACTORS):
        for threshold_index, threshold in enumerate(THRESHOLDS):
            for scope_index, (scope, bins) in enumerate(SCOPES.items()):
                bins_array = np.asarray(bins, dtype=int)
                values = cluster_raw[:, factor_index, bins_array] - control_raw[:, factor_index, bins_array]
                nominal = (cluster_nominal[:, factor_index] + control_nominal[:, factor_index])
                nominal = nominal[:, bins_array][:, :, bins_array]
                coverage_ok = ((cluster_coverage[:, factor_index, bins_array] >= threshold) &
                               (control_coverage[:, factor_index, bins_array] >= threshold))
                selected = np.all(coverage_ok, axis=1) & np.isfinite(values).all(axis=1)
                selected &= np.isfinite(nominal).all(axis=(1, 2))
                selected_indices, labels, cluster_cells, control_cells, blocks, largest_block = pair_block_labels(
                    cluster_records, control_rows_for_targets(validation_rows, target_to_control), selected)
                info, matrices = covariance_diagnostics(
                    values[selected], labels, nominal[selected],
                    seed=20261016 + factor_index * 10000 + threshold_index * 100 + scope_index,
                    replicas=2000,
                )
                key = f"k{factor}_c{int(threshold*100)}_r{scope.replace('-', '_')}"
                for name, matrix in matrices.items():
                    covariance_archive[f"{key}__{name}"] = matrix
                covariance_archive[f"{key}__selected_indices"] = selected_indices
                cohorts = [validation_rows[i]["preparation_cohort"] for i in selected_indices]
                bootstrap_ratio = info.get("iid_bootstrap_diagonal_ratio")
                folds = info.get("folds", [])
                pilot_gate = bool(
                    info.get("covariance_sample_size_gate") and
                    info.get("rank") == len(bins) and bootstrap_ratio and
                    np.isfinite(bootstrap_ratio).all() and
                    np.max(np.abs(np.asarray(bootstrap_ratio) - 1)) <= 0.10 and
                    len(folds) == 5 and info.get("cv_fraction_within_1p96", 0) >= 0.80
                )
                full_rows.append(dict(
                    name=key, kron_factor=factor, coverage_floor=threshold, radial_scope_kpc=scope,
                    bins=",".join(str(i) for i in bins), n_pairs=int(selected.sum()),
                    n_area_qualified=sum(value == "area_qualified" for value in cohorts),
                    n_mask_stress=sum(value != "area_qualified" for value in cohorts),
                    n_cluster_sky_cells=cluster_cells, n_control_sky_cells=control_cells,
                    n_joint_spatial_blocks=blocks, largest_joint_block=largest_block,
                    n_bins=len(bins), rank=info.get("rank"), condition_number=info.get("condition_number"),
                    minimum_eigenvalue=info.get("minimum_eigenvalue"),
                    sample_size_gate=info.get("covariance_sample_size_gate", False),
                    pilot_covariance_diagnostic_gate=pilot_gate,
                    iid_bootstrap_diagonal_ratio=json.dumps(bootstrap_ratio),
                    cv_fraction_within_1p96=info.get("cv_fraction_within_1p96"),
                    n_complete_spatial_folds=len(folds),
                    empirical_to_nominal_variance_ratio=json.dumps(info.get("empirical_to_nominal_variance_ratio")),
                    mean_pair_difference=json.dumps(
                        matrices["mean"].tolist() if "mean" in matrices else None),
                    independent_stack_standard_error=json.dumps(
                        np.sqrt(np.diag(matrices["independent_stack_covariance"])).tolist()
                        if "independent_stack_covariance" in matrices else None),
                    joint_block_standard_error=json.dumps(
                        np.sqrt(np.diag(matrices["spatial_bootstrap_covariance"])).tolist()
                        if "spatial_bootstrap_covariance" in matrices and
                        np.isfinite(matrices["spatial_bootstrap_covariance"]).all() else None),
                    status=info["status"], science_measurement=False,
                ))

    args.work.mkdir(parents=True, exist_ok=False)
    support_path = args.work / "support_counts.tsv"
    covariance_path = args.work / "pilot_covariances.npz"
    recovery_output = args.work / "injection_recovery_sensitivity.tsv"
    write_tsv(support_path, full_rows)
    np.savez_compressed(covariance_path, **covariance_archive)
    recovery_rows = recovery_sensitivity(recovery_path)
    write_tsv(recovery_output, recovery_rows)
    program_path = Path(__file__).resolve()
    summary = dict(
        status="complete", validation_pairs=len(validation_rows), matched_controls=len(control_ids),
        thresholds=list(THRESHOLDS), radial_scopes_kpc=SCOPES,
        covariance_configs=len(full_rows), injection_sensitivity_rows=len(recovery_rows),
        purpose="Exploratory paired-support sensitivity using existing 50 validation cluster/control pairs.",
        limitations=[
            "The 50 validation targets are a selected pilot cohort, not the full 800-cluster science sample.",
            "Coverage floors below 50% and alternate Kron factors are sensitivity cases, not approved cuts.",
            "The 5-degree joint blocks link pairs sharing either target or control sky cells; residual spatial dependence remains possible.",
            "Paired profile scatter includes cluster diversity and control/background variation; it is not pure pixel noise.",
            "Injection recovery is a coadd-level smooth-signal test and does not validate PSF wings, upstream sky subtraction, or unmasked-source contamination.",
        ],
        primary_definition_n=next(row["n_pairs"] for row in full_rows
                                  if row["kron_factor"] == 3.5 and row["coverage_floor"] == 0.5
                                  and row["radial_scope_kpc"] == "50-500"),
        science_measurement=False, sky_correction_approved=False, common_pair_covariance_approved=False,
        source_sha256=input_hashes,
        output_sha256={str(path): sha256(path) for path in (support_path, covariance_path, recovery_output)},
        script_sha256={str(program_path): sha256(program_path)},
    )
    write_json(args.work / "summary.json", summary)
    print(f"Pair-support screen complete: {len(full_rows)} covariance configurations; "
          f"{len(recovery_rows)} injection sensitivity groups", flush=True)


def control_rows_for_targets(target_rows, target_to_control):
    return [target_to_control[str(row["id"])] for row in target_rows]


if __name__ == "__main__":
    main()
