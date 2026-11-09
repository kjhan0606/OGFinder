#!/usr/bin/env python3
"""Audit selection into the i2 ICL–satellite–random-subtracted-WL overlap.

Uses the existing hash-pinned 68-lens ICL/satellite table and WL-eligibility
audit. It makes no network request and does not read image pixels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PAIRS = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
OVERLAP_DIR = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03"
ELIGIBILITY = OVERLAP_DIR / "supported_pair_wl_eligibility.tsv"
OVERLAP_SUMMARY = OVERLAP_DIR / "summary.json"
DEFAULT_OUTPUT = DATA / "icl_wl_overlap_selection_i2_2026-11-06_02"
EXPECTED = {
    PAIRS: "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da",
    ELIGIBILITY: "68630ea51174de5feeb3bb8571da55beae4b8816dfc55ffe0d32d21479ad4ae9",
    OVERLAP_SUMMARY: "98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12",
}
REPLICATES = 20000
SEED = 20261108


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def standardized_difference(values: np.ndarray, selected: np.ndarray) -> float:
    a, b = values[selected], values[~selected]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = math.sqrt(((len(a) - 1) * np.var(a, ddof=1)
                        + (len(b) - 1) * np.var(b, ddof=1))
                       / (len(a) + len(b) - 2))
    return float((np.mean(a) - np.mean(b)) / pooled) if pooled > 0 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=PAIRS)
    parser.add_argument("--eligibility", type=Path, default=ELIGIBILITY)
    parser.add_argument("--overlap-summary", type=Path, default=OVERLAP_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"Refusing to overwrite existing output: {args.output}")
    if args.replicates < 1000:
        parser.error("At least 1,000 restricted-block permutations are required")

    paths = {"pairs": args.pairs, "eligibility": args.eligibility,
             "overlap_summary": args.overlap_summary}
    hashes = {key: sha256(path) for key, path in paths.items()}
    for key, path in paths.items():
        pinned = {"pairs": PAIRS, "eligibility": ELIGIBILITY,
                  "overlap_summary": OVERLAP_SUMMARY}[key]
        if path.resolve() == pinned.resolve() and hashes[key] != EXPECTED[path.resolve()]:
            raise ValueError(f"Pinned input hash mismatch: {path}")

    audit = json.loads(args.overlap_summary.read_text(encoding="utf-8"))
    if (audit.get("schema") != "astrafek-icl-wl-sample-overlap-audit-v1"
            or audit.get("status") != "complete_exploratory"
            or audit.get("science_measurement") is not False
            or audit.get("network_used") is not False
            or audit.get("eligibility_table_sha256") != hashes["eligibility"]):
        raise ValueError("Overlap audit status or eligibility provenance changed")
    expected = audit["joint_icl_satellite_support_by_filter"]["i2"]
    pair_rows = [r for r in read_tsv(args.pairs) if r["filter"] == "i2"]
    eligibility_rows = [r for r in read_tsv(args.eligibility) if r["filter"] == "i2"]
    pair_by_id = {r["cluster_id"]: r for r in pair_rows}
    elig_by_id = {r["cluster_id"]: r for r in eligibility_rows}
    if len(pair_by_id) != 68 or len(elig_by_id) != 68 or set(pair_by_id) != set(elig_by_id):
        raise ValueError("Expected exact one-to-one join of the 68 i2 supported pairs")
    if any(r.get("science_measurement", "False").lower() != "false" for r in pair_rows):
        raise ValueError("Input support table contains science_measurement promotion")
    common_ids = set(expected["cluster_ids_in_common_support"])
    if len(common_ids) != 27:
        raise ValueError(f"Expected 27 verified common-support IDs, got {len(common_ids)}")

    rows = []
    for cluster_id in sorted(pair_by_id, key=int):
        pair, eligible = pair_by_id[cluster_id], elig_by_id[cluster_id]
        for key, pkey, ekey in (
            ("redshift", "redshift", "redshift"),
            ("richness", "richness", "richness"),
            ("icl_delta_200_300", "icl_delta_200_300", "icl_delta_200_300"),
            ("satellite_delta_50_75", "satellite_delta_50_75", "satellite_delta_50_75"),
        ):
            if not math.isclose(float(pair[pkey]), float(eligible[ekey]), rel_tol=1e-10, abs_tol=1e-9):
                raise ValueError(f"Input mismatch for {cluster_id}/{key}")
        rows.append({**pair, **{f"wl_{k}": v for k, v in eligible.items()
                                if k not in ("cluster_id", "filter")}})

    flags = ("wl_matched_sample", "target_shape_rows",
             "positive_random_shape_count", "random_subtracted_wl_common_support")
    stage_counts = {flag: sum(r[f"wl_{flag}"] == "True" for r in rows)
                    for flag in flags}
    if (stage_counts["wl_matched_sample"] != expected["n_in_wl_matched_620"]
            or stage_counts["target_shape_rows"] != expected["n_with_target_shape_rows"]
            or stage_counts["positive_random_shape_count"] != expected["n_with_positive_random_shape_counts"]
            or stage_counts["random_subtracted_wl_common_support"]
            != expected["n_in_verified_random_subtracted_common_support"]):
        raise ValueError(f"Eligibility stage counts differ from the audited summary: {stage_counts}")
    selected = np.asarray([r["wl_random_subtracted_wl_common_support"] == "True" for r in rows])
    if set(r["cluster_id"] for r, take in zip(rows, selected) if take) != common_ids:
        raise ValueError("Eligibility flags do not reproduce the pinned 27-lens common sample")

    blocks: dict[str, list[int]] = {}
    for i, row in enumerate(rows):
        blocks.setdefault(row["connected_5deg_block"], []).append(i)
    mixed = [np.asarray(ix, dtype=int) for ix in blocks.values()
             if 0 < int(selected[ix].sum()) < len(ix)]
    if len(mixed) < 5:
        raise ValueError(f"Too few mixed spatial blocks for restricted selection audit: {len(mixed)}")

    variables = {
        "redshift": np.asarray([float(r["redshift"]) for r in rows]),
        "log_richness": np.log([float(r["richness"]) for r in rows]),
        "icl_delta_200_300": np.asarray([float(r["icl_delta_200_300"]) for r in rows]),
        "satellite_delta_50_75": np.asarray([float(r["satellite_delta_50_75"]) for r in rows]),
    }
    # Shape-pair count is a direct component of the WL support rule. Record it
    # as the selection mechanism, but do not mix this tautological variable
    # into the family of physical/covariate balance tests.
    target_pair_rows = np.asarray([float(r["wl_target_pair_rows"]) for r in rows])
    log_target_pair_rows = np.log1p(target_pair_rows)
    pair_count_mechanism = {
        "selected_median_target_pairs": float(np.median(target_pair_rows[selected])),
        "excluded_median_target_pairs": float(np.median(target_pair_rows[~selected])),
        "selected_mean_log1p_target_pairs": float(np.mean(log_target_pair_rows[selected])),
        "excluded_mean_log1p_target_pairs": float(np.mean(log_target_pair_rows[~selected])),
        "log1p_target_pair_rows_standardized_mean_difference":
            standardized_difference(log_target_pair_rows, selected),
        "interpretation": "expected structural association because positive target shape support is part of the common-WL eligibility gate",
    }
    observed = {name: standardized_difference(values, selected)
                for name, values in variables.items()}
    if any(not math.isfinite(v) for v in observed.values()):
        raise ValueError("Non-finite observed selection balance statistic")

    rng = np.random.default_rng(SEED)
    max_abs = np.empty(args.replicates, dtype=float)
    raw_abs = {name: np.empty(args.replicates, dtype=float) for name in variables}
    for b in range(args.replicates):
        permuted = selected.copy()
        for ix in mixed:
            n_selected = int(selected[ix].sum())
            permuted[ix] = False
            permuted[rng.choice(ix, size=n_selected, replace=False)] = True
        stats = {name: standardized_difference(values, permuted)
                 for name, values in variables.items()}
        for name, value in stats.items():
            raw_abs[name][b] = abs(value)
        max_abs[b] = max(abs(value) for value in stats.values())

    output_rows = []
    for name, values in variables.items():
        within = (1 + sum(raw_abs[name] >= abs(observed[name]))) / (args.replicates + 1)
        family = (1 + sum(max_abs >= abs(observed[name]))) / (args.replicates + 1)
        sel = values[selected]
        exc = values[~selected]
        output_rows.append({
            "covariate": name, "n_common_wl": len(sel), "n_no_common_wl": len(exc),
            "selected_mean": float(np.mean(sel)), "excluded_mean": float(np.mean(exc)),
            "standardized_mean_difference": observed[name],
            "restricted_block_p": within, "familywise_maxT_p_4": family,
            "mixed_blocks": len(mixed),
        })

    args.output.mkdir(parents=True, exist_ok=False)
    balance_path = args.output / "selection_balance.tsv"
    with balance_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    block_path = args.output / "block_support.tsv"
    with block_path.open("w", newline="", encoding="utf-8") as stream:
        fields = ["block", "n_i2_pairs", "n_common_wl", "n_without_common_wl"]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for block, ix in sorted(blocks.items()):
            n_sel = int(selected[ix].sum())
            writer.writerow({"block": block, "n_i2_pairs": len(ix),
                             "n_common_wl": n_sel,
                             "n_without_common_wl": len(ix) - n_sel})
    summary = {
        "schema": "astrafek-icl-wl-overlap-selection-audit-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "hsc_query_submitted": False,
        "parent_i2_satellite_supported_n": len(rows),
        "common_wl_supported_n": int(selected.sum()),
        "no_common_wl_support_n": int((~selected).sum()),
        "stage_counts": stage_counts,
        "selection_rule": "random_subtracted_wl_common_support from pinned eligibility audit",
        "spatial_blocks": len(blocks), "mixed_selection_blocks": len(mixed),
        "largest_block_n": max(map(len, blocks.values())),
        "covariates": list(variables),
        "permutation": {
            "method": "permute common-WL-support labels within saved connected 5-degree blocks, preserving each block's selected count",
            "replicates": args.replicates, "seed": SEED,
            "familywise_method": "max absolute standardized mean difference across four prespecified physical/covariate balance variables; target shape-pair count is reported separately as a structural support gate",
        },
        "support_mechanism": pair_count_mechanism,
        "selection_ids": sorted(common_ids, key=int),
        "excluded_ids": sorted(set(pair_by_id) - common_ids, key=int),
        "input_sha256": {str(path): hashes[key] for key, path in paths.items()},
        "output_sha256": {balance_path.name: sha256(balance_path),
                          block_path.name: sha256(block_path)},
        "limitations": [
            "selection balance within the selected 68-lens satellite-supported subset does not establish missing-at-random",
            "only mixed 5-degree blocks contribute to the restricted permutation; the block partition is approximate",
            "small selected sample limits power and cannot support a population correction",
            "this audit does not recover missing shape rows or add HSC data",
        ],
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "stage_counts": stage_counts, "mixed_blocks": len(mixed),
                      "familywise_results": output_rows}, indent=2))


if __name__ == "__main__":
    main()
