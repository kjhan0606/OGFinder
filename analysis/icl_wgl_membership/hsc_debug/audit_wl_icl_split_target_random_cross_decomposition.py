#!/usr/bin/env python3
"""Separate target and random contributions to the i-band split cross profile.

This deterministic downstream audit uses the saved source-region stack
summary. It verifies that the regional high-minus-low random-subtracted
profile equals the target difference minus the random-aperture difference,
then writes all endpoint terms with provenance. It reads no shape rows or
images and makes no HSC request.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCE_REGION = DATA / "wl_icl_split_source_region_cross_null_2026-11-06_01/summary.json"
DEFAULT_OUTPUT = DATA / "wl_icl_split_target_random_cross_decomposition_2026-11-06_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def delta(high: dict, low: dict, key: str) -> float:
    return float(high[key]) - float(low[key])


def run(output: Path = DEFAULT_OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite target/random decomposition: {output}")
    source = json.loads(SOURCE_REGION.read_text(encoding="utf-8"))
    if (source.get("schema") != "astrafek-wl-icl-split-source-region-cross-null-v1"
            or source.get("status") != "complete_descriptive_diagnostic"
            or source.get("science_measurement") is not False):
        raise ValueError("unexpected source-region diagnostic schema or status")

    endpoints = source["estimator_counts_and_region_profiles"]
    regional_saved = source["regional_random_subtracted_high_minus_low"]
    output_regions = {}
    for region, saved_rows in regional_saved.items():
        term_rows = {}
        for endpoint in ("target", "random"):
            low = endpoints[endpoint]["by_group_region"]["low_icl"][region]["profile"]
            high = endpoints[endpoint]["by_group_region"]["high_icl"][region]["profile"]
            if len(low) != len(high) or len(low) != len(saved_rows):
                raise ValueError(f"profile length mismatch in region {region}/{endpoint}")
            term_rows[endpoint] = (low, high)

        rows = []
        for index, saved in enumerate(saved_rows):
            low_target, high_target = term_rows["target"][0][index], term_rows["target"][1][index]
            low_random, high_random = term_rows["random"][0][index], term_rows["random"][1][index]
            row = {
                "r_min_mpc": float(saved["r_min_mpc"]),
                "r_max_mpc": float(saved["r_max_mpc"]),
                "target_high_minus_low_tangential_msun_mpc2": delta(
                    high_target, low_target, "tangential_msun_mpc2"),
                "random_high_minus_low_tangential_msun_mpc2": delta(
                    high_random, low_random, "tangential_msun_mpc2"),
                "target_high_minus_low_cross_msun_mpc2": delta(
                    high_target, low_target, "cross_msun_mpc2"),
                "random_high_minus_low_cross_msun_mpc2": delta(
                    high_random, low_random, "cross_msun_mpc2"),
            }
            for component, signal in (("tangential", "tangential_msun_mpc2"),
                                      ("cross", "cross_msun_mpc2")):
                target_delta = row[f"target_high_minus_low_{component}_msun_mpc2"]
                random_delta = row[f"random_high_minus_low_{component}_msun_mpc2"]
                reconstructed = target_delta - random_delta
                saved_key = f"high_minus_low_{component}_random_subtracted_msun_mpc2"
                if abs(reconstructed - float(saved[saved_key])) > max(1e-5, abs(reconstructed) * 1e-12):
                    raise ValueError(f"endpoint decomposition mismatch for {region}/{index}/{component}")
                row[f"random_subtracted_high_minus_low_{component}_msun_mpc2"] = reconstructed
            row["target_cross_used_pairs_low_high"] = [
                int(low_target["used_pairs"]), int(high_target["used_pairs"]),
            ]
            row["random_cross_used_pairs_low_high"] = [
                int(low_random["used_pairs"]), int(high_random["used_pairs"]),
            ]
            rows.append(row)
        output_regions[str(region)] = rows

    result = {
        "schema": "astrafek-wl-icl-split-target-random-cross-decomposition-v1",
        "status": "complete_descriptive_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "shape_rows_read": False,
        "images_read": False,
        "method": "exact endpoint arithmetic on the region-specific HSC Y3 profiles; high-minus-low random-subtracted signal = target high-minus-low minus random-center high-minus-low",
        "units": "M_sun Mpc^-2",
        "regions": output_regions,
        "interpretation": "Endpoint decomposition identifies which saved point estimates contribute to the cross profile; it provides no covariance, uncertainty, or physical/instrumental cause.",
        "provenance": {
            "source_region_summary": str(SOURCE_REGION),
            "source_region_summary_sha256": sha256(SOURCE_REGION),
            "source_region_runner_sha256": source["analysis_script_sha256"],
        },
    }
    output.mkdir(parents=True)
    report = output / "summary.json"
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    report.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                 allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(report)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "cross_target_delta_1e12": {
            region: [round(row["target_high_minus_low_cross_msun_mpc2"] / 1e12, 3)
                     for row in rows]
            for region, rows in output_regions.items()
        },
        "cross_random_delta_1e12": {
            region: [round(row["random_high_minus_low_cross_msun_mpc2"] / 1e12, 3)
                     for row in rows]
            for region, rows in output_regions.items()
        },
        "cross_subtracted_1e12": {
            region: [round(row["random_subtracted_high_minus_low_cross_msun_mpc2"] / 1e12, 3)
                     for row in rows]
            for region, rows in output_regions.items()
        },
        "summary": str(report),
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
