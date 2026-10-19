#!/usr/bin/env python3
"""Test sensitivity of the fixed ICL–satellite sample to separate pixel support."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_satellite_icl_spatial_blocks import coordinate_labels, load_verified_rows
from verify_icl_independent_support import verify as verify_independent_support

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_ICL = ROOT / "icl_annulus_support_2026-10-18_01"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"
DEFAULT_WORKS = {
    2.5: ROOT / "icl_2d_independent_support_k2p5_2026-10-19_02",
    3.5: ROOT / "icl_2d_independent_support_k3p5_2026-10-19_02",
    5.0: ROOT / "icl_2d_independent_support_k5p0_2026-10-19_02",
}
DEFAULT_OUTPUT = ROOT / "icl_satellite_independent_support_association_2026-10-20_01"
PAIR_FILE = "paired_independent_support_profiles.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def independent_annulus(profile: list[dict], lo_kpc: int = 200,
                        hi_kpc: int = 300) -> dict:
    """Combine 25 kpc annuli with pixel-count weighting on each frame separately."""
    selected = [row for row in profile
                if lo_kpc <= float(row["r_min_kpc"])
                and float(row["r_max_kpc"]) <= hi_kpc]
    expected_edges = list(range(lo_kpc, hi_kpc, 25))
    if ([int(float(row["r_min_kpc"])) for row in selected] != expected_edges
            or any(int(float(row["r_max_kpc"])) != int(float(row["r_min_kpc"])) + 25
                   for row in selected)):
        raise ValueError("independent profile does not contain the exact 25 kpc 200–300 bins")
    nt = sum(int(row["n_target_pixels"]) for row in selected)
    nc = sum(int(row["n_control_pixels"]) for row in selected)
    geom = sum(int(row["geometric_pixels"]) for row in selected)
    if nt <= 0 or nc <= 0 or geom <= 0:
        raise ValueError("independent target/control support is empty in 200–300 kpc")
    target = sum(float(row["target_mean_njy_arcsec2"]) * int(row["n_target_pixels"])
                 for row in selected) / nt
    control = sum(float(row["control_mean_njy_arcsec2"]) * int(row["n_control_pixels"])
                  for row in selected) / nc
    if not all(math.isfinite(x) for x in (target, control)):
        raise ValueError("non-finite independent annulus mean")
    return {
        "target_mean_njy_arcsec2": target,
        "control_mean_njy_arcsec2": control,
        "delta_njy_arcsec2": target - control,
        "n_target_pixels": nt, "n_control_pixels": nc,
        "geometric_pixels": geom,
        "target_coverage": nt / geom, "control_coverage": nc / geom,
    }


def partial_pearson(rows: list[dict], x_key: str, y_key: str) -> float:
    x = np.asarray([float(row[x_key]) for row in rows])
    y = np.asarray([float(row[y_key]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    if len(rows) < 5 or not all(np.all(np.isfinite(a)) for a in (x, y, z, richness)):
        raise ValueError("partial correlation needs at least five finite rows")
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("redshift/richness design is rank deficient")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if denominator <= 0:
        raise ValueError("partial correlation is undefined")
    return float(rx @ ry / denominator)


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = (start + stop - 1) / 2.0
        start = stop
    return ranks


def spearman_rho(rows: list[dict], x_key: str, y_key: str) -> float:
    x = np.asarray([float(row[x_key]) for row in rows])
    y = np.asarray([float(row[y_key]) for row in rows])
    rx, ry = average_ranks(x), average_ranks(y)
    rx -= rx.mean(); ry -= ry.mean()
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if denominator <= 0:
        raise ValueError("Spearman correlation is undefined")
    return float(rx @ ry / denominator)


def paired_block_bootstrap(rows: list[dict], labels: np.ndarray, replicas: int,
                           seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    if len(labels) != len(rows):
        raise ValueError("one connected block label is required per row")
    groups = np.unique(labels)
    if len(groups) < 6:
        raise ValueError("at least six connected spatial blocks are required")
    indices = {group: np.flatnonzero(labels == group) for group in groups}
    baseline = partial_pearson(rows, "icl_delta_200_300", "satellite_delta_50_75")
    alternate = partial_pearson(rows, "icl_independent_delta_200_300", "satellite_delta_50_75")
    rng = np.random.default_rng(seed)
    base_draws, alternate_draws, paired_differences = [], [], []
    invalid = 0
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        sample_indices = np.concatenate([indices[group] for group in chosen])
        sample = [rows[index] for index in sample_indices]
        try:
            base_draw = partial_pearson(sample, "icl_delta_200_300", "satellite_delta_50_75")
            alt_draw = partial_pearson(sample, "icl_independent_delta_200_300", "satellite_delta_50_75")
        except ValueError:
            invalid += 1
            continue
        base_draws.append(base_draw)
        alternate_draws.append(alt_draw)
        paired_differences.append(alt_draw - base_draw)
    if not paired_differences:
        raise ValueError("no valid paired block-bootstrap replicates")
    q = (2.5, 97.5)
    sizes = [len(indices[group]) for group in groups]
    return {
        "n_blocks": int(len(groups)), "largest_block_pairs": int(max(sizes)),
        "baseline_partial_r": baseline, "independent_support_partial_r": alternate,
        "partial_r_difference": alternate - baseline,
        "baseline_block_ci95": np.percentile(base_draws, q).tolist(),
        "independent_support_block_ci95": np.percentile(alternate_draws, q).tolist(),
        "paired_difference_block_ci95": np.percentile(paired_differences, q).tolist(),
        "baseline_spearman_rho": spearman_rho(rows, "icl_delta_200_300", "satellite_delta_50_75"),
        "independent_support_spearman_rho": spearman_rho(
            rows, "icl_independent_delta_200_300", "satellite_delta_50_75"),
        "replicas_requested": replicas, "replicas_valid": len(paired_differences),
        "replicas_invalid": invalid, "seed": seed,
    }


def load_independent_pairs(work: Path) -> tuple[dict, dict]:
    verification = verify_independent_support(work)
    if verification.get("status") != "verified" or verification.get("science_measurement") is not False:
        raise ValueError(f"independent-support run is not verified: {work}")
    summary = json.loads((work / "summary.json").read_text(encoding="utf-8"))
    pair_path = work / PAIR_FILE
    rows = read_tsv(pair_path)
    if sha256(pair_path) != summary["outputs"]["pair_profiles_sha256"]:
        raise ValueError(f"independent pair profiles fail saved checksum: {work}")
    result = {}
    for row in rows:
        key = (row["target_id"], row["control_id"], row["filter"])
        if key in result:
            raise ValueError(f"duplicate independent-support pair {key}")
        result[key] = {"profile": json.loads(row["profile_json"]),
                       "z": float(row["z"]), "ra_deg": float(row["ra_deg"]),
                       "dec_deg": float(row["dec_deg"])}
    return summary, result


def combine_rows(base_rows: list[dict], profiles: dict) -> tuple[list[dict], list[dict]]:
    joined, pair_records = [], []
    for base in base_rows:
        key = (str(base["cluster_id"]), str(base["control_id"]), str(base["filter"]))
        pair_input = profiles.get(key)
        if pair_input is None:
            raise ValueError(f"fixed association row missing independent profile: {key}")
        if not math.isclose(float(base["redshift"]), pair_input["z"], rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"redshift differs between fixed sample and independent profile: {key}")
        profile = pair_input["profile"]
        annulus = independent_annulus(profile)
        if not math.isfinite(float(base["richness"])):
            raise ValueError(f"invalid association metadata for {key}")
        row = dict(base)
        row["icl_independent_delta_200_300"] = annulus["delta_njy_arcsec2"]
        joined.append(row)
        pair_records.append({
            "cluster_id": key[0], "control_id": key[1], "filter": key[2],
            "icl_common_delta_200_300": float(base["icl_delta_200_300"]),
            "icl_independent_delta_200_300": annulus["delta_njy_arcsec2"],
            "satellite_delta_50_75": float(base["satellite_delta_50_75"]),
            "redshift": float(base["redshift"]), "richness": float(base["richness"]),
            "target_coverage_independent": annulus["target_coverage"],
            "control_coverage_independent": annulus["control_coverage"],
            "n_target_pixels_independent": annulus["n_target_pixels"],
            "n_control_pixels_independent": annulus["n_control_pixels"],
            "science_measurement": False,
        })
    if len(joined) != len(base_rows):
        raise ValueError("independent-support join changed the fixed association sample")
    return joined, pair_records


def analyze(association_dir: Path, icl_dir: Path, extraction_summary: Path,
            independent_works: dict[float, Path], output_dir: Path,
            replicas: int = 5000, seed: int = 20261021) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 connected-block replicas are required")
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_dir}")
    base_rows, baseline_provenance = load_verified_rows(
        association_dir, icl_dir, extraction_summary)
    expected_n = len(base_rows)
    if expected_n != 94:
        raise ValueError(f"expected fixed 94-row association sample, found {expected_n}")
    if set(independent_works) != {2.5, 3.5, 5.0}:
        raise ValueError("exactly K=2.5, 3.5, and 5.0 inputs are required")

    metric_rows, pair_rows, input_hashes, code_hashes = [], [], {}, {}
    input_paths = [Path(association_dir) / "summary.json",
                   Path(association_dir) / "matched_cluster_values.tsv",
                   Path(icl_dir) / "summary.json",
                   Path(icl_dir) / "annulus_pair_profiles.tsv", Path(extraction_summary)]
    input_hashes.update({str(Path(path).resolve()): sha256(Path(path).resolve())
                         for path in input_paths})
    source_files = [Path(__file__).resolve(), Path(coordinate_labels.__code__.co_filename).resolve(),
                    Path(verify_independent_support.__code__.co_filename).resolve()]
    code_hashes.update({str(path): sha256(path) for path in source_files})
    selection_defs = (
        ("all_supported_pairs", lambda row: True),
        ("center_eligible", lambda row: str(row["center_eligible"]).lower() == "true"),
    )
    for index, mask in enumerate((2.5, 3.5, 5.0)):
        work = Path(independent_works[mask]).resolve()
        independent_summary, profiles = load_independent_pairs(work)
        if float(independent_summary["mask_factor"]) != mask:
            raise ValueError(f"unexpected mask factor in {work}")
        input_hashes[str(work / "summary.json")] = sha256(work / "summary.json")
        input_hashes[str(work / "independent_verification.json")] = sha256(work / "independent_verification.json")
        input_hashes[str(work / PAIR_FILE)] = sha256(work / PAIR_FILE)
        code_hashes.update({path: digest for path, digest in independent_summary["scripts_sha256"].items()})
        rows, these_pairs = combine_rows(base_rows, profiles)
        pair_rows.extend({"mask_factor": mask, **row} for row in these_pairs)
        for band_index, band in enumerate(("i2", "i")):
            band_rows = [row for row in rows if row["filter"] == band]
            for selection_index, (selection, predicate) in enumerate(selection_defs):
                selected = [row for row in band_rows if predicate(row)]
                if len(selected) < 8:
                    raise ValueError(f"too few rows for {mask}/{band}/{selection}: {len(selected)}")
                labels = coordinate_labels(selected)
                result = paired_block_bootstrap(
                    selected, labels, replicas,
                    seed + index * 100 + band_index * 10 + selection_index)
                result.update({"mask_factor": mask, "filter": band,
                               "selection": selection, "n_pairs": len(selected),
                               "science_measurement": False})
                metric_rows.append(result)
                print(f"Independent-support association K={mask:g} {band}/{selection}: "
                      f"pairs={len(selected)} blocks={result['n_blocks']} "
                      f"delta_r={result['partial_r_difference']:+.4f} "
                      f"valid_draws={result['replicas_valid']}/{replicas}", flush=True)

    output_dir.mkdir(parents=True, exist_ok=False)
    pair_path = output_dir / "fixed_sample_pair_values.tsv"
    metric_path = output_dir / "association_sensitivity.tsv"
    write_tsv(pair_path, pair_rows)
    write_tsv(metric_path, metric_rows)
    summary = {
        "schema": "astrafek-icl-satellite-independent-support-association-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "sample": "fixed 94 cluster/filter rows from the pre-existing K2.5 per-annulus support association; no re-selection",
            "baseline": "existing common-support K2.5 200–300 kpc ICL delta",
            "alternate": "200–300 kpc ICL target/control means recombined from four 25 kpc bins with each frame's own valid pixels, weighted by valid pixel counts",
            "association": "partial Pearson residualized on redshift and richness; Spearman is an unadjusted companion",
            "uncertainty": "paired connected target/control 5-degree block bootstrap; each draw uses the same blocks for baseline and alternate estimators",
            "replicas": replicas, "seed": seed,
        },
        "n_fixed_rows": expected_n,
        "n_by_filter": {band: sum(row["filter"] == band for row in base_rows) for band in ("i", "i2")},
        "baseline_provenance": baseline_provenance,
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__,
                        "platform": platform.platform()},
        "independent_support_runs": {str(mask): {
            "work": str(Path(independent_works[mask]).resolve()),
            "summary_sha256": sha256(Path(independent_works[mask]) / "summary.json"),
            "verification_sha256": sha256(Path(independent_works[mask]) / "independent_verification.json"),
        } for mask in (2.5, 3.5, 5.0)},
        "input_sha256": input_hashes, "code_sha256": code_hashes,
        "outputs": {
            pair_path.name: sha256(pair_path), metric_path.name: sha256(metric_path),
        },
        "limitations": [
            "The 94-row association sample was selected by the common-support K2.5 estimator; this is a fixed-sample estimator sensitivity, not a population reanalysis.",
            "The alternate estimator uses different valid target and control pixels; azimuthal symmetry and matched populations are assumptions.",
            "The paired block interval describes the change between estimators within this selected sample; it does not include mask selection, extended PSF, sky transfer, extinction, or full coadd covariance uncertainty.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "n_fixed_rows": expected_n,
                      "science_measurement": False,
                      "n_metric_rows": len(metric_rows)}, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl-support", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--k2p5-work", type=Path, default=DEFAULT_WORKS[2.5])
    parser.add_argument("--k3p5-work", type=Path, default=DEFAULT_WORKS[3.5])
    parser.add_argument("--k5p0-work", type=Path, default=DEFAULT_WORKS[5.0])
    parser.add_argument("--work", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    independent_works = {2.5: args.k2p5_work, 3.5: args.k3p5_work, 5.0: args.k5p0_work}
    analyze(args.association, args.icl_support, args.extraction_summary,
            independent_works, args.work, args.replicas, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
