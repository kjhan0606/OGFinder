#!/usr/bin/env python3
"""Compare verified 2-D ICL stacks across target/control Kron-mask radii."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_RUNS = {
    2.5: DEFAULT_SAMPLE / "icl_2d_matched_stack_2026-10-19_02",
    3.5: DEFAULT_SAMPLE / "icl_2d_matched_stack_2026-10-19_k3p5",
    5.0: DEFAULT_SAMPLE / "icl_2d_matched_stack_2026-10-19_k5p0",
}
DEFAULT_OUTPUT = DEFAULT_SAMPLE / "icl_2d_mask_sensitivity_assessment_2026-10-19_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_run(work: Path, expected_factor: float) -> tuple[dict, list[dict]]:
    summary_path = work / "summary.json"
    verification_path = work / "independent_verification.json"
    evaluation_path = work / "evaluation.json"
    profile_path = work / "broad_radial_profiles.tsv"
    for path in (summary_path, verification_path, evaluation_path, profile_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing completed sensitivity input: {path}")
    summary, verification, evaluation = map(read_json, (summary_path, verification_path, evaluation_path))
    factor = float(summary.get("mask_factor", 2.5))
    if factor != expected_factor:
        raise ValueError(f"{work}: expected mask factor {expected_factor}, found {factor}")
    if summary.get("status") != "complete_exploratory" or summary.get("science_measurement") is not False:
        raise ValueError(f"{work}: stack is incomplete or science gate is not false")
    if verification.get("status") != "verified" or verification.get("science_measurement") is not False:
        raise ValueError(f"{work}: independent verification did not pass")
    if evaluation.get("status") != "completed_exploratory" or evaluation.get("science_measurement") is not False:
        raise ValueError(f"{work}: radial evaluation is incomplete")
    if verification.get("n_eligible_pairs") != summary.get("n_pairs_stacked"):
        raise ValueError(f"{work}: pair counts disagree between summary and verifier")
    rows = read_tsv(profile_path)
    if len(rows) != 12:
        raise ValueError(f"{work}: expected six radial bins in each of two filters")
    for row in rows:
        if row.get("mean_delta_njy_arcsec2") in (None, ""):
            raise ValueError(f"{work}: missing radial delta")
    metadata = {
        "mask_factor": factor, "work": str(work.resolve()),
        "summary_sha256": sha256(summary_path),
        "verification_sha256": sha256(verification_path),
        "evaluation_sha256": sha256(evaluation_path),
        "profile_sha256": sha256(profile_path),
        "n_pairs": int(summary["n_pairs_stacked"]),
        "n_by_filter": summary["n_by_filter"],
    }
    return metadata, rows


def assemble_comparison(runs: dict[float, list[dict]]) -> list[dict]:
    factors = (2.5, 3.5, 5.0)
    lookups = {}
    for factor in factors:
        lookup = {}
        for row in runs[factor]:
            key = (row["filter"], int(row["r_min_kpc"]), int(row["r_max_kpc"]))
            if key in lookup:
                raise ValueError(f"duplicate radial bin {key} for Kron {factor}")
            lookup[key] = row
        if len(lookup) != 12:
            raise ValueError(f"Kron {factor}: expected 12 unique filter/radius rows")
        lookups[factor] = lookup
    keys = set(lookups[2.5])
    if any(set(lookups[f]) != keys for f in factors[1:]):
        raise ValueError("filter/radial-bin grids differ between mask factors")
    output = []
    for key in sorted(keys):
        base = lookups[2.5][key]
        row = {"filter": key[0], "r_min_kpc": key[1], "r_max_kpc": key[2],
               "baseline_k2p5_n_pairs": int(base["n_pairs"]),
               "baseline_k2p5_median_coverage": float(base["median_pair_common_coverage"]),
               "baseline_k2p5_delta_njy_arcsec2": float(base["mean_delta_njy_arcsec2"]),
               "baseline_k2p5_pair_sem_njy_arcsec2": float(base["pair_sem_delta_njy_arcsec2"])}
        for factor in factors[1:]:
            item = lookups[factor][key]
            suffix = str(factor).replace(".", "p")
            delta = float(item["mean_delta_njy_arcsec2"])
            row[f"k{suffix}_n_pairs"] = int(item["n_pairs"])
            row[f"k{suffix}_median_coverage"] = float(item["median_pair_common_coverage"])
            row[f"k{suffix}_delta_njy_arcsec2"] = delta
            row[f"k{suffix}_pair_sem_njy_arcsec2"] = float(item["pair_sem_delta_njy_arcsec2"])
            row[f"k{suffix}_delta_minus_k2p5_njy_arcsec2"] = delta - row["baseline_k2p5_delta_njy_arcsec2"]
        output.append(row)
    return output


def plot_comparison(rows: list[dict], factor_metadata: list[dict], output: Path) -> dict:
    figures = {}
    factors = ((2.5, "k2p5"), (3.5, "k3p5"), (5.0, "k5p0"))
    styles = (("#1b6ca8", "o"), ("#c05a27", "s"), ("#4c8c45", "^") )
    for band in ("i", "i2"):
        selected = sorted((r for r in rows if r["filter"] == band), key=lambda r: r["r_min_kpc"])
        x = np.asarray([(r["r_min_kpc"] + r["r_max_kpc"]) / 2 for r in selected], dtype=float)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
        for factor, stem, (color, marker) in zip((f[0] for f in factors), (f[1] for f in factors), styles):
            y = np.asarray([r["baseline_k2p5_delta_njy_arcsec2"] if factor == 2.5 else r[f"{stem}_delta_njy_arcsec2"] for r in selected])
            e = np.asarray([r["baseline_k2p5_pair_sem_njy_arcsec2"] if factor == 2.5 else r[f"{stem}_pair_sem_njy_arcsec2"] for r in selected])
            n = np.asarray([r["baseline_k2p5_n_pairs"] if factor == 2.5 else r[f"{stem}_n_pairs"] for r in selected])
            coverage = np.asarray([r["baseline_k2p5_median_coverage"] if factor == 2.5 else r[f"{stem}_median_coverage"] for r in selected])
            axes[0].errorbar(x, y, yerr=e, color=color, marker=marker, capsize=2, label=f"Kron {factor:g}")
            axes[1].plot(x, n, color=color, marker=marker, label=f"Kron {factor:g}")
            axes[2].plot(x, coverage, color=color, marker=marker, label=f"Kron {factor:g}")
        axes[0].axhline(0, color="black", linewidth=.8)
        axes[0].set_ylabel("target − control (nJy/arcsec²)")
        axes[0].set_title("Equal-pair radial mean ± pair SEM")
        axes[1].set_ylabel("pairs with annular support")
        axes[1].set_title("Radius-dependent sample size")
        axes[2].set_ylabel("median common-area coverage")
        axes[2].set_title("Among pairs contributing")
        for ax in axes:
            ax.set_xlabel("proper radius (kpc)")
            ax.grid(alpha=.2)
        axes[0].legend()
        path = output / f"mask_sensitivity_{band}.png"
        fig.suptitle(f"HSC-{band} 2-D stack mask-radius sensitivity (exploratory)")
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures[band] = str(path)
    return figures


def run(work_paths: dict[float, Path], output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    metadata, runs = [], {}
    for factor in (2.5, 3.5, 5.0):
        item, rows = read_run(work_paths[factor], factor)
        metadata.append(item)
        runs[factor] = rows
    rows = assemble_comparison(runs)
    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "mask_sensitivity_comparison.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)
    figures = plot_comparison(rows, metadata, output)
    result = {
        "schema": "astrafek-hsc-icl-2d-mask-sensitivity-v1",
        "status": "completed_exploratory", "mask_factors": [2.5, 3.5, 5.0],
        "runs": metadata, "n_comparison_rows": len(rows),
        "comparison_table": str(table_path), "comparison_table_sha256": sha256(table_path),
        "figures": figures, "science_measurement": False,
        "interpretation_limit": "Pair SEM is descriptive only; changing masks alter support and the contributing sample with radius. This comparison is not a detection test or a full covariance-aware uncertainty analysis.",
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-k2p5", type=Path, default=DEFAULT_RUNS[2.5])
    parser.add_argument("--work-k3p5", type=Path, default=DEFAULT_RUNS[3.5])
    parser.add_argument("--work-k5p0", type=Path, default=DEFAULT_RUNS[5.0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    run({2.5: args.work_k2p5.resolve(), 3.5: args.work_k3p5.resolve(),
         5.0: args.work_k5p0.resolve()}, args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
