#!/usr/bin/env python3
"""Estimate spatial-block uncertainty for verified matched-control ICL stacks.

This resamples connected 5-degree sky cells shared by target or control
centers. It reports the pair bootstrap beside the spatial-block bootstrap and
keeps both exploratory; neither is a full image-level covariance model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_RUNS = {
    2.5: ROOT / "icl_2d_matched_stack_2026-10-19_02",
    3.5: ROOT / "icl_2d_matched_stack_2026-10-19_k3p5",
    5.0: ROOT / "icl_2d_matched_stack_2026-10-19_k5p0",
}
DEFAULT_OUTPUT = ROOT / "icl_2d_spatial_block_bootstrap_2026-10-19_01"
RADII = ((50, 100), (100, 200), (200, 300), (300, 400), (400, 500))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def spatial_cell(ra: float, dec: float, width: float) -> str:
    if not np.isfinite(ra) or not 0 <= ra < 360:
        raise ValueError(f"invalid RA: {ra}")
    if not np.isfinite(dec) or not -90 <= dec <= 90:
        raise ValueError(f"invalid Dec: {dec}")
    if not np.isfinite(width) or width <= 0:
        raise ValueError("cell width must be positive")
    return f"{int(ra // width)}_{int((dec + 90.0) // width)}"


def connected_pair_blocks(targets: list[dict], controls: list[dict], width: float = 5.0) -> np.ndarray:
    """Connect pairs sharing any target/control sky cell, regardless of role."""
    if len(targets) != len(controls) or not targets:
        raise ValueError("target/control coordinate lists must have equal non-zero lengths")
    parent = list(range(len(targets)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        root_l, root_r = find(left), find(right)
        if root_l != root_r:
            parent[max(root_l, root_r)] = min(root_l, root_r)

    owner = {}
    for index, (target, control) in enumerate(zip(targets, controls)):
        cells = {
            spatial_cell(float(target["ra_deg"]), float(target["dec_deg"]), width),
            spatial_cell(float(control["ra_deg"]), float(control["dec_deg"]), width),
        }
        for cell in cells:
            if cell in owner:
                union(index, owner[cell])
            else:
                owner[cell] = index
    roots = {root: label for label, root in enumerate(sorted({find(i) for i in range(len(parent))}))}
    return np.asarray([f"sky_block_{roots[find(i)]}" for i in range(len(parent))])


def profile_matrix(rows: list[dict[str, str]], pair_keys: list[tuple[str, str, str]],
                   band: str) -> np.ndarray:
    key_index = {key: index for index, key in enumerate(pair_keys)}
    matrix = np.full((len(pair_keys), len(RADII)), np.nan, dtype=float)
    seen = set()
    for row in rows:
        if row["filter"] != band:
            continue
        key = (row["target_id"], row["control_id"], row["filter"])
        if key not in key_index:
            raise ValueError(f"unexpected pair profile row: {key}")
        radius = (int(row["r_min_kpc"]), int(row["r_max_kpc"]))
        if radius not in RADII:
            continue
        slot = (key, radius)
        if slot in seen:
            raise ValueError(f"duplicate pair/radius profile: {slot}")
        seen.add(slot)
        value = row.get("delta_njy_arcsec2", "")
        matrix[key_index[key], RADII.index(radius)] = float(value) if value not in ("", "None") else np.nan
    expected = sum(key[2] == band for key in pair_keys) * len(RADII)
    if len(seen) != expected:
        raise ValueError(f"{band}: expected {expected} outer-annulus pair rows, found {len(seen)}")
    return matrix


def summarize_resamples(values: np.ndarray, labels: np.ndarray, *, replicas: int, seed: int,
                        cell_width_deg: float = 5.0) -> dict:
    """Return pair/block bootstrap distributions for a matrix of pair profiles."""
    if values.ndim != 2 or values.shape[0] != len(labels):
        raise ValueError("profile matrix and block labels have incompatible shapes")
    if replicas < 100:
        raise ValueError("at least 100 bootstrap replicas are required")
    n_pairs, n_radii = values.shape
    usable = np.isfinite(values)
    counts = usable.sum(axis=0)
    if np.any(counts == 0):
        raise ValueError("an annulus has no supported pairs")
    point = np.nansum(values, axis=0) / counts
    pair_sd = np.nanstd(values, axis=0, ddof=1)
    pair_sem = pair_sd / np.sqrt(counts)
    groups = np.unique(labels)
    if len(groups) < 6:
        raise ValueError(f"only {len(groups)} connected spatial blocks; at least 6 required")

    rng = np.random.default_rng(seed)
    pair_draws = rng.integers(0, n_pairs, size=(replicas, n_pairs))
    pair_boot = np.full((replicas, n_radii), np.nan)
    for radial_index in range(n_radii):
        sampled = values[pair_draws, radial_index]
        valid = np.isfinite(sampled)
        den = valid.sum(axis=1)
        ok = den > 0
        pair_boot[ok, radial_index] = np.nansum(sampled[ok], axis=1) / den[ok]

    block_sum = np.zeros((len(groups), n_radii), dtype=float)
    block_n = np.zeros((len(groups), n_radii), dtype=float)
    for index, group in enumerate(groups):
        rows = labels == group
        block_sum[index] = np.nansum(values[rows], axis=0)
        block_n[index] = np.isfinite(values[rows]).sum(axis=0)
    block_draws = rng.multinomial(len(groups), np.full(len(groups), 1 / len(groups)), size=replicas)
    numerator = block_draws @ block_sum
    denominator = block_draws @ block_n
    block_boot = np.full((replicas, n_radii), np.nan)
    np.divide(numerator, denominator, out=block_boot, where=denominator > 0)

    def intervals(draws):
        result = []
        for index in range(n_radii):
            valid = draws[:, index][np.isfinite(draws[:, index])]
            result.append({
                "valid_replicas": int(valid.size),
                "ci95": np.percentile(valid, [2.5, 97.5]).tolist() if valid.size else None,
                "bootstrap_sd": float(valid.std(ddof=1)) if valid.size > 1 else None,
            })
        complete = draws[np.isfinite(draws).all(axis=1)]
        covariance = (np.cov(complete, rowvar=False, ddof=1).tolist()
                      if len(complete) > 1 else None)
        return result, covariance, int(len(complete))

    pair_intervals, pair_covariance, pair_complete = intervals(pair_boot)
    block_intervals, block_covariance, block_complete = intervals(block_boot)
    return {
        "n_pairs_by_annulus": counts.astype(int).tolist(),
        "point_mean_njy_arcsec2": point.tolist(),
        "pair_sem_njy_arcsec2": pair_sem.tolist(),
        "n_connected_blocks": int(len(groups)),
        "largest_block_pairs": int(max(np.count_nonzero(labels == group) for group in groups)),
        "replicas_requested": int(replicas),
        "pair_bootstrap": {"intervals": pair_intervals, "complete_profile_replicas": pair_complete,
                           "covariance_njy2_arcsec4": pair_covariance},
        "spatial_block_bootstrap": {"method": f"connected {cell_width_deg:g}-degree RA/Dec cells; target/control roles share cells",
                                    "intervals": block_intervals,
                                    "complete_profile_replicas": block_complete,
                                    "covariance_njy2_arcsec4": block_covariance},
    }


def load_run(work: Path, factor: float, target_manifest: dict, control_manifest: dict,
             expected_pairs: list[tuple[str, str, str]]):
    summary_path = work / "summary.json"
    verification_path = work / "independent_verification.json"
    evaluation_path = work / "evaluation.json"
    profiles_path = work / "broad_radial_pair_profiles.tsv"
    summary, verification, evaluation = map(read_json, (summary_path, verification_path, evaluation_path))
    if summary.get("status") != "complete_exploratory" or summary.get("science_measurement") is not False:
        raise ValueError(f"{work}: stack incomplete or science status changed")
    if float(summary.get("mask_factor", 2.5)) != factor:
        raise ValueError(f"{work}: expected Kron mask factor {factor}")
    for path in (ROOT / "sample_manifest.json", ROOT / "control_download_01/download_manifest.json"):
        if summary.get("input_sha256", {}).get(str(path)) != sha256(path):
            raise ValueError(f"{work}: target/control coordinate manifest changed: {path}")
    for band in ("i", "i2"):
        stack_path = Path(summary["outputs"]["maps"][band]["file"])
        with fits.open(stack_path, memmap=False, checksum=True) as hdus:
            header_factor = hdus[0].header.get("MASKFAC")
            if header_factor is None:
                # The original verified K=2.5 run predates MASKFAC in FITS headers.
                mask_text = summary.get("stack_definition", {}).get("mask", "")
                if factor != 2.5 or "MASK_K2P5" not in mask_text:
                    raise ValueError(f"{work}: missing explicit mask-factor metadata")
            elif float(header_factor) != factor:
                raise ValueError(f"{work}: {band} FITS header has the wrong mask factor")
            if hdus[0].header.get("SCIVAL") is not False:
                raise ValueError(f"{work}: {band} FITS has an unexpected science gate")
    if verification.get("status") != "verified" or verification.get("science_measurement") is not False:
        raise ValueError(f"{work}: independent verification is not passed")
    if evaluation.get("status") != "completed_exploratory" or evaluation.get("science_measurement") is not False:
        raise ValueError(f"{work}: radial evaluation is not complete")
    rows = read_tsv(profiles_path)
    pair_keys = [(row["target_id"], row["control_id"], row["filter"]) for row in rows
                 if (row["r_min_kpc"], row["r_max_kpc"]) == ("50", "100")]
    if len(pair_keys) != summary.get("n_pairs_stacked") or len(set(pair_keys)) != len(pair_keys):
        raise ValueError(f"{work}: pair-profile identities are duplicated or incomplete")
    if len(rows) != len(pair_keys) * 6:
        raise ValueError(f"{work}: expected six annular profiles for every target/control pair")
    if expected_pairs is not None and pair_keys != expected_pairs:
        raise ValueError(f"{work}: pair identities/order differ from the K=2.5 run")
    target_coords = []
    control_coords = []
    for target_id, control_id, _ in pair_keys:
        if target_id not in target_manifest or control_id not in control_manifest:
            raise ValueError(f"missing sky coordinate for target/control pair {target_id}/{control_id}")
        target_coords.append({"ra_deg": target_manifest[target_id]["ra_deg"],
                              "dec_deg": target_manifest[target_id]["dec_deg"]})
        control_coords.append({"ra_deg": control_manifest[control_id]["ra_deg"],
                               "dec_deg": control_manifest[control_id]["dec_deg"]})
    return rows, pair_keys, target_coords, control_coords, {
        "mask_factor": factor, "work": str(work.resolve()),
        "summary_sha256": sha256(summary_path),
        "verification_sha256": sha256(verification_path),
        "evaluation_sha256": sha256(evaluation_path),
        "profile_sha256": sha256(profiles_path),
        "n_pairs": summary["n_pairs_stacked"],
    }


def run(work_paths: dict[float, Path], output: Path, replicas: int = 5000, seed: int = 20261019,
        cell_width_deg: float = 5.0) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    targets_path = ROOT / "sample_manifest.json"
    controls_path = ROOT / "control_download_01/download_manifest.json"
    targets_doc, controls_doc = read_json(targets_path), read_json(controls_path)
    target_manifest = {str(row["id"]): row for row in targets_doc["clusters"]}
    control_manifest = {str(row["id"]): row for row in controls_doc["clusters"]}
    if len(target_manifest) != 800 or len(control_manifest) != 798:
        raise ValueError("expected the hash-audited 800 target and 798 control manifests")

    metadata = []
    base_rows = None
    base_pair_keys = None
    blocks = None
    results = {}
    for factor in (2.5, 3.5, 5.0):
        rows, pair_keys, target_coords, control_coords, info = load_run(
            work_paths[factor], factor, target_manifest, control_manifest,
            base_pair_keys)
        if base_pair_keys is None:
            base_rows, base_pair_keys = rows, pair_keys
            blocks = connected_pair_blocks(target_coords, control_coords, cell_width_deg)
        else:
            for band in ("i", "i2"):
                a = [(r["target_id"], r["control_id"], r["filter"]) for r in base_rows
                     if r["filter"] == band and (r["r_min_kpc"], r["r_max_kpc"]) == ("50", "100")]
                b = [(r["target_id"], r["control_id"], r["filter"]) for r in rows
                     if r["filter"] == band and (r["r_min_kpc"], r["r_max_kpc"]) == ("50", "100")]
                if a != b:
                    raise ValueError(f"Kron {factor}: pair support differs from K=2.5 for {band}")
        metadata.append(info)
        results[factor] = {}
        for band_index, band in enumerate(("i", "i2")):
            matrix = profile_matrix(rows, pair_keys, band)
            labels = blocks[np.asarray([key[2] == band for key in pair_keys])]
            values = matrix[np.asarray([key[2] == band for key in pair_keys])]
            results[factor][band] = summarize_resamples(
                values, labels, replicas=replicas,
                seed=seed + int(factor * 10) + band_index,
                cell_width_deg=cell_width_deg)

    output.mkdir(parents=True, exist_ok=False)
    rows_out = []
    radii = [{"r_min_kpc": lo, "r_max_kpc": hi} for lo, hi in RADII]
    for factor in (2.5, 3.5, 5.0):
        for band in ("i", "i2"):
            item = results[factor][band]
            for index, radius in enumerate(radii):
                pair_ci = item["pair_bootstrap"]["intervals"][index]["ci95"]
                block_ci = item["spatial_block_bootstrap"]["intervals"][index]["ci95"]
                rows_out.append({
                    "mask_factor": factor, "filter": band, **radius,
                    "n_pairs": item["n_pairs_by_annulus"][index],
                    "n_connected_blocks": item["n_connected_blocks"],
                    "largest_block_pairs": item["largest_block_pairs"],
                    "mean_delta_njy_arcsec2": item["point_mean_njy_arcsec2"][index],
                    "pair_sem_njy_arcsec2": item["pair_sem_njy_arcsec2"][index],
                    "pair_bootstrap_lo95_njy_arcsec2": pair_ci[0] if pair_ci else None,
                    "pair_bootstrap_hi95_njy_arcsec2": pair_ci[1] if pair_ci else None,
                    "block_bootstrap_lo95_njy_arcsec2": block_ci[0] if block_ci else None,
                    "block_bootstrap_hi95_njy_arcsec2": block_ci[1] if block_ci else None,
                    "block_bootstrap_sd_njy_arcsec2": item["spatial_block_bootstrap"]["intervals"][index]["bootstrap_sd"],
                    "block_valid_replicas": item["spatial_block_bootstrap"]["intervals"][index]["valid_replicas"],
                })

    table_path = output / "spatial_bootstrap_results.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

    figures = {}
    colors = {2.5: "#1b6ca8", 3.5: "#c05a27", 5.0: "#4c8c45"}
    x = np.asarray([(lo + hi) / 2 for lo, hi in RADII], dtype=float)
    for band in ("i", "i2"):
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        ax.axhline(0, color="black", linewidth=.8)
        for factor in (2.5, 3.5, 5.0):
            item = results[factor][band]
            y = np.asarray(item["point_mean_njy_arcsec2"])
            interval = item["spatial_block_bootstrap"]["intervals"]
            lo = np.asarray([row["ci95"][0] for row in interval])
            hi = np.asarray([row["ci95"][1] for row in interval])
            ax.plot(x, y, marker="o", color=colors[factor], label=f"Kron {factor:g}")
            ax.fill_between(x, lo, hi, color=colors[factor], alpha=.15)
        ax.set_xlabel("proper radius (kpc)")
        ax.set_ylabel("target − matched control (nJy/arcsec²)")
        ax.set_title(f"HSC-{band} ICL profile, connected spatial-block bootstrap")
        ax.legend()
        ax.grid(alpha=.2)
        figure_path = output / f"spatial_bootstrap_{band}.png"
        fig.savefig(figure_path, dpi=160)
        plt.close(fig)
        figures[band] = str(figure_path)

    summary = {
        "schema": "astrafek-hsc-icl-2d-spatial-block-bootstrap-v1",
        "status": "complete_exploratory",
        "mask_factors": [2.5, 3.5, 5.0],
        "radial_bins_kpc": [list(pair) for pair in RADII],
        "n_target_control_pairs": len(base_pair_keys),
        "n_by_filter": {band: sum(key[2] == band for key in base_pair_keys) for band in ("i", "i2")},
        "n_connected_blocks": int(len(np.unique(blocks))),
        "largest_block_pairs": int(max(np.count_nonzero(blocks == group) for group in np.unique(blocks))),
        "cell_width_equatorial_deg": float(cell_width_deg),
        "block_construction": f"transitive connected components of pairs sharing any {cell_width_deg:g}-degree RA/Dec cell among target or matched-control centers; target/control roles share cells",
        "replicas": int(replicas), "seed": int(seed),
        "runs": metadata, "results": results,
        "inputs": {
            "target_manifest": str(targets_path), "target_manifest_sha256": sha256(targets_path),
            "control_manifest": str(controls_path), "control_manifest_sha256": sha256(controls_path),
        },
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "table": str(table_path), "table_sha256": sha256(table_path),
        "figures": figures,
        "science_measurement": False,
        "limitations": [
            "5-degree rectangular RA/Dec cells are an approximate spatial resampling unit and may combine disconnected survey fields",
            "connected components can be large and unequal; pair and block intervals are sensitivity estimates, not detection significance",
            "radial support and contributing pairs change with mask factor and radius",
            "this resampling does not include image-level pixel covariance, PSF, sky transfer, extinction, calibration, or astrophysical covariance",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "n_pairs": summary["n_target_control_pairs"],
                      "n_connected_blocks": summary["n_connected_blocks"],
                      "largest_block_pairs": summary["largest_block_pairs"],
                      "table_sha256": summary["table_sha256"],
                      "science_measurement": False}, indent=2))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-k2p5", type=Path, default=DEFAULT_RUNS[2.5])
    parser.add_argument("--work-k3p5", type=Path, default=DEFAULT_RUNS[3.5])
    parser.add_argument("--work-k5p0", type=Path, default=DEFAULT_RUNS[5.0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    args = parser.parse_args(argv)
    run({2.5: args.work_k2p5.resolve(), 3.5: args.work_k3p5.resolve(),
         5.0: args.work_k5p0.resolve()}, args.output.resolve(), args.replicas,
        args.seed, args.cell_width_deg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
