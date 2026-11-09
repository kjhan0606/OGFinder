#!/usr/bin/env python3
"""Test radial ICL--satellite association on a fixed, complete profile sample.

The input is the saved 94-pair ICL-annulus/satellite-outer-support intersection.
For each filter, clusters are retained only when all seven satellite radial
bins have a valid profile. Every radial relation therefore uses the same
clusters and connected target/control sky blocks. A joint 14-test restricted
Webb wild-cluster bootstrap controls familywise error across filters and bins.
All outputs remain exploratory; this does not correct the photometry.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = DATA / "icl_satellite_supported_94_association_2026-11-05_01"
DEFAULT_SAT = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUT = DATA / "icl_satellite_fixed_profile_2026-11-06_01"
BANDS = ("i", "i2")
RADIAL = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
WEBB = np.asarray([-math.sqrt(1.5), -1., -math.sqrt(.5), math.sqrt(.5), 1., math.sqrt(1.5)])


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        w.writeheader()
        w.writerows(rows)


def cluster_t(beta: np.ndarray, resid: np.ndarray, x: np.ndarray,
              bread: np.ndarray, labels: np.ndarray) -> float:
    groups = np.unique(labels)
    n, k = x.shape
    if n <= k or len(groups) < 2:
        return math.nan
    score = np.vstack([x[labels == group].T @ resid[labels == group] for group in groups])
    correction = (len(groups) / (len(groups) - 1)) * ((n - 1) / (n - k))
    variance = correction * (bread @ (score.T @ score) @ bread)[-1, -1]
    return float(beta[-1] / math.sqrt(variance)) if variance > 0 else math.nan


def prepare(rows: list[dict], band: str, radial: str, min_blocks: int) -> dict:
    y = np.asarray([row["satellite_by_radius"][radial] for row in rows], dtype=float)
    icl = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    logrich = np.log(np.asarray([float(row["richness"]) for row in rows]))
    zsd, rsd, xsd = float(z.std()), float(logrich.std()), float(icl.std())
    if min(zsd, rsd, xsd) <= 0:
        raise ValueError(f"zero-variance covariate in {band}/{radial}")
    z, logrich, x_signal = (z-z.mean())/zsd, (logrich-logrich.mean())/rsd, (icl-icl.mean())/xsd
    nuisance = np.column_stack((np.ones(len(rows)), z, logrich))
    design = np.column_stack((nuisance, x_signal))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError(f"rank-deficient full model in {band}/{radial}")
    blocks = np.asarray([row["block"] for row in rows], dtype=str)
    unique = np.unique(blocks)
    if len(unique) < min_blocks:
        raise ValueError(f"only {len(unique)} blocks in {band}/{radial}")
    bread = np.linalg.inv(design.T @ design)
    null_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_fit = nuisance @ (null_bread @ nuisance.T @ y)
    null_resid = y - null_fit
    leverage = np.sum((nuisance @ null_bread) * nuisance, axis=1)
    adjusted = null_resid / np.sqrt(np.maximum(1-leverage, .05))
    beta = bread @ design.T @ y
    residual = y - design @ beta
    tvalue = cluster_t(beta, residual, design, bread, blocks)
    if not math.isfinite(tvalue):
        raise ValueError(f"invalid CR1 statistic in {band}/{radial}")
    xresid = icl - nuisance @ (null_bread @ nuisance.T @ icl)
    denom = np.linalg.norm(xresid) * np.linalg.norm(null_resid)
    partial_r = float(xresid @ null_resid / denom) if denom > 0 else math.nan
    return {"band": band, "radial": radial, "rows": rows, "y": y,
            "nuisance": nuisance, "design": design, "bread": bread,
            "blocks": blocks, "unique_blocks": unique, "null_fit": null_fit,
            "adjusted_null_resid": adjusted, "beta": beta, "residual": residual,
            "t": tvalue, "slope": float(beta[-1]), "partial_r": partial_r,
            "n": len(rows), "nblocks": len(unique), "largest_block": max(np.count_nonzero(blocks==b) for b in unique),
            "se": abs(float(beta[-1])/tvalue) if tvalue else math.nan}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-dir", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SAT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261110)
    parser.add_argument("--minimum-blocks", type=int, default=6)
    args = parser.parse_args()
    if args.replicas < 1000 or args.minimum_blocks < 6:
        parser.error("replicas must be >=1000 and minimum-blocks >=6")
    assoc_dir, sat_dir, out = args.association_dir.resolve(), args.satellite_dir.resolve(), args.output.resolve()
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    assoc_path, assoc_sum_path = assoc_dir/"supported_pairs.tsv", assoc_dir/"summary.json"
    sat_path, sat_sum_path = sat_dir/"cluster_profiles.tsv", sat_dir/"summary.json"
    assoc_sum, sat_sum = (json.loads(p.read_text(encoding="utf-8")) for p in (assoc_sum_path, sat_sum_path))
    if (assoc_sum.get("status") != "complete_diagnostic" or assoc_sum.get("science_measurement") is not False
            or assoc_sum.get("n_pairs") != 94 or assoc_sum.get("output_sha256", {}).get(str(assoc_path)) != sha(assoc_path)
            or sat_sum.get("status") != "complete" or sat_sum.get("science_measurement") is not False
            or sat_sum.get("output_sha256", {}).get(str(sat_path)) != sha(sat_path)):
        raise ValueError("input summaries/status/hash do not match the pinned exploratory products")
    association = read_tsv(assoc_path)
    profiles = [r for r in read_tsv(sat_path) if r.get("cohort") == "full_sample" and r.get("radial_bin_kpc") in RADIAL]
    if len(association) != 94 or len({(r["cluster_id"], r["filter"]) for r in association}) != 94:
        raise ValueError("expected 94 unique supported association pairs")
    profile_by_key = {}
    for r in profiles:
        key = (str(r["cluster_id"]), r["filter"], r["radial_bin_kpc"])
        if key in profile_by_key:
            raise ValueError(f"duplicate profile row {key}")
        if int(r["n_satellites"]) > 0:
            profile_by_key[key] = float(r["mean_delta_sb_njy_arcsec2"])

    from evaluate_prestack_pair_support import pair_block_labels
    target = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])} for r in association]
    control = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])} for r in association]
    indices, labels, *_ = pair_block_labels(target, control, np.ones(len(association), dtype=bool))
    if len(indices) != len(association):
        raise ValueError("spatial block assignment changed input row order")
    block_by_key = {(r["cluster_id"], r["filter"]): str(label) for r, label in zip(association, labels)}

    samples, counts = {}, {}
    for band in BANDS:
        band_assoc = [r for r in association if r["filter"] == band]
        support_by_bin = [{(r["cluster_id"], band) for r in band_assoc
                           if (r["cluster_id"], band, radial) in profile_by_key} for radial in RADIAL]
        common = set.intersection(*support_by_bin)
        if len(common) < 10:
            raise ValueError(f"fixed seven-bin sample unexpectedly small for {band}: {len(common)}")
        rows = []
        for a in band_assoc:
            key = (a["cluster_id"], band)
            if key not in common:
                continue
            rows.append({**a, "block": block_by_key[key],
                         "satellite_by_radius": {radial: profile_by_key[(key[0], band, radial)] for radial in RADIAL}})
        samples[band] = rows
        counts[band] = {"n_fixed_pairs": len(rows), "n_global_connected_blocks": len({r["block"] for r in rows}),
                        "n_available_by_radius_before_intersection": {
                            radial: len(support_by_bin[i]) for i, radial in enumerate(RADIAL)}}

    groups = {(band, radial): prepare(samples[band], band, radial, args.minimum_blocks)
              for band in BANDS for radial in RADIAL}
    global_blocks = sorted(set(block_by_key.values()))
    rng = np.random.default_rng(args.seed)
    draws = {key: np.empty(args.replicas, dtype=float) for key in groups}
    max_t = np.empty(args.replicas, dtype=float)
    for i in range(args.replicas):
        weights = dict(zip(global_blocks, rng.choice(WEBB, size=len(global_blocks), replace=True)))
        maximum = 0.
        for key, g in groups.items():
            w = np.asarray([weights[b] for b in g["blocks"]])
            ystar = g["null_fit"] + g["adjusted_null_resid"]*w
            beta_star = g["bread"] @ g["design"].T @ ystar
            resid_star = ystar - g["design"] @ beta_star
            tstar = cluster_t(beta_star, resid_star, g["design"], g["bread"], g["blocks"])
            if not math.isfinite(tstar):
                raise ValueError(f"invalid bootstrap t for {key}")
            draws[key][i] = tstar
            maximum = max(maximum, abs(tstar))
        max_t[i] = maximum
    raw = {key: (1+np.count_nonzero(np.abs(draws[key]) >= abs(g["t"]))) /(args.replicas+1)
           for key, g in groups.items()}
    ordered = sorted(raw, key=raw.get)
    holm, running = {}, 0.
    for rank, key in enumerate(ordered):
        running = max(running, (len(ordered)-rank)*raw[key])
        holm[key] = min(1., running)

    records = []
    for key, g in groups.items():
        critical = float(student_t.ppf(.975, g["nblocks"]-1))
        records.append({"filter": key[0], "radial_bin_kpc": key[1],
                        "n_fixed_clusters": g["n"], "n_connected_blocks": g["nblocks"],
                        "largest_block_n": g["largest_block"],
                        "partial_r_z_logrichness": g["partial_r"],
                        "slope_njy_arcsec2_per_sd_icl": g["slope"], "cr1_se": g["se"],
                        "cr1_ci95_low": g["slope"]-critical*g["se"],
                        "cr1_ci95_high": g["slope"]+critical*g["se"], "t_cr1": g["t"],
                        "p_wild_raw": raw[key], "p_holm_14": holm[key],
                        "p_maxT_14": (1+np.count_nonzero(max_t>=abs(g["t"])))/(args.replicas+1)})

    out.mkdir(parents=True, exist_ok=False)
    table = out/"fixed_sample_radial_relations.tsv"
    write_tsv(table, records)
    null_path = out/"maxT_null.tsv"
    write_tsv(null_path, [{"replica": i+1, "max_abs_t_14": float(v)} for i,v in enumerate(max_t)])
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 7.0), sharex=True, sharey=False, constrained_layout=True)
    centers = np.asarray([7.5, 12.5, 17.5, 25., 35., 45., 62.5])
    halfwidth = np.asarray([2.5,2.5,2.5,5.,5.,5.,12.5])
    for ax, band in zip(axes, BANDS):
        band_rows = [r for r in records if r["filter"] == band]
        y = np.asarray([r["slope_njy_arcsec2_per_sd_icl"] for r in band_rows])
        lo = np.asarray([r["cr1_ci95_low"] for r in band_rows])
        hi = np.asarray([r["cr1_ci95_high"] for r in band_rows])
        ax.axhline(0., color="0.35", lw=1, zorder=0)
        ax.errorbar(centers, y, yerr=np.vstack((y-lo, hi-y)), xerr=halfwidth,
                    fmt="o", color="#235789", ecolor="#607d9b", capsize=3, ms=5)
        ax.set_title(f"{band}: fixed N={counts[band]['n_fixed_pairs']} clusters",
                     loc="left", fontsize=10)
        ax.set_ylabel("satellite excess slope\n[nJy arcsec⁻² per 1σ ICL]")
        ax.grid(axis="y", alpha=.2)
        extent = float(np.nanmax(np.abs(np.concatenate((lo, hi)))))
        ax.set_ylim(-1.08*extent, 1.08*extent)
    axes[-1].set_xlabel("Satellite-centered radius [proper kpc]")
    fig.savefig(out/"icl_satellite_fixed_sample_radial_slopes.png", dpi=180)
    fig.savefig(out/"icl_satellite_fixed_sample_radial_slopes.pdf")
    plt.close(fig)
    output_hashes = {p.name: sha(p) for p in (table, null_path, out/"icl_satellite_fixed_sample_radial_slopes.png", out/"icl_satellite_fixed_sample_radial_slopes.pdf")}
    summary = {"schema": "astrafek-icl-satellite-fixed-profile-v1", "status": "complete_exploratory",
               "science_measurement": False, "hsc_network_used": False, "new_image_download": False,
               "method": "same complete seven-bin cluster/filter set within each filter; satellite excess ~ ICL contrast + standardized redshift + standardized log richness; CR1 by joint connected 5-degree target/control block; restricted HC2 Webb six-point wild-cluster bootstrap-t; Holm and single-step maxT across all 14 filter/radius tests",
               "radial_bins_kpc": list(RADIAL), "replicas": args.replicas, "seed": args.seed,
               "fixed_support": counts,
               "input_sha256": {str(p): sha(p) for p in (assoc_path, assoc_sum_path, sat_path, sat_sum_path)},
               "runner_sha256": sha(Path(__file__).resolve()), "output_sha256": output_hashes,
               "limitations": ["94-pair parent is selected by prior ICL and outer-satellite support; complete seven-bin intersection is smaller and not representative of all 800 clusters",
                              "connected equal-area 5-degree cells approximate HSC survey fields",
                              "CR1 intervals and bootstrap do not include source Poisson, coadd weight, sky-transfer, extended-PSF, completeness, or injection-response uncertainty",
                              "14-test familywise correction does not account for earlier exploratory analyses",
                              "fixed-support profile does not imply a fixed underlying satellite population if source masks or member detections vary by radius",
                              "all results remain exploratory"]}
    summary_path = out/"summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "fixed_support": counts,
                      "minimum_maxT_p": min(r["p_maxT_14"] for r in records),
                      "outputs": output_hashes}, indent=2))


if __name__ == "__main__":
    main()
