#!/usr/bin/env python3
"""Fit continuous ICL--weak-lensing profiles on the verified i2 common sample.

Uses only previously approved local HSC shape and random-center exports. For
each supported lens, it measures a target profile minus the pooled profile of
that lens's matched random centers, then regresses the seven radial bins on a
continuous ICL contrast while controlling redshift and log richness. A shared
equal-area delete-one jackknife removes both lens/center and source rows in a
cell. This is an exploratory association analysis, not a calibrated detection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import sha
from fit_wl_random_stack_nfw_mass import DATA, load_pairs_and_provenance
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import (
    delete_block_covariance, equal_area_block_ids,
)
from run_wl_random_shape_stack import verify_and_load_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc, stack_s19a_pairs,
    tomographic_nz_lensing_efficiency,
)
from analyze_wl_random_subtracted_icl_split import (
    ICL_SUPPORT, ICL_WL_OVERLAP, RANDOM_PLAN, RANDOM_ROWS, RANDOM_SUMMARY,
    SAMPLE_MANIFEST, SOURCE_BINS, RADIAL_EDGES_MPC, active_lens_ids, subset,
)

NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_continuous_icl_relation_i2_2026-11-06_01"
FILTER = "i2"
BLOCK_WIDTH_DEG = 5.0
EXPECTED_ICL_SHA = "86e6d453df847ba4c9316c4fdc188565b1bc292bb41b0c5ba8ad310864f02aa5"
EXPECTED_OVERLAP_SHA = "98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12"
EXPECTED_RANDOM_SUMMARY_SHA = "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8"


def sha256(path: Path) -> str:
    return sha(Path(path))


def load_icl_values():
    if sha256(ICL_SUPPORT) != EXPECTED_ICL_SHA:
        raise ValueError("ICL support table changed from the audited input")
    values = {}
    with ICL_SUPPORT.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["filter"] != FILTER or not row["delta_sb_200_300"]:
                continue
            if float(row["coverage_200_300"]) < 0.40:
                continue
            lens_id = str(row["id"])
            if lens_id in values:
                raise ValueError(f"duplicate ICL row: {lens_id}")
            values[lens_id] = float(row["delta_sb_200_300"])
    return values


def get_profile_vector(pairs, beta, cosmology):
    if len(pairs["lens_id"]) == 0:
        return np.full(2 * (len(RADIAL_EDGES_MPC) - 1), np.nan), [0] * (len(RADIAL_EDGES_MPC) - 1)
    result = stack_s19a_pairs(
        pairs, None, None, None, RADIAL_EDGES_MPC, cosmology,
        _pair_lensing_efficiency=beta,
    )
    vals = []
    counts = []
    for row in result["profile"]:
        tangent = row["delta_sigma_weak_msun_mpc2"]
        cross = row["cross_delta_sigma_weak_msun_mpc2"]
        vals.extend((float(tangent) if tangent is not None else np.nan,
                     float(cross) if cross is not None else np.nan))
        counts.append(int(row["n_pairs"]))
    return np.asarray(vals), counts


def design_matrix(ids, meta, icl_values, centers):
    z = np.asarray([float(meta[lens]["z"]) for lens in ids])
    log_rich = np.log([float(meta[lens]["richness"]) for lens in ids])
    icl = np.asarray([icl_values[lens] for lens in ids])
    z0, zscale = float(np.mean(z)), float(np.std(z, ddof=1))
    r0, rscale = float(np.mean(log_rich)), float(np.std(log_rich, ddof=1))
    i0, iscale = float(np.mean(icl)), float(np.std(icl, ddof=1))
    if min(zscale, rscale, iscale) <= 0:
        raise ValueError("cannot standardize constant covariate or ICL metric")
    x = np.column_stack((
        np.ones(len(ids)), (icl - i0) / iscale,
        (z - z0) / zscale, (log_rich - r0) / rscale,
    ))
    return x, {"icl_mean_njy_arcsec2": i0, "icl_sd_njy_arcsec2": iscale,
               "redshift_mean": z0, "redshift_sd": zscale,
               "log_richness_mean": r0, "log_richness_sd": rscale}


def fit_ols_hc3(x, y):
    if len(y) <= x.shape[1] + 1 or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("regression design is rank deficient or too small")
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ (x.T @ y)
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    adjusted = residual / np.maximum(1.0 - leverage, 1e-8)
    meat = x.T @ ((adjusted * adjusted)[:, None] * x)
    cov = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    dof = len(y) - x.shape[1]
    p = 2.0 * student_t.sf(np.abs(beta / np.maximum(se, 1e-300)), dof)
    return beta, se, p


def run(output: Path, *, block_width_deg=BLOCK_WIDTH_DEG):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite continuous ICL/WL output: {output}")
    if sha256(ICL_WL_OVERLAP) != EXPECTED_OVERLAP_SHA:
        raise ValueError("ICL/WL overlap audit changed")
    if sha256(RANDOM_SUMMARY) != EXPECTED_RANDOM_SUMMARY_SHA:
        raise ValueError("random-subtracted stack summary changed")
    overlap = json.loads(ICL_WL_OVERLAP.read_text(encoding="utf-8"))
    if overlap.get("science_measurement") is not False or overlap.get("status") != "complete_exploratory":
        raise ValueError("ICL/WL overlap audit is not complete exploratory evidence")

    target, random_doc, _ = load_pairs_and_provenance()
    random, random_provenance = verify_and_load_pairs(RANDOM_ROWS, RANDOM_PLAN, masked_value=0)
    center_to_lens = {str(k): str(v) for k, v in random_doc["input"]["random_center_matched_lens_ids"].items()}
    expected_center_counts = random_doc["input"]["random_center_expected_source_counts"]
    target_lens_ids = {str(x) for x in np.unique(target["lens_id"])}
    if not target_lens_ids <= set(center_to_lens.values()):
        raise ValueError("a target lens is absent from the approved random-center mapping")

    keep_t = np.isin(target["source_zbin"], SOURCE_BINS)
    keep_r = np.isin(random["source_zbin"], SOURCE_BINS)
    target = subset(target, keep_t)
    random = subset(random, keep_r)
    z_grid, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    target_beta = tomographic_nz_lensing_efficiency(
        target["lens_z"], target["source_zbin"], z_grid, nz_bins, cosmology,
    )
    random_beta = tomographic_nz_lensing_efficiency(
        random["lens_z"], random["source_zbin"], z_grid, nz_bins, cosmology,
    )
    target_active = active_lens_ids(target, z_grid, nz_bins, cosmology)
    random_active_centers = active_lens_ids(random, z_grid, nz_bins, cosmology)
    random_active_lenses = {center_to_lens[c] for c in random_active_centers if c in center_to_lens}
    approved = set(map(str, overlap["icl_profile_support_by_filter"][FILTER]["icl_supported_ids_in_common_wl_support"]))
    icl_values = load_icl_values()
    ids = sorted(target_active & random_active_lenses & approved & set(icl_values))
    expected_ids = set(map(str, overlap["icl_profile_support_by_filter"][FILTER]["icl_supported_ids_in_common_wl_support"]))
    if len(ids) != 70 or set(ids) != expected_ids:
        raise ValueError(f"expected the exact 70-lens i2 ICL/WL common sample, found {len(ids)}")

    manifest = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    meta = {str(row["id"]): row for row in manifest["clusters"]}
    x, scaling = design_matrix(ids, meta, icl_values, None)

    target_by_lens = {lens: subset(target, np.asarray(target["lens_id"]).astype(str) == lens) for lens in ids}
    target_beta_by_lens = {
        lens: target_beta[np.asarray(target["lens_id"]).astype(str) == lens] for lens in ids
    }
    centers_by_lens = {
        lens: sorted(center for center, owner in center_to_lens.items() if owner == lens)
        for lens in ids
    }
    if any(len(centers_by_lens[lens]) != 5 for lens in ids):
        raise ValueError("each common-support lens must map to exactly five random centers")
    random_by_lens = {
        lens: subset(random, np.isin(np.asarray(random["lens_id"]).astype(str), centers_by_lens[lens]))
        for lens in ids
    }
    random_beta_by_lens = {
        lens: random_beta[np.isin(np.asarray(random["lens_id"]).astype(str), centers_by_lens[lens])]
        for lens in ids
    }
    if any(not len(random_by_lens[lens]["lens_id"]) for lens in ids):
        raise ValueError("a common lens has no random-center source rows in bins 2--4")

    endpoints = [target_by_lens[lens] for lens in ids] + [random_by_lens[lens] for lens in ids]
    endpoint_betas = [target_beta_by_lens[lens] for lens in ids] + [random_beta_by_lens[lens] for lens in ids]
    endpoint_cells = []
    for pairs in endpoints:
        lens_cells = equal_area_block_ids(pairs["lens_ra_deg"], pairs["lens_dec_deg"], cell_width_deg=block_width_deg)
        source_cells = equal_area_block_ids(pairs["source_ra_deg"], pairs["source_dec_deg"], cell_width_deg=block_width_deg)
        endpoint_cells.append((lens_cells, source_cells))
    block_ids = np.unique(np.concatenate([v for pair in endpoint_cells for v in pair]))
    if len(block_ids) < 6:
        raise ValueError(f"too few occupied shared sky cells for slope jackknife: {len(block_ids)}")

    def profile_matrix(deleted_block=None):
        outputs = []
        valid_lenses = []
        count_rows = []
        for i, lens in enumerate(ids):
            vectors = []
            endpoint_counts = []
            for endpoint_index in (i, len(ids) + i):
                pairs = endpoints[endpoint_index]
                lens_cells, source_cells = endpoint_cells[endpoint_index]
                if deleted_block is None:
                    keep = np.ones(len(pairs["lens_id"]), dtype=bool)
                else:
                    keep = (lens_cells != deleted_block) & (source_cells != deleted_block)
                selected = subset(pairs, keep)
                value = get_profile_vector(selected, endpoint_betas[endpoint_index][keep], cosmology)
                vectors.append(value[0])
                endpoint_counts.append(value[1])
            outputs.append(vectors[0] - vectors[1])
            valid_lenses.append(lens)
            count_rows.append(endpoint_counts)
        if not outputs:
            raise ValueError("no lenses remain with full radial support")
        return valid_lenses, np.asarray(outputs), count_rows

    full_ids, y, pair_counts = profile_matrix()
    betas, hc3, nominal_p = [], [], []
    for column in range(y.shape[1]):
        valid = np.isfinite(y[:, column])
        b, e, p = fit_ols_hc3(x[valid], y[valid, column])
        betas.append(b)
        hc3.append(e)
        nominal_p.append(p)
    betas = np.asarray(betas)
    hc3 = np.asarray(hc3)
    nominal_p = np.asarray(nominal_p)

    replicate_slopes = []
    deletion_info = []
    for block in block_ids:
        rep_ids, rep_y, _ = profile_matrix(int(block))
        id_index = [ids.index(lens) for lens in rep_ids]
        rep_beta = []
        for column in range(rep_y.shape[1]):
            valid = np.isfinite(rep_y[:, column])
            rep_beta.append(fit_ols_hc3(x[np.asarray(id_index)[valid]], rep_y[valid, column])[0][1])
        replicate_slopes.append(rep_beta)
        deletion_info.append({
            "block_id": int(block),
            "n_lenses_with_any_finite_profile": int(np.count_nonzero(np.any(np.isfinite(rep_y), axis=1))),
            "n_lenses_with_profile_by_bin": [
                int(np.count_nonzero(np.isfinite(rep_y[:, 2 * radial_bin])))
                for radial_bin in range(len(RADIAL_EDGES_MPC) - 1)
            ],
        })
    replicate_slopes = np.asarray(replicate_slopes)
    slope_cov = delete_block_covariance(replicate_slopes)
    slope_jk = np.sqrt(np.maximum(np.diag(slope_cov), 0.0))

    regression_rows = []
    for index in range(len(RADIAL_EDGES_MPC) - 1):
        for component, offset in (("tangential", 0), ("cross", 1)):
            col = 2 * index + offset
            regression_rows.append({
                "bin": index, "r_min_mpc": float(RADIAL_EDGES_MPC[index]),
                "r_max_mpc": float(RADIAL_EDGES_MPC[index + 1]), "component": component,
                "icl_slope_per_1sd_msun_mpc2": float(betas[col, 1]),
                "hc3_standard_error_msun_mpc2": float(hc3[col, 1]),
                "spatial_jackknife_error_msun_mpc2": float(slope_jk[col]),
                "hc3_nominal_two_sided_p": float(nominal_p[col][1]),
                "intercept_msun_mpc2": float(betas[col, 0]),
                "redshift_control_slope_per_1sd_msun_mpc2": float(betas[col, 2]),
                "log_richness_control_slope_per_1sd_msun_mpc2": float(betas[col, 3]),
                "n_lenses": int(np.count_nonzero(np.isfinite(y[:, col]))),
            })

    output.mkdir(parents=True, exist_ok=False)
    with (output / "per_lens_profiles.tsv").open("w", encoding="utf-8", newline="") as stream:
        fields = ["lens_id", "icl_delta_sb_200_300_njy_arcsec2", "redshift", "richness"]
        for i in range(len(RADIAL_EDGES_MPC) - 1):
            fields.extend((f"tangential_bin{i}_msun_mpc2", f"cross_bin{i}_msun_mpc2",
                           f"target_pairs_bin{i}", f"random_pairs_bin{i}"))
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row_i, lens in enumerate(full_ids):
            row = {"lens_id": lens, "icl_delta_sb_200_300_njy_arcsec2": icl_values[lens],
                   "redshift": meta[lens]["z"], "richness": meta[lens]["richness"]}
            for i in range(len(RADIAL_EDGES_MPC) - 1):
                row[f"tangential_bin{i}_msun_mpc2"] = y[row_i, 2*i]
                row[f"cross_bin{i}_msun_mpc2"] = y[row_i, 2*i+1]
                row[f"target_pairs_bin{i}"] = pair_counts[row_i][0][i]
                row[f"random_pairs_bin{i}"] = pair_counts[row_i][1][i]
            writer.writerow(row)
    with (output / "regression.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(regression_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(regression_rows)

    summary = {
        "schema": "astrafek-hsc-continuous-icl-wl-relation-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used": False,
        "new_data_requested": False,
        "filter": FILTER,
        "n_common_support_lenses": len(full_ids),
        "icl_metric": "target-minus-control 200-300 kpc surface-brightness contrast; continuous standardized predictor",
        "dependent_variable": "per-lens target DeltaSigma minus pooled matched-random-center DeltaSigma; HSC Y3 calibrated n(z), source bins 2-4",
        "model": "OLS at each radial bin/component: profile = intercept + beta_ICL*zscore(ICL) + beta_z*zscore(redshift) + beta_lnrich*zscore(ln richness); HC3 and shared equal-area delete-one sky-cell slope covariance",
        "covariate_scaling": scaling,
        "radial_edges_proper_mpc": RADIAL_EDGES_MPC.tolist(),
        "number_of_parameters": 4,
        "regression_rows": regression_rows,
        "slope_covariance_order": "tangential bin0, cross bin0, ..., tangential bin6, cross bin6",
        "icl_slope_spatial_jackknife_covariance": slope_cov.tolist(),
        "n_spatial_delete_one_blocks": int(len(block_ids)),
        "delete_block_details": deletion_info,
        "random_centers_per_lens": 5,
        "source_bin_counts_per_lens": {"target": "stored in per_lens_profiles.tsv", "random": "stored in per_lens_profiles.tsv"},
        "provenance": {
            "icl_support_table_sha256": sha256(ICL_SUPPORT),
            "icl_wl_overlap_summary_sha256": sha256(ICL_WL_OVERLAP),
            "random_stack_summary_sha256": sha256(RANDOM_SUMMARY),
            "random_row_plan_sha256": random_provenance["random_row_plan_sha256"],
            "random_batch_sha256": random_provenance["random_row_batch_sha256"],
            "hsc_y3_nz_sacc_sha256": sha256(NZ_SACC),
            "sample_manifest_sha256": sha256(SAMPLE_MANIFEST),
        },
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "nominal HC3 p-values are conditional model diagnostics and are not calibrated detection significances",
            "the approximate equal-area grid is not a validated contiguous HSC field jackknife",
            "sample includes only 70 i2 lenses with both ICL annulus support and usable target/random shape rows; it does not represent all 800 clusters",
            "ICL values retain unresolved injected-source Poisson, coadd-weight, sky-transfer, completeness, and extended-PSF-wing uncertainties",
            "weak-lensing systematics including photo-z posterior uncertainty, shear calibration, reduced shear, magnification, centering, and spatially varying depth/mask/PSF are not fully propagated",
            "linear conditional association does not establish causal direction or a physical ICL-to-halo-mass relation",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "n_lenses": len(full_ids),
                      "n_spatial_blocks": len(block_ids),
                      "first_tangential_slopes": [x["icl_slope_per_1sd_msun_mpc2"] for x in regression_rows if x["component"] == "tangential"]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-width-deg", type=float, default=BLOCK_WIDTH_DEG)
    args = parser.parse_args()
    run(args.output, block_width_deg=args.block_width_deg)


if __name__ == "__main__":
    main()
