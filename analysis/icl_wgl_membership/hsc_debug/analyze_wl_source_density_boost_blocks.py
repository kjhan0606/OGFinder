#!/usr/bin/env python3
"""Recompute WL source-density boosts and joint boost-DeltaSigma covariance.

Reads the already approved, checksum-verified local HSC source-row exports.
It makes no HSC request and does not alter the saved DeltaSigma stack. The
output remains a diagnostic because cluster/random aperture matching and
photo-z dilution calibration are not yet science-validated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[2]))

from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from run_wl_random_shape_stack import verify_and_load_pairs as load_random_pairs
from run_wl_shape_stack import verify_and_load_pairs as load_cluster_pairs
from s19a_shape_stack import (
    C_KMS,
    G_MPC,
    _to_tangential,
    load_hsc_y3_nz_sacc,
    pair_geometry,
    tomographic_nz_lensing_efficiency,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RANDOM_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_ROWS = RANDOM_PLAN / "results"
CLUSTER_ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
CLUSTER_PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
CLUSTER_PLAN = DATA / "wl_shape_rows_plan_2026-10-18_04"
RANDOM_CENTERS = DATA / "wl_random_centers_lens_reuse_2026-11-03_01/random_centers.json"
SAMPLE_MANIFEST = DATA / "sample_manifest.json"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
STACK_5 = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
STACK_2P5 = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_02p5/summary.json"
MASKED_BMODE_VALUE = 0
RADIAL_EDGES = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0], dtype=float)
METRICS = (
    "pair_count", "shape_weight", "sigma_crit_weight", "sigma_crit_weight_rms2",
    "sigma_crit_weight_mag_edge", "sigma_crit_weight_resolution_edge",
    "response_weight_without_responsivity", "response_lensing_weight_without_responsivity",
    "sigma_crit_weight_psf_e1", "sigma_crit_weight_psf_e2",
    "ds_w_beta_1pm_e_t", "ds_w_beta_1pm_e_x",
    "ds_w_beta_1pm_c_t", "ds_w_beta_1pm_c_x",
    "ds_w_beta_1pm2", "ds_w_beta_1pm2_cos2phi", "ds_w_beta_1pm2_sin2phi",
    "ds_w_beta2_1pm2",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reduce_pairs(pairs: dict, keep_lenses: set[str], center_to_lens: dict[str, str] | None,
                 z_nz: np.ndarray, nz_by_bin: np.ndarray, cosmology: Cosmology,
                 width: float) -> dict:
    if center_to_lens is None:
        keep = np.isin(pairs["lens_id"].astype(str), list(keep_lenses))
    else:
        keep_centers = {center for center, lens in center_to_lens.items() if lens in keep_lenses}
        keep = np.isin(pairs["lens_id"].astype(str), list(keep_centers))
    indices = np.flatnonzero(keep)
    if not len(indices):
        raise ValueError("selected matched sample has no source rows")
    subset = {key: values[indices] for key, values in pairs.items()}
    beta = tomographic_nz_lensing_efficiency(
        subset["lens_z"], subset["source_zbin"], z_nz, nz_by_bin, cosmology)
    sigma_crit_inv = ((4.0 * np.pi * G_MPC / (C_KMS * C_KMS))
                      * np.asarray(cosmology.angular_diameter(subset["lens_z"]), dtype=float)
                      * beta)
    radius, phi = pair_geometry(subset, cosmology)
    radial = np.searchsorted(RADIAL_EDGES, radius, side="right") - 1
    radial[radius == RADIAL_EDGES[-1]] = len(RADIAL_EDGES) - 2
    in_range = ((radius >= RADIAL_EDGES[0]) & (radius <= RADIAL_EDGES[-1])
                & (radial >= 0) & (radial < len(RADIAL_EDGES) - 1))
    shape = ((subset["sigma_e"] >= 0.0) & (subset["sigma_e"] <= 0.4)
             & (np.hypot(subset["e1"], subset["e2"]) < 2.0)
             & (subset["mag_a10"] <= 25.5) & (subset["resolution"] >= 0.3)
             & subset["b_mode_pass"] & (beta > 0.0))
    use = in_range & shape
    lens_ra = subset["lens_ra_deg"][use]
    lens_dec = subset["lens_dec_deg"][use]
    source_ra = subset["source_ra_deg"][use]
    source_dec = subset["source_dec_deg"][use]
    radial_bin = radial[use]
    beta_k = sigma_crit_inv[use]
    phi = phi[use]
    weight = subset["weight"][use]
    rms = subset["rms_e"][use]
    m = subset["m"][use]
    mag = subset["mag_a10"][use]
    resolution = subset["resolution"][use]
    e_t, e_x = _to_tangential(phi, subset["e1"][use], subset["e2"][use])
    c_t, c_x = _to_tangential(phi, subset["c1"][use], subset["c2"][use])
    psf_trace = subset["psf_shape11"][use] + subset["psf_shape22"][use]
    if np.any(psf_trace <= 0):
        raise ValueError("selected HSC PSF moment trace must be positive")
    psf_e1 = (subset["psf_shape11"][use] - subset["psf_shape22"][use]) / psf_trace
    psf_e2 = 2.0 * subset["psf_shape12"][use] / psf_trace
    all_endpoints = np.concatenate((lens_ra, source_ra))
    all_decs = np.concatenate((lens_dec, source_dec))
    raw_blocks = equal_area_block_ids(all_endpoints, all_decs, cell_width_deg=width)
    raw_lens_blocks, raw_source_blocks = raw_blocks[:len(radial_bin)], raw_blocks[len(radial_bin):]
    unique_blocks = np.unique(raw_blocks)
    lens_block = np.searchsorted(unique_blocks, raw_lens_blocks)
    source_block = np.searchsorted(unique_blocks, raw_source_blocks)

    n_radial = len(RADIAL_EDGES) - 1
    n_blocks = len(unique_blocks)
    flat = ((radial_bin * n_blocks + lens_block) * n_blocks + source_block)
    shape_weight = weight
    q = weight * beta_k * beta_k
    values = {
        "pair_count": np.ones(len(weight), dtype=float),
        "shape_weight": shape_weight,
        "sigma_crit_weight": q,
        "sigma_crit_weight_rms2": q * rms * rms,
        "sigma_crit_weight_mag_edge": q * (mag >= (25.5 - 0.025)),
        "sigma_crit_weight_resolution_edge": q * (resolution <= (0.3 + 0.01)),
        "response_weight_without_responsivity": weight * (1.0 + m) ** 2,
        "response_lensing_weight_without_responsivity": weight * (1.0 + m) ** 2 * beta_k * beta_k,
        "sigma_crit_weight_psf_e1": q * psf_e1,
        "sigma_crit_weight_psf_e2": q * psf_e2,
        "ds_w_beta_1pm_e_t": weight * beta_k * (1.0 + m) * e_t,
        "ds_w_beta_1pm_e_x": weight * beta_k * (1.0 + m) * e_x,
        "ds_w_beta_1pm_c_t": weight * beta_k * (1.0 + m) * c_t,
        "ds_w_beta_1pm_c_x": weight * beta_k * (1.0 + m) * c_x,
        "ds_w_beta_1pm2": weight * beta_k * (1.0 + m) ** 2,
        "ds_w_beta_1pm2_cos2phi": weight * beta_k * (1.0 + m) ** 2 * np.cos(2.0 * phi),
        "ds_w_beta_1pm2_sin2phi": weight * beta_k * (1.0 + m) ** 2 * np.sin(2.0 * phi),
        "ds_w_beta2_1pm2": weight * beta_k * beta_k * (1.0 + m) ** 2,
    }
    tensor = np.zeros((n_radial, n_blocks, n_blocks, len(METRICS)), dtype=np.float64)
    for metric_index, metric_name in enumerate(METRICS):
        tensor[..., metric_index] = np.bincount(
            flat, weights=values[metric_name], minlength=n_radial * n_blocks * n_blocks
        ).reshape(n_radial, n_blocks, n_blocks)

    endpoint_lenses = np.asarray(pairs["lens_id"][indices], dtype=str)
    unique_center_rows = {}
    for ident, ra, dec in zip(endpoint_lenses, pairs["lens_ra_deg"][indices], pairs["lens_dec_deg"][indices]):
        unique_center_rows.setdefault(ident, (float(ra), float(dec)))
    center_coords = np.asarray([unique_center_rows[ident] for ident in sorted(unique_center_rows)], dtype=float)
    total = tensor.sum(axis=(1, 2))
    return {
        "tensor": tensor, "blocks": unique_blocks,
        "selected_pair_counts": total[:, METRICS.index("pair_count")].astype(int),
        "full_metrics": total,
        "n_input_rows_in_matched_sample": int(len(indices)),
        "n_selected_rows": int(use.sum()),
        "endpoint_center_coords": center_coords,
        "radial_bins_mpc": [[float(a), float(b)] for a, b in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])],
    }


def source_centers_for_sample(random_center_doc: dict, center_to_lens: dict[str, str],
                              common_lenses: set[str]) -> list[dict]:
    selected = [row for row in random_center_doc["centers"]
                if str(row["random_center_id"]) in center_to_lens
                and center_to_lens[str(row["random_center_id"])] in common_lenses]
    if len(selected) != 5 * len(common_lenses):
        raise ValueError("random-center positions do not contain exactly five centers per common lens")
    return selected


def delta_sigma_from_totals(total: np.ndarray) -> np.ndarray:
    """Evaluate the calibrated HSC S19A weak-shear estimator from block sums."""
    q = total[:, METRICS.index("sigma_crit_weight")]
    q_rms2 = total[:, METRICS.index("sigma_crit_weight_rms2")]
    q_mag = total[:, METRICS.index("sigma_crit_weight_mag_edge")]
    q_res = total[:, METRICS.index("sigma_crit_weight_resolution_edge")]
    q_psf_e1 = total[:, METRICS.index("sigma_crit_weight_psf_e1")]
    q_psf_e2 = total[:, METRICS.index("sigma_crit_weight_psf_e2")]
    if np.any(q <= 0):
        raise ValueError("every jackknife radial bin must retain positive SigmaCrit weight")

    responsivity = 1.0 - q_rms2 / q
    p_mag = q_mag / q / 0.025
    p_res = q_res / q / 0.01
    m_sel = -0.059 * p_mag + 0.019 * p_res
    a_sel = 0.0064 * p_mag + 0.0063 * p_res
    c1_sel = a_sel * q_psf_e1 / q
    c2_sel = a_sel * q_psf_e2 / q
    if np.any(responsivity <= 0) or np.any(1.0 + m_sel <= 0):
        raise ValueError("jackknife radial bin has a non-positive S19A calibrated response")

    e_t = total[:, METRICS.index("ds_w_beta_1pm_e_t")]
    e_x = total[:, METRICS.index("ds_w_beta_1pm_e_x")]
    c_t = total[:, METRICS.index("ds_w_beta_1pm_c_t")]
    c_x = total[:, METRICS.index("ds_w_beta_1pm_c_x")]
    cos_sum = total[:, METRICS.index("ds_w_beta_1pm2_cos2phi")]
    sin_sum = total[:, METRICS.index("ds_w_beta_1pm2_sin2phi")]
    denominator = (4.0 * responsivity**2 * (1.0 + m_sel)**2
                   * total[:, METRICS.index("ds_w_beta2_1pm2")])
    selection_c_t = -c1_sel * cos_sum - c2_sel * sin_sum
    selection_c_x = c1_sel * sin_sum - c2_sel * cos_sum
    numerator_t = (
        2.0 * responsivity * (1.0 + m_sel) * e_t
        - 4.0 * responsivity**2 * (1.0 + m_sel) * c_t
        - 4.0 * responsivity**2 * (1.0 + m_sel) * selection_c_t
    )
    numerator_x = (
        2.0 * responsivity * (1.0 + m_sel) * e_x
        - 4.0 * responsivity**2 * (1.0 + m_sel) * c_x
        - 4.0 * responsivity**2 * (1.0 + m_sel) * selection_c_x
    )
    if np.any(denominator <= 0) or np.any(~np.isfinite(denominator)):
        raise ValueError("jackknife radial bin has a non-positive DeltaSigma denominator")
    result = np.empty((len(q), 2), dtype=float)
    result[:, 0] = numerator_t / denominator
    result[:, 1] = numerator_x / denominator
    if not np.all(np.isfinite(result)):
        raise ValueError("sufficient-statistic DeltaSigma profile is non-finite")
    return result


def delta_sigma_leave_one_block(jackknife: dict) -> tuple[np.ndarray, np.ndarray]:
    tensor = jackknife["tensor"]
    total = tensor.sum(axis=(1, 2))
    full = delta_sigma_from_totals(total)
    replicas = []
    for block in range(tensor.shape[1]):
        kept = (total - tensor[:, block, :, :].sum(axis=1)
                - tensor[:, :, block, :].sum(axis=1) + tensor[:, block, block, :])
        replicas.append(delta_sigma_from_totals(kept))
    return full, np.asarray(replicas)


def jackknife_covariance(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] < 3 or not np.all(np.isfinite(values)):
        raise ValueError("jackknife replicas must be a finite [>=3 blocks, parameters] array")
    centered = values - values.mean(axis=0, keepdims=True)
    return ((values.shape[0] - 1.0) / values.shape[0]) * centered.T @ centered


def totals_and_replicas(tensor: np.ndarray, center_blocks: np.ndarray,
                        n_centers: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    total = tensor.sum(axis=(1, 2))
    n_radial = tensor.shape[0]
    n_blocks = tensor.shape[1]
    block_center_counts = np.bincount(center_blocks, minlength=n_blocks).astype(float)
    metric_replicas = {name: np.full((n_blocks, n_radial), np.nan) for name in
                       ("count", "shape_weight", "sigma_crit_weight", "estimator_weight")}
    for block in range(n_blocks):
        kept = total - tensor[:, block, :, :].sum(axis=1) - tensor[:, :, block, :].sum(axis=1) + tensor[:, block, block, :]
        remaining = n_centers - block_center_counts[block]
        if remaining <= 0:
            continue
        count = kept[:, 0]
        shape_w = kept[:, 1]
        q, q_rms2, q_mag, q_res = (kept[:, i] for i in (2, 3, 4, 5))
        response_q = kept[:, 7]
        response = np.divide(q - q_rms2, q, out=np.full(n_radial, np.nan), where=q > 0)
        p_mag = np.divide(q_mag, q * 0.025, out=np.full(n_radial, np.nan), where=q > 0)
        p_res = np.divide(q_res, q * 0.01, out=np.full(n_radial, np.nan), where=q > 0)
        msel = -0.059 * p_mag + 0.019 * p_res
        estimator_w = response_q * (2.0 * response * (1.0 + msel)) ** 2
        metric_replicas["count"][block] = count / remaining
        metric_replicas["shape_weight"][block] = shape_w / remaining
        metric_replicas["sigma_crit_weight"][block] = q / remaining
        metric_replicas["estimator_weight"][block] = estimator_w / remaining
    return total, metric_replicas


def jackknife_summary(lens: dict, random: dict, n_blocks: int) -> tuple[dict, dict]:
    lens_total, lens_rep = totals_and_replicas(
        lens["tensor"], lens["cluster_center_blocks"], lens["n_lens_centers"])
    random_total, random_rep = totals_and_replicas(
        random["tensor"], random["random_center_blocks"], random["n_random_centers"])
    out = {}
    replica_columns = []
    parameter_order = []
    for metric, name in (("count", "count_boost"), ("shape_weight", "shape_weight_boost"),
                         ("sigma_crit_weight", "sigma_crit_weight_boost"),
                         ("estimator_weight", "estimator_weight_boost")):
        lens_total_metric = lens_total[:, METRICS.index({
            "count": "pair_count", "shape_weight": "shape_weight",
            "sigma_crit_weight": "sigma_crit_weight",
            "estimator_weight": "response_lensing_weight_without_responsivity",
        }[metric])]
        if metric == "estimator_weight":
            q = lens_total[:, METRICS.index("sigma_crit_weight")]
            # Recalculate the exact full-sample S19A estimator weight after
            # deriving responsivity and selection-bias terms in each bin.
            q_rms2 = lens_total[:, METRICS.index("sigma_crit_weight_rms2")]
            q_mag = lens_total[:, METRICS.index("sigma_crit_weight_mag_edge")]
            q_res = lens_total[:, METRICS.index("sigma_crit_weight_resolution_edge")]
            resp = np.divide(q - q_rms2, q, out=np.full(len(q), np.nan), where=q > 0)
            p_mag = np.divide(q_mag, q * 0.025, out=np.full(len(q), np.nan), where=q > 0)
            p_res = np.divide(q_res, q * 0.01, out=np.full(len(q), np.nan), where=q > 0)
            msel = -0.059 * p_mag + 0.019 * p_res
            lens_total_metric = lens_total[:, METRICS.index("response_lensing_weight_without_responsivity")] * (2.0 * resp * (1.0 + msel)) ** 2
            q = random_total[:, METRICS.index("sigma_crit_weight")]
            q_rms2 = random_total[:, METRICS.index("sigma_crit_weight_rms2")]
            q_mag = random_total[:, METRICS.index("sigma_crit_weight_mag_edge")]
            q_res = random_total[:, METRICS.index("sigma_crit_weight_resolution_edge")]
            resp_r = np.divide(q - q_rms2, q, out=np.full(len(q), np.nan), where=q > 0)
            p_mag_r = np.divide(q_mag, q * 0.025, out=np.full(len(q), np.nan), where=q > 0)
            p_res_r = np.divide(q_res, q * 0.01, out=np.full(len(q), np.nan), where=q > 0)
            msel_r = -0.059 * p_mag_r + 0.019 * p_res_r
            random_total_metric = random_total[:, METRICS.index("response_lensing_weight_without_responsivity")] * (2.0 * resp_r * (1.0 + msel_r)) ** 2
        else:
            random_total_metric = random_total[:, METRICS.index({
                "count": "pair_count", "shape_weight": "shape_weight",
                "sigma_crit_weight": "sigma_crit_weight",
            }[metric])]
        point = (lens_total_metric / lens["n_lens_centers"]) / (random_total_metric / random["n_random_centers"])
        lens_draws = lens_rep[metric]
        random_draws = random_rep[metric]
        draws = lens_draws / random_draws
        valid = np.isfinite(draws)
        if not valid.all():
            raise ValueError(f"{name} has non-finite leave-one-block replicas; covariance is undefined")
        n_valid = valid.sum(axis=0)
        mean = np.nanmean(draws, axis=0)
        se = np.sqrt(np.divide(
            (n_valid - 1.0) * np.nansum((draws - mean) ** 2, axis=0), n_valid,
            out=np.full(len(point), np.nan), where=n_valid > 1))
        out[name] = [{
            "point": float(point[i]), "delete_one_block_sigma": float(se[i]),
            "delete_one_block_min": float(np.nanmin(draws[:, i])),
            "delete_one_block_max": float(np.nanmax(draws[:, i])),
            "valid_blocks": int(n_valid[i]),
        } for i in range(len(point))]
        replica_columns.append(draws)
        parameter_order.extend({"metric": name, "radial_bin": i} for i in range(len(point)))

    joint_draws = np.concatenate(replica_columns, axis=1)
    centered = joint_draws - joint_draws.mean(axis=0, keepdims=True)
    covariance = ((n_blocks - 1.0) / n_blocks) * (centered.T @ centered)
    expected_variance = np.concatenate([
        np.square([row["delete_one_block_sigma"] for row in out[name]])
        for name in ("count_boost", "shape_weight_boost", "sigma_crit_weight_boost",
                     "estimator_weight_boost")
    ])
    if not np.allclose(np.diag(covariance), expected_variance, rtol=2e-10, atol=1e-14):
        raise ValueError("joint boost covariance diagonal does not reproduce per-bin jackknife errors")
    sigma = np.sqrt(np.diag(covariance))
    denominator = np.outer(sigma, sigma)
    correlation = np.divide(covariance, denominator, out=np.zeros_like(covariance),
                            where=denominator > 0)
    joint = {
        "method": "paired delete-one equal-area sky-cell jackknife covariance across all boost metrics and radial bins",
        "n_blocks": int(n_blocks),
        "parameter_order": parameter_order,
        "covariance": covariance.tolist(),
        "correlation": correlation.tolist(),
    }
    return out, joint, joint_draws


def run(output: Path, bmode_masked_value: int = MASKED_BMODE_VALUE) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_25 = json.loads(STACK_2P5.read_text(encoding="utf-8"))
    summary_5 = json.loads(STACK_5.read_text(encoding="utf-8"))
    if any(s.get("science_measurement") is not False for s in (summary_25, summary_5)):
        raise ValueError("input stacks must retain science_measurement=false")
    for item in (summary_25, summary_5):
        source = item.get("input", {})
        if (source.get("n_z_sacc_sha256") != sha256(NZ_SACC)
                or source.get("b_mode_masked_value") != MASKED_BMODE_VALUE
                or source.get("cosmology") != {"H0_km_s_Mpc": 70.0, "Omega_m": 0.3}):
            raise ValueError("saved stack n(z), B-mode polarity, or cosmology differs from this audit")
    random_pairs, random_prov = load_random_pairs(RANDOM_ROWS, RANDOM_PLAN, bmode_masked_value)
    lens_pairs, lens_prov = load_cluster_pairs(
        CLUSTER_ROWS, CLUSTER_PREFLIGHT, CLUSTER_PLAN, bmode_masked_value)
    matched_lens_ids = set(random_prov["matched_lens_ids"])
    lens_ids_with_rows = set(lens_pairs["lens_id"].astype(str))
    random_support_lens_ids = {
        random_prov["random_center_matched_lens_ids"][center]
        for center, count in random_prov["random_center_expected_source_counts"].items()
        if int(count) > 0
    }
    common_lens_ids = matched_lens_ids & lens_ids_with_rows & random_support_lens_ids
    matched = summary_25["input"]["matched_sample"]
    if len(common_lens_ids) != int(matched["n_common_shape_support_lenses"]):
        raise ValueError("row-derived common lens sample does not match the saved 435-lens stack")
    if int(matched["n_random_centers_in_common_support"]) != 5 * len(common_lens_ids):
        raise ValueError("saved random-center count is not five per common lens")
    common_sha = hashlib.sha256("\n".join(sorted(common_lens_ids)).encode("utf-8")).hexdigest()
    if common_sha != matched.get("common_lens_ids_sha256"):
        raise ValueError("row-derived common lenses do not match the pinned 435-lens sample")
    if (random_prov["random_row_plan_sha256"] != summary_25["input"].get("random_row_plan_sha256")
            or lens_prov["row_plan_sha256"] != summary_25["input"].get("lens_row_plan_sha256")):
        raise ValueError("verified row exports differ from the saved matched-stack plans")

    center_doc = json.loads(RANDOM_CENTERS.read_text(encoding="utf-8"))
    centers = source_centers_for_sample(
        center_doc, random_prov["random_center_matched_lens_ids"], common_lens_ids)
    z_nz, nz_by_bin = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmo = Cosmology(H0=70.0, Om=0.3)
    lens_data = {}
    random_data = {}
    for width, saved in ((5.0, summary_5), (2.5, summary_25)):
        lens = reduce_pairs(lens_pairs, common_lens_ids, None, z_nz, nz_by_bin, cosmo, width)
        center_to_lens = random_prov["random_center_matched_lens_ids"]
        random = reduce_pairs(random_pairs, common_lens_ids, center_to_lens,
                              z_nz, nz_by_bin, cosmo, width)
        expected_profile = saved["paired_spatial_jackknife"]["profile"]
        lens_counts = [int(row["n_lens_pairs"]) for row in expected_profile]
        random_counts = [int(row["n_random_pairs"]) for row in expected_profile]
        if lens["selected_pair_counts"].tolist() != lens_counts:
            raise ValueError(f"recomputed cluster source counts disagree at {width:g}-degree scale")
        if random["selected_pair_counts"].tolist() != random_counts:
            raise ValueError(f"recomputed random source counts disagree at {width:g}-degree scale")

        lens_centers = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))["clusters"]
        lens_by_id = {str(row["id"]): row for row in lens_centers}
        lens_positions = [lens_by_id[ident] for ident in sorted(common_lens_ids)]
        lens_center_blocks = equal_area_block_ids(
            [row["ra_deg"] for row in lens_positions],
            [row["dec_deg"] for row in lens_positions], cell_width_deg=width)
        random_positions = np.asarray([[row["ra_deg"], row["dec_deg"]] for row in centers], dtype=float)
        random_center_blocks = equal_area_block_ids(
            random_positions[:, 0], random_positions[:, 1], cell_width_deg=width)

        # Rebuild a common endpoint-block index, including the 129 random
        # centers with zero source rows so the density denominator is spatially
        # deleted with the same planned center sample.
        pair_block_ids = np.unique(np.concatenate((
            lens["blocks"], random["blocks"], lens_center_blocks, random_center_blocks)))
        expanded_tensors = []
        for dataset in (lens, random):
            old_to_new = np.searchsorted(pair_block_ids, dataset["blocks"])
            expanded = np.zeros((len(RADIAL_EDGES) - 1, len(pair_block_ids), len(pair_block_ids), len(METRICS)))
            for i in range(len(old_to_new)):
                for j in range(len(old_to_new)):
                    expanded[:, old_to_new[i], old_to_new[j], :] = dataset["tensor"][:, i, j, :]
            expanded_tensors.append(expanded)
        lens_center_blocks = np.searchsorted(pair_block_ids, lens_center_blocks)
        random_center_blocks = np.searchsorted(pair_block_ids, random_center_blocks)
        lens_jk = {"tensor": expanded_tensors[0], "cluster_center_blocks": lens_center_blocks,
                   "n_lens_centers": len(common_lens_ids), "block_ids": pair_block_ids}
        random_jk = {"tensor": expanded_tensors[1], "random_center_blocks": random_center_blocks,
                     "n_random_centers": len(centers), "block_ids": pair_block_ids}
        lens_data[str(width)] = {
            "diagnostic": lens, "jackknife": lens_jk,
            "n_endpoint_blocks": len(pair_block_ids),
            "saved_stack": saved,
        }
        random_data[str(width)] = {"diagnostic": random, "jackknife": random_jk}

    profile = []
    summaries = {}
    for width in (5.0, 2.5):
        key = str(width)
        lens = lens_data[key]["diagnostic"]
        random = random_data[key]["diagnostic"]
        lens_jk = lens_data[key]["jackknife"]
        random_jk = random_data[key]["jackknife"]
        metrics, boost_covariance, boost_replicas = jackknife_summary(
            lens_jk, random_jk, lens_data[key]["n_endpoint_blocks"])
        lens_ds_full, lens_ds_replicas = delta_sigma_leave_one_block(lens_jk)
        random_ds_full, random_ds_replicas = delta_sigma_leave_one_block(random_jk)
        ds_full = lens_ds_full - random_ds_full
        ds_replicas = lens_ds_replicas - random_ds_replicas
        saved_jk = lens_data[key]["saved_stack"]["paired_spatial_jackknife"]
        expected_ds = np.asarray([
            [row["delta_sigma_random_subtracted_msun_mpc2"],
             row["cross_delta_sigma_random_subtracted_msun_mpc2"]]
            for row in saved_jk["profile"]
        ], dtype=float)
        if not np.allclose(ds_full, expected_ds, rtol=5e-9, atol=1e-8):
            raise ValueError(f"sufficient-statistic DeltaSigma point differs from saved {width:g}-degree stack")
        block_ids = np.asarray(lens_jk["block_ids"], dtype=np.int64)
        saved_block_ids = [int(row["block_id"]) for row in saved_jk["block_deletions"]]
        if len(saved_block_ids) != int(saved_jk["n_occupied_blocks"]):
            raise ValueError("saved DeltaSigma block deletion list is inconsistent")
        block_index = {int(block): index for index, block in enumerate(block_ids)}
        if not set(saved_block_ids).issubset(block_index):
            raise ValueError("saved DeltaSigma jackknife blocks differ from boost endpoint blocks")
        saved_indices = np.asarray([block_index[block] for block in saved_block_ids], dtype=int)
        saved_block_covariance = jackknife_covariance(ds_replicas[saved_indices].reshape(len(saved_indices), -1))
        expected_block_covariance = np.asarray(saved_jk["covariance"], dtype=float)
        if expected_block_covariance.shape != saved_block_covariance.shape or not np.allclose(
                saved_block_covariance, expected_block_covariance, rtol=2e-6, atol=1e-8):
            raise ValueError(f"sufficient-statistic jackknife covariance differs from saved {width:g}-degree stack")
        extra_indices = np.asarray([
            index for index, block in enumerate(block_ids) if int(block) not in set(saved_block_ids)
        ], dtype=int)
        if any(np.any(lens_jk["tensor"][:, index, :, :] != 0)
               or np.any(random_jk["tensor"][:, index, :, :] != 0)
               or np.any(lens_jk["tensor"][:, :, index, :] != 0)
               or np.any(random_jk["tensor"][:, :, index, :] != 0)
               for index in extra_indices):
            raise ValueError("additional density-denominator blocks unexpectedly contain shear pairs")

        combined_replicas = np.concatenate((boost_replicas, ds_replicas.reshape(len(ds_replicas), -1)), axis=1)
        combined_covariance = jackknife_covariance(combined_replicas)
        combined_sigma = np.sqrt(np.diag(combined_covariance))
        combined_correlation = combined_covariance / np.outer(combined_sigma, combined_sigma)
        ds_parameter_order = [
            {"metric": metric, "radial_bin": radial_bin}
            for radial_bin in range(len(RADIAL_EDGES) - 1)
            for metric in ("delta_sigma_random_subtracted_msun_mpc2",
                           "cross_delta_sigma_random_subtracted_msun_mpc2")
        ]
        joint_covariance = {
            "method": "paired delete-one equal-area sky-cell jackknife of boost and saved-stack-calibrated DeltaSigma sufficient statistics",
            "n_blocks": int(len(block_ids)),
            "n_saved_stack_pair_endpoint_blocks": int(len(saved_block_ids)),
            "additional_center_only_blocks": [int(block_ids[index]) for index in extra_indices],
            "saved_stack_point_reproduced": True,
            "saved_stack_covariance_reproduced_on_pair_endpoint_blocks": True,
            "parameter_order": boost_covariance["parameter_order"] + ds_parameter_order,
            "covariance": combined_covariance.tolist(),
            "correlation": combined_correlation.tolist(),
        }
        # Diagnostic only: exact paired leave-one-block replicas of the
        # candidate product B_W(R) * [DeltaSigma_lens - DeltaSigma_random].
        # The ratio is not a validated dilution correction, so the saved
        # profile and science flags remain unchanged.
        estimator_weight_columns = [
            index for index, item in enumerate(boost_covariance["parameter_order"])
            if item["metric"] == "estimator_weight_boost"
        ]
        if len(estimator_weight_columns) != len(RADIAL_EDGES) - 1:
            raise ValueError("expected one estimator-weight boost per radial bin")
        boost_replicas_by_bin = boost_replicas[:, estimator_weight_columns]
        boosted_replicas = ds_replicas * boost_replicas_by_bin[:, :, None]
        boost_points = np.asarray([item["point"] for item in metrics["estimator_weight_boost"]])
        boosted_full = ds_full * boost_points[:, None]
        boosted_covariance = jackknife_covariance(
            boosted_replicas.reshape(len(boosted_replicas), -1))
        boosted_sigma = np.sqrt(np.diag(boosted_covariance)).reshape(-1, 2)
        boosted_profile = []
        for radial_bin, (lo, hi) in enumerate(zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])):
            boosted_profile.append({
                "bin": radial_bin, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                "estimator_weight_boost": float(boost_points[radial_bin]),
                "delta_sigma_sensitivity_msun_mpc2": float(boosted_full[radial_bin, 0]),
                "cross_delta_sigma_sensitivity_msun_mpc2": float(boosted_full[radial_bin, 1]),
                "delta_sigma_sensitivity_jackknife_sigma_msun_mpc2": float(boosted_sigma[radial_bin, 0]),
                "cross_sensitivity_jackknife_sigma_msun_mpc2": float(boosted_sigma[radial_bin, 1]),
                "fractional_tangential_shift": float(boost_points[radial_bin] - 1.0),
            })
        tangential_indices = np.arange(0, 2 * (len(RADIAL_EDGES) - 1), 2)
        cross_indices = tangential_indices + 1
        boosted_profile_covariance = {
            "method": "exact paired delete-one jackknife of (per-bin response-adjusted estimator-weight ratio × random-subtracted DeltaSigma)",
            "n_blocks": int(len(boosted_replicas)),
            "parameter_order": ds_parameter_order,
            "profile": boosted_profile,
            "covariance": boosted_covariance.tolist(),
            "conditional_chi_square_vs_zero": {
                "tangential": float(boosted_full[:, 0] @ np.linalg.solve(
                    boosted_covariance[np.ix_(tangential_indices, tangential_indices)],
                    boosted_full[:, 0])),
                "cross": float(boosted_full[:, 1] @ np.linalg.solve(
                    boosted_covariance[np.ix_(cross_indices, cross_indices)],
                    boosted_full[:, 1])),
                "dof": int(len(RADIAL_EDGES) - 1),
                "note": "conditional diagnostic only; not a corrected science result",
            },
        }
        summaries[key] = {
            "n_endpoint_blocks": lens_data[key]["n_endpoint_blocks"],
            "n_lenses": lens_jk["n_lens_centers"],
            "n_random_centers": random_jk["n_random_centers"],
            "cluster_input_rows": lens["n_input_rows_in_matched_sample"],
            "random_input_rows": random["n_input_rows_in_matched_sample"],
            "cluster_selected_rows": lens["n_selected_rows"],
            "random_selected_rows": random["n_selected_rows"],
            "boosts_by_metric": metrics,
            "joint_boost_jackknife_covariance": boost_covariance,
            "delta_sigma_random_subtracted_full_sample": ds_full.tolist(),
            "delta_sigma_random_subtracted_jackknife_covariance": jackknife_covariance(
                ds_replicas.reshape(len(ds_replicas), -1)).tolist(),
            "joint_boost_delta_sigma_jackknife_covariance": joint_covariance,
            "unvalidated_boost_profile_sensitivity": boosted_profile_covariance,
        }
    for i, edges in enumerate(zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])):
        row = {"bin": i, "r_min_mpc": float(edges[0]), "r_max_mpc": float(edges[1])}
        for width in ("5.0", "2.5"):
            for metric, values in summaries[width]["boosts_by_metric"].items():
                row[f"{metric}_{width.replace('.', 'p')}deg"] = values[i]["point"]
                row[f"{metric}_{width.replace('.', 'p')}deg_jackknife_sigma"] = values[i]["delete_one_block_sigma"]
                row[f"{metric}_{width.replace('.', 'p')}deg_delete_range"] = (
                    f"{values[i]['delete_one_block_min']:.8g}:{values[i]['delete_one_block_max']:.8g}")
            sensitivity = summaries[width]["unvalidated_boost_profile_sensitivity"]["profile"][i]
            for field in ("estimator_weight_boost", "delta_sigma_sensitivity_msun_mpc2",
                          "cross_delta_sigma_sensitivity_msun_mpc2",
                          "delta_sigma_sensitivity_jackknife_sigma_msun_mpc2",
                          "cross_sensitivity_jackknife_sigma_msun_mpc2",
                          "fractional_tangential_shift"):
                row[f"{field}_{width.replace('.', 'p')}deg"] = sensitivity[field]
        profile.append(row)

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "source_density_boost_blocks.tsv"
    import csv
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(profile[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(profile)
    inputs = (STACK_5, STACK_2P5, RANDOM_CENTERS, SAMPLE_MANIFEST, NZ_SACC,
              RANDOM_PLAN / "random_shape_row_plan.json",
              RANDOM_ROWS / "row_summary.json", RANDOM_ROWS / "row_status.json",
              CLUSTER_PLAN / "row_export_plan.json",
              CLUSTER_ROWS / "row_summary.json", CLUSTER_ROWS / "row_status.json")
    summary = {
        "schema": "astrafek-wl-source-density-boost-joint-covariance-v2",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "boost_correction_applied": False,
        "matched_sample": {
            "n_common_lenses": len(common_lens_ids), "n_random_centers": len(centers),
            "random_centers_per_lens": 5,
            "zero_shape_row_random_centers": sum(
                int(count) == 0 and random_prov["random_center_matched_lens_ids"][center] in common_lens_ids
                for center, count in random_prov["random_center_expected_source_counts"].items()),
            "common_lens_ids_sha256_from_stack": matched["common_lens_ids_sha256"],
        },
        "row_export_plans": {
            "random": random_prov["random_row_plan_sha256"],
            "cluster": lens_prov["row_plan_sha256"],
            "random_row_batches": len(random_prov["random_row_batch_sha256"]),
            "cluster_row_batches": len(lens_prov["row_batch_sha256"]),
        },
        "estimator": "count, shape-weight, SigmaCritInv-squared weight, and exact S19A response-adjusted estimator-weight density ratios per planned lens/random center; SigmaCritInv includes angular-diameter distance",
        "block_method": "paired delete-one equal-area RA/sin(Dec) cells; remove pairs if either lens/random endpoint or source endpoint is in the deleted block; delete the corresponding center positions in the density denominator",
        "n_z_sacc_sha256": sha256(NZ_SACC),
        "results_by_block_width_deg": summaries,
        "profile_rows": profile,
        "profile": str(table_path.resolve()),
        "profile_sha256": sha256(table_path),
        "input_sha256": {str(path.resolve()): sha256(path) for path in inputs},
        "data_requested": False,
        "hsc_network_used": False,
        "images_read": False,
        "limitations": [
            "Boost ratios remain diagnostics until cluster and random apertures are validated to have matched full-aperture depth, mask, PSF, and sky response.",
            "The HSC Y3 tomographic n(z) provides fiducial lensing efficiency; photo-z dilution/calibration uncertainty is not propagated.",
            "A density ratio does not by itself prove cluster-member contamination; magnification, depth, masking, or residual selection differences can also change source counts.",
            "The candidate estimator-weight-boosted profile is an unvalidated sensitivity and does not prove the excess weight is member dilution rather than magnification, depth, mask, PSF, sky, or source-selection effects.",
            "No correction was applied to the saved random-subtracted DeltaSigma science profile; boost_correction_applied remains false and science_measurement remains false.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--b-mode-masked-value", type=int, choices=(0, 1), default=MASKED_BMODE_VALUE)
    args = parser.parse_args()
    result = run(args.output, args.b_mode_masked_value)
    print(json.dumps({
        "status": result["status"], "science_measurement": result["science_measurement"],
        "matched_sample": result["matched_sample"],
        "boosts": [{"radius_mpc": [r["r_min_mpc"], r["r_max_mpc"]],
                    "count_2p5deg": r["count_boost_2p5deg"],
                    "estimator_weight_2p5deg": r["estimator_weight_boost_2p5deg"]}
                   for r in result["profile_rows"]],
    }, indent=2))


if __name__ == "__main__":
    main()
