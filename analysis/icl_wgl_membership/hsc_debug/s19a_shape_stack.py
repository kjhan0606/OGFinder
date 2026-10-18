"""Calibrated pair stack for the HSC S19A re-Gaussianization shape catalog.

This is an analysis adapter, not a data downloader.  It expects exact object IDs,
the S19A shape columns, and the full photo-z PDF table matched by object_id.
Selection-bias corrections are a required, caller-supplied calibration per radial
bin because they depend on the selected source sample.  The returned DeltaSigma
uses the weak-shear approximation; it is not a reduced-shear mass fit or a final
science measurement.
"""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np


C_KMS = 299792.458
G_MPC = 4.30091e-9  # Mpc (km/s)^2 / Msun
REQUIRED_PAIR_COLUMNS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "weight", "rms_e",
    "m", "c1", "c2", "b_mode_pass",
)


def _trapz(y, x, axis=-1):
    fn = getattr(np, "trapezoid", None) or np.trapz
    return fn(y, x, axis=axis)


def _ids(values, name):
    """Keep catalog IDs exact; float IDs may already have lost int64 precision."""
    result = []
    for value in values:
        if isinstance(value, (float, np.floating)):
            raise ValueError(f"{name} contains a floating-point ID; preserve IDs as integers or strings")
        if isinstance(value, (int, np.integer)):
            result.append(str(int(value)))
        elif isinstance(value, str) and value:
            result.append(value)
        else:
            raise ValueError(f"{name} contains an empty or unsupported ID")
    return np.asarray(result, dtype=object)


def pdf_lensing_efficiency(lens_redshift, pair_object_ids, photoz_object_ids,
                           photoz_redshift, photoz_pdf, cosmology):
    """Return pair-wise ``<D_ls / D_s>`` from full source-redshift PDFs.

    ``photoz_pdf`` has one row per unique ``photoz_object_ids`` entry and a
    shared redshift grid. PDFs are normalized here. The cosmology object must
    provide ``comoving(z)`` in Mpc and must cover the full redshift grid.
    Sources at or in front of each lens contribute zero efficiency.
    """
    zl = np.asarray(lens_redshift, dtype=float)
    pair_ids = _ids(pair_object_ids, "pair_object_ids")
    pz_ids = _ids(photoz_object_ids, "photoz_object_ids")
    z = np.asarray(photoz_redshift, dtype=float)
    pdf = np.asarray(photoz_pdf, dtype=float)
    if zl.ndim != 1 or len(zl) != len(pair_ids) or not np.all(np.isfinite(zl)) or np.any(zl < 0):
        raise ValueError("lens_redshift must be a finite non-negative value per pair")
    if z.ndim != 1 or len(z) < 2 or not np.all(np.isfinite(z)) or np.any(z < 0) or np.any(np.diff(z) <= 0):
        raise ValueError("photoz_redshift must be a strictly increasing, non-negative grid")
    if pdf.shape != (len(pz_ids), len(z)) or not np.all(np.isfinite(pdf)) or np.any(pdf < 0):
        raise ValueError("photoz_pdf must be a finite non-negative [object_id, redshift] array")
    if len(set(pz_ids.tolist())) != len(pz_ids):
        raise ValueError("photoz_object_ids contains duplicates")
    lookup = {object_id: row for row, object_id in enumerate(pz_ids)}
    try:
        pdf_rows = np.fromiter((lookup[object_id] for object_id in pair_ids), dtype=np.int64, count=len(pair_ids))
    except KeyError as error:
        raise ValueError(f"Missing full photo-z PDF for object_id {error.args[0]}") from None
    norms = _trapz(pdf, z, axis=1)
    if not np.all(np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("every source photo-z PDF must have positive finite integral")
    pdf = pdf / norms[:, None]
    chi_grid = np.asarray(cosmology.comoving(z), dtype=float)
    chi_l = np.asarray(cosmology.comoving(zl), dtype=float)
    if not np.all(np.isfinite(chi_grid)) or not np.all(np.isfinite(chi_l)):
        raise ValueError("cosmology does not cover the lens and source redshift ranges")
    beta = np.zeros(len(zl), dtype=float)
    # Chunk to bound memory for a large lens-source pair table.
    chunk_size = 65536
    for start in range(0, len(zl), chunk_size):
        stop = min(start + chunk_size, len(zl))
        z_l = zl[start:stop, None]
        chi_l_chunk = chi_l[start:stop, None]
        chi_s = chi_grid[None, :]
        efficiency = np.divide(
            chi_s - chi_l_chunk, chi_s,
            out=np.zeros((stop - start, len(z)), dtype=float),
            where=(z[None, :] > z_l) & (chi_s > 0),
        )
        efficiency = np.clip(efficiency, 0.0, 1.0)
        beta[start:stop] = _trapz(pdf[pdf_rows[start:stop]] * efficiency, z, axis=1)
    return beta


def pair_geometry(pairs, cosmology):
    """Return proper projected radius (Mpc) and position angle (east to north)."""
    ra_l = np.deg2rad(np.asarray(pairs["lens_ra_deg"], dtype=float))
    dec_l = np.deg2rad(np.asarray(pairs["lens_dec_deg"], dtype=float))
    ra_s = np.deg2rad(np.asarray(pairs["source_ra_deg"], dtype=float))
    dec_s = np.deg2rad(np.asarray(pairs["source_dec_deg"], dtype=float))
    zl = np.asarray(pairs["lens_z"], dtype=float)
    dra = (ra_s - ra_l + math.pi) % (2.0 * math.pi) - math.pi
    sin_l, cos_l = np.sin(dec_l), np.cos(dec_l)
    sin_s, cos_s = np.sin(dec_s), np.cos(dec_s)
    denominator = sin_l * sin_s + cos_l * cos_s * np.cos(dra)
    if np.any(denominator <= 0):
        raise ValueError("source lies outside the lens-centred tangent hemisphere")
    x_east = cos_s * np.sin(dra) / denominator
    y_north = (cos_l * sin_s - sin_l * cos_s * np.cos(dra)) / denominator
    theta = np.arctan(np.hypot(x_east, y_north))
    phi = np.arctan2(y_north, x_east)
    d_l = np.asarray(cosmology.angular_diameter(zl), dtype=float)
    if not np.all(np.isfinite(d_l)) or np.any(d_l <= 0):
        raise ValueError("lens angular-diameter distances must be positive and finite")
    return theta * d_l, phi


def _to_tangential(phi, e1, e2):
    c2, s2 = np.cos(2.0 * phi), np.sin(2.0 * phi)
    e_t = -e1 * c2 - e2 * s2
    e_x = e1 * s2 - e2 * c2
    return e_t, e_x


def _check_columns(pairs):
    missing = [key for key in REQUIRED_PAIR_COLUMNS if key not in pairs]
    if missing:
        raise ValueError("missing S19A pair columns: " + ", ".join(missing))
    arrays = {key: np.asarray(value) for key, value in pairs.items()}
    lengths = {len(value) for value in arrays.values() if np.ndim(value) > 0}
    if len(lengths) != 1:
        raise ValueError("pair columns must be one-dimensional arrays with equal lengths")
    n = next(iter(lengths), 0)
    if any(np.asarray(value).ndim != 1 for value in arrays.values()):
        raise ValueError("pair columns must be one-dimensional")
    arrays["lens_id"] = _ids(arrays["lens_id"], "lens_id")
    arrays["object_id"] = _ids(arrays["object_id"], "object_id")
    if len(set(zip(arrays["lens_id"].tolist(), arrays["object_id"].tolist()))) != n:
        raise ValueError("duplicate (lens_id, object_id) pairs would double count sources")
    numeric = [key for key in REQUIRED_PAIR_COLUMNS if key not in ("lens_id", "object_id", "b_mode_pass")]
    for key in numeric:
        arrays[key] = np.asarray(arrays[key], dtype=float)
        if not np.all(np.isfinite(arrays[key])):
            raise ValueError(f"{key} contains non-finite values")
    if np.any(arrays["lens_z"] < 0) or np.any(np.abs(arrays["lens_dec_deg"]) > 90) or np.any(np.abs(arrays["source_dec_deg"]) > 90):
        raise ValueError("invalid lens/source sky coordinates or redshifts")
    if np.any((arrays["lens_ra_deg"] < 0) | (arrays["lens_ra_deg"] >= 360)) or np.any((arrays["source_ra_deg"] < 0) | (arrays["source_ra_deg"] >= 360)):
        raise ValueError("RA must be in [0, 360) degrees")
    if np.any(arrays["weight"] <= 0) or np.any(arrays["rms_e"] < 0) or np.any(arrays["rms_e"] >= 1):
        raise ValueError("S19A weight must be positive and rms_e must be in [0, 1)")
    if np.any(arrays["e1"] ** 2 + arrays["e2"] ** 2 >= 1.0 + 1e-10):
        raise ValueError("ellipticity magnitude must be smaller than one")
    if np.any(arrays["m"] <= -1):
        raise ValueError("per-source multiplicative bias must be greater than -1")
    flag = np.asarray(arrays["b_mode_pass"])
    if flag.dtype.kind not in "biu" or np.any((flag != 0) & (flag != 1)):
        raise ValueError("b_mode_pass must be a decoded boolean mask (0 or 1)")
    arrays["b_mode_pass"] = flag.astype(bool)
    return arrays, n


def stack_s19a_pairs(pairs, photoz_object_ids, photoz_redshift, photoz_pdf,
                     radial_edges_mpc, cosmology, selection_corrections,
                     calibration_reference):
    """Make HSC-calibrated tangential/cross and weak-approximation DeltaSigma bins.

    ``pairs`` contains lens-source rows. A source may appear with multiple lenses,
    but a (lens_id, object_id) pair must be unique. ``b_mode_pass`` is a decoded
    boolean indicating that the source passes the documented HSC B-mode region
    mask. ``selection_corrections`` is one mapping per radial bin with ``m_sel``,
    ``c1_sel``, and ``c2_sel`` in shear units, calibrated for that exact selected
    ensemble. The HSC weighted responsivity is recomputed in each radial bin.

    Returned shape-noise errors assume independent lens-source pairs. Repeated
    source IDs, field covariance, calibration uncertainty, random-lens subtraction,
    and reduced-shear corrections are not included; the product always has
    ``science_measurement=False``.
    """
    arr, n_input = _check_columns(pairs)
    edges = np.asarray(radial_edges_mpc, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError("radial_edges_mpc must be finite, positive, and strictly increasing")
    if len(selection_corrections) != len(edges) - 1:
        raise ValueError("one selection-bias calibration is required per radial bin")
    if not isinstance(calibration_reference, str) or not calibration_reference.strip():
        raise ValueError("calibration_reference must identify the source of selection-bias corrections")
    selection = []
    for index, correction in enumerate(selection_corrections):
        required = ("m_sel", "c1_sel", "c2_sel")
        if not isinstance(correction, Mapping) or any(key not in correction for key in required):
            raise ValueError(f"selection correction {index} needs m_sel, c1_sel, and c2_sel")
        item = {key: float(correction[key]) for key in required}
        if not all(np.isfinite(value) for value in item.values()) or item["m_sel"] <= -1:
            raise ValueError(f"selection correction {index} is invalid")
        selection.append(item)

    beta = pdf_lensing_efficiency(arr["lens_z"], arr["object_id"], photoz_object_ids,
                                  photoz_redshift, photoz_pdf, cosmology)
    radius, phi = pair_geometry(arr, cosmology)
    bmode = arr["b_mode_pass"]
    behind = beta > 0
    in_range = (radius >= edges[0]) & (radius <= edges[-1])
    usable = bmode & behind & in_range
    counts = {
        "input_pairs": int(n_input),
        "b_mode_rejected": int(np.count_nonzero(~bmode)),
        "foreground_or_zero_efficiency": int(np.count_nonzero(bmode & ~behind)),
        "outside_radial_range": int(np.count_nonzero(bmode & behind & ~in_range)),
        "used_pairs": int(np.count_nonzero(usable)),
        "unique_used_lenses": int(len(set(arr["lens_id"][usable].tolist()))),
        "unique_used_sources": int(len(set(arr["object_id"][usable].tolist()))),
    }
    chi_factor = 4.0 * math.pi * G_MPC / (C_KMS * C_KMS)
    d_l = np.asarray(cosmology.angular_diameter(arr["lens_z"]), dtype=float)
    sigma_crit_inv = chi_factor * d_l * beta  # Mpc^2 / Msun
    e_t, e_x = _to_tangential(phi, arr["e1"], arr["e2"])
    c_t, c_x = _to_tangential(phi, arr["c1"], arr["c2"])
    profile = []
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        in_bin = usable & (radius >= lo) & ((radius < hi) if index < len(edges) - 2 else (radius <= hi))
        selection_bin = selection[index]
        if not np.any(in_bin):
            profile.append({
                "bin": index, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                "n_pairs": 0, "n_lenses": 0, "n_sources": 0,
                "responsivity": None, "mean_gt": None, "mean_gx": None,
                "delta_sigma_weak_msun_mpc2": None,
                "delta_sigma_shape_noise_msun_mpc2": None,
                "cross_delta_sigma_weak_msun_mpc2": None,
                "cross_shape_noise_msun_mpc2": None,
            })
            continue
        w = arr["weight"][in_bin]
        rms = arr["rms_e"][in_bin]
        responsivity = float(1.0 - np.sum(w * rms * rms) / np.sum(w))
        if not np.isfinite(responsivity) or responsivity <= 0:
            raise ValueError(f"radial bin {index} has non-positive HSC shear responsivity")
        m = arr["m"][in_bin]
        msel = selection_bin["m_sel"]
        if np.any((1.0 + m) * (1.0 + msel) <= 0):
            raise ValueError(f"radial bin {index} has a non-positive calibrated response")
        gt = (e_t[in_bin] / (2.0 * responsivity) - c_t[in_bin]) / (1.0 + m)
        gx = (e_x[in_bin] / (2.0 * responsivity) - c_x[in_bin]) / (1.0 + m)
        csel_t, csel_x = _to_tangential(
            phi[in_bin], np.full(np.count_nonzero(in_bin), selection_bin["c1_sel"]),
            np.full(np.count_nonzero(in_bin), selection_bin["c2_sel"]),
        )
        gt = (gt - csel_t) / (1.0 + msel)
        gx = (gx - csel_x) / (1.0 + msel)
        # The catalog weight is inverse variance for distortion e; propagate the
        # responsivity and multiplicative corrections into the shear weight.
        w_shear = w * (2.0 * responsivity * (1.0 + m) * (1.0 + msel)) ** 2
        mean_weight = float(np.sum(w_shear))
        mean_gt = float(np.sum(w_shear * gt) / mean_weight)
        mean_gx = float(np.sum(w_shear * gx) / mean_weight)
        beta_k = sigma_crit_inv[in_bin]
        denominator = float(np.sum(w_shear * beta_k * beta_k))
        if not np.isfinite(denominator) or denominator <= 0:
            raise ValueError(f"radial bin {index} has no positive lensing efficiency")
        ds = float(np.sum(w_shear * beta_k * gt) / denominator)
        ds_x = float(np.sum(w_shear * beta_k * gx) / denominator)
        sigma_ds = float(math.sqrt(1.0 / denominator))
        profile.append({
            "bin": index, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
            "n_pairs": int(np.count_nonzero(in_bin)),
            "n_lenses": int(len(set(arr["lens_id"][in_bin].tolist()))),
            "n_sources": int(len(set(arr["object_id"][in_bin].tolist()))),
            "responsivity": responsivity, "mean_gt": mean_gt, "mean_gx": mean_gx,
            "mean_gt_shape_noise": float(math.sqrt(1.0 / mean_weight)),
            "mean_gx_shape_noise": float(math.sqrt(1.0 / mean_weight)),
            "delta_sigma_weak_msun_mpc2": ds,
            "delta_sigma_shape_noise_msun_mpc2": sigma_ds,
            "cross_delta_sigma_weak_msun_mpc2": ds_x,
            "cross_shape_noise_msun_mpc2": sigma_ds,
            "effective_pairs": float(np.sum(w_shear) ** 2 / np.sum(w_shear * w_shear)),
        })
    return {
        "schema": "hsc-s19a-shape-stack-v1",
        "status": "exploratory",
        "science_measurement": False,
        "counts": counts,
        "profile": profile,
        "calibration": {
            "shape_catalog": "S19A PDR3 re-Gaussianization",
            "responsivity": "per-radial-bin weighted HSC responsivity from derived_rms_e",
            "per_object_bias": "derived_shear_bias_m/c1/c2",
            "selection_bias_reference": calibration_reference,
            "photoz_dilution": "full P(z) integrated as <D_ls/D_s> for each lens-source pair",
            "b_mode_mask": "input b_mode_pass must be decoded from the HSC b_mode_mask field",
            "delta_sigma": "weak-shear approximation; reduced shear is not modeled",
            "uncertainty": "independent-pair shape noise only; not a field covariance",
            "null_tests": "random-lens subtraction and spatial block covariance not included",
        },
    }
