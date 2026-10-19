"""Calibrated pair stack for the HSC S19A re-Gaussianization shape catalog.

This is an analysis adapter, not a data downloader. It expects exact object IDs,
the S19A shape columns, and either full photo-z PDFs matched by object_id or
HSC Y3 tomographic source-bin assignments with calibrated n(z). Selection-bias
corrections follow the HSC Y3/S19A edge-density prescription and are computed
per radial bin. The returned DeltaSigma uses the weak-shear approximation; it is
not a reduced-shear mass fit or a final science measurement.
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np


C_KMS = 299792.458
G_MPC = 4.30091e-9  # Mpc (km/s)^2 / Msun
REQUIRED_PAIR_COLUMNS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "weight", "rms_e",
    "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_pass",
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


def tomographic_nz_lensing_efficiency(lens_redshift, source_redshift_bin,
                                     nz_redshift, nz_by_bin, cosmology):
    """Return pair-wise ``<D_ls/D_s>`` from calibrated tomographic ``n(z)``.

    ``source_redshift_bin`` uses HSC ``hsc_y3_zbin`` convention: 0 is not in
    the four Y3 source bins; bins 1..K map to rows of ``nz_by_bin``. Each row
    of ``nz_by_bin`` is a calibrated source-redshift distribution on the common
    ``nz_redshift`` grid. Bin 0 receives zero lensing efficiency and is
    reported separately by the caller.
    """
    zl = np.asarray(lens_redshift, dtype=float)
    raw_bins = np.asarray(source_redshift_bin)
    z = np.asarray(nz_redshift, dtype=float)
    nz = np.asarray(nz_by_bin, dtype=float)
    if zl.ndim != 1 or not np.all(np.isfinite(zl)) or np.any(zl < 0):
        raise ValueError("lens_redshift must be a finite non-negative vector")
    if raw_bins.ndim != 1 or len(raw_bins) != len(zl):
        raise ValueError("source_redshift_bin must have one value per lens-source pair")
    try:
        numeric_bins = raw_bins.astype(float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("source_redshift_bin must contain integer HSC Y3 bins") from None
    if not np.all(np.isfinite(numeric_bins)):
        raise ValueError("source_redshift_bin must contain finite integer HSC Y3 bins")
    bins = numeric_bins.astype(int)
    if not np.array_equal(numeric_bins, bins) or np.any(bins < 0):
        raise ValueError("source_redshift_bin must contain non-negative integer bins")
    if z.ndim != 1 or len(z) < 2 or not np.all(np.isfinite(z)) or np.any(z < 0) or np.any(np.diff(z) <= 0):
        raise ValueError("nz_redshift must be a strictly increasing, non-negative grid")
    if nz.ndim != 2 or nz.shape[1] != len(z) or nz.shape[0] < 1:
        raise ValueError("nz_by_bin must have shape [tomographic bin, redshift]")
    if not np.all(np.isfinite(nz)) or np.any(nz < 0):
        raise ValueError("nz_by_bin must be finite and non-negative")
    if np.any(bins > nz.shape[0]):
        raise ValueError("source_redshift_bin exceeds the supplied n(z) bins")
    norms = _trapz(nz, z, axis=1)
    if not np.all(np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("every calibrated n(z) must have positive finite integral")
    normalized = nz / norms[:, None]
    chi_grid = np.asarray(cosmology.comoving(z), dtype=float)
    chi_l = np.asarray(cosmology.comoving(zl), dtype=float)
    if not np.all(np.isfinite(chi_grid)) or not np.all(np.isfinite(chi_l)):
        raise ValueError("cosmology does not cover the lens and n(z) redshift ranges")
    beta = np.zeros(len(zl), dtype=float)
    chunk_size = 65536
    for start in range(0, len(zl), chunk_size):
        stop = min(start + chunk_size, len(zl))
        z_l = zl[start:stop, None]
        chi_l_chunk = chi_l[start:stop, None]
        for bin_index in range(1, nz.shape[0] + 1):
            selected = bins[start:stop] == bin_index
            if not np.any(selected):
                continue
            efficiency = np.divide(
                chi_grid[None, :] - chi_l_chunk[selected], chi_grid[None, :],
                out=np.zeros((int(np.count_nonzero(selected)), len(z)), dtype=float),
                where=(z[None, :] > z_l[selected]) & (chi_grid[None, :] > 0),
            )
            beta[start:stop][selected] = _trapz(
                normalized[bin_index - 1][None, :] * np.clip(efficiency, 0.0, 1.0),
                z, axis=1,
            )
    return beta


def load_hsc_y3_nz_sacc(path):
    """Load the four released HSC Y3 ``wl_0``..``wl_3`` n(z) tracers.

    The SACC file is an explicit local input; this function performs no network
    access. It returns the shared redshift grid and an ``[4, z]`` n(z) array in
    SACC tracer order. HSC ``hsc_y3_zbin`` values 1..4 map to these rows 0..3;
    bin 0 remains excluded by :func:`tomographic_nz_lensing_efficiency`.
    The tracer n(z) curves are fiducial distributions, not posterior samples,
    so their uncertainty still has to be propagated separately.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"HSC Y3 n(z) SACC file does not exist: {path}")
    try:
        import sacc
    except ImportError:
        # Astropy is already needed by the HSC image pipeline. The official
        # HSC release is a FITS SACC file, so read only its NZ tracer tables if
        # the higher-level SACC package is unavailable on the analysis host.
        try:
            from astropy.io import fits
        except ImportError as error:
            raise RuntimeError(
                "reading HSC Y3 n(z) requires either 'sacc' or 'astropy'") from error
        tracer_tables = {}
        with fits.open(path, memmap=True) as hdus:
            for hdu in hdus:
                header = hdu.header
                if (header.get("SACCTYPE") != "tracer"
                        or header.get("SACCCLSS") != "NZ"):
                    continue
                name = header.get("SACCNAME")
                if not name or hdu.data is None:
                    continue
                columns = {str(column).lower() for column in (hdu.data.names or [])}
                if not {"z", "nz"}.issubset(columns):
                    continue
                z = np.asarray(hdu.data["z"], dtype=float).copy()
                nz = np.asarray(hdu.data["nz"], dtype=float).copy()
                if name in tracer_tables:
                    old_z, old_nz = tracer_tables[name]
                    if not (np.array_equal(z, old_z) and np.array_equal(nz, old_nz)):
                        raise ValueError(f"SACC contains conflicting duplicate n(z) tracer {name}")
                else:
                    tracer_tables[name] = (z, nz)

        class FitsTracerView:
            def __getitem__(self, name):
                z, nz = tracer_tables[name]
                return type("NZTracer", (), {"z": z, "nz": nz})()

        dataset = type("SaccFitsView", (), {"tracers": FitsTracerView()})()
    else:
        dataset = sacc.Sacc.load_fits(str(path))
    redshift = None
    distributions = []
    for index in range(4):
        name = f"wl_{index}"
        try:
            tracer = dataset.tracers[name]
            z = np.asarray(tracer.z, dtype=float)
            nz = np.asarray(tracer.nz, dtype=float)
        except (KeyError, AttributeError, TypeError) as error:
            raise ValueError(f"SACC file is missing a valid n(z) tracer {name}") from error
        if z.ndim != 1 or len(z) < 2 or not np.all(np.isfinite(z)) or np.any(z < 0) or np.any(np.diff(z) <= 0):
            raise ValueError(f"SACC tracer {name} has an invalid redshift grid")
        if nz.shape != z.shape or not np.all(np.isfinite(nz)) or np.any(nz < 0):
            raise ValueError(f"SACC tracer {name} has an invalid n(z) distribution")
        norm = float(_trapz(nz, z))
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError(f"SACC tracer {name} has zero-integral n(z)")
        if redshift is None:
            redshift = z.copy()
        elif not np.array_equal(z, redshift):
            raise ValueError("HSC Y3 SACC n(z) tracers must share one redshift grid")
        distributions.append(nz.copy())
    return redshift, np.vstack(distributions)


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


def hsc_selection_bias(weight, mag_a10, resolution,
                       psf_shape11, psf_shape22, psf_shape12):
    """Estimate HSC edge-selection bias for one selected lens-source sample.

    This implements the published HSC Y3/S19A calibration prescription for the
    resolution (R2 > 0.3) and 1-arcsec aperture magnitude (mag_A < 25.5) cuts.
    ``weight`` should match the lensing statistic; this stack uses
    ``shape_weight * Sigma_crit^{-2}``. The PSF ellipticity is formed from the
    SDSS PSF second moments. Calibration coefficients and their quoted errors
    follow Li et al.'s released HSC calibration utility.
    """
    w = np.asarray(weight, dtype=float)
    mag = np.asarray(mag_a10, dtype=float)
    res = np.asarray(resolution, dtype=float)
    mxx = np.asarray(psf_shape11, dtype=float)
    myy = np.asarray(psf_shape22, dtype=float)
    mxy = np.asarray(psf_shape12, dtype=float)
    arrays = (w, mag, res, mxx, myy, mxy)
    if any(value.ndim != 1 for value in arrays) or len({len(value) for value in arrays}) != 1:
        raise ValueError("selection-bias inputs must be equal-length one-dimensional arrays")
    if len(w) == 0 or any(not np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("selection-bias inputs must be non-empty and finite")
    if np.any(w < 0) or np.sum(w) <= 0:
        raise ValueError("selection-bias weights must be non-negative with positive sum")
    trace = mxx + myy
    if np.any(trace <= 0):
        raise ValueError("PSF moment traces must be positive")
    e1_psf = (mxx - myy) / trace
    e2_psf = 2.0 * mxy / trace
    w_sum = float(np.sum(w))
    mag_width, res_width = 0.025, 0.01
    p_edge_mag = float(np.sum(w[mag >= 25.5 - mag_width]) / w_sum / mag_width)
    p_edge_res = float(np.sum(w[res <= 0.3 + res_width]) / w_sum / res_width)
    m_sel = -0.059 * p_edge_mag + 0.019 * p_edge_res
    a_sel = 0.0064 * p_edge_mag + 0.0063 * p_edge_res
    m_sel_err = math.sqrt((0.0089 * p_edge_mag) ** 2 + (0.0013 * p_edge_res) ** 2)
    a_sel_err = math.sqrt((0.0034 * p_edge_mag) ** 2 + (0.0009 * p_edge_res) ** 2)
    mean_psf_e1 = float(np.sum(w * e1_psf) / w_sum)
    mean_psf_e2 = float(np.sum(w * e2_psf) / w_sum)
    return {
        "m_sel": float(m_sel),
        "a_sel": float(a_sel),
        "c1_sel": float(a_sel * mean_psf_e1),
        "c2_sel": float(a_sel * mean_psf_e2),
        "m_sel_err": float(m_sel_err),
        "a_sel_err": float(a_sel_err),
        "c1_sel_err": float(a_sel_err * abs(mean_psf_e1)),
        "c2_sel_err": float(a_sel_err * abs(mean_psf_e2)),
        "p_edge_mag": p_edge_mag,
        "p_edge_resolution": p_edge_res,
        "mean_psf_e1": mean_psf_e1,
        "mean_psf_e2": mean_psf_e2,
    }


def decode_hsc_b_mode_mask(raw_mask, *, masked_value):
    """Decode a binary HSC ``b_mode_mask`` only with an explicit polarity.

    The HSC PDR3 page says this column identifies the region to remove, but does
    not document whether stored 1 or 0 denotes that region. Callers must verify
    the schema/data convention and pass that value as ``masked_value``; this
    function intentionally does not guess. Values outside {0, 1} are rejected.
    The returned boolean is ``True`` for rows to retain in the lensing stack.
    """
    values = np.asarray(raw_mask)
    if values.ndim != 1:
        raise ValueError("raw b_mode_mask must be a one-dimensional vector")
    if masked_value not in (0, 1, False, True):
        raise ValueError("masked_value must be the verified binary value 0 or 1")
    try:
        numeric = values.astype(float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("raw b_mode_mask must contain binary numeric values") from None
    if not np.all(np.isfinite(numeric)) or np.any((numeric != 0) & (numeric != 1)):
        raise ValueError("raw b_mode_mask is not binary; verify its catalog definition")
    return numeric != int(masked_value)


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
    if np.any(arrays["m"] <= -1):
        raise ValueError("per-source multiplicative bias must be greater than -1")
    flag = np.asarray(arrays["b_mode_pass"])
    if flag.dtype.kind not in "biu" or np.any((flag != 0) & (flag != 1)):
        raise ValueError("b_mode_pass must be a decoded boolean mask (0 or 1)")
    arrays["b_mode_pass"] = flag.astype(bool)
    return arrays, n


def stack_s19a_pairs(pairs, photoz_object_ids, photoz_redshift, photoz_pdf,
                     radial_edges_mpc, cosmology, *, apply_b_mode_mask=True,
                     _pair_lensing_efficiency=None):
    """Make HSC-calibrated tangential/cross and weak-approximation DeltaSigma bins.

    ``pairs`` contains lens-source rows. A source may appear with multiple lenses,
    but a (lens_id, object_id) pair must be unique. ``b_mode_pass`` is a decoded
    boolean indicating that the source passes the documented HSC B-mode region
    mask. ``apply_b_mode_mask=False`` is an explicit sensitivity mode: it keeps
    mask-failing rows while retaining their rejection count in the diagnostics.
    The edge-selection correction is computed per radial bin using the
    weighted number density near the HSC R2 and aperture-magnitude cuts, with
    pair weights ``shape_weight * Sigma_crit^{-2}``. The HSC weighted
    responsivity is recomputed in each radial bin.

    Returned shape-noise errors assume independent lens-source pairs. Repeated
    source IDs, field covariance, calibration uncertainty, random-lens subtraction,
    and reduced-shear corrections are not included; the product always has
    ``science_measurement=False``.
    """
    arr, n_input = _check_columns(pairs)
    edges = np.asarray(radial_edges_mpc, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError("radial_edges_mpc must be finite, positive, and strictly increasing")

    if _pair_lensing_efficiency is None:
        beta = pdf_lensing_efficiency(arr["lens_z"], arr["object_id"], photoz_object_ids,
                                      photoz_redshift, photoz_pdf, cosmology)
    else:
        beta = np.asarray(_pair_lensing_efficiency, dtype=float)
        if beta.shape != (n_input,) or not np.all(np.isfinite(beta)) or np.any((beta < 0) | (beta > 1)):
            raise ValueError("precomputed pair lensing efficiencies must be finite [0, 1] values")
    radius, phi = pair_geometry(arr, cosmology)
    bmode = arr["b_mode_pass"]
    ellipticity_amplitude = np.hypot(arr["e1"], arr["e2"])
    shape_selection = ((arr["sigma_e"] >= 0.0) & (arr["sigma_e"] <= 0.4)
                       & (ellipticity_amplitude < 2.0))
    source_selection = (shape_selection & (arr["mag_a10"] <= 25.5)
                        & (arr["resolution"] >= 0.3))
    mask_selection = bmode if apply_b_mode_mask else np.ones(len(bmode), dtype=bool)
    behind = beta > 0
    in_range = (radius >= edges[0]) & (radius <= edges[-1])
    usable = mask_selection & source_selection & behind & in_range
    counts = {
        "input_pairs": int(n_input),
        "hsc_source_selection_rejected": int(np.count_nonzero(~source_selection)),
        "hsc_shape_quality_rejected": int(np.count_nonzero(~shape_selection)),
        "ellipticity_amplitude_ge_2_rejected": int(np.count_nonzero(ellipticity_amplitude >= 2.0)),
        "sigma_e_outside_0_0p4_rejected": int(np.count_nonzero(
            (arr["sigma_e"] < 0.0) | (arr["sigma_e"] > 0.4))),
        "aperture_magnitude_over_25p5_rejected": int(np.count_nonzero(arr["mag_a10"] > 25.5)),
        "resolution_below_0p3_rejected": int(np.count_nonzero(arr["resolution"] < 0.3)),
        "b_mode_mask_applied": bool(apply_b_mode_mask),
        "b_mode_mask_rejected": int(np.count_nonzero(~bmode)),
        "b_mode_rejected": int(np.count_nonzero(~bmode)),
        "foreground_or_zero_efficiency": int(np.count_nonzero(mask_selection & source_selection & ~behind)),
        "outside_radial_range": int(np.count_nonzero(mask_selection & source_selection & behind & ~in_range)),
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
        if not np.any(in_bin):
            profile.append({
                "bin": index, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                "n_pairs": 0, "n_lenses": 0, "n_sources": 0,
                "responsivity": None, "mean_gt": None, "mean_gx": None,
                "delta_sigma_weak_msun_mpc2": None,
                "delta_sigma_shape_noise_msun_mpc2": None,
                "cross_delta_sigma_weak_msun_mpc2": None,
                "cross_shape_noise_msun_mpc2": None,
                "selection_bias": None,
            })
            continue
        w = arr["weight"][in_bin]
        beta_k = sigma_crit_inv[in_bin]
        lensing_weight = w * beta_k * beta_k
        rms = arr["rms_e"][in_bin]
        responsivity = float(1.0 - np.sum(lensing_weight * rms * rms) / np.sum(lensing_weight))
        if not np.isfinite(responsivity) or responsivity <= 0:
            raise ValueError(f"radial bin {index} has non-positive HSC shear responsivity")
        m = arr["m"][in_bin]
        selection_weight = lensing_weight
        selection_bin = hsc_selection_bias(
            selection_weight,
            arr["mag_a10"][in_bin],
            arr["resolution"][in_bin],
            arr["psf_shape11"][in_bin],
            arr["psf_shape22"][in_bin],
            arr["psf_shape12"][in_bin],
        )
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
            "selection_bias": selection_bin,
        })
    return {
        "schema": "hsc-s19a-shape-stack-v1",
        "status": "exploratory",
        "science_measurement": False,
        "counts": counts,
        "profile": profile,
        "calibration": {
            "shape_catalog": "S19A PDR3 re-Gaussianization",
            "responsivity": "per-radial-bin HSC responsivity from derived_rms_e weighted by shape_weight*SigmaCritInv^2",
            "per_object_bias": "derived_shear_bias_m/c1/c2",
            "selection_bias_method": "HSC edge-density prescription for R2>0.3 and magA10<25.5; weights=shape_weight*SigmaCritInv^2",
            "source_selection": "HSC Y3 cuts |e|<2, 0<=sigma_e<=0.4, R2>=0.3, and magA10<=25.5 applied before stacking",
            "photoz_geometry": "full P(z) integrated to mean SigmaCritInv for each lens-source pair",
            "photoz_bias_calibration": "not applied; PDF bias/dilution factor needs external calibration sample",
            "b_mode_mask": "input b_mode_pass must be decoded with verified stored-value polarity; decode_hsc_b_mode_mask requires masked_value",
            "delta_sigma": "weak-shear approximation; reduced shear is not modeled",
            "uncertainty": "independent-pair shape noise only; not a field covariance",
            "null_tests": "random-lens subtraction and spatial block covariance not included",
        },
    }


def stack_s19a_tomographic_pairs(pairs, nz_redshift, nz_by_bin,
                                  radial_edges_mpc, cosmology,
                                  *, apply_b_mode_mask=True):
    """Stack using HSC Y3 calibrated tomographic ``n(z)`` distributions.

    This path consumes ``source_zbin`` from the minimal shape-row export and
    does not require per-object full-PDF rows. The calibrated ``n(z)`` must
    match the HSC Y3 shape-source selection. HSC zbin 0 is excluded; bins 1..K
    use corresponding rows in ``nz_by_bin``. Redshift-distribution uncertainty
    still needs propagation before a science result.
    """
    if "source_zbin" not in pairs:
        raise ValueError("tomographic n(z) stacking requires source_zbin")
    arr, n_input = _check_columns(pairs)
    raw_bins = np.asarray(pairs["source_zbin"])
    try:
        numeric_bins = raw_bins.astype(float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("source_zbin must contain integer HSC Y3 bins") from None
    if not np.all(np.isfinite(numeric_bins)):
        raise ValueError("source_zbin must contain finite integer HSC Y3 bins")
    source_bins = numeric_bins.astype(int)
    if not np.array_equal(numeric_bins, source_bins):
        raise ValueError("source_zbin must contain integer HSC Y3 bins")
    source_bin_by_id = {}
    for object_id, source_bin in zip(arr["object_id"], source_bins):
        previous = source_bin_by_id.setdefault(str(object_id), int(source_bin))
        if previous != int(source_bin):
            raise ValueError("a repeated object_id has inconsistent source_zbin assignments")
    beta = tomographic_nz_lensing_efficiency(
        arr["lens_z"], source_bins, nz_redshift, nz_by_bin, cosmology,
    )
    result = stack_s19a_pairs(
        pairs, None, None, None, radial_edges_mpc, cosmology,
        apply_b_mode_mask=apply_b_mode_mask,
        _pair_lensing_efficiency=beta,
    )
    result["counts"]["unassigned_source_zbin"] = int(np.count_nonzero(source_bins == 0))
    result["calibration"]["photoz_geometry"] = (
        "HSC Y3 tomographic n(z) integrated to pair-wise mean SigmaCritInv"
    )
    result["calibration"]["photoz_bias_calibration"] = (
        "calibrated n(z) applied; posterior/systematic uncertainty not propagated"
    )
    result["calibration"]["photoz_source_bins"] = (
        "hsc_y3_zbin 0 excluded; 1..K mapped to the supplied calibrated n(z) rows"
    )
    return result


def subtract_random_lens_profile(lens_result, random_result):
    """Subtract a matched random-lens profile from the cluster-lens profile.

    This returns point estimates only. It intentionally does not combine the
    two independent-pair shape-noise errors: overlapping sources and common
    survey modes require a paired spatial jackknife or equivalent covariance.
    """
    for name, result in (("lens_result", lens_result), ("random_result", random_result)):
        if result.get("schema") != "hsc-s19a-shape-stack-v1":
            raise ValueError(f"{name} has an unexpected S19A stack schema")
        if result.get("science_measurement") is not False:
            raise ValueError(f"{name} has an unexpected science gate")
    lens_profile, random_profile = lens_result.get("profile"), random_result.get("profile")
    if not isinstance(lens_profile, list) or not isinstance(random_profile, list):
        raise ValueError("both inputs must contain radial profile lists")
    if len(lens_profile) != len(random_profile):
        raise ValueError("lens and random profiles have different radial-bin counts")
    combined = []
    for lens_row, random_row in zip(lens_profile, random_profile):
        for key in ("bin", "r_min_mpc", "r_max_mpc"):
            if lens_row.get(key) != random_row.get(key):
                raise ValueError("lens and random radial-bin definitions differ")
        row = {key: lens_row[key] for key in ("bin", "r_min_mpc", "r_max_mpc")}
        for signal, output in (("delta_sigma_weak_msun_mpc2", "delta_sigma_random_subtracted_msun_mpc2"),
                               ("cross_delta_sigma_weak_msun_mpc2", "cross_delta_sigma_random_subtracted_msun_mpc2")):
            lens_value, random_value = lens_row.get(signal), random_row.get(signal)
            row[output] = (float(lens_value) - float(random_value)
                           if lens_value is not None and random_value is not None else None)
        row["n_lens_pairs"] = int(lens_row.get("n_pairs", 0))
        row["n_random_pairs"] = int(random_row.get("n_pairs", 0))
        row["covariance_status"] = "not_estimated; paired spatial jackknife required"
        combined.append(row)
    return {
        "schema": "hsc-s19a-random-subtracted-profile-v1",
        "status": "exploratory",
        "science_measurement": False,
        "method": "cluster-lens minus matched-random-lens DeltaSigma and cross profile",
        "photoz_dilution_correction_applied": False,
        "covariance": "not estimated; requires matched spatial-block jackknife including source overlap",
        "profile": combined,
    }


def equal_area_source_blocks(source_ra_deg, source_dec_deg, *, n_ra=4, n_sin_dec=4):
    """Assign source positions to deterministic equal-area lon/sin(lat) cells.

    This is a reproducible fallback when survey-region labels are unavailable.
    It is not a substitute for validated contiguous survey jackknife regions.
    Identical sky positions always receive the same block, including across
    cluster-lens and random-lens pair tables.
    """
    ra = np.asarray(source_ra_deg, dtype=float)
    dec = np.asarray(source_dec_deg, dtype=float)
    if (ra.ndim != 1 or dec.shape != ra.shape or not np.all(np.isfinite(ra))
            or not np.all(np.isfinite(dec)) or np.any((ra < 0) | (ra >= 360))
            or np.any(np.abs(dec) > 90)):
        raise ValueError("source positions must be finite RA/Dec vectors in degrees")
    if n_ra < 1 or n_sin_dec < 1:
        raise ValueError("equal-area grid dimensions must be positive")
    ra_bin = np.minimum((ra / 360.0 * n_ra).astype(int), n_ra - 1)
    sin_dec = np.sin(np.deg2rad(dec))
    dec_bin = np.minimum(((sin_dec + 1.0) * 0.5 * n_sin_dec).astype(int), n_sin_dec - 1)
    return dec_bin * n_ra + ra_bin


def paired_tomographic_spatial_jackknife(
        lens_pairs, random_pairs, nz_redshift, nz_by_bin, radial_edges_mpc,
        cosmology, *, n_ra_blocks=4, n_sin_dec_blocks=4, apply_b_mode_mask=True):
    """Recompute both stacks after deleting the same source-sky block.

    Blocks are defined by equal-area RA/sin(Dec) cells. Every leave-one-block
    replicate reruns the full tomographic estimator, including radial-bin
    responsivity and edge-selection calibration, before random subtraction.
    The returned covariance covers both DeltaSigma and its cross component.
    """
    lens, n_lens = _check_columns(lens_pairs)
    random, n_random = _check_columns(random_pairs)
    if n_lens == 0 or n_random == 0:
        raise ValueError("lens and random pair tables must both be non-empty")
    lens_blocks = equal_area_source_blocks(
        lens["source_ra_deg"], lens["source_dec_deg"],
        n_ra=n_ra_blocks, n_sin_dec=n_sin_dec_blocks,
    )
    random_blocks = equal_area_source_blocks(
        random["source_ra_deg"], random["source_dec_deg"],
        n_ra=n_ra_blocks, n_sin_dec=n_sin_dec_blocks,
    )
    block_ids = np.unique(np.concatenate((lens_blocks, random_blocks)))
    if len(block_ids) < 3:
        raise ValueError("paired spatial jackknife requires at least three occupied sky blocks")

    def subset(pairs, keep):
        return {key: np.asarray(value)[keep] for key, value in pairs.items()}

    def profile_vector(pairs):
        result = stack_s19a_tomographic_pairs(
            pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
            apply_b_mode_mask=apply_b_mode_mask,
        )
        values = []
        for row in result["profile"]:
            ds = row["delta_sigma_weak_msun_mpc2"]
            cross = row["cross_delta_sigma_weak_msun_mpc2"]
            if ds is None or cross is None or not np.isfinite(ds) or not np.isfinite(cross):
                raise ValueError("every radial bin must remain populated in every jackknife replicate")
            values.extend((float(ds), float(cross)))
        return result, np.asarray(values, dtype=float)

    lens_full_result, lens_full = profile_vector(lens)
    random_full_result, random_full = profile_vector(random)
    lens_replicas = []
    random_replicas = []
    block_counts = []
    for block in block_ids:
        keep_lens = lens_blocks != block
        keep_random = random_blocks != block
        _, lens_replica = profile_vector(subset(lens, keep_lens))
        _, random_replica = profile_vector(subset(random, keep_random))
        lens_replicas.append(lens_replica)
        random_replicas.append(random_replica)
        block_counts.append({
            "block_id": int(block),
            "lens_pairs_deleted": int(np.count_nonzero(~keep_lens)),
            "random_pairs_deleted": int(np.count_nonzero(~keep_random)),
        })
    covariance = paired_jackknife_covariance(
        np.asarray(lens_replicas), np.asarray(random_replicas),
    )
    full_delta = lens_full - random_full
    errors = covariance["standard_error"]
    profile = []
    for index, (lens_row, random_row) in enumerate(zip(
            lens_full_result["profile"], random_full_result["profile"])):
        profile.append({
            "bin": index,
            "r_min_mpc": lens_row["r_min_mpc"],
            "r_max_mpc": lens_row["r_max_mpc"],
            "delta_sigma_random_subtracted_msun_mpc2": float(full_delta[2 * index]),
            "cross_delta_sigma_random_subtracted_msun_mpc2": float(full_delta[2 * index + 1]),
            "delta_sigma_jackknife_error_msun_mpc2": float(errors[2 * index]),
            "cross_jackknife_error_msun_mpc2": float(errors[2 * index + 1]),
            "n_lens_pairs": lens_row["n_pairs"],
            "n_random_pairs": random_row["n_pairs"],
        })
    return {
        "schema": "hsc-s19a-paired-spatial-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "block_method": "equal-area RA/sin(Dec) grid assigned from source coordinates",
        "block_grid": {"n_ra": int(n_ra_blocks), "n_sin_dec": int(n_sin_dec_blocks)},
        "n_occupied_blocks": int(len(block_ids)),
        "block_deletions": block_counts,
        "profile_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "profile": profile,
        "covariance": covariance["covariance"].tolist(),
        "jackknife_mean_diagnostic": covariance["mean_random_subtracted"].tolist(),
        "limitations": [
            "equal-area grid cells are a fallback and may combine disconnected survey regions",
            "survey-region jackknife, photo-z calibration uncertainty, B-mode polarity, and reduced shear remain independent gates",
            "this covariance does not by itself authorize a science detection",
        ],
    }


def paired_jackknife_covariance(lens_leave_one_block, random_leave_one_block):
    """Covariance of random-subtracted profiles from matched delete-one blocks.

    Inputs are ``[block, radial_bin]`` arrays produced by rerunning the full
    calibrated estimator after deleting the same sky block from the cluster and
    random pair samples. Shared source IDs must be assigned to the same sky
    block in both samples before these arrays are produced.
    """
    lens = np.asarray(lens_leave_one_block, dtype=float)
    random = np.asarray(random_leave_one_block, dtype=float)
    if lens.ndim != 2 or random.shape != lens.shape or lens.shape[0] < 3 or lens.shape[1] < 1:
        raise ValueError("matched lens/random jackknife arrays must share shape [>=3 blocks, >=1 bins]")
    if not np.all(np.isfinite(lens)) or not np.all(np.isfinite(random)):
        raise ValueError("jackknife profile values must be finite")
    delta = lens - random
    mean = np.mean(delta, axis=0)
    centered = delta - mean
    covariance = (lens.shape[0] - 1.0) / lens.shape[0] * centered.T @ centered
    return {
        "n_blocks": int(lens.shape[0]),
        "mean_random_subtracted": mean,
        "covariance": covariance,
        "standard_error": np.sqrt(np.diag(covariance)),
    }
