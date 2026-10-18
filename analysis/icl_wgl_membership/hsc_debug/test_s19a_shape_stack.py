import math

import numpy as np
import pytest

from ogfkit import wlcluster as wl
from s19a_shape_stack import (
    G_MPC,
    C_KMS,
    _to_tangential,
    hsc_selection_bias,
    pair_geometry,
    pdf_lensing_efficiency,
    stack_s19a_pairs,
)


@pytest.fixture(scope="module")
def cosmo():
    return wl.Cosmology(H0=70.0, Om=0.3)


def test_full_pdf_efficiency_normalizes_and_drops_foreground(cosmo):
    ids = np.array(["source-a", "source-b"], dtype=object)
    z = np.linspace(0.0, 2.0, 401)
    pdf = np.vstack([
        np.exp(-0.5 * ((z - 1.0) / 0.12) ** 2),
        np.where(z < 0.3, np.exp(-0.5 * ((z - 0.2) / 0.04) ** 2), 0.0),
    ])
    beta = pdf_lensing_efficiency(
        [0.3, 0.3], ids, ids, z, pdf * np.array([[3.0], [0.25]]), cosmo,
    )
    assert 0.0 < beta[0] < 1.0
    assert beta[1] == 0.0
    assert beta[0] == pytest.approx(
        pdf_lensing_efficiency([0.3], ["source-a"], ["source-a"], z, pdf[:1], cosmo)[0]
    )


def test_pair_geometry_and_tangential_sign_follow_hsc_sky_axes(cosmo):
    pairs = {
        "lens_ra_deg": np.array([150.0, 150.0]),
        "lens_dec_deg": np.array([0.0, 0.0]),
        "lens_z": np.array([0.3, 0.3]),
        "source_ra_deg": np.array([150.01, 150.0]),
        "source_dec_deg": np.array([0.0, 0.01]),
    }
    radius, phi = pair_geometry(pairs, cosmo)
    assert radius[0] > 0 and radius[1] > 0
    assert phi[0] == pytest.approx(0.0, abs=1e-7)  # east is +x
    assert phi[1] == pytest.approx(math.pi / 2.0, abs=1e-7)  # north is +y
    gt, gx = _to_tangential(phi, np.array([-0.1, 0.1]), np.zeros(2))
    assert np.allclose(gt, 0.1)
    assert np.allclose(gx, 0.0, atol=1e-12)


def test_hsc_edge_selection_correction_uses_pair_weights_and_psf_ellipticity():
    weight = np.ones(4)
    mag_a10 = np.array([25.49, 25.48, 24.0, 24.5])
    resolution = np.array([0.305, 0.31, 0.4, 0.5])
    e1_psf = np.array([0.02, 0.04, 0.06, 0.08])
    e2_psf = np.array([-0.04, -0.02, 0.02, 0.04])
    result = hsc_selection_bias(
        weight, mag_a10, resolution,
        1.0 + e1_psf, 1.0 - e1_psf, e2_psf,
    )
    p_mag = 2.0 / 4.0 / 0.025
    p_res = 2.0 / 4.0 / 0.01
    expected_a = 0.0064 * p_mag + 0.0063 * p_res
    assert result["p_edge_mag"] == pytest.approx(p_mag)
    assert result["p_edge_resolution"] == pytest.approx(p_res)
    assert result["m_sel"] == pytest.approx(-0.059 * p_mag + 0.019 * p_res)
    assert result["a_sel"] == pytest.approx(expected_a)
    assert result["c1_sel"] == pytest.approx(expected_a * np.mean(e1_psf))
    assert result["c2_sel"] == pytest.approx(expected_a * np.mean(e2_psf))


def _synthetic_catalog(cosmo):
    rng = np.random.default_rng(5931)
    z_lens = 0.3
    ra_lens, dec_lens = 150.0, 2.0
    radial_edges = np.array([0.4, 0.7, 1.1, 1.5, 2.1])
    radius = np.concatenate([
        rng.uniform(lo + 0.03, hi - 0.03, 100)
        for lo, hi in zip(radial_edges[:-1], radial_edges[1:])
    ])
    phi = rng.uniform(-math.pi, math.pi, len(radius))
    d_l = float(cosmo.angular_diameter(z_lens))
    theta = radius / d_l
    ra_source = (ra_lens + np.rad2deg(theta * np.cos(phi) / np.cos(np.deg2rad(dec_lens)))) % 360.0
    dec_source = dec_lens + np.rad2deg(theta * np.sin(phi))
    n = len(radius)
    lens_ids = np.full(n, "lens-1", dtype=object)
    source_ids = np.array([f"source-{i}" for i in range(n)], dtype=object)
    weights = rng.uniform(0.5, 3.0, n)
    rms = rng.uniform(0.25, 0.42, n)
    m = rng.uniform(-0.03, 0.02, n)
    c1 = rng.normal(0.0, 0.001, n)
    c2 = rng.normal(0.0, 0.001, n)
    mag_a10 = rng.uniform(24.1, 25.5, n)
    resolution = rng.uniform(0.3, 0.7, n)
    psf_e1 = rng.normal(0.015, 0.004, n)
    psf_e2 = rng.normal(-0.006, 0.004, n)
    bmode = np.ones(n, dtype=bool)
    bmode[250] = False
    lens_z = np.full(n, z_lens)
    redshift_grid = np.linspace(0.0, 2.0, 401)
    source_pdf = np.exp(-0.5 * ((redshift_grid[None, :] - 1.0) / 0.16) ** 2)
    source_pdf = np.repeat(source_pdf, n, axis=0)
    pairs = {
        "lens_id": lens_ids,
        "object_id": source_ids,
        "lens_ra_deg": np.full(n, ra_lens),
        "lens_dec_deg": np.full(n, dec_lens),
        "lens_z": lens_z,
        "source_ra_deg": ra_source,
        "source_dec_deg": dec_source,
        "e1": np.zeros(n),
        "e2": np.zeros(n),
        "weight": weights,
        "rms_e": rms,
        "m": m,
        "c1": c1,
        "c2": c2,
        "mag_a10": mag_a10,
        "resolution": resolution,
        "psf_shape11": 1.0 + psf_e1,
        "psf_shape22": 1.0 - psf_e1,
        "psf_shape12": psf_e2,
        "b_mode_pass": bmode,
    }
    actual_radius, actual_phi = pair_geometry(pairs, cosmo)
    bin_index = np.searchsorted(radial_edges, actual_radius, side="right") - 1
    bin_index = np.clip(bin_index, 0, len(radial_edges) - 2)
    beta = pdf_lensing_efficiency(lens_z, source_ids, source_ids, redshift_grid, source_pdf, cosmo)
    sigma_crit_inv = 4.0 * math.pi * G_MPC / C_KMS**2 * d_l * beta
    selection = []
    for i in range(len(radial_edges) - 1):
        selected = (bin_index == i) & bmode
        selection.append(hsc_selection_bias(
            weights[selected] * sigma_crit_inv[selected] ** 2,
            mag_a10[selected], resolution[selected],
            pairs["psf_shape11"][selected], pairs["psf_shape22"][selected],
            pairs["psf_shape12"][selected],
        ))
    responsivity = np.zeros(len(radial_edges) - 1)
    for i in range(len(responsivity)):
        selected = (bin_index == i) & bmode
        pair_weight = weights[selected] * sigma_crit_inv[selected] ** 2
        responsivity[i] = 1.0 - np.sum(pair_weight * rms[selected] ** 2) / np.sum(pair_weight)
    true_delta_sigma = 2.0e14
    gamma_t = true_delta_sigma * sigma_crit_inv
    c_ind_t, c_ind_x = _to_tangential(actual_phi, c1, c2)
    m_sel = np.array([selection[i]["m_sel"] for i in bin_index])
    csel_1 = np.array([selection[i]["c1_sel"] for i in bin_index])
    csel_2 = np.array([selection[i]["c2_sel"] for i in bin_index])
    csel_t, csel_x = _to_tangential(actual_phi, csel_1, csel_2)
    response_by_pair = responsivity[bin_index]
    observed_t = 2.0 * response_by_pair * ((1.0 + m) * ((1.0 + m_sel) * gamma_t + csel_t) + c_ind_t)
    observed_x = 2.0 * response_by_pair * ((1.0 + m) * csel_x + c_ind_x)
    cos2, sin2 = np.cos(2.0 * actual_phi), np.sin(2.0 * actual_phi)
    pairs["e1"] = -observed_t * cos2 + observed_x * sin2
    pairs["e2"] = -observed_t * sin2 - observed_x * cos2
    return pairs, source_ids, redshift_grid, source_pdf, radial_edges, true_delta_sigma


def test_synthetic_hsc_calibration_and_pdf_beta_recover_injected_delta_sigma(cosmo):
    pairs, source_ids, z, pdf, edges, injected = _synthetic_catalog(cosmo)
    result = stack_s19a_pairs(
        pairs, source_ids, z, pdf, edges, cosmo,
    )
    assert result["science_measurement"] is False
    assert result["counts"]["input_pairs"] == 400
    assert result["counts"]["b_mode_rejected"] == 1
    assert len(result["profile"]) == 4
    for row in result["profile"]:
        assert row["n_pairs"] > 90
        assert row["delta_sigma_weak_msun_mpc2"] == pytest.approx(injected, rel=2e-10)
        assert row["mean_gx"] == pytest.approx(0.0, abs=2e-12)
        assert row["delta_sigma_shape_noise_msun_mpc2"] > 0


def test_missing_calibration_pdf_or_decoded_bmode_is_rejected(cosmo):
    pairs, source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    with pytest.raises(ValueError, match="Missing full photo-z PDF"):
        stack_s19a_pairs(pairs, source_ids[:-1], z, pdf[:-1], edges, cosmo)
    bad = dict(pairs)
    bad["b_mode_pass"] = np.full(len(source_ids), 2, dtype=int)
    with pytest.raises(ValueError, match="decoded boolean mask"):
        stack_s19a_pairs(bad, source_ids, z, pdf, edges, cosmo)


def test_duplicate_pair_and_lossy_float_ids_are_rejected(cosmo):
    pairs, source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    duplicate = {key: np.append(value, value[0]) for key, value in pairs.items()}
    with pytest.raises(ValueError, match=r"duplicate \(lens_id, object_id\)"):
        stack_s19a_pairs(duplicate, np.append(source_ids, source_ids[0]), z,
                         np.vstack([pdf, pdf[0]]), edges, cosmo)
    lossy = dict(pairs)
    lossy["object_id"] = np.arange(len(source_ids), dtype=float)
    with pytest.raises(ValueError, match="floating-point ID"):
        stack_s19a_pairs(lossy, source_ids, z, pdf, edges, cosmo)
