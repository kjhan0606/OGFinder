import math
from types import SimpleNamespace

import numpy as np
import pytest

from ogfkit import wlcluster as wl
from s19a_shape_stack import (
    G_MPC,
    C_KMS,
    _to_tangential,
    decode_hsc_b_mode_mask,
    hsc_selection_bias,
    load_hsc_y3_nz_sacc,
    pair_geometry,
    pdf_lensing_efficiency,
    stack_s19a_pairs,
    stack_s19a_tomographic_pairs,
    paired_jackknife_covariance,
    equal_area_source_blocks,
    paired_tomographic_spatial_jackknife,
    subtract_random_lens_profile,
    tomographic_nz_lensing_efficiency,
)
from run_wl_bmode_polarity_sensitivity import calculate_sensitivity


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


def test_b_mode_mask_decoder_requires_explicit_catalog_polarity():
    raw = np.array([0, 1, 0, 1])
    assert decode_hsc_b_mode_mask(raw, masked_value=1).tolist() == [True, False, True, False]
    assert decode_hsc_b_mode_mask(raw, masked_value=0).tolist() == [False, True, False, True]
    with pytest.raises(ValueError, match="not binary"):
        decode_hsc_b_mode_mask([0, 2], masked_value=1)
    with pytest.raises(ValueError, match="verified binary value"):
        decode_hsc_b_mode_mask([0, 1], masked_value=2)


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
        "sigma_e": rng.uniform(0.02, 0.30, n),
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
    assert result["counts"]["b_mode_mask_applied"] is True
    assert result["counts"]["b_mode_mask_rejected"] == 1
    assert len(result["profile"]) == 4
    for row in result["profile"]:
        assert row["n_pairs"] > 90
        assert row["delta_sigma_weak_msun_mpc2"] == pytest.approx(injected, rel=2e-10)
        assert row["mean_gx"] == pytest.approx(0.0, abs=2e-12)
        assert row["delta_sigma_shape_noise_msun_mpc2"] > 0


def test_stack_applies_hsc_resolution_and_aperture_magnitude_cuts(cosmo):
    pairs, source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    pairs["mag_a10"][0] = 25.5001
    pairs["resolution"][1] = 0.2999
    result = stack_s19a_pairs(pairs, source_ids, z, pdf, edges, cosmo)
    assert result["counts"]["hsc_source_selection_rejected"] == 2
    assert result["counts"]["used_pairs"] == 397  # two source cuts plus one B-mode row
    assert result["calibration"]["source_selection"] == (
        "HSC Y3 cuts |e|<2, 0<=sigma_e<=0.4, R2>=0.3, and magA10<=25.5 applied before stacking"
    )


def test_hsc_shape_quality_uses_published_ellipticity_and_sigma_cuts(cosmo):
    pairs, source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    pairs["e1"][0] = 1.2  # HSC catalog |e| may exceed one; published cut is |e|<2.
    pairs["e1"][1] = 2.0
    pairs["sigma_e"][2] = 0.4001
    result = stack_s19a_pairs(pairs, source_ids, z, pdf, edges, cosmo)
    assert result["counts"]["ellipticity_amplitude_ge_2_rejected"] == 1
    assert result["counts"]["sigma_e_outside_0_0p4_rejected"] == 1
    assert result["counts"]["hsc_shape_quality_rejected"] == 2
    assert result["counts"]["used_pairs"] == 397  # two shape cuts plus one B-mode row


def test_bmode_mask_can_be_disabled_for_explicit_sensitivity(cosmo):
    pairs, source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    result = stack_s19a_pairs(
        pairs, source_ids, z, pdf, edges, cosmo, apply_b_mode_mask=False,
    )
    assert result["counts"]["b_mode_mask_applied"] is False
    assert result["counts"]["b_mode_mask_rejected"] == 1
    assert result["counts"]["used_pairs"] == 400
    assert all(row["n_pairs"] > 90 for row in result["profile"])


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


def test_tomographic_nz_efficiency_matches_pdf_and_excludes_bin_zero(cosmo):
    z = np.linspace(0.0, 2.0, 401)
    nz = np.exp(-0.5 * ((z - 1.0) / 0.16) ** 2)
    bins = np.array([1, 2, 3, 4, 0])
    zl = np.full(len(bins), 0.3)
    beta = tomographic_nz_lensing_efficiency(
        zl, bins, z, np.tile(nz, (4, 1)), cosmo,
    )
    one_pdf_beta = pdf_lensing_efficiency(
        np.full(4, 0.3), ["s1", "s2", "s3", "s4"],
        ["s1", "s2", "s3", "s4"], z, np.tile(nz, (4, 1)), cosmo,
    )
    assert beta[:4] == pytest.approx(one_pdf_beta)
    assert beta[4] == 0.0
    with pytest.raises(ValueError, match="exceeds the supplied n\\(z\\) bins"):
        tomographic_nz_lensing_efficiency([0.3], [5], z, np.tile(nz, (4, 1)), cosmo)


def test_bmode_polarity_sensitivity_keeps_both_hypotheses(cosmo):
    pairs, _source_ids, z_pdf, _pdf, edges, _injected = _synthetic_catalog(cosmo)
    pairs["source_zbin"] = np.ones(len(pairs["object_id"]), dtype=int)
    pairs["b_mode_mask_raw"] = np.arange(len(pairs["object_id"]), dtype=np.int8) % 2
    nz = np.exp(-0.5 * ((z_pdf - 1.0) / 0.16) ** 2)
    scenarios, differences = calculate_sensitivity(
        pairs, z_pdf, np.tile(nz, (4, 1)), edges, cosmo,
    )
    exclude_zero = scenarios["exclude_raw_0_hypothesis"]["stack"]
    exclude_one = scenarios["exclude_raw_1_hypothesis"]["stack"]

    assert exclude_zero["counts"]["used_pairs"] == 200
    assert exclude_one["counts"]["used_pairs"] == 200
    assert exclude_zero["counts"]["used_pairs"] + exclude_one["counts"]["used_pairs"] == 400
    assert len(differences) == len(edges) - 1
    assert all(row["delta_sigma_exclude_0_minus_exclude_1_msun_mpc2"] is not None
               for row in differences)
    assert exclude_zero["science_measurement"] is False
    assert exclude_one["science_measurement"] is False


def test_load_hsc_y3_nz_sacc_returns_four_ordered_tracers(tmp_path, monkeypatch):
    z = np.array([0.0, 0.5, 1.0])
    tracers = {
        f"wl_{index}": SimpleNamespace(z=z, nz=np.array([0.0, index + 1.0, 0.0]))
        for index in range(4)
    }
    fake_sacc = SimpleNamespace(Sacc=SimpleNamespace(
        load_fits=lambda filename: SimpleNamespace(tracers=tracers),
    ))
    monkeypatch.setitem(__import__("sys").modules, "sacc", fake_sacc)
    path = tmp_path / "hsc_y3.sacc"
    path.write_bytes(b"fixture")

    loaded_z, loaded_nz = load_hsc_y3_nz_sacc(path)

    assert np.array_equal(loaded_z, z)
    assert loaded_nz.shape == (4, 3)
    assert loaded_nz[:, 1] == pytest.approx([1.0, 2.0, 3.0, 4.0])


def test_load_hsc_y3_nz_sacc_rejects_missing_or_misaligned_tracer(tmp_path, monkeypatch):
    path = tmp_path / "hsc_y3.sacc"
    path.write_bytes(b"fixture")
    good_z = np.array([0.0, 0.5, 1.0])

    missing = {f"wl_{index}": SimpleNamespace(z=good_z, nz=np.ones(3)) for index in range(3)}
    monkeypatch.setitem(__import__("sys").modules, "sacc", SimpleNamespace(
        Sacc=SimpleNamespace(load_fits=lambda _filename: SimpleNamespace(tracers=missing)),
    ))
    with pytest.raises(ValueError, match=r"missing a valid n\(z\) tracer wl_3"):
        load_hsc_y3_nz_sacc(path)

    shifted = dict(missing)
    shifted["wl_3"] = SimpleNamespace(z=np.array([0.0, 0.6, 1.0]), nz=np.ones(3))
    monkeypatch.setitem(__import__("sys").modules, "sacc", SimpleNamespace(
        Sacc=SimpleNamespace(load_fits=lambda _filename: SimpleNamespace(tracers=shifted)),
    ))
    with pytest.raises(ValueError, match="share one redshift grid"):
        load_hsc_y3_nz_sacc(path)


def test_load_hsc_y3_nz_sacc_astropy_fallback_deduplicates_exact_tables(tmp_path, monkeypatch):
    fits = pytest.importorskip("astropy.io.fits")
    monkeypatch.delitem(__import__("sys").modules, "sacc", raising=False)
    z = np.array([0.0, 0.5, 1.0])
    hdus = [fits.PrimaryHDU()]
    for index in range(4):
        nz = np.array([0.0, float(index + 1), 0.0])
        for _ in range(2):  # Official HSC file repeats each n(z) tracer exactly.
            table = fits.BinTableHDU.from_columns([
                fits.Column(name="z", format="D", array=z),
                fits.Column(name="nz", format="D", array=nz),
            ])
            table.header["SACCTYPE"] = "tracer"
            table.header["SACCCLSS"] = "NZ"
            table.header["SACCNAME"] = f"wl_{index}"
            hdus.append(table)
    path = tmp_path / "hsc_y3.sacc"
    fits.HDUList(hdus).writeto(path)

    loaded_z, loaded_nz = load_hsc_y3_nz_sacc(path)

    assert np.array_equal(loaded_z, z)
    assert loaded_nz[:, 1] == pytest.approx([1.0, 2.0, 3.0, 4.0])


def test_load_hsc_y3_nz_sacc_astropy_fallback_rejects_conflicting_duplicates(tmp_path, monkeypatch):
    fits = pytest.importorskip("astropy.io.fits")
    monkeypatch.delitem(__import__("sys").modules, "sacc", raising=False)
    z = np.array([0.0, 0.5, 1.0])
    hdus = [fits.PrimaryHDU()]
    for index in range(4):
        values = [0.0, float(index + 1), 0.0]
        if index == 2:
            duplicate_values = [0.0, 99.0, 0.0]
        else:
            duplicate_values = values
        for nz in (values, duplicate_values):
            table = fits.BinTableHDU.from_columns([
                fits.Column(name="z", format="D", array=z),
                fits.Column(name="nz", format="D", array=nz),
            ])
            table.header["SACCTYPE"] = "tracer"
            table.header["SACCCLSS"] = "NZ"
            table.header["SACCNAME"] = f"wl_{index}"
            hdus.append(table)
    path = tmp_path / "hsc_y3_conflict.sacc"
    fits.HDUList(hdus).writeto(path)

    with pytest.raises(ValueError, match=r"conflicting duplicate n\(z\) tracer wl_2"):
        load_hsc_y3_nz_sacc(path)


def test_tomographic_stack_uses_exported_source_zbin_and_recovers_signal(cosmo):
    pairs, source_ids, z, pdf, edges, injected = _synthetic_catalog(cosmo)
    # CSV catalogues commonly arrive as strings, so accept exact integer text.
    pairs["source_zbin"] = np.ones(len(source_ids), dtype=str)
    # HSC zbin 0 is already the one B-mode-rejected row in this fixture.
    pairs["source_zbin"][250] = 0
    nz_by_bin = np.tile(pdf[0], (4, 1))
    result = stack_s19a_tomographic_pairs(pairs, z, nz_by_bin, edges, cosmo)
    assert result["science_measurement"] is False
    assert result["counts"]["unassigned_source_zbin"] == 1
    assert result["counts"]["used_pairs"] == 399
    for row in result["profile"]:
        assert row["delta_sigma_weak_msun_mpc2"] == pytest.approx(injected, rel=2e-10)
    assert result["calibration"]["photoz_geometry"] == (
        "HSC Y3 tomographic n(z) integrated to pair-wise mean SigmaCritInv"
    )


def test_tomographic_stack_rejects_inconsistent_source_bin(cosmo):
    pairs, _source_ids, z, pdf, edges, _ = _synthetic_catalog(cosmo)
    pairs["source_zbin"] = np.ones(len(pairs["object_id"]), dtype=int)
    pairs["lens_id"][1] = "lens-2"
    pairs["object_id"][1] = pairs["object_id"][0]
    pairs["source_zbin"][1] = 2
    with pytest.raises(ValueError, match="inconsistent source_zbin"):
        stack_s19a_tomographic_pairs(pairs, z, np.tile(pdf[0], (4, 1)), edges, cosmo)


def test_random_lens_profile_subtracts_signal_but_does_not_fake_errors():
    def result(ds, cross):
        return {
            "schema": "hsc-s19a-shape-stack-v1", "science_measurement": False,
            "profile": [{
                "bin": 0, "r_min_mpc": 0.4, "r_max_mpc": 0.7,
                "n_pairs": 100, "delta_sigma_weak_msun_mpc2": ds,
                "cross_delta_sigma_weak_msun_mpc2": cross,
                "delta_sigma_shape_noise_msun_mpc2": 0.2,
            }],
        }
    subtracted = subtract_random_lens_profile(result(15.0, 2.0), result(3.0, 1.0))
    row = subtracted["profile"][0]
    assert row["delta_sigma_random_subtracted_msun_mpc2"] == pytest.approx(12.0)
    assert row["cross_delta_sigma_random_subtracted_msun_mpc2"] == pytest.approx(1.0)
    assert "not_estimated" in row["covariance_status"]
    assert not any("error" in key for key in row)


def test_paired_jackknife_covariance_uses_subtracted_replicas():
    lens = np.array([[10.0], [12.0], [8.0]])
    random = np.array([[2.0], [1.0], [3.0]])
    result = paired_jackknife_covariance(lens, random)
    assert result["n_blocks"] == 3
    assert result["mean_random_subtracted"] == pytest.approx([8.0])
    assert np.allclose(result["covariance"], [[12.0]])
    assert result["standard_error"] == pytest.approx([math.sqrt(12.0)])
    with pytest.raises(ValueError, match="share shape"):
        paired_jackknife_covariance(lens, random[:2])


def test_equal_area_source_blocks_are_deterministic_across_ra_wrap():
    blocks = equal_area_source_blocks(
        [359.9, 0.1, 180.1], [0.0, 0.0, 0.0], n_ra=4, n_sin_dec=2,
    )
    assert blocks.tolist() == [7, 4, 6]
    assert equal_area_source_blocks([0.1], [0.0], n_ra=4, n_sin_dec=2)[0] == blocks[1]


def test_paired_tomographic_jackknife_recomputes_both_stacks(cosmo):
    pairs, source_ids, z, pdf, edges, _injected = _synthetic_catalog(cosmo)
    lens_pairs = {key: value.copy() for key, value in pairs.items()}
    random_pairs = {key: value.copy() for key, value in pairs.items()}
    lens_pairs["source_zbin"] = np.ones(len(source_ids), dtype=int)
    random_pairs["source_zbin"] = np.ones(len(source_ids), dtype=int)
    random_pairs["lens_id"][:] = "random-lens-1"
    nz_by_bin = np.tile(pdf[0], (4, 1))

    result = paired_tomographic_spatial_jackknife(
        lens_pairs, random_pairs, z, nz_by_bin, edges, cosmo,
        n_ra_blocks=7200, n_sin_dec_blocks=7200,
    )

    assert result["science_measurement"] is False
    assert result["n_occupied_blocks"] >= 3
    assert len(result["profile"]) == len(edges) - 1
    assert np.asarray(result["covariance"]).shape == (2 * (len(edges) - 1),) * 2
    assert all(row["delta_sigma_random_subtracted_msun_mpc2"] == pytest.approx(0.0, abs=1e-8)
               for row in result["profile"])
    assert all(row["delta_sigma_jackknife_error_msun_mpc2"] >= 0 for row in result["profile"])
