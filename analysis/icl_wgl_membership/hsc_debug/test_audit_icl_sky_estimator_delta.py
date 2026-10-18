import numpy as np

from audit_icl_sky_estimator_delta import paired_block_contrast, partial_pearson


def make_rows(plane_equals_pixel=True):
    rows = []
    for i in range(48):
        z = 0.15 + 0.006 * i
        richness = 18 + ((i * 7) % 23)
        satellite = np.sin(i * 0.57) + 0.3 * z - 0.01 * richness
        pixel = satellite + np.cos(i * 0.31) * 0.4 + 0.2 * z
        plane = pixel if plane_equals_pixel else pixel + np.sin(i * 0.19) * 0.8
        rows.append({
            "icl_delta_pixel_sfd_njy_arcsec2": pixel,
            "icl_delta_pixel_sfd_after_plane_njy_arcsec2": plane,
            "satellite_delta_center_sfd_njy_arcsec2": satellite,
            "z": z,
            "richness": richness,
        })
    labels = np.repeat([f"block_{i}" for i in range(8)], 6)
    return rows, labels


def test_partial_correlation_removes_linear_covariates():
    z = np.linspace(0, 1, 30)
    richness = np.linspace(1, 3, 30) ** 2
    residual = np.sin(np.arange(30))
    x = 4 * z - 2 * richness + residual
    y = -z + 0.3 * richness + 2 * residual
    assert np.isclose(partial_pearson(x, y, z, richness), 1.0)


def test_paired_bootstrap_is_deterministic_and_zero_for_identical_estimator():
    rows, labels = make_rows()
    first = paired_block_contrast(rows, labels, replicas=1200, seed=19)
    second = paired_block_contrast(rows, labels, replicas=1200, seed=19)
    assert first == second
    assert first["point_delta_plane_minus_pixel"] == 0
    assert first["paired_block_bootstrap_delta_ci95"] == [0.0, 0.0]
    assert first["valid_paired_replicas"] == 1200


def test_paired_bootstrap_reports_same_fixed_pair_count():
    rows, labels = make_rows(plane_equals_pixel=False)
    result = paired_block_contrast(rows, labels, replicas=1200, seed=20)
    assert result["n_pairs"] == len(rows)
    assert result["n_connected_5deg_blocks"] == len(set(labels))
    assert result["valid_paired_replicas"] > 1000
    assert result["point_delta_plane_minus_pixel"] != 0
