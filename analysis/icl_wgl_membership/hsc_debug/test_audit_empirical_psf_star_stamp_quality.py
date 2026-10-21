import numpy as np

from audit_empirical_psf_star_stamp_quality import neighbor_summary, precompute_rings, ring_support


def test_precomputed_rings_are_centered_and_annuli_are_disjoint():
    geometry = precompute_rings((401, 401), 200.0, 200.0)

    assert geometry["core_ring"].shape == (401, 401)
    assert 430 < np.count_nonzero(geometry["core_ring"]) < 445
    inner = geometry["wing_48_64"]["ring"]
    outer = geometry["wing_64_96"]["ring"]
    assert not np.any(inner & outer)
    assert np.all(geometry["wing_48_64"]["sector"] == geometry["wing_64_96"]["sector"])


def test_ring_support_counts_open_pixels_and_azimuthal_sectors():
    mask = np.zeros((101, 101), dtype=bool)
    mask[:, 50:] = True
    result = ring_support(mask, 50.0, 50.0, 10.0, 20.0)

    assert 0.45 < result["available_fraction"] < 0.55
    assert result["sectors_at_least_25pct_open"] == 4


def test_neighbor_summary_excludes_target_and_counts_strict_radii():
    rows = [
        {"NUMBER": "1", "X_IMAGE": "51", "Y_IMAGE": "51", "FLUX_ISO": "1000"},
        {"NUMBER": "2", "X_IMAGE": "56", "Y_IMAGE": "51", "FLUX_ISO": "100"},
        {"NUMBER": "3", "X_IMAGE": "62", "Y_IMAGE": "51", "FLUX_ISO": "50"},
        {"NUMBER": "4", "X_IMAGE": "151", "Y_IMAGE": "51", "FLUX_ISO": "500"},
    ]
    result = neighbor_summary(rows, 1, 50.0, 50.0, 1000.0)

    assert result["nearest_neighbor_distance_pix"] == 5.0
    assert result["n_neighbors_within_12pix"] == 2
    assert result["n_neighbors_within_96pix"] == 2
    assert result["max_neighbor_to_target_flux_iso_ratio_within_96pix"] == 0.1
