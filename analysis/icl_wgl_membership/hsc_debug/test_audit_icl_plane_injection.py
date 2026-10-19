import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM

from audit_icl_plane_injection import (
    SHAPES,
    block_bootstrap_mean,
    measure_transfer,
    select_shape_rows,
)


def test_block_bootstrap_is_deterministic_for_constant_paired_response():
    values = np.full(36, 0.42)
    labels = np.repeat([f"block_{i}" for i in range(9)], 4)
    first = block_bootstrap_mean(values, labels, replicas=1200, seed=38)
    assert first == block_bootstrap_mean(values, labels, replicas=1200, seed=38)
    assert first["mean"] == 0.42
    assert first["ci95"] == [0.42, 0.42]


def test_block_bootstrap_requires_spatial_blocks():
    with pytest.raises(ValueError, match="six connected"):
        block_bootstrap_mean(np.ones(5), np.array(["a", "b", "c", "d", "e"]),
                             replicas=1000, seed=4)


def test_injection_recovers_exact_raw_annulus_mean_on_synthetic_frame():
    shape = (2200, 2200)
    yy, xx = np.indices(shape, dtype=np.float32)
    image = 0.0003 * xx - 0.0002 * yy
    variance = np.ones(shape, dtype=np.float32)
    mask = np.zeros(shape, dtype=bool)
    info = {"center_pixel_xy": [1100.0, 1100.0], "pixel_scale_arcsec": 1.0}
    result = measure_transfer(
        image, variance, mask, info, conversion=1.0, z=0.3,
        scale_kpc=100.0, axis_ratio=1.0,
        cosmology=FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0), amplitude=1.0)
    assert result["annulus_pixels"] > 100
    assert result["n_plane_tiles_baseline"] >= 10
    assert result["raw_recovery_fraction"] == pytest.approx(1.0, abs=2e-6)
    assert np.isfinite(result["plane_corrected_recovery_fraction"])


def test_shape_group_keeps_all_pairs_without_mixing_injection_models():
    rows = []
    for band, count in (("i2", 3), ("i", 2)):
        for pair in range(count):
            for scale, ratio in SHAPES:
                rows.append({"cluster_id": f"{band}_{pair}",
                             "control_id": f"ctrl_{band}_{pair}",
                             "filter": band, "scale_kpc": scale,
                             "axis_ratio": ratio})
    selected = select_shape_rows(rows, "i2", 100, 0.65, expected_pairs=3)
    assert len(selected) == 3
    assert {r["cluster_id"] for r in selected} == {"i2_0", "i2_1", "i2_2"}
    incomplete = [r for r in rows if not (r["cluster_id"] == "i2_0"
                   and r["scale_kpc"] == 100.0 and r["axis_ratio"] == 0.65)]
    with pytest.raises(ValueError, match="fixed paired sample"):
        select_shape_rows(incomplete, "i2", 100, 0.65, expected_pairs=3)


def test_extended_500_kpc_profiles_are_explicit_stress_cases():
    assert (500.0, 1.0) in SHAPES
    assert (500.0, 0.65) in SHAPES
    assert len(SHAPES) == 6
