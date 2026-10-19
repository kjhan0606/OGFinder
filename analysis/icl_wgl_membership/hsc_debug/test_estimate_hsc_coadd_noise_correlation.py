import numpy as np

from estimate_hsc_coadd_noise_correlation import (
    _matched_control_source_mask,
    lag_correlations,
    load_matched_controls,
)


def test_lag_correlation_recovers_correlated_noise_and_unit_diagonal():
    rng = np.random.default_rng(20261019)
    white = rng.normal(size=(128, 128))
    image = (4 * white + np.roll(white, 1, axis=1) + np.roll(white, -1, axis=1)) / np.sqrt(18)
    variance = np.ones_like(image)
    valid = np.ones_like(image, dtype=bool)
    result = lag_correlations(image, variance, valid, max_lag=2)
    assert abs(result[(0, 0)][0] - 1.0) < 0.10
    assert result[(0, 1)][0] > 0.15
    assert result[(0, -1)][0] > 0.15


def test_lag_correlation_ignores_masked_extreme_pixels():
    rng = np.random.default_rng(22)
    image = rng.normal(size=(64, 64))
    variance = np.ones_like(image)
    valid = np.ones_like(image, dtype=bool)
    valid[:8] = False
    image[:8] = 1e6
    result = lag_correlations(image, variance, valid, max_lag=1)
    assert abs(result[(0, 0)][0] - 1.0) < 0.15
    assert abs(result[(0, 1)][0]) < 0.15


def test_matched_controls_follow_target_id_and_filter(tmp_path):
    import hashlib
    import json

    clusters, completed, targets = [], [], []
    (tmp_path / "images").mkdir()
    for index in range(57):
        ident = str(1000 + index)
        filt = "i" if index % 2 else "i2"
        image_name = f"control_{ident}.fits"
        image_path = tmp_path / "images" / image_name
        image_path.write_bytes(f"image-{ident}".encode())
        digest = hashlib.sha256(image_path.read_bytes()).hexdigest()
        clusters.append({
            "id": f"control-{ident}", "target_cluster_id": ident,
            "coadd_header_filter": filt, "image": image_name,
            "bad_mask": {"exclude_bits": [1, 2]},
            "bright_star_mask": {"exclude_bits": [512]},
        })
        completed.append({"id": f"control-{ident}", "file": image_name, "sha256": digest,
                          "review": {"pixel_scale_arcsec": 0.168}})
        targets.append({"cluster_id": ident, "filter": filt})
    manifest_path = tmp_path / "download_manifest.json"
    status_path = tmp_path / "download_status.json"
    manifest_path.write_text(json.dumps({"frame_role": "matched_control", "science_measurement": False,
                                        "clusters": clusters}))
    status_path.write_text(json.dumps({"status": "complete", "n_expected": 57,
                                      "completed": completed}))

    rows = load_matched_controls({"frame_headers": targets}, manifest_path, status_path)

    assert len(rows) == 57
    assert rows[0]["role"] == "matched_control"
    assert rows[0]["filter"] == "i2"
    assert rows[0]["mask_exclude_bits"] == [1, 2, 512]


def test_control_source_mask_uses_same_kron_mask_builder():
    rng = np.random.default_rng(5)
    image = rng.normal(0.0, 0.5, size=(128, 128)).astype(np.float32)
    yy, xx = np.indices(image.shape)
    image += (20.0 * np.exp(-((xx - 63.2) ** 2 + (yy - 65.1) ** 2) / 8.0)).astype(np.float32)
    variance = np.ones_like(image)
    image_mask = np.zeros_like(image, dtype=np.int64)

    mask = _matched_control_source_mask(image, variance, image_mask, 0.168)

    assert mask.shape == image.shape
    assert mask[65, 63]
    assert not mask[8, 8]

