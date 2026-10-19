import numpy as np

from estimate_hsc_coadd_noise_correlation import lag_correlations, load_matched_controls


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
    for index in range(57):
        ident = str(1000 + index)
        filt = "i" if index % 2 else "i2"
        image_name = f"control_{ident}.fits"
        image_path = tmp_path / image_name
        image_path.write_bytes(f"image-{ident}".encode())
        digest = hashlib.sha256(image_path.read_bytes()).hexdigest()
        clusters.append({
            "id": f"control-{ident}", "target_cluster_id": ident,
            "coadd_header_filter": filt, "image": image_name,
            "bad_mask": {"exclude_bits": [1, 2]},
            "bright_star_mask": {"exclude_bits": [512]},
        })
        completed.append({"id": f"control-{ident}", "file": image_name, "sha256": digest})
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

