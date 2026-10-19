import numpy as np

from diagnose_empirical_psf import objects_from_rows


def test_missing_unrelated_detection_morphology_uses_finite_mask_fallbacks():
    rows = [{
        "NUMBER": "1", "X_IMAGE": "10", "Y_IMAGE": "20",
        "A_IMAGE": "1.2", "B_IMAGE": "", "THETA_RAD": "nan",
        "KRON_RADIUS": "0", "KRON_FLAGS": "64",
    }]

    objects, kron, flags = objects_from_rows(rows, x0=4, y0=5)

    np.testing.assert_allclose([objects[key][0] for key in ("x", "y", "a", "b", "theta")],
                               (5, 14, 1.2, 0.5, 0))
    assert kron.tolist() == [0]
    assert flags.tolist() == [64]
    assert np.isfinite(objects.view(float)).all()


def test_missing_position_is_kept_outside_the_stamp_for_segmentation_masking():
    rows = [{
        "NUMBER": "1", "X_IMAGE": "", "Y_IMAGE": "nan",
        "A_IMAGE": "", "B_IMAGE": "", "THETA_RAD": "",
        "KRON_RADIUS": "", "KRON_FLAGS": "",
    }]

    objects, kron, flags = objects_from_rows(rows, x0=100, y0=200)

    assert objects["x"][0] < -1e5
    assert objects["y"][0] < -1e5
    assert np.isfinite(objects.view(float)).all()
    assert kron.tolist() == [1.2]
    assert flags.tolist() == [1]
