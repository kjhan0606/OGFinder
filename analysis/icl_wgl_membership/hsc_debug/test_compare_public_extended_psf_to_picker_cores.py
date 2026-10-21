from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compare_public_extended_psf_to_picker_cores import measure_public_core


def test_public_core_crop_is_normalized_only_for_shape_metrics():
    y, x = np.indices((101, 101), dtype=float)
    image = np.exp(-((x - 50) ** 2 + (y - 50) ** 2) / (2 * 3.0 ** 2))
    image /= image.sum()
    metrics, crop_fraction = measure_public_core(image, side=43)
    assert 0 < crop_fraction < 1
    assert metrics["signed_sum"] == pytest.approx(1.0)
    assert metrics["encircled_energy_r50_arcsec"] > 0
    assert metrics["shape_yx"] == [43, 43]


@pytest.mark.parametrize("image", [np.ones((10, 10)), np.full((11, 11), np.nan), np.zeros((11, 11))])
def test_public_core_rejects_invalid_image(image):
    with pytest.raises(ValueError):
        measure_public_core(image, side=5)


def test_public_core_rejects_even_crop():
    with pytest.raises(ValueError):
        measure_public_core(np.ones((11, 11)), side=4)
