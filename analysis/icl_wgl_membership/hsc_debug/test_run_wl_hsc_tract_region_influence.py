#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path

import numpy as np

from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
    unit_vectors,
)


class TractRegionInfluenceTests(unittest.TestCase):
    def test_connected_components_cross_ra_zero(self):
        centers = unit_vectors([359.5, 0.5, 20.0], [0.0, 0.0, 0.0])
        labels = tract_region_labels(centers, adjacency_deg=2.0)
        self.assertEqual(labels[0], labels[1])
        self.assertNotEqual(labels[1], labels[2])

    def test_nearest_tract_assigns_region(self):
        centers = unit_vectors([359.5, 0.5, 20.0], [0.0, 0.0, 0.0])
        regions = tract_region_labels(centers, adjacency_deg=2.0)
        result = nearest_region_labels([0.0, 20.1], [0.0, 0.0], centers, regions)
        self.assertEqual(result[0], regions[0])
        self.assertEqual(result[1], regions[2])

    def test_patch_qa_groups_patch_centers_by_tract(self):
        text = (
            "# skymap_id,tract,patch,ra,dec\n"
            "82780004,8278,0,359.0,0.0\n"
            "82780004,8278,1,1.0,0.0\n"
            "82790004,8279,0,10.0,2.0\n"
            "82790004,8279,1,10.2,2.0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "patch_qa.csv"
            path.write_text(text, encoding="utf-8")
            keys, centers = tract_centers_from_patch_qa(path)
        self.assertEqual(keys, ["8278", "8279"])
        self.assertEqual(centers.shape, (2, 3))
        np.testing.assert_allclose(np.linalg.norm(centers, axis=1), 1.0)
        angle = float(np.rad2deg(np.arctan2(centers[0, 1], centers[0, 0])) % 360)
        self.assertLess(min(angle, 360.0 - angle), 1e-6)

    def test_invalid_coordinates_and_threshold_fail(self):
        with self.assertRaises(ValueError):
            unit_vectors([1.0], [91.0])
        with self.assertRaises(ValueError):
            tract_region_labels(np.eye(3), adjacency_deg=0)


if __name__ == "__main__":
    unittest.main()
