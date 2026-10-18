import unittest

import numpy as np

from measure_satellite_outer_profiles import (detection_arrays, make_stamp_mask,
                                             match_members, profile_at)


class SatelliteOuterProfileHelpersTest(unittest.TestCase):
    def test_member_matching_rejects_ambiguous_and_reused_sources(self):
        detections = [
            {"NUMBER": 1, "x": 10.0, "y": 10.0},
            {"NUMBER": 2, "x": 20.0, "y": 10.0},
            {"NUMBER": 3, "x": 30.0, "y": 10.0},
            {"NUMBER": 4, "x": 30.5, "y": 10.0},
            {"NUMBER": 5, "x": 40.0, "y": 10.0},
        ]
        members = [
            {"INDEX": "0", "X_IMAGE": "11.1", "Y_IMAGE": "11.0", "IN_FRAME": "True"},
            {"INDEX": "1", "X_IMAGE": "21.1", "Y_IMAGE": "11.0", "IN_FRAME": "True"},
            {"INDEX": "2", "X_IMAGE": "30.1", "Y_IMAGE": "10.0", "IN_FRAME": "True"},
            {"INDEX": "3", "X_IMAGE": "41.0", "Y_IMAGE": "11.0", "IN_FRAME": "True"},
            {"INDEX": "4", "X_IMAGE": "41.1", "Y_IMAGE": "11.0", "IN_FRAME": "True"},
        ]
        matches, stats = match_members(members, detections, pixel_scale=0.168)
        self.assertEqual([row["member_index"] for row in matches], [0, 1])
        self.assertEqual(stats["ambiguous"], 1)
        self.assertEqual(stats["duplicate_assignments"], 1)

    def test_target_light_is_preserved_but_neighbor_and_bad_pixels_are_masked(self):
        shape = (100, 100)
        base = np.zeros(shape, dtype=bool)
        base[50, 50] = True
        segmentation = np.zeros(shape, dtype=np.int32)
        segmentation[48:53, 48:53] = 1
        segmentation[48:53, 63:68] = 2
        detections = {
            "number": np.array([1, 2]), "x": np.array([50., 65.]),
            "y": np.array([50., 50.]), "a": np.array([2., 2.]),
            "b": np.array([2., 2.]), "theta": np.array([0., 0.]),
            "kron_scale": np.array([2.5, 2.5]), "extent": np.array([5., 5.]),
        }
        target_mask = make_stamp_mask(base, segmentation, detections, None,
                                      (50., 50.), 10, 10, 1, .168, 40, 2.5)
        self.assertTrue(target_mask[50 - 10, 50 - 10])  # Original bad pixel remains masked.
        self.assertFalse(target_mask[49 - 10, 49 - 10])  # Target source pixels are retained.
        self.assertTrue(target_mask[49 - 10, 64 - 10])  # Neighbor source pixels are excluded.
        control_mask = make_stamp_mask(base, segmentation, detections, None,
                                       (50., 50.), 10, 10, None, .168, 40, 2.5)
        self.assertTrue(control_mask[49 - 10, 49 - 10])
        self.assertTrue(control_mask[49 - 10, 64 - 10])

    def test_constant_surface_brightness_is_recovered_in_each_annulus(self):
        shape = (180, 180)
        image = np.ones(shape, dtype=np.float32)
        variance = np.ones(shape, dtype=np.float32)
        base = np.zeros(shape, dtype=bool)
        segmentation = np.zeros(shape, dtype=np.int32)
        segmentation[88:93, 88:93] = 1
        detections = detection_arrays([{
            "NUMBER": 1, "x": 90., "y": 90., "A_IMAGE": 2., "B_IMAGE": 2.,
            "THETA_RAD": 0., "KRON_RADIUS": 1., "KRON_FLAGS": 0,
        }])
        result = profile_at(image, variance, base, segmentation, detections, None,
                            (90., 90.), .168, 1., 2., 1, .6)
        self.assertIsNotNone(result)
        self.assertTrue(result["passes_coverage"])
        np.testing.assert_allclose(result["sb"], 2., rtol=0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
