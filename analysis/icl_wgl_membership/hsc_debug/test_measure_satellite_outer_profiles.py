import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from measure_satellite_outer_profiles import (cohort_name, detection_arrays,
                                             make_stamp_mask, match_members,
                                             parse_detections, profile_at, resolve_image_paths,
                                             summarize)


class SatelliteOuterProfileHelpersTest(unittest.TestCase):
    def test_full_manifest_defaults_to_an_explicit_full_sample_cohort(self):
        self.assertEqual(cohort_name({"preparation_cohort": None}), "full_sample")
        self.assertEqual(cohort_name({"preparation_cohort": "mask_stress"}), "mask_stress")

    def test_relative_images_resolve_from_manifest_image_directory(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "images").mkdir()
            (root / "images" / "cluster_1.fits").write_bytes(b"fits")
            manifest = root / "download_manifest.json"
            frames = resolve_image_paths([{"id": "1", "image": "cluster_1.fits"}], manifest)
            self.assertEqual(frames[0]["image"], str((root / "images" / "cluster_1.fits").resolve()))
            self.assertEqual(frames[0]["id"], "1")

    def test_empty_profile_result_still_writes_reviewable_tables(self):
        with TemporaryDirectory() as temporary:
            output = Path(temporary)
            result = summarize([], output, replicas=100, seed=10)
            self.assertEqual(result, {})
            self.assertEqual((output / "cluster_profiles.tsv").read_text().splitlines()[0],
                             "cluster_id\tfilter\tcohort\tradial_bin_kpc\tn_satellites\t"
                             "mean_target_sb_njy_arcsec2\tmean_control_sb_njy_arcsec2\t"
                             "mean_delta_sb_njy_arcsec2")
            self.assertTrue((output / "filter_profiles.tsv").is_file())

    def test_flagged_detection_with_blank_minor_axis_uses_finite_fallback(self):
        parsed = parse_detections([{
            "NUMBER": "3", "X_IMAGE": "20.0", "Y_IMAGE": "30.0",
            "A_IMAGE": "0.57", "B_IMAGE": "", "THETA_RAD": "0.8",
            "KRON_RADIUS": "0.0", "KRON_FLAGS": "64",
        }])
        arrays = detection_arrays(parsed)
        self.assertEqual(len(parsed), 1)
        self.assertEqual(arrays["b"][0], 1.0)
        self.assertEqual(arrays["kron_scale"][0], 3.0)
        self.assertFalse(arrays["shape_valid"][0])
        self.assertGreaterEqual(arrays["extent"][0], 2.0 / .168)
        self.assertTrue(np.isfinite(arrays["extent"]).all())

    def test_partial_mask_flag_keeps_measured_kron_radius(self):
        parsed = parse_detections([{
            "NUMBER": "3", "X_IMAGE": "20.0", "Y_IMAGE": "30.0",
            "A_IMAGE": "2.0", "B_IMAGE": "1.2", "THETA_RAD": "0.8",
            "KRON_RADIUS": "1.7", "KRON_FLAGS": "32",
        }])
        arrays = detection_arrays(parsed)
        self.assertTrue(arrays["shape_valid"][0])
        self.assertAlmostEqual(arrays["kron_scale"][0], 2.5 * 1.7)

    def test_invalid_neighbor_shape_uses_two_arcsec_circular_mask(self):
        shape = (100, 100)
        segmentation = np.zeros(shape, dtype=np.int32)
        segmentation[48:53, 48:53] = 1
        segmentation[48:53, 63:68] = 2
        detections = {
            "number": np.array([1, 2]), "x": np.array([50., 65.]),
            "y": np.array([50., 50.]), "a": np.array([2., 1.]),
            "b": np.array([2., 1.]), "theta": np.array([0., 0.]),
            "kron_scale": np.array([2.5, 3.]), "extent": np.array([5., 12.]),
            "shape_valid": np.array([True, False]),
        }
        mask = make_stamp_mask(np.zeros(shape, bool), segmentation, detections, None,
                               (50., 50.), 10, 10, 1, .168, 40, 2.5)
        self.assertTrue(mask[50 - 10, 73 - 10])

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
            "shape_valid": np.array([True, True]),
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
