import unittest

import numpy as np

from stack_empirical_psf_outer_field_wings import radial_profile, stack_pairs


class OuterFieldWingTests(unittest.TestCase):
    def test_radial_profile_recovers_background_subtracted_core_normalized_bins(self):
        yy, xx = np.indices((501, 501), dtype=float)
        radius = np.hypot(xx - 250.2, yy - 249.7)
        image = 4.0 + 1000.0 * np.exp(-0.5 * (radius / 2.2) ** 2)
        variance = np.ones_like(image)
        mask = np.zeros_like(image, dtype=bool)
        rows = radial_profile(image, variance, mask, (250.2, 249.7), 1.0)
        first = next(row for row in rows if row["sky_r_min_pix"] == 120 and row["r_min_kpc"] == 0)
        outer = next(row for row in rows if row["sky_r_min_pix"] == 120 and row["r_min_kpc"] == 275)
        self.assertEqual(first["status"], "ok")
        self.assertGreater(first["target_core_flux"], 0)
        self.assertAlmostEqual(first["background"], 4.0, places=8)
        self.assertEqual(outer["status"], "ok")
        self.assertLess(abs(outer["profile_per_target_core_per_pixel"]), 1e-8)

    def test_missing_azimuth_sector_blocks_radial_bin(self):
        yy, xx = np.indices((501, 501), dtype=float)
        dx, dy = xx - 250, yy - 250
        radius = np.hypot(dx, dy)
        image = 1.0 + 100.0 * np.exp(-radius**2 / 8.0)
        variance = np.ones_like(image)
        mask = np.zeros_like(image, dtype=bool)
        angle = np.arctan2(dy, dx)
        target_ring = (radius >= 275) & (radius < 300) & (angle >= -np.pi) & (angle < -3 * np.pi / 4)
        mask[target_ring] = True
        rows = radial_profile(image, variance, mask, (250, 250), 1.0)
        outer = next(row for row in rows if row["sky_r_min_pix"] == 120 and row["r_min_kpc"] == 275)
        self.assertEqual(outer["status"], "blocked_radial_coverage")
        self.assertLess(outer["n_covered_sectors"], 8)

    def test_control_requires_matching_target_core_flux_per_sky_annulus(self):
        yy, xx = np.indices((501, 501), dtype=float)
        radius = np.hypot(xx - 250, yy - 250)
        image = 2.0 + 100.0 * np.exp(-radius**2 / 8.0)
        variance = np.ones_like(image)
        mask = np.zeros_like(image, dtype=bool)
        rows = radial_profile(image, variance, mask, (250, 250), 1.0,
                              core_flux={(160.0, 190.0): 25.0})
        first = next(row for row in rows if row["sky_r_min_pix"] == 120 and row["r_min_kpc"] == 25)
        second = next(row for row in rows if row["sky_r_min_pix"] == 160 and row["r_min_kpc"] == 25)
        self.assertEqual(first["status"], "blocked_target_core_flux")
        self.assertIsNone(first["profile_per_target_core_per_pixel"])
        self.assertEqual(second["status"], "ok")
        self.assertEqual(second["target_core_flux"], 25.0)

    def test_frame_paired_stack_resamples_frames_and_is_reproducible(self):
        rows = []
        for candidate, target, control in (("a", 3.0, 1.0), ("b", 7.0, 2.0), ("c", 11.0, 5.0)):
            for role, value in (("target", target), ("control", control)):
                rows.append({"candidate_id": candidate, "role": role, "status": "ok",
                             "profile_per_target_core_per_pixel": value,
                             "sky_r_min_pix": 120, "sky_r_max_pix": 150,
                             "r_min_kpc": 0.0, "r_max_kpc": 25.0})
        first = stack_pairs(rows, n_bootstrap=500, seed=7)
        second = stack_pairs(rows, n_bootstrap=500, seed=7)
        self.assertEqual(first, second)
        self.assertEqual(first[0]["n_paired_frames"], 3)
        self.assertAlmostEqual(first[0]["paired_difference_frame_equal_mean"], 13 / 3)
        self.assertEqual(first[0]["bootstrap_unit"], "matched HSC target/control frame pair")


if __name__ == "__main__":
    unittest.main()
