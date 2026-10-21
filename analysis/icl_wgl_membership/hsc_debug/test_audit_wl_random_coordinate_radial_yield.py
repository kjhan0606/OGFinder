import unittest

from audit_wl_random_coordinate_radial_yield import (
    separation_arcmin,
    spherical_annulus_area_arcmin2,
)


class RandomCoordinateRadialYieldTests(unittest.TestCase):
    def test_separation_handles_ra_wrap(self):
        self.assertAlmostEqual(separation_arcmin(359.9, 0.0, 0.1, 0.0), 12.0, places=8)

    def test_spherical_annulus_area_is_positive_and_monotone(self):
        inner = spherical_annulus_area_arcmin2(10.0, 15.0)
        wider = spherical_annulus_area_arcmin2(10.0, 20.0)
        self.assertGreater(inner, 0.0)
        self.assertGreater(wider, inner)
        with self.assertRaises(ValueError):
            spherical_annulus_area_arcmin2(15.0, 10.0)


if __name__ == "__main__":
    unittest.main()
