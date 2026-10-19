#!/usr/bin/env python3
import unittest

from audit_satellite_icl_perbin_block_scale import partial_pearson as baseline_partial
from audit_satellite_icl_region_fixed_effects import partial_correlation


class RegionFixedEffectsTests(unittest.TestCase):
    def test_base_model_reproduces_existing_partial_correlation(self):
        rows = []
        for i in range(24):
            x = i * i + (i % 4) * 0.2
            rows.append({
                "icl_delta_200_300": x,
                "satellite_delta_50_75": 0.7 * x + ((i * 7) % 9),
                "satellite_delta": 0.7 * x + ((i * 7) % 9),
                "redshift": 0.3 + i * 0.01,
                "richness": 15.0 + ((i * 11) % 17),
                "target_region": i % 3,
                "control_region": (i + 1) % 3,
            })
        got = partial_correlation(rows)
        want = baseline_partial(rows)
        self.assertAlmostEqual(got["partial_r"], want, places=14)
        self.assertEqual(got["residual_degrees_of_freedom"], len(rows) - 3)

    def test_region_dummies_remove_shared_field_offset(self):
        rows = []
        for region in range(3):
            for i in range(20):
                noise_x = ((i * 7 + region * 3) % 17) - 8
                noise_y = ((i * 11 + region * 5) % 19) - 9
                rows.append({
                    "icl_delta_200_300": 30.0 * region + noise_x,
                    "satellite_delta_50_75": 20.0 * region + noise_y,
                    "redshift": 0.3 + 0.001 * (i + region),
                    "richness": 15.0 + ((i * 13 + region) % 23),
                    "target_region": region,
                    "control_region": region,
                })
        base = partial_correlation(rows)["partial_r"]
        adjusted = partial_correlation(rows, region_effects=("target_region",))["partial_r"]
        self.assertGreater(base, 0.8)
        self.assertLess(abs(adjusted), abs(base))

    def test_identical_endpoint_regions_are_counted_once(self):
        rows = [{
            "icl_delta_200_300": float(i * i + i % 3),
            "satellite_delta_50_75": float(0.4 * i * i + (i * 7) % 11),
            "redshift": 0.2 + (i % 9) * 0.03,
            "richness": 10 + (i * 5) % 17,
            "target_region": i % 3,
            "control_region": i % 3,
        } for i in range(36)]
        target = partial_correlation(rows, region_effects=("target_region",))
        both = partial_correlation(rows, region_effects=("target_region", "control_region"))
        self.assertAlmostEqual(target["partial_r"], both["partial_r"], places=14)
        self.assertEqual(both["redundant_region_terms"], ["control_region[1]", "control_region[2]"])

    def test_rank_deficient_region_design_fails_explicitly(self):
        rows = [{
            "icl_delta_200_300": float(i * i),
            "satellite_delta_50_75": float(i * 3 + (i % 2)),
            "redshift": 0.3 + i * 0.01,
            "richness": 10.0 + i * 2,
            "target_region": 0,
            "control_region": 0,
        } for i in range(8)]
        with self.assertRaises(ValueError):
            partial_correlation(rows, region_effects=("target_region",))


if __name__ == "__main__":
    unittest.main()
