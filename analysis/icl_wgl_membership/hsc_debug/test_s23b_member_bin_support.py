"""Regression checks for member-bin decisions and nonexclusive failure causes."""
import json
import unittest
import numpy as np
from measure_s23b_camira_member_outer_profiles import bin_support


def profile(coverage):
    values = np.array(coverage, dtype=float)
    return dict(coverage=values, sb=np.where(values > 0, 3., np.nan),
                n_pixels=(values * 100).astype(int))


class BinSupportTests(unittest.TestCase):
    def test_decisions_match_legacy_rule_over_boundary_cases(self):
        # Independent legacy decision expression, including exact threshold and
        # bins that pass despite a different bin failing in the same profile.
        for cov in (0., .599999, .6, .600001, 1.):
            for n in range(5):
                target = profile([cov, 0.])
                controls = [profile([.6, 0.])] * n + [profile([.59, 1.])]
                expected = cov >= .6 and n >= 3
                result = bin_support(target, controls, 0, len(controls), .6, .6, 3)
                self.assertEqual(result['retained'], expected)
                self.assertEqual(result['n_valid_controls'], n)
                self.assertEqual(result['failure_reason'] == 'retained', expected)

    def test_overlapping_target_and_control_area_failures(self):
        result = bin_support(profile([.2]), [profile([.4])]*4, 0, 4, .6, .6, 3)
        self.assertEqual(set(result['failure_reason'].split('|')),
                         {'target_coverage_below_floor', 'control_bin_coverage_below_floor'})
        self.assertFalse(result['retained'])

    def test_control_pool_and_stamp_failures_are_distinct(self):
        target = profile([.8])
        a = bin_support(target, [profile([.9])]*2, 0, 2, .6, .6, 3)
        b = bin_support(target, [profile([.9])]*2, 0, 4, .6, .6, 3)
        self.assertEqual(a['failure_reason'], 'too_few_accepted_control_centers')
        self.assertEqual(b['failure_reason'], 'too_few_available_control_profiles')

    def test_unmatched_is_unmeasured(self):
        result = bin_support(None, [], 0, 0, .6, .6, 3, matched=False)
        self.assertEqual(result['failure_reason'], 'no_unique_sep_match')
        self.assertIsNone(result['target_coverage'])
        self.assertIsNone(result['n_valid_controls'])
        self.assertIsNone(result['control_count_pass'])

    def test_zero_pixel_flux_is_null_in_json(self):
        result = bin_support(profile([0.]), [profile([0.])], 0, 1, .6, .6, 3)
        self.assertIsNone(result['target_sb_njy_arcsec2'])
        control = json.loads(result['control_bin_measurements_json'])[0]
        self.assertIsNone(control['sb_njy_arcsec2'])
        json.dumps(result, allow_nan=False)

    def test_missing_target_is_recorded_with_control_support(self):
        result = bin_support(None, [profile([.8])]*3, 0, 3, .6, .6, 3)
        self.assertEqual(result['failure_reason'], 'target_profile_unavailable')
        self.assertTrue(result['control_count_pass'])
        self.assertFalse(result['retained'])


if __name__ == '__main__':
    unittest.main()
