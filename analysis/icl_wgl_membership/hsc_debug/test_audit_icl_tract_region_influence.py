#!/usr/bin/env python3
import unittest

import numpy as np

from audit_icl_tract_region_influence import leave_one_region_out_means


class ICLTractRegionInfluenceTests(unittest.TestCase):
    def test_equal_pair_means_and_region_deletions(self):
        full, deletions = leave_one_region_out_means(
            [1.0, 3.0, 10.0, 14.0], [0, 0, 1, 2],
        )
        self.assertEqual(full, 7.0)
        self.assertEqual(deletions, [(0, 2, 12.0), (1, 3, 6.0), (2, 3, 4.666666666666667)])

    def test_rejects_invalid_vectors_and_single_region(self):
        with self.assertRaises(ValueError):
            leave_one_region_out_means([1.0, float("nan")], [0, 1])
        with self.assertRaises(ValueError):
            leave_one_region_out_means([1.0, 2.0], [0, 0])
        with self.assertRaises(ValueError):
            leave_one_region_out_means([1.0], [0])


if __name__ == "__main__":
    unittest.main()
