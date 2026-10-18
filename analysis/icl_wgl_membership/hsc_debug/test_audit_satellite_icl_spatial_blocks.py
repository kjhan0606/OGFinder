import unittest

import numpy as np

from audit_satellite_icl_spatial_blocks import bootstrap, coordinate_labels, spearman_rho


def rows_for_test():
    rows = []
    target_ra = [1.0, 1.2, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    control_ra = [140.0, 160.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0]
    x = [2.0, -1.0, 5.0, 1.0, 7.0, 3.0, 8.0, 4.0]
    y = [3.0, 0.0, 6.0, 4.0, 9.0, 2.0, 10.0, 5.0]
    for i in range(8):
        rows.append({
            "icl_delta_200_300": str(x[i]),
            "satellite_delta_50_75": str(y[i]),
            "redshift": str(0.25 + 0.013 * i + 0.001 * i * i),
            "richness": str(12 + (i * 7) % 13 + 0.2 * i),
            "target_ra_deg": str(target_ra[i]),
            "target_dec_deg": str(-2.0 + 0.03 * i),
            "control_ra_deg": str(control_ra[i]),
            "control_dec_deg": str(1.0 + 0.02 * i),
        })
    return rows


class SpatialBlockAuditTests(unittest.TestCase):
    def test_coordinates_connect_pairs_sharing_either_sky_cell(self):
        labels = coordinate_labels(rows_for_test())
        self.assertEqual(labels[0], labels[1])
        self.assertEqual(len(np.unique(labels)), 7)

    def test_spearman_rank_handles_monotonic_data(self):
        rows = rows_for_test()
        for index, row in enumerate(rows):
            row["icl_delta_200_300"] = str(index)
            row["satellite_delta_50_75"] = str(index * 3)
        self.assertAlmostEqual(spearman_rho(rows, "satellite_delta_50_75"), 1.0)

    def test_spatial_bootstrap_is_reproducible_and_retains_spatial_method(self):
        rows = rows_for_test()
        labels = coordinate_labels(rows)
        first = bootstrap(rows, "satellite_delta_50_75", labels, 1000, 20261020)
        second = bootstrap(rows, "satellite_delta_50_75", labels, 1000, 20261020)
        self.assertEqual(first, second)
        self.assertEqual(first["method"], "connected_target_control_5deg_block_bootstrap")
        self.assertEqual(first["n_blocks"], 7)
        self.assertGreaterEqual(first["replicas_valid"], 950)
        self.assertEqual(first["replicas_valid"] + first["replicas_invalid"], 1000)


if __name__ == "__main__":
    unittest.main()
