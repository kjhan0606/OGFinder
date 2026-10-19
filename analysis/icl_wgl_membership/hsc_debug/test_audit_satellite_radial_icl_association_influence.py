import unittest

from audit_satellite_radial_icl_association_influence import analyze_group, summarize_deletions


def rows(n=12):
    return [{
        "cluster_id": str(i), "filter": "i", "radial_bin_kpc": "5-10",
        "icl_raw": float(i), "outer_baseline": float(i % 3) - 1,
        "icl_outer_subtracted": float(i) - (float(i % 3) - 1),
        "satellite": float(i * i + 1), "n_satellites": i + 1,
        "redshift": 0.15 + 0.002 * i, "richness": 8.0 + i,
        "ra_deg": 4.0 + 20.0 * i, "dec_deg": -2.0 + 0.1 * i,
        "control_ra_deg": 5.0 + 20.0 * i, "control_dec_deg": -2.0 + 0.1 * i,
        "inner_coverage": 0.2 + 0.001 * i, "outer_coverage": 0.3 + 0.001 * i,
    } for i in range(1, n + 1)]


class RadialInfluenceAuditTests(unittest.TestCase):
    def test_deletes_each_block_and_cluster_once(self):
        summary, block_rows, cluster_rows = analyze_group(rows(), "i", "5-10")
        self.assertEqual(summary["n_clusters"], 12)
        self.assertEqual(summary["n_connected_blocks"], 12)
        self.assertEqual(len(block_rows), 12)
        self.assertEqual(len(cluster_rows), 12)
        self.assertIn("icl_outer_subtracted", summary["delete_one_block"]["metrics"])

    def test_summary_counts_sign_reversals(self):
        point = {"icl_raw": {"spearman_rho": 0.2},
                 "icl_outer_subtracted": {"spearman_rho": -0.2}}
        deleted = [
            {"icl_raw": {"spearman_rho": 0.1}, "icl_outer_subtracted": {"spearman_rho": 0.1}},
            {"icl_raw": {"spearman_rho": -0.1}, "icl_outer_subtracted": {"spearman_rho": -0.3}},
        ]
        result = summarize_deletions(point, deleted, "unit")
        self.assertEqual(result["metrics"]["icl_raw"]["spearman_rho"]["n_sign_reversals"], 1)
        self.assertFalse(result["metrics"]["icl_raw"]["spearman_rho"]["all_same_sign_as_full_sample"])


if __name__ == "__main__":
    unittest.main()
