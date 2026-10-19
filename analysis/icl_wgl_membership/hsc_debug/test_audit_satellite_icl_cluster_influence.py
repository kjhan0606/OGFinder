#!/usr/bin/env python3
import unittest

from audit_satellite_icl_cluster_influence import leave_one_cluster_out, partial_spearman


def make_rows(n=12):
    return [{
        "cluster_id": str(i),
        "icl_delta_200_300": float(i * i + (i % 3) * 0.3),
        "satellite_delta_50_75": float(0.5 * i * i + (i * 7) % 11),
        "redshift": 0.2 + (i % 7) * 0.03,
        "richness": 10 + (i * 5) % 19,
    } for i in range(n)]


class ClusterInfluenceTests(unittest.TestCase):
    def test_leave_one_out_covers_each_cluster_and_preserves_sample_size(self):
        result = leave_one_cluster_out(make_rows())
        self.assertEqual(result["n_clusters"], 12)
        self.assertEqual(len(result["records"]), 12)
        self.assertEqual({row["n_remaining"] for row in result["records"]}, {11})
        self.assertEqual({row["removed_cluster_id"] for row in result["records"]},
                         {str(i) for i in range(12)})

    def test_duplicate_cluster_ids_are_rejected(self):
        rows = make_rows()
        rows[1]["cluster_id"] = rows[0]["cluster_id"]
        with self.assertRaises(ValueError):
            leave_one_cluster_out(rows)

    def test_too_few_rows_are_rejected(self):
        with self.assertRaises(ValueError):
            leave_one_cluster_out(make_rows(5))

    def test_partial_spearman_is_invariant_to_monotone_rescaling(self):
        rows = make_rows()
        original = partial_spearman(rows)
        transformed = [{**row,
                        "icl_delta_200_300": row["icl_delta_200_300"] ** 3,
                        "satellite_delta_50_75": row["satellite_delta_50_75"] ** 3,
                        "redshift": row["redshift"] ** 3,
                        "richness": row["richness"] ** 3}
                       for row in rows]
        self.assertAlmostEqual(original, partial_spearman(transformed), places=14)


if __name__ == "__main__":
    unittest.main()
