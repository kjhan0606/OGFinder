#!/usr/bin/env python3
import unittest

from audit_satellite_icl_tract_region_influence import (
    endpoint_touch_counts,
    leave_one_region_out,
)


class SatelliteICLTractInfluenceTests(unittest.TestCase):
    def test_deletes_pairs_touching_region_and_recomputes_statistics(self):
        rows = []
        for i in range(18):
            icl = float(i * i + (i % 3) * 0.7)
            rows.append({
                "cluster_id": str(i),
                "icl_delta_200_300": icl,
                "satellite_delta": 0.4 * icl + ((i * 7) % 5) - 2.0,
                "redshift": 0.3 + i * 0.01,
                "richness": 15.0 + ((i * 11) % 13),
                "target_region": i % 3,
                "control_region": (i + 1) % 3,
            })
        result = leave_one_region_out(rows)
        self.assertEqual(result["n_pairs"], 18)
        self.assertEqual(len(result["leave_one_region_out"]), 3)
        for deletion in result["leave_one_region_out"]:
            expected_omitted = sum(
                row["target_region"] == deletion["region_id"]
                or row["control_region"] == deletion["region_id"]
                for row in rows
            )
            self.assertEqual(deletion["n_pairs_omitted"], expected_omitted)
            self.assertEqual(deletion["n_pairs_retained"] + expected_omitted, len(rows))
            self.assertEqual(len(deletion["omitted_cluster_ids"]), expected_omitted)
            self.assertEqual(deletion["status"], "complete_exploratory")

    def test_empty_sample_rejected(self):
        with self.assertRaises(ValueError):
            leave_one_region_out([])

    def test_endpoint_counts_are_pairwise_unions(self):
        rows = [
            {"filter": "i", "target_region": 0, "control_region": 1},
            {"filter": "i", "target_region": 1, "control_region": 1},
            {"filter": "i2", "target_region": 2, "control_region": 0},
        ]
        self.assertEqual(endpoint_touch_counts(rows), {
            "i": {"0": 1, "1": 2},
            "i2": {"0": 1, "2": 1},
        })


if __name__ == "__main__":
    unittest.main()
