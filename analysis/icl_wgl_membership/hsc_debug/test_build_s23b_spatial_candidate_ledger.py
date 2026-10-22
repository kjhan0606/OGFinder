import unittest

from build_s23b_spatial_candidate_ledger import (
    angular_sep_arcsec, build_ledger, sky_cell,
)


class SpatialCandidateLedgerTests(unittest.TestCase):
    @staticmethod
    def fixture():
        parent = []
        for i in range(800):
            parent.append({
                "id": str(i),
                "ra_deg": str((i % 60) * 5 + 0.25),
                "dec_deg": str(-87.5 + (i // 60) * 10),
            })
        association = []
        for i in range(26):
            association.append({"cluster_id": str(i), "filter": "i"})
        for i in range(100, 168):
            association.append({"cluster_id": str(i), "filter": "i2"})
        candidates = [{
            "id": "0",  # Deliberate numeric collision with a different-release ID.
            "name": "S23B-candidate",
            "ra_deg": 181.0,
            "dec_deg": 44.0,
            "z": 0.4,
            "richness": 30.0,
        }]
        return candidates, parent, association

    def test_release_id_collision_is_not_a_coordinate_match(self):
        candidates, parent, association = self.fixture()
        rows, stats = build_ledger(candidates, parent, association, {"0"})
        self.assertEqual(rows[0]["candidate_key"], "S23B:0")
        self.assertFalse(rows[0]["coordinate_match_parent_within_1arcsec"])
        self.assertEqual(stats["parent_coordinate_matches_within_1arcsec"], 0)
        self.assertTrue(rows[0]["center_inside_existing_control_cutout"])

    def test_spatial_cell_matches_covariance_helper(self):
        self.assertEqual(sky_cell(11.117649, 1.836322), "2_18")
        self.assertAlmostEqual(angular_sep_arcsec(359.999, 0., 0.001, 0.), 7.2, places=3)


if __name__ == "__main__":
    unittest.main()
