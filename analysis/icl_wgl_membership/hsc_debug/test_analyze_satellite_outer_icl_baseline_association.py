import unittest

from analyze_satellite_outer_icl_baseline_association import (
    SATELLITE_RADIAL_BINS, expand_common_radial_profiles, join_rows, metrics,
)


def fixture():
    clusters = {"1": {"z": 0.2, "ra_deg": 12.0, "dec_deg": -1.0, "richness": 20.0}}
    controls = {"1": {"control_id": "c1", "ra_deg": 12.1, "dec_deg": -1.0}}
    icl = [{"target_id": "1", "control_id": "c1", "filter": "i",
            "z": "0.2", "ra_deg": "12.0", "dec_deg": "-1.0"}]
    outer = [{"filter": "i", "r_min_kpc": "200", "r_max_kpc": "300",
              "target_id": "1", "control_id": "c1", "inner_delta_njy_arcsec2": "4",
              "outer_400_500_delta_njy_arcsec2": "1",
              "outer_subtracted_delta_njy_arcsec2": "3", "inner_coverage": "0.5",
              "outer_coverage": "0.5"}]
    satellite = [{"cluster_id": "1", "filter": "i", "cohort": "full_sample",
                  "radial_bin_kpc": "50-75", "mean_delta_sb_njy_arcsec2": "2",
                  "n_satellites": "3"}]
    return icl, outer, satellite, clusters, controls


class AssociationJoinTests(unittest.TestCase):
    def test_joins_same_pair_and_subtracts_its_outer_baseline(self):
        rows = join_rows(*fixture())["i"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["icl_raw"], 4.0)
        self.assertEqual(rows[0]["outer_baseline"], 1.0)
        self.assertEqual(rows[0]["icl_outer_subtracted"], 3.0)
        self.assertEqual(rows[0]["satellite"], 2.0)

    def test_other_radial_bins_and_filters_are_not_joined(self):
        icl, outer, satellite, clusters, controls = fixture()
        satellite.append({**satellite[0], "filter": "i2"})
        outer.append({**outer[0], "r_min_kpc": "100", "r_max_kpc": "200"})
        joined = join_rows(icl, outer, satellite, clusters, controls)
        self.assertEqual(len(joined["i"]), 1)
        self.assertEqual(joined["i2"], [])

    def test_rejects_inconsistent_subtraction_or_pair_mapping(self):
        icl, outer, satellite, clusters, controls = fixture()
        outer[0]["outer_subtracted_delta_njy_arcsec2"] = "2.5"
        with self.assertRaisesRegex(ValueError, "subtraction mismatch"):
            join_rows(icl, outer, satellite, clusters, controls)

    def test_rejects_outer_baseline_below_inherited_support_cut(self):
        icl, outer, satellite, clusters, controls = fixture()
        outer[0]["outer_coverage"] = "0.049"
        with self.assertRaisesRegex(ValueError, "5% support cut"):
            join_rows(icl, outer, satellite, clusters, controls)

    def test_partial_metrics_report_raw_and_corrected_separately(self):
        rows = [
            {"icl_raw": float(i), "icl_outer_subtracted": float(i * 2),
             "satellite": float(i + (i % 2)), "redshift": 0.2 + i * 0.01,
             "richness": 10.0 + i, "n_satellites": i,
             "inner_coverage": 0.2 + i * 0.01, "outer_coverage": 0.3 + i * 0.01}
            for i in range(1, 8)
        ]
        result = metrics(rows)
        self.assertIn("spearman_rho", result["icl_raw"])
        self.assertIn("spearman_rho", result["icl_outer_subtracted"])
        self.assertIn("partial_pearson_r_control_z_richness_satellite_count_area_support",
                      result["icl_raw"])

    def test_radial_expansion_keeps_one_fixed_cluster_set(self):
        base_rows = {"i": [{"cluster_id": "1", "filter": "i", "icl_raw": 4.0}], "i2": []}
        profiles = [{"cluster_id": "1", "filter": "i", "cohort": "full_sample",
                     "radial_bin_kpc": radial_bin,
                     "mean_delta_sb_njy_arcsec2": str(index),
                     "n_satellites": str(index + 1)}
                    for index, radial_bin in enumerate(SATELLITE_RADIAL_BINS)]
        expanded = expand_common_radial_profiles(base_rows, profiles)
        self.assertEqual(len(expanded["i"]), 7)
        self.assertEqual({row["cluster_id"] for row in expanded["i"]}, {"1"})
        self.assertEqual([row["radial_bin_kpc"] for row in expanded["i"]],
                         list(SATELLITE_RADIAL_BINS))
        self.assertEqual(expanded["i2"], [])


if __name__ == "__main__":
    unittest.main()
