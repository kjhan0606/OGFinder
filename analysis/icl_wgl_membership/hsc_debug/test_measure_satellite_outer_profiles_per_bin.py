import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

from measure_satellite_outer_profiles_per_bin import summarize_per_bin


class PerAnnulusSatelliteSummaryTest(unittest.TestCase):
    def test_bins_keep_independent_support_and_weight_clusters_equally(self):
        rows = [
            {"cluster_id": "A", "filter": "i2", "cohort": "full_sample",
             "member_index": "1", "radial_bin_kpc": "5-10",
             "target_sb_njy_arcsec2": 3., "control_sb_njy_arcsec2": 3.,
             "delta_sb_njy_arcsec2": 0.},
            {"cluster_id": "A", "filter": "i2", "cohort": "full_sample",
             "member_index": "2", "radial_bin_kpc": "5-10",
             "target_sb_njy_arcsec2": 23., "control_sb_njy_arcsec2": 3.,
             "delta_sb_njy_arcsec2": 20.},
            {"cluster_id": "B", "filter": "i2", "cohort": "full_sample",
             "member_index": "3", "radial_bin_kpc": "5-10",
             "target_sb_njy_arcsec2": 103., "control_sb_njy_arcsec2": 3.,
             "delta_sb_njy_arcsec2": 100.},
            {"cluster_id": "A", "filter": "i2", "cohort": "full_sample",
             "member_index": "1", "radial_bin_kpc": "50-75",
             "target_sb_njy_arcsec2": 13., "control_sb_njy_arcsec2": 3.,
             "delta_sb_njy_arcsec2": 10.},
        ]
        with TemporaryDirectory() as temporary:
            summary = summarize_per_bin(rows, Path(temporary), replicas=100, seed=7)

        inner = summary["full_sample"]["i2"]["bins"]["5-10"]
        outer = summary["full_sample"]["i2"]["bins"]["50-75"]
        self.assertEqual(inner["n_clusters"], 2)
        self.assertEqual(inner["n_satellites"], 3)
        self.assertAlmostEqual(inner["mean_delta"], 55.)
        self.assertEqual(outer["n_clusters"], 1)
        self.assertEqual(outer["n_satellites"], 1)
        self.assertAlmostEqual(outer["mean_delta"], 10.)

    def test_empty_input_writes_valid_empty_summary_tables(self):
        with TemporaryDirectory() as temporary:
            summary = summarize_per_bin([], Path(temporary), replicas=100, seed=7)
            self.assertEqual(summary, {})
            self.assertEqual(len((Path(temporary) / "cluster_profiles.tsv").read_text().splitlines()), 1)
            self.assertEqual(len((Path(temporary) / "filter_profiles.tsv").read_text().splitlines()), 1)


if __name__ == "__main__":
    unittest.main()
