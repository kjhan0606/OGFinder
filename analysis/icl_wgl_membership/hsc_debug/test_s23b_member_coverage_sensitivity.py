"""Check fixed-candidate support and paired outcome differences independently."""
import json
import unittest
import numpy as np
from analyze_s23b_member_coverage_sensitivity import candidate_profiles, cluster_profiles, fit


def row(member, target_cov, control_coverages, target=20.):
    return dict(target_id='A',filter='i',member_index=str(member),radial_bin_kpc='50-75',
        sep_matched='True',target_profile_available='True',target_coverage=str(target_cov),
        target_sb_njy_arcsec2=str(target),control_bin_measurements_json=json.dumps([
            dict(coverage=c,sb_njy_arcsec2=v) for c,v in control_coverages]))


class FixedSupportTests(unittest.TestCase):
    def test_nested_members_and_control_mean_changes(self):
        controls=[(.9,1.),(.9,2.),(.9,3.),(.5,10.)]
        pool=candidate_profiles([row(1,.9,controls),row(2,.7,controls,target=50.)])
        self.assertEqual(len(pool[.4]),2)
        self.assertEqual(len(pool[.8]),1)
        self.assertAlmostEqual(pool[.4][('A','i','1')]['delta'],16.)
        self.assertAlmostEqual(pool[.8][('A','i','1')]['delta'],18.)
        fixed=cluster_profiles(pool,'common_candidates')
        self.assertEqual([fixed[f][('A','i')]['n_member_positions_50_75'] for f in (.4,.6,.8)],[1,1,1])
        self.assertAlmostEqual(fixed[.4][('A','i')]['member_outer_delta_njy_arcsec2'],16.)
        self.assertAlmostEqual(fixed[.8][('A','i')]['member_outer_delta_njy_arcsec2'],18.)
        self.assertEqual(fixed[.4][('A','i')]['median_member_control_coverage_50_75'],.9)
        native=cluster_profiles(pool,'native')
        self.assertEqual(native[.4][('A','i')]['n_member_positions_50_75'],2)

    def test_common_clusters_does_not_mean_common_members(self):
        pool=candidate_profiles([row(1,.9,[(.9,1.)]*3),row(2,.5,[(.9,2.)]*3)])
        common=cluster_profiles(pool,'common_clusters')
        self.assertEqual(common[.4][('A','i')]['n_member_positions_50_75'],2)
        self.assertEqual(common[.8][('A','i')]['n_member_positions_50_75'],1)

    def test_duplicate_candidate_is_rejected(self):
        r=row(1,.9,[(.9,1.)]*3)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            candidate_profiles([r,r])

    def test_paired_difference_slope_uses_same_design(self):
        rng=np.random.default_rng(17)
        rows=[]
        for i in range(40):
            x=float(rng.normal())
            rows.append(dict(target_id=str(i),ra_deg=8.*i,dec_deg=10.,control_ra_deg=8.*i+1.,control_dec_deg=10.,
                member_outer_delta_njy_arcsec2=2*x,icl_contrast_njy_arcsec2=x,richness=20.+rng.uniform(0,20),
                z_cl=.3+rng.uniform(0,.2),min_icl_common_support_fraction=rng.uniform(.2,.5),
                median_member_target_coverage_50_75=.9,median_member_control_coverage_50_75=.9))
        regular=fit(rows,2.5,0.,100,17)
        paired=fit(rows,2.5,0.,100,17,difference=[.5*r['icl_contrast_njy_arcsec2'] for r in rows])
        self.assertEqual(regular['status'],'complete')
        self.assertAlmostEqual(regular['slope'],2.)
        self.assertAlmostEqual(paired['slope'],.5)


if __name__=='__main__': unittest.main()
