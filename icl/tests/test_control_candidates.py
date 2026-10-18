from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from prepare_control_candidates import xyz, angle, exclusion_counts, center_ok, qa_ok, BAD_FLAGS, select_nonoverlap


def test_spherical_geometry_handles_wrap_and_pole():
    np.testing.assert_allclose(angle(xyz([359.9], [0]), xyz([.1], [0])), [.2], atol=1e-10)
    assert angle(xyz([0], [89.9]), xyz([180], [89.9]))[0] < .201


def test_whole_footprint_and_non_nearest_large_cluster_are_excluded():
    points = xyz([0], [0]); clusters = xyz([.15, .3], [0, 0])
    radius = np.array([.02, .25])
    assert exclusion_counts(points, np.array([.1]), clusters, radius, 1.)[0] == 1
    assert exclusion_counts(points, np.array([0.]), clusters, radius, 1.)[0] == 0
    assert exclusion_counts(points, np.array([.1]), clusters, radius, 2.)[0] == 1


def test_catalogue_unknown_flags_and_nonfinite_quality_fail_closed():
    row = {k:'True' for k in ('isprimary', 'detect_ispatchinner', 'detect_istractinner')}
    row.update({'i_'+k:'False' for k in BAD_FLAGS}); row.update(i_inputcount_value='4', i_pix_variance='.01')
    assert center_ok(row)
    row['i_pixelflags_crcenter'] = ''; assert not center_ok(row)
    q = dict(icalexp_exists='True', imag_psf_depth='26', iseeing='.6')
    assert qa_ok(q)
    q['iseeing'] = 'nan'; assert not qa_ok(q)


def test_assignment_avoids_shared_id_and_overlapping_different_controls():
    def row(cid, oid, ra, score):
        return dict(target_cluster_id=cid, object_id=oid, ra=ra, dec=0., footprint_cap_deg=.1,
                    metadata_match=True, match_score=score, **{'overlap_1.0':0})
    rows = [row('a','same',0,0), row('b','same',0,0), row('b','near',.1,1), row('b','far',.3,2)]
    chosen = select_nonoverlap(rows, 1.)
    assert [(r['target_cluster_id'],r['object_id']) for r in chosen] == [('a','same'),('b','far')]
    rows[-1]['overlap_1.0'] = 1
    assert len(select_nonoverlap(rows, 1.)) == 1
