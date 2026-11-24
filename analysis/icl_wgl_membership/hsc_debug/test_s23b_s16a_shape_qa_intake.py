from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_s23b_s16a_shape_qa_pilot as qa


def test_radial_edges_and_final_cone_boundary_account_for_every_source():
    radius=np.array([0,.099,.1,.2,.4,.7,1,1.5,2,3,3.000001])
    assert qa.radial_bins(radius).tolist()==[0,0,1,2,3,4,5,6,7,7,7]
    with pytest.raises(ValueError):qa.radial_bins([3.01])
    with pytest.raises(ValueError):qa.radial_bins([float('nan')])


def test_effective_shape_weight_number_is_scaled_and_empty_support_is_explicit():
    assert qa.weight_diagnostics([])=={'sum_shape_weight':0.,'effective_number_shape_weight':0.}
    small=qa.weight_diagnostics([1,2,3])
    large=qa.weight_diagnostics(np.array([1,2,3])*1e200)
    assert small['effective_number_shape_weight']==pytest.approx(36/14)
    assert large['effective_number_shape_weight']==pytest.approx(36/14)
    assert qa.weight_diagnostics([10,10,10])['effective_number_shape_weight']==pytest.approx(3)
    for bad in ([0],[-1],[float('inf')]):
        with pytest.raises(ValueError):qa.weight_diagnostics(bad)


def test_angular_separation_handles_ra_wrap_without_unit_conversion():
    assert np.rad2deg(qa.separations_radians([.1],[0],359.9,0))[0]==pytest.approx(.2)


def test_known_zero_random_prevents_common_support_even_when_other_four_are_full():
    centers=[dict(center_id='L',matched_lens_id='A',role='lens',icl_filter='i')]
    centers += [dict(center_id=f'R{i}',matched_lens_id='A',role='random',icl_filter='i') for i in range(5)]
    profiles=[]
    for c in centers:
        for k in range(len(qa.RADIAL_EDGES)-1):
            profiles.append(dict(center_id=c['center_id'],radial_bin=k,
                raw_count=0 if c['center_id']=='R0' else 10,
                numeric_positive_weight_count=0 if c['center_id']=='R0' else 9))
    result=qa.common_support(profiles,centers)
    assert len(result)==8
    assert all(r['n_randoms_with_raw_rows']==4 and not r['raw_common_support'] for r in result)
    assert all(not r['numeric_weight_common_support'] for r in result)


def test_numeric_failure_does_not_remove_raw_radial_support():
    centers=[dict(center_id='L',matched_lens_id='A',role='lens',icl_filter='i')]
    centers += [dict(center_id=f'R{i}',matched_lens_id='A',role='random',icl_filter='i') for i in range(5)]
    profiles=[dict(center_id=c['center_id'],radial_bin=k,raw_count=1,
                   numeric_positive_weight_count=0 if c['center_id']=='L' else 1)
              for c in centers for k in range(8)]
    result=qa.common_support(profiles,centers)
    assert all(r['raw_common_support'] and not r['numeric_weight_common_support'] for r in result)


def test_source_reuse_and_large_integer_identity_are_kept_separate():
    ids=['L','R','B']
    objects=['9007199254740992','9007199254740992','9007199254740993']
    values=np.ones((3,10));values[1]=values[0]
    mapping={'L':dict(matched_lens_id='A',role='lens'),
             'R':dict(matched_lens_id='A',role='random'),
             'B':dict(matched_lens_id='B',role='lens')}
    result=qa.reuse_audit(ids,objects,values,mapping)
    assert result['unique_objects']==2
    assert result['objects_in_multiple_centers']==1
    assert result['objects_in_own_lens_and_random']==1
    assert result['inconsistent_source_properties']==0
    values[1,2]=2
    assert qa.reuse_audit(ids,objects,values,mapping)['inconsistent_source_properties']==1
    with pytest.raises(ValueError,match='duplicate'):
        qa.reuse_audit(['L','L','B'],objects,values,mapping)


def test_weight_identity_diagnoses_mismatch_without_selection_or_repair():
    result=qa.weight_identity([.1,.1,float('nan')],[.3,.3,.3],[10,20,10])
    assert result['checked_rows']==2 and result['excluded_rows']==1
    assert result['above_tolerance']==1
    assert result['used_as_selection_cut'] is False
