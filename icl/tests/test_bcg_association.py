from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from review_bcg_association import peak_candidate, stable_candidate
from refine_bcg_member_masks import satellite_members


def fixture():
    objects=np.zeros(2,dtype=[(key,'f8') for key in ('x','y','xpeak','ypeak','xcpeak','ycpeak')])
    for key in objects.dtype.names:objects[key]=[20,80]
    objects['x'][0]=35  # A shifted moment centroid must not move the peak.
    seg=np.zeros((100,100),dtype=int);seg[15:40,15:40]=1;seg[75:85,75:85]=2
    return objects,seg,np.zeros(seg.shape,bool)


def test_shifted_centroid_and_stable_peak():
    objects,seg,bad=fixture()
    assert peak_candidate(objects,seg,bad,(20,20),.168)==(0,'peak_segment_candidate')
    assert stable_candidate(objects[0],objects[0],.168)
    other=objects[0].copy();other['xpeak']+=3
    assert not stable_candidate(objects[0],other,.168)


def test_invalid_center_and_segment_disagreement_block():
    objects,seg,bad=fixture();bad[20,20]=True
    assert peak_candidate(objects,seg,bad,(20,20),.168)[1]=='center_invalid'
    bad[:]=False;seg[20,20]=2
    assert peak_candidate(objects,seg,bad,(20,20),.168)[1]=='central_segment_disagrees'


def test_multiple_or_filtered_peak_disagreement_block():
    objects,seg,bad=fixture()
    objects[1]=objects[0]
    assert peak_candidate(objects,seg,bad,(20,20),.168)[1]=='multiple_peaks'
    objects['xcpeak']=90
    assert peak_candidate(objects,seg,bad,(20,20),.168)[1]=='no_peak'


def test_bcg_member_exclusion_requires_unique_confirmed_position():
    members=np.array([[20.,20.],[40.,40.]])
    kept,n=satellite_members(members,(20,20),.168,True)
    assert n==1
    np.testing.assert_array_equal(kept,[[40.,40.]])
    assert satellite_members(members,(20,20),.168,False)[1]==0
    assert satellite_members(np.array([[20.,20.],[21.,20.]]),(20,20),.168,True)[1]==0


def test_member_fallback_excludes_bcg_but_preserves_neighbor_overlap():
    import prepare_detection_psf as prep
    objects=np.array([(35.,20.,2.,2.,0.)],dtype=[(k,'f8') for k in ('x','y','a','b','theta')])
    seg=np.zeros((100,100),dtype=int);seg[15:40,15:40]=1;base=np.zeros(seg.shape,bool)
    members=np.array([[20.,20.]])
    kept,_=satellite_members(members,(20,20),.168,True)
    masks,_=prep.source_masks(objects,seg,np.ones(1),np.zeros(1),base,0,kept,.168)
    assert not masks[0][20,20]
    # A real nearby member remains excluded even when its mask overlaps the BCG.
    kept,_=satellite_members(np.array([[20.,20.],[20.,30.]]),(20,20),.168,True)
    masks,_=prep.source_masks(objects,seg,np.ones(1),np.zeros(1),base,0,kept,.168)
    assert masks[0][20,20]
