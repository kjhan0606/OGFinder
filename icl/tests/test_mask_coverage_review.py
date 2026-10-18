from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from review_mask_coverage import coverage,center_causes


def test_coverage_and_sectors_reveal_directional_loss():
    y,x=np.indices((101,101));dx=x-50;dy=y-50;r=np.hypot(dx,dy);angle=np.mod(np.arctan2(dy,dx),2*np.pi)
    mask=dx>=0
    rows=coverage(mask,r,angle,[10,20,30])
    for row in rows:
        assert .45<row['grid_coverage']<.5
        assert row['min_sector_coverage']==0
        assert row['sectors_below_half']==4
    assert coverage(np.zeros_like(mask),r,angle,[10,20])[0]['grid_coverage']==1
    assert coverage(np.ones_like(mask),r,angle,[10,20])[0]['grid_coverage']==0


def test_center_attribution_preserves_overlapping_reasons():
    objects=np.array([(20.,20.,2.,2.,0.),(23.,20.,2.,2.,0.)],dtype=[(k,'f8') for k in ('x','y','a','b','theta')])
    seg=np.zeros((50,50),int);seg[20,20]=1;base=np.zeros(seg.shape,bool);base[20,20]=True
    out=center_causes(objects,seg,np.ones(2),np.zeros(2),base,0,np.array([[20.,30.]]),(20.,20.),.168,2.5)
    assert out['base_center'] and not out['source_segment_center']
    assert out['neighbor_ellipse_numbers']==[2]
    assert out['fallback_members_covering_center']==1
