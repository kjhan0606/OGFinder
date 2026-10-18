import importlib.util
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'
sys.path.insert(0,str(HERE))
spec=importlib.util.spec_from_file_location('prepare_detection_psf',HERE/'prepare_detection_psf.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def field():
    y,x=np.indices((256,256));image=np.random.default_rng(7).normal(100,1,(256,256))
    for cx,cy,amp in [(128,128,200),(65,65,150),(190,190,120)]:
        image+=amp*np.exp(-((x-cx)**2+(y-cy)**2)/8)
    return image.astype(np.float32)


def test_detection_preserves_input_and_masks_keep_defects_on_bcg():
    image=field();original=image.copy();bad=np.zeros(image.shape,bool);variance=np.ones(image.shape,np.float32)
    objects,seg,kron,flags,bcg,snr,info=m.detect(image,variance,bad,(128,128),np.empty((0,2)),.168)
    assert bcg>=0 and len(objects)>=3
    np.testing.assert_array_equal(image,original)
    masks,_=m.source_masks(objects,seg,kron,flags,bad,bcg,np.empty((0,2)),.168)
    assert not masks[0][128,128] and masks[0][65,65]
    assert np.all(~masks[0]|masks[1]) and np.all(~masks[1]|masks[2])
    bright=bad.copy();bright[128,128]=True
    masks,_=m.source_masks(objects,seg,kron,flags,bright,bcg,np.empty((0,2)),.168)
    assert all(mask[128,128] for mask in masks)


def test_no_matching_bcg_is_not_silently_kept_and_member_fallback_works():
    image=field();bad=np.zeros(image.shape,bool);members=np.array([[30.,200.]])
    objects,seg,kron,flags,bcg,_,info=m.detect(image,np.ones_like(image),bad,(10,10),members,.168)
    assert bcg==-1 and info['bcg_matches']==0
    masks,count=m.source_masks(objects,seg,kron,flags,bad,bcg,members,.168)
    assert count==1 and masks[0][200,30] and masks[0][128,128]


def test_psf_candidates_reject_members_and_masked_cores_without_fallback():
    image=field();bad=np.zeros(image.shape,bool);variance=np.ones_like(image)
    objects,seg,kron,flags,bcg,snr,_=m.detect(image,variance,bad,(128,128),np.empty((0,2)),.168)
    rows,stamps,profiles,mode=m.psf_candidates(objects,seg,image,variance,bad,snr,bcg,np.empty((0,2)),.168)
    assert rows and len(stamps)==len(rows) and profiles
    for star in rows:
        assert star['NUMBER']!=bcg+1
        assert np.isfinite(stamps[star['STAMP_INDEX']][50,50])
    fully_masked=np.ones(image.shape,bool)
    rows,stamps,_,_=m.psf_candidates(objects,seg,image,variance,fully_masked,snr,bcg,np.empty((0,2)),.168)
    assert rows==[] and stamps==[]
    members=np.column_stack([objects['x'],objects['y']])
    rows,stamps,_,_=m.psf_candidates(objects,seg,image,variance,bad,snr,bcg,members,.168)
    assert rows==[] and stamps==[]
