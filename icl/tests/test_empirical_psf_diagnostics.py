from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from diagnose_empirical_psf import measure, coordinates


def data():
    r,_=coordinates((401,401),200.2,199.7)
    return r,1000*np.exp(-r*r/8)+7,np.ones(r.shape),np.zeros(r.shape,bool)


def test_constant_sky_offset_invariance_and_finite_core_normalization():
    r,image,var,mask=data()
    s,p=measure(image,var,mask,200.2,199.7,(120,150))
    t,q=measure(image+20,var,mask,200.2,199.7,(120,150))
    assert s['status']==t['status']=='diagnostic_only'
    np.testing.assert_allclose(s['core_flux'],np.sum((image-7)[r<12]))
    np.testing.assert_allclose([a['mean_per_core_flux'] for a in p],[a['mean_per_core_flux'] for a in q],atol=1e-15)
    np.testing.assert_allclose(sum(a['observed_ring_flux_per_core_flux'] for a in p if a['r_max_pix']<=12),1.)


def test_masked_core_and_incomplete_background_fail():
    r,image,var,mask=data();mask[r<3]=True
    s,p=measure(image,var,mask,200.2,199.7,(120,150))
    assert s['status']=='blocked_core_mask' and p==[]
    mask[:]=False;mask[:200]=True
    s,p=measure(image,var,mask,200.2,199.7,(120,150))
    assert s['status']=='blocked_background_coverage' and p==[]


def test_background_curvature_changes_signed_wings_without_clipping():
    r,image,var,mask=data();image+=r*r*1e-4
    s,p=measure(image,var,mask,200.2,199.7,(120,150))
    t,q=measure(image,var,mask,200.2,199.7,(160,190))
    assert s['status']==t['status']=='diagnostic_only'
    assert t['background']>s['background'] and t['core_flux']<s['core_flux']
    assert p[-1]['mean_per_core_flux']<0 and q[-1]['mean_per_core_flux']<p[-1]['mean_per_core_flux']


def test_masked_contaminant_does_not_enter_profile_or_background():
    r,image,var,mask=data();mask[210:230,225:245]=True
    s,p=measure(image,var,mask,200.2,199.7,(120,150))
    image[mask]+=1e7
    t,q=measure(image,var,mask,200.2,199.7,(120,150))
    assert s==t and p==q
