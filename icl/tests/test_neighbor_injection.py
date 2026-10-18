from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from inject_neighbor_recovery import profile,fit_operator,estimate,design
from summarize_neighbor_recovery import transfer


def test_profiles_conserve_finite_stamp_flux_and_psf_wing_changes_shape():
    core=profile((81,81),(40.3,40.1),3)
    wing=profile((81,81),(40.3,40.1),3,wing=True)
    assert abs(core.sum()-1)<1e-12 and abs(wing.sum()-1)<1e-12
    assert wing[40,40]<core[40,40] and np.all(core>=0)


def test_known_template_amplitudes_recovered_without_noise():
    core,diffuse,neighbors,halos,sky=design((101,101),(50,50),np.array([[65.,52.]]),1.)
    cols=[core,diffuse]+neighbors+halos+sky;truth=np.array([4000.,1500.,2000.,800.,.03,.02,-.01])
    image=sum(t*c for t,c in zip(truth,cols));before=image.copy()
    mask=np.zeros(image.shape,bool);mask[20:25,:]=True
    op,status=fit_operator(cols,np.ones(image.shape),mask)
    assert status=='ok';np.testing.assert_allclose(estimate(op,image),truth,atol=1e-7)
    np.testing.assert_array_equal(image,before)
    assert not op['valid'][22,50]


def test_degenerate_or_empty_fit_is_blocked():
    col=np.ones((20,20));var=col.copy();mask=np.zeros(col.shape,bool)
    assert fit_operator([col,col],var,mask)[1]=='ill_conditioned'
    assert fit_operator([col],var,~mask)[1]=='too_few_pixels'


def test_conditional_covariance_scales_with_variance_and_preserves_sign():
    x=np.indices((20,20))[1]/20;cols=[np.ones(x.shape),x]
    op,_=fit_operator(cols,np.ones(x.shape),np.zeros(x.shape,bool))
    noisier,_=fit_operator(cols,np.full(x.shape,4.),np.zeros(x.shape,bool))
    np.testing.assert_allclose(noisier['covariance'],4*op['covariance'])
    np.testing.assert_allclose(estimate(op,-2+3*x),[-2,3],atol=1e-12)


def test_transfer_recovers_cross_talk_and_rejects_survivor_only_grid():
    rows=[]
    for d in [1500,4500]:
        for h in [0,1000,2000]:
            rows.append(dict(status='ok',true_diffuse_flux=str(d),true_neighbor_halo_flux=str(h),
                diffuse_estimate=str(100+.8*d+.2*h),neighbor_halo_estimate=str(-20+.1*d+.9*h)))
    tr=transfer(rows)
    np.testing.assert_allclose([tr['diffuse_gain'],tr['halo_to_diffuse_leakage'],tr['diffuse_to_halo_leakage'],tr['halo_gain']],[.8,.2,.1,.9],atol=1e-12)
    rows[0]['status']='blocked';assert transfer(rows) is None
