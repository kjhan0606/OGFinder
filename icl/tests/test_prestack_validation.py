from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from prestack_core import geometry, measure, diffuse_model, covariance_diagnostics


def test_planar_sky_is_removed_and_shared_sky_fit_covariance_is_propagated():
    shape = (256, 256); center = (127.4, 128.2); kpp = 8.
    geom = geometry(shape, center, kpp)
    image = np.broadcast_to(2.+.3*geom[0]/1000-.5*geom[1]/1000, shape).copy()
    variance = np.ones(shape); mask = np.zeros(shape, bool)
    mask[:50, :100] = True; image[mask] = 1e12
    result = measure(image, variance, mask, geom, center, kpp)
    np.testing.assert_allclose(result['corrected'], 0, atol=1e-12)
    np.testing.assert_allclose(result['plane'], [2, .3, -.5], atol=1e-12)
    assert result['nominal_corrected'][0, 1] != 0
    assert np.linalg.eigvalsh(result['nominal_corrected']).min() > 0
    assert np.all(np.diag(result['nominal_corrected']) > np.diag(result['nominal_raw']))


def test_fixed_operator_linearity_and_extended_halo_self_subtraction():
    shape = (256, 256); center = (128., 128.); kpp = 8.
    geom = geometry(shape, center, kpp); var = np.ones(shape); mask = np.zeros(shape, bool)
    before = np.random.default_rng(3).normal(size=shape)
    model = diffuse_model(geom[0], geom[1], 250, .65, 3.)
    a = measure(before, var, mask, geom, center, kpp)
    b = measure(before+model, var, mask, geom, center, kpp)
    t = measure(model, var, mask, geom, center, kpp)
    for key in ('raw', 'corrected'):
        np.testing.assert_allclose(b[key]-a[key], t[key], atol=1e-12)
    assert np.all(t['corrected'] < t['raw'])
    assert t['corrected'][-1]/t['raw'][-1] < .9


def test_covariance_resampling_preserves_cross_bins_and_heterogeneous_nominal_errors():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(200, 3)); x[:, 1] += x[:, 0]*2
    labels = np.array([str(i//5) for i in range(200)])
    nominal = np.array([np.eye(3)*(1+i) for i in range(200)])
    summary, arrays = covariance_diagnostics(x, labels, nominal, replicas=4000)
    assert summary['rank'] == 3 and summary['covariance_sample_size_gate']
    np.testing.assert_allclose(arrays['independent_stack_covariance'], np.cov(x.T)/200, atol=1e-14)
    np.testing.assert_allclose(arrays['nominal_independent_stack_covariance'], nominal.sum(axis=0)/200**2)
    assert arrays['sample_covariance'][0, 1] > 1
    assert max(abs(np.array(summary['iid_bootstrap_diagonal_ratio'])-1)) < .1
    assert np.linalg.eigvalsh(arrays['spatial_bootstrap_covariance']).min() > 0


def test_missing_sky_support_and_insufficient_covariance_are_explicit():
    image = np.ones((128, 128)); geom = geometry(image.shape, (64,64), 2.)
    result = measure(image, image, np.zeros(image.shape, bool), geom, (64,64), 2.)
    assert result['n_tiles'] == 0 and np.isnan(result['corrected']).all()
    info, arrays = covariance_diagnostics(np.ones((1, 6)), ['one'], np.ones((1,6,6)))
    assert info['status'] == 'insufficient_complete_fields' and arrays == {}


def test_nominal_plane_propagation_predicts_independent_noise_monte_carlo():
    shape=(256,256);center=(127.4,128.2);kpp=8.
    geom=geometry(shape,center,kpp);var=np.ones(shape);mask=np.zeros(shape,bool)
    mask[80:120,100:140]=True
    reference=measure(np.zeros(shape),var,mask,geom,center,kpp)['nominal_corrected']
    rng=np.random.default_rng(19)
    samples=np.array([measure(rng.normal(size=shape),var,mask,geom,center,kpp)['corrected'] for _ in range(500)])
    whitened=np.linalg.solve(np.linalg.cholesky(reference),samples.T).T
    eigenvalues=np.linalg.eigvalsh(np.cov(whitened.T))
    assert eigenvalues.min()>.7 and eigenvalues.max()<1.3


def test_spatial_bootstrap_preserves_shared_cell_offsets():
    rng=np.random.default_rng(32);cells=rng.normal(size=(40,2))
    values=np.repeat(cells,10,axis=0)+rng.normal(0,.01,(400,2))
    labels=np.repeat(np.arange(40),10)
    _,matrices=covariance_diagnostics(values,labels,np.tile(np.eye(2),(400,1,1)),replicas=1500)
    ratio=np.diag(matrices['spatial_bootstrap_covariance'])/np.diag(matrices['iid_bootstrap_covariance'])
    assert np.all(ratio>7) and np.all(ratio<13)
