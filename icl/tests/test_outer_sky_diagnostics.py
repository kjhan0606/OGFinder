from pathlib import Path
import sys
import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from diagnose_outer_sky import tile_diagnostics


def test_plane_recovery_and_masked_extreme_pixel_does_not_leak():
    y,x = np.indices((128,128)); image = .2+.03*x-.02*y
    mask = np.zeros(image.shape,bool); mask[45,45] = True; image[45,45] = 1e20
    tiles,lags,result = tile_diagnostics(image,np.ones(image.shape),mask,np.ones(image.shape,bool),(0,0),1.,32)
    assert len(tiles)==16 and result['status']=='diagnostic_only'
    np.testing.assert_allclose([result['plane_constant'],result['plane_x_per_mpc'],result['plane_y_per_mpc']], [.2,30,-20], atol=1e-12)
    assert result['plane_residual_mad'] < 1e-12
    assert max(t['raw_mean'] for t in tiles) < 5


def test_independent_noise_and_smoothed_noise_have_distinct_pair_statistics():
    noise = np.random.default_rng(41).normal(size=(256,256)); var=np.ones(noise.shape)
    mask=np.zeros(noise.shape,bool); region=~mask
    _,white,_=tile_diagnostics(noise,var,mask,region,(128,128),1.,64)
    _,smooth,_=tile_diagnostics(gaussian_filter(noise,1.5),var,mask,region,(128,128),1.,64)
    assert all(abs(r['correlation']) < .04 for r in white if r['lag_pixels']==1)
    assert all(r['correlation'] > .8 for r in smooth if r['lag_pixels']==1)


def test_fully_masked_field_is_reported_without_fabricated_measurement():
    image=np.ones((64,64)); mask=np.ones(image.shape,bool)
    tiles,lags,result=tile_diagnostics(image,image,mask,mask,(32,32),1.)
    assert not tiles and result['status']=='insufficient_tiles'
    assert result['median_tile_mean'] is None and result['outer_valid_fraction']==0
    assert all(r['n_pairs']==0 and r['correlation'] is None for r in lags)


def test_nominal_aperture_error_uses_only_unmasked_variance():
    image=np.ones((64,64)); var=np.full(image.shape,4.); mask=np.zeros(image.shape,bool)
    mask[:8,:]=True; var[mask]=1e20
    tiles,_,_=tile_diagnostics(image,var,mask,np.ones(image.shape,bool),(32,32),1.)
    assert len(tiles)==1 and tiles[0]['n_valid']==56*64
    np.testing.assert_allclose(tiles[0]['nominal_mean_sigma'],2/np.sqrt(56*64))
