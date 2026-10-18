from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from prepare_extended_psf import radial_summary
from compare_extended_psf import radial_template,forward_case
from diagnose_empirical_psf import coordinates
from inventory_psf_support import full_stamp,source_path


def test_reference_requires_centered_normalized_finite_image():
    image=np.zeros((401,401));image[200,200]=1
    rows,info=radial_summary(image)
    assert info['finite_grid_sum']==1 and info['center_pixel_xy']==[200,200]
    with pytest.raises(ValueError):radial_summary(image*2)
    with pytest.raises(ValueError):radial_summary(np.roll(image,1,axis=0))
    image[0,0]=np.nan
    with pytest.raises(ValueError):radial_summary(image)


def test_signed_reference_pixels_are_retained_and_crop_not_renormalized():
    image=np.zeros((401,401));image[200,200]=1.01;image[205,200]=-.01
    rows,info=radial_summary(image)
    assert info['negative_pixels']==1 and any(r['mean']<0 for r in rows)
    model=radial_template(rows,(31,31),15.2,14.8)
    assert model.min()<0 and not np.isclose(model.sum(),1)
    with pytest.raises(ValueError):radial_template(rows,(801,801),400,400)


def test_differential_injection_and_reference_sky_loss_with_real_background():
    r,_=coordinates((401,401),200.1,200.3)
    model=(1+r*r/9)**-1.6;model/=model.sum()
    variance=np.ones(r.shape);mask=np.zeros(r.shape,bool)
    image=200000*model+3+.0001*r*r
    measured,profiles,added=forward_case(image,variance,mask,model,200.1,200.3,(120,150))
    assert measured['status']=='diagnostic_only' and len(added)==33
    assert 0<measured['reference_core_transfer']<1
    assert 0<profiles[-1]['reference_sky_removal_fraction']<1
    assert max(r['absolute_difference'] for r in added)<1e-10


def test_incomplete_sky_does_not_produce_comparison_or_injection():
    r,_=coordinates((401,401),200,200);model=np.exp(-r*r/10);model/=model.sum()
    mask=r>100
    measured,profiles,added=forward_case(1e6*model,np.ones(r.shape),mask,model,200,200,(120,150))
    assert measured['status']=='blocked_background_coverage' and profiles==added==[]


def test_native_stamp_support_includes_last_pixel_and_rejects_off_edge():
    assert full_stamp(200,200,(401,401))
    assert full_stamp(199.6,200.4,(401,401))
    assert not full_stamp(199.4,200,(401,401))
    assert not full_stamp(200,200.6,(401,401))
    assert not full_stamp(np.nan,200,(401,401))


def test_downloader_relative_image_names_resolve_under_images():
    manifest=Path('/data/run/download_manifest.json')
    assert source_path(manifest,'cluster_7.fits')==Path('/data/run/images/cluster_7.fits')
    assert source_path(manifest,'/data/old/cluster_7.fits')==Path('/data/old/cluster_7.fits')
    with pytest.raises(ValueError):source_path(manifest,'../cluster_7.fits')
