from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from connect_control_patches import parse_patches, contains, Geometry, transport_manifest
from prepare_control_candidates import xyz
from download_pilot import estimate_bytes, build_manifest_row


def polygon(tract=1):
    coords=[(359.9,-.1),(.1,-.1),(.1,.1),(359.9,.1),(359.9,-.1)]
    return f'Tract: {tract} Patch: 0,0 Center (RA, Dec): (0, 0)\n'+'\n'.join(
        f'Tract: {tract} Patch: 0,0 Corner{i} (RA, Dec): ({ra}, {dec})' for i,(ra,dec) in enumerate(coords))


def test_official_parser_wrap_boundary_and_hemisphere():
    patch=parse_patches(polygon())['10000']
    assert contains(patch,xyz([0],[0])[0])
    assert not contains(patch,xyz([180],[0])[0])
    assert not contains(patch,xyz([.2],[0])[0])
    assert contains(patch,xyz([.1],[0])[0])


def test_missing_or_duplicate_vertices_fail_closed():
    with pytest.raises(ValueError): parse_patches(polygon().split('Corner4')[0])
    with pytest.raises(ValueError): parse_patches(polygon()+'\n'+polygon())
    with pytest.raises(ValueError): parse_patches('not a geometry table')


def test_tract_overlap_is_preserved_instead_of_selecting_first():
    geometry=Geometry(parse_patches(polygon(1)+'\n'+polygon(2)))
    assert set(geometry.locate([0],[0])[0])=={'10000','20000'}
    assert geometry.locate([1],[0])==[[]]


def test_acquisition_contract_keeps_control_identity_and_target_scale():
    row=dict(target_cluster_id='7',object_id='12345',ra=3.,dec=2.,assigned_z=.4,
             requested_half_width_arcsec=200.,target_patch_ids='10000',skymap_id='10000',ebv_sfd=.03)
    plan=transport_manifest([row],dict(h0=70.,omega_m=.3,requested_outer_radius_kpc=1000.))
    frame=plan['clusters'][0]
    assert plan['schema']=='hsc-bcg-diffuse-sample-v1' and plan['bcg_policy']=='none'
    assert frame['frame_role']=='matched_control' and frame['z']==.4 and frame['ra_deg']==3.
    assert 'image' not in frame and not plan['science_measurement']
    assert estimate_bytes(plan['clusters'])[0]>0
    review=dict(mask_hdu=2,science_hdu=1,variance_hdu=3,coadd_header_filter='i',input_unit='instrumental_flux/pixel',ab_zeropoint=27.)
    saved=build_manifest_row(frame,'mock.fits',review)
    assert saved['target_cluster_id']=='7' and saved['bcg_policy']=='none' and saved['frame_role']=='matched_control'
    assert saved['calibration']['a_i_mag'] is None
