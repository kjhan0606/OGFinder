"""The named-data gate counts caller-supplied rows and does not open a sky measurement."""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from icl.named_gate import FAILURE_MODES, named_data_gate  # noqa: E402


def _box(ra_min, ra_max, dec_min=0.0, dec_max=10.0):
    return {'ra_min': ra_min, 'ra_max': ra_max, 'dec_min': dec_min, 'dec_max': dec_max}


def _image(identity, sky, ra_min, ra_max, product='deepCoadd', version='HSC-SSP PDR3'):
    return {
        'id': identity,
        'version': version,
        'product': product,
        'sky': sky,
        'footprint': _box(ra_min, ra_max),
    }


def _cluster(identity, ra, dec, z, version='CAMIRA'):
    return {'id': identity, 'version': version, 'ra': ra, 'dec': dec, 'z': z}


def _shape(identity, ra, dec, z, quality=1, crossmatch=True, **extra):
    row = {
        'id': identity,
        'ra': ra,
        'dec': dec,
        'z': z,
        'spec_quality': quality,
        'crossmatch_ok': crossmatch,
    }
    row.update(extra)
    return row


def _shapes(rows, ra_min=0.0, ra_max=8.0, version='S19A', quality=True, crossmatch=True):
    return {
        'version': version,
        'footprint': _box(ra_min, ra_max),
        'rows': rows,
        'spec_quality_supplied': quality,
        'crossmatch_supplied': crossmatch,
    }


def _catalogues(rows=None, clusters=None, images=None, **shapes_kw):
    if images is None:
        images = [
            _image('global', 'global', 0.0, 10.0),
            _image('local', 'local', 20.0, 30.0),
            _image('calib', 'global', 40.0, 50.0, product='calexp'),
        ]
    if clusters is None:
        clusters = [
            _cluster('A', 5.0, 5.0, 0.25),
            _cluster('E', 6.0, 6.0, 0.5),
            _cluster('B', 25.0, 5.0, 0.25),
            _cluster('C', 15.0, 5.0, 0.25),
            _cluster('D', 5.0, 5.0, 0.25, version='redMaPPer'),
            _cluster('F', 9.0, 5.0, 0.25),
        ]
    if rows is None:
        # Dyadic redshifts keep the cluster window and the lens gap exact in binary.
        rows = [
            _shape('s1', 5.1, 5.1, 1.0),
            _shape('s2', 5.2, 5.2, 0.3125),
            _shape('s3', 5.3, 5.0, 0.0),
            _shape('s4', 5.4, 5.0, float('nan'), z_phot=0.9),
            _shape('s5', 5.5, 5.0, 1.0, quality=0),
            _shape('s6', 5.6, 5.0, 1.0, crossmatch=False),
            _shape('s7', 4.0, 4.0, 0.5),
            _shape('s8', 25.0, 5.0, 1.0),
            _shape('s9', 5.7, 5.2, 0.75, quality=True),
            _shape('s10', 8.0, 5.0, 1.0),
            _shape('s11', 4.2, 4.2, 1.0, z_is_spec=False),
            _shape('s12', 9.0, 5.0, 1.0),
            _shape('s13', 4.5, 4.5, 1.0, star=True),
        ]
    return images, clusters, _shapes(rows, **shapes_kw)


def test_empty_gate_is_not_a_sky_measurement():
    gate = named_data_gate()
    assert gate['status'] == 'empty'
    assert gate['n_lens_sources'] is None
    assert gate['sky_measurement'] is False
    assert gate['science_measurement'] is False
    assert gate['named_data_mass'] is False
    assert gate['survey_calibration'] is False
    assert gate['observed_sentence'] is None
    assert gate['observed_sentence_allowed'] is False
    assert gate['combination_matched'] is False
    assert gate['checklist_filled'] is False
    assert gate['design_scale']['detection'] is False
    assert 'not an individual-image point-source limit' in gate['design_scale']['note']
    assert 'before the spectroscopic cut' in gate['design_scale']['note']
    assert gate['frontier_fields']['use'] == 'centre mask and blending example'
    assert gate['frontier_fields']['stands_in_for_wide_radius_icl_mean'] is False
    assert gate['frontier_fields']['stands_in_for_representative_nfw'] is False
    assert gate['shape_calibration']['method'] == 're-Gaussianization'
    assert gate['shape_calibration']['same_as_measure_stamps'] is False
    roles = {item['data']: item['role'] for item in gate['other_data']}
    assert 'low-redshift ICL comparison' in roles.values()
    assert any('redMaPPer' in name for name in roles)
    assert any('KiDS' in name for name in roles)
    assert any('DES DR2' in name for name in roles)
    for key in FAILURE_MODES:
        assert gate['failure_modes'][key]['filled'] is False
        assert gate['failure_modes'][key]['statement']
    assert gate['reinforcement']['bins'] == {
        'brightness': None,
        'colour': None,
        'surface_brightness': None,
        'radius': None,
    }
    assert gate['reinforcement']['n_rows'] is None


def test_common_footprint_counts_only_spectroscopic_lens_sources():
    images, clusters, shapes = _catalogues()
    gate = named_data_gate(images, clusters, shapes, dz=0.125, z_gap=0.25)
    assert gate['status'] == 'counted'
    assert gate['combination_matched'] is True
    assert gate['n_lens_sources'] == 3
    assert gate['sky_measurement'] is False
    assert gate['science_measurement'] is False
    assert gate['observed_sentence'] is None
    assert gate['n_shapes_in_common_footprint'] == 11
    assert gate['n_shapes_outside'] == 2
    assert gate['n_rejected_quality'] == 1
    assert gate['n_rejected_crossmatch'] == 1
    assert gate['n_rejected_no_spec_z'] == 2
    assert gate['n_rejected_membership'] == 4
    assert (gate['n_lens_sources'] + gate['n_rejected_quality'] + gate['n_rejected_crossmatch']
            + gate['n_rejected_no_spec_z'] + gate['n_rejected_membership']
            == gate['n_shapes_in_common_footprint'])
    by_id = {row['id']: row for row in gate['per_cluster']}
    assert by_id['A']['source_ids'] == ['s1', 's9', 's10']
    assert by_id['E']['source_ids'] == ['s1', 's10']
    assert gate['images_accepted'] == ['global']
    reasons = {row['id']: row['reason'] for row in gate['images_rejected']}
    assert 'local-sky' in reasons['local']
    assert 'deepCoadd' in reasons['calib']
    excluded = {row['id']: row['reason'] for row in gate['clusters_excluded']}
    assert 'local-sky' in excluded['B']
    assert 'global-sky' in excluded['C']
    assert 'CAMIRA' in excluded['D']
    assert 'S19A' in excluded['F']
    assert 'not a sky measurement' in gate['note']
    assert 'science_measurement stays false' in gate['note']


def test_local_sky_alone_does_not_form_the_common_footprint():
    images = [_image('local', 'local-sky', 0.0, 10.0)]
    clusters = [_cluster('B', 5.0, 5.0, 0.25)]
    shapes = _shapes([_shape('s1', 5.0, 5.0, 1.0)], ra_max=10.0)
    gate = named_data_gate(images, clusters, shapes, dz=0.125, z_gap=0.25)
    assert gate['status'] == 'no_common_footprint'
    assert gate['n_lens_sources'] is None
    assert gate['combination_matched'] is False
    assert gate['science_measurement'] is False
    assert 'local-sky' in gate['clusters_excluded'][0]['reason']


def test_window_boundary_stays_a_lens_source_when_it_is_behind_the_gap():
    # |z - z_cl| = dz is outside the satellite window. z_gap is narrower than dz,
    # so that boundary redshift is still behind the lens.
    images = [_image('global', 'global-sky', 0.0, 10.0)]
    clusters = [_cluster('A', 5.0, 5.0, 0.25)]
    rows = [
        _shape('boundary', 5.0, 5.0, 0.5),
        _shape('satellite', 5.1, 5.0, 0.4375),
    ]
    gate = named_data_gate(images, clusters, _shapes(rows, ra_max=10.0), dz=0.25, z_gap=0.125)
    assert gate['n_lens_sources'] == 1
    assert gate['per_cluster'][0]['source_ids'] == ['boundary']
    assert gate['n_rejected_membership'] == 1


def test_a_finite_redshift_without_a_quality_column_is_not_the_lens_count():
    images = [_image('global', 'global', 0.0, 10.0)]
    clusters = [_cluster('A', 5.0, 5.0, 0.25)]
    rows = [_shape('s1', 5.0, 5.0, 1.0), _shape('s2', 5.2, 5.0, 0.3125)]
    shapes = _shapes(rows, ra_max=10.0, quality=False, crossmatch=False)
    gate = named_data_gate(images, clusters, shapes, dz=0.125, z_gap=0.25)
    assert gate['status'] == 'untrusted_redshift'
    assert gate['n_lens_sources'] is None
    assert gate['n_untrusted_finite_z'] == 1
    assert gate['science_measurement'] is False
    assert gate['n_rejected_quality'] is None
    assert 'quality flag' in gate['note']


def test_a_missing_crossmatch_column_does_not_drop_the_row():
    images = [_image('global', 'global', 0.0, 10.0)]
    clusters = [_cluster('A', 5.0, 5.0, 0.25)]
    rows = [_shape('s1', 5.0, 5.0, 1.0, crossmatch=False)]
    open_match = named_data_gate(
        images, clusters, _shapes(rows, ra_max=10.0, crossmatch=False), dz=0.125, z_gap=0.25)
    closed_match = named_data_gate(
        images, clusters, _shapes(rows, ra_max=10.0, crossmatch=True), dz=0.125, z_gap=0.25)
    assert open_match['n_lens_sources'] == 1
    assert open_match['n_rejected_crossmatch'] is None
    assert closed_match['n_lens_sources'] == 0
    assert closed_match['n_rejected_crossmatch'] == 1


def test_wrong_shape_version_is_not_counted_as_zero():
    images, clusters, shapes = _catalogues(version='S16A')
    gate = named_data_gate(images, clusters, shapes, dz=0.125)
    assert gate['status'] == 'version_mismatch'
    assert gate['n_lens_sources'] is None
    assert gate['combination_matched'] is False
    assert gate['science_measurement'] is False


def test_wrapped_right_ascension_joins_the_footprint():
    images = [_image('global', 'global', 350.0, 10.0)]
    images[0]['footprint']['dec_min'] = -1.0
    images[0]['footprint']['dec_max'] = 1.0
    clusters = [_cluster('A', 359.0, 0.0, 0.25), _cluster('C', 20.0, 0.0, 0.25)]
    box = {'ra_min': 350.0, 'ra_max': 10.0, 'dec_min': -1.0, 'dec_max': 1.0}
    shapes = {
        'version': 'S19A',
        'footprint': box,
        'spec_quality_supplied': True,
        'crossmatch_supplied': True,
        'rows': [
            _shape('near', 1.0, 0.0, 1.0),
            _shape('far', 20.0, 0.0, 1.0),
        ],
    }
    gate = named_data_gate(images, clusters, shapes, dz=0.125, z_gap=0.25)
    assert gate['n_lens_sources'] == 1
    assert gate['per_cluster'][0]['source_ids'] == ['near']
    assert gate['clusters_excluded'][0]['id'] == 'C'


def test_filled_failure_notes_still_leave_the_science_flag_down():
    notes = {key: 'checked on the supplied rows' for key in FAILURE_MODES}
    images, clusters, shapes = _catalogues()
    gate = named_data_gate(images, clusters, shapes, dz=0.125, z_gap=0.25, failure_notes=notes)
    assert gate['checklist_filled'] is True
    assert gate['n_lens_sources'] == 3
    assert gate['science_measurement'] is False
    assert gate['named_data_mass'] is False
    assert gate['survey_calibration'] is False
    assert gate['observed_sentence'] is None
    assert gate['observed_sentence_allowed'] is False
    with pytest.raises(ValueError):
        named_data_gate(failure_notes={'sky_oversubtraction': True})
    with pytest.raises(ValueError):
        named_data_gate(failure_notes={'not_a_mode': 'text'})


def test_reinforcement_bins_count_only_supplied_edges():
    bare = named_data_gate(reinforcement={'survey': 'DESI DR1', 'rows': [
        {'brightness': 19.0, 'colour': 1.0, 'surface_brightness': 26.0, 'radius': 0.4,
         'in_intersection': True},
    ]})
    assert bare['reinforcement']['n_rows'] == 1
    assert bare['reinforcement']['bins']['brightness'] is None
    assert bare['sky_measurement'] is False
    rows = [
        {'brightness': 19.0, 'colour': 0.5, 'surface_brightness': 25.0, 'radius': 0.2,
         'in_intersection': True},
        {'brightness': 20.0, 'colour': 1.5, 'surface_brightness': 27.0, 'radius': 0.8,
         'in_intersection': True},
        {'brightness': 19.5, 'colour': 0.5, 'surface_brightness': 25.0, 'radius': 0.2,
         'in_intersection': False},
    ]
    edges = {
        'brightness': [18.0, 20.0, 22.0],
        'colour': [0.0, 1.0, 2.0],
        'surface_brightness': [24.0, 26.0, 28.0],
        'radius': [0.0, 0.5, 1.0],
    }
    gate = named_data_gate(reinforcement={'survey': 'SDSS', 'rows': rows, 'edges': edges})
    table = gate['reinforcement']
    assert table['n_rows'] == 3
    assert table['n_in_intersection'] == 2
    assert table['n_outside'] == 1
    assert table['sky_measurement'] is False
    assert [item['n'] for item in table['bins']['brightness']] == [1, 1]
    assert table['bins']['brightness'][1]['low'] == 20.0
    assert 'not a survey measurement' in table['note']
    with pytest.raises(ValueError):
        named_data_gate(reinforcement={'survey': '2dF', 'rows': []})


def test_partial_catalogues_and_bad_windows_raise():
    images, clusters, shapes = _catalogues()
    with pytest.raises(ValueError):
        named_data_gate(images, clusters, dz=0.1)
    with pytest.raises(ValueError):
        named_data_gate(images, clusters, shapes)
    with pytest.raises(ValueError):
        named_data_gate(images, clusters, shapes, dz=0.0)
    with pytest.raises(ValueError):
        named_data_gate(images, clusters, shapes, dz=0.1, z_gap=-0.1)
    with pytest.raises(ValueError):
        named_data_gate(dz=0.1)
    repeated = _shapes([_shape('s1', 1.0, 1.0, 1.0), _shape('s1', 2.0, 2.0, 1.0)], ra_max=10.0)
    with pytest.raises(ValueError):
        named_data_gate([_image('global', 'global', 0.0, 10.0)], [_cluster('A', 5.0, 5.0, 0.25)],
                        repeated, dz=0.1)
