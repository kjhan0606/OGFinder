"""Named-data gate in front of an observed science sentence.

The first combination is HSC-SSP PDR3 images, a CAMIRA cluster list, and the
S19A public shape catalogue. `named_data_gate` matches those versions and the
caller-supplied footprints, then counts lens sources that survive the current
spectroscopic limit. An ICL image is a `deepCoadd` frame whose sky is global.
A local-sky frame stays out of the common footprint even when a cluster and a
shape row lie inside it.

The spectroscopic limit is the one already used for weak lensing. A finite
redshift is not a quality flag, so the count stays empty unless the caller
marks the quality column as supplied. A bad cross-match drops out when that
column is supplied. `classify_spectroscopic` labels `|z - z_cl| < dz` as a
satellite and a redshift outside that window as background, foreground included.
`background_for_lensing` then keeps a background row whose redshift is finite
and behind the lens by `z_gap`. The boundary `|z - z_cl| = dz` is outside the
cluster window. A photometric redshift is not copied into `z`. The radial window
and the shear cuts inside `select_sources` are a later cut.

This function does not download or reduce a survey. When the versions and the
common footprint match, `n_lens_sources` counts the rows the caller passed.
`science_measurement` stays false. `observed_sentence` stays empty. A filled
failure-mode checklist does not make the count a named-data mass or a survey
calibration. S19A shapes are a re-Gaussianization catalogue, not the
metacalibration in `measure_stamps`.
"""

import numpy as np

from icl.members import background_for_lensing, classify_spectroscopic


IMAGES = 'HSC-SSP PDR3'
CLUSTERS = 'CAMIRA'
SHAPES = 'S19A'
ICL_PRODUCT = 'deepCoadd'
GLOBAL_SKY = ('global', 'global-sky')
LOCAL_SKY = ('local', 'local-sky')
REINFORCEMENT_SURVEYS = ('SDSS', 'DESI DR1')
INTERSECTION_AXES = ('brightness', 'colour', 'surface_brightness', 'radius')

FAILURE_MODES = (
    'sky_oversubtraction',
    'psf_wings_after_u',
    'galactic_cirrus',
    'blending',
    'member_dilution',
    'photoz_nz',
    'centring',
    'selection_depth',
)

_FAILURE_STATEMENTS = {
    'sky_oversubtraction': (
        'Sky over-subtraction is a common systematic. Adding clusters does not remove it.'),
    'psf_wings_after_u': (
        'BCG and star PSF wings take different angular sizes after the u transform.'),
    'galactic_cirrus': 'Galactic cirrus can couple to extinction and selection.',
    'blending': (
        'Blending changes the satellite mask, the light denominator, shapes, and colours together.'),
    'member_dilution': (
        'Member dilution weakens tangential shear. '
        'A source-density boost does not by itself separate obscuration from magnification.'),
    'photoz_nz': (
        'Photo-z and n(z) errors tilt lensing efficiency and the mass normalisation. '
        'Cluster-neighbour contamination can differ from an empty-field correction.'),
    'centring': 'Centring error blurs the inner lens and the ICL profile.',
    'selection_depth': (
        'Selection and depth can average different populations at different radius, mass, and redshift '
        'and can create a false trend.'),
}

_OTHER_DATA = (
    {
        'data': 'HSC-SSP PDR3 images, CAMIRA, S19A shapes and photo-z PDFs',
        'role': 'first combination for ICL and lensing on one footprint',
        'constraint': (
            'ICL images are deepCoadd global-sky frames. '
            'S19A shapes are re-Gaussianization, a different calibration from measure_stamps.'),
    },
    {
        'data': 'SDSS corrected images and DR8 redMaPPer',
        'role': 'low-redshift ICL comparison',
        'constraint': (
            'Sky treatment and the outer light lost around bright galaxies stay in the comparison. '
            'Ordinary photometric shapes are not a shear catalogue.'),
    },
    {
        'data': 'DES DR2 images and Y3 metacalibration products',
        'role': 'southern ICL and a separately calibrated shear',
        'constraint': (
            'The imaging and the Y3 lensing products are different products. '
            'Selection, mask, response, and redshift calibration have to be matched again.'),
    },
    {
        'data': 'KiDS DR4/4.1 and KiDS-1000 shear',
        'role': 'lensing alternative and images of overlapping clusters',
        'constraint': (
            'ICL preservation needs its own check. '
            'lensfit is not the metacalibration in measure_stamps.'),
    },
)

_FRONTIER_FIELDS = {
    'use': 'centre mask and blending example',
    'sample': 'six clusters, including Abell S1063, from MAST',
    'stands_in_for_wide_radius_icl_mean': False,
    'stands_in_for_representative_nfw': False,
}

_DESIGN_SCALE = {
    'detection': False,
    'note': (
        'These are design ranges, not detections. An ICL mean profile starts from 10^2 to 10^3 clusters, '
        'and from hundreds in each mass and redshift bin. Zibetti et al. 2005 and Zhang et al. 2019 '
        'discussed outer profiles around 30-32 mag arcsec^-2 with samples of hundreds. That level is '
        'not an individual-image point-source limit. An intermediate-redshift sample starts from hundreds '
        'to thousands of square degrees. A cutout keeps the target u range and an outer region that '
        'constrains the sky. Survey area and the effective area of each cluster are reported separately. '
        'The lensing design range is roughly 24-26 AB in several bands, with good seeing, and a few to '
        'tens of shape sources per square arcminute after quality cuts. That density is the scale before '
        'the spectroscopic cut. The figure of merit is sigma(g_bar) ~ sigma_e/sqrt(N_eff) and the '
        'lensing efficiency of each source.'),
}

_SHAPE_CALIBRATION = {
    'catalogue': SHAPES,
    'method': 're-Gaussianization',
    'same_as_measure_stamps': False,
    'note': (
        'S19A shapes are a re-Gaussianization catalogue. They are not the metacalibration in '
        'measure_stamps, a KiDS lensfit calibration, or a DES Y3 metacalibration.'),
}

_COUNT_NOTE = (
    'n_lens_sources counts supplied rows in the matched footprint that pass the spectroscopic lens cut. '
    'The count is not a sky measurement.')
_CLOSED_NOTE = (
    'science_measurement stays false. This function does not read a survey release, and a filled '
    'failure-mode checklist is not a named-data mass or a survey calibration.')


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _finite(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError('%s must be finite' % name)
    if not np.isfinite(number):
        raise ValueError('%s must be finite' % name)
    return number


def _box(value, name):
    if not isinstance(value, dict):
        raise ValueError('%s footprint must give ra_min, ra_max, dec_min, and dec_max in degrees' % name)
    edges = []
    for key in ('ra_min', 'ra_max', 'dec_min', 'dec_max'):
        if key not in value:
            raise ValueError('%s footprint needs %s' % (name, key))
        edges.append(_finite(value[key], '%s %s' % (name, key)))
    ra_min, ra_max, dec_min, dec_max = edges
    if dec_min < -90.0 or dec_max > 90.0 or dec_min > dec_max:
        raise ValueError('%s declination must lie in [-90, 90] with dec_min <= dec_max' % name)
    if ra_min == ra_max or ra_min < 0.0 or ra_max > 360.0 or ra_min > 360.0 or ra_max < 0.0:
        raise ValueError('%s right ascension must lie in [0, 360] and the box must have a width' % name)
    return (ra_min, ra_max, dec_min, dec_max)


def _boxes(record, name):
    if 'footprints' in record and record['footprints'] is not None:
        raw = record['footprints']
        if not isinstance(raw, (list, tuple)) or len(raw) == 0:
            raise ValueError('%s needs at least one footprint' % name)
        return [_box(item, name) for item in raw]
    if 'footprint' not in record or record['footprint'] is None:
        raise ValueError('%s needs a footprint' % name)
    return [_box(record['footprint'], name)]


def _inside(ra, dec, box):
    ra_min, ra_max, dec_min, dec_max = box
    if dec < dec_min or dec > dec_max:
        return False
    if ra_min < ra_max:
        return ra_min <= ra <= ra_max
    return ra >= ra_min or ra <= ra_max


def _inside_any(ra, dec, boxes):
    return any(_inside(ra, dec, box) for box in boxes)


def _identifier(record, index):
    if record.get('id') is None:
        return str(index)
    return str(record['id'])


def _failure_report(notes):
    notes = {} if notes is None else dict(notes)
    unknown = [key for key in notes if key not in FAILURE_MODES]
    if unknown:
        raise ValueError('unknown failure-mode key: %s' % ', '.join(sorted(unknown)))
    report = {}
    for key in FAILURE_MODES:
        caller = notes.get(key)
        if caller is None:
            filled = False
            stored = None
        elif not _text(caller):
            raise ValueError('failure-mode %s needs a non-empty note of what was checked' % key)
        else:
            filled = True
            stored = str(caller).strip()
        report[key] = {
            'statement': _FAILURE_STATEMENTS[key],
            'caller': stored,
            'filled': filled,
        }
    return report, all(item['filled'] for item in report.values())


def _axis_bins(values, edges, axis):
    edges = [_finite(value, '%s edge' % axis) for value in edges]
    if len(edges) < 2:
        raise ValueError('%s needs a low edge and a high edge' % axis)
    if any(high <= low for low, high in zip(edges, edges[1:])):
        raise ValueError('%s edges must increase' % axis)
    counts = [0 for _ in range(len(edges) - 1)]
    unbinned = 0
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = np.nan
        if not np.isfinite(number) or number < edges[0] or number > edges[-1]:
            unbinned += 1
            continue
        for index in range(len(edges) - 1):
            last = index == len(edges) - 2
            below_high = number <= edges[index + 1] if last else number < edges[index + 1]
            if number >= edges[index] and below_high:
                counts[index] += 1
                break
    bins = [
        {'low': edges[index], 'high': edges[index + 1], 'n': counts[index]}
        for index in range(len(counts))]
    return bins, unbinned


def _reinforcement(spec):
    base = {
        'surveys': list(REINFORCEMENT_SURVEYS),
        'survey': None,
        'bins': {axis: None for axis in INTERSECTION_AXES},
        'n_rows': None,
        'n_in_intersection': None,
        'n_outside': None,
        'n_unbinned': None,
        'sky_measurement': False,
        'note': (
            'SDSS and DESI DR1 are spectroscopic reinforcement candidates. '
            'Target selection does not represent every background galaxy in a deep image. '
            'The brightness, colour, surface-brightness, and radius dependence stays empty until '
            'the caller supplies rows and bin edges. Bin counts then count those rows and are not '
            'a survey measurement.'),
    }
    if spec is None:
        return base
    if not isinstance(spec, dict):
        raise ValueError('reinforcement must describe one of SDSS or DESI DR1')
    survey = spec.get('survey')
    if survey not in REINFORCEMENT_SURVEYS:
        raise ValueError('reinforcement survey must be SDSS or DESI DR1')
    rows = spec.get('rows') or []
    if not isinstance(rows, (list, tuple)):
        raise ValueError('reinforcement rows must be a list')
    base['survey'] = survey
    base['n_rows'] = len(rows)
    inside = []
    outside = 0
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError('each reinforcement row must be a record')
        if row.get('in_intersection') is True:
            inside.append(row)
        else:
            outside += 1
    base['n_in_intersection'] = len(inside)
    base['n_outside'] = outside
    edges = spec.get('edges')
    if edges is None:
        base['note'] += ' Bin edges were not supplied, so no dependence was counted.'
        return base
    if not isinstance(edges, dict):
        raise ValueError('reinforcement edges must name brightness, colour, surface_brightness, and radius')
    missing = [axis for axis in INTERSECTION_AXES if axis not in edges or edges[axis] is None]
    if missing:
        raise ValueError('reinforcement edges need %s' % ', '.join(missing))
    unbinned = {}
    bins = {}
    for axis in INTERSECTION_AXES:
        values = []
        for row in inside:
            if axis not in row:
                raise ValueError('reinforcement row needs %s' % axis)
            values.append(row[axis])
        axis_bins, axis_unbinned = _axis_bins(values, edges[axis], axis)
        bins[axis] = axis_bins
        unbinned[axis] = axis_unbinned
    base['bins'] = bins
    base['n_unbinned'] = unbinned
    return base


def _empty_result(failure, checklist_filled, reinforcement):
    return {
        'combination': '%s + %s + %s' % (IMAGES, CLUSTERS, SHAPES),
        'status': 'empty',
        'combination_matched': False,
        'sky_measurement': False,
        'science_measurement': False,
        'named_data_mass': False,
        'survey_calibration': False,
        'observed_sentence': None,
        'observed_sentence_allowed': False,
        'n_lens_sources': None,
        'n_untrusted_finite_z': None,
        'n_accepted_images': None,
        'n_accepted_clusters': None,
        'n_shapes': None,
        'n_shapes_in_common_footprint': None,
        'n_shapes_outside': None,
        'n_rejected_quality': None,
        'n_rejected_crossmatch': None,
        'n_rejected_no_spec_z': None,
        'n_rejected_membership': None,
        'per_cluster': [],
        'clusters_excluded': [],
        'images_accepted': [],
        'images_rejected': [],
        'spec_quality_supplied': False,
        'crossmatch_supplied': False,
        'dz': None,
        'z_gap': None,
        'checklist_filled': checklist_filled,
        'failure_modes': failure,
        'reinforcement': reinforcement,
        'other_data': [dict(item) for item in _OTHER_DATA],
        'frontier_fields': dict(_FRONTIER_FIELDS),
        'design_scale': dict(_DESIGN_SCALE),
        'shape_calibration': dict(_SHAPE_CALIBRATION),
        'note': (
            'No catalogue was supplied. ' + _CLOSED_NOTE),
    }


def _sky_kind(value):
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in GLOBAL_SKY:
        return 'global'
    if text in LOCAL_SKY:
        return 'local'
    return None


def _accept_image(record, index):
    identity = _identifier(record, index)
    version = record.get('version')
    product = record.get('product')
    sky = _sky_kind(record.get('sky'))
    if version != IMAGES:
        return None, {'id': identity, 'reason': 'version is not HSC-SSP PDR3'}
    if product != ICL_PRODUCT:
        return None, {'id': identity, 'reason': 'product is not deepCoadd'}
    if sky == 'local':
        return None, {
            'id': identity,
            'reason': 'local-sky is not an ICL input, including when the catalogue matches the frame',
        }
    if sky != 'global':
        return None, {'id': identity, 'reason': 'sky was not named global'}
    try:
        boxes = _boxes(record, 'image %s' % identity)
    except ValueError as exc:
        raise ValueError(str(exc))
    return {'id': identity, 'boxes': boxes}, None


def _quality_passes(row, minimum):
    if 'spec_quality' not in row or row['spec_quality'] is None:
        return False
    quality = row['spec_quality']
    if isinstance(quality, (bool, np.bool_)):
        return bool(quality)
    try:
        number = float(quality)
    except (TypeError, ValueError):
        raise ValueError('spec_quality must be boolean or numeric')
    return bool(np.isfinite(number) and number >= minimum)


def _spec_redshift(row):
    """Return the spectroscopic redshift, or None when the row does not have one.

    A photometric value, including one copied by the caller into a flag that says
    the column is not spectroscopic, does not become the lensing redshift.
    """
    if 'z_is_spec' in row and row['z_is_spec'] is False:
        return None
    if 'z' not in row or row['z'] is None:
        return None
    try:
        redshift = float(row['z'])
    except (TypeError, ValueError):
        return None
    if not np.isfinite(redshift):
        return None
    return redshift


def _behind_lens(redshifts, stars, z_cluster, dz, z_gap):
    if redshifts.size == 0:
        return np.zeros(0, bool)
    labels = classify_spectroscopic(redshifts, z_cluster, dz, star=stars)
    return np.asarray(background_for_lensing(labels, redshifts, z_cluster, z_gap), bool)


def _row_gate(row, quality_supplied, crossmatch_supplied, minimum):
    """Classify one shape row before the per-cluster membership cut.

    Returns 'quality', 'crossmatch', 'no_spec_z', or a finite redshift.
    """
    if quality_supplied and not _quality_passes(row, minimum):
        return 'quality'
    if crossmatch_supplied and row.get('crossmatch_ok') is not True:
        return 'crossmatch'
    redshift = _spec_redshift(row)
    if redshift is None:
        return 'no_spec_z'
    return redshift


def named_data_gate(images=None, clusters=None, shapes=None, dz=None, z_gap=0.1,
                    min_spec_quality=1.0, failure_notes=None, reinforcement=None):
    """Match the first named-data combination and count surviving lens sources.

    Pass all three catalogues together, or pass none. `images` is a list of
    records with version, product, sky, and a footprint in degrees. `clusters`
    is a list of records with version, ra, dec, and z. `shapes` is one record
    with version, a footprint, and rows. Footprint edges are inclusive. A box
    with ra_min greater than ra_max wraps through 0 degrees.

    `spec_quality_supplied` on the shape catalogue defaults to false. While it
    is false, a finite redshift is counted only as `n_untrusted_finite_z` and
    `n_lens_sources` stays empty. `z_gap` defaults to 0.1, the same default as
    `select_sources`. `dz` is the caller's spectroscopic window and has no
    default. Failure-mode notes are text that says what was checked. A boolean
    does not fill an item.
    """
    failure, checklist_filled = _failure_report(failure_notes)
    reinforcement_report = _reinforcement(reinforcement)
    supplied = [images is not None, clusters is not None, shapes is not None]
    if not any(supplied):
        if dz is not None:
            raise ValueError('dz is unused when no catalogue is supplied')
        return _empty_result(failure, checklist_filled, reinforcement_report)
    if not all(supplied):
        raise ValueError('the first combination needs images, clusters, and shapes together')
    if not isinstance(images, (list, tuple)) or not isinstance(clusters, (list, tuple)):
        raise ValueError('images and clusters must be lists of records')
    if not isinstance(shapes, dict) or 'rows' not in shapes:
        raise ValueError('shapes must name a version, a footprint, and rows')
    dz = _finite(dz, 'dz')
    if dz <= 0:
        raise ValueError('dz must be positive')
    z_gap = _finite(z_gap, 'z_gap')
    if z_gap < 0:
        raise ValueError('z_gap must be non-negative')
    minimum = _finite(min_spec_quality, 'min_spec_quality')

    accepted_images = []
    rejected_images = []
    for index, record in enumerate(images):
        if not isinstance(record, dict):
            raise ValueError('each image must be a record')
        accepted, rejected = _accept_image(record, index)
        if accepted is None:
            rejected_images.append(rejected)
        else:
            accepted_images.append(accepted)

    shape_version = shapes.get('version')
    shape_boxes = _boxes(shapes, 'shapes') if shape_version == SHAPES else []
    rows = shapes.get('rows')
    if not isinstance(rows, (list, tuple)):
        raise ValueError('shape rows must be a list')
    quality_supplied = bool(shapes.get('spec_quality_supplied', False))
    crossmatch_supplied = bool(shapes.get('crossmatch_supplied', False))

    seen_ids = set()
    shape_rows = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError('each shape row must be a record')
        identity = _identifier(row, index)
        if identity in seen_ids:
            raise ValueError('shape id is repeated')
        seen_ids.add(identity)
        if 'ra' not in row or 'dec' not in row:
            raise ValueError('shape row %s needs ra and dec' % identity)
        ra = _finite(row['ra'], 'ra')
        dec = _finite(row['dec'], 'dec')
        if ra < 0.0 or ra > 360.0 or dec < -90.0 or dec > 90.0:
            raise ValueError('shape row %s coordinates are outside the sphere' % identity)
        shape_rows.append({
            'id': identity,
            'ra': ra,
            'dec': dec,
            'row': row,
            'star': bool(row.get('star', False)),
        })

    local_boxes = []
    for raw in images:
        if not isinstance(raw, dict):
            continue
        if _sky_kind(raw.get('sky')) != 'local' or raw.get('version') != IMAGES or raw.get('product') != ICL_PRODUCT:
            continue
        try:
            local_boxes.extend(_boxes(raw, 'image'))
        except ValueError:
            continue

    accepted_clusters = []
    excluded_clusters = []
    seen_clusters = set()
    for index, record in enumerate(clusters):
        if not isinstance(record, dict):
            raise ValueError('each cluster must be a record')
        identity = _identifier(record, index)
        if identity in seen_clusters:
            raise ValueError('cluster id is repeated')
        seen_clusters.add(identity)
        if record.get('version') != CLUSTERS:
            excluded_clusters.append({'id': identity, 'reason': 'version is not CAMIRA'})
            continue
        if 'ra' not in record or 'dec' not in record or 'z' not in record:
            raise ValueError('cluster %s needs ra, dec, and z' % identity)
        ra = _finite(record['ra'], 'cluster ra')
        dec = _finite(record['dec'], 'cluster dec')
        redshift = _finite(record['z'], 'cluster z')
        if ra < 0.0 or ra > 360.0 or dec < -90.0 or dec > 90.0:
            raise ValueError('cluster %s coordinates are outside the sphere' % identity)
        if redshift < 0.0:
            raise ValueError('cluster redshift must be non-negative')
        covering = [image for image in accepted_images if _inside_any(ra, dec, image['boxes'])]
        if shape_version != SHAPES:
            excluded_clusters.append({'id': identity, 'reason': 'shape catalogue version is not S19A'})
            continue
        if not covering and _inside_any(ra, dec, local_boxes):
            excluded_clusters.append({
                'id': identity,
                'reason': 'The only covering image is local-sky, which is not an ICL input.',
            })
            continue
        if not covering:
            excluded_clusters.append({
                'id': identity,
                'reason': 'outside the HSC-SSP PDR3 deepCoadd global-sky footprint',
            })
            continue
        if not _inside_any(ra, dec, shape_boxes):
            excluded_clusters.append({'id': identity, 'reason': 'outside the S19A footprint'})
            continue
        accepted_clusters.append({
            'id': identity,
            'ra': ra,
            'dec': dec,
            'z': redshift,
            'images': covering,
        })

    result = _empty_result(failure, checklist_filled, reinforcement_report)
    result['dz'] = dz
    result['z_gap'] = z_gap
    result['spec_quality_supplied'] = quality_supplied
    result['crossmatch_supplied'] = crossmatch_supplied
    result['n_shapes'] = len(shape_rows)
    result['n_accepted_images'] = len(accepted_images)
    result['n_accepted_clusters'] = len(accepted_clusters)
    result['images_accepted'] = [image['id'] for image in accepted_images]
    result['images_rejected'] = rejected_images
    result['clusters_excluded'] = excluded_clusters
    result['note'] = _CLOSED_NOTE

    if shape_version != SHAPES:
        result['status'] = 'version_mismatch'
        result['note'] = 'The shape catalogue version is not S19A. ' + _CLOSED_NOTE
        return result
    if not accepted_clusters:
        result['status'] = 'no_common_footprint'
        result['note'] = (
            'The common footprint of HSC-SSP PDR3 deepCoadd global-sky, CAMIRA, and S19A is empty. '
            + _CLOSED_NOTE)
        return result

    result['combination_matched'] = True
    per_cluster_ids = {cluster['id']: [] for cluster in accepted_clusters}
    survivors = []
    outside = 0
    rejected = {'quality': 0, 'crossmatch': 0, 'no_spec_z': 0, 'membership': 0}
    untrusted = []
    in_common = 0
    for source in shape_rows:
        if not _inside_any(source['ra'], source['dec'], shape_boxes):
            outside += 1
            continue
        clusters_here = []
        for cluster in accepted_clusters:
            covered = any(_inside_any(source['ra'], source['dec'], image['boxes']) for image in cluster['images'])
            if covered:
                clusters_here.append(cluster)
        if not clusters_here:
            outside += 1
            continue
        in_common += 1
        decision = _row_gate(source['row'], quality_supplied, crossmatch_supplied, minimum)
        if decision in rejected:
            rejected[decision] += 1
            continue
        matched = []
        for cluster in clusters_here:
            passed = _behind_lens(
                np.array([decision], float), np.array([source['star']], bool),
                cluster['z'], dz, z_gap)
            if bool(passed[0]):
                matched.append(cluster['id'])
        if not matched:
            rejected['membership'] += 1
            continue
        for identity in matched:
            per_cluster_ids[identity].append(source['id'])
        if quality_supplied:
            survivors.append(source['id'])
        else:
            untrusted.append(source['id'])

    result['n_shapes_outside'] = outside
    result['n_shapes_in_common_footprint'] = in_common
    result['n_rejected_quality'] = rejected['quality'] if quality_supplied else None
    result['n_rejected_crossmatch'] = rejected['crossmatch'] if crossmatch_supplied else None
    result['n_rejected_no_spec_z'] = rejected['no_spec_z']
    result['n_rejected_membership'] = rejected['membership']
    result['per_cluster'] = [
        {
            'id': cluster['id'],
            'z': cluster['z'],
            'n_lens_sources': len(per_cluster_ids[cluster['id']]) if quality_supplied else None,
            'n_untrusted_finite_z': None if quality_supplied else len(per_cluster_ids[cluster['id']]),
            'source_ids': list(per_cluster_ids[cluster['id']]),
        }
        for cluster in accepted_clusters]
    if quality_supplied:
        result['status'] = 'counted'
        result['n_lens_sources'] = len(survivors)
        result['n_untrusted_finite_z'] = None
        result['note'] = _COUNT_NOTE + ' ' + _CLOSED_NOTE
    else:
        result['status'] = 'untrusted_redshift'
        result['n_lens_sources'] = None
        result['n_untrusted_finite_z'] = len(untrusted)
        result['note'] = (
            'A finite redshift was not treated as a spectroscopic quality flag. '
            'n_untrusted_finite_z counts supplied rows and is not a sky measurement. '
            + _CLOSED_NOTE)
    return result
