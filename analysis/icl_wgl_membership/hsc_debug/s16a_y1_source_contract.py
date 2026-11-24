"""Offline S16A/Y1 Table-4 input checks. No lensing estimator or network.

This is a diagnostic reconstruction of the published object cuts, not a
complete weak-lensing selection (area/mask and redshift calibration are separate).
The pilot's ``sigma_e`` is calibrated; ``raw_sigma`` below must never alias it.
"""
from __future__ import annotations

import math
import re

PAPER = 'https://arxiv.org/html/1705.06745v3#A2'
HSC = 'https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/'
COSMOS = 'https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/'
BANDS = ('g', 'r', 'z', 'y')
BAD_FLAGS = (
    'deblend_skipped', 'bad_centroid', 'centroid_sdss_failed', 'pixel_edge',
    'pixel_interpolated_center', 'pixel_saturated_center', 'pixel_cr_center',
    'pixel_bad', 'pixel_suspect_center', 'pixel_clipped_any', 'shape_failed',
)
# Explicit public PDR1 database identifiers, rather than paper shorthand.
EXPRESSIONS = {
    'object_id': 'r.object_id',
    'is_primary': 'm.idetect_is_primary',
    'deblend_skipped': 'm.ideblend_skipped',
    'bad_centroid': 'm.iflags_badcentroid',
    'centroid_sdss_failed': 'a.icentroid_sdss_flags',
    'pixel_edge': 'm.iflags_pixel_edge',
    'pixel_interpolated_center': 'm.iflags_pixel_interpolated_center',
    'pixel_saturated_center': 'm.iflags_pixel_saturated_center',
    'pixel_cr_center': 'm.iflags_pixel_cr_center',
    'pixel_bad': 'm.iflags_pixel_bad',
    'pixel_suspect_center': 'm.iflags_pixel_suspect_center',
    'pixel_clipped_any': 'm.iflags_pixel_clipped_any',
    'shape_failed': 'a.ishape_hsm_regauss_flags',
    'raw_sigma': 'a.ishape_hsm_regauss_sigma',
    'resolution': 'a.ishape_hsm_regauss_resolution',
    'extendedness': 'm.iclassification_extendedness',
    'blendedness': 'm.iblendedness_abs_flux',
    'i_flux': 'm.icmodel_flux', 'i_flux_err': 'm.icmodel_flux_err',
    'i_mag': 'm.icmodel_mag', 'a_i': 'f.a_i',
    **{key: f'm.{band}cmodel_{suffix}'
       for band in BANDS for key, suffix in
       ((f'{band}_flux', 'flux'), (f'{band}_flux_err', 'flux_err'))},
    'photoz_best': 'z.photoz_best', 'photoz_risk_best': 'z.photoz_risk_best',
}
TABLES = {'m': 'meas', 'a': 'meas2', 'f': 'forced', 'z': 'photoz_ephor_ab'}
NUMERIC_INPUTS = ('raw_sigma', 'resolution', 'extendedness', 'blendedness',
                  'i_flux', 'i_flux_err', 'i_mag', 'a_i',
                  *(f'{b}_{suffix}' for b in BANDS for suffix in ('flux', 'flux_err')))
REQUIRED = ('e1', 'e2', 'is_primary', *BAD_FLAGS, *NUMERIC_INPUTS)
BLEND_MAX = 10 ** (-0.375)


def object_id(value):
    """Canonical positive signed-bigint ID; floating-point IDs are forbidden."""
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ValueError('object_id must be an exact integer string or integer')
    text = str(value)
    if not re.fullmatch(r'[1-9][0-9]*', text) or int(text) > 2**63 - 1:
        raise ValueError('invalid signed-bigint object_id')
    return text


def finite(value):
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (ValueError, TypeError):
        return None
    return number if math.isfinite(number) else None


def boolean(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return {'true': True, 't': True, 'false': False, 'f': False}.get(value.strip().lower())
    return None  # No Python truthiness: in particular bool('false') is True.


def state(value):
    return 'unknown' if value is None else ('pass' if value else 'fail')


def select_row(row):
    """Three-valued, diagnostic object cuts; unknown input never passes.

    For the two-of-four detection rule, two known passing bands suffice even
    if another is missing. A remaining ambiguity is reported, not filled with 0.
    No m, derived-sigma, point-z or Y3 selection is implicitly imposed.
    """
    c = {'primary': boolean(row.get('is_primary'))}
    c.update({key: None if (value := boolean(row.get(key))) is None else not value
              for key in BAD_FLAGS})
    n = {key: finite(row.get(key)) for key in (*NUMERIC_INPUTS, 'e1', 'e2')}

    def check(keys, predicate):
        return None if any(n[k] is None for k in keys) else bool(predicate())

    def snr(band, threshold):
        keys = (f'{band}_flux', f'{band}_flux_err')
        if any(n[k] is None for k in keys) or n[keys[1]] <= 0:
            return None
        return n[keys[0]] / n[keys[1]] >= threshold

    c['i_snr_ge_10'] = snr('i', 10)
    c['resolution_ge_0p3'] = check(['resolution'], lambda: n['resolution'] >= .3)
    c['distortion_lt_2'] = check(['e1', 'e2'], lambda: math.hypot(n['e1'], n['e2']) < 2)
    c['raw_sigma_0_to_0p4'] = check(['raw_sigma'], lambda: 0 <= n['raw_sigma'] <= .4)
    c['extended'] = check(['extendedness'], lambda: n['extendedness'] != 0)
    c['dereddened_i_le_24p5'] = check(['i_mag', 'a_i'], lambda: n['i_mag'] - n['a_i'] <= 24.5)
    c['blend_lt_limit'] = check(['blendedness'], lambda: n['blendedness'] < BLEND_MAX)
    bands = {b: snr(b, 5) for b in BANDS}
    passed = sum(v is True for v in bands.values())
    unknown = sum(v is None for v in bands.values())
    c['two_other_bands_ge_5'] = True if passed >= 2 else (False if passed + unknown < 2 else None)
    decision = False if any(v is False for v in c.values()) else (None if any(v is None for v in c.values()) else True)
    invalid = [key for key in REQUIRED if
               (boolean(row.get(key)) is None if key in ('is_primary', *BAD_FLAGS)
                else finite(row.get(key)) is None)]
    return dict(table4_diagnostic=state(decision), criteria={k: state(v) for k, v in c.items()},
                other_band_snr={k: state(v) for k, v in bands.items()},
                missing_or_invalid_inputs=invalid, science_measurement=False)


def merge_supplement(pilot_rows, supplement_rows):
    """Exact object join, preserving every center pair and original property.

    Entirely missing rows, extras, duplicates or overwritten properties abort.
    NULL cells remain NULL for diagnostic selection after the complete join.
    """
    pilot_rows = list(pilot_rows)
    expected = {object_id(row['object_id']) for row in pilot_rows}
    supplied = {}
    for row in supplement_rows:
        key = object_id(row['object_id'])
        if key not in expected or key in supplied:
            raise ValueError('unknown or duplicate supplement object_id')
        supplied[key] = row
    if set(supplied) != expected:
        raise ValueError('missing supplement object_id')
    pairs, output = set(), []
    for row in pilot_rows:
        key = object_id(row['object_id'])
        pair = (row['center_id'], key)
        if pair in pairs:
            raise ValueError('duplicate pilot center/object pair')
        pairs.add(pair)
        add = supplied[key]
        if (set(row) & set(add)) - {'object_id'}:
            raise ValueError('supplement would overwrite pilot properties')
        output.append({**row, **add})
    return output


def contract():
    return dict(
        schema='astrafek-s16a-y1-input-contract-v1', purpose='offline_table4_diagnostics_only',
        references=dict(table4=PAPER, release=HSC, cosmos=COSMOS),
        catalog='pdr1.s16a_wide', additional_expressions=EXPRESSIONS,
        existing_fields=['center_id', 'object_id', 'source_ra', 'source_dec', 'e1', 'e2',
                         'sigma_e', 'rms_e', 'weight', 'm', 'c1', 'c2'],
        required_diagnostic_inputs=list(REQUIRED),
        thresholds=dict(i_snr_min=10, resolution_min=.3, distortion_max_exclusive=2,
                        raw_sigma_min=0, raw_sigma_max=.4, dereddened_i_max=24.5,
                        blend_max_exclusive=BLEND_MAX, other_band_snr_min=5, other_bands_required=2),
        null_policy='Three-valued pass/fail/unknown. Never substitute derived sigma for raw sigma or fill missing inputs.',
        photometry='Unforced meas cmodel in all five bands; meas2 centroid and raw HSM; forced Galactic extinction only.',
        inherited_selection='Release documentation states weaklensing table is a weak-lensing-flag/Sirius selected subset. This does not supply a usable-area map for random centers.',
        flag_caveat='Table-4 cosmic-ray/bad flags were ineffective in the original pipeline; false does not establish uncontaminated pixels.',
        photoz=dict(method='Ephor_AB', purpose='point-estimate completeness/risk QA only',
                    background_cut=None, pdfs_requested=False, calibrated_beta_available=False,
                    justification='Public S16A-matched S17A COSMOS page lists Ephor_ab. This establishes a possible matched calibration route, not validated accuracy or available calibration files.'),
        calibration=dict(weight='Keep catalog inverse variance and its matched m/c calibration; no alternative weights without recalibration.',
                         sigma='raw_sigma controls the Table-4 cut; existing sigma_e is calibrated noise used in the weight identity.',
                         ensemble='Use weighted responsivity and ensemble m, rotating e and c consistently. Never divide each object by 1+m.',
                         m_tail='No automatic per-object m>-1 cut; diagnose tails against raw selection first.',
                         redshift='Same method and identical cuts on matched COSMOS, combining SOM and shape weights; point-z alone is not a dilution correction.',
                         selection='Y1-specific selection bias, PSF terms and redshift-dependent responsivity need validation; do not import S19A/Y3 coefficients or n(z).'),
        unresolved_measurement_gates=['usable S16A area and Sirius-mask geometry',
                                     'validated Y1 selection-bias and PSF null tests',
                                     'matched source-z/dilution calibration and lens-source efficiency',
                                     'redshift-dependent responsivity correction',
                                     'radial source support after identical cuts',
                                     'random overlap, spatial covariance and field influence'],
        science_measurement=False,
    )
