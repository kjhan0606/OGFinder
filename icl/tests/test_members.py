"""Spectroscopic membership window for ICL masks and weak-lensing sources."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from icl.members import (  # noqa: E402
    BCG, BACKGROUND, SATELLITE, STAR, UNKNOWN, background_for_lensing,
    classify_members, classify_spectroscopic, diffuse_stack_masked, fraction_kept,
    object_mask, photometric_redshift_labels)
from ogfkit.wlcluster import Cosmology, select_sources  # noqa: E402


def test_window_splits_satellites_from_background_and_leaves_gaps_unknown():
    z = np.array([0.30, 0.33, 0.35, 0.50, np.nan, 0.29])
    labels = classify_spectroscopic(z, 0.30, 0.02)
    assert list(labels) == [SATELLITE, BACKGROUND, BACKGROUND, BACKGROUND, UNKNOWN, SATELLITE]
    # The boundary itself is outside the window.
    assert classify_spectroscopic([0.32], 0.30, 0.02)[0] == BACKGROUND
    assert classify_spectroscopic([0.319], 0.30, 0.02)[0] == SATELLITE


def test_star_and_chosen_centre_win_over_the_redshift():
    z = np.array([0.30, 0.30, 0.80, np.nan])
    star = np.array([False, True, False, True])
    labels = classify_spectroscopic(z, 0.30, 0.05, bcg_index=0, star=star)
    assert list(labels) == [BCG, STAR, BACKGROUND, STAR]
    # A centre with no redshift is still the BCG.
    assert classify_spectroscopic([np.nan], 0.30, 0.05, bcg_index=0)[0] == BCG


def test_diffuse_mask_drops_resolved_sources_and_the_fraction_keeps_members():
    labels = np.array([BCG, SATELLITE, BACKGROUND, STAR, UNKNOWN])
    assert list(diffuse_stack_masked(labels)) == [False, True, True, True, True]
    assert list(fraction_kept(labels)) == [True, True, False, False, False]
    seg = np.array([[0, 1, 2], [3, 4, 5]])
    diffuse = object_mask(seg, labels, 'diffuse')
    assert diffuse.tolist() == [[False, False, True], [True, True, True]]
    fraction = object_mask(seg, labels, 'nonmember')
    assert fraction.tolist() == [[False, False, False], [True, True, True]]


def test_lensing_keeps_only_background_behind_the_gap():
    labels = np.array([SATELLITE, BACKGROUND, BACKGROUND, UNKNOWN])
    z = np.array([0.31, 0.35, 1.0, 1.2])
    keep = background_for_lensing(labels, z, 0.30, 0.1)
    assert list(keep) == [False, False, True, False]
    cosmo = Cosmology()
    x = np.zeros(4)
    y = np.array([100.0, 200.0, 300.0, 400.0])
    g1 = np.full(4, 0.02)
    g2 = np.zeros(4)
    w = np.ones(4)
    chosen = select_sources(x, y, g1, g2, z, w, 0.30, cosmo, z_gap=0.1, rmin=0.01, rmax=20.0,
                            member=labels)
    assert list(chosen['index']) == [2]
    untouched = select_sources(x, y, g1, g2, z, w, 0.30, cosmo, z_gap=0.1, rmin=0.01, rmax=20.0)
    assert list(untouched['index']) == [2, 3]


def test_rejects_a_bad_window_or_a_short_table():
    with pytest.raises(ValueError):
        classify_spectroscopic([0.3], 0.3, 0.0)
    with pytest.raises(ValueError):
        classify_spectroscopic([0.3], -0.1, 0.05)
    with pytest.raises(ValueError):
        classify_spectroscopic([0.3, 0.4], 0.3, 0.05, bcg_index=2)
    with pytest.raises(ValueError):
        classify_spectroscopic([0.3, 0.4], 0.3, 0.05, star=[True])
    with pytest.raises(ValueError):
        object_mask(np.array([[1, 9]]), ['bcg'], 'diffuse')
    with pytest.raises(ValueError):
        object_mask(np.zeros((1, 1), int), ['bcg'], 'other')
    cosmo = Cosmology()
    with pytest.raises(ValueError):
        select_sources([0.0], [10.0], [0.0], [0.0], [1.0], [1.0], 0.3, cosmo,
                       rmin=0.01, rmax=5.0, member=['background', 'background'])


def _line(mag):
    mag = np.asarray(mag, float)
    return 1.10 - 0.03 * (mag - 20.0)


def test_point_photoz_uses_one_window_and_treats_overlap_as_unknown():
    z_phot = np.array([0.30, 0.80, 0.10, 0.50, np.nan])
    labels = photometric_redshift_labels(z_phot, (0.28, 0.34), (0.50, 1.20))
    assert list(labels) == [SATELLITE, BACKGROUND, UNKNOWN, BACKGROUND, UNKNOWN]
    # The upper edge is outside the half-open window. Overlap is unknown.
    assert photometric_redshift_labels([0.34], (0.28, 0.34), (0.50, 1.20))[0] == UNKNOWN
    assert photometric_redshift_labels([0.40], (0.28, 0.55), (0.40, 1.20))[0] == UNKNOWN


def test_pz_decides_a_row_and_an_unusable_row_uses_the_point_estimate():
    grid = np.linspace(0.0, 2.0, 401)
    peak = np.exp(-0.5 * ((grid - 0.30) / 0.008) ** 2)
    far = np.exp(-0.5 * ((grid - 1.00) / 0.008) ** 2)
    split = np.exp(-0.5 * ((grid - 0.30) / 0.008) ** 2) + np.exp(-0.5 * ((grid - 1.00) / 0.008) ** 2)
    pz = np.vstack([peak, far, split, np.full(grid.size, np.nan)])
    labels = photometric_redshift_labels(
        [np.nan, np.nan, np.nan, 0.30], (0.20, 0.40), (0.70, 1.40),
        pz=pz, z_grid=grid, p_min=0.8)
    assert list(labels) == [SATELLITE, BACKGROUND, UNKNOWN, SATELLITE]


def test_red_sequence_adds_satellites_and_leaves_the_spectroscopic_label():
    mag = np.array([18.0, 19.0, 20.0, 21.0, 22.0, 22.5, 20.0, 21.0, 21.0, 23.8, 21.0, 18.5])
    color = _line(mag)
    color[6] = _line(20.0) - 0.50
    z = np.array([0.30, 0.31, 0.29, 0.30, 0.31, 0.29, 0.30, 0.80, np.nan, np.nan, np.nan, np.nan])
    # Rows 0-5 sit on the line with a cluster redshift. Row 6 is a blue satellite.
    # Row 7 is a spectroscopic background on the line. Rows 8-11 have no redshift:
    # on the line, fainter than the limit, just outside the width, and a star.
    color[10] = _line(21.0) + 0.15
    out = classify_members(
        z, 0.30, 0.02, bcg_index=0, star=[False] * 11 + [True],
        color=color, mag=mag, rs_width=0.15, mag_limit=23.0, min_satellites=4)
    labels = list(out['labels'])
    assert labels[:7] == [BCG, SATELLITE, SATELLITE, SATELLITE, SATELLITE, SATELLITE, SATELLITE]
    assert labels[7:] == [BACKGROUND, SATELLITE, UNKNOWN, UNKNOWN, STAR]
    fit = out['red_sequence']
    assert fit['ok'] and fit['n'] == 4
    assert fit['a'] + fit['b'] * (21.0 - fit['m0']) == pytest.approx(_line(21.0), abs=1e-8)
    # A spectroscopic redshift is kept when the photo-z and the colour disagree.
    kept = classify_members(
        [0.80], 0.30, 0.02, z_phot=[0.30], cluster_window=(0.28, 0.34),
        background_window=(0.6, 1.2), color=[_line(20.0)], mag=[20.0],
        rs_width=0.08, mag_limit=23.0, min_satellites=2)
    assert kept['labels'][0] == BACKGROUND


def test_photoz_fills_only_a_missing_redshift_and_pz_can_make_a_satellite():
    grid = np.linspace(0.0, 2.0, 401)
    peak = np.exp(-0.5 * ((grid - 0.31) / 0.008) ** 2)
    z = np.array([0.50, np.nan, np.nan])
    z_phot = np.array([0.30, 0.90, np.nan])
    pz = np.vstack([np.full(grid.size, np.nan), np.full(grid.size, np.nan), peak])
    color = np.array([_line(20.0), _line(20.0), _line(20.0) - 0.4])
    out = classify_members(
        z, 0.30, 0.02, z_phot=z_phot, cluster_window=(0.28, 0.34),
        background_window=(0.70, 1.30), pz=pz, z_grid=grid, p_min=0.8,
        color=color, mag=np.array([20.0, 20.0, 20.0]), rs_width=0.08,
        mag_limit=23.0, min_satellites=2)
    # No spectroscopic satellite, so the red sequence does not relabel the blue P(z) row.
    assert out['red_sequence']['ok'] is False
    assert list(out['labels']) == [BACKGROUND, BACKGROUND, SATELLITE]


def test_photometric_arguments_are_rejected_when_incomplete():
    with pytest.raises(ValueError):
        photometric_redshift_labels([0.3], (0.4, 0.2), (0.5, 1.0))
    with pytest.raises(ValueError):
        photometric_redshift_labels(pz=np.ones((1, 4)), z_grid=np.linspace(0, 1, 4),
                                    cluster_window=(0.2, 0.4), background_window=(0.6, 1.0))
    with pytest.raises(ValueError):
        classify_members([0.3, np.nan], 0.3, 0.02, color=[1.0, 1.0], mag=[20.0, 21.0])
    with pytest.raises(ValueError):
        classify_members([np.nan], 0.3, 0.02, z_phot=[0.3], cluster_window=(0.2, 0.4),
                         background_window=(0.5, 1.0), color=[1.0, 1.1], mag=[20.0],
                         rs_width=0.1, mag_limit=23.0, min_satellites=2)
