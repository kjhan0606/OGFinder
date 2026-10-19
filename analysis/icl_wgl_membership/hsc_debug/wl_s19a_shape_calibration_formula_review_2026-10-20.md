# HSC S19A shear-calibration formula review — 2026-10-20

## Scope

This is a source-level audit of Astrafek's existing HSC S19A estimator. It
uses the HSC Y3 shear-catalog paper and the calibration utility linked from the
official HSC release page. It does not access the CAS, transfer data, or
recompute a profile.

## Formula correspondence

The Astrafek implementation in
[`s19a_shape_stack.py`](s19a_shape_stack.py) matches the HSC prescriptions in
the following points:

- The weighted responsivity is `R = 1 - sum(w e_RMS^2)/sum(w)`, recomputed
  within each projected-radius bin using the lensing pair weight
  `shape_weight * Sigma_crit^{-2}`.
- Per-object multiplicative and additive calibrations use the released S19A
  `derived_shear_bias_m`, `derived_shear_bias_c1`, and
  `derived_shear_bias_c2` fields.
- The magnitude and resolution edge-density terms use the HSC bands
  `mag_A >= 25.5 - 0.025` and `R2 <= 0.3 + 0.01`, normalized by their band
  widths and by the selected-sample weight sum. Astrafek's coefficients for
  `m_sel`, `a_sel`, and their quoted errors match the linked HSC
  `get_sel_bias` utility. The additive correction is projected from
  `a_sel * <e_PSF>` using the same pair weights.
- The applied source-quality cuts `|e|<2`, `0<=sigma_e<=0.4`, `R2>=0.3`, and
  `mag_A<=25.5` match Table 2 of the HSC Y3 shear-catalog paper.

The relevant code is in
[`hsc_selection_bias`](s19a_shape_stack.py) and
[`stack_s19a_pairs`](s19a_shape_stack.py). These checks validate formula
correspondence; they do not validate the new DESI redshift distributions,
random-lens subtraction, reduced shear, field covariance, or a measured lens
signal. `science_measurement=false` remains required.

One endpoint detail is recorded for future audits: the linked utility's
`get_wl_cuts` helper writes strict comparisons for `sigma_e` at 0 and 0.4,
whereas Table 2 states inclusive bounds. Astrafek follows the paper table and
the S19A catalog documentation. The catalog is selected through HSC's
weak-lensing table, and no rows were downloaded or reselected in this
source-level review; endpoint counts have not been independently recomputed
from the row files.

## Sources

- Li et al., [The three-year shear catalog of the Subaru Hyper Suprime-Cam SSP
  Survey](https://arxiv.org/abs/2107.00136), especially the responsivity
  equation and Table 2 cuts.
- Official HSC [S19A Shape Catalog (PDR3) documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).
- HSC-linked [S19A calibration utility](https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py),
  `get_sel_bias` and `get_wl_cuts`.

## Scientific gate

The estimator implementation passes this formula correspondence review. The
observed source-bin consistency test still fails in inner radial bins and has
not been repeated with an updated calibrated `n(z)`. The existing cluster-only
profile is not promoted to a physical detection until the approved random
controls, calibrated redshift-distribution uncertainty, paired source-overlap
covariance, and remaining shear systematics are resolved.
