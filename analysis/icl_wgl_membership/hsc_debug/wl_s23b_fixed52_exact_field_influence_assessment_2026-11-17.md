# Fixed-52 weak-lensing influence from official HSC footprint fields

## Method and validation

The fixed 52 lens sets and their five planned random centers per lens are the
same approved S19A shape-row and HSC Y3 n(z) inputs used by the preceding exact
tract calculation. Existing official PDR3 Wide patch polygons assign each
active lens/control and source endpoint. Tracts are then grouped into their
parent official footprint field. A field deletion removes an individual pair
if **either endpoint** lies in that field; calibration, random-center
averaging and the requirement of at least four positive control centers in
each radial bin are recalculated. This is a leave-field-out *influence*
calculation, not a survey covariance estimator. The four local geometry files
were checked against the preceding exact-tract hashes. All three full-sample
reference estimators agree exactly with the prior exact-tract outputs, and all
output file hashes were independently checked. No CAS, network, new download,
or image read was involved. `science_measurement=false`.

## Results

The sample occupies three fields: AUTUMN (21 lens centers), HECTOMAP (3) and
SPRING (28). The equal-random-center full-sample tangential profile is 0.7726
in the 1.5–2.0 proper-Mpc bin, in units of 10^13 solar masses per Mpc².

| Omitted field | Remaining supported lens sets | 1.5–2.0 Mpc profile | Shift from full sample |
|---|---:|---:|---:|
| AUTUMN | 31 | 1.5422 | +0.7696 |
| HECTOMAP | 49 | 0.6684 | −0.1042 |
| SPRING | 24 | −0.0087 | −0.7813 |

The two large-field exclusions shift that bin by approximately +99.6% and
−101.1% of the full value, respectively. This confirms the previously observed
large-scale sensitivity using official patch-derived field labels instead of
nearest-patch-center connected components. It also incorporates both lens and
source endpoint deletions. In the inner 0.1–0.2 Mpc bin the full tangential
profile is 24.2044; field-exclusion shifts are +1.6311 (AUTUMN), −1.2835
(HECTOMAP), and +0.5136 (SPRING) in the same units. The SPRING-exclusion inner
cross-component shift is −10.6462 in the same units. These are raw influence
shifts, not independent error bars or significance estimates.

Each deletion retains at least four control centers for every included lens set
and radius, but the common supported set changes. Consequently these shifts
combine field dependence and the effect of changing the lens sample. The three
occupied fields cannot supply a stable full survey covariance. The sample is
also selected for support; n(z) uncertainty, source reuse, PSF/mask matching,
shear-selection response, contamination/boost and other systematics remain.
Do not interpret this as a mass fit or physical detection.

- Runner: `analyze_s23b_fixed52_tract_region_influence.py --exact-group field`
- Grouping tests: `test_s23b_wl_exact_field_grouping.py` (plus S19A estimator tests)
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_exact_field_influence_2026-11-17_01/`
- Summary SHA-256: `422afe04d11cd10f5bb047f02dae35d0990a1890615ab7b5c49e3abb716c7d33`
- The original tract-level result and the field-level run remain distinct.
