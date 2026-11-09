# Random-subtracted weak-lensing influence by connected HSC tract region

## Scope and method

This offline analysis re-evaluates the existing HSC S19A cluster-minus-random
stack after deleting each approximate connected HSC tract-center component.
The analysis uses the already approved 435-lens common-support sample, 2,175
matched random centers, 2,075,062 target lens-source pairs, and 8,905,892
random-center-source pairs. The calibrated S19A estimator, B-mode mask,
tomographic `n(z)`, radial bins, and cosmology match the saved stack.

Astrafek groups the 840 tracts in the saved HSC patch-QA table by proximity of
spherical-mean tract centers at a 2-degree adjacency threshold. The graph has
four connected components. Astrafek assigns each lens, random center, and
source to its nearest tract center and deletes a pair when either endpoint
belongs to the component under deletion. The analysis recomputes the complete
calibrated estimator for both the target and random samples after every
deletion.

The tract-center graph approximates the survey footprint. The graph does not
use tract polygons or official field labels. The leave-one-region-out results
measure regional influence. Three occupied components do not support a stable
full covariance estimate.

## Results

The full paired estimator reproduces all seven saved random-subtracted
tangential and cross profile values exactly at the stored precision. The
random-subtracted tangential profile is positive and declines from
`16.005 × 10^13` to `1.679 × 10^13 M☉ Mpc⁻²` between 0.1–0.2 and 2–3 proper
Mpc. Removing one connected region changes the tangential profile as follows.

| Proper radius (Mpc) | Full ΔΣ (10¹³ M☉ Mpc⁻²) | Region 0 removed | Region 1 removed | Region 2 removed |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | 16.005 | −28.6% | +12.1% | +3.7% |
| 0.2–0.4 | 11.175 | −10.0% | +4.3% | +1.3% |
| 0.4–0.7 | 7.326 | −26.8% | +10.9% | +3.9% |
| 0.7–1.0 | 5.172 | −18.7% | +10.1% | +1.0% |
| 1.0–1.5 | 3.936 | −6.7% | +1.6% | +1.9% |
| 1.5–2.0 | 2.677 | +0.7% | −0.2% | −0.2% |
| 2.0–3.0 | 1.679 | −11.3% | +5.0% | +1.5% |

Region 0 contains 261 of 435 common-support lenses, or 60.0% of the sample.
Deleting region 0 shifts the 0.4–0.7 Mpc signal by 3.29 times the corresponding
5-degree jackknife diagonal error. The 0.1–0.2 Mpc shift is 3.01 times that
error. Region 1 contains 141 lenses and produces a maximum 1.33-error shift.
Region 2 contains 33 lenses and produces a maximum 0.48-error shift. The
regions therefore have strongly unequal statistical weight.

An offline topology check repeated the tract graph at adjacency thresholds of
1.6°, 2.0°, 2.4°, and 3.0°. All four settings produce the same exact partition
of the 840 tract IDs and retain the same common-lens counts of 261, 141, and
33, plus the same random-center counts of 1,305, 705, and 165. The large
leave-one-region-out shifts therefore do not arise from the tested adjacency
threshold choice. The check does not validate the approximate nearest-center
assignment or replace polygon-defined field geometry.

## Interpretation

The positive declining point profile is reproducible, but the innermost and
intermediate-radius amplitudes depend strongly on the dominant connected
region. The 5-degree equal-area jackknife does not capture this influence
adequately. The three occupied components also provide too few independent
regions for a defensible region-jackknife covariance. The lensing amplitude
and its formal significance therefore remain exploratory. `science_measurement`
remains `false`.

The result does not establish a halo-mass measurement or a cluster detection.
The approximate tract assignment, random-center full-aperture matching,
photo-z calibration uncertainty, reduced shear, magnification, shear
calibration, and PSF/selection systematics remain unresolved. A calibrated
survey-region covariance requires validated field geometry and adequate
independent-region support.

## Reproducibility

- Astrafek runner: `analyze_wl_random_tract_region_influence.py`
- Runner SHA-256: `1446b31779fb9b15c88e825afe135d8aaab621e60cdfd2a52d2c8cb271a7b27e`
- Result summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_region_influence_2026-11-06_01/summary.json`
- Summary SHA-256: `ae47953cb7af1dc89070c8b486e522dfa9ac088fd9761ed7cb1f6ce43494d343`
- Reference paired stack: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json`
- Reference stack SHA-256: `83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8`
- Plan hashes, catalogue hashes, estimator settings, and component-level row counts are stored in the result summary.
- The analysis uses existing HSC exports only. No CAS query or download was performed.
- Partition audit runner: `audit_wl_random_tract_partition_sensitivity.py`, SHA-256 `c7bb5872aef6575150926ce20db3aeaa45ab00098b3dd3ab55b9f75d30467c6b`.
- Partition audit summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_partition_sensitivity_2026-11-06_01/summary.json`, SHA-256 `da77b4f3d9db6af7b514259156b342ced14e61874ffaab564b2e9ab061ce1ebb`.
