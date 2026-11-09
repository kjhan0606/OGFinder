# ICL--satellite association: evidence status

This report consolidates saved exploratory diagnostics for the fixed
ICL--satellite association sample. It does not reread image pixels or
submit HSC queries. All component summaries retain
`science_measurement=false`; the numerical intervals below are
diagnostic and do not establish a detection or causal relation.

## Main robustness comparisons

The first table uses the same pair sample and z/richness-adjusted
partial-correlation estimator for the raw and finite-lag pixel-noise
results. The shell-exclusion values are a paired aperture-overlap
sensitivity: supported satellite 50--75 kpc annulus pixels are removed
from the target ICL 200--300 kpc ring; matched-control ICL and satellite
measurements remain fixed.

| Filter | Pairs | Raw partial r (95% block interval) | Finite-lag pixel-noise corrected r (95% interval) | After shell-pixel exclusion (95% interval) | Paired change, excluded minus raw (95% interval) |
|---|---:|---:|---:|---:|---:|
| i | 21 | 0.120 [-0.227, 0.522] (median 0.130) | 0.122 [-0.232, 0.555] (median 0.133) | 0.082 [-0.240, 0.518] (median 0.101) | -0.037 [-0.084, -0.003] (median -0.036) |
| i2 | 53 | 0.234 [0.016, 0.416] (median 0.227) | 0.249 [0.018, 0.449] (median 0.242) | 0.193 [-0.028, 0.379] (median 0.191) | -0.041 [-0.064, -0.021] (median -0.039) |

The second table treats the measured 600--900 kpc target-minus-control
plane intercept as a residual-sky diagnostic covariate. This annulus
can contain real cluster light and is not subtracted from the ICL
measurement.

| Filter | Pairs | Baseline partial r (z, richness) | With outer-residual covariate | Paired change (95% block interval) |
|---|---:|---:|---:|---:|
| i | 21 | 0.112 | 0.107 | [-0.137, 0.106] |
| i2 | 53 | 0.239 | 0.252 | [-0.034, 0.111] |

## Selection and sample-size limits

A separate support-selection audit has 224 ICL-supported rows; 94 also
have the required satellite-profile support (26 in i, 68 in i2). The
correlation tables above use the narrower 74-pair sample (21 in i, 53
in i2) retained by their diagnostic inputs. These counts describe
different analysis gates and should not be interchanged.

The supported rows are selected: in i2, association rows differ from
the non-associated ICL-supported rows by +1.055 SD in target image
coverage, +0.318 SD in redshift, and -0.325 SD in richness. The audit
describes selection; it does not correct it or make the sample
representative of all 800 parent clusters.

The following are design scenarios, not observed-effect claims. They
assume 80% power, two covariates, and conservative Bonferroni control
for 168 hypotheses. Parent-cluster projections hold the observed i2
support yield fixed; design effects are sensitivity assumptions.

| Assumed partial r | Design effect | Required i2 usable rows | Projected parent clusters |
|---:|---:|---:|---:|
| 0.200 | 1.000 | 536 | 6,306 |
| 0.200 | 1.500 | 804 | 9,459 |
| 0.293 | 1.000 | 323 | 3,800 |
| 0.293 | 1.500 | 485 | 5,706 |

The offline S23B catalog audit gives at most 1,333 additional catalog
candidates (2,133 parent
clusters total), projecting about 181.305 i2 support rows
if the current yield transferred. This remains below every displayed
family-wise planning target. It is a catalog-only upper bound: candidate
PDR3 images, masks, depth, PSF, footprint, and matched controls are not
verified, and no HSC query or new image download is represented here.
The r=0.293 row is a post-selection exploratory estimate and is not a
confirmatory effect-size assumption.

## Interpretation and limits

- The current association sample has 74 pairs (21 in i and 53 in i2).
- The i2 point estimate is positive in these diagnostics, but the
  supported-shell exclusion attenuates it; the paired change interval
  is an exploratory sensitivity with only 28 spatial blocks. The i-band
  estimate is much less precise (21 pairs; 13 blocks).
- The finite 9x9 lag-kernel correction covers only the supplied
  delivered-variance noise model. It is not a correction for large-scale
  sky modes, PSF wings, selection, completeness, extinction, or calibration.
- The shell-pixel test is not a deblending or PSF correction. It removes
  supported pixels from the ICL aperture while keeping the satellite
  estimator and matched-control ICL fixed.
- The outer annulus is a residual diagnostic, not clean sky; no flux was
  subtracted. Its paired change interval includes zero in both filters.
- The coadd-header audit found no gain or exposure metadata in the 1,598
  target/control images. A physically based injected-source Poisson
  variance and completeness correction therefore remains unavailable
  from the currently saved coadds and variance planes alone.
- These results support a tentative i2 association worth following up,
  not a confirmed ICL--satellite relation. Connected HSC fields, PSF
  wings, sky transfer, source completeness, and noise covariance need
  stronger validation before a science claim.

## Provenance

Bootstrap intervals use 5,000 replicates in their component analyses.
Source-summary SHA-256 values:

- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_sky_association_reblocked_2026-11-03_01/summary.json`: `89c0f39d72f76a8f4864a68f8ca4a09fa05fe0d047f937c888055c81b5c8c983`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_poisson_noise_input_audit_2026-11-04_01/summary.json`: `f75e704de3ae8dbdfc6c2f52700bebb9407e13532d0f1dd10eb6b71990fd533c`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_expansion_feasibility_2026-10-21_02/summary.json`: `e482ebd640230126b5d7d7e4438b78609b22899728765e512cc43b6bd21ab4e5`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sample_size_plan_2026-10-21_01/sample_size_plan.json`: `df3d6d06adecaa310ab94df94bbebd1441262087520c16385927f285846f0b88`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_association_selection_2026-10-20_03/summary.json`: `c3c552aa1a4ca1d47f52cef2c90120d52f96250ed21419ad278eca03aa016fc7`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_sky_covariance_2026-11-03_01/summary.json`: `947d9ecd8bb728e5a3bb1ec193f88687a893752403c9c7a4d409a3e6752888e2`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_overlap_sensitivity_2026-11-03_01/summary.json`: `733433f2495ef4cf1592d7885a9845de57c726b6b64a2fdf93ca91514a0260e9`
- `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_pixel_noise_bias_2026-11-03_01/summary.json`: `268c5ee1a420e9cd97c479f8fe5ea8000bd3afaac41c7469ff3a1028bfed8396`
- report builder `/home/kjhan/BACKUP/Astrafek/OGFinder/analysis/icl_wgl_membership/hsc_debug/build_icl_association_evidence_report.py`: `fcd1fe96029a141c882242715ff0e3d12ece6261a119568cbc4171e1885f9bf6`

This consolidation reads saved JSON summaries only; it makes no network
request, does not alter image data, and leaves the science gate false.
