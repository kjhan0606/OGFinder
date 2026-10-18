# ICL outer-annulus center-quality sensitivity

Evaluation date: 2026-10-18 (Asia/Seoul)<br>
Status: descriptive, post hoc; `science_measurement=false`

## Method

Starting from the 224 K2.5/40% 200–300 kpc target-control pairs, I compared the
full sample with a post hoc center-eligible subset. A cluster is center-eligible
only when the automatic center check reports one unique candidate and the K2.5
mask does not cover the central pixel. The complementary group fails at least
one of those two checks. The estimate is the unweighted mean target-minus-control
surface brightness within each filter. The 95% intervals resample connected
5-degree sky blocks linking target and control positions (10,000 replicates).
No p-values or selection correction were computed.

## Results

| Filter | Group | Clusters | Joint blocks | Mean excess [spatial-block 95% interval] (nJy/arcsec²) | Leave-one-block-out mean range |
|---|---|---:|---:|---:|---:|
| i | All annulus-qualified | 71 | 18 | 3.23 [2.08, 4.23] | 2.86–3.49 |
| i | Center-eligible | 22 | 11 | 4.53 [2.31, 6.60] | 3.90–5.05 |
| i | Not center-eligible | 49 | 17 | 2.64 [1.17, 4.00] | 2.18–2.91 |
| i2 | All annulus-qualified | 153 | 40 | 3.79 [2.68, 4.82] | 3.60–4.04 |
| i2 | Center-eligible | 34 | 22 | 1.99 [−0.85, 5.03] | 1.32–3.17 |
| i2 | Not center-eligible | 119 | 38 | 4.30 [2.89, 5.58] | 4.02–4.62 |

The i estimate is higher in the center-eligible subset, but its interval overlaps
the complementary estimate. In i2 the center-eligible subset has a lower estimate
and its interval includes zero, while the complementary subset remains positive.
This filter-dependent pattern shows that a single full-sample outer excess does
not establish robustness to center quality. It does not establish that the
complementary group's excess is spurious: the flags are automatic diagnostics,
and these are post hoc subsets of the same images.

The leave-one-block-out ranges show that no single connected block reverses the
sign of any group's point estimate. They do not provide another uncertainty
interval or remove shared measurement systematics. In particular, the i2
center-eligible range remains positive while its bootstrap interval includes zero,
because the interval reflects spatial resampling rather than single-block deletion.

## Limits and use

This is a sensitivity diagnostic for planning manual center review. It does not
validate BCG identity, remove ICL, certify the sky or PSF transfer, model the
coadd pixel covariance, or correct for center-selection effects. The filter
profiles remain separate, and the result is not a science detection. The next
useful check is blinded visual review of the center and residual-source masks,
followed by a predeclared center-quality cut and a repeat on untouched validation
regions.

Reproduction code: [`audit_icl_center_quality.py`](audit_icl_center_quality.py).
Input and output hashes are recorded in
`/gpfs/kjhan/Astrafak/data/sample_800_2026-10-13/icl_center_quality_sensitivity_2026-10-18_02/summary.json`.
