# S23B member outer-profile support: two-floor follow-up

The earlier three-floor analysis found no HSC-i candidates meeting the 80%
coverage requirement. This follow-up uses the previously saved, audited
50–75 proper-kpc candidate-bin diagnostics to investigate that loss and
compare 40% versus 60% coverage on shared support. No images were read and no
network request or new data transfer was made. `science_measurement=false`.

## Coverage ceiling

For each uniquely SEP-matched candidate, the joint ceiling is the smaller
of the target coverage and the third-highest coverage among its accepted
rotated control profiles; fewer than three available controls gives zero.
Passing a common floor requires both coverages to reach that floor.

| Filter | Matched candidates | Target ≥0.8 | Third control ≥0.8 | Both ≥0.8 | Median joint ceiling | Maximum joint ceiling |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| i | 2,223 | 31 | 17 | 0 | 0.1611 | 0.7892 |
| i2 | 8,730 | 434 | 406 | 113 | 0.2638 | 0.9253 |

These counts are conditional on unique source matching. Target and control
coverage can each exceed 80% in i, but not simultaneously for any candidate.
The recorded coverage alone does not identify which mask components caused
the missing area. Coverage is not a satellite completeness estimate.

## 40% and 60% fixed support

This comparison was selected **after** observing the empty 80% i intersection.
It is a post hoc sensitivity analysis. It reuses the original adjusted model,
5,000 block-bootstrap draws, six connected sky-block layouts, and the 60%
reference for frozen coverage covariates in the fixed-sample models. It is not
a new detection test. Candidate light is a target-minus-rotated-control mean;
ICL is a relative inner-minus-outer contrast.

For i, the same-cluster sample contains 62 clusters. These contain 367
retained candidates at 40% and 160 at 60%. The same-candidate sample fixes
all 160 candidates in those 62 clusters. On the 5° phase-zero sky layout:

| Fixed support | 40% slope | 60% slope | Paired change 40% − 60% | Marginal 95% block interval for change |
| --- | ---: | ---: | ---: | ---: |
| Same clusters | 0.0568 | 0.0857 | −0.0288 | [−0.2523, 0.1674] |
| Same candidates | 0.0737 | 0.0857 | −0.0120 | [−0.0695, 0.0264] |

All 36 i slope intervals and 12 paired-change intervals across the six
layouts include zero. The broad common-cluster paired interval admits
substantial changes; these results do not prove that coverage selection has
no effect or that a physical association is absent. Existing membership,
PSF wings, sky transfer, Poisson and coadd covariance gates remain open.

## Reproduction and provenance

Runners:

- `analysis/icl_wgl_membership/hsc_debug/audit_s23b_outer_support_ceiling.py`
- `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_member_two_floor_sensitivity.py`

The two-floor wrapper explicitly sets floors before using the unchanged
three-floor estimator and writes its own invocation provenance. Its output
hashes were verified. Four existing estimator tests passed, including nested
member support and a known paired slope difference.

Results:

- Ceiling: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_support_ceiling_2026-11-17_01/summary.json`, SHA-256 `573bd016e4aabec3387f42e7e53670f4a15417c8142cf62f456167a9b7ad1420`.
- Two-floor: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_two_floor_sensitivity_2026-11-17_01/summary.json`, SHA-256 `363c05231943d6d79d28912a59803a720c9579f86fb1a8bcb87b654c9929a311`.

Both use the completed 979-frame diagnostic root
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_03/`.
The S16A weak-lensing random-center query remains a separate approval-gated
request. No commit or push was performed.
