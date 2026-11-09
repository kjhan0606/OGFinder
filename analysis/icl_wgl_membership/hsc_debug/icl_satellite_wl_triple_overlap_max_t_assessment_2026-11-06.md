# Familywise ICL–satellite–weak-lensing association audit (2026-11-06)

## Question

The existing triple-overlap report found two pointwise 95% spatial-block
intervals excluding zero among 29 reported associations. This analysis asks
whether either pattern remains after accounting for the full 29-association
family, using the same selected 27-lens `i2` common-support sample and existing
HSC-derived tables only.

## Method

For every relationship, the point estimate is the same partial Pearson
correlation as in the original report. Tests use CR1 cluster-studentized OLS
slopes with the same conditional covariates. The null is generated with a
restricted Webb six-point wild-cluster bootstrap (10,000 replicates). Each
replicate uses one shared Webb multiplier per connected 5-degree target/control
sky block across all 29 relationships, preserving dependence among WL
components, radial bins, and the two predictors. The maximum absolute
studentized statistic across the family supplies the max-T adjusted p-value.
The shared partition has 20 connected blocks and at most two lenses per block.
These null-based wild-cluster p-values are a complementary inference method,
not a multiplicity-adjusted version of the original percentile-bootstrap
intervals; differences between their pointwise conclusions can occur.

## Findings

| Relationship | Lenses | Partial r | Wild-cluster p | 29-way max-T p |
|---|---:|---:|---:|---:|
| WL cross vs satellite contrast, 0.4–0.7 proper Mpc | 25 | −0.517 | 0.043 | 0.672 |
| WL tangential vs satellite contrast, 2.0–3.0 proper Mpc | 27 | +0.412 | 0.070 | 0.923 |
| Direct ICL vs satellite contrast | 27 | +0.067 | 0.460 | 1.000 |

No relationship survives the familywise adjustment. The original pointwise
patterns are not evidence for a resolved ICL–satellite–halo relation after
accounting for the 29 tested associations. This does not establish the absence
of a physical relation: the common sample contains only 27 selected lenses and
the block count is small.

All 29 partial-correlation point estimates reproduce the original metrics
exactly (maximum absolute numerical difference 0). The bootstrap is exploratory
and inherits approximate 5-degree sky blocks, the selected support, and
unresolved WL photo-z, shear calibration, magnification, mask, PSF, and
covariance systematics. It does not repair ICL sky transfer, source-Poisson,
PSF-wing, or image-covariance uncertainties. Keep `science_measurement=false`.

## Reproducibility

- Runner: `adjust_icl_satellite_wl_triple_overlap_max_t.py`, SHA-256
  `162c90585875d0811f6921407c77ec72837f6ecfb5fd9e571bf2c5bbde5e156d`.
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_maxt_2026-11-06_01/`.
- `summary.json` SHA-256: `d5e2d3aa05be32b07be5fa62b8de467aea5314411daba49cb733444850126b0e`.
- `max_t_metrics.tsv` SHA-256: `8f58ec651bf09d6734aba344a77ef52613c0351b68e4d10f8b65f03333dde59f`.
- Source lens table, metric table, and summary hashes are recorded in the output
  `summary.json`; no HSC request or download occurred.
