# S23B member mask: paired cluster-level comparison

## Scope and method

This offline analysis uses the previously declared 96-frame spatial sample
(48 i, 48 i2), saved per-cluster mask fractions, and the pinned relative ICL
profile tables. It reads no image pixels and makes no network request. Each
pair is **one cluster's pooled member-centered target fraction and pooled
rotated-control fraction**. The controls share the image and are not treated
as independent clusters or as one-to-one matched stamps.

The estimated difference is control minus target, equally weighted across
paired clusters. Each 5° sky-cell bootstrap draw carries every paired cluster
in a drawn cell together and is shared across components and contexts within
that filter. The reported 95% intervals use 5,000 draws over 33 sampled i
cells or 48 sampled i2 cells. This deliberately sky-balanced selection is
not population weighted; intervals describe these sampled cells, not an
unbiased survey-population interval. Mask footprints overlap and cannot be
added. All results are exploratory with `science_measurement=false`.

## Paired mask differences

| Filter | Component | Mean control − target | Marginal 95% sky-cell interval |
| --- | --- | ---: | ---: |
| i | Selected quality-bit footprint | +9.56 percentage points | [5.43, 14.33] |
| i | Bright-object footprint | −0.77 pp | [−5.06, 3.46] |
| i | Neighbor-source-shape footprint | −4.94 pp | [−7.38, −2.60] |
| i | Unmasked area | −1.35 pp | [−3.73, 0.93] |
| i2 | Selected quality-bit footprint | +3.18 pp | [0.56, 6.47] |
| i2 | Bright-object footprint | +2.08 pp | [−2.40, 6.45] |
| i2 | Neighbor-source-shape footprint | −3.39 pp | [−5.78, −1.15] |
| i2 | Unmasked area | −0.36 pp | [−3.35, 2.49] |

These are marginal intervals for several correlated components, without
multiplicity adjustment. They show that the target/control quality-bit and
neighbor-shape differences deserve explicit modeling in this selected sample;
they do not quantify resulting photometric bias or justify removing flags.
The unmasked-area intervals contain zero despite differing component footprints,
consistent with overlap and compensation among components. This does not
prove the total selection is balanced in the parent population.

## Exploratory ICL and support associations

The same paired differences were compared with relative inner-minus-outer ICL
contrast and the minimum common ICL pixel-support fraction using Spearman rank
correlations and the same sky-cell draws. These are unadjusted, marginal and
post-selection summaries, without covariate adjustment or multiple-comparison
control. Among 24 component/filter/context correlations, one marginal interval
excludes zero: i neighbor-shape difference versus relative ICL contrast has
rho 0.338 [0.009, 0.611]. It is a hypothesis-generating result only and is not
evidence for a physical ICL–satellite relation or a selected primary finding.
The other 23 marginal intervals include zero. Full results are in the table.

The mask geometry, relative ICL contrast and member profiles share images;
this can create measurement coupling. The comparison does not remove such
coupling or solve PSF, sky transfer, source Poisson, coadd covariance or
membership uncertainty. No causal attribution or calibrated science claim is
made.

## Verification and provenance

The analysis joins the exact 96 planned frames, verifies their sky-cell IDs,
and pins the cluster summary, sample plan and two ICL context tables. All 48
clusters per filter have both roles. A synthetic tied-rank comparison matched
SciPy's Spearman statistic, a constant rank vector was rejected, and a
synthetic bootstrap confirmed that all clusters in each drawn sky cell remain
together. Output file hashes were independently checked.

- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_member_mask_paired.py`.
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_paired_2026-11-17_01/`.
- `summary.json` SHA-256: `32d0cb99149a5234bc630c369703923f45646ec2fe41f6f0776ac7144597d33c`.
- `paired_cluster_fractions.tsv` SHA-256: `bd8646be837d978f8bb6ab88ba1c31bc30c80e4b6a0d94462a55f06a3318624f`.
- `paired_mask_results.tsv` SHA-256: `51bb1e370c77d11719bd9cd5cb4b389832df31e73f6a47e69e4c69c59a2f1ffe`.

A useful next gate is a fixed-sample analysis of how observed mask differences
change the measurement response under injected signals, with PSF convolution
and source Poisson handled explicitly. Such calibration must preserve target
and control selection and must not infer causation from the rank correlations.
The separate S16A random-center count query remains approval-gated. No commit
or push was performed.
