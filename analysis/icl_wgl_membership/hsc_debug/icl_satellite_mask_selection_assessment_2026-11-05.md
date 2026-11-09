# Mask-driven selection in the ICL–satellite sample

This offline analysis asks whether candidates with usable 50–75 kpc annuli are
selected differently by cluster-centric radius, SEP detection S/N, cluster
redshift, or richness. It uses 2,691 candidates in 224 cluster/filter frames
with ICL 200–300 kpc support. Quartile contrasts are computed within filter;
95% intervals resample whole cluster IDs (5,000 replicates). No image pixels or
HSC service were accessed. This is a mask-selection diagnostic, not an ICL
measurement; `science_measurement=false`.

## Results

The share passing the 40% outer-annulus target-coverage floor increases from
the lowest to highest cluster-centric-radius quartile by `+0.265` in `i`
(95% cluster-bootstrap interval `[+0.154,+0.374]`) and `+0.267` in `i2`
(`[+0.189,+0.348]`). The retained-profile fraction, which also requires at
least three valid azimuthal controls, increases by `+0.227` in `i`
(`[+0.106,+0.351]`) and `+0.206` in `i2` (`[+0.127,+0.282]`). Thus the
available outer-annulus measurements are preferentially drawn from satellites
farther from the cluster center in both filters.

The highest-richness quartile has a lower coverage-pass fraction than the
lowest-richness quartile in `i` by `−0.237` (`[−0.341,−0.119]`); the retained
profile contrast is `−0.199` (`[−0.295,−0.092]`). The corresponding `i2`
contrasts are smaller and their intervals include zero: coverage `−0.077`
(`[−0.164,+0.007]`) and retained profiles `−0.073` (`[−0.170,+0.024]`).

Detection S/N has no evident `i` coverage contrast (`0.000`,
`[−0.099,+0.097]`). In `i2`, the highest-minus-lowest S/N quartile contrast is
`−0.075` (`[−0.150, approximately 0]`) for both coverage passage and retained
profiles. Redshift contrasts are consistent with zero in both filters. These
are separate, unadjusted quartile comparisons; they do not isolate a causal
effect of any one variable.

An adjusted descriptive logistic model gives the same main pattern. Predictors
are standardized within filter; S/N and richness use `log1p(max(S/N,0))` and
`ln(richness)`. Intervals use a cluster-ID sandwich covariance and a `t`
critical value (`G−1` degrees of freedom). All four models converged with
gradient infinity norms below `7.3×10⁻⁷`; information-matrix condition numbers
are 2.49–2.70.

| Outcome and filter | Radius OR per SD (95% CI) | Log-richness OR per SD (95% CI) | Other adjusted predictor |
| --- | --- | --- | --- |
| Coverage pass, `i` | 1.592 (1.317–1.925) | 0.665 (0.561–0.787) | log S/N and redshift intervals include 1 |
| Coverage pass, `i2` | 1.556 (1.379–1.755) | 0.816 (0.726–0.916) | log S/N 0.871 (0.778–0.975); redshift interval includes 1 |
| Profile retained, `i` | 1.548 (1.246–1.924) | 0.642 (0.509–0.808) | log S/N and redshift intervals include 1 |
| Profile retained, `i2` | 1.448 (1.277–1.641) | 0.807 (0.693–0.940) | log S/N 0.876 (0.781–0.984); redshift interval includes 1 |

The radial and richness associations remain after this model adjustment. The
model is a linear-logit description, not a validated propensity model; the
intervals do not include injection, PSF, sky, or calibration uncertainty and
are not corrected for the several reported coefficients. They must not be
used as inverse-probability weights yet.

## Interpretation

The earlier layer audit found neighboring-source masks to be the main cause
of low 50–75 kpc coverage. This candidate-level analysis shows that the
resulting profile support depends strongly on cluster-centric radius in both
filters and also on richness in `i`. The measured ICL–satellite association
therefore faces a non-random radial support cut. A correlation computed only
from retained profiles can mix a physical radial trend with the changing
availability of measurements.

Do not correct this by simply removing source masks or lowering the floor.
The mask protects the diffuse-light estimate from neighbor flux, and the
coverage diagnostic is not an injection-calibrated selection function. The
next analysis should (1) fit a multivariable candidate-selection model with
cluster-block uncertainty, (2) test its calibration against injected/recovered
satellites and diffuse profiles, then (3) compare association estimates under
predeclared radial-support and completeness weighting. Any corrected estimate
must carry the selection uncertainty and remain separate from the current
raw-supported-sample result.

## Reproducibility

Astrafek script: [`analyze_icl_satellite_mask_selection.py`](analyze_icl_satellite_mask_selection.py),
SHA-256 `0dc191bcfcafe7a97af3527177203e25aeae3ff9727379136834e69579f7b919`.
Adjusted model script:
[`analyze_icl_satellite_mask_selection_adjusted.py`](analyze_icl_satellite_mask_selection_adjusted.py),
SHA-256 `b87356251be01ef120ceac75043ded9002f178ef2e08b87d51fc2091925fde07`.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_mask_selection_2026-11-05_01`

- Candidate input table SHA-256: `347de5a6c12f065227687a9127f9954904f912288413af2b980fd3a61f6f20e8`
- Summary SHA-256: `04aefeace05f58b5a5035ecc7d43de07d6dfeebba2f2399684a33c7314623876`
- Bootstrap seed: `20261105`; 5,000 cluster-block replicates.
- Adjusted model summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_mask_selection_adjusted_2026-11-05_02/summary.json`
- Adjusted summary SHA-256: `9fe70bf1da239c05a4bf72d6675918e1864b7186a458dcf0a7ec6f0ea74800af`.
