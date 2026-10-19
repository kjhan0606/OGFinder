# Satellite radial ICL association: family-wise intervals — 2026-10-20

## Result

The earlier seven-bin radial comparisons used ordinary 95% bootstrap intervals
and did not account for searching across radial bins. This follow-up resampled
the same connected 5-degree target/control sky blocks jointly across all seven
bins, separately for each filter. It calculated simultaneous max-T 95%
intervals for raw ICL association, outer-subtracted ICL association, and the
paired change between them. The three association families and the two filters
were kept separate.

**All 42 simultaneous intervals include zero.** The marginal i-band result at
5–10 kpc for outer-subtracted ICL (`rho=0.212`) no longer excludes zero after
family-wise adjustment (`[-0.054, 0.479]`). The largest positive estimates
remain exploratory and do not establish a radial ICL–satellite association.

### i filter (`n=77` clusters, 21 connected blocks)

| Satellite radius (kpc) | Raw ICL rho [simultaneous 95% interval] | Outer-subtracted ICL rho [simultaneous 95% interval] | Change: corrected − raw rho [simultaneous 95% interval] |
|---:|---:|---:|---:|
| 5–10 | 0.069 [−0.276, 0.414] | 0.212 [−0.054, 0.479] | 0.144 [−0.104, 0.391] |
| 10–15 | 0.031 [−0.320, 0.382] | 0.202 [−0.057, 0.461] | 0.171 [−0.077, 0.418] |
| 15–20 | −0.014 [−0.353, 0.327] | 0.144 [−0.096, 0.384] | 0.158 [−0.111, 0.427] |
| 20–30 | −0.017 [−0.346, 0.313] | 0.145 [−0.079, 0.369] | 0.162 [−0.100, 0.423] |
| 30–40 | −0.035 [−0.331, 0.261] | 0.125 [−0.114, 0.365] | 0.160 [−0.114, 0.434] |
| 40–50 | −0.030 [−0.305, 0.246] | 0.015 [−0.238, 0.268] | 0.045 [−0.246, 0.335] |
| 50–75 | 0.001 [−0.317, 0.318] | 0.074 [−0.208, 0.357] | 0.074 [−0.206, 0.354] |

### i2 filter (`n=183` clusters, 44 connected blocks)

| Satellite radius (kpc) | Raw ICL rho [simultaneous 95% interval] | Outer-subtracted ICL rho [simultaneous 95% interval] | Change: corrected − raw rho [simultaneous 95% interval] |
|---:|---:|---:|---:|
| 5–10 | 0.053 [−0.140, 0.246] | 0.011 [−0.149, 0.171] | −0.042 [−0.221, 0.138] |
| 10–15 | 0.049 [−0.149, 0.247] | 0.026 [−0.147, 0.199] | −0.023 [−0.209, 0.163] |
| 15–20 | 0.080 [−0.127, 0.287] | 0.041 [−0.138, 0.220] | −0.038 [−0.248, 0.171] |
| 20–30 | 0.102 [−0.098, 0.302] | 0.026 [−0.141, 0.193] | −0.076 [−0.287, 0.135] |
| 30–40 | 0.126 [−0.053, 0.306] | 0.074 [−0.097, 0.245] | −0.053 [−0.266, 0.161] |
| 40–50 | 0.089 [−0.108, 0.286] | 0.067 [−0.124, 0.258] | −0.022 [−0.225, 0.183] |
| 50–75 | 0.051 [−0.115, 0.218] | 0.086 [−0.102, 0.273] | 0.034 [−0.165, 0.234] |

Intervals are descriptive uncertainty bands from 5,000 connected-block
bootstrap replicates, not p-values. Resampling uses the same blocks jointly
across the seven radial bins within each filter; ties and duplicate blocks are
retained as produced by the bootstrap. The max-T critical value is estimated
from the largest absolute standardized deviation across the seven bins.

## Interpretation and limits

The family-wise adjustment removes the apparent marginal i-band 5–10 kpc
exclusion of zero. Across both bands, the adjusted intervals remain broad,
including for the paired corrected-minus-raw comparison. The results support
no detection claim and do not show that outer subtraction recovers a physical
ICL component.

The 400–500 kpc annulus is not verified blank sky. Satellites differ between
radial bins even though the cluster support is fixed. The block partition is a
chosen spatial dependence approximation and does not include pixel/coadd
covariance. PSF wings, sky transfer, completeness, extinction and photometric
systematics also remain unresolved. `science_measurement=false` remains in
force.

## Reproducibility

The analysis uses the frozen cluster-level table from
[`satellite_radial_icl_baseline_association_assessment_2026-10-20.md`](satellite_radial_icl_baseline_association_assessment_2026-10-20.md)
and is implemented in
[`assess_satellite_radial_icl_association_familywise.py`](assess_satellite_radial_icl_association_familywise.py).
Three focused tests pass in
[`test_assess_satellite_radial_icl_association_familywise.py`](test_assess_satellite_radial_icl_association_familywise.py).
The CPU run used `lageunha`, with 5,000 replicates per filter and fixed seeds
20261020 and 20271020. It made no HSC request.

The machine-readable output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_radial_icl_familywise_2026-10-20_01/summary.json`.
Its SHA-256 is
`f824325eefb9ddaa1792112d8e645ca6226a642ac0b7cfc98a0eb6e691a57de8`.
The output records source and code hashes, bootstrap settings and all
simultaneous intervals.
