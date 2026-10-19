# Spatial-block uncertainty for the ICL scale-transfer test — 2026-10-20

## Result

The four-scale ICL injection grid was reanalysed with 5,000 paired connected
5-degree target/control block-bootstrap replicates per filter. The same blocks
were drawn jointly across the 100, 200, 300, and 500 kpc profiles. The primary
statistic is the median per-frame ratio of plane-corrected to raw ICL response;
the paired change is the median across frames of the scale-specific ratio
minus the same frame's 100 kpc ratio. The reported 95% percentile intervals
are spatial-block bootstrap intervals, not p-values.

| Filter | Scale length (kpc) | Median retention [95% block interval] | Paired median change from 100 kpc [95% block interval] |
|---|---:|---:|---:|
| i | 100 | 1.007 [0.954, 1.017] | reference |
| i | 200 | 0.912 [0.889, 0.929] | −0.064 [−0.123, −0.057] |
| i | 300 | 0.807 [0.782, 0.831] | −0.178 [−0.225, −0.158] |
| i | 500 | 0.644 [0.633, 0.656] | −0.349 [−0.407, −0.329] |
| i2 | 100 | 0.997 [0.982, 1.007] | reference |
| i2 | 200 | 0.918 [0.906, 0.929] | −0.075 [−0.089, −0.060] |
| i2 | 300 | 0.812 [0.803, 0.817] | −0.182 [−0.198, −0.166] |
| i2 | 500 | 0.633 [0.624, 0.643] | −0.362 [−0.386, −0.341] |

The denominator-free median response difference (`plane-corrected − raw`, in
response per unit injected amplitude) shows the same progression. For scale
lengths 100/200/300/500 kpc its medians and 95% block intervals are `0.002`
`[−0.008, 0.004]`, `−0.041 [−0.052, −0.028]`, `−0.114 [−0.122, −0.100]`, and
`−0.259 [−0.273, −0.248]` in `i`; for `i2` they are `−0.001 [−0.003, 0.002]`,
`−0.040 [−0.042, −0.034]`, `−0.112 [−0.116, −0.106]`, and
`−0.270 [−0.276, −0.265]`.

All paired median retention changes from 100 kpc are below zero in these
bootstrap intervals. This shows that, within this selected injection sample,
the downstream plane estimator's attenuation increases with diffuse scale
length after position-dependent PSF-core convolution and mask redetection.
It does not validate a correction to observed ICL profiles.

## Denominator audit and limits

One `i2` frame (cluster 4687) has a positive 100 kpc raw response of only
`0.00338`; its corrected response is `−0.01707`, giving a ratio of `−5.06`.
Including that ratio drives the arithmetic mean retention at 100 kpc to `0.855`
while the median is `0.997`. The median is therefore the primary ratio
statistic; the machine-readable output records the denominator and frame.
One `i` frame has nonpositive raw response at one or more scales and is
excluded from the common positive-response ratio sample, leaving 12/13 frames.
All 44 `i2` frames remain.

After this response-support selection, `i` contains 12 frames in 9 connected
blocks (largest block: 3 frames); `i2` contains 44 frames in 24 blocks (largest
block: 5). The bootstrap treats each connected block as the resampling unit.
The small `i` sample and approximate 5-degree spatial partition limit the
interval calibration. These intervals do not include pixel/coadd covariance,
HSC upstream sky processing, extended PSF wings, photon noise, extinction, or
passband uncertainty. No HSC query or download was made; `science_measurement`
remains false.

## Reproduction and provenance

The Astrafek implementation is
[`assess_icl_psf_diffuse_scale_block_bootstrap.py`](assess_icl_psf_diffuse_scale_block_bootstrap.py),
with three focused tests in
[`test_assess_icl_psf_diffuse_scale_block_bootstrap.py`](test_assess_icl_psf_diffuse_scale_block_bootstrap.py).
It verifies all four scale-run hashes and audits, joins the same 57 target and
matched-control positions, and retains identical block draws across scales.
The CPU calculation ran on `lageunha` with 5,000 replicates per filter.

The authoritative result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_scale_block_bootstrap_2026-10-20_02/summary.json`
(SHA-256 `d87f9f794296a1744de8cb4e92b6d6a317904d51f73e649a2aa49e6760dd1248`).
The earlier `_01` output is preserved; it reported mean ratios only and is
superseded because the near-zero raw-response denominator made its `i2`
100 kpc mean unstable.
