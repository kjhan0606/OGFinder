# ICL–satellite outer-annulus coverage loss

This offline audit separates missing 50–75 kpc annulus area into image edge,
base quality mask, SEP segmentation, neighboring-source masks, and fallback
CAMIRA-member circles. It uses 2,691 satellite candidates in the 224
cluster/filter frames with ICL 200–300 kpc support. It reads saved masks and
catalogue rows, but no science-image pixel values and makes no HSC request.
`science_measurement=false`.

## Results

Using unrounded pixel coverage, 1,171/2,691 candidates (43.5%) pass the 40%
outer-annulus floor. The sequential first-stage-at-which-coverage-falls-below-
40% counts are:

| First limiting layer | Candidates | Fraction of floor failures |
| --- | ---: | ---: |
| Image footprint | 0 | 0.0% |
| Base quality mask | 465 | 30.6% |
| SEP segmentation mask | 81 | 5.3% |
| Neighbor detection masks (Kron ellipses or invalid-shape circles) | 973 | 64.0% |
| Fallback unmatched-member circles | 1 | 0.1% |

There are 1,520 candidates below the floor. The source-related stages account
for 1,055 (69.4%) of those failures; the dominant term is neighboring-source
masking. In the `i` filter, 318/844 candidates pass and 360 fail first at the
neighbor detection masks. In `i2`, 853/1,847 pass and 613 fail first at that
stage. After target coverage and the minimum of three valid azimuthal controls,
919 profiles remain, consistent with the 919 saved 50–75 kpc profiles; these
span 204 of the 224 ICL-supported frames.

The mean fraction of expected annulus pixels removed at each sequential layer
is about 0.0003 by the footprint, 0.287 by the base quality mask, 0.097 by
segmentation, 0.268 by detection-neighbor masks, and 0.00009 by fallback
member circles. These are ordered, overlapping mask effects; they are not
independent error terms. Candidate failures are classified at the first stage
crossing the floor, so the table does not imply that only that mask affected a
candidate.

The initial audit stopped because a float64 geometry implementation differed
from the measurement's float32 pixel-center arithmetic at annulus boundaries.
The audit now uses the measurement arithmetic. Four candidates whose rounded
diagnostic displayed `0.400` have reconstructed coverage just below 0.4; final
threshold counts therefore use unrounded pixel values. One-candidate parity,
monotonic layer coverage, all 224 frame hashes, and all 2,691 candidate
reconstructions passed.

## Interpretation and next analysis

Image boundaries do not explain the loss in this 50–75 kpc annulus. Neighbor
detection masks are the main measured coverage bottleneck, followed by the
base quality mask. This is evidence of mask-driven selection, not evidence that
those masks are too conservative: removing a mask could admit real contaminant
light. The result is a selection/completeness diagnostic, not a flux-bias
correction or an ICL–satellite detection.

Next, quantify how the mask-driven pass fraction varies with satellite
cluster-centric radius, detection brightness/SNR, cluster redshift and
richness, with cluster-level resampling. Then use injection/recovery to test
whether any proposed mask adjustment recovers diffuse signal without admitting
neighbor flux. Keep the 40% diagnostic floor and the established 60% baseline
separate; do not reinterpret either as a corrected science measurement.

## Reproducibility

Astrafek script: [`audit_icl_satellite_outer_coverage_loss.py`](audit_icl_satellite_outer_coverage_loss.py),
SHA-256 `0d0ae7fd6e06c9a0e5f8e6acfc2dbd3324b7c6628a84899efefc24b3b735fc65`.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_outer_coverage_loss_2026-11-05_03`

- Candidate table SHA-256: `3d8d6fe8cb86ac5d54caa56d947eeffcf8212eac09efa3ea8adaec8714ea16bb`
- Frame table SHA-256: `2ff9dd4079379e965f86c4325f439f97694fc04d14e6ed33dd6ce8c0e890fab0`
- Summary SHA-256: `42788e0c1084fbfbe18c9b368ecee6944adf85fa0441f7f434942200034704d4`
