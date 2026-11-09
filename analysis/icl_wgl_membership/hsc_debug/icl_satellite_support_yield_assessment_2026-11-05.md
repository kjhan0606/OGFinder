# ICL–satellite support yield in the 800-parent sample

This offline audit joins the saved 200–300 kpc ICL support table to each
annulus-specific satellite cluster profile. It reads no image pixels and
requests no HSC data. It quantifies where the current parent sample loses
cluster/filter pairs; it does not change the measurement cuts or the
`science_measurement=false` gate.

## ICL and satellite radial support

There are 224 ICL-supported cluster/filter pairs: 71 in `i` and 153 in `i2`.
The number that also has a satellite profile decreases with satellite radius:

| Satellite annulus (kpc) | `i` intersection | `i2` intersection | Total | Center-eligible total |
|---|---:|---:|---:|---:|
| 5–10 | 62 | 136 | 198 | 52 |
| 10–15 | 60 | 133 | 193 | 50 |
| 15–20 | 52 | 125 | 177 | 45 |
| 20–30 | 48 | 109 | 157 | 41 |
| 30–40 | 44 | 94 | 138 | 39 |
| 40–50 | 32 | 86 | 118 | 31 |
| 50–75 | 26 | 68 | 94 | 23 |

At the primary 50–75 kpc outer annulus, the intersection is 26 `i` and 68
`i2` pairs. Only 8 `i` and 15 `i2` pairs are center-eligible under the saved
Kron-2.5 central mask. Thus the unmasked-center subset is especially small.
Inner satellite bins have more support, but they also have greater exposure
to central-light, ICL-overlap, and PSF systematics; their larger counts do not
replace the outer-structure measurement.

## Candidate retention

The per-annulus pipeline records 9,468 satellite candidates from the existing
800-parent images. The saved status counts are:

| Status | Candidates | Fraction |
|---|---:|---:|
| No target annulus reaches the 60% coverage floor | 6,383 | 67.4% |
| Target annulus support exists, but fewer than three azimuthal controls pass | 1,306 | 13.8% |
| At least one annulus retained | 1,768 | 18.7% |
| Invalid target shape | 11 | 0.1% |

The first status is assigned before control positions are generated: no target
annulus has coverage at or above 0.60. Of these candidates, 994 have at least
one annulus with coverage at or above 0.40, so a looser cut could admit some
additional measurements. That would change the estimator and could amplify
mask-edge and background bias; this audit does not recommend relaxing it. The
current tables do not distinguish whether the lost target coverage comes
from image boundaries, source/neighbor masks, or detection geometry.

The 50–75 kpc `i2` intersection is 68 pairs, or 8.5% of the 800 parent IDs.
If that selected-sample yield persisted, 323–485 usable `i2` pairs would
require roughly 3,800–5,706 parent clusters. This is only a planning
extrapolation from the exploratory `r=0.293` point estimate; it is not a
confirmatory effect-size assumption or a claim that more data alone will
remove the existing selection and image-systematic biases.

## Next analysis implication

Before enlarging the HSC sample, investigate the 60% target-annulus coverage
loss using the saved masks, member/detection geometry, and cutout WCS. If the
dominant cause is a cutout edge, a separately planned wider cutout may recover
support. If masks or source confusion dominate, a wider cutout will not fix
the measurement. Any altered mask, control placement, or coverage threshold
requires a separately versioned estimator, a paired-support comparison, and
new noise/completeness/PSF checks. Do not merge different support samples into
a single radial trend.

## Reproduction and provenance

Astrafek script: [`audit_icl_satellite_support_yield.py`](audit_icl_satellite_support_yield.py)

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_icl_satellite_support_yield.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_support_yield_2026-11-05_01
```

Archived output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_support_yield_2026-11-05_01/summary.json`.
The summary includes SHA-256 hashes for both upstream support products, the
candidate diagnostics, and both generated tables. The output status is
`complete_offline_diagnostic`, with `hsc_network_used=false` and
`science_measurement=false`.
