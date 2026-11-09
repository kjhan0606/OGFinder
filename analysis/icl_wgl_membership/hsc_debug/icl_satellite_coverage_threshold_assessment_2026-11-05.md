# ICL–satellite association sensitivity to satellite coverage floor

The 800 existing target images were remeasured with the annulus-specific
satellite pipeline at a 40% target/control coverage floor, retaining the
minimum of three azimuthal controls. This is a deliberate selection and
estimator sensitivity; the 60% floor remains the baseline. No images were
downloaded, and all outputs remain `science_measurement=false`.

## Support recovered

For the ICL 200–300 kpc and satellite 50–75 kpc intersection, the 60% run has
94 cluster/filter pairs (26 in `i`, 68 in `i2`), of which 23 are center
eligible under the existing central-mask flag. The 40% run has 204 pairs (64
in `i`, 140 in `i2`), of which 52 are center eligible. This is a 2.17× increase
in outer-annulus pairs and a 2.26× increase in center-eligible pairs. These are
cluster/filter pairs; the filters are disjoint in this particular supported
table.

| Satellite annulus (kpc) | `i` 40% support | `i2` 40% support | Total pairs | Center eligible |
|---|---:|---:|---:|---:|
| 5–10 | 66 | 144 | 210 | 54 |
| 10–15 | 66 | 142 | 208 | 54 |
| 15–20 | 66 | 144 | 210 | 55 |
| 20–30 | 67 | 144 | 211 | 55 |
| 30–40 | 65 | 142 | 207 | 53 |
| 40–50 | 64 | 144 | 208 | 53 |
| 50–75 | 64 | 140 | 204 | 52 |

The profiles remain annulus-specific, and the comparison does not make all
radial bins a common sample.

## Association comparison at 50–75 kpc

The estimator is the redshift- and richness-adjusted partial Pearson
correlation between ICL target-minus-control surface brightness at 200–300 kpc
and the cluster-mean satellite target-minus-azimuthal-control surface
brightness at 50–75 kpc. Intervals use 2,000 five-degree connected-block
bootstrap replicates.

| Filter | Floor | Pairs | Partial r (95% block interval) | Center eligible |
|---|---:|---:|---:|---:|
| `i` | 60% | 26 | 0.079 `[-0.187, 0.512]` | 8 |
| `i` | 40% | 64 | -0.151 `[-0.333, 0.162]` | 22 |
| `i2` | 60% | 68 | 0.293 `[0.028, 0.500]` | 15 |
| `i2` | 40% | 140 | 0.249 `[-0.002, 0.454]` | 30 |

![Support and association sensitivity](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_coverage_threshold_2026-11-05_02/coverage_threshold_radial_association.png)

The paired 40%-minus-60% correlation intervals are `[-0.725, 0.401]` in `i`
and `[-0.252, 0.169]` in `i2`; both include zero. The i2 point estimate
decreases from 0.293 to 0.249 when the expanded low-coverage sample is
included, while its interval reaches zero. This sensitivity does not
establish that the association is absent or present; it shows that the
current tentative i2 result is not stable enough to promote to a science
claim. Only 30 of the 140 low-floor i2 pairs are center eligible.

Restricting to the existing center-eligible flag further reduces the sample to
22 `i` and 30 `i2` pairs at the 40% floor. Their partial correlations are
`-0.187 [-0.508, 0.349]` and `0.242 [-0.024, 0.684]`, respectively. At the 60%
floor, only 8 `i` and 15 `i2` center-eligible pairs remain; the `i` interval
spans `[-1,1]`, and the `i2` interval is `[-0.081, 0.785]`. All include zero,
showing that center-eligible support is currently too small for a precise
association test.

For the 94 pairs supported at both floors, the mean 50–75 kpc satellite
surface-brightness shift (40% minus 60%) is `-1.236 nJy arcsec^-2` in `i` and
`-0.287 nJy arcsec^-2` in `i2`. Their paired block-bootstrap 95% intervals are
`[-5.349, 2.694]` and `[-1.211, 0.756]`, respectively. The change occurs
because the 40% estimator admits additional azimuthal controls; the altered
control mean is not just a sample-size change.

## Decision and limitations

The relaxed floor more than doubles the outer-annulus support and gives a
useful stability check, but it also admits partly masked annuli and changes
the control baseline. The i2 estimate remains positive but its 40% interval
includes zero, and the paired floor-change interval includes zero in each
filter. Do not replace the 60% result with the larger 40% sample or combine
the two floors. Keep both as an explicit systematic sensitivity until
coverage-dependent flux bias, mask edges, PSF wings, sky transfer, completeness,
and source-Poisson noise are calibrated. The center-eligible outer sample
remains only 52 pairs and needs a dedicated central-light/PSF analysis.

## Reproduction and provenance

The 40% full-sample rerun used the existing `measure_satellite_outer_profiles_per_bin.py`
script on the 800 verified images. Astrafek's
[`analyze_icl_satellite_coverage_threshold.py`](analyze_icl_satellite_coverage_threshold.py)
performs the pair joins, correlations, and paired spatial bootstrap. Run it with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_coverage_threshold.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_coverage_threshold_2026-11-05_02 \
  --replicas 2000 --seed 20261105
```

Archived comparison summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_coverage_threshold_2026-11-05_02/summary.json`.
The 40% satellite-profile summary is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_outer_profile_perbin_cov040_2026-11-05_01/summary.json`.
The comparison records input and output SHA-256 hashes, used no HSC network,
and leaves the science gate false.

Figure SHA-256: `c2641f029a72066fc279b890918eb63e638519bd6f8edadec931ad08e5f451ac`.
The plotting script is
[`plot_icl_satellite_coverage_threshold.py`](plot_icl_satellite_coverage_threshold.py).
