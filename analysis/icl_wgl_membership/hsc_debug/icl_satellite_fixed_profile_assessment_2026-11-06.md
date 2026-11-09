# Fixed-support radial ICL–satellite relation (2026-11-06)

## Design

The earlier radial relation analysis used the available profile rows at each
radius, so the contributing clusters changed with radius. This analysis fixes
the sample first: within each filter, a cluster must belong to the saved
94-pair ICL/satellite support table and have valid satellite excess in all
seven bins from 5 to 75 proper kpc. Every radial coefficient then uses that
same filter-specific cluster set.

The model is satellite surface-brightness excess versus the 200–300 kpc ICL
contrast, standardized redshift, and standardized log richness. Inference uses
CR1 uncertainty by connected 5-degree target/control sky blocks and a
restricted HC2 Webb wild-cluster bootstrap. Holm and single-step max-T
corrections cover all 14 filter/radius comparisons.

| Filter | Fixed clusters | Connected blocks | ICL slope at 50–75 kpc (95% CR1 interval) | Raw bootstrap p | 14-test max-T p |
|---|---:|---:|---:|---:|---:|
| `i` | 21 | 13 | +0.99 [−2.31, +4.29] | 0.508 | 1.000 |
| `i2` | 53 | 28 | +1.27 [+0.29, +2.24] | 0.041 | 0.410 |

The outer `i2` pointwise interval is positive, but its max-T-adjusted p-value is
0.410; across the full family, the smallest max-T p-value is 0.237. The `i`
inner slopes are negative and their first three pointwise intervals exclude
zero, but none survives the familywise correction. Thus the same-cluster
radial comparison does not establish an ICL–satellite association at any
radius. It removes changing cluster membership across the seven bins, but not
member-detection or mask changes within a cluster.

## Artifacts and limits

- Runner: [`analyze_icl_satellite_fixed_profile.py`](analyze_icl_satellite_fixed_profile.py)
- Runner SHA-256: `b19d18de0b54e4918463b55232ac260e6b06e109e585475ed42511cbdcde187c`
- Result summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_fixed_profile_2026-11-06_02/summary.json`
- Summary SHA-256: `08e60ec574a624f4c177c7fbbc66c9ee93a57a4072cbea1755acf23aaeda113d`
- Results table SHA-256: `23e7d9a26d2f070787f52ce011917b67e8ad13ddfde410738f3dc84d196f9efb`
- Figure: `icl_satellite_fixed_sample_radial_slopes.png` and `.pdf` in the result directory.
- Bootstrap: 10,000 replicates, seed `20261110`; no HSC network or new image download.

The 94-pair parent was selected by earlier support gates. Requiring all radial
bins leaves 21 `i` and 53 `i2` clusters, so this is not representative of all
800 images. Five-degree connected cells approximate HSC survey fields. The
errors omit source Poisson, full coadd-weight, sky-transfer, extended-PSF,
completeness, and injection-response uncertainty. The familywise adjustment
does not cover the broader history of exploratory analyses. Keep
`science_measurement=false`.
