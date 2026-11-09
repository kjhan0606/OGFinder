# Fixed-sample high/low ICL satellite profiles (2026-11-06)

## Design and results

Using the same complete seven-bin profile sample as the fixed-slope analysis,
Astrafek split clusters at the within-filter median 200–300 kpc ICL contrast.
The groups are fixed across radius: `i` has 10 low/11 high ICL clusters, and
`i2` has 26 low/27 high. At each radius, the adjusted high-minus-low satellite
excess comes from an OLS model controlling standardized redshift and log
richness. CR1 uncertainty uses joint target/control connected 5-degree blocks;
restricted HC2 Webb bootstrap and max-T correction cover all 14 tests.

| Filter | Adjusted high–low at 50–75 kpc (95% CR1 interval) | Raw bootstrap p | 14-test max-T p |
|---|---:|---:|---:|
| `i` | +4.61 [−3.05, +12.27] nJy arcsec⁻² | 0.232 | 0.896 |
| `i2` | +1.35 [−1.23, +3.93] nJy arcsec⁻² | 0.309 | 0.967 |

The smallest adjusted p across the entire family is 0.452. The positive outer
`i2` continuous slope therefore does not appear as a resolved high-versus-low
profile contrast. No radial bin passes the familywise correction. This result
uses the same selected complete-profile support as the continuous-slope test,
so it is a complementary estimand, not an independent replication.

## Reproducibility and limits

- Runner: [`analyze_icl_satellite_fixed_profile_split.py`](analyze_icl_satellite_fixed_profile_split.py)
- Runner SHA-256: `79a54e88e03d1e2827701d276264a368e4da149a2d5b74afd4fb2aa7bd382d30`
- Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_fixed_profile_split_2026-11-06_02/summary.json`
- Summary SHA-256: `49403d38eb267f61bee9812c7b2a24694a2fa45b5463fb2db60c13cbd4d1f909`
- Results table SHA-256: `ca54827661f5ab1723cdd511685db0300333620aa3d87d27adfbc266c9cf30f6`
- Figure SHA-256: `58177773b0e7dfc95fc6a6845aaf5c17510104a7473fec08d32caf450a077034`
- Bootstrap: 10,000 replicates, seed `20261111`; no HSC network or image download.

The complete-profile sample retains only 21 `i` and 53 `i2` pairs from a
selected 94-pair support intersection. Five-degree connected cells approximate
survey fields. The errors omit source Poisson, coadd weights, sky transfer,
extended PSF, completeness, and injection-response uncertainties. Keep
`science_measurement=false`.
