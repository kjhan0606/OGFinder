# ICL–satellite association: block-width sensitivity of the 94-pair support sample

## Question and method

The existing ICL–satellite association table has 94 pairs with both a 200–300
kpc ICL profile and a 50–75 kpc satellite profile: 26 in `i`, 68 in `i2`.
This Astrafek diagnostic re-estimates the raw target-minus-control partial
Pearson correlation while controlling linearly for redshift and CAMIRA
richness. It compares the full 94-pair support subset with the nested earlier
74-pair diagnostic subset at 2.5°, 5°, and 10° equal-area cell widths. Pairs
whose target or matched-control endpoint shares a cell are joined into the
same connected block. Each bootstrap resamples those blocks with replacement;
the 94-minus-74 change uses identical 94-sample block resamples for both
nested scopes. Each interval uses 5,000 draws.

No image pixels or new HSC data were read or requested. The analysis retains
`science_measurement=false`.

## Results

Values are partial `r` with percentile 95% block-bootstrap intervals.

| Filter | Pairs | 2.5° blocks | 5° blocks | 10° blocks |
|---|---:|---:|---:|---:|
| `i`, prior subset | 21 | 0.120 `[-0.404, 0.659]` | 0.120 `[-0.218, 0.527]` | 0.120 `[-0.547, 0.533]` |
| `i`, full support | 26 | 0.079 `[-0.348, 0.585]` | 0.079 `[-0.190, 0.495]` | 0.079 `[-0.490, 0.499]` |
| `i2`, prior subset | 53 | 0.234 `[-0.011, 0.399]` | 0.234 `[0.017, 0.418]` | 0.234 `[0.004, 0.393]` |
| `i2`, full support | 68 | 0.293 `[0.002, 0.504]` | 0.293 `[0.048, 0.498]` | 0.293 `[0.056, 0.496]` |

The point estimates are invariant to block width, while uncertainty changes.
The `i` intervals are broad and cross zero at every width. The full-support
`i2` interval is only just above zero at 2.5° and remains positive at 5° and
10°. The prior 53-pair `i2` interval crosses zero at 2.5° but not at the two
coarser widths. This dependence on sample scope and block construction does
not establish a robust population signal.

The paired change from 74 to 94 pairs is `−0.041` for `i` and `+0.059` for
`i2`; every paired interval includes zero:

| Filter | 2.5° paired Δr CI | 5° paired Δr CI | 10° paired Δr CI |
|---|---:|---:|---:|
| `i` | `[-0.253, 0.136]` | `[-0.240, 0.089]` | `[-0.201, 0.152]` |
| `i2` | `[-0.191, 0.351]` | `[-0.187, 0.337]` | `[-0.203, 0.367]` |

![Block-width sensitivity of the 74- and 94-pair partial correlations](</gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_supported_94_blockwidth_2026-11-05_01/icl_satellite_94_blockwidth_clean.png>)

## Interpretation and limits

The positive `i2` point estimate merits follow-up, but the 94 pairs are a
selected subset of the 800-cluster image parent and show measurable support
selection. The intervals describe this finite selected sample under
approximate equal-area blocks; these blocks are not validated contiguous HSC
footprint regions. The paired scope-change intervals do not show that adding
the 20 pairs changes the association. The point estimate is unweighted and
does not correct for image covariance, injected-source Poisson noise, PSF
wings, sky transfer, extinction, source completeness, or the probability of
satellite-profile support. No causal interpretation or detection claim is
made.

The next scientific gate remains broader injection/recovery and image-noise
validation, with support selection modeled against the 800-cluster parent.
The existing 2,005-position PSF request is still a separate unsubmitted plan
and requires the user's explicit approval before HSC access.

## Reproducibility

- Runner: `analyze_icl_satellite_supported_94_blockwidth.py`, SHA-256
  `dfd845145c3e5c752aa9bf8ab0e0f77d7e99ed0bcd8364b197de51812388945b`.
- Plotter: `plot_icl_satellite_supported_94_blockwidth.py`, SHA-256
  `5115445514c12ff47d32edbcc4cb54f863c632bf83b436ff6c2b44ceb86e83cf`.
- Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_supported_94_blockwidth_2026-11-05_01/summary.json`, SHA-256
  `7aad81051d46b88804db2923cf0de38e85e62384fa8d82455839b827e96c7227`.
- Figure: `icl_satellite_94_blockwidth_clean.png`, SHA-256
  `b01765d7a6aab05bf3e60e5ea9f518f65317bad7c128b8cb476890fb7ec45a93`.
