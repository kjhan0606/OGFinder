# ICL–satellite–weak-lensing overlap selection audit (2026-11-06)

## Scope

The selected ICL–satellite–WL result uses 27 i2 lenses, while the audited
ICL–satellite support sample has 68 i2 systems. This audit measures where that
attrition occurs and compares the selected and excluded systems using the
already saved pair, eligibility, and sky-block tables. It uses no new HSC data.

## Support attrition and sky coverage

| Eligibility stage | i2 lenses |
|---|---:|
| ICL–satellite supported | 68 |
| In the 620-lens matched WL parent | 60 |
| With target shape rows | 35 |
| With positive random-shape counts | 31 |
| Verified random-subtracted WL common support | 27 |

The 68 lenses occupy 33 saved connected 5-degree blocks. Only 9 blocks contain
both common-support and no-common-support lenses. Those mixed blocks contain
10 of the 27 selected lenses and 16 of the 41 excluded lenses. The remaining
17 selected and 25 excluded lenses lie in blocks that have no counterpart
group. This is a spatial positivity failure for transporting the 27-lens WL
overlap to the complete 68-lens ICL–satellite sample. A fitted propensity or
large inverse-probability weights cannot supply the missing cross-group
comparisons in those blocks.

## Covariate balance

The table reports the selected-minus-excluded standardized mean difference.
Labels were permuted only within each mixed sky block, preserving each block's
observed selection count. The four physical/covariate balance variables form
one max-T family (20,000 permutations). Source-pair count is shown separately
because it is a direct ingredient of WL eligibility.

| Variable | Selected mean | Excluded mean | Standardized difference | Restricted-block p | Four-variable max-T p |
|---|---:|---:|---:|---:|---:|
| Redshift | 0.457 | 0.452 | +0.070 | 0.677 | 1.000 |
| ln(richness) | 2.958 | 3.021 | −0.321 | 0.069 | 0.843 |
| ICL 200–300 kpc contrast | 3.314 | 4.705 | −0.173 | 0.799 | 0.998 |
| Satellite 50–75 kpc contrast | 1.227 | 2.914 | −0.257 | 0.911 | 0.952 |

The max-T analysis does not detect covariate imbalance, but only nine mixed
blocks inform the restricted permutation; it has little power and does not
establish missing-at-random. Target shape-pair counts differ sharply by
construction: the selected median is 4,342 versus zero for excluded lenses.
This confirms that source support is a major part of the selection mechanism,
not an independent physical result.

## Implication

The 27-lens triple-overlap correlations remain descriptive for the selected
footprint. Do not inverse-probability-weight that result to represent all 68
i2 ICL–satellite systems. A broader WL conclusion requires additional valid
shape/random support or a design with independent field coverage. Keep
science_measurement=false.

## Reproducibility

- Runner: audit_icl_wl_overlap_selection.py, SHA-256
  b2d474d5ecc234e3b3fcb9bea21739f72af720d895d730017c1dfba111895a94.
- Results: /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/.
- summary.json SHA-256:
  fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98.
- selection_balance.tsv SHA-256:
  abb280237413cddde93727ebde5e0b901191b7a4fcc1bc7210f0b31affcd0250.
- block_support.tsv SHA-256:
  07fececd4c438d20cc926f3016e93674aaccea8feeef2232cd92402cfd6d988c.
- All input table and overlap-audit hashes are stored in summary.json.
