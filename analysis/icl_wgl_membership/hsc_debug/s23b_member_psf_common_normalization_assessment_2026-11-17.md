# Member PSF-core response with common input normalization

## Purpose and reproducibility

The prior paired comparison divides each variant's recovered flux by **its own**
nominal 50–75 proper-kpc input. PSF convolution changes that denominator, so a
small difference between the two normalized responses can hide a change in the
recovered flux. This offline analysis uses the original completed 427-member ×
three-amplitude ledgers and the exact both-supported selection. It verifies
pinned summaries, ledger hashes, the paired keys, baseline replay, support
counts and stored flux/denominator arithmetic. It reads no images and makes no
network requests. All 1,218 matched member-amplitude pairs supported by both
variants enter the corresponding conditional samples. `science_measurement=false`.

For each matched row, with recovered increments F0 (unconvolved), F1 (PSF
core), and nominal annular inputs D0 and D1, the exact decomposition is

`F1/D1 - F0/D0 = (F1-F0)/D0 + (F1/D1 - F1/D0)`.

The first term uses the same original input D0 for both variants; the second
is the change due to using D1 for the convolved variant. This is an algebraic
accounting of these simulations, not an observational PSF correction.

## Results at the 1.0 nJy arcsec^-2 injection amplitude

The 95% intervals resample the same target sky cells with replacement and are
marginal, selected-sample, exploratory intervals. Numbers are medians of
within-member quantities. Medians of the two decomposition terms need not add
to the median total, even though the identity holds for every row.

| Quantity | i: 58 pairs, 20 cells | i2: 351 pairs, 61 cells |
|---|---:|---:|
| Fractional nominal-input change D1/D0 − 1 | +0.00923 [0.00866, 0.00973] | +0.01014 [0.00966, 0.01057] |
| Recovered-flux change / common D0 | +0.00784 [0.00547, 0.00976] | +0.00906 [0.00828, 0.01001] |
| Denominator contribution | −0.00742 [−0.00866, −0.00383] | −0.00766 [−0.00839, −0.00695] |
| Original different-denominator ratio change | −0.00146 [−0.00581, 0.00196] | +0.00200 [0.00083, 0.00329] |

Thus the changed nominal denominator offsets much of the common-denominator
recovery change in this particular experiment. A near-zero original ratio
change does not establish an absence of PSF effects. The observed differences
are small on this chosen normalization, but no general bound follows for other
profiles, radii or samples. Selection changed in several cases; excluding those
cases prevents completeness or population inference. The i/i2 labels are kept
separate and should not be treated as independent survey realizations.

Only finite local PSF cores and one circular exponential halo are tested.
Extended PSF wings, source Poisson, sky transfer, coadd weights, correlated
noise, membership and sample selection remain unresolved. Intervals do not
correct for multiple comparisons. No flux or science measurement was promoted.

- Analysis: `analyze_s23b_member_psf_common_normalization.py`
- Unit tests: `test_s23b_member_psf_common_normalization.py`
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_psf_common_normalization_2026-11-17_01/`
- Summary SHA-256: `cfb3d8fc5422fad350939ec762b5e58533efc731209eebedfa2af784069e27a8`
- Source and comparison summary hashes are embedded in the output. The full
  per-row decomposition and sky-cell bootstrap results are retained.
